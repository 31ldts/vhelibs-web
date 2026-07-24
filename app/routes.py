# -*- coding: utf-8 -*-
"""
Flask routes for VHELIBS web.

Endpoints:
  GET  /                        – serve the SPA
  POST /api/analyse             – start an analysis job
  GET  /api/status/<job>        – poll job progress / results
  GET  /api/edm/<pdbid>         – full density map, downloaded/cached on demand.
                                   Default: CCP4 2Fo-Fc map from RCSB/EBI.
                                   ?source=pdb_redo: PDB-REDO's own .map file
                                   (map-maker service) — used instead of
                                   density-box for PDB-REDO analyses, since
                                   EBI's density server only has maps for the
                                   original (non-PDB-REDO) entry.
  GET  /api/density-box/<pdbid> – cropped density chunk (BinaryCIF) for a given region,
                                   proxied + cached from EBI's density server.
                                   Only valid for the original PDB entry, not
                                   PDB-REDO — see /api/edm above for that case.
  GET  /api/density-mask/<pdbid>/<region> – pre-masked CCP4 map for one region
                                   (ligand/binding_site/residues_to_examine),
                                   cropped + masked server-side via
                                   core.density_mask/gemmi so the 3D viewer
                                   renders one isosurface per region instead
                                   of one per atom.
  POST /api/cache/clear         – delete every file under the on-disk cache
                                   (structures, validation stats, density
                                   maps, masked-map cache, UniProt lookups).
                                   Nothing in the cache is ever removed
                                   automatically otherwise.
"""
import os
import re
import shutil
import threading
import uuid
import logging
import concurrent.futures
from functools import wraps
from urllib.parse import urlencode

import requests
from flask import (
    Blueprint, render_template, request, jsonify, current_app,
    Response, abort, send_file,
)

import core.pdb_utils as pdb_utils
import core.eds_utils as eds_utils
import core.pdb_redo_utils as pdb_redo_utils
import core.cofactors as cofactors
import core.http_cache as http_cache
import core.density_mask as density_mask
from core.rsr_core import AnalysisConfig, analyse_pdbids

logger = logging.getLogger(__name__)
bp = Blueprint("main", __name__)

# In-memory job store  {job_id: {"status": ..., "progress": ..., "results": ...}}
_jobs = {}
_jobs_lock = threading.Lock()

# Shared background executor used to fire-and-forget the default-radius
# density-mask precompute for every ligand result as soon as a job
# finishes (see _prefetch_density_masks below and
# core.density_mask.prefetch_default_masks), so that by the time a user
# opens the 3D viewer the default view is usually already cached. Kept
# here (web layer) rather than inside core.rsr_core deliberately: that
# module is exercised directly by tests/test_rsr_core.py against real
# core.eds_utils/core.pdb_redo_utils calls, and firing background network
# I/O from inside it would make those tests do real network calls. This
# executor, and the whole prefetch step, is best-effort — never allowed
# to affect the job result itself.
_density_prefetch_executor = concurrent.futures.ThreadPoolExecutor(
    max_workers=2, thread_name_prefix="density-mask-prefetch"
)

# Canonical UniProt ID pattern (6-char, or 10-char introduced in
# later releases), e.g. "P00734", "A0A0A0MRZ7". Anything that doesn't match
# this is treated as a plain PDB ID.
UNIPROT_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$",
    re.IGNORECASE,
)

# UniProt "entry name" / mnemonic, e.g. "PPARG_HUMAN", "1433B_HUMAN": a
# 1-5 character protein/gene code, an underscore, then a 1-5 character
# organism mnemonic (https://www.uniprot.org/help/entry_name). These are
# not stable identifiers the way accessions are (they can be reassigned if
# a gene/organism is reclassified), so tokens matching this are resolved
# down to their current primary accession in _resolve_uniprot_entry_name
# before anything else touches them — accession is what actually gets
# cached/used by core.pdb_utils.get_pdbids_for_uniprot.
UNIPROT_ENTRY_NAME_RE = re.compile(r"^[A-Z0-9]{1,5}_[A-Z0-9]{1,5}$", re.IGNORECASE)

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def _resolve_uniprot_entry_name(entry_name):
    """
    Resolve a UniProt entry name / mnemonic (e.g. "PPARG_HUMAN") to its
    current primary accession (e.g. "P37231") via the UniProt REST search
    API, caching the (small) JSON response to disk the same way the other
    core.*_utils modules do via core.http_cache.

    Returns the primary accession as a string, or None if no UniProt entry
    matches ``entry_name``.
    """
    entry_name = entry_name.upper()
    cache_dir = http_cache.entry_cache_dir(pdb_utils.CACHEDIR, "_uniprot_names")
    cache_path = os.path.join(cache_dir, entry_name + ".json")

    qs = urlencode({"query": f"id:{entry_name}", "fields": "accession", "format": "json"})
    url = f"{UNIPROT_SEARCH_URL}?{qs}"

    data = http_cache.fetch_json(url, cache_path)
    if not data:
        return None
    results = data.get("results") or []
    if not results:
        return None
    return results[0].get("primaryAccession")


def _with_cache_dir(view):
    """Point core.pdb_utils' shared cache dir at the app's configured one.

    Every route that ends up calling into core.pdb_utils/eds_utils/
    pdb_redo_utils (which all read the module-level CACHEDIR set by
    core.pdb_utils.set_cache_dir) needs this done first. Applied as a
    decorator instead of repeating the same call in each view.
    """
    @wraps(view)
    def wrapper(*args, **kwargs):
        pdb_utils.set_cache_dir(current_app.config["CACHE_DIR"])
        return view(*args, **kwargs)
    return wrapper


def _expand_ids(tokens):
    """
    Expand any UniProt ID or entry names/mnemonics found in
    *tokens* into the PDB IDs of the structures that reference them
    (core.pdb_utils.get_pdbids_for_uniprot). Entry names (e.g.
    "PPARG_HUMAN") are first resolved to their primary accession via
    _resolve_uniprot_entry_name; accessions (e.g. "P37231") are used as-is.
    Plain PDB IDs pass through unchanged. De-duplicates while preserving
    order of first appearance.

    Returns (pdbids, origin_map, unresolved):
      pdbids      – flat list of PDB IDs to analyse
      origin_map  – {pdbid.lower(): uniprot_accession_or_None}, so results
                    can be tagged with the UniProt code they came from
      unresolved  – UniProt-looking tokens that returned no PDB entries
                    (or, for entry names, no matching UniProt ID)
    """
    pdbids = []
    origin_map = {}
    unresolved = []
    seen = set()

    for token in tokens:
        uniprot_id = None
        if UNIPROT_RE.match(token):
            uniprot_id = token.upper()
        elif UNIPROT_ENTRY_NAME_RE.match(token):
            uniprot_id = _resolve_uniprot_entry_name(token)
            if not uniprot_id:
                unresolved.append(token.upper())
                continue

        if uniprot_id:
            related = pdb_utils.get_pdbids_for_uniprot(uniprot_id)
            if not related:
                unresolved.append(uniprot_id)
                continue
            for pid in related:
                key = pid.lower()
                if key not in seen:
                    seen.add(key)
                    pdbids.append(pid)
                origin_map.setdefault(key, uniprot_id)
        else:
            key = token.lower()
            if key not in seen:
                seen.add(key)
                pdbids.append(token)
            origin_map.setdefault(key, None)

    return pdbids, origin_map, unresolved


def _prefetch_density_masks(pdbid, result):
    """Fire-and-forget precompute of the default-radius masked density
    maps for every ligand in one analysis result (see
    core.density_mask.prefetch_default_masks).
    """
    for ligand in (result or {}).get("ligands") or []:
        boxes = ligand.get("density_boxes")
        atoms = ligand.get("density_atoms")
        if not boxes or not atoms:
            continue
        try:
            density_mask.prefetch_default_masks(
                pdbid, boxes, atoms, executor=_density_prefetch_executor,
            )
        except Exception:
            logger.exception("Density-mask prefetch failed for %s (non-fatal)", pdbid)


def _run_job(job_id, pdbids, cfg, origin_map=None):
    # The actual concurrent analysis (thread pool, per-entry error
    # handling, order preservation) lives in core.rsr_core.analyse_pdbids
    # — see that function's docstring for why a thread pool is safe here.
    # This wrapper only adds the two things that are specific to *this*
    # web job: live progress updates for /api/status polling, and tagging
    # each result with the UniProt ID (if any) that produced it.
    origin_map = origin_map or {}

    def _on_progress(completed, total):
        # analyse_pdbids calls this from its own background thread (the
        # one _run_job itself runs in), never concurrently with itself,
        # so only the shared `_jobs` dict needs locking here — it's also
        # read from request-handling threads via /api/status/<job_id>.
        with _jobs_lock:
            _jobs[job_id]["progress"] = completed
            _jobs[job_id]["total"] = total

    results = analyse_pdbids(pdbids, cfg, on_progress=_on_progress)

    for pdbid, res in zip(pdbids, results):
        # Tag with the UniProt ID that produced this PDB ID, if any,
        # so the frontend can group/label results accordingly.
        res["uniprot"] = origin_map.get(pdbid.strip().lower())
        # Deliberately done here rather than inside core.rsr_core, so the
        # analysis engine itself stays free of this web-layer side effect
        _prefetch_density_masks(pdbid, res)

    with _jobs_lock:
        _jobs[job_id]["status"] = "done"
        _jobs[job_id]["results"] = results


@bp.route("/")
def index():
    return render_template("index.html")


@bp.route("/api/blacklist")
def blacklist_defaults():
    """
    Return the built-in ligand-blacklist/metal entries (core.cofactors), so
    the Analysis tab can render them as togglable checkboxes instead of
    them being invisible/hardcoded. The response never reflects any single
    user's customization — per-request overrides (unchecked entries, custom
    additions, or a fully replaced list from an uploaded file) are sent back
    to the server with each /api/analyse call instead, see that route.
    """
    return jsonify({"entries": cofactors.get_default_entries()})


@bp.route("/api/blacklist/parse", methods=["POST"])
def blacklist_parse():
    """
    Parse an uploaded blacklist file (sent as raw text) and return it as
    structured entries, so the frontend can show a preview ("this file
    defines N blacklist + M metal entries") before the user commits to
    replacing their current list with it.
    """
    data = request.get_json(force=True, silent=True) or {}
    text = data.get("text", "")
    if not isinstance(text, str) or not text.strip():
        return jsonify({"error": "The file appears to be empty."}), 400

    metals_dict, blacklist_dict = cofactors.parse_uploaded_list(text)
    if not metals_dict and not blacklist_dict:
        return jsonify({"error": "No valid entries found in the file."}), 400

    entries = (
        [{"code": c, "name": n, "category": "blacklist"} for c, n in blacklist_dict.items()]
        + [{"code": c, "name": n, "category": "metal"} for c, n in metals_dict.items()]
    )
    entries.sort(key=lambda e: (e["category"], e["code"]))
    return jsonify({
        "entries": entries,
        "metals": metals_dict,
        "ligand_blacklist": blacklist_dict,
    })


@bp.route("/api/analyse", methods=["POST"])
@_with_cache_dir  # UniProt resolution below also caches to disk, so set this first
def analyse():
    data = request.get_json(force=True, silent=True) or {}

    raw_ids = data.get("pdbids", "")
    if isinstance(raw_ids, list):
        tokens = [p.strip() for p in raw_ids if p.strip()]
    else:
        tokens = [p.strip() for p in str(raw_ids).replace(",", "\n").split() if p.strip()]

    if not tokens:
        return jsonify({"error": "No PDB IDs or UniProt IDs provided"}), 400

    pdbids, origin_map, unresolved = _expand_ids(tokens)

    if not pdbids:
        msg = "No valid PDB IDs found."
        if unresolved:
            msg += " Could not find any PDB entries for UniProt ID(s): " + ", ".join(unresolved)
        return jsonify({"error": msg}), 400

    # Build config from request params
    def _f(key, default):
        v = data.get(key)
        return float(v) if v not in (None, "") else default

    def _i(key, default):
        v = data.get(key)
        return int(v) if v not in (None, "") else default

    def _b(key):
        return bool(data.get(key, False))

    # Per-request blacklist customization. Built fresh for this job only
    # — never mutates core.cofactors' shared module-level dicts — so 
    # concurrent jobs from other users/tabs are unaffected.
    blacklist_cfg = data.get("blacklist") or {}
    effective_metals, effective_ligand_blacklist = cofactors.build_effective_lists(
        disabled_codes=blacklist_cfg.get("disabled"),
        custom_entries=blacklist_cfg.get("custom"),
        replace=blacklist_cfg.get("replace"),
    )

    cfg = AnalysisConfig(
        rsr_upper=_f("rsr_upper", 0.4),
        rsr_lower=_f("rsr_lower", 0.24),
        rscc_min=_f("rscc_min", 0.9),
        rfree_max=_f("rfree_max", 1.0),
        occupancy_min=_f("occupancy_min", 1.0),
        tolerance=_i("tolerance", 2),
        distance=_f("distance", 4.5),
        pdb_redo=_b("use_pdb_redo"),
        check_owab=_b("check_owab"),
        owab_max=_f("owab_max", 50.0),
        check_resolution=_b("check_resolution"),
        resolution_max=_f("resolution_max", 3.5),
        use_rdiff=_b("use_rdiff"),
        rdiff_max=_f("rdiff_max", 0.05),
        use_dpi=_b("use_dpi"),
        dpi_max=_f("dpi_max", 0.42),
        use_cache=bool(data.get("use_cache", True)),
        metals=effective_metals,
        ligand_blacklist=effective_ligand_blacklist,
    )

    job_id = str(uuid.uuid4())
    with _jobs_lock:
        _jobs[job_id] = {"status": "running", "progress": 0, "total": len(pdbids), "results": None}

    t = threading.Thread(target=_run_job, args=(job_id, pdbids, cfg, origin_map), daemon=True)
    t.start()

    response = {"job_id": job_id, "total": len(pdbids)}
    if unresolved:
        response["warnings"] = [
            f"No PDB entries found for UniProt ID {u}" for u in unresolved
        ]
    return jsonify(response)


@bp.route("/api/status/<job_id>")
def status(job_id):
    with _jobs_lock:
        job = _jobs.get(job_id)
    if not job:
        return jsonify({"error": "Unknown job"}), 404
    return jsonify({
        "status": job["status"],
        "progress": job["progress"],
        "total": job["total"],
        "results": job["results"],
    })


def _dir_size_and_count(path):
    """Total size (bytes) and file count of everything under `path`."""
    total_size = 0
    total_files = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            fp = os.path.join(root, name)
            try:
                total_size += os.path.getsize(fp)
                total_files += 1
            except OSError:
                pass
    return total_size, total_files


@bp.route("/api/cache/clear", methods=["POST"])
@_with_cache_dir
def clear_cache():
    """
    Delete everything under the server's on-disk cache: downloaded mmCIF
    structures, validation stats, electron-density maps, masked-map
    cache, UniProt lookups, etc. — every file ever written via
    core.http_cache or core.density_mask into CACHEDIR.

    Nothing in the cache is ever removed automatically, so this is the
    only way to reclaim that disk space.

    Individual entries that fail to delete (e.g. a permissions error, or
    a file locked by another process) are skipped and reported in
    ``errors`` rather than aborting the whole cleanup.
    """
    cache_dir = current_app.config["CACHE_DIR"]
    removed_files = 0
    freed_bytes = 0
    errors = []

    if os.path.isdir(cache_dir):
        for entry in os.listdir(cache_dir):
            path = os.path.join(cache_dir, entry)
            try:
                if os.path.isdir(path) and not os.path.islink(path):
                    size, count = _dir_size_and_count(path)
                    shutil.rmtree(path)
                    freed_bytes += size
                    removed_files += count
                else:
                    freed_bytes += os.path.getsize(path)
                    os.remove(path)
                    removed_files += 1
            except OSError as exc:
                logger.warning("Could not remove cache entry %s: %s", path, exc)
                errors.append(entry)

    # Recreate an empty cache dir so subsequent requests don't have to
    # special-case a missing directory.
    os.makedirs(cache_dir, exist_ok=True)

    logger.info(
        "Cache cleared: %d file(s), %d bytes freed from %s (%d error(s))",
        removed_files, freed_bytes, cache_dir, len(errors),
    )
    return jsonify({
        "removed_files": removed_files,
        "freed_bytes": freed_bytes,
        "errors": errors,
    })


# ---------------------------------------------------------------------------
# Electron density map (EDM) endpoints
# ---------------------------------------------------------------------------

@bp.route("/api/edm/<pdbid>")
@_with_cache_dir
def edm(pdbid):
    """
    Serve an electron density map for pdbid, downloading it if needed.

    By default this serves the full CCP4 2Fo-Fc map from the same source
    used by /api/density-box (RCSB/EBI), for entries analysed against the
    original PDB. When the frontend is displaying a PDB-REDO analysis
    (?source=pdb_redo), density-box's EBI proxy does NOT apply — PDB-REDO
    re-refines the model, so its density map differs from the original
    entry's. In that case we serve PDB-REDO's own precomputed density map
    instead (a CCP4-style .map file from the map-maker service, NOT the
    same file/format as the default .ccp4 map below).
    """
    pdbid = pdbid.lower()

    if request.args.get("source") == "pdb_redo":
        mapfile = pdb_redo_utils.get_EDM(pdbid)
        if not mapfile or not os.path.isfile(mapfile):
            abort(404, description=f"No PDB-REDO electron density map available for {pdbid}")
        return send_file(
            mapfile,
            mimetype="application/octet-stream",
            as_attachment=False,
            download_name=f"{pdbid}_final.map",
            conditional=True,  # enables Range requests / 304s for large maps
        )

    mapfile, _sigma = eds_utils.get_edm(pdbid)
    if not mapfile or not os.path.isfile(mapfile):
        abort(404, description=f"No electron density map available for {pdbid}")
    return send_file(
        mapfile,
        mimetype="application/octet-stream",
        as_attachment=False,
        download_name=f"{pdbid}.ccp4",
        conditional=True,  # enables Range requests / 304s for large maps
    )


@bp.route("/api/edm-exists/<pdbid>")
@_with_cache_dir
def edm_exists(pdbid):
    """
    Lightweight existence check for the full 2Fo-Fc map served by /api/edm
    (default RCSB/EBI source only — see eds_utils.edm_exists). Used by the
    Results-export feature to report whether a density map is available
    for a structure without actually downloading it.

    Note: this does NOT cover PDB-REDO (?source=pdb_redo) maps — those are
    generated on demand by PDB-REDO's map-maker service, and a successful
    PDB-REDO analysis already implies its map is obtainable, so the
    frontend infers that case itself instead of calling this endpoint.
    """
    pdbid = pdbid.lower()
    use_cache = request.args.get("use_cache", "1") != "0"
    exists = eds_utils.edm_exists(pdbid, use_cache=use_cache)
    return jsonify({"pdbid": pdbid, "exists": bool(exists)})


def _parse_xyz(raw, name):
    try:
        parts = [float(v) for v in raw.split(",")]
    except (AttributeError, ValueError):
        parts = []
    if len(parts) != 3:
        abort(400, description=f"'{name}' must be 3 comma-separated numbers, e.g. 1.0,2.0,3.0")
    return parts


def _parse_atom_list(raw):
    """'x1,y1,z1;x2,y2,z2;...' -> [{"center": [x,y,z]}, ...]. Aborts 400
    on malformed input, same convention as _parse_xyz above."""
    if not raw:
        return []
    atoms = []
    for i, triplet in enumerate(raw.split(";")):
        triplet = triplet.strip()
        if not triplet:
            continue
        atoms.append({"center": _parse_xyz(triplet, f"atoms[{i}]")})
    return atoms


@bp.route("/api/density-mask/<pdbid>/<region>")
@_with_cache_dir
def density_mask_region(pdbid, region):
    """
    Serve a pre-masked, per-region CCP4 density map:
    the source map is cropped to the region's padded box and every voxel
    further than `radius` from every atom in the region is zeroed out,
    server-side, so the 3D viewer only has to render a single isosurface
    for this region.

    Query params:
      min=minx,miny,minz    (required) region bounding box, e.g. from
                             a result's density_boxes[region]["min"]
      max=maxx,maxy,maxz    (required) same, density_boxes[region]["max"]
      atoms=x,y,z;x,y,z;... (required) atom centers, from a result's
                             density_atoms[region]
      radius=1.6             (optional) atom-mask radius in Å; quantized
                             and cached server-side, see
                             core.density_mask.quantize_radius
      source=pdb|pdb_redo    (optional, default "pdb")
    """
    pdbid = pdbid.lower()
    if region not in density_mask.VALID_REGIONS:
        return jsonify({"error": f"Unknown density region '{region}'"}), 404

    source = request.args.get("source", "pdb")
    if source not in ("pdb", "pdb_redo"):
        return jsonify({"error": f"Unknown source '{source}'"}), 400

    box = {
        "min": _parse_xyz(request.args.get("min", ""), "min"),
        "max": _parse_xyz(request.args.get("max", ""), "max"),
    }
    atoms = _parse_atom_list(request.args.get("atoms", ""))
    if not atoms:
        return jsonify({"error": "'atoms' must be a non-empty ';'-separated list of x,y,z triples"}), 400

    radius = request.args.get("radius", 1.6)

    mapfile = density_mask.get_masked_region_map(
        pdbid, region, box, atoms, radius=radius, source=source, use_cache=True,
    )
    if not mapfile or not os.path.isfile(mapfile):
        return jsonify({"error": f"No density available for {pdbid}/{region} ({source})"}), 404

    return send_file(
        mapfile,
        mimetype="application/octet-stream",
        as_attachment=False,
        download_name=f"{pdbid}_{region}.ccp4",
        conditional=True,
    )


@bp.route("/api/density-box/<pdbid>")
@_with_cache_dir
def density_box(pdbid):
    """
    Proxy a region query to EBI's density server and stream back the
    BinaryCIF chunk. Query params:
      min=minx,miny,minz   (required)
      max=maxx,maxy,maxz   (required)
      detail=0-6           (optional, default 3)
    """
    pdbid = pdbid.lower()
    box_min = _parse_xyz(request.args.get("min", ""), "min")
    box_max = _parse_xyz(request.args.get("max", ""), "max")
    detail = request.args.get("detail", "3")

    cachedir = os.path.join(pdb_utils.CACHEDIR, pdbid, "density-boxes")
    os.makedirs(cachedir, exist_ok=True)
    cachekey = "{}_{}_{}.bcif".format(
        "-".join(f"{v:.2f}" for v in box_min),
        "-".join(f"{v:.2f}" for v in box_max),
        detail,
    )
    cachepath = os.path.join(cachedir, cachekey)

    if not os.path.isfile(cachepath):
        url = eds_utils.edm_box_url(
            pdbid, {"min": box_min, "max": box_max}, detail=detail
        )
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
        except Exception as exc:
            logger.error("density-box fetch failed for %s: %s", pdbid, exc)
            abort(502, description="Could not fetch density data from EBI")
        with open(cachepath, "wb") as fh:
            fh.write(r.content)

    with open(cachepath, "rb") as fh:
        data = fh.read()
    return Response(data, mimetype="application/octet-stream")
