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
"""
import os
import re
import threading
import uuid
import logging
from functools import wraps

import requests
from flask import (
    Blueprint, render_template, request, jsonify, current_app,
    Response, abort, send_file,
)

import core.pdb_utils as pdb_utils
import core.eds_utils as eds_utils
import core.pdb_redo_utils as pdb_redo_utils
from core.rsr_core import AnalysisConfig, analyse_pdbids

logger = logging.getLogger(__name__)
bp = Blueprint("main", __name__)

# In-memory job store  {job_id: {"status": ..., "progress": ..., "results": ...}}
_jobs = {}
_jobs_lock = threading.Lock()

# Canonical UniProt accession pattern (6-char, or 10-char introduced in
# later releases), e.g. "P00734", "A0A0A0MRZ7". Anything that doesn't match
# this is treated as a plain PDB ID.
UNIPROT_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$",
    re.IGNORECASE,
)


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
    Expand any UniProt accessions found in *tokens* into the PDB IDs of the
    structures that reference them (core.pdb_utils.get_pdbids_for_uniprot).
    Plain PDB IDs pass through unchanged. De-duplicates while preserving
    order of first appearance.

    Returns (pdbids, origin_map, unresolved):
      pdbids      – flat list of PDB IDs to analyse
      origin_map  – {pdbid.lower(): uniprot_accession_or_None}, so results
                    can be tagged with the UniProt code they came from
      unresolved  – UniProt-looking tokens that returned no PDB entries
    """
    pdbids = []
    origin_map = {}
    unresolved = []
    seen = set()

    for token in tokens:
        if UNIPROT_RE.match(token):
            uniprot_id = token.upper()
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


def _run_job(job_id, pdbids, cfg, origin_map=None):
    # The actual concurrent analysis (thread pool, per-entry error
    # handling, order preservation) lives in core.rsr_core.analyse_pdbids
    # — see that function's docstring for why a thread pool is safe here.
    # This wrapper only adds the two things that are specific to *this*
    # web job: live progress updates for /api/status polling, and tagging
    # each result with the UniProt accession (if any) that produced it.
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
        # Tag with the UniProt accession that produced this PDB ID, if any,
        # so the frontend can group/label results accordingly.
        res["uniprot"] = origin_map.get(pdbid.strip().lower())

    with _jobs_lock:
        _jobs[job_id]["status"] = "done"
        _jobs[job_id]["results"] = results


@bp.route("/")
def index():
    return render_template("index.html")


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
        return jsonify({"error": "No PDB IDs or UniProt accessions provided"}), 400

    pdbids, origin_map, unresolved = _expand_ids(tokens)

    if not pdbids:
        msg = "No valid PDB IDs found."
        if unresolved:
            msg += " Could not find any PDB entries for UniProt accession(s): " + ", ".join(unresolved)
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
    )

    job_id = str(uuid.uuid4())
    with _jobs_lock:
        _jobs[job_id] = {"status": "running", "progress": 0, "total": len(pdbids), "results": None}

    t = threading.Thread(target=_run_job, args=(job_id, pdbids, cfg, origin_map), daemon=True)
    t.start()

    response = {"job_id": job_id, "total": len(pdbids)}
    if unresolved:
        response["warnings"] = [
            f"No PDB entries found for UniProt accession {u}" for u in unresolved
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


def _parse_xyz(raw, name):
    try:
        parts = [float(v) for v in raw.split(",")]
    except (AttributeError, ValueError):
        parts = []
    if len(parts) != 3:
        abort(400, description=f"'{name}' must be 3 comma-separated numbers, e.g. 1.0,2.0,3.0")
    return parts


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
