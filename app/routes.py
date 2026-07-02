# -*- coding: utf-8 -*-
"""
Flask routes for VHELIBS web.

Endpoints:
  GET  /                        – serve the SPA
  POST /api/analyse             – start an analysis job
  GET  /api/status/<job>        – poll job progress / results
  GET  /api/edm/<pdbid>         – full 2Fo-Fc map (CCP4), downloaded/cached on demand
  GET  /api/density-box/<pdbid> – cropped density chunk (BinaryCIF) for a given region,
                                   proxied + cached from EBI's density server
"""
import os
import threading
import uuid
import logging
import requests
from flask import (
    Blueprint, render_template, request, jsonify, current_app,
    Response, abort, send_file,
)

import core.pdb_utils as pdb_utils
import core.eds_utils as eds_utils
from core.rsr_core import AnalysisConfig, analyse_pdbids

logger = logging.getLogger(__name__)
bp = Blueprint("main", __name__)

# In-memory job store  {job_id: {"status": ..., "progress": ..., "results": ...}}
_jobs = {}
_jobs_lock = threading.Lock()


def _run_job(job_id, pdbids, cfg):
    total = len(pdbids)
    results = []
    for i, pdbid in enumerate(pdbids, 1):
        from core.rsr_core import parse_binding_site
        try:
            res = parse_binding_site(pdbid.strip().lower(), cfg)
        except Exception as exc:
            logger.exception("Error analysing %s", pdbid)
            res = {"pdbid": pdbid, "error": str(exc)}
        results.append(res)
        with _jobs_lock:
            _jobs[job_id]["progress"] = i
            _jobs[job_id]["total"] = total

    with _jobs_lock:
        _jobs[job_id]["status"] = "done"
        _jobs[job_id]["results"] = results


@bp.route("/")
def index():
    return render_template("index.html")


@bp.route("/api/analyse", methods=["POST"])
def analyse():
    data = request.get_json(force=True, silent=True) or {}

    raw_ids = data.get("pdbids", "")
    if isinstance(raw_ids, list):
        pdbids = [p.strip() for p in raw_ids if p.strip()]
    else:
        pdbids = [p.strip() for p in str(raw_ids).replace(",", "\n").split() if p.strip()]

    if not pdbids:
        return jsonify({"error": "No PDB IDs provided"}), 400

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
    )

    # Set cache dir from app config
    pdb_utils.set_cache_dir(current_app.config["CACHE_DIR"])

    job_id = str(uuid.uuid4())
    with _jobs_lock:
        _jobs[job_id] = {"status": "running", "progress": 0, "total": len(pdbids), "results": None}

    t = threading.Thread(target=_run_job, args=(job_id, pdbids, cfg), daemon=True)
    t.start()

    return jsonify({"job_id": job_id, "total": len(pdbids)})


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
def edm(pdbid):
    """Serve the full CCP4 2Fo-Fc map for pdbid, downloading it if needed."""
    pdb_utils.set_cache_dir(current_app.config["CACHE_DIR"])
    pdbid = pdbid.lower()
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
def density_box(pdbid):
    """
    Proxy a region query to EBI's density server and stream back the
    BinaryCIF chunk. Query params:
      min=minx,miny,minz   (required)
      max=maxx,maxy,maxz   (required)
      detail=0-6           (optional, default 3)
    """
    pdb_utils.set_cache_dir(current_app.config["CACHE_DIR"])
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
