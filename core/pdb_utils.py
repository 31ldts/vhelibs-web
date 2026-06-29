# -*- coding: utf-8 -*-
#
#   Copyright 2010-2024 Adrià Cereto Massagué
#   Migrated to web version.
#   Changes: uses requests with proper SSL/timeout handling, no Jython/Cython deps,
#            CACHEDIR configurable via set_cache_dir().
#
import os
import json
import logging
import tempfile

import requests

logger = logging.getLogger(__name__)

PDBbase = "https://files.rcsb.org/download/{}.cif.gz"
PDBREDObase_full = "https://pdb-redo.eu/db/{pdbid}/{pdbid}_final.cif"
QUERY_TPL = "https://data.rcsb.org/rest/v1/core/entry/{}"

CACHEDIR = os.path.join(tempfile.gettempdir(), "vhelibs_cache")


def set_cache_dir(path):
    global CACHEDIR
    CACHEDIR = path
    os.makedirs(CACHEDIR, exist_ok=True)


def _download(url, dest_path, retries=3):
    """Download *url* to *dest_path* with retry logic. Returns True on success."""
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=30, verify=True)
            r.raise_for_status()
            with open(dest_path, "wb") as fh:
                fh.write(r.content)
            return True
        except Exception as exc:
            logger.warning("Download attempt %d/%d failed for %s: %s", attempt, retries, url, exc)
    return False


def get_custom_report(pdbid):
    """
    Fetch structure metadata from RCSB REST API.
    Returns {PDBID: rowdict} on success, or {} if the entry is unusable
    (e.g. NMR / cryo-EM with no refinement block, or network failure).

    Robust against missing/null fields: every value has a safe fallback so a
    KeyError on one field never kills the whole entry.
    """
    pdbid = pdbid.upper()
    url = QUERY_TPL.format(pdbid.lower())
    cachedir = os.path.join(CACHEDIR, pdbid.lower())
    os.makedirs(cachedir, exist_ok=True)
    cache_path = os.path.join(cachedir, "pdb_stats.json")

    # ── Load from cache or network ─────────────────────────────────────────
    rawdict = None
    if os.path.isfile(cache_path) and os.path.getsize(cache_path) > 0:
        logger.debug("Loading cached stats: %s", cache_path)
        try:
            with open(cache_path, "rt") as fh:
                rawdict = json.load(fh)
        except Exception:
            rawdict = None

    if rawdict is None:
        logger.info("Fetching %s", url)
        try:
            r = requests.get(url, timeout=30, verify=True)
            r.raise_for_status()
            rawdict = r.json()
            with open(cache_path, "wt") as fh:
                json.dump(rawdict, fh)
        except Exception as exc:
            logger.error("Could not fetch report for %s: %s", pdbid, exc)
            return {}

    # ── Extract fields with safe fallbacks ────────────────────────────────
    def _g(d, *keys, default=None):
        """Safely get a nested key from a dict; returns *default* if any step is missing/None."""
        for k in keys:
            if not isinstance(d, dict):
                return default
            d = d.get(k)
            if d is None:
                return default
        return d

    try:
        info    = rawdict.get("rcsb_entry_info") or {}
        refine  = (rawdict.get("refine") or [{}])[0]  # first refinement block or empty dict
        cell    = rawdict.get("cell") or {}

        method = info.get("experimental_method", "")

        # Resolution: try ls_dres_high first, then ls_dres_low (some entries only have one),
        # then the entry-level resolution from rcsb_entry_info.
        resolution = (
            refine.get("ls_dres_high")
            or refine.get("ls_dres_low")
            or info.get("resolution_combined", [None])[0]
            or 0.0
        )

        rowdict = {
            "experimentalTechnique": method,
            "rFree":                 refine.get("ls_rfactor_rfree")  or 9999,
            "rWork":                 refine.get("ls_rfactor_rwork")  or 9999,
            "refinementResolution":  float(resolution) if resolution else 0.0,
            "nreflections":          refine.get("ls_number_reflns_rfree") or 0,
            "unitCellAngleAlpha":    cell.get("angle_alpha")  or 0.0,
            "unitCellAngleBeta":     cell.get("angle_beta")   or 0.0,
            "unitCellAngleGamma":    cell.get("angle_gamma")  or 0.0,
            "lengthOfUnitCellLatticeA": cell.get("length_a") or 0.0,
            "lengthOfUnitCellLatticeB": cell.get("length_b") or 0.0,
            "lengthOfUnitCellLatticeC": cell.get("length_c") or 0.0,
        }

        # If there is genuinely no refinement data (pure NMR, etc.) the rFree
        # will be 9999 — the caller (rsr_core) will reject the entry gracefully.
        logger.debug("Stats for %s: rFree=%.4f res=%.2f Å method=%s",
                     pdbid, rowdict["rFree"], rowdict["refinementResolution"], method)
        return {pdbid: rowdict}

    except Exception as exc:
        logger.error("Error parsing stats for %s: %s", pdbid, exc)
        return {}


def get_pdb_file(pdbcode, pdb_redo=False):
    """Download the mmCIF file for *pdbcode*. Returns local file path or empty string."""
    pdbcode_lower = pdbcode.lower()
    os.makedirs(CACHEDIR, exist_ok=True)
    if not pdb_redo:
        url = PDBbase.format(pdbcode_lower)
        filename = os.path.join(CACHEDIR, pdbcode.upper() + ".cif.gz")
    else:
        url = PDBREDObase_full.format(pdbid=pdbcode_lower)
        filename = os.path.join(CACHEDIR, os.path.basename(url))

    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        logger.debug("Using cached file: %s", filename)
        return filename

    logger.info("Downloading %s → %s", url, filename)
    if _download(url, filename):
        return filename

    logger.error("Could not download %s", url)
    return ""
