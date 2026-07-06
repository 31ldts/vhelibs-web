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

# RCSB Search API — used to resolve a UniProt accession to every PDB entry
# whose polymer entities reference it (see get_pdbids_for_uniprot below).
UNIPROT_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

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

        # Resolution: try ls_d_res_high first, then ls_d_res_low (some entries only have one),
        # then the entry-level resolution from rcsb_entry_info.
        resolution = (
            refine.get("ls_d_res_high")
            or refine.get("ls_d_res_low")
            or info.get("resolution_combined", [None])[0]
            or 0.0
        )

        rowdict = {
            "experimentalTechnique": method,
            "rFree":                 refine.get("ls_R_factor_R_free")  or 9999,
            "rWork":                 refine.get("ls_R_factor_R_work")  or 9999,
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


def get_pdbids_for_uniprot(uniprot_id, max_results=200):
    """
    Resolve a UniProt accession (e.g. "P00734") to the PDB entries whose
    polymer entities reference it, via the RCSB Search API.

    Returns a list of upper-case 4-character PDB IDs (empty on no hits or
    on network/parse failure — callers should treat that as "nothing
    found" rather than a hard error).
    """
    uniprot_id = uniprot_id.strip().upper()
    cachedir = os.path.join(CACHEDIR, "uniprot")
    os.makedirs(cachedir, exist_ok=True)
    cache_path = os.path.join(cachedir, f"{uniprot_id}.json")

    if os.path.isfile(cache_path) and os.path.getsize(cache_path) > 0:
        try:
            with open(cache_path, "rt") as fh:
                return json.load(fh)
        except Exception:
            pass  # fall through and re-fetch on a corrupt cache entry

    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": (
                    "rcsb_polymer_entity_container_identifiers."
                    "reference_sequence_identifiers.database_accession"
                ),
                "operator": "exact_match",
                "value": uniprot_id,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_results}},
    }

    try:
        r = requests.post(UNIPROT_SEARCH_URL, json=query, timeout=30)
        if r.status_code == 204:
            # RCSB's convention for "search executed fine, zero hits"
            pdbids = []
        else:
            r.raise_for_status()
            data = r.json()
            pdbids = sorted({
                hit["identifier"].upper()
                for hit in data.get("result_set", [])
                if hit.get("identifier")
            })
    except Exception as exc:
        logger.error("UniProt->PDB lookup failed for %s: %s", uniprot_id, exc)
        return []

    try:
        with open(cache_path, "wt") as fh:
            json.dump(pdbids, fh)
    except Exception:
        pass

    return pdbids


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
