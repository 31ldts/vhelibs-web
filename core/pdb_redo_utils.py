# -*- coding: utf-8 -*-
#
#   Copyright 2013-2024 Adrià Cereto Massagué
#   Migrated to web version.
#   Changes: replaced urllib with requests; removed Java/Jython locale hacks;
#            removed CSV-based legacy fallback; no Cython deps.
#
import os
import json
import logging

import requests

import core.pdb_utils as pdb_utils
from core.pdb_atom import format_reskey

logger = logging.getLogger(__name__)

PDB_REDO_ED_DATA_URL = "https://pdb-redo.eu/db/{pdbid}/{pdbid}_final.json"
PDB_REDO_EDM_URL = "https://pdb-redo.eu/db/{pdbid}/{pdbid}_final.mtz"
ALLDATA_URL = "https://pdb-redo.eu/db/{pdbid}/data.json"


def _download(url, dest_path, retries=3):
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=60, verify=True)
            r.raise_for_status()
            with open(dest_path, "wb") as fh:
                fh.write(r.content)
            return True
        except Exception as exc:
            logger.warning("Attempt %d/%d failed for %s: %s", attempt, retries, url, exc)
    return False


def get_ED_data(pdbid):
    """
    Fetch per-residue RSR/RSCC from PDB-REDO for *pdbid*.
    Returns edd_dict {residue_key: {"RSR": float, "RSCC": float}} or None on failure.
    """
    pdbid = pdbid.lower()
    downloaddir = os.path.join(pdb_utils.CACHEDIR, pdbid)
    os.makedirs(downloaddir, exist_ok=True)
    url = PDB_REDO_ED_DATA_URL.format(pdbid=pdbid)
    filename = os.path.join(downloaddir, f"{pdbid}_final.json")

    if not (os.path.isfile(filename) and os.path.getsize(filename) > 0):
        logger.info("Downloading %s", url)
        if not _download(url, filename):
            logger.error("Unable to download %s", url)
            return None

    try:
        with open(filename, "rt") as fh:
            ed_data = json.load(fh)
    except Exception as exc:
        logger.error("Could not parse %s: %s", filename, exc)
        return None

    edd_dict = {}
    for comp in ed_data:
        residue = format_reskey(
            comp["pdb"]["compID"],
            comp["pdb"]["strandID"],
            comp["pdb"]["seqNum"],
        )
        edd_dict[residue] = {
            "RSR": float(comp.get("RSR") or 100),
            "RSCC": float(comp.get("RSCCS") or 0),
        }
    return edd_dict


def get_pdbredo_data(pdbid):
    """
    Fetch overall structure statistics from PDB-REDO for *pdbid*.
    Returns a rowdict compatible with pdb_utils.get_custom_report() or None.
    """
    pdbid = pdbid.lower()
    cachedir = os.path.join(pdb_utils.CACHEDIR, pdbid)
    os.makedirs(cachedir, exist_ok=True)
    url = ALLDATA_URL.format(pdbid=pdbid)
    cache_path = os.path.join(cachedir, "data.json")

    rawdict = None
    if os.path.isfile(cache_path) and os.path.getsize(cache_path) > 0:
        logger.debug("Loading cached PDB-REDO data: %s", cache_path)
        try:
            with open(cache_path, "rt") as fh:
                rawdict = json.load(fh)
        except Exception:
            rawdict = None

    if rawdict is None:
        logger.info("Fetching %s", url)
        try:
            r = requests.get(url, timeout=60, verify=True)
            r.raise_for_status()
            rawdict = r.json()
            with open(cache_path, "wt") as fh:
                json.dump(rawdict, fh)
        except Exception as exc:
            logger.error("Could not fetch PDB-REDO data for %s: %s", pdbid, exc)
            return None

    try:
        props = rawdict["properties"]
        return {
            "experimentalTechnique": props.get("EXPTYP"),
            "rFree": props.get("RFFIN", 9999),
            "rWork": props.get("RFIN", 9999),
            "refinementResolution": props.get("RESOLUTION", 0),
            "unitCellAngleAlpha": props.get("ALPHA", 0),
            "unitCellAngleBeta": props.get("BETA", 0),
            "unitCellAngleGamma": props.get("GAMMA", 0),
            "lengthOfUnitCellLatticeA": props.get("AAXIS", 0),
            "lengthOfUnitCellLatticeB": props.get("BAXIS", 0),
            "lengthOfUnitCellLatticeC": props.get("CAXIS", 0),
            "nreflections": props.get("NREFCNT", 0),
        }
    except Exception as exc:
        logger.error("Error parsing PDB-REDO properties for %s: %s", pdbid, exc)
        return None


def get_EDM(pdbid):
    """Download the PDB-REDO MTZ map file. Returns local path or None."""
    pdbid = pdbid.lower()
    downloaddir = os.path.join(pdb_utils.CACHEDIR, pdbid)
    os.makedirs(downloaddir, exist_ok=True)
    url = PDB_REDO_EDM_URL.format(pdbid=pdbid)
    filename = os.path.join(downloaddir, f"{pdbid}_final.mtz")
    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        return filename
    logger.info("Downloading %s", url)
    return filename if _download(url, filename) else None
