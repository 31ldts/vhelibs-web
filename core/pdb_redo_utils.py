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
    """Download a URL to a local file, retrying on failure.

    Attempts to fetch ``url`` and write its raw content to ``dest_path``,
    retrying up to ``retries`` times if a request fails.

    Args:
        url (str): URL to download.
        dest_path (str): Local filesystem path where the downloaded
            content will be written.
        retries (int, optional): Maximum number of attempts before giving
            up. Defaults to ``3``.

    Returns:
        bool: ``True`` if the download succeeded and the file was written,
        ``False`` if all attempts failed.

    Raises:
        None: Request and I/O errors are caught internally on each
            attempt and logged as warnings; the function returns
            ``False`` instead of propagating exceptions.
    """
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
    """Fetch per-residue RSR/RSCC electron-density statistics from PDB-REDO.

    Downloads (and caches on disk) the PDB-REDO ``*_final.json`` file for
    the given entry, then parses it into a per-residue dictionary of
    density-fit statistics.

    Args:
        pdbid (str): PDB identifier of the structure to fetch. It is
            lower-cased before use.

    Returns:
        dict or None: A dictionary mapping residue key (as produced by
        :func:`core.pdb_atom.format_reskey`) to a dict with keys
        ``"RSR"`` and ``"RSCC"`` (floats). Returns ``None`` if the file
        could not be downloaded or parsed.

    Raises:
        None: Download and parsing errors are caught internally and
            reported via the module logger; the function returns ``None``
            instead of propagating exceptions.
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
    """Fetch overall structure refinement statistics from PDB-REDO.

    Retrieves (using a local cache when available) the PDB-REDO
    ``data.json`` file for the given entry and extracts a set of
    structure-level refinement properties from it.

    Args:
        pdbid (str): PDB identifier of the structure to fetch. It is
            lower-cased before use.

    Returns:
        dict or None: A dictionary compatible with
        ``pdb_utils.get_custom_report()``, containing the keys
        ``"experimentalTechnique"``, ``"rFree"``, ``"rWork"``,
        ``"refinementResolution"``, ``"unitCellAngleAlpha"``,
        ``"unitCellAngleBeta"``, ``"unitCellAngleGamma"``,
        ``"lengthOfUnitCellLatticeA"``, ``"lengthOfUnitCellLatticeB"``,
        ``"lengthOfUnitCellLatticeC"``, and ``"nreflections"``. Returns
        ``None`` if the data could not be fetched or parsed.

    Raises:
        None: Fetch, cache, and parsing errors are caught internally and
            reported via the module logger; the function returns ``None``
            instead of propagating exceptions.
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
    """Download (and cache) the PDB-REDO MTZ electron density map file.

    If a valid cached copy of the MTZ file already exists locally, the
    download is skipped and the cached path is returned directly.

    Args:
        pdbid (str): PDB identifier of the structure to fetch. It is
            lower-cased before use.

    Returns:
        str or None: Local filesystem path to the downloaded (or cached)
        MTZ file, or ``None`` if the download failed.

    Raises:
        None: Download errors are caught internally (via
            :func:`_download`) and reported via the module logger; the
            function returns ``None`` instead of propagating exceptions.
    """
    pdbid = pdbid.lower()
    downloaddir = os.path.join(pdb_utils.CACHEDIR, pdbid)
    os.makedirs(downloaddir, exist_ok=True)
    url = PDB_REDO_EDM_URL.format(pdbid=pdbid)
    filename = os.path.join(downloaddir, f"{pdbid}_final.mtz")
    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        return filename
    logger.info("Downloading %s", url)
    return filename if _download(url, filename) else None
