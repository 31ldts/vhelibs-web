# -*- coding: utf-8 -*-
#
#   Copyright 2013-2024 Adrià Cereto Massagué
#   Migrated to web version.
#   Changes: replaced urllib with requests; removed Java/Jython locale hacks;
#            removed CSV-based legacy fallback; no Cython deps.
#   Refactor: download/cache-on-disk logic now delegates to core.http_cache
#   (was duplicated near-verbatim across this module, pdb_utils.py and
#   eds_utils.py); split JSON parsing in get_ED_data and get_pdbredo_data
#   into private helpers.
#
import os
import logging

import core.pdb_utils as pdb_utils
import core.http_cache as http_cache
from core.pdb_atom import format_reskey

logger = logging.getLogger(__name__)

PDB_REDO_ED_DATA_URL = "https://pdb-redo.eu/db/{pdbid}/{pdbid}_final.json"
PDB_REDO_EDM_URL = "https://pdb-redo.eu/db/{pdbid}/{pdbid}_final.mtz"
ALLDATA_URL = "https://pdb-redo.eu/db/{pdbid}/data.json"

# Timeouts here are intentionally longer than pdb_utils'/eds_utils' defaults:
# PDB-REDO's server has historically been slower to respond than RCSB/EBI.
_DOWNLOAD_TIMEOUT = 60


def _parse_ed_data(ed_data):
    """Convert PDB-REDO's raw ``*_final.json`` payload into a per-residue dict.

    Args:
        ed_data (list): Parsed JSON payload as returned by the PDB-REDO
            ``{pdbid}_final.json`` endpoint: a list of per-component
            dictionaries.

    Returns:
        dict: Mapping of residue key (as produced by
        :func:`core.pdb_atom.format_reskey`) to a dict with keys
        ``"RSR"`` and ``"RSCC"`` (floats).

    Raises:
        Exception: Propagates any error from indexing into a malformed
            entry; callers are expected to catch this.
    """
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
    downloaddir = http_cache.entry_cache_dir(pdb_utils.CACHEDIR, pdbid)
    url = PDB_REDO_ED_DATA_URL.format(pdbid=pdbid)
    filename = os.path.join(downloaddir, f"{pdbid}_final.json")

    if not http_cache.download_if_missing(url, filename, timeout=_DOWNLOAD_TIMEOUT):
        logger.error("Unable to download %s", url)
        return None

    ed_data = http_cache.load_cached_json(filename)
    if ed_data is None:
        logger.error("Could not parse %s", filename)
        return None

    try:
        return _parse_ed_data(ed_data)
    except Exception as exc:
        logger.error("Could not parse %s: %s", filename, exc)
        return None


def _extract_pdbredo_props(rawdict):
    """Extract structure-level refinement properties from PDB-REDO's ``data.json``.

    Args:
        rawdict (dict): Parsed JSON payload from the PDB-REDO ``data.json``
            endpoint (see ``ALLDATA_URL``).

    Returns:
        dict: Dictionary compatible with ``pdb_utils.get_custom_report()``
        (see :func:`get_pdbredo_data` for the full key list).

    Raises:
        KeyError: If ``rawdict`` has no ``"properties"`` key, i.e. it isn't
            shaped like a PDB-REDO ``data.json`` payload; callers are
            expected to catch this.
    """
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
    cachedir = http_cache.entry_cache_dir(pdb_utils.CACHEDIR, pdbid)
    url = ALLDATA_URL.format(pdbid=pdbid)
    cache_path = os.path.join(cachedir, "data.json")

    rawdict = http_cache.fetch_json(url, cache_path, timeout=_DOWNLOAD_TIMEOUT)
    if rawdict is None:
        logger.error("Could not fetch PDB-REDO data for %s", pdbid)
        return None

    try:
        return _extract_pdbredo_props(rawdict)
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
            :func:`core.http_cache.download_if_missing`) and reported via
            the module logger; the function returns ``None`` instead of
            propagating exceptions.
    """
    pdbid = pdbid.lower()
    downloaddir = http_cache.entry_cache_dir(pdb_utils.CACHEDIR, pdbid)
    url = PDB_REDO_EDM_URL.format(pdbid=pdbid)
    filename = os.path.join(downloaddir, f"{pdbid}_final.mtz")

    if http_cache.download_if_missing(url, filename, timeout=_DOWNLOAD_TIMEOUT):
        return filename

    logger.error("Could not download %s", url)
    return None
