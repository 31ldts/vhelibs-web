# -*- coding: utf-8 -*-
#
#   Copyright 2010-2024 Adrià Cereto Massagué
#   Migrated to web version.
#   Changes: uses requests with proper SSL/timeout handling, no Jython/Cython deps,
#            CACHEDIR configurable via set_cache_dir().
#   Refactor: download/cache-on-disk logic now delegates to core.http_cache
#   (was duplicated near-verbatim across this module, eds_utils.py and
#   pdb_redo_utils.py); dropped the unused ``_g`` nested helper in
#   get_custom_report and split its field-extraction into
#   ``_extract_report_fields`` for readability.
#
import os
import logging
import tempfile

import requests

import core.http_cache as http_cache

logger = logging.getLogger(__name__)

PDBbase = "https://files.rcsb.org/download/{}.cif.gz"
PDBREDObase_full = "https://pdb-redo.eu/db/{pdbid}/{pdbid}_final.cif"
QUERY_TPL = "https://data.rcsb.org/rest/v1/core/entry/{}"

# RCSB Search API — used to resolve a UniProt ID to every PDB entry
# whose polymer entities reference it (see get_pdbids_for_uniprot below).
UNIPROT_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

CACHEDIR = os.path.join(tempfile.gettempdir(), "vhelibs_cache")


def set_cache_dir(path):
    """Set the directory used to cache downloaded files.

    Updates the module-level ``CACHEDIR`` global to point at ``path`` and
    ensures the directory exists, creating it (and any missing parent
    directories) if necessary.

    Args:
        path (str): Filesystem path to use as the new cache directory.

    Returns:
        None: This function does not return a value; it mutates the
        module-level ``CACHEDIR`` global in place.

    Raises:
        OSError: If the directory cannot be created (e.g. due to
            insufficient permissions).
    """
    global CACHEDIR
    CACHEDIR = path
    os.makedirs(CACHEDIR, exist_ok=True)


def _get_nested(d, *keys, default=None):
    """Safely retrieve a nested value from a dict of dicts.

    Walks through ``d`` following ``keys`` in order, returning ``default``
    as soon as any intermediate value is not a dict or is ``None``.

    Args:
        d (dict): Dictionary (possibly nested) to traverse.
        *keys: Sequence of keys to look up successively.
        default: Value to return if any step of the traversal is missing
            or ``None``. Defaults to ``None``.

    Returns:
        Any: The value found at the end of the key path, or ``default`` if
        the path could not be fully resolved.

    Raises:
        None
    """
    for key in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(key)
        if d is None:
            return default
    return d


def _extract_report_fields(rawdict):
    """Extract refinement/unit-cell statistics from a raw RCSB entry payload.

    Args:
        rawdict (dict): Parsed JSON response from the RCSB entry REST API
            (see ``QUERY_TPL``).

    Returns:
        dict: Dictionary of refinement/unit-cell fields with safe
        fallbacks (see :func:`get_custom_report` for the full key list).

    Raises:
        Exception: Propagates any error encountered while indexing into
            ``rawdict`` (e.g. if it isn't shaped like an RCSB entry
            payload at all); callers are expected to catch this.
    """
    info = rawdict.get("rcsb_entry_info") or {}
    refine = (rawdict.get("refine") or [{}])[0]  # first refinement block or empty dict
    cell = rawdict.get("cell") or {}

    method = info.get("experimental_method", "")

    # Resolution: try ls_d_res_high first, then ls_d_res_low (some entries only have one),
    # then the entry-level resolution from rcsb_entry_info.
    resolution = (
        refine.get("ls_d_res_high")
        or refine.get("ls_d_res_low")
        or info.get("resolution_combined", [None])[0]
        or 0.0
    )

    return {
        "experimentalTechnique": method,
        "rFree":                 refine.get("ls_R_factor_R_free") or 9999,
        "rWork":                 refine.get("ls_R_factor_R_work") or 9999,
        "refinementResolution":  float(resolution) if resolution else 0.0,
        "nreflections":          refine.get("ls_number_reflns_rfree") or 0,
        "unitCellAngleAlpha":    cell.get("angle_alpha") or 0.0,
        "unitCellAngleBeta":     cell.get("angle_beta") or 0.0,
        "unitCellAngleGamma":    cell.get("angle_gamma") or 0.0,
        "lengthOfUnitCellLatticeA": cell.get("length_a") or 0.0,
        "lengthOfUnitCellLatticeB": cell.get("length_b") or 0.0,
        "lengthOfUnitCellLatticeC": cell.get("length_c") or 0.0,
    }


def get_custom_report(pdbid, use_cache=True):
    """Fetch structure refinement/metadata statistics from the RCSB REST API.

    Retrieves (using a local cache when available and ``use_cache`` is
    ``True``) the RCSB entry data for the given PDB entry and extracts a
    set of refinement and unit-cell statistics from it. Every extracted
    value has a safe fallback, so a missing or null field never causes the
    whole entry to fail.

    Args:
        pdbid (str): PDB identifier of the structure to fetch. Case is
            normalized internally (upper-cased for the returned dict key,
            lower-cased for the API request URL).
        use_cache (bool, optional): Whether a cached response for this
            entry may be reused instead of re-querying the RCSB API.
            Defaults to ``True``.

    Returns:
        dict: A dictionary of the form ``{PDBID: rowdict}`` where
        ``rowdict`` contains the keys ``"experimentalTechnique"``,
        ``"rFree"``, ``"rWork"``, ``"refinementResolution"``,
        ``"nreflections"``, ``"unitCellAngleAlpha"``,
        ``"unitCellAngleBeta"``, ``"unitCellAngleGamma"``,
        ``"lengthOfUnitCellLatticeA"``, ``"lengthOfUnitCellLatticeB"``,
        and ``"lengthOfUnitCellLatticeC"``. Returns an empty dict ``{}``
        if the entry is unusable (e.g. NMR/cryo-EM with no refinement
        block) or if the data could not be fetched or parsed.

    Raises:
        None: Fetch, cache, and parsing errors are caught internally and
            reported via the module logger; the function returns ``{}``
            instead of propagating exceptions.
    """
    pdbid = pdbid.upper()
    url = QUERY_TPL.format(pdbid.lower())
    cachedir = http_cache.entry_cache_dir(CACHEDIR, pdbid)
    cache_path = os.path.join(cachedir, "pdb_stats.json")

    rawdict = http_cache.fetch_json(url, cache_path, timeout=30, use_cache=use_cache)
    if rawdict is None:
        return {}

    try:
        rowdict = _extract_report_fields(rawdict)
        # If there is genuinely no refinement data (pure NMR, etc.) the rFree
        # will be 9999 — the caller (rsr_core) will reject the entry gracefully.
        logger.debug("Stats for %s: rFree=%.4f res=%.2f Å method=%s",
                     pdbid, rowdict["rFree"], rowdict["refinementResolution"],
                     rowdict["experimentalTechnique"])
        return {pdbid: rowdict}
    except Exception as exc:
        logger.error("Error parsing stats for %s: %s", pdbid, exc)
        return {}


def _build_uniprot_query(uniprot_id, max_results):
    """Build the RCSB Search API request body for a UniProt->PDB lookup.

    Args:
        uniprot_id (str): UniProt ID to search for (exact match).
        max_results (int): Maximum number of PDB entries to request.

    Returns:
        dict: JSON-serializable request body for ``UNIPROT_SEARCH_URL``.

    Raises:
        None
    """
    return {
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


def get_pdbids_for_uniprot(uniprot_id, max_results=200):
    """Resolve a UniProt ID to the PDB entries that reference it.

    Queries the RCSB Search API to find every PDB entry whose polymer
    entities reference the given UniProt ID (e.g. ``"P00734"``),
    caching results on disk for subsequent calls.

    Args:
        uniprot_id (str): UniProt ID to resolve. Leading/trailing
            whitespace is stripped and the value is upper-cased before
            use.
        max_results (int, optional): Maximum number of PDB entries to
            request from the search API. Defaults to ``200``.

    Returns:
        list: A sorted list of upper-case 4-character PDB IDs referencing
        the given UniProt ID. Empty if there are no hits or if a
        network/parse failure occurs; callers should treat an empty list
        as "nothing found" rather than a hard error.

    Raises:
        None: Request and parsing errors are caught internally and
            reported via the module logger; the function returns ``[]``
            instead of propagating exceptions.
    """
    uniprot_id = uniprot_id.strip().upper()
    cachedir = os.path.join(CACHEDIR, "uniprot")
    os.makedirs(cachedir, exist_ok=True)
    cache_path = os.path.join(cachedir, f"{uniprot_id}.json")

    cached = http_cache.load_cached_json(cache_path)
    if cached is not None:
        return cached

    query = _build_uniprot_query(uniprot_id, max_results)
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

    http_cache.save_json(cache_path, pdbids)
    return pdbids


def get_pdb_file(pdbcode, pdb_redo=False, use_cache=True):
    """Download (and cache) the mmCIF structure file for a PDB entry.

    Fetches the mmCIF file either from the standard RCSB file archive or,
    if requested, from PDB-REDO's re-refined coordinate set. If a valid
    cached copy already exists locally and ``use_cache`` is ``True``, the
    download is skipped.

    Args:
        pdbcode (str): PDB identifier of the structure to fetch.
        pdb_redo (bool, optional): If ``True``, download the PDB-REDO
            re-refined ``.cif`` file instead of the standard RCSB
            ``.cif.gz`` file. Defaults to ``False``.
        use_cache (bool, optional): Whether an existing cached copy of the
            file may be reused instead of re-downloading it. Defaults to
            ``True``.

    Returns:
        str: Local filesystem path to the downloaded (or cached) mmCIF
        file, or an empty string ``""`` if the download failed.

    Raises:
        None: Download errors are caught internally (via
            :func:`core.http_cache.download_if_missing`) and reported via
            the module logger; the function returns ``""`` instead of
            propagating exceptions.
    """
    pdbcode_lower = pdbcode.lower()
    os.makedirs(CACHEDIR, exist_ok=True)
    if not pdb_redo:
        url = PDBbase.format(pdbcode_lower)
        filename = os.path.join(CACHEDIR, pdbcode.upper() + ".cif.gz")
    else:
        url = PDBREDObase_full.format(pdbid=pdbcode_lower)
        filename = os.path.join(CACHEDIR, os.path.basename(url))

    if http_cache.download_if_missing(url, filename, timeout=30, retries=3, use_cache=use_cache):
        logger.debug("Using file: %s", filename)
        return filename

    logger.error("Could not download %s", url)
    return ""
