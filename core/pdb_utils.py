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
            r = requests.get(url, timeout=30, verify=True)
            r.raise_for_status()
            with open(dest_path, "wb") as fh:
                fh.write(r.content)
            return True
        except Exception as exc:
            logger.warning("Download attempt %d/%d failed for %s: %s", attempt, retries, url, exc)
    return False


def get_custom_report(pdbid):
    """Fetch structure refinement/metadata statistics from the RCSB REST API.

    Retrieves (using a local cache when available) the RCSB entry data for
    the given PDB entry and extracts a set of refinement and unit-cell
    statistics from it. Every extracted value has a safe fallback, so a
    missing or null field never causes the whole entry to fail.

    Args:
        pdbid (str): PDB identifier of the structure to fetch. Case is
            normalized internally (upper-cased for the returned dict key,
            lower-cased for the API request URL).

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
        """Safely retrieve a nested value from a dict of dicts.

        Walks through ``d`` following ``keys`` in order, returning
        ``default`` as soon as any intermediate value is not a dict or is
        ``None``.

        Args:
            d (dict): Dictionary (possibly nested) to traverse.
            *keys: Sequence of keys to look up successively.
            default: Value to return if any step of the traversal is
                missing or ``None``. Defaults to ``None``.

        Returns:
            Any: The value found at the end of the key path, or
            ``default`` if the path could not be fully resolved.

        Raises:
            None
        """
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
    """Resolve a UniProt accession to the PDB entries that reference it.

    Queries the RCSB Search API to find every PDB entry whose polymer
    entities reference the given UniProt accession (e.g. ``"P00734"``),
    caching results on disk for subsequent calls.

    Args:
        uniprot_id (str): UniProt accession to resolve. Leading/trailing
            whitespace is stripped and the value is upper-cased before
            use.
        max_results (int, optional): Maximum number of PDB entries to
            request from the search API. Defaults to ``200``.

    Returns:
        list: A sorted list of upper-case 4-character PDB IDs referencing
        the given UniProt accession. Empty if there are no hits or if a
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
    """Download (and cache) the mmCIF structure file for a PDB entry.

    Fetches the mmCIF file either from the standard RCSB file archive or,
    if requested, from PDB-REDO's re-refined coordinate set. If a valid
    cached copy already exists locally, the download is skipped.

    Args:
        pdbcode (str): PDB identifier of the structure to fetch.
        pdb_redo (bool, optional): If ``True``, download the PDB-REDO
            re-refined ``.cif`` file instead of the standard RCSB
            ``.cif.gz`` file. Defaults to ``False``.

    Returns:
        str: Local filesystem path to the downloaded (or cached) mmCIF
        file, or an empty string ``""`` if the download failed.

    Raises:
        None: Download errors are caught internally (via
            :func:`_download`) and reported via the module logger; the
            function returns ``""`` instead of propagating exceptions.
    """
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
