# -*- coding: utf-8 -*-
#
#   Copyright 2010-2024 Adrià Cereto Massagué
#   Migrated to web version.
#   Changes: replaced urllib with requests; proper SSL/timeout; no Cython deps.
#   Added: get_edm() to fetch the full electron density map (CCP4 format),
#   which EDS_parser.py used to fetch from the now-defunct
#   edmaps.rcsb.org/maps/{}_2fofc.dsn6 endpoint.
#   Refactor: get_edm/get_EDS now share their single-attempt,
#   404-aware download logic via core.http_cache.download_or_404 (this was
#   duplicated, with only minor variations, across this module,
#   pdb_utils.py and pdb_redo_utils.py); the XML parsing in get_EDS was
#   split into the private _parse_validation_xml helper.
#
import os
import logging
import xml.etree.ElementTree as ET

import core.pdb_utils as pdb_utils
import core.http_cache as http_cache
from core.pdb_atom import format_reskey

logger = logging.getLogger(__name__)

EDSTATS_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/{}_validation.xml"

# Full 2Fo-Fc map in CCP4 format. Confirmed working endpoint (see the
# official Mol*/MolViewSpec documentation examples, which use this same
# URL pattern against 1tqn). Mol* can parse .ccp4/.map files client-side
# via `parse({format: 'map'})` or `plugin.builders.volume.parseVolume`.
EDM_URL = "https://www.ebi.ac.uk/pdbe/entry-files/{}.ccp4"

# Per-region density "box" endpoint (EBI's density VolumeServer). Given a
# bounding box (in Angstrom, orthogonal coordinates) it streams back only
# the relevant chunk of the map as BinaryCIF, with two channels: 2FO-FC
# and FO-FC. This is what the 3D viewer uses for segmented density
# (ligand / binding site / residues-to-examine), since it avoids
# downloading and masking the whole map in the browser.
EDM_BOX_URL = (
    "https://www.ebi.ac.uk/pdbe/densities/x-ray/{pdbid}/box/"
    "{minx},{miny},{minz}/{maxx},{maxy},{maxz}?detail={detail}"
)


def edm_box_url(pdbid, box, detail=3):
    """Build the EBI density-server "box" query URL for a given bounding box.

    Constructs a URL against the EBI density VolumeServer's "box" endpoint,
    which streams back only the density data within the specified
    bounding box rather than the full electron density map.

    Args:
        pdbid (str): PDB identifier of the structure. It is lower-cased
            before being inserted into the URL.
        box (dict): Bounding box in Angstrom, orthogonal coordinates, with
            the shape ``{"min": [x, y, z], "max": [x, y, z]}``.
        detail (int, optional): Level of detail for the returned density,
            ranging from ``0`` (coarsest) to ``6`` (finest). Defaults to
            ``3``, a good default for ligand-sized regions.

    Returns:
        str: The fully formatted density-server "box" query URL.

    Raises:
        KeyError: If ``box`` does not contain both ``"min"`` and ``"max"``
            keys.
        ValueError: If ``box["min"]`` or ``box["max"]`` cannot be unpacked
            into exactly three coordinate values.
    """
    (minx, miny, minz) = box["min"]
    (maxx, maxy, maxz) = box["max"]
    return EDM_BOX_URL.format(
        pdbid=pdbid.lower(),
        minx=minx, miny=miny, minz=minz,
        maxx=maxx, maxy=maxy, maxz=maxz,
        detail=detail,
    )


def get_edm(pdbid, use_cache=True):
    """Download and cache the full 2Fo-Fc electron density map for a PDB entry.

    Fetches the CCP4-format electron density map from the EBI entry-files
    endpoint, storing it under the shared PDB cache directory. If a valid
    cached copy already exists and ``use_cache`` is ``True``, the download
    is skipped.

    Note:
        `sigma` is a fixed placeholder (``1.0``) here: unlike the old
        ``.dsn6`` format, CCP4 maps carry their own header statistics
        (mean/rms), so Mol* / any CCP4-aware client computes contour
        levels directly from the file instead of needing a
        separately-tracked sigma value. It is kept in the return
        signature for backwards compatibility with callers of the old
        ``EDS_parser.get_EDM()``.

    Args:
        pdbid (str): PDB identifier of the structure to fetch. It is
            lower-cased before use.
        use_cache (bool, optional): Whether to reuse an existing cached
            map file instead of re-downloading it. Defaults to ``True``.

    Returns:
        tuple: A 2-tuple ``(filepath, sigma)``:

            - filepath (str or None): Path to the downloaded (or cached)
              CCP4 map file, or ``None`` if the download failed or no map
              is available.
            - sigma (float or None): Fixed placeholder value of ``1.0``
              on success, or ``None`` on failure.

    Raises:
        None: All request and I/O errors are caught internally and
            reported via the module logger; the function returns
            ``(None, None)`` instead of propagating exceptions.
    """
    pdbid = pdbid.lower()
    downloaddir = http_cache.entry_cache_dir(pdb_utils.CACHEDIR, pdbid)
    mapfile = os.path.join(downloaddir, f"{pdbid}.ccp4")

    if use_cache and os.path.isfile(mapfile) and os.path.getsize(mapfile) > 0:
        return mapfile, 1.0

    url = EDM_URL.format(pdbid)
    logger.info("Downloading EDM %s", url)
    result = http_cache.download_or_404(url, mapfile, timeout=120, stream=True)

    if result is None:
        logger.warning("No EDM available for %s", pdbid)
        return None, None
    if result is False:
        return None, None
    return mapfile, 1.0


def _parse_validation_xml(xml_path):
    """Parse a PDBe validation XML file into a per-residue statistics dict.

    Args:
        xml_path (str): Path to a previously downloaded
            ``{pdbid}_validation.xml`` file.

    Returns:
        dict: Mapping of residue key (as produced by
        :func:`core.pdb_atom.format_reskey`) to a dict with keys
        ``"RSR"``, ``"RSCC"``, ``"OWAB"``, ``"RSRZ"``, and ``"occupancy"``
        (all floats).

    Raises:
        xml.etree.ElementTree.ParseError: If ``xml_path`` is not
            well-formed XML.
        OSError: If ``xml_path`` cannot be read.
    """
    edd_dict = {}
    tree = ET.parse(xml_path)
    for res in tree.findall("ModelledSubgroup"):
        residue = format_reskey(
            res.get("resname", ""),
            res.get("chain", ""),
            res.get("resnum", "") or "0",
            res.get("icode", ""),
        )
        edd_dict[residue] = {
            "RSR": float(res.get("rsr") or 100),
            "RSCC": float(res.get("rscc") or 0),
            "OWAB": float(res.get("owab") or 1000),
            "RSRZ": float(res.get("rsrz") or 9999),
            "occupancy": float(res.get("avgoccu") or 0),
        }
    return edd_dict


def get_EDS(pdbid):
    """Fetch EDS (Electron Density Server) validation statistics from PDBe.

    Downloads (and caches on disk) the PDBe validation XML file for the
    given PDB entry, then parses per-residue electron-density validation
    metrics (RSR, RSCC, OWAB, RSRZ, occupancy) from it.

    Args:
        pdbid (str): PDB identifier of the structure to fetch. It is
            lower-cased before use.

    Returns:
        tuple: A 2-tuple ``(pdbdict, edd_dict)``:

            - pdbdict (dict): Mapping ``{pdbid: status}`` where ``status``
              is ``True`` on success, ``False`` if no validation data is
              available (e.g. HTTP 404 or missing file), or a string
              containing the error message if an exception occurred.
            - edd_dict (dict): Mapping of residue key (as produced by
              :func:`core.pdb_atom.format_reskey`) to a dict with keys
              ``"RSR"``, ``"RSCC"``, ``"OWAB"``, ``"RSRZ"``, and
              ``"occupancy"`` (all floats). Empty if no data could be
              parsed.

    Raises:
        None: All request, I/O, and XML-parsing errors are caught
            internally and reported via the module logger and the
            returned ``pdbdict`` status string; the function does not
            propagate exceptions.
    """
    pdbid = pdbid.lower()
    pdbdict = {pdbid: None}

    url = EDSTATS_URL.format(pdbid)
    downloaddir = http_cache.entry_cache_dir(pdb_utils.CACHEDIR, pdbid)
    stat_path = os.path.join(downloaddir, f"{pdbid}_validation.xml")

    try:
        if not os.path.isfile(stat_path):
            logger.info("Downloading %s", url)
            result = http_cache.download_or_404(url, stat_path, timeout=60)
            if not result:
                # None (404) and False (request error) both mean "no data".
                pdbdict[pdbid] = False
                return pdbdict, {}

        if not os.path.isfile(stat_path):
            pdbdict[pdbid] = False
            return pdbdict, {}

        edd_dict = _parse_validation_xml(stat_path)
        pdbdict[pdbid] = True
        return pdbdict, edd_dict

    except Exception as exc:
        logger.error("EDS fetch error for %s: %s", pdbid, exc)
        pdbdict[pdbid] = str(exc)
        return pdbdict, {}


# Clearer alias used by newer call sites; kept get_EDS for backwards
# compatibility with rsr_core.py and any other existing callers.
get_validation_data = get_EDS
