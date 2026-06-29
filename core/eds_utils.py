# -*- coding: utf-8 -*-
#
#   Copyright 2010-2024 Adrià Cereto Massagué
#   Migrated to web version.
#   Changes: replaced urllib with requests; proper SSL/timeout; no Cython deps.
#
import os
import logging
import xml.etree.ElementTree as ET

import requests

import core.pdb_utils as pdb_utils

logger = logging.getLogger(__name__)

EDSTATS_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/{}_validation.xml"


def get_EDS(pdbid):
    """
    Fetch EDS validation statistics for *pdbid* from PDBe.

    Returns (pdbdict, edd_dict) where:
      pdbdict = {pdbid: True | False | reason_string}
      edd_dict = {residue_key: {"RSR": float, "RSCC": float, "OWAB": float, ...}}
    """
    pdbid = pdbid.lower()
    pdbdict = {pdbid: None}
    edd_dict = {}

    url = EDSTATS_URL.format(pdbid)
    downloaddir = os.path.join(pdb_utils.CACHEDIR, pdbid)
    os.makedirs(downloaddir, exist_ok=True)
    stat_path = os.path.join(downloaddir, f"{pdbid}_validation.xml")

    try:
        if not os.path.isfile(stat_path):
            logger.info("Downloading %s", url)
            r = requests.get(url, timeout=60, verify=True)
            if r.status_code == 404:
                pdbdict[pdbid] = False
                return pdbdict, edd_dict
            r.raise_for_status()
            with open(stat_path, "wb") as fh:
                fh.write(r.content)

        if not os.path.isfile(stat_path):
            pdbdict[pdbid] = False
            return pdbdict, edd_dict

        tree = ET.parse(stat_path)
        for res in tree.findall("ModelledSubgroup"):
            resname = list("   ")
            for i, c in enumerate(res.get("resname", "")[::-1]):
                resname[2 - i] = c
            resname = "".join(resname)

            resnum = list("    ")
            for i, c in enumerate(res.get("resnum", "")[::-1]):
                resnum[3 - i] = c
            resnum = "".join(resnum)

            residue = "{} {}{}{}".format(
                resname,
                res.get("chain", ""),
                resnum,
                res.get("icode", ""),
            ).strip()

            resdict = {
                "RSR": float(res.get("rsr") or 100),
                "RSCC": float(res.get("rscc") or 0),
                "OWAB": float(res.get("owab") or 1000),
                "RSRZ": float(res.get("rsrz") or 9999),
                "occupancy": float(res.get("avgoccu") or 0),
            }
            edd_dict[residue] = resdict

        pdbdict[pdbid] = True

    except Exception as exc:
        logger.error("EDS fetch error for %s: %s", pdbid, exc)
        pdbdict[pdbid] = str(exc)

    return pdbdict, edd_dict
