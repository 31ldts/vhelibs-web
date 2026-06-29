# -*- coding: utf-8 -*-
#
#   Copyright 2011-2024 Adrià Cereto Massagué
#   Migrated to web version.
#
#   Major changes vs rsr_analysis.py:
#     - Removed CLI argument parsing (argparse).
#     - Removed multiprocessing / multithreading pool (not needed for Flask).
#     - Removed CSV/file output; returns dicts/lists suitable for JSON serialisation.
#     - Removed Java/Jython/Cython imports; uses pure-Python PdbAtom.
#     - Uses refactored core.* modules instead of top-level modules.
#     - Replaced print() with logging.
#     - Fixed SSL context hack removed (requests handles SSL properly).
#
import math
import logging

import gemmi

from core.pdb_atom import PdbAtom
import core.pdb_utils as pdb_utils
import core.eds_utils as eds_utils
import core.pdb_redo_utils as pdb_redo_utils
import core.cofactors as cofactors

logger = logging.getLogger(__name__)


class AnalysisConfig:
    """All tunable parameters for a single analysis run."""

    def __init__(
        self,
        rsr_upper=0.4,
        rsr_lower=0.24,
        rscc_min=0.9,
        rfree_max=1.0,
        occupancy_min=1.0,
        tolerance=2,
        distance=4.5,
        check_owab=False,
        owab_max=50.0,
        check_resolution=False,
        resolution_max=3.5,
        use_rdiff=False,
        rdiff_max=0.05,
        use_dpi=False,
        dpi_max=0.42,
        pdb_redo=False,
    ):
        self.rsr_upper = rsr_upper
        self.rsr_lower = rsr_lower
        self.rscc_min = rscc_min
        self.rfree_max = rfree_max
        self.occupancy_min = occupancy_min
        self.tolerance = tolerance
        self.inner_distance = distance ** 2
        self.check_owab = check_owab
        self.owab_max = owab_max
        self.check_resolution = check_resolution
        self.resolution_max = resolution_max
        self.use_rdiff = use_rdiff
        self.rdiff_max = rdiff_max
        self.use_dpi = use_dpi
        self.dpi_max = dpi_max
        self.pdb_redo = pdb_redo


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _average_occ(residue_atoms):
    return sum(a.occupancy for a in residue_atoms) / len(residue_atoms)


def _dpi(a, b, c, alpha, beta, gamma, natoms, reflections, rfree):
    cosa = math.cos(math.radians(alpha))
    cosb = math.cos(math.radians(beta))
    cosg = math.cos(math.radians(gamma))
    V = a * b * c * math.sqrt(max(0.0, 1 - cosa**2 - cosb**2 - cosg**2 + 2*cosa*cosb*cosg))
    if V <= 0 or reflections <= 0:
        return float("nan")
    return 1.28 * (natoms ** 0.5) * (V ** (1/3)) * (reflections ** (-5/6)) * rfree


# ---------------------------------------------------------------------------
# mmCIF parsing
# ---------------------------------------------------------------------------

def _fmt_reskey(comp_id, asym_id, seq_id):
    """Build a fixed-width residue key identical to the original code's format."""
    pos = str(seq_id)
    while len(pos) < 4:
        pos = " " + pos
    return "{} {}{}".format(comp_id, asym_id, pos)


def parse_mmcif_file(mmciffilepath, pdbid, inner_distance):
    """
    Parse an mmCIF file (plain or gzip) using gemmi and return the same
    five-tuple as the original pdbx-based implementation:
        (natoms, res_atom_dict, ligand_res_atom_dict, notligands, links)
    or a 1-tuple (error_string,) on failure.
    """
    natoms = 0
    res_atom_dict = {}
    ligand_res_atom_dict = {}
    notligands = {}
    links = []

    logger.debug("Reading %s", mmciffilepath)
    try:
        structure = gemmi.read_structure(mmciffilepath)
    except Exception as exc:
        return (f"Could not parse mmCIF file {mmciffilepath}: {exc}",)

    # ── Atom loop ──────────────────────────────────────────────────────────
    for model in structure:
        for chain in model:
            for residue in chain:
                comp_id  = residue.name          # auth_comp_id  (e.g. "ATP")
                asym_id  = chain.name            # auth_asym_id  (e.g. "A")
                seq_id   = residue.seqid.num     # auth_seq_id   (integer)
                res_key  = _fmt_reskey(comp_id, asym_id, seq_id)

                is_hetatm = residue.entity_type == gemmi.EntityType.NonPolymer or \
                            residue.entity_type == gemmi.EntityType.Water or \
                            not residue.entity_type  # ligands / waters / unknowns

                for atom in residue:
                    alt = atom.altloc or "."
                    atom_dict = {
                        "auth_atom_id":  atom.name,
                        "label_alt_id":  alt,
                        "auth_comp_id":  comp_id,
                        "auth_asym_id":  asym_id,
                        "auth_seq_id":   str(seq_id),
                        "Cartn_x":       str(atom.pos.x),
                        "Cartn_y":       str(atom.pos.y),
                        "Cartn_z":       str(atom.pos.z),
                        "occupancy":     str(atom.occ),
                        "B_iso_or_equiv": str(atom.b_iso),
                        "type_symbol":   atom.element.name,
                        "group_PDB":     "HETATM" if is_hetatm else "ATOM",
                        "id":            str(atom.serial),
                    }

                    pdb_atom = PdbAtom(atom_dict)
                    natoms += pdb_atom.occupancy

                    if not is_hetatm:
                        # Protein / nucleic acid atom
                        if inner_distance:
                            res_atom_dict.setdefault(res_key, set()).add(pdb_atom)
                    else:
                        if comp_id == "HOH":
                            continue
                        if comp_id in cofactors.ligand_blacklist or comp_id in cofactors.metals:
                            res_atom_dict.setdefault(res_key, set()).add(pdb_atom)
                            notligands[res_key] = "Blacklisted ligand"
                        else:
                            ligand_res_atom_dict.setdefault(res_key, set()).add(pdb_atom)

    # ── struct_conn (covalent links) ───────────────────────────────────────
    # gemmi exposes connections on the structure object
    for conn in structure.connections:
        if conn.type not in (gemmi.ConnectionType.Covale, gemmi.ConnectionType.Disulf,
                             gemmi.ConnectionType.MetalC):
            continue
        p1 = conn.partner1
        p2 = conn.partner2

        res1 = _fmt_reskey(p1.res_id.name, p1.chain_name, p1.res_id.seqid.num)
        res2 = _fmt_reskey(p2.res_id.name, p2.chain_name, p2.res_id.seqid.num)

        # gemmi doesn't give us a pre-computed distance from the connection record;
        # use 1714 (unknown) as the original code did for "?" values.
        links.append((res1, res2, 1714.0))

    return natoms, res_atom_dict, ligand_res_atom_dict, notligands, links


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classificate_residue(residue, residue_dict, struc_dict, good_rsr, dubious_rsr, bad_rsr, cfg):
    score = 0
    reason = None

    if not residue_dict:
        score += 1000
        reason = f"No data for {residue}"
    else:
        if cfg.rscc_min > residue_dict["RSCC"]:
            score += 1
        if cfg.check_owab:
            owab = residue_dict.get("OWAB", 0)
            if not (1 < owab < cfg.owab_max):
                score += 1
        occ = residue_dict["occupancy"]
        if occ > 1.00:
            score += 1000
            reason = "Occupancy above 1"
        elif occ < cfg.occupancy_min:
            score += 1
        rsr = residue_dict["RSR"]
        if rsr > cfg.rsr_lower:
            score += 1
            if rsr > cfg.rsr_upper:
                score += 1

    if struc_dict:
        rfree = struc_dict.get("rFree", 9999)
        if rfree > cfg.rfree_max:
            score += 1
        if rfree < 0:
            score += 1000
            reason = f"No rFree data for {residue}"
        if cfg.check_resolution:
            resolution = struc_dict.get("Resolution", 10)
            if resolution > cfg.resolution_max:
                score += 1
            if resolution == 10:
                score += 1000
                reason = f"No resolution data for {residue}"
        if cfg.use_rdiff:
            rdiff = struc_dict.get("Rdiff", float("nan"))
            if math.isnan(rdiff):
                score += 1000
                reason = f"No reliable rFree/rWork data for {residue}"
            elif rdiff > cfg.rdiff_max:
                score += 1
        if cfg.use_dpi:
            dpi_val = struc_dict.get("DPI", -1)
            if math.isnan(dpi_val) or dpi_val <= 0:
                score += 1000
                reason = f"No reliable structural data for {residue}"
            elif dpi_val >= cfg.dpi_max:
                score += 1
    else:
        if cfg.use_dpi or cfg.use_rdiff or cfg.check_resolution:
            score += 1000
            reason = f"No structural data for {residue}"

    if score == 0:
        good_rsr.add(residue)
    elif score > cfg.tolerance:
        bad_rsr.add(residue)
    else:
        dubious_rsr.add(residue)

    return score, reason


def validate(residues, good_rsr, bad_rsr, dubious_rsr):
    if residues <= good_rsr:
        return "Good"
    if residues & bad_rsr:
        return "Bad"
    if residues & dubious_rsr:
        return "Dubious"
    return "Dubious"


# ---------------------------------------------------------------------------
# Ligand grouping
# ---------------------------------------------------------------------------

def group_ligands(ligand_residues, links):
    ligands = []
    linked_ligand_res = set()
    ligand_links = []
    for res1, res2, blen in links:
        if res1 in ligand_residues and res2 in ligand_residues:
            ligand_links.append((res1, res2, blen))
            linked_ligand_res.update([res1, res2])

    while ligand_links:
        for res1, res2, blen in ligand_links:
            for ligand in ligands:
                if res1 in ligand:
                    ligand.add(res2)
                    ligand_links.remove((res1, res2, blen))
                    break
                elif res2 in ligand:
                    ligand.add(res1)
                    ligand_links.remove((res1, res2, blen))
                    break
            else:
                ligands.append({res1, res2})
                ligand_links.remove((res1, res2, blen))
                break

    for lres in ligand_residues:
        if not any(lres in lig for lig in ligands):
            ligands.append({lres})

    # Merge overlapping sets
    merged = True
    while merged:
        merged = False
        for i, lig1 in enumerate(ligands):
            for j, lig2 in enumerate(ligands):
                if i >= j:
                    continue
                if lig1 & lig2:
                    ligands[i] = lig1 | lig2
                    ligands.pop(j)
                    merged = True
                    break
            if merged:
                break

    return ligands


# ---------------------------------------------------------------------------
# Binding-site extraction
# ---------------------------------------------------------------------------

def get_binding_site(ligand, ligand_score, good_rsr, bad_rsr, dubious_rsr,
                     pdbid, res_atom_dict, ligands, ligand_res_atom_dict,
                     edd_dict, struc_dict, notligands, cfg):
    inner_binding_site = set()
    for ligandres in ligand:
        if ligandres in notligands:
            return [notligands[ligandres]]
        for res, atoms in res_atom_dict.items():
            for atom in atoms:
                for ligandatom in ligand_res_atom_dict[ligandres]:
                    dist = atom | ligandatom
                    if dist <= cfg.inner_distance:
                        if dist < 2.1:
                            hetid = ligandres[:3].strip()
                            if hetid in cofactors.ligand_blacklist:
                                reason = "Covalently bound to a blacklisted ligand"
                                notligands[ligandres] = reason
                                return [reason]
                            elif hetid in cofactors.metals:
                                reason = "Covalently bound to the sequence"
                                notligands[ligandres] = reason
                                return [reason]
                        inner_binding_site.add(atom.residue)
                        break
        for other_lig in ligands:
            if other_lig == ligand:
                continue
            for lres in other_lig:
                for latom in ligand_res_atom_dict[lres]:
                    for ligandatom in ligand_res_atom_dict[ligandres]:
                        if (latom | ligandatom) <= cfg.inner_distance:
                            inner_binding_site.add(lres)
                            break

    bad_occupancy = [lr for lr in ligand
                     if edd_dict.get(lr, {"occupancy": 0})["occupancy"] < 1]
    bs_score = 0
    for res in inner_binding_site:
        resatoms = res_atom_dict.get(res)
        residue_dict = edd_dict.get(res)
        if not (residue_dict and resatoms):
            bad_occupancy.append(res)
        elif residue_dict.get("occupancy", 1) < 1:
            bad_occupancy.append(res)
        score, _ = classificate_residue(res, residue_dict, struc_dict,
                                        good_rsr, dubious_rsr, bad_rsr, cfg)
        bs_score = max(bs_score, score)

    rte = inner_binding_site | ligand - good_rsr
    ligandgood = validate(ligand, good_rsr, bad_rsr, dubious_rsr)
    bsgood = validate(inner_binding_site, good_rsr, bad_rsr, dubious_rsr)
    return ligand, inner_binding_site, rte, ligandgood, bsgood, bad_occupancy, ligand_score, bs_score


# ---------------------------------------------------------------------------
# Main per-PDB entry point
# ---------------------------------------------------------------------------

def parse_binding_site(pdbid, cfg=None):
    """
    Analyse a single PDB entry.

    Returns a dict:
      On success:
        {"pdbid": ..., "ligands": [...], "rejected": {...}, "struc_dict": {...}}
      On failure:
        {"pdbid": ..., "error": "reason"}
    """
    if cfg is None:
        cfg = AnalysisConfig()

    pdbid = pdbid.lower()
    logger.info("Analysing %s", pdbid)

    # --- Fetch structure-level stats ---
    if cfg.pdb_redo:
        pdbid_stats = pdb_redo_utils.get_pdbredo_data(pdbid)
        if not pdbid_stats:
            return {"pdbid": pdbid, "error": "Not in PDB_REDO"}
        edd_dict = pdb_redo_utils.get_ED_data(pdbid)
        if not edd_dict:
            return {"pdbid": pdbid, "error": "No PDB-REDO ED data available"}
    else:
        report = pdb_utils.get_custom_report(pdbid)
        # report may be empty for NMR/cryo-EM — continue with empty stats;
        # rFree will be 9999, adding to score, but analysis still runs.
        pdbid_stats = report.get(pdbid.upper(), {}) if report else {}
        if not pdbid_stats:
            logger.warning("%s: no refinement stats; proceeding without R-factor data", pdbid)
        pdbdict, edd_dict = eds_utils.get_EDS(pdbid)
        if not edd_dict:
            return {"pdbid": pdbid, "error": "No EDM/validation data available (may not be an X-ray entry)"}

    struc_dict = {
        "rFree": pdbid_stats.get("rFree", float("nan")),
        "rWork": pdbid_stats.get("rWork", float("nan")),
    }
    if cfg.check_resolution:
        struc_dict["Resolution"] = pdbid_stats.get("refinementResolution", 10)
    if cfg.use_rdiff:
        struc_dict["Rdiff"] = struc_dict["rFree"] - struc_dict["rWork"]

    # --- Download & parse mmCIF ---
    pdbfilepath = pdb_utils.get_pdb_file(pdbid.upper(), cfg.pdb_redo)
    if not pdbfilepath:
        return {"pdbid": pdbid, "error": "Unable to load PDBx/mmCIF model"}

    parsed = parse_mmcif_file(pdbfilepath, pdbid, cfg.inner_distance)
    if len(parsed) == 1:
        return {"pdbid": pdbid, "error": str(parsed[0])}

    natoms, res_atom_dict, ligand_res_atom_dict, notligands, links = parsed

    if cfg.use_dpi:
        rfl = pdbid_stats.get("nreflections", 0)
        if rfl:
            struc_dict["DPI"] = _dpi(
                pdbid_stats.get("lengthOfUnitCellLatticeA", 0),
                pdbid_stats.get("lengthOfUnitCellLatticeB", 0),
                pdbid_stats.get("lengthOfUnitCellLatticeC", 0),
                pdbid_stats.get("unitCellAngleAlpha", 0),
                pdbid_stats.get("unitCellAngleBeta", 0),
                pdbid_stats.get("unitCellAngleGamma", 0),
                natoms, rfl, struc_dict["rFree"],
            )
        else:
            struc_dict["DPI"] = float("nan")

    # --- Fill occupancy gaps ---
    bad_res = set()
    for residue, rdict in edd_dict.items():
        residue = residue.strip()
        if "occupancy" not in rdict:
            atoms = res_atom_dict.get(residue) or ligand_res_atom_dict.get(residue)
            if atoms:
                rdict["occupancy"] = _average_occ(atoms)
            else:
                bad_res.add(residue)
    edd_dict = {k.strip(): v for k, v in edd_dict.items() if k.strip() not in bad_res}

    # --- Prune covalently bound ligands ---
    all_links_parsed = False
    while not all_links_parsed:
        for res1, res2, blen in list(links):
            checklink = sum([res1 in res_atom_dict, res2 in res_atom_dict])
            if checklink == 2:
                links.remove((res1, res2, blen))
                break
            sres = ligres = None
            if res1 in res_atom_dict:
                sres, ligres = res1, res2
            elif res2 in res_atom_dict:
                sres, ligres = res2, res1
            if checklink == 1:
                if blen and blen >= 2.1:
                    continue
                if (res1[:3].strip() in cofactors.metals) or (res2[:3].strip() in cofactors.metals):
                    continue
                if ligres:
                    if (res1[:3].strip() in cofactors.ligand_blacklist) or \
                       (res2[:3].strip() in cofactors.ligand_blacklist):
                        notligands[ligres] = "Covalently bound to a blacklisted ligand"
                    else:
                        notligands[ligres] = "Covalently bound to the sequence"
                    links.remove((res1, res2, blen))
                    if ligres in ligand_res_atom_dict:
                        res_atom_dict[ligres] = ligand_res_atom_dict.pop(ligres)
                    break
        else:
            all_links_parsed = True

    if not ligand_res_atom_dict:
        return {"pdbid": pdbid, "error": "No ligands found"}

    good_rsr, dubious_rsr, bad_rsr = set(), set(), set()
    ligands = group_ligands(ligand_res_atom_dict.keys(), links)

    # Score individual ligand residues
    ligand_scores = []
    for ligand in ligands:
        ligand_score = 0
        for res in list(ligand):
            residue_dict = edd_dict.get(res)
            score, reason = classificate_residue(
                res, residue_dict, struc_dict, good_rsr, dubious_rsr, bad_rsr, cfg)
            if reason and score >= 1000:
                notligands[res] = reason
                ligand.discard(res)
            ligand_score = max(ligand_score, score)
        ligand_scores.append(ligand_score)

    ligand_bs_list = []
    for ligand, ligand_score in zip(ligands, ligand_scores):
        if not ligand:
            continue
        bs = get_binding_site(
            ligand, ligand_score, good_rsr, bad_rsr, dubious_rsr,
            pdbid, res_atom_dict, ligands, ligand_res_atom_dict,
            edd_dict, struc_dict, notligands, cfg)
        if len(bs) == 1:
            for lr in ligand:
                res_atom_dict[lr] = ligand_res_atom_dict.pop(lr, set())
                notligands.setdefault(lr, bs[0])
            continue
        ligand_bs_list.append(bs)

    # --- Serialise to JSON-safe structures ---
    result_ligands = []
    source = "PDB_REDO" if cfg.pdb_redo else "PDB"
    for data in ligand_bs_list:
        ligandresidues, binding_site, rte, ligandgood, bsgood, bad_occupancy, lig_score, bs_score = data
        if not ligandresidues:
            continue
        result_ligands.append({
            "ligand_residues": sorted(ligandresidues),
            "binding_site_residues": sorted(binding_site),
            "residues_to_examine": sorted(rte),
            "ligand_quality": ligandgood,
            "binding_site_quality": bsgood,
            "source": source,
            "ligand_score": lig_score,
            "binding_site_score": bs_score,
            "low_occupancy": sorted(bad_occupancy),
        })

    safe_struc = {k: (v if not (isinstance(v, float) and math.isnan(v)) else None)
                  for k, v in struc_dict.items()}

    return {
        "pdbid": pdbid,
        "ligands": result_ligands,
        "rejected": {k: str(v) for k, v in notligands.items()},
        "struc_dict": safe_struc,
    }


def analyse_pdbids(pdbids, cfg=None):
    """
    Analyse a list of PDB IDs sequentially and return a list of result dicts.
    The web layer may call this in a background thread.
    """
    if cfg is None:
        cfg = AnalysisConfig()
    results = []
    for pdbid in pdbids:
        try:
            result = parse_binding_site(pdbid.strip().lower(), cfg)
        except Exception as exc:
            logger.exception("Unexpected error analysing %s", pdbid)
            result = {"pdbid": pdbid, "error": str(exc)}
        results.append(result)
    return results
