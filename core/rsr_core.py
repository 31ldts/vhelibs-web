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
import concurrent.futures

import gemmi

from core.pdb_atom import PdbAtom, format_reskey
import core.pdb_utils as pdb_utils
import core.eds_utils as eds_utils
import core.pdb_redo_utils as pdb_redo_utils
import core.cofactors as cofactors

logger = logging.getLogger(__name__)


class AnalysisConfig:
    """All tunable parameters for a single analysis run.

    Attributes:
        rsr_upper (float): Upper RSR (real-space R-factor) threshold above
            which a residue accrues an extra penalty point.
        rsr_lower (float): Lower RSR threshold above which a residue
            accrues a penalty point.
        rscc_min (float): Minimum acceptable RSCC (real-space correlation
            coefficient); values below this accrue a penalty point.
        rfree_max (float): Maximum acceptable structure-level R-free
            value; values above this accrue a penalty point.
        occupancy_min (float): Minimum acceptable per-residue occupancy;
            values below this (but not above 1.0) accrue a penalty point.
        tolerance (int): Maximum cumulative score a residue may have while
            still being classified as "dubious" rather than "bad".
        inner_distance (float): Squared contact distance (Angstrom²) used
            to determine binding-site membership, derived from the
            ``distance`` constructor argument.
        check_owab (bool): Whether to factor per-atom OWAB (occupancy-
            weighted average B-factor) into residue scoring.
        owab_max (float): Maximum acceptable OWAB value when
            ``check_owab`` is enabled.
        check_resolution (bool): Whether to factor structure resolution
            into scoring.
        resolution_max (float): Maximum acceptable resolution (Angstrom)
            when ``check_resolution`` is enabled.
        use_rdiff (bool): Whether to factor the R-free/R-work difference
            into scoring.
        rdiff_max (float): Maximum acceptable R-free/R-work difference
            when ``use_rdiff`` is enabled.
        use_dpi (bool): Whether to factor the estimated coordinate
            precision (DPI) into scoring.
        dpi_max (float): Maximum acceptable DPI value when ``use_dpi`` is
            enabled.
        pdb_redo (bool): Whether to source structural/refinement data from
            PDB-REDO instead of the primary PDB/RCSB archive.
    """

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
        """Initialize an AnalysisConfig with the given tunable parameters.

        Args:
            rsr_upper (float, optional): Upper RSR threshold. Defaults to
                ``0.4``.
            rsr_lower (float, optional): Lower RSR threshold. Defaults to
                ``0.24``.
            rscc_min (float, optional): Minimum acceptable RSCC. Defaults
                to ``0.9``.
            rfree_max (float, optional): Maximum acceptable R-free.
                Defaults to ``1.0``.
            occupancy_min (float, optional): Minimum acceptable
                occupancy. Defaults to ``1.0``.
            tolerance (int, optional): Score tolerance separating
                "dubious" from "bad" residues. Defaults to ``2``.
            distance (float, optional): Contact distance in Angstrom used
                to compute ``inner_distance`` (its square). Defaults to
                ``4.5``.
            check_owab (bool, optional): Whether to check OWAB. Defaults
                to ``False``.
            owab_max (float, optional): Maximum acceptable OWAB. Defaults
                to ``50.0``.
            check_resolution (bool, optional): Whether to check
                resolution. Defaults to ``False``.
            resolution_max (float, optional): Maximum acceptable
                resolution. Defaults to ``3.5``.
            use_rdiff (bool, optional): Whether to use the R-free/R-work
                difference. Defaults to ``False``.
            rdiff_max (float, optional): Maximum acceptable R-free/R-work
                difference. Defaults to ``0.05``.
            use_dpi (bool, optional): Whether to use DPI. Defaults to
                ``False``.
            dpi_max (float, optional): Maximum acceptable DPI. Defaults to
                ``0.42``.
            pdb_redo (bool, optional): Whether to source data from
                PDB-REDO. Defaults to ``False``.

        Returns:
            None

        Raises:
            None
        """
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
    """Compute the average occupancy over a collection of atoms.

    Args:
        residue_atoms (iterable): Iterable of :class:`core.pdb_atom.PdbAtom`
            objects (typically all atoms of a single residue).

    Returns:
        float: The mean of the ``occupancy`` attribute across
        ``residue_atoms``.

    Raises:
        ZeroDivisionError: If ``residue_atoms`` is empty.
    """
    return sum(a.occupancy for a in residue_atoms) / len(residue_atoms)


# ---------------------------------------------------------------------------
# Density-map bounding boxes
# ---------------------------------------------------------------------------
#
# For each region of interest (ligand, binding site, residues-to-examine)
# we compute an axis-aligned bounding box over the atoms of that region,
# padded by `padding` Angstrom. The frontend uses this box to ask EBI's
# density server (see core.eds_utils.edm_box_url) for just the relevant
# chunk of the 2Fo-Fc/Fo-Fc map, instead of downloading and masking the
# whole map client-side. `padding` defaults to 2.1 Å, the same distance
# rsr_analysis.py/rsr_core.py already use elsewhere to detect covalent
# contacts, and a sensible "shell" around the region for density display.

DEFAULT_BOX_PADDING = 2.1


def _atoms_for_residues(residues, res_atom_dict, ligand_res_atom_dict):
    """Collect all atoms belonging to a set of residue keys.

    Looks up each residue key in both the protein/cofactor atom dict and
    the ligand atom dict, combining the results into a single flat list.

    Args:
        residues (iterable): Iterable of residue key strings (as produced
            by :func:`core.pdb_atom.format_reskey`) to collect atoms for.
        res_atom_dict (dict): Mapping of residue key to a collection of
            :class:`core.pdb_atom.PdbAtom` for protein/nucleic-acid and
            blacklisted/metal atoms.
        ligand_res_atom_dict (dict): Mapping of residue key to a
            collection of :class:`core.pdb_atom.PdbAtom` for ligand
            atoms.

    Returns:
        list: A flat list of :class:`core.pdb_atom.PdbAtom` objects
        belonging to any of the given residue keys.

    Raises:
        None
    """
    atoms = []
    for res in residues:
        atoms.extend(res_atom_dict.get(res, ()))
        atoms.extend(ligand_res_atom_dict.get(res, ()))
    return atoms


def _bbox(atoms, padding=DEFAULT_BOX_PADDING):
    """Compute an axis-aligned, padded bounding box over a set of atoms.

    Args:
        atoms (iterable): Iterable of :class:`core.pdb_atom.PdbAtom`
            objects to compute the bounding box over.
        padding (float, optional): Amount, in Angstrom, to expand the box
            on every side beyond the atoms' extreme coordinates. Defaults
            to :data:`DEFAULT_BOX_PADDING`.

    Returns:
        dict or None: A dict of the form
        ``{"min": [x, y, z], "max": [x, y, z]}`` describing the padded
        bounding box, or ``None`` if ``atoms`` is empty.

    Raises:
        None
    """
    atoms = list(atoms)
    if not atoms:
        return None
    xs, ys, zs = zip(*(a.xyz for a in atoms))
    return {
        "min": [min(xs) - padding, min(ys) - padding, min(zs) - padding],
        "max": [max(xs) + padding, max(ys) + padding, max(zs) + padding],
    }


def residues_bbox(residues, res_atom_dict, ligand_res_atom_dict, padding=DEFAULT_BOX_PADDING):
    """Compute a padded bounding box covering all atoms of given residues.

    This is only used to size the *download window* sent to the density
    server (which only supports box/cell queries, not per-atom masking) —
    it is intentionally loose, since the whole point is to have enough
    margin to render a good isosurface. It should NOT be used directly to
    decide what density is *shown*, since a box covering e.g. a binding
    site will inevitably also cover the ligand and unrelated solvent; see
    `residue_atom_centers` for the per-atom masks used for that.

    Args:
        residues (iterable): Iterable of residue key strings to include in
            the bounding box.
        res_atom_dict (dict): Mapping of residue key to a collection of
            :class:`core.pdb_atom.PdbAtom` for protein/nucleic-acid and
            blacklisted/metal atoms.
        ligand_res_atom_dict (dict): Mapping of residue key to a
            collection of :class:`core.pdb_atom.PdbAtom` for ligand
            atoms.
        padding (float, optional): Amount, in Angstrom, to expand the box
            on every side. Defaults to :data:`DEFAULT_BOX_PADDING`.

    Returns:
        dict or None: A dict of the form
        ``{"min": [x, y, z], "max": [x, y, z]}`` describing the padded
        bounding box, or ``None`` if no atoms are found for the given
        residues.

    Raises:
        None
    """
    atoms = _atoms_for_residues(residues, res_atom_dict, ligand_res_atom_dict)
    return _bbox(atoms, padding)


def residue_atom_centers(residues, res_atom_dict, ligand_res_atom_dict):
    """Build a flat list of per-atom coordinates for the given residue keys.

    Used for true per-atom density masking in the viewer: the frontend
    clips the isosurface to a small sphere around *each* atom (radius
    controlled by a UI slider), instead of one region-wide box or a single
    per-residue sphere. Since the ligand/binding-site/residues-to-examine
    atom lists are already known from the analysis, this lets the density
    shown for each layer trace the actual atoms of that layer, rather than
    everything inside a loose enclosing volume.

    Args:
        residues (iterable): Iterable of residue key strings to collect
            atom centers for.
        res_atom_dict (dict): Mapping of residue key to a collection of
            :class:`core.pdb_atom.PdbAtom` for protein/nucleic-acid and
            blacklisted/metal atoms.
        ligand_res_atom_dict (dict): Mapping of residue key to a
            collection of :class:`core.pdb_atom.PdbAtom` for ligand
            atoms.

    Returns:
        list: A list of dicts, one per atom, each of the form
        ``{"residue": res, "center": [x, y, z]}``.

    Raises:
        None
    """
    atoms_out = []
    for res in residues:
        atoms = res_atom_dict.get(res) or ligand_res_atom_dict.get(res) or ()
        for a in atoms:
            atoms_out.append({"residue": res, "center": list(a.xyz)})
    return atoms_out


def _dpi(a, b, c, alpha, beta, gamma, natoms, reflections, rfree):
    """Estimate the diffraction-component precision indicator (DPI).

    Computes the DPI using the unit-cell dimensions, number of atoms,
    number of reflections, and R-free value, following the standard DPI
    formula based on unit-cell volume.

    Args:
        a (float): Unit cell lattice length A, in Angstrom.
        b (float): Unit cell lattice length B, in Angstrom.
        c (float): Unit cell lattice length C, in Angstrom.
        alpha (float): Unit cell angle alpha, in degrees.
        beta (float): Unit cell angle beta, in degrees.
        gamma (float): Unit cell angle gamma, in degrees.
        natoms (float): Number of atoms (occupancy-weighted) in the
            structure.
        reflections (float): Number of reflections used in refinement.
        rfree (float): The structure's R-free value.

    Returns:
        float: The estimated DPI value, or ``float("nan")`` if the unit
        cell volume is non-positive or ``reflections`` is non-positive.

    Raises:
        None
    """
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
    """Build a residue key.

    Thin wrapper kept for call-site compatibility — see
    core.pdb_atom.format_reskey for the canonical (shared) implementation
    and why asym_id/seq_id must always be joined with an explicit space.

    Args:
        comp_id (str): Component (residue) identifier, e.g. the residue
            name.
        asym_id (str): Author-assigned asymmetric unit (chain) identifier.
        seq_id (int or str): Author-assigned sequence number of the
            residue.

    Returns:
        str: The formatted, canonical residue key, as produced by
        :func:`core.pdb_atom.format_reskey`.

    Raises:
        None
    """
    return format_reskey(comp_id, asym_id, seq_id)


def parse_mmcif_file(mmciffilepath, pdbid, inner_distance):
    """Parse an mmCIF file into atom and covalent-link dictionaries.

    Parses an mmCIF file (plain or gzip) using gemmi, classifying each
    atom as protein/nucleic-acid, blacklisted ligand/metal, or ligand,
    and extracting covalent/disulfide/metal-coordination connections.

    Args:
        mmciffilepath (str): Path to the mmCIF file to parse (may be
            gzip-compressed).
        pdbid (str): PDB identifier of the structure, used for logging.
        inner_distance (float): If falsy (e.g. ``0``), protein/nucleic-acid
            atoms are not added to ``res_atom_dict``, skipping binding-site
            distance calculations for non-ligand residues.

    Returns:
        tuple: On success, a 5-tuple
        ``(natoms, res_atom_dict, ligand_res_atom_dict, notligands, links)``:

            - natoms (float): Total occupancy-weighted atom count.
            - res_atom_dict (dict): Mapping of residue key to a set of
              :class:`core.pdb_atom.PdbAtom` for protein/nucleic-acid and
              blacklisted/metal atoms.
            - ligand_res_atom_dict (dict): Mapping of residue key to a set
              of :class:`core.pdb_atom.PdbAtom` for ligand atoms.
            - notligands (dict): Mapping of residue key to a reason string
              for residues excluded from ligand consideration (e.g.
              blacklisted).
            - links (list): List of ``(res1, res2, bond_length)`` tuples
              describing covalent/disulfide/metal-coordination
              connections.

        On failure, a 1-tuple ``(error_string,)`` describing the parsing
        error.

    Raises:
        None: Parsing errors from gemmi are caught internally and
            returned as a 1-tuple error message instead of propagating.
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
    """Score a residue's electron-density fit and structural quality.

    Accumulates a penalty score based on RSCC, occupancy, RSR, and
    (optionally) OWAB from ``residue_dict``, plus structure-level R-free,
    resolution, R-diff, and DPI checks from ``struc_dict`` as enabled by
    ``cfg``. Based on the final score, the residue is added to exactly one
    of ``good_rsr``, ``dubious_rsr``, or ``bad_rsr``.

    Args:
        residue (str): Residue key of the residue being scored, used only
            for building reason strings.
        residue_dict (dict or None): Per-residue validation stats with
            keys such as ``"RSCC"``, ``"RSR"``, ``"occupancy"``, and
            optionally ``"OWAB"``. If falsy, the residue is scored as
            having no data.
        struc_dict (dict or None): Structure-level stats with keys such as
            ``"rFree"``, ``"Resolution"``, ``"Rdiff"``, and ``"DPI"``. If
            falsy, structure-dependent checks enabled in ``cfg`` add a
            large penalty.
        good_rsr (set): Set of residue keys classified as "good"; updated
            in place.
        dubious_rsr (set): Set of residue keys classified as "dubious";
            updated in place.
        bad_rsr (set): Set of residue keys classified as "bad"; updated in
            place.
        cfg (AnalysisConfig): Analysis configuration providing thresholds
            and which optional checks to perform.

    Returns:
        tuple: A 2-tuple ``(score, reason)`` where ``score`` (int) is the
        cumulative penalty score and ``reason`` (str or None) is a
        human-readable explanation set when a severe (>=1000) penalty was
        applied, or ``None`` otherwise.

    Raises:
        KeyError: If ``residue_dict`` is truthy but missing the
            ``"occupancy"`` or ``"RSR"`` key.
    """
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
    """Determine an overall quality label for a set of residues.

    A group of residues is labeled "Good" only if every residue is in
    ``good_rsr``; otherwise it is "Bad" if any residue is in ``bad_rsr``,
    "Dubious" if any residue is in ``dubious_rsr``, and "Dubious" as a
    final fallback.

    Args:
        residues (set): Set of residue keys to evaluate.
        good_rsr (set): Set of residue keys classified as "good".
        bad_rsr (set): Set of residue keys classified as "bad".
        dubious_rsr (set): Set of residue keys classified as "dubious".

    Returns:
        str: One of ``"Good"``, ``"Bad"``, or ``"Dubious"``.

    Raises:
        None
    """
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
    """Group ligand residue keys into connected ligand components.

    Residues connected by a covalent link that are both ligand residues
    are merged into the same group. Any remaining residue not linked to
    another becomes its own singleton group. Overlapping groups are then
    iteratively merged until none overlap.

    Args:
        ligand_residues (iterable): Iterable of ligand residue key
            strings.
        links (list): List of ``(res1, res2, bond_length)`` tuples
            describing covalent connections between residues, as returned
            by :func:`parse_mmcif_file`.

    Returns:
        list: A list of sets, each set containing the residue keys that
        belong to the same connected ligand component.

    Raises:
        None
    """
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
    """Compute the binding site and quality assessment for a ligand group.

    Identifies all non-ligand residues (and residues from other ligand
    groups) within contact distance of the given ligand's atoms, checking
    for disqualifying covalent bonds to blacklisted ligands or metals
    along the way. Scores the binding site residues and determines overall
    quality labels for both the ligand and its binding site.

    Args:
        ligand (set): Set of residue keys belonging to the ligand group
            being analyzed.
        ligand_score (int): Precomputed maximum per-residue score for this
            ligand (from the caller's scoring pass).
        good_rsr (set): Set of residue keys classified as "good"; may be
            updated via :func:`classificate_residue`.
        bad_rsr (set): Set of residue keys classified as "bad"; may be
            updated via :func:`classificate_residue`.
        dubious_rsr (set): Set of residue keys classified as "dubious";
            may be updated via :func:`classificate_residue`.
        pdbid (str): PDB identifier of the structure (currently unused
            within the function body but kept for call-site
            compatibility).
        res_atom_dict (dict): Mapping of residue key to a collection of
            :class:`core.pdb_atom.PdbAtom` for protein/nucleic-acid and
            blacklisted/metal atoms.
        ligands (list): List of all ligand groups (sets of residue keys)
            in the structure, used to detect proximity to other ligands.
        ligand_res_atom_dict (dict): Mapping of residue key to a
            collection of :class:`core.pdb_atom.PdbAtom` for ligand atoms.
        edd_dict (dict): Mapping of residue key to per-residue validation
            stats (RSR, RSCC, occupancy, etc.).
        struc_dict (dict): Structure-level stats (rFree, resolution,
            etc.) used for scoring.
        notligands (dict): Mapping of residue key to disqualification
            reason strings; updated in place if a covalent disqualifying
            bond is found.
        cfg (AnalysisConfig): Analysis configuration providing thresholds
            and which optional checks to perform.

    Returns:
        list or tuple: If the ligand is disqualified due to a covalent
        bond to a blacklisted ligand or metal, a 1-element list
        ``[reason_string]``. Otherwise, an 8-tuple:

            - ligand (set): The input ligand residue set (possibly
              unchanged).
            - inner_binding_site (set): Residue keys within contact
              distance of the ligand.
            - rte (set): "Residues to examine" — the union of the binding
              site and ligand residues not already classified as good.
            - ligandgood (str): Quality label for the ligand ("Good",
              "Bad", or "Dubious").
            - bsgood (str): Quality label for the binding site.
            - bad_occupancy (list): Residue keys with occupancy below 1 or
              missing data.
            - ligand_score (int): The input ``ligand_score``, passed
              through unchanged.
            - bs_score (int): Maximum per-residue score among the binding
              site residues.

    Raises:
        KeyError: If a ligand residue key is missing from
            ``ligand_res_atom_dict``.
    """
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

    rte = (inner_binding_site | ligand) - good_rsr
    ligandgood = validate(ligand, good_rsr, bad_rsr, dubious_rsr)
    bsgood = validate(inner_binding_site, good_rsr, bad_rsr, dubious_rsr)
    return ligand, inner_binding_site, rte, ligandgood, bsgood, bad_occupancy, ligand_score, bs_score


# ---------------------------------------------------------------------------
# Main per-PDB entry point
# ---------------------------------------------------------------------------

def parse_binding_site(pdbid, cfg=None):
    """Analyse a single PDB entry for ligand binding-site quality.

    Fetches structure-level refinement stats and per-residue validation
    data (from either PDB/RCSB or PDB-REDO, per ``cfg.pdb_redo``),
    downloads and parses the mmCIF model, then classifies each ligand and
    its binding site by electron-density fit quality, producing
    JSON-serializable results including bounding boxes and per-atom
    coordinates for 3D density display.

    Args:
        pdbid (str): PDB identifier of the structure to analyse. Case is
            normalized internally.
        cfg (AnalysisConfig, optional): Analysis configuration. If
            ``None``, a default :class:`AnalysisConfig` is used.

    Returns:
        dict: On success, a dict with keys:

            - ``"pdbid"`` (str): The analysed PDB identifier.
            - ``"ligands"`` (list): List of per-ligand result dicts, each
              containing ligand/binding-site residues, quality labels,
              scores, density boxes, and density atom coordinates.
            - ``"rejected"`` (dict): Mapping of residue key to the reason
              it was excluded from ligand consideration.
            - ``"struc_dict"`` (dict): Structure-level statistics with
              NaN values replaced by ``None`` for JSON safety.

        On failure, a dict of the form
        ``{"pdbid": pdbid, "error": "reason"}``.

    Raises:
        None: Fetch and parsing errors are caught internally (or via the
            helper functions called) and surfaced through the returned
            ``"error"`` key rather than propagating exceptions.
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

    # Snapshot of every residue belonging to *any* ligand group in this
    # structure, taken after the per-residue pruning above. Used below to
    # keep each ligand's binding-site/residues-to-examine lists — and the
    # density boxes/atoms derived from them — free of atoms belonging to a
    # *different* ligand. Without this, a second ligand sitting close
    # enough to be pulled into get_binding_site()'s cross-ligand check
    # would end up rendered in the 3D viewer even though it isn't the
    # ligand actually under study for that entry.
    all_ligand_keys = set().union(*ligands) if ligands else set()

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

        # Residues belonging to any *other* ligand in this structure (e.g. a
        # second binding site close enough to have been pulled into this
        # one's binding site by get_binding_site's cross-ligand check).
        # These still count towards this ligand's scoring above, but are
        # excluded from everything sent below for 3D display, so the viewer
        # only ever shows the ligand actually under study for this entry —
        # any other ligand present in the structure stays hidden.
        other_ligand_residues = all_ligand_keys - set(ligandresidues)
        display_binding_site = [r for r in binding_site if r not in other_ligand_residues]
        display_rte = [r for r in rte if r not in other_ligand_residues]
        display_bad_occupancy = [r for r in bad_occupancy if r not in other_ligand_residues]

        result_ligands.append({
            "ligand_residues": sorted(ligandresidues),
            "binding_site_residues": sorted(display_binding_site),
            "residues_to_examine": sorted(display_rte),
            "ligand_quality": ligandgood,
            "binding_site_quality": bsgood,
            "source": source,
            "ligand_score": lig_score,
            "binding_site_score": bs_score,
            "low_occupancy": sorted(display_bad_occupancy),
            # Other ligand(s) present in this structure that are NOT shown
            # in the 3D viewer for this entry. Purely informational for the
            # UI (e.g. "2 other ligand(s) hidden") — not used for scoring.
            "other_ligands": sorted(other_ligand_residues),
            # Padded bounding boxes for on-demand, segmented density
            # display in the 3D viewer (see core.eds_utils.edm_box_url).
            # None if a region has no atoms (shouldn't normally happen for
            # ligand_residues, but binding_site/rte can theoretically be
            # empty for a solvent-exposed ligand with distance=0).
            "density_boxes": {
                "ligand": residues_bbox(ligandresidues, res_atom_dict, ligand_res_atom_dict),
                "binding_site": residues_bbox(display_binding_site, res_atom_dict, ligand_res_atom_dict),
                "residues_to_examine": residues_bbox(display_rte, res_atom_dict, ligand_res_atom_dict),
            },
            # Per-atom coordinates used to actually differentiate the
            # density shown for each region (see residue_atom_centers
            # docstring) — the boxes above only size the download window.
            "density_atoms": {
                "ligand": residue_atom_centers(ligandresidues, res_atom_dict, ligand_res_atom_dict),
                "binding_site": residue_atom_centers(display_binding_site, res_atom_dict, ligand_res_atom_dict),
                "residues_to_examine": residue_atom_centers(display_rte, res_atom_dict, ligand_res_atom_dict),
            },
        })

    safe_struc = {k: (v if not (isinstance(v, float) and math.isnan(v)) else None)
                  for k, v in struc_dict.items()}

    return {
        "pdbid": pdbid,
        "ligands": result_ligands,
        "rejected": {k: str(v) for k, v in notligands.items()},
        "struc_dict": safe_struc,
    }



# Default cap on how many PDB entries are analysed concurrently. Each entry's
# runtime is dominated by network I/O (mmCIF + validation/stats downloads
# against RCSB/EBI/PDB-REDO), which releases the GIL while waiting, so a
# modest thread pool gives a large wall-clock speedup for multi-entry jobs.
# Kept fairly conservative to avoid tripping rate-limiting (HTTP 429) on the
# upstream APIs, which don't publish official concurrency limits.
DEFAULT_ANALYSE_WORKERS = 8


def analyse_pdbids(pdbids, cfg=None, max_workers=DEFAULT_ANALYSE_WORKERS):
    """Analyse a list of PDB IDs concurrently.

    Calls :func:`parse_binding_site` for each PDB ID via a thread pool
    (since each call is I/O-bound, dominated by network downloads rather
    than CPU), catching and logging any unexpected exception per entry so
    that one failing entry does not abort the whole batch. The web layer
    may call this in a background thread; results preserve the input
    order regardless of completion order.

    Args:
        pdbids (iterable): Iterable of PDB identifier strings to analyse.
            Each is stripped of whitespace and lower-cased before
            analysis.
        cfg (AnalysisConfig, optional): Analysis configuration shared
            across all entries. If ``None``, a default
            :class:`AnalysisConfig` is used. ``AnalysisConfig`` instances
            are read-only during analysis, so sharing one across threads
            is safe.
        max_workers (int, optional): Maximum number of PDB entries
            analysed concurrently. Defaults to
            :data:`DEFAULT_ANALYSE_WORKERS`. Each PDB ID writes only to
            its own per-entry cache directory (see
            :func:`core.http_cache.entry_cache_dir`), so entries never
            contend with each other on disk.

    Returns:
        list: A list of result dicts, one per input PDB ID, **in the same
        order as the input** (not completion order), in the same format
        returned by :func:`parse_binding_site` (either a success dict or
        an ``{"pdbid": ..., "error": ...}`` dict).

    Raises:
        None: Unexpected exceptions from analysing an individual entry are
            caught, logged, and converted into an error result for that
            entry rather than propagating.
    """
    if cfg is None:
        cfg = AnalysisConfig()
    pdbids = list(pdbids)
    results = [None] * len(pdbids)
    if not pdbids:
        return results

    def _analyse_one(pdbid):
        try:
            return parse_binding_site(pdbid.strip().lower(), cfg)
        except Exception as exc:
            logger.exception("Unexpected error analysing %s", pdbid)
            return {"pdbid": pdbid, "error": str(exc)}

    workers = min(max_workers, len(pdbids))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_index = {
            executor.submit(_analyse_one, pdbid): i
            for i, pdbid in enumerate(pdbids)
        }
        for future in concurrent.futures.as_completed(future_to_index):
            results[future_to_index[future]] = future.result()

    return results
