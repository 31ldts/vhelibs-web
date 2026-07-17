# -*- coding: utf-8 -*-
"""
Unit tests for core.rsr_core.

External dependencies mocked/isolated:
  - gemmi.read_structure: monkeypatched to return synthetic gemmi.Structure
    objects built in-memory (via the gemmi Python API), so no real mmCIF
    files or downloads are needed.
  - core.pdb_utils / core.eds_utils / core.pdb_redo_utils: patched at the
    call sites used by core.rsr_core.parse_binding_site.
  - core.rsr_core.parse_mmcif_file itself is patched for the
    parse_binding_site/analyse_pdbids orchestration tests, since that
    parsing logic is already covered directly by TestParseMmcifFile.
"""
import math
from unittest.mock import patch

import gemmi
import pytest

import core.rsr_core as rsr_core
from core.pdb_atom import format_reskey


# ---------------------------------------------------------------------------
# gemmi structure-building helpers (test fixtures, not production code)
# ---------------------------------------------------------------------------

def make_gemmi_atom(name, x, y, z, occ=1.0, b_iso=20.0, element="C", serial=1):
    a = gemmi.Atom()
    a.name = name
    a.pos = gemmi.Position(x, y, z)
    a.occ = occ
    a.b_iso = b_iso
    a.element = gemmi.Element(element)
    a.serial = serial
    a.altloc = "\0"
    return a


def make_gemmi_residue(name, seqnum, entity_type, atoms):
    r = gemmi.Residue()
    r.name = name
    r.seqid = gemmi.SeqId(seqnum, " ")
    r.entity_type = entity_type
    for a in atoms:
        r.add_atom(a)
    return r


def make_gemmi_structure(chains):
    """chains: dict[chain_name] -> list[gemmi.Residue]"""
    model = gemmi.Model("1")
    for chain_name, residues in chains.items():
        chain = gemmi.Chain(chain_name)
        for res in residues:
            chain.add_residue(res)
        model.add_chain(chain)
    structure = gemmi.Structure()
    structure.add_model(model)
    return structure


def make_connection(conn_type, chain1, resname1, seqnum1, atom1,
                     chain2, resname2, seqnum2, atom2):
    conn = gemmi.Connection()
    conn.type = conn_type
    conn.partner1 = gemmi.AtomAddress(
        chain1, gemmi.SeqId(seqnum1, " "), resname1, atom1)
    conn.partner2 = gemmi.AtomAddress(
        chain2, gemmi.SeqId(seqnum2, " "), resname2, atom2)
    return conn


# ---------------------------------------------------------------------------
# AnalysisConfig
# ---------------------------------------------------------------------------

class TestAnalysisConfig:
    def test_defaults(self):
        cfg = rsr_core.AnalysisConfig()
        assert cfg.rsr_upper == 0.4
        assert cfg.rsr_lower == 0.24
        assert cfg.rscc_min == 0.9
        assert cfg.rfree_max == 1.0
        assert cfg.occupancy_min == 1.0
        assert cfg.tolerance == 2
        assert cfg.inner_distance == pytest.approx(4.5 ** 2)
        assert cfg.check_owab is False
        assert cfg.owab_max == 50.0
        assert cfg.check_resolution is False
        assert cfg.resolution_max == 3.5
        assert cfg.use_rdiff is False
        assert cfg.rdiff_max == 0.05
        assert cfg.use_dpi is False
        assert cfg.dpi_max == 0.42
        assert cfg.pdb_redo is False

    def test_distance_is_squared_into_inner_distance(self):
        cfg = rsr_core.AnalysisConfig(distance=2.0)
        assert cfg.inner_distance == 4.0

    def test_custom_values_are_stored_verbatim(self):
        cfg = rsr_core.AnalysisConfig(
            rsr_upper=0.5, rsr_lower=0.3, rscc_min=0.8, rfree_max=0.3,
            occupancy_min=0.9, tolerance=1, check_owab=True, owab_max=40,
            check_resolution=True, resolution_max=2.5, use_rdiff=True,
            rdiff_max=0.02, use_dpi=True, dpi_max=0.3, pdb_redo=True,
        )
        assert cfg.rsr_upper == 0.5
        assert cfg.check_owab is True
        assert cfg.owab_max == 40
        assert cfg.check_resolution is True
        assert cfg.use_rdiff is True
        assert cfg.use_dpi is True
        assert cfg.pdb_redo is True


# ---------------------------------------------------------------------------
# _average_occ
# ---------------------------------------------------------------------------

class FakeAtom:
    """Minimal stand-in exposing only `.occupancy`, for testing _average_occ
    without depending on full PdbAtom construction."""
    def __init__(self, occupancy):
        self.occupancy = occupancy


class TestAverageOcc:
    def test_average_of_multiple_atoms(self):
        atoms = [FakeAtom(1.0), FakeAtom(0.5), FakeAtom(0.5)]
        assert rsr_core._average_occ(atoms) == pytest.approx(2.0 / 3)

    def test_single_atom(self):
        assert rsr_core._average_occ([FakeAtom(0.75)]) == 0.75

    def test_empty_list_raises_zero_division(self):
        with pytest.raises(ZeroDivisionError):
            rsr_core._average_occ([])


# ---------------------------------------------------------------------------
# _dpi
# ---------------------------------------------------------------------------

class TestDpi:
    def test_cubic_cell_matches_manual_formula(self):
        # Orthogonal cell (all angles 90deg) -> V = a*b*c exactly.
        a, b, c = 10.0, 10.0, 10.0
        natoms, reflections, rfree = 100, 5000, 0.2
        result = rsr_core._dpi(a, b, c, 90, 90, 90, natoms, reflections, rfree)
        V = a * b * c
        expected = 1.28 * (natoms ** 0.5) * (V ** (1 / 3)) * (reflections ** (-5 / 6)) * rfree
        assert result == pytest.approx(expected)

    def test_zero_reflections_returns_nan(self):
        assert math.isnan(rsr_core._dpi(10, 10, 10, 90, 90, 90, 100, 0, 0.2))

    def test_degenerate_cell_returns_nan(self):
        # angles that make the volume term negative under sqrt are clamped
        # to 0 by max(0.0, ...), which yields V == 0 -> nan (division by
        # zero guarded by `V <= 0`).
        result = rsr_core._dpi(10, 10, 10, 0, 0, 0, 100, 1000, 0.2)
        assert math.isnan(result)

    def test_higher_rfree_increases_dpi(self):
        low = rsr_core._dpi(10, 10, 10, 90, 90, 90, 100, 5000, 0.1)
        high = rsr_core._dpi(10, 10, 10, 90, 90, 90, 100, 5000, 0.3)
        assert high > low


# ---------------------------------------------------------------------------
# _fmt_reskey (thin wrapper around core.pdb_atom.format_reskey)
# ---------------------------------------------------------------------------

class TestFmtReskey:
    def test_matches_shared_format_reskey(self):
        assert rsr_core._fmt_reskey("REA", "A", "200") == format_reskey("REA", "A", "200")

    def test_no_icode_argument_available(self):
        # rsr_core._fmt_reskey never forwards an icode, unlike
        # core.pdb_atom.format_reskey directly.
        assert rsr_core._fmt_reskey("REA", "A", 200) == "REA A  200"


# ---------------------------------------------------------------------------
# _bbox / residues_bbox / _atoms_for_residues / residue_atom_centers
# ---------------------------------------------------------------------------

class FakeXyzAtom:
    def __init__(self, xyz):
        self.xyz = xyz


class TestBoundingBoxHelpers:
    def test_bbox_of_single_atom_applies_padding(self):
        atoms = [FakeXyzAtom((1.0, 2.0, 3.0))]
        box = rsr_core._bbox(atoms, padding=2.0)
        assert box == {
            "min": [-1.0, 0.0, 1.0],
            "max": [3.0, 4.0, 5.0],
        }

    def test_bbox_of_multiple_atoms_covers_extremes(self):
        atoms = [FakeXyzAtom((0, 0, 0)), FakeXyzAtom((10, -5, 2))]
        box = rsr_core._bbox(atoms, padding=0.0)
        assert box == {"min": [0, -5, 0], "max": [10, 0, 2]}

    def test_bbox_of_empty_atoms_is_none(self):
        assert rsr_core._bbox([]) is None

    def test_default_padding_constant(self):
        assert rsr_core.DEFAULT_BOX_PADDING == 2.1

    def test_atoms_for_residues_merges_both_dicts(self):
        res_atom_dict = {"A": {FakeXyzAtom((0, 0, 0))}}
        ligand_res_atom_dict = {"B": {FakeXyzAtom((1, 1, 1))}}
        atoms = rsr_core._atoms_for_residues(
            ["A", "B", "C"], res_atom_dict, ligand_res_atom_dict)
        assert len(atoms) == 2

    def test_residues_bbox_end_to_end(self):
        res_atom_dict = {"A": {FakeXyzAtom((0, 0, 0))}}
        ligand_res_atom_dict = {"B": {FakeXyzAtom((4, 0, 0))}}
        box = rsr_core.residues_bbox(["A", "B"], res_atom_dict, ligand_res_atom_dict, padding=1.0)
        assert box["min"] == [-1.0, -1.0, -1.0]
        assert box["max"] == [5.0, 1.0, 1.0]

    def test_residues_bbox_returns_none_for_unknown_residues(self):
        assert rsr_core.residues_bbox(["ZZZ"], {}, {}) is None

    def test_residue_atom_centers_flattens_all_atoms(self):
        res_atom_dict = {"A": {FakeXyzAtom((0, 0, 0))}}
        ligand_res_atom_dict = {"B": {FakeXyzAtom((1, 1, 1)), FakeXyzAtom((2, 2, 2))}}
        centers = rsr_core.residue_atom_centers(["A", "B"], res_atom_dict, ligand_res_atom_dict)
        assert len(centers) == 3
        residues_seen = {c["residue"] for c in centers}
        assert residues_seen == {"A", "B"}
        for c in centers:
            assert "center" in c and isinstance(c["center"], list)

    def test_residue_atom_centers_unknown_residue_contributes_nothing(self):
        centers = rsr_core.residue_atom_centers(["ZZZ"], {}, {})
        assert centers == []

    def test_residue_atom_centers_prefers_res_atom_dict_over_ligand_dict(self):
        # `res_atom_dict.get(res) or ligand_res_atom_dict.get(res)` — if a
        # residue key exists (non-empty) in both, res_atom_dict wins.
        res_atom_dict = {"A": {FakeXyzAtom((9, 9, 9))}}
        ligand_res_atom_dict = {"A": {FakeXyzAtom((0, 0, 0))}}
        centers = rsr_core.residue_atom_centers(["A"], res_atom_dict, ligand_res_atom_dict)
        assert centers == [{"residue": "A", "center": [9, 9, 9]}]


# ---------------------------------------------------------------------------
# parse_mmcif_file
# ---------------------------------------------------------------------------

class TestParseMmcifFile:
    def test_protein_ligand_water_metal_classification(self):
        protein_res = make_gemmi_residue(
            "ALA", 1, gemmi.EntityType.Polymer,
            [make_gemmi_atom("CA", 1, 2, 3, serial=1),
             make_gemmi_atom("N", 0, 1, 2, serial=2)])
        ligand_res = make_gemmi_residue(
            "REA", 200, gemmi.EntityType.NonPolymer,
            [make_gemmi_atom("C1", 5, 5, 5, serial=3)])
        water_res = make_gemmi_residue(
            "HOH", 301, gemmi.EntityType.Water,
            [make_gemmi_atom("O", 9, 9, 9, serial=4)])
        # ZN is present in core.cofactors.metals
        metal_res = make_gemmi_residue(
            "ZN", 401, gemmi.EntityType.NonPolymer,
            [make_gemmi_atom("ZN", 8, 8, 8, element="Zn", serial=5)])

        structure = make_gemmi_structure(
            {"A": [protein_res, ligand_res, water_res, metal_res]})

        with patch("core.rsr_core.gemmi.read_structure", return_value=structure):
            result = rsr_core.parse_mmcif_file("/fake/path.cif", "1cbs", inner_distance=20.25)

        assert len(result) == 6
        natoms, res_atom_dict, ligand_res_atom_dict, notligands, links, res_names = result

        protein_key = format_reskey("ALA", "A", 1)
        ligand_key = format_reskey("REA", "A", 200)
        metal_key = format_reskey("ZN", "A", 401)

        assert protein_key in res_atom_dict
        assert len(res_atom_dict[protein_key]) == 2
        assert ligand_key in ligand_res_atom_dict
        assert metal_key in res_atom_dict  # metals go to res_atom_dict, not ligand dict
        assert notligands[metal_key] == "Blacklisted ligand"
        # Water is skipped from both dicts entirely.
        water_key = format_reskey("HOH", "A", 301)
        assert water_key not in res_atom_dict
        assert water_key not in ligand_res_atom_dict
        # natoms accumulates occupancy for every atom, water included
        # (the water `continue` happens after the natoms increment).
        assert natoms == pytest.approx(5.0)  # 5 atoms, occupancy 1.0 each
        assert links == []
        # /fake/path.cif doesn't exist on disk, so the separate
        # _chem_comp-name read (gemmi.cif.read, not mocked here) fails and
        # is swallowed -> no names available, not a crash.
        assert res_names == {}

    def test_zero_inner_distance_skips_protein_atoms(self):
        protein_res = make_gemmi_residue(
            "ALA", 1, gemmi.EntityType.Polymer, [make_gemmi_atom("CA", 1, 2, 3)])
        structure = make_gemmi_structure({"A": [protein_res]})

        with patch("core.rsr_core.gemmi.read_structure", return_value=structure):
            _, res_atom_dict, _, _, _, _ = rsr_core.parse_mmcif_file(
                "/fake/path.cif", "1cbs", inner_distance=0)

        assert res_atom_dict == {}

    def test_nonzero_inner_distance_keeps_protein_atoms(self):
        protein_res = make_gemmi_residue(
            "ALA", 1, gemmi.EntityType.Polymer, [make_gemmi_atom("CA", 1, 2, 3)])
        structure = make_gemmi_structure({"A": [protein_res]})

        with patch("core.rsr_core.gemmi.read_structure", return_value=structure):
            _, res_atom_dict, _, _, _, _ = rsr_core.parse_mmcif_file(
                "/fake/path.cif", "1cbs", inner_distance=20.25)

        assert len(res_atom_dict) == 1

    def test_non_blacklisted_ligand_goes_to_ligand_dict_only(self):
        ligand_res = make_gemmi_residue(
            "REA", 200, gemmi.EntityType.NonPolymer, [make_gemmi_atom("C1", 5, 5, 5)])
        structure = make_gemmi_structure({"A": [ligand_res]})

        with patch("core.rsr_core.gemmi.read_structure", return_value=structure):
            _, res_atom_dict, ligand_res_atom_dict, notligands, _, _ = rsr_core.parse_mmcif_file(
                "/fake/path.cif", "1cbs", inner_distance=20.25)

        key = format_reskey("REA", "A", 200)
        assert key in ligand_res_atom_dict
        assert key not in res_atom_dict
        assert key not in notligands

    def test_blacklisted_nonmetal_ligand_goes_to_res_atom_dict(self):
        # SO4 is in core.cofactors.ligand_blacklist but not in metals.
        so4_res = make_gemmi_residue(
            "SO4", 50, gemmi.EntityType.NonPolymer, [make_gemmi_atom("S", 1, 1, 1)])
        structure = make_gemmi_structure({"A": [so4_res]})

        with patch("core.rsr_core.gemmi.read_structure", return_value=structure):
            _, res_atom_dict, ligand_res_atom_dict, notligands, _, _ = rsr_core.parse_mmcif_file(
                "/fake/path.cif", "1cbs", inner_distance=20.25)

        key = format_reskey("SO4", "A", 50)
        assert key in res_atom_dict
        assert key not in ligand_res_atom_dict
        assert notligands[key] == "Blacklisted ligand"

    def test_covale_and_disulf_and_metalc_connections_are_kept(self):
        structure = make_gemmi_structure({"A": []})
        structure.connections.append(make_connection(
            gemmi.ConnectionType.Covale, "A", "ALA", 1, "CA", "A", "REA", 200, "C1"))
        structure.connections.append(make_connection(
            gemmi.ConnectionType.Disulf, "A", "CYS", 5, "SG", "A", "CYS", 20, "SG"))
        structure.connections.append(make_connection(
            gemmi.ConnectionType.MetalC, "A", "HIS", 8, "NE2", "A", "ZN", 401, "ZN"))

        with patch("core.rsr_core.gemmi.read_structure", return_value=structure):
            _, _, _, _, links, _ = rsr_core.parse_mmcif_file(
                "/fake/path.cif", "1cbs", inner_distance=20.25)

        assert len(links) == 3
        for res1, res2, blen in links:
            assert blen == 1714.0

    def test_hydrog_and_unknown_connections_are_dropped(self):
        structure = make_gemmi_structure({"A": []})
        structure.connections.append(make_connection(
            gemmi.ConnectionType.Hydrog, "A", "SER", 1, "OG", "A", "ASP", 2, "OD1"))
        structure.connections.append(make_connection(
            gemmi.ConnectionType.Unknown, "A", "SER", 1, "OG", "A", "ASP", 2, "OD1"))

        with patch("core.rsr_core.gemmi.read_structure", return_value=structure):
            _, _, _, _, links, _ = rsr_core.parse_mmcif_file(
                "/fake/path.cif", "1cbs", inner_distance=20.25)

        assert links == []

    def test_link_residue_keys_use_format_reskey(self):
        structure = make_gemmi_structure({"A": []})
        structure.connections.append(make_connection(
            gemmi.ConnectionType.Covale, "A", "ALA", 1, "CA", "B", "REA", 200, "C1"))

        with patch("core.rsr_core.gemmi.read_structure", return_value=structure):
            _, _, _, _, links, _ = rsr_core.parse_mmcif_file(
                "/fake/path.cif", "1cbs", inner_distance=20.25)

        res1, res2, _ = links[0]
        assert res1 == format_reskey("ALA", "A", 1)
        assert res2 == format_reskey("REA", "B", 200)

    def test_read_failure_returns_single_element_error_tuple(self):
        with patch("core.rsr_core.gemmi.read_structure", side_effect=RuntimeError("bad file")):
            result = rsr_core.parse_mmcif_file("/fake/path.cif", "1cbs", inner_distance=20.25)

        assert len(result) == 1
        assert "Could not parse mmCIF file" in result[0]
        assert "bad file" in result[0]


# ---------------------------------------------------------------------------
# classificate_residue
# ---------------------------------------------------------------------------

def base_residue_dict(**overrides):
    d = {"RSR": 0.10, "RSCC": 0.95, "occupancy": 1.0}
    d.update(overrides)
    return d


def new_sets():
    return set(), set(), set()  # good, dubious, bad


class TestClassificateResidue:
    def test_perfect_residue_is_good(self):
        cfg = rsr_core.AnalysisConfig()
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score == 0
        assert reason is None
        assert "R1" in good

    def test_missing_residue_dict_is_bad_with_reason(self):
        cfg = rsr_core.AnalysisConfig()
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", None, {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score == 1000
        assert reason == "No data for R1"
        assert "R1" in bad

    def test_empty_residue_dict_is_also_treated_as_missing(self):
        cfg = rsr_core.AnalysisConfig()
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", {}, {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score == 1000
        assert reason == "No data for R1"

    def test_moderate_rsr_gives_dubious(self):
        cfg = rsr_core.AnalysisConfig()  # rsr_lower=0.24, rsr_upper=0.4, tolerance=2
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(RSR=0.30), {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score == 1
        assert "R1" in dubious

    def test_high_rsr_and_low_rscc_gives_bad(self):
        cfg = rsr_core.AnalysisConfig()
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(RSR=0.5, RSCC=0.5), {"rFree": 0.2}, good, dubious, bad, cfg)
        # RSCC fail (+1) + RSR>upper (+2) = 3 > tolerance(2)
        assert score == 3
        assert "R1" in bad

    def test_occupancy_above_one_is_bad_with_specific_reason(self):
        cfg = rsr_core.AnalysisConfig()
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(occupancy=1.5), {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score >= 1000
        assert reason == "Occupancy above 1"
        assert "R1" in bad

    def test_low_occupancy_adds_one_point(self):
        cfg = rsr_core.AnalysisConfig(occupancy_min=1.0)
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(occupancy=0.8), {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score == 1
        assert "R1" in dubious

    def test_rfree_above_max_adds_one_point(self):
        cfg = rsr_core.AnalysisConfig(rfree_max=0.25)
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": 0.3}, good, dubious, bad, cfg)
        assert score == 1

    def test_negative_rfree_is_bad_with_specific_reason(self):
        cfg = rsr_core.AnalysisConfig()
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": -1}, good, dubious, bad, cfg)
        assert reason == "No rFree data for R1"
        assert "R1" in bad

    def test_empty_struc_dict_with_no_advanced_checks_adds_nothing(self):
        cfg = rsr_core.AnalysisConfig()
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {}, good, dubious, bad, cfg)
        assert score == 0
        assert reason is None

    def test_empty_struc_dict_with_advanced_check_enabled_is_bad(self):
        cfg = rsr_core.AnalysisConfig(use_dpi=True)
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {}, good, dubious, bad, cfg)
        assert score >= 1000
        assert reason == "No structural data for R1"
        assert "R1" in bad

    def test_check_owab_out_of_range_adds_one_point(self):
        cfg = rsr_core.AnalysisConfig(check_owab=True, owab_max=50)
        good, dubious, bad = new_sets()
        score, _ = rsr_core.classificate_residue(
            "R1", base_residue_dict(OWAB=60), {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score == 1

    def test_check_owab_in_range_adds_nothing(self):
        cfg = rsr_core.AnalysisConfig(check_owab=True, owab_max=50)
        good, dubious, bad = new_sets()
        score, _ = rsr_core.classificate_residue(
            "R1", base_residue_dict(OWAB=25), {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score == 0

    def test_owab_of_exactly_one_is_out_of_range(self):
        # Condition is `not (1 < owab < max)`, so owab == 1 fails (not <1<1).
        cfg = rsr_core.AnalysisConfig(check_owab=True, owab_max=50)
        good, dubious, bad = new_sets()
        score, _ = rsr_core.classificate_residue(
            "R1", base_residue_dict(OWAB=1), {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score == 1

    def test_check_resolution_above_max_adds_one_point(self):
        cfg = rsr_core.AnalysisConfig(check_resolution=True, resolution_max=3.0)
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": 0.2, "Resolution": 3.5}, good, dubious, bad, cfg)
        assert score == 1
        assert reason is None

    def test_check_resolution_missing_value_uses_sentinel_and_is_bad(self):
        cfg = rsr_core.AnalysisConfig(check_resolution=True, resolution_max=3.0)
        good, dubious, bad = new_sets()
        # struc_dict has no "Resolution" key -> .get default 10 -> triggers
        # both the ">max" (+1) and the "==10 missing sentinel" (+1000) branches.
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score >= 1000
        assert reason == "No resolution data for R1"

    def test_use_rdiff_nan_is_bad_with_specific_reason(self):
        cfg = rsr_core.AnalysisConfig(use_rdiff=True)
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": 0.2, "Rdiff": float("nan")},
            good, dubious, bad, cfg)
        assert score >= 1000
        assert reason == "No reliable rFree/rWork data for R1"

    def test_use_rdiff_missing_key_defaults_to_nan_and_is_bad(self):
        cfg = rsr_core.AnalysisConfig(use_rdiff=True)
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score >= 1000
        assert reason == "No reliable rFree/rWork data for R1"

    def test_use_rdiff_above_max_adds_one_point(self):
        cfg = rsr_core.AnalysisConfig(use_rdiff=True, rdiff_max=0.05)
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": 0.2, "Rdiff": 0.10}, good, dubious, bad, cfg)
        assert score == 1
        assert reason is None

    def test_use_dpi_nan_is_bad(self):
        cfg = rsr_core.AnalysisConfig(use_dpi=True)
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": 0.2, "DPI": float("nan")},
            good, dubious, bad, cfg)
        assert score >= 1000
        assert reason == "No reliable structural data for R1"

    def test_use_dpi_missing_key_defaults_to_minus_one_and_is_bad(self):
        cfg = rsr_core.AnalysisConfig(use_dpi=True)
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": 0.2}, good, dubious, bad, cfg)
        assert score >= 1000
        assert reason == "No reliable structural data for R1"

    def test_use_dpi_above_max_adds_one_point(self):
        cfg = rsr_core.AnalysisConfig(use_dpi=True, dpi_max=0.3)
        good, dubious, bad = new_sets()
        score, reason = rsr_core.classificate_residue(
            "R1", base_residue_dict(), {"rFree": 0.2, "DPI": 0.35}, good, dubious, bad, cfg)
        assert score == 1
        assert reason is None

    def test_score_exactly_at_tolerance_is_dubious_not_bad(self):
        cfg = rsr_core.AnalysisConfig(tolerance=2)
        good, dubious, bad = new_sets()
        # RSR moderate (+1) + low occupancy (+1) = score 2 == tolerance
        score, _ = rsr_core.classificate_residue(
            "R1", base_residue_dict(RSR=0.30, occupancy=0.5), {"rFree": 0.2},
            good, dubious, bad, cfg)
        assert score == 2
        assert "R1" in dubious
        assert "R1" not in bad


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------

class TestValidate:
    def test_all_residues_good(self):
        good, bad, dubious = {"A", "B"}, set(), set()
        assert rsr_core.validate({"A", "B"}, good, bad, dubious) == "Good"

    def test_subset_of_good_is_good(self):
        good, bad, dubious = {"A", "B", "C"}, set(), set()
        assert rsr_core.validate({"A"}, good, bad, dubious) == "Good"

    def test_any_bad_residue_makes_it_bad(self):
        good, bad, dubious = {"A"}, {"B"}, set()
        assert rsr_core.validate({"A", "B"}, good, bad, dubious) == "Bad"

    def test_bad_takes_priority_over_dubious(self):
        good, bad, dubious = set(), {"A"}, {"B"}
        assert rsr_core.validate({"A", "B"}, good, bad, dubious) == "Bad"

    def test_any_dubious_residue_without_bad_is_dubious(self):
        good, bad, dubious = {"A"}, set(), {"B"}
        assert rsr_core.validate({"A", "B"}, good, bad, dubious) == "Dubious"

    def test_unclassified_residue_defaults_to_dubious(self):
        # Residue not present in good/bad/dubious at all.
        good, bad, dubious = set(), set(), set()
        assert rsr_core.validate({"Z"}, good, bad, dubious) == "Dubious"

    def test_empty_residues_set_is_trivially_good(self):
        # empty set <= any set is True, so this hits the "Good" branch.
        good, bad, dubious = set(), set(), set()
        assert rsr_core.validate(set(), good, bad, dubious) == "Good"


# ---------------------------------------------------------------------------
# group_ligands
# ---------------------------------------------------------------------------

class TestGroupLigands:
    def test_single_unlinked_ligand(self):
        result = rsr_core.group_ligands({"A"}, [])
        assert result == [{"A"}]

    def test_multiple_unlinked_ligands_stay_separate(self):
        result = rsr_core.group_ligands({"A", "B"}, [])
        assert {"A"} in result and {"B"} in result
        assert len(result) == 2

    def test_two_linked_residues_are_merged(self):
        result = rsr_core.group_ligands({"A", "B"}, [("A", "B", 1.5)])
        assert result == [{"A", "B"}]

    def test_chain_of_links_merges_transitively(self):
        result = rsr_core.group_ligands(
            {"A", "B", "C"}, [("A", "B", 1.5), ("B", "C", 1.5)])
        assert result == [{"A", "B", "C"}]

    def test_links_to_residues_outside_ligand_set_are_ignored(self):
        result = rsr_core.group_ligands({"A"}, [("A", "X", 1.5)])
        assert result == [{"A"}]

    def test_unrelated_groups_stay_separate_and_linked_group_merges(self):
        result = rsr_core.group_ligands(
            {"A", "B", "C", "D"}, [("A", "B", 1.5), ("C", "D", 1.5)])
        result_sets = sorted(result, key=lambda s: sorted(s))
        assert result_sets == [{"A", "B"}, {"C", "D"}]

    def test_two_pairs_bridged_by_third_link_merge_into_one_group(self):
        # (A,B) and (C,D) form independent pairs, but (B,C) bridges them —
        # the trailing "merge overlapping sets" pass must combine all four.
        result = rsr_core.group_ligands(
            {"A", "B", "C", "D"},
            [("A", "B", 1.5), ("C", "D", 1.5), ("B", "C", 1.5)])
        assert len(result) == 1
        assert result[0] == {"A", "B", "C", "D"}


# ---------------------------------------------------------------------------
# get_binding_site
# ---------------------------------------------------------------------------

from core.pdb_atom import PdbAtom  # noqa: E402  (local import, keeps top imports tidy)


def make_real_atom(comp_id, asym_id, seq_id, x, y, z, occ=1.0, name="X1"):
    return PdbAtom({
        "auth_atom_id": name,
        "auth_comp_id": comp_id,
        "auth_asym_id": asym_id,
        "auth_seq_id": seq_id,
        "pdbx_PDB_ins_code": "",
        "Cartn_x": x, "Cartn_y": y, "Cartn_z": z,
        "occupancy": occ,
        "label_alt_id": ".",
    })


class TestGetBindingSite:
    def test_happy_path_single_residue_binding_site(self):
        cfg = rsr_core.AnalysisConfig()  # inner_distance = 4.5**2 = 20.25
        ligand_key = format_reskey("REA", "A", 200)
        bs_key = format_reskey("ALA", "A", 1)

        ligand_atom = make_real_atom("REA", "A", 200, 0, 0, 0)
        bs_atom = make_real_atom("ALA", "A", 1, 3, 0, 0)  # dist^2 = 9

        ligand_res_atom_dict = {ligand_key: {ligand_atom}}
        res_atom_dict = {bs_key: {bs_atom}}
        edd_dict = {
            ligand_key: {"RSR": 0.1, "RSCC": 0.95, "occupancy": 1.0},
            bs_key: {"RSR": 0.1, "RSCC": 0.95, "occupancy": 1.0},
        }
        struc_dict = {"rFree": 0.2}
        notligands = {}
        good_rsr, dubious_rsr, bad_rsr = set(), set(), set()
        ligand = {ligand_key}
        ligands = [ligand]

        # Mirrors the pipeline: ligand residues are classified before
        # get_binding_site() is called.
        rsr_core.classificate_residue(
            ligand_key, edd_dict[ligand_key], struc_dict,
            good_rsr, dubious_rsr, bad_rsr, cfg)

        result = rsr_core.get_binding_site(
            ligand, 0, good_rsr, bad_rsr, dubious_rsr, "1cbs",
            res_atom_dict, ligands, ligand_res_atom_dict,
            edd_dict, struc_dict, notligands, cfg)

        (ligandresidues, binding_site, rte, ligandgood, bsgood,
         bad_occupancy, lig_score, bs_score) = result

        assert binding_site == {bs_key}
        assert bad_occupancy == []
        assert ligandgood == "Good"
        assert bsgood == "Good"
        assert bs_score == 0
        assert rte == set()

    def test_distant_protein_residue_excluded(self):
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("REA", "A", 200)
        far_key = format_reskey("GLY", "A", 2)

        ligand_atom = make_real_atom("REA", "A", 200, 0, 0, 0)
        far_atom = make_real_atom("GLY", "A", 2, 100, 100, 100)

        ligand_res_atom_dict = {ligand_key: {ligand_atom}}
        res_atom_dict = {far_key: {far_atom}}
        edd_dict = {ligand_key: {"RSR": 0.1, "RSCC": 0.95, "occupancy": 1.0}}
        struc_dict = {"rFree": 0.2}
        good_rsr, dubious_rsr, bad_rsr = set(), set(), set()
        ligand = {ligand_key}

        result = rsr_core.get_binding_site(
            ligand, 0, good_rsr, bad_rsr, dubious_rsr, "1cbs",
            res_atom_dict, [ligand], ligand_res_atom_dict,
            edd_dict, struc_dict, {}, cfg)

        binding_site = result[1]
        assert binding_site == set()

    def test_ligand_already_in_notligands_short_circuits(self):
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("REA", "A", 200)
        notligands = {ligand_key: "Some earlier rejection reason"}

        result = rsr_core.get_binding_site(
            {ligand_key}, 0, set(), set(), set(), "1cbs",
            {}, [{ligand_key}], {ligand_key: set()},
            {}, {}, notligands, cfg)

        assert result == ["Some earlier rejection reason"]

    def test_close_contact_to_metal_ligand_rejects_as_covalent(self):
        cfg = rsr_core.AnalysisConfig()
        # "ZN" is in core.cofactors.metals.
        ligand_key = format_reskey("ZN", "A", 401)
        bs_key = format_reskey("HIS", "A", 8)

        ligand_atom = make_real_atom("ZN", "A", 401, 0, 0, 0)
        bs_atom = make_real_atom("HIS", "A", 8, 0, 0, 0)  # dist == 0 < 2.1

        result = rsr_core.get_binding_site(
            {ligand_key}, 0, set(), set(), set(), "1cbs",
            {bs_key: {bs_atom}}, [{ligand_key}], {ligand_key: {ligand_atom}},
            {}, {}, {}, cfg)

        assert result == ["Covalently bound to the sequence"]

    def test_close_contact_to_blacklisted_ligand_rejects_with_blacklist_reason(self):
        cfg = rsr_core.AnalysisConfig()
        # "SO4" is in core.cofactors.ligand_blacklist but not in metals.
        ligand_key = format_reskey("SO4", "A", 50)
        bs_key = format_reskey("HIS", "A", 8)

        ligand_atom = make_real_atom("SO4", "A", 50, 0, 0, 0)
        bs_atom = make_real_atom("HIS", "A", 8, 0, 0, 0)

        result = rsr_core.get_binding_site(
            {ligand_key}, 0, set(), set(), set(), "1cbs",
            {bs_key: {bs_atom}}, [{ligand_key}], {ligand_key: {ligand_atom}},
            {}, {}, {}, cfg)

        assert result == ["Covalently bound to a blacklisted ligand"]

    def test_cross_ligand_contact_pulls_in_other_ligand_residue(self):
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("REA", "A", 200)
        other_key = format_reskey("EDO", "A", 500)

        ligand_atom = make_real_atom("REA", "A", 200, 0, 0, 0)
        other_atom = make_real_atom("EDO", "A", 500, 2, 0, 0)  # dist^2 = 4

        ligand_res_atom_dict = {
            ligand_key: {ligand_atom},
            other_key: {other_atom},
        }
        ligand_group = {ligand_key}
        other_group = {other_key}

        result = rsr_core.get_binding_site(
            ligand_group, 0, set(), set(), set(), "1cbs",
            {}, [ligand_group, other_group], ligand_res_atom_dict,
            {}, {}, {}, cfg)

        binding_site = result[1]
        assert other_key in binding_site

    def test_low_occupancy_ligand_residue_flagged(self):
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("REA", "A", 200)
        edd_dict = {ligand_key: {"occupancy": 0.5}}

        result = rsr_core.get_binding_site(
            {ligand_key}, 0, set(), set(), set(), "1cbs",
            {}, [{ligand_key}], {ligand_key: set()},
            edd_dict, {}, {}, cfg)

        bad_occupancy = result[5]
        assert ligand_key in bad_occupancy

    def test_binding_site_residue_missing_validation_data_flagged_bad_occupancy(self):
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("REA", "A", 200)
        bs_key = format_reskey("ALA", "A", 1)
        ligand_atom = make_real_atom("REA", "A", 200, 0, 0, 0)
        bs_atom = make_real_atom("ALA", "A", 1, 1, 0, 0)

        result = rsr_core.get_binding_site(
            {ligand_key}, 0, set(), set(), set(), "1cbs",
            {bs_key: {bs_atom}}, [{ligand_key}], {ligand_key: {ligand_atom}},
            {},  # no edd_dict entry at all for bs_key
            {}, {}, cfg)

        bad_occupancy = result[5]
        assert bs_key in bad_occupancy

    def test_rte_ignores_binding_site_quality_when_ligand_is_good(self):
        # Regression test documenting current operator-precedence behaviour:
        # `rte = inner_binding_site | ligand - good_rsr` evaluates as
        # `inner_binding_site | (ligand - good_rsr)` (`-` binds tighter than
        # `|`), so whenever the ligand itself is fully "good", rte is
        # always exactly `inner_binding_site`, even for a "Bad" binding site.
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("REA", "A", 200)
        bs_key = format_reskey("ALA", "A", 1)
        ligand_atom = make_real_atom("REA", "A", 200, 0, 0, 0)
        bs_atom = make_real_atom("ALA", "A", 1, 1, 0, 0)

        good_rsr = {ligand_key}  # ligand pre-classified as good
        edd_dict = {bs_key: {"RSR": 0.9, "RSCC": 0.1, "occupancy": 1.0}}  # very bad BS

        result = rsr_core.get_binding_site(
            {ligand_key}, 0, good_rsr, set(), set(), "1cbs",
            {bs_key: {bs_atom}}, [{ligand_key}], {ligand_key: {ligand_atom}},
            edd_dict, {"rFree": 0.2}, {}, cfg)

        _, binding_site, rte, _, bsgood, _, _, _ = result
        assert bsgood == "Bad"
        assert rte == binding_site  # NOT empty, despite bsgood == "Bad"


# ---------------------------------------------------------------------------
# parse_binding_site
# ---------------------------------------------------------------------------

class TestParseBindingSite:
    def test_pdb_redo_missing_stats_returns_error(self):
        cfg = rsr_core.AnalysisConfig(pdb_redo=True)
        with patch("core.rsr_core.pdb_redo_utils.get_pdbredo_data", return_value=None):
            result = rsr_core.parse_binding_site("1cbs", cfg)
        assert result == {"pdbid": "1cbs", "error": "Not in PDB_REDO"}

    def test_pdb_redo_missing_ed_data_returns_error(self):
        cfg = rsr_core.AnalysisConfig(pdb_redo=True)
        with patch("core.rsr_core.pdb_redo_utils.get_pdbredo_data",
                   return_value={"rFree": 0.2, "rWork": 0.18}), \
             patch("core.rsr_core.pdb_redo_utils.get_ED_data", return_value=None):
            result = rsr_core.parse_binding_site("1cbs", cfg)
        assert result == {"pdbid": "1cbs", "error": "No PDB-REDO ED data available"}

    def test_default_path_missing_validation_data_returns_error(self):
        cfg = rsr_core.AnalysisConfig()
        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value={}), \
             patch("core.rsr_core.eds_utils.get_EDS", return_value=({}, {})):
            result = rsr_core.parse_binding_site("1cbs", cfg)
        assert result == {
            "pdbid": "1cbs",
            "error": "No EDM/validation data available (may not be an X-ray entry)",
        }

    def test_pdb_file_download_failure_returns_error(self):
        cfg = rsr_core.AnalysisConfig()
        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value={}), \
             patch("core.rsr_core.eds_utils.get_EDS",
                   return_value=({"1cbs": True}, {"X": {"RSR": 0.1}})), \
             patch("core.rsr_core.pdb_utils.get_pdb_file", return_value=""):
            result = rsr_core.parse_binding_site("1cbs", cfg)
        assert result == {"pdbid": "1cbs", "error": "Unable to load PDBx/mmCIF model"}

    def test_mmcif_parse_failure_propagates_error_message(self):
        cfg = rsr_core.AnalysisConfig()
        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value={}), \
             patch("core.rsr_core.eds_utils.get_EDS",
                   return_value=({"1cbs": True}, {"X": {"RSR": 0.1}})), \
             patch("core.rsr_core.pdb_utils.get_pdb_file", return_value="/tmp/f.cif"), \
             patch("core.rsr_core.parse_mmcif_file",
                   return_value=("boom: could not parse",)):
            result = rsr_core.parse_binding_site("1cbs", cfg)
        assert result == {"pdbid": "1cbs", "error": "boom: could not parse"}

    def test_no_ligands_found_returns_error(self):
        cfg = rsr_core.AnalysisConfig()
        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value={}), \
             patch("core.rsr_core.eds_utils.get_EDS",
                   return_value=({"1cbs": True}, {"X": {"RSR": 0.1, "RSCC": 0.9, "occupancy": 1.0}})), \
             patch("core.rsr_core.pdb_utils.get_pdb_file", return_value="/tmp/f.cif"), \
             patch("core.rsr_core.parse_mmcif_file",
                   return_value=(10.0, {}, {}, {}, [], {})):  # empty ligand_res_atom_dict
            result = rsr_core.parse_binding_site("1cbs", cfg)
        assert result == {"pdbid": "1cbs", "error": "No ligands found"}

    def test_successful_minimal_analysis(self):
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("REA", "A", 200)
        bs_key = format_reskey("ALA", "A", 1)
        ligand_atom = make_real_atom("REA", "A", 200, 0, 0, 0)
        bs_atom = make_real_atom("ALA", "A", 1, 1, 0, 0)

        edd_dict = {
            ligand_key: {"RSR": 0.1, "RSCC": 0.95, "occupancy": 1.0},
            bs_key: {"RSR": 0.1, "RSCC": 0.95, "occupancy": 1.0},
        }
        parsed = (
            10.0,
            {bs_key: {bs_atom}},
            {ligand_key: {ligand_atom}},
            {},
            [],
            {ligand_key: "RETINOIC ACID"},
        )
        rcsb_report = {"1CBS": {"rFree": 0.2, "rWork": 0.18}}

        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value=rcsb_report), \
             patch("core.rsr_core.eds_utils.get_EDS", return_value=({"1cbs": True}, edd_dict)), \
             patch("core.rsr_core.pdb_utils.get_pdb_file", return_value="/tmp/f.cif"), \
             patch("core.rsr_core.parse_mmcif_file", return_value=parsed):
            result = rsr_core.parse_binding_site("1cbs", cfg)

        assert result["pdbid"] == "1cbs"
        assert "error" not in result
        assert len(result["ligands"]) == 1
        lig = result["ligands"][0]
        assert lig["ligand_residues"] == [ligand_key]
        assert lig["ligand_names"] == ["RETINOIC ACID"]
        assert lig["binding_site_residues"] == [bs_key]
        assert lig["ligand_quality"] == "Good"
        assert lig["binding_site_quality"] == "Good"
        assert lig["source"] == "PDB"
        assert result["struc_dict"]["rFree"] == 0.2
        assert result["rejected"] == {}

    def test_struc_dict_nan_values_become_none(self):
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("REA", "A", 200)
        ligand_atom = make_real_atom("REA", "A", 200, 0, 0, 0)
        edd_dict = {ligand_key: {"RSR": 0.1, "RSCC": 0.95, "occupancy": 1.0}}
        parsed = (5.0, {}, {ligand_key: {ligand_atom}}, {}, [], {})

        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value={}), \
             patch("core.rsr_core.eds_utils.get_EDS", return_value=({"1cbs": True}, edd_dict)), \
             patch("core.rsr_core.pdb_utils.get_pdb_file", return_value="/tmp/f.cif"), \
             patch("core.rsr_core.parse_mmcif_file", return_value=parsed):
            result = rsr_core.parse_binding_site("1cbs", cfg)

        # No RCSB stats at all -> rFree/rWork default to nan -> serialised as None.
        assert result["struc_dict"]["rFree"] is None
        assert result["struc_dict"]["rWork"] is None

    def test_occupancy_gap_filled_from_atom_average(self):
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("REA", "A", 200)
        ligand_atom1 = make_real_atom("REA", "A", 200, 0, 0, 0, occ=1.0, name="C1")
        ligand_atom2 = make_real_atom("REA", "A", 200, 0.1, 0, 0, occ=0.5, name="C2")
        # RSR/RSCC given, but no "occupancy" key -> should be filled from
        # the average occupancy of the ligand's own atoms (0.75).
        edd_dict = {ligand_key: {"RSR": 0.1, "RSCC": 0.95}}
        parsed = (5.0, {}, {ligand_key: {ligand_atom1, ligand_atom2}}, {}, [], {})

        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value={}), \
             patch("core.rsr_core.eds_utils.get_EDS", return_value=({"1cbs": True}, edd_dict)), \
             patch("core.rsr_core.pdb_utils.get_pdb_file", return_value="/tmp/f.cif"), \
             patch("core.rsr_core.parse_mmcif_file", return_value=parsed):
            result = rsr_core.parse_binding_site("1cbs", cfg)

        # occupancy 0.75 < occupancy_min(1.0) -> ligand flagged low_occupancy
        assert result["ligands"][0]["low_occupancy"] == [ligand_key]

    def test_residue_with_no_occupancy_and_no_atoms_is_dropped_from_edd(self):
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("REA", "A", 200)
        ligand_atom = make_real_atom("REA", "A", 200, 0, 0, 0)
        ghost_key = format_reskey("XXX", "A", 999)
        # "ghost_key" has no occupancy and no matching atoms anywhere ->
        # dropped from edd_dict entirely by the occupancy-gap-filling step.
        edd_dict = {
            ligand_key: {"RSR": 0.1, "RSCC": 0.95, "occupancy": 1.0},
            ghost_key: {"RSR": 0.2, "RSCC": 0.90},
        }
        parsed = (5.0, {}, {ligand_key: {ligand_atom}}, {}, [], {})

        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value={}), \
             patch("core.rsr_core.eds_utils.get_EDS", return_value=({"1cbs": True}, edd_dict)), \
             patch("core.rsr_core.pdb_utils.get_pdb_file", return_value="/tmp/f.cif"), \
             patch("core.rsr_core.parse_mmcif_file", return_value=parsed) as mock_parse:
            result = rsr_core.parse_binding_site("1cbs", cfg)

        assert result["pdbid"] == "1cbs"
        # No direct way to inspect edd_dict post-hoc, but the call must not
        # raise (KeyError/ZeroDivisionError would surface as an "error" key
        # coming from analyse_pdbids, not from this function directly, so
        # here we simply assert the happy path still completed).
        assert "error" not in result

    def test_use_dpi_populates_dpi_field(self):
        cfg = rsr_core.AnalysisConfig(use_dpi=True)
        ligand_key = format_reskey("REA", "A", 200)
        ligand_atom = make_real_atom("REA", "A", 200, 0, 0, 0)
        edd_dict = {ligand_key: {"RSR": 0.1, "RSCC": 0.95, "occupancy": 1.0}}
        parsed = (100.0, {}, {ligand_key: {ligand_atom}}, {}, [], {})
        rcsb_report = {
            "1CBS": {
                "rFree": 0.2, "rWork": 0.18, "nreflections": 5000,
                "lengthOfUnitCellLatticeA": 40.0, "lengthOfUnitCellLatticeB": 40.0,
                "lengthOfUnitCellLatticeC": 40.0, "unitCellAngleAlpha": 90.0,
                "unitCellAngleBeta": 90.0, "unitCellAngleGamma": 90.0,
            }
        }

        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value=rcsb_report), \
             patch("core.rsr_core.eds_utils.get_EDS", return_value=({"1cbs": True}, edd_dict)), \
             patch("core.rsr_core.pdb_utils.get_pdb_file", return_value="/tmp/f.cif"), \
             patch("core.rsr_core.parse_mmcif_file", return_value=parsed):
            result = rsr_core.parse_binding_site("1cbs", cfg)

        assert result["struc_dict"]["DPI"] is not None
        assert result["struc_dict"]["DPI"] > 0

    def test_use_dpi_zero_reflections_gives_none_dpi(self):
        cfg = rsr_core.AnalysisConfig(use_dpi=True)
        ligand_key = format_reskey("REA", "A", 200)
        ligand_atom = make_real_atom("REA", "A", 200, 0, 0, 0)
        edd_dict = {ligand_key: {"RSR": 0.1, "RSCC": 0.95, "occupancy": 1.0}}
        parsed = (100.0, {}, {ligand_key: {ligand_atom}}, {}, [], {})
        rcsb_report = {"1CBS": {"rFree": 0.2, "rWork": 0.18, "nreflections": 0}}

        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value=rcsb_report), \
             patch("core.rsr_core.eds_utils.get_EDS", return_value=({"1cbs": True}, edd_dict)), \
             patch("core.rsr_core.pdb_utils.get_pdb_file", return_value="/tmp/f.cif"), \
             patch("core.rsr_core.parse_mmcif_file", return_value=parsed):
            result = rsr_core.parse_binding_site("1cbs", cfg)

        assert result["struc_dict"]["DPI"] is None  # NaN serialised to None

    def test_covalently_pruned_ligand_is_moved_to_rejected(self):
        cfg = rsr_core.AnalysisConfig()
        ligand_key = format_reskey("EDO", "A", 500)  # not blacklisted, not metal
        bs_key = format_reskey("ALA", "A", 1)
        ligand_atom = make_real_atom("EDO", "A", 500, 0, 0, 0)
        bs_atom = make_real_atom("ALA", "A", 1, 1, 0, 0)

        # edd_dict is deliberately non-empty (so the earlier "no validation
        # data at all" guard doesn't trigger) but has no entry that matters
        # to the pruning logic itself.
        edd_dict = {"UNUSED X    9": {"RSR": 1.0}}
        # A short (<2.1) covalent-type link between the ligand and a
        # protein residue -> pruning step should move EDO into
        # res_atom_dict/notligands instead of treating it as a ligand.
        links = [(bs_key, ligand_key, 1.5)]
        parsed = (
            5.0,
            {bs_key: {bs_atom}},
            {ligand_key: {ligand_atom}},
            {},
            links,
            {},
        )

        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value={}), \
             patch("core.rsr_core.eds_utils.get_EDS", return_value=({"1cbs": True}, edd_dict)), \
             patch("core.rsr_core.pdb_utils.get_pdb_file", return_value="/tmp/f.cif"), \
             patch("core.rsr_core.parse_mmcif_file", return_value=parsed):
            result = rsr_core.parse_binding_site("1cbs", cfg)

        assert result == {"pdbid": "1cbs", "error": "No ligands found"}

    def test_lowercases_and_uses_pdbid_in_result(self):
        cfg = rsr_core.AnalysisConfig()
        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value={}), \
             patch("core.rsr_core.eds_utils.get_EDS", return_value=({}, {})):
            result = rsr_core.parse_binding_site("1CBS", cfg)
        assert result["pdbid"] == "1cbs"

    def test_default_cfg_is_created_when_none_passed(self):
        with patch("core.rsr_core.pdb_utils.get_custom_report", return_value={}), \
             patch("core.rsr_core.eds_utils.get_EDS", return_value=({}, {})):
            # Should not raise despite cfg=None.
            result = rsr_core.parse_binding_site("1cbs", None)
        assert result["pdbid"] == "1cbs"


# ---------------------------------------------------------------------------
# analyse_pdbids
# ---------------------------------------------------------------------------

class TestAnalysePdbids:
    def test_calls_parse_binding_site_for_each_id_stripped_and_lowercased(self):
        with patch("core.rsr_core.parse_binding_site") as mock_parse:
            mock_parse.side_effect = lambda pdbid, cfg: {"pdbid": pdbid}
            results = rsr_core.analyse_pdbids([" 1CBS ", "3DZU"])

        assert [c.args[0] for c in mock_parse.call_args_list] == ["1cbs", "3dzu"]
        assert results == [{"pdbid": "1cbs"}, {"pdbid": "3dzu"}]

    def test_uses_default_config_when_none_given(self):
        with patch("core.rsr_core.parse_binding_site") as mock_parse:
            mock_parse.return_value = {"pdbid": "1cbs"}
            rsr_core.analyse_pdbids(["1cbs"])
        cfg_arg = mock_parse.call_args[0][1]
        assert isinstance(cfg_arg, rsr_core.AnalysisConfig)

    def test_passes_through_given_config(self):
        cfg = rsr_core.AnalysisConfig(tolerance=5)
        with patch("core.rsr_core.parse_binding_site") as mock_parse:
            mock_parse.return_value = {"pdbid": "1cbs"}
            rsr_core.analyse_pdbids(["1cbs"], cfg)
        assert mock_parse.call_args[0][1] is cfg

    def test_one_failure_does_not_abort_the_batch(self):
        def fake_parse(pdbid, cfg):
            if pdbid == "bad1":
                raise RuntimeError("kaboom")
            return {"pdbid": pdbid, "ligands": []}

        with patch("core.rsr_core.parse_binding_site", side_effect=fake_parse):
            results = rsr_core.analyse_pdbids(["1cbs", "bad1", "3dzu"])

        assert len(results) == 3
        assert results[0]["pdbid"] == "1cbs"
        assert results[1] == {"pdbid": "bad1", "error": "kaboom"}
        assert results[2]["pdbid"] == "3dzu"

    def test_empty_list_returns_empty_results(self):
        assert rsr_core.analyse_pdbids([]) == []

    def test_on_progress_called_once_per_entry_with_final_total(self):
        calls = []

        def on_progress(completed, total):
            calls.append((completed, total))

        with patch("core.rsr_core.parse_binding_site") as mock_parse:
            mock_parse.side_effect = lambda pdbid, cfg: {"pdbid": pdbid}
            rsr_core.analyse_pdbids(["1cbs", "3dzu", "2xyz"], on_progress=on_progress)

        assert len(calls) == 3
        # completed values are exactly 1..N, seen exactly once each,
        # regardless of the order background workers finish in.
        assert sorted(c for c, _ in calls) == [1, 2, 3]
        # total is always the full batch size on every call.
        assert all(t == 3 for _, t in calls)

    def test_on_progress_still_called_when_an_entry_fails(self):
        calls = []

        def fake_parse(pdbid, cfg):
            if pdbid == "bad1":
                raise RuntimeError("kaboom")
            return {"pdbid": pdbid}

        with patch("core.rsr_core.parse_binding_site", side_effect=fake_parse):
            rsr_core.analyse_pdbids(
                ["1cbs", "bad1"], on_progress=lambda completed, total: calls.append((completed, total))
            )

        assert len(calls) == 2
        assert sorted(c for c, _ in calls) == [1, 2]

    def test_on_progress_not_required(self):
        # Default is None; must not raise when omitted.
        with patch("core.rsr_core.parse_binding_site") as mock_parse:
            mock_parse.return_value = {"pdbid": "1cbs"}
            results = rsr_core.analyse_pdbids(["1cbs"])
        assert results == [{"pdbid": "1cbs"}]

    def test_on_progress_not_called_for_empty_input(self):
        calls = []
        rsr_core.analyse_pdbids([], on_progress=lambda c, t: calls.append((c, t)))
        assert calls == []
