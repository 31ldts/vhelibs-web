# -*- coding: utf-8 -*-
"""
Unit tests for core.pdb_atom.

Pure in-memory logic, no external dependencies to mock.
"""
import pytest

from core.pdb_atom import PdbAtom, format_reskey


def make_atom_dict(**overrides):
    base = {
        "auth_atom_id": "CA",
        "auth_comp_id": "ALA",
        "auth_asym_id": "A",
        "auth_seq_id": "42",
        "pdbx_PDB_ins_code": "",
        "Cartn_x": "1.0",
        "Cartn_y": "2.0",
        "Cartn_z": "3.0",
        "occupancy": "1.0",
        "label_alt_id": ".",
    }
    base.update(overrides)
    return base


class TestFormatReskey:
    def test_pads_comp_id_to_three_chars(self):
        assert format_reskey("CA", "A", "1") == " CA A    1"

    def test_no_padding_needed_for_three_char_comp_id(self):
        assert format_reskey("ATP", "A", "1") == "ATP A    1"

    def test_pads_seq_id_to_four_chars(self):
        assert format_reskey("REA", "A", "200") == "REA A  200"

    def test_no_padding_needed_for_four_digit_seq_id(self):
        # This is the exact regression this function's docstring describes:
        # a literal space always separates asym_id and seq_id, even when
        # seq_id already has 4+ digits.
        assert format_reskey("REA", "A", "1000") == "REA A 1000"

    def test_five_digit_seq_id_is_not_truncated(self):
        assert format_reskey("REA", "A", "10000") == "REA A 10000"

    def test_icode_is_appended_verbatim(self):
        assert format_reskey("REA", "A", "200", "B") == "REA A  200B"

    def test_icode_defaults_to_empty_string_when_falsy(self):
        assert format_reskey("REA", "A", "200", "") == "REA A  200"
        assert format_reskey("REA", "A", "200", None) == "REA A  200"

    def test_accepts_integer_seq_id(self):
        assert format_reskey("REA", "A", 200) == "REA A  200"

    def test_single_char_asym_id(self):
        key = format_reskey("HOH", "B", "5")
        assert key.split(" ")[1] == "B"


class TestPdbAtomInit:
    def test_basic_field_extraction(self):
        atom = PdbAtom(make_atom_dict())
        assert atom.name == "CA"
        assert atom.hetid == "ALA"
        assert atom.xyz == (1.0, 2.0, 3.0)
        assert atom.occupancy == 1.0
        assert atom.variant == "."

    def test_residue_uses_format_reskey(self):
        atom = PdbAtom(make_atom_dict(
            auth_comp_id="REA", auth_asym_id="A", auth_seq_id="200"))
        assert atom.residue == format_reskey("REA", "A", "200")

    def test_residue_includes_insertion_code(self):
        atom = PdbAtom(make_atom_dict(pdbx_PDB_ins_code="B"))
        assert atom.residue.endswith("B")

    def test_missing_ins_code_key_defaults_to_empty(self):
        d = make_atom_dict()
        del d["pdbx_PDB_ins_code"]
        atom = PdbAtom(d)
        assert atom.residue == format_reskey("ALA", "A", "42", "")

    def test_missing_required_key_raises_keyerror(self):
        d = make_atom_dict()
        del d["auth_atom_id"]
        with pytest.raises(KeyError):
            PdbAtom(d)

    def test_coordinates_and_occupancy_are_cast_to_float(self):
        atom = PdbAtom(make_atom_dict(
            Cartn_x="1.5", Cartn_y="-2.25", Cartn_z="0", occupancy="0.5"))
        assert atom.xyz == (1.5, -2.25, 0.0)
        assert atom.occupancy == 0.5
        assert isinstance(atom.occupancy, float)

    def test_invalid_numeric_field_raises_valueerror(self):
        with pytest.raises(ValueError):
            PdbAtom(make_atom_dict(occupancy="not-a-number"))

    def test_has_no_dict_due_to_slots(self):
        atom = PdbAtom(make_atom_dict())
        with pytest.raises(AttributeError):
            atom.__dict__


class TestPdbAtomOr:
    def test_squared_distance_between_identical_points_is_zero(self):
        a1 = PdbAtom(make_atom_dict(Cartn_x="1", Cartn_y="1", Cartn_z="1"))
        a2 = PdbAtom(make_atom_dict(Cartn_x="1", Cartn_y="1", Cartn_z="1"))
        assert (a1 | a2) == 0.0

    def test_squared_distance_matches_manual_computation(self):
        a1 = PdbAtom(make_atom_dict(Cartn_x="0", Cartn_y="0", Cartn_z="0"))
        a2 = PdbAtom(make_atom_dict(Cartn_x="3", Cartn_y="4", Cartn_z="0"))
        # squared distance = 3^2 + 4^2 + 0^2 = 25
        assert (a1 | a2) == 25.0

    def test_or_is_symmetric(self):
        a1 = PdbAtom(make_atom_dict(Cartn_x="0", Cartn_y="0", Cartn_z="0"))
        a2 = PdbAtom(make_atom_dict(Cartn_x="1", Cartn_y="2", Cartn_z="2"))
        assert (a1 | a2) == (a2 | a1)


class TestPdbAtomHashAndEq:
    def test_equal_atoms_are_equal_and_share_hash(self):
        a1 = PdbAtom(make_atom_dict())
        a2 = PdbAtom(make_atom_dict())
        assert a1 == a2
        assert hash(a1) == hash(a2)

    def test_different_name_makes_atoms_unequal(self):
        a1 = PdbAtom(make_atom_dict(auth_atom_id="CA"))
        a2 = PdbAtom(make_atom_dict(auth_atom_id="CB"))
        assert a1 != a2

    def test_different_xyz_makes_atoms_unequal(self):
        a1 = PdbAtom(make_atom_dict(Cartn_x="1.0"))
        a2 = PdbAtom(make_atom_dict(Cartn_x="2.0"))
        assert a1 != a2

    def test_equality_ignores_occupancy_and_variant(self):
        # __eq__/__hash__ only look at (residue, name, xyz) — occupancy and
        # variant are not part of the identity, by current implementation.
        a1 = PdbAtom(make_atom_dict(occupancy="1.0", label_alt_id="A"))
        a2 = PdbAtom(make_atom_dict(occupancy="0.5", label_alt_id="B"))
        assert a1 == a2
        assert hash(a1) == hash(a2)

    def test_atoms_are_usable_in_a_set(self):
        a1 = PdbAtom(make_atom_dict())
        a2 = PdbAtom(make_atom_dict())  # duplicate of a1
        a3 = PdbAtom(make_atom_dict(auth_atom_id="CB"))
        s = {a1, a2, a3}
        assert len(s) == 2

    def test_eq_accesses_attributes_directly_no_type_guard(self):
        # __eq__ does not check isinstance(other, PdbAtom); comparing
        # against an object lacking the expected attributes raises
        # AttributeError rather than returning False/NotImplemented.
        atom = PdbAtom(make_atom_dict())
        with pytest.raises(AttributeError):
            atom == object()
