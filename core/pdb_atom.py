# -*- coding: utf-8 -*-
#
#   Copyright 2010-2024 Adrià Cereto Massagué
#   Migrated to web version – pure Python fallback is now the primary implementation.
#


class PdbAtom:
    """Represents an atom from a PDB/mmCIF file."""

    __slots__ = ('name', 'residue', 'hetid', 'xyz', 'occupancy', 'variant')

    def __init__(self, atom_dict):
        self.name = atom_dict["auth_atom_id"]
        pos = atom_dict["auth_seq_id"]
        while len(pos) < 4:
            pos = " " + pos
        self.residue = "{} {}{}".format(
            atom_dict["auth_comp_id"], atom_dict["auth_asym_id"], pos)
        self.hetid = atom_dict["auth_comp_id"]
        self.xyz = (
            float(atom_dict["Cartn_x"]),
            float(atom_dict["Cartn_y"]),
            float(atom_dict["Cartn_z"]),
        )
        self.occupancy = float(atom_dict["occupancy"])
        self.variant = atom_dict["label_alt_id"]

    def __or__(self, other):
        """Return squared distance between two atoms."""
        sx, sy, sz = self.xyz
        ox, oy, oz = other.xyz
        return (sx - ox) ** 2 + (sy - oy) ** 2 + (sz - oz) ** 2

    def __hash__(self):
        return hash((self.residue, self.name, self.xyz))

    def __eq__(self, other):
        return (self.residue, self.name, self.xyz) == (other.residue, other.name, other.xyz)
