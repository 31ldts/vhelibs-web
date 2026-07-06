# -*- coding: utf-8 -*-
#
#   Copyright 2010-2024 Adrià Cereto Massagué
#   Migrated to web version – pure Python fallback is now the primary implementation.
#


def format_reskey(comp_id, asym_id, seq_id, icode=""):
    """
    Canonical, unambiguous residue key: "{comp_id:>3} {asym_id} {seq_id:>4}{icode}".

    This is the single source of truth for residue-key formatting, used by
    this module, core.rsr_core, core.eds_utils and core.pdb_redo_utils.
    They must all agree byte-for-byte, since these keys are compared by
    plain string equality to cross-reference validation stats (edd_dict)
    against parsed atoms (res_atom_dict / ligand_res_atom_dict).

    IMPORTANT: asym_id and seq_id are always joined with a literal space,
    regardless of how many digits seq_id has. A previous version relied on
    seq_id's own left-padding (e.g. " 200") to *look* separated from
    asym_id, which silently broke for any 4+-digit residue number (no
    padding left to add), fusing asym_id and seq_id into a single token
    (e.g. "A1000"). That corrupted every downstream string-based parser
    that assumed a fixed 3-field layout — most visibly, the frontend's
    residueSelector(), which would then fail to build a selector for that
    residue, fall back to selecting *every* ligand in the structure, and
    show ligands that should have stayed hidden.
    """
    comp_id = str(comp_id)
    while len(comp_id) < 3:
        comp_id = " " + comp_id
    seq_id = str(seq_id)
    while len(seq_id) < 4:
        seq_id = " " + seq_id
    return "{} {} {}{}".format(comp_id, asym_id, seq_id, icode or "")


class PdbAtom:
    """Represents an atom from a PDB/mmCIF file."""

    __slots__ = ('name', 'residue', 'hetid', 'xyz', 'occupancy', 'variant')

    def __init__(self, atom_dict):
        self.name = atom_dict["auth_atom_id"]
        self.residue = format_reskey(
            atom_dict["auth_comp_id"],
            atom_dict["auth_asym_id"],
            atom_dict["auth_seq_id"],
            atom_dict.get("pdbx_PDB_ins_code", ""),
        )
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
