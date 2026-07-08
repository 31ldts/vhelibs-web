# -*- coding: utf-8 -*-
#
#   Copyright 2010-2024 Adrià Cereto Massagué
#   Migrated to web version
#


def format_reskey(comp_id, asym_id, seq_id, icode=""):
    """Build a canonical, unambiguous residue key.

    Produces a key of the form ``"{comp_id:>3} {asym_id} {seq_id:>4}{icode}"``.

    This is the single source of truth for residue-key formatting, used by
    this module, core.rsr_core, core.eds_utils and core.pdb_redo_utils.
    They must all agree byte-for-byte, since these keys are compared by
    plain string equality to cross-reference validation stats (edd_dict)
    against parsed atoms (res_atom_dict / ligand_res_atom_dict).

    Args:
        comp_id (str): Component (residue) identifier, e.g. the residue
            name. Left-padded with spaces to a minimum width of 3
            characters.
        asym_id (str): Author-assigned asymmetric unit (chain) identifier.
        seq_id (int or str): Author-assigned sequence number of the
            residue. Left-padded with spaces to a minimum width of 4
            characters.
        icode (str, optional): PDB insertion code. Defaults to ``""``.

    Returns:
        str: The formatted, canonical residue key.

    Raises:
        None
    """
    comp_id = str(comp_id)
    while len(comp_id) < 3:
        comp_id = " " + comp_id
    seq_id = str(seq_id)
    while len(seq_id) < 4:
        seq_id = " " + seq_id
    return "{} {} {}{}".format(comp_id, asym_id, seq_id, icode or "")


class PdbAtom:
    """Represents an atom from a PDB/mmCIF file.

    Attributes:
        name (str): Author-assigned atom name.
        residue (str): Canonical residue key for the atom's parent
            residue, as produced by :func:`format_reskey`.
        hetid (str): Author-assigned component (residue) identifier.
        xyz (tuple): 3-tuple of floats ``(x, y, z)`` giving the atom's
            Cartesian coordinates.
        occupancy (float): Crystallographic occupancy of the atom.
        variant (str): Alternate location (altloc) identifier.
    """

    __slots__ = ('name', 'residue', 'hetid', 'xyz', 'occupancy', 'variant')

    def __init__(self, atom_dict):
        """Initialize a PdbAtom from a raw mmCIF atom record.

        Args:
            atom_dict (dict): Dictionary of mmCIF ``_atom_site`` fields for
                a single atom. Must contain the keys ``"auth_atom_id"``,
                ``"auth_comp_id"``, ``"auth_asym_id"``, ``"auth_seq_id"``,
                ``"Cartn_x"``, ``"Cartn_y"``, ``"Cartn_z"``,
                ``"occupancy"``, and ``"label_alt_id"``. The key
                ``"pdbx_PDB_ins_code"`` is optional.

        Returns:
            None

        Raises:
            KeyError: If a required key is missing from ``atom_dict``.
            ValueError: If ``Cartn_x``, ``Cartn_y``, ``Cartn_z``, or
                ``occupancy`` cannot be converted to ``float``.
        """
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
        """Compute the squared Euclidean distance to another atom.

        Args:
            other (PdbAtom): The other atom to measure the distance to.

        Returns:
            float: The squared distance between this atom's and ``other``'s
            Cartesian coordinates.

        Raises:
            AttributeError: If ``other`` does not have an ``xyz``
                attribute.
        """
        sx, sy, sz = self.xyz
        ox, oy, oz = other.xyz
        return (sx - ox) ** 2 + (sy - oy) ** 2 + (sz - oz) ** 2

    def __hash__(self):
        """Compute a hash based on residue, atom name, and coordinates.

        Returns:
            int: Hash value derived from the tuple
            ``(self.residue, self.name, self.xyz)``.

        Raises:
            None
        """
        return hash((self.residue, self.name, self.xyz))

    def __eq__(self, other):
        """Compare two atoms for equality.

        Two atoms are considered equal if they share the same residue key,
        atom name, and Cartesian coordinates.

        Args:
            other (PdbAtom): The other atom to compare against.

        Returns:
            bool: ``True`` if both atoms have the same
            ``(residue, name, xyz)`` tuple, ``False`` otherwise.

        Raises:
            AttributeError: If ``other`` does not have ``residue``,
                ``name``, or ``xyz`` attributes.
        """
        return (self.residue, self.name, self.xyz) == (other.residue, other.name, other.xyz)
