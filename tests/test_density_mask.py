# -*- coding: utf-8 -*-
"""
Unit tests for core.density_mask.

External dependencies mocked/isolated:
  - core.eds_utils.get_edm / core.pdb_redo_utils.get_EDM: patched to return
    a small synthetic CCP4 map built with the real gemmi API (see
    make_synthetic_ccp4_map below) instead of downloading anything.
  - core.pdb_utils.CACHEDIR: monkeypatched to tmp_path, same convention as
    tests/test_routes.py's `client` fixture, so nothing touches a real
    on-disk cache directory.
  - gemmi itself is used for real (not mocked): the whole point of this
    module is what it does with gemmi (crop/mask/write), so faking that
    out would test nothing. See tests/test_rsr_core.py for the project's
    existing precedent of exercising real gemmi against small synthetic
    structures/maps built in-memory.
"""
import os
import time
import concurrent.futures
from unittest.mock import patch

import gemmi
import numpy as np
import pytest

import core.density_mask as density_mask
import core.pdb_utils as pdb_utils


# ---------------------------------------------------------------------------
# gemmi map-building helpers (test fixtures, not production code)
# ---------------------------------------------------------------------------

def make_synthetic_ccp4_map(path, side=40, spacing=1.0, seed=0):
    """Write a small synthetic P1 CCP4 map to `path` and return `path`.

    `side` Angstrom cubic cell at `spacing` Angstrom/voxel, filled with a
    smooth, non-trivial, deterministic pattern (never all-zero/all-one) so
    tests can tell masking actually changed something.
    """
    n = int(round(side / spacing))
    cell = gemmi.UnitCell(side, side, side, 90, 90, 90)
    grid = gemmi.FloatGrid(n, n, n)
    grid.set_unit_cell(cell)
    grid.spacegroup = gemmi.SpaceGroup("P1")

    arr = np.array(grid, copy=False)
    xs, ys, zs = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    rng = np.random.default_rng(seed)
    arr[:] = np.sin(xs / 3.0) + np.cos(ys / 4.0) + zs * 0.01 + rng.normal(0, 0.01, xs.shape)

    m = gemmi.Ccp4Map()
    m.grid = grid
    m.update_ccp4_header(2, True)
    m.write_ccp4_map(str(path))
    return str(path)


def read_nonzero_voxel_count(path):
    m = gemmi.read_ccp4_map(path)
    return int((np.array(m.grid, copy=False) != 0).sum())


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def cache_dir(tmp_path, monkeypatch):
    """Redirect core.pdb_utils.CACHEDIR to tmp_path for every test — same
    convention as tests/test_routes.py's `client` fixture."""
    monkeypatch.setattr(pdb_utils, "CACHEDIR", str(tmp_path))
    yield tmp_path


@pytest.fixture
def synthetic_map(tmp_path):
    return make_synthetic_ccp4_map(tmp_path / "source.ccp4")


@pytest.fixture
def region_box_atoms():
    box = {"min": [10.0, 10.0, 10.0], "max": [15.0, 15.0, 15.0]}
    atoms = [
        {"residue": "A/LIG`501/", "center": [12.0, 12.0, 12.0]},
        {"residue": "A/LIG`501/", "center": [13.5, 12.5, 13.0]},
    ]
    return box, atoms


# ---------------------------------------------------------------------------
# quantize_radius
# ---------------------------------------------------------------------------

class TestQuantizeRadius:
    @pytest.mark.parametrize("raw,expected", [
        (1.6, 1.5),
        (1.63, 1.75),
        (1.5, 1.5),
        (1.75, 1.75),
        (0.0, density_mask.RADIUS_MIN),
        (-3.0, density_mask.RADIUS_MIN),
        (0.1, density_mask.RADIUS_MIN),
        (999.0, density_mask.RADIUS_MAX),
    ])
    def test_snaps_to_step_and_clamps(self, raw, expected):
        assert density_mask.quantize_radius(raw) == expected

    def test_non_numeric_input_falls_back_to_minimum(self):
        assert density_mask.quantize_radius("not-a-number") == density_mask.RADIUS_MIN
        assert density_mask.quantize_radius(None) == density_mask.RADIUS_MIN

    def test_result_is_a_multiple_of_the_step(self):
        for raw in (0.4, 1.1, 2.37, 4.99):
            q = density_mask.quantize_radius(raw)
            # round() to sidestep float noise like 1.7500000000000002
            assert round(q / density_mask.RADIUS_STEP, 6) % 1 == 0


# ---------------------------------------------------------------------------
# get_masked_region_map — input validation
# ---------------------------------------------------------------------------

class TestGetMaskedRegionMapValidation:
    def test_unknown_region_returns_none(self, region_box_atoms):
        box, atoms = region_box_atoms
        assert density_mask.get_masked_region_map("1cbs", "not_a_region", box, atoms) is None

    def test_missing_box_returns_none(self, region_box_atoms):
        _box, atoms = region_box_atoms
        assert density_mask.get_masked_region_map("1cbs", "ligand", {}, atoms) is None
        assert density_mask.get_masked_region_map("1cbs", "ligand", None, atoms) is None

    def test_box_missing_min_or_max_returns_none(self, region_box_atoms):
        box, atoms = region_box_atoms
        assert density_mask.get_masked_region_map(
            "1cbs", "ligand", {"min": box["min"]}, atoms) is None
        assert density_mask.get_masked_region_map(
            "1cbs", "ligand", {"max": box["max"]}, atoms) is None

    def test_empty_atoms_returns_none(self, region_box_atoms):
        box, _atoms = region_box_atoms
        assert density_mask.get_masked_region_map("1cbs", "ligand", box, []) is None
        assert density_mask.get_masked_region_map("1cbs", "ligand", box, None) is None

    @patch("core.eds_utils.get_edm")
    def test_invalid_input_never_touches_the_source_map(self, mock_get_edm, region_box_atoms):
        # Validation must short-circuit before any download/gemmi work.
        density_mask.get_masked_region_map("1cbs", "ligand", {}, [])
        mock_get_edm.assert_not_called()

    @patch("core.eds_utils.get_edm", return_value=(None, None))
    def test_no_source_map_available_returns_none(self, mock_get_edm, region_box_atoms):
        box, atoms = region_box_atoms
        assert density_mask.get_masked_region_map("1cbs", "ligand", box, atoms) is None

    @patch("core.pdb_redo_utils.get_EDM", return_value=None)
    def test_no_pdb_redo_map_available_returns_none(self, mock_get_edm, region_box_atoms):
        box, atoms = region_box_atoms
        result = density_mask.get_masked_region_map(
            "1cbs", "ligand", box, atoms, source="pdb_redo")
        assert result is None


# ---------------------------------------------------------------------------
# get_masked_region_map — real gemmi masking behaviour
# ---------------------------------------------------------------------------

class TestGetMaskedRegionMapMasking:
    @patch("core.eds_utils.get_edm")
    def test_produces_a_smaller_masked_file(self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, atoms = region_box_atoms

        out = density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=1.6)

        assert out is not None
        assert os.path.isfile(out)
        assert os.path.getsize(out) < os.path.getsize(synthetic_map)

    @patch("core.eds_utils.get_edm")
    def test_masked_file_has_nonzero_but_not_all_voxels(self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, atoms = region_box_atoms

        out = density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=1.6)

        m = gemmi.read_ccp4_map(out)
        arr = np.array(m.grid, copy=False)
        nonzero = int((arr != 0).sum())
        assert 0 < nonzero < arr.size, (
            "masked output should keep some density (near the atoms) and "
            "zero out the rest, not be all-zero or fully unmasked"
        )

    @patch("core.eds_utils.get_edm")
    def test_larger_radius_masks_more_voxels(self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, atoms = region_box_atoms

        small = density_mask.get_masked_region_map(
            "1cbs", "ligand", box, atoms, radius=0.5, use_cache=False)
        large = density_mask.get_masked_region_map(
            "1cbs", "ligand", box, atoms, radius=3.0, use_cache=False)

        assert read_nonzero_voxel_count(small) < read_nonzero_voxel_count(large)

    @patch("core.eds_utils.get_edm")
    def test_output_is_a_valid_ccp4_map_matching_the_cell(self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, atoms = region_box_atoms

        out = density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=1.6)

        source_cell = gemmi.read_ccp4_map(synthetic_map).grid.unit_cell
        out_cell = gemmi.read_ccp4_map(out).grid.unit_cell
        # Cropping shrinks the stored extent, not the crystallographic
        # unit cell itself.
        assert (source_cell.a, source_cell.b, source_cell.c) == (out_cell.a, out_cell.b, out_cell.c)

    @patch("core.pdb_redo_utils.get_EDM")
    def test_pdb_redo_source_is_masked_too(self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = synthetic_map  # bare path, not a tuple — see get_EDM's real signature
        box, atoms = region_box_atoms

        out = density_mask.get_masked_region_map(
            "1cbs", "ligand", box, atoms, radius=1.6, source="pdb_redo")

        assert out is not None and os.path.isfile(out)
        mock_get_edm.assert_called_once_with("1cbs", use_cache=True)

    @patch("core.pdb_redo_utils.get_EDM")
    def test_pdb_redo_map_is_cached_after_first_download(
            self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = synthetic_map
        box, atoms = region_box_atoms

        density_mask.get_masked_region_map(
            "1cbs", "ligand", box, atoms, radius=1.6, source="pdb_redo")
        density_mask.get_masked_region_map(
            "1cbs", "ligand", box, atoms, radius=1.6, source="pdb_redo")

        # The masked-region cache alone already prevents a 2nd call for
        # identical box/atoms/radius — get_EDM's own on-disk cache is
        # exercised directly in tests/test_pdb_redo_utils.py.
        mock_get_edm.assert_called_once()


# ---------------------------------------------------------------------------
# get_masked_region_map — caching
# ---------------------------------------------------------------------------

class TestGetMaskedRegionMapCaching:
    @patch("core.eds_utils.get_edm")
    def test_second_call_reuses_cache_without_hitting_get_edm_again(
            self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, atoms = region_box_atoms

        out1 = density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=1.6)
        out2 = density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=1.6)

        assert out1 == out2
        mock_get_edm.assert_called_once()

    @patch("core.eds_utils.get_edm")
    def test_use_cache_false_recomputes(self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, atoms = region_box_atoms

        density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=1.6, use_cache=True)
        density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=1.6, use_cache=False)

        assert mock_get_edm.call_count == 2

    @patch("core.eds_utils.get_edm")
    def test_distinct_ligands_in_the_same_region_do_not_collide(
            self, mock_get_edm, synthetic_map):
        # Regression test: a PDB entry with two different ligands produces
        # two different box/atoms for the SAME region name ("ligand").
        # The cache key must disambiguate them, or the second ligand's
        # viewer request would silently be served the first ligand's map.
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box_a = {"min": [10.0, 10.0, 10.0], "max": [15.0, 15.0, 15.0]}
        atoms_a = [{"center": [12.0, 12.0, 12.0]}]
        box_b = {"min": [20.0, 20.0, 20.0], "max": [24.0, 24.0, 24.0]}
        atoms_b = [{"center": [22.0, 22.0, 22.0]}]

        out_a = density_mask.get_masked_region_map("1cbs", "ligand", box_a, atoms_a, radius=1.6)
        out_b = density_mask.get_masked_region_map("1cbs", "ligand", box_b, atoms_b, radius=1.6)

        assert out_a != out_b
        assert os.path.isfile(out_a) and os.path.isfile(out_b)

    @patch("core.eds_utils.get_edm")
    def test_same_ligand_requested_twice_hits_the_same_cache_file(
            self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, atoms = region_box_atoms

        out1 = density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=1.6)
        # Same coordinates, freshly-constructed dicts/lists (not the same
        # objects) — the cache key must be based on the *values*.
        box_copy = {"min": list(box["min"]), "max": list(box["max"])}
        atoms_copy = [{"center": list(a["center"])} for a in atoms]
        out2 = density_mask.get_masked_region_map("1cbs", "ligand", box_copy, atoms_copy, radius=1.6)

        assert out1 == out2
        mock_get_edm.assert_called_once()

    @patch("core.eds_utils.get_edm")
    def test_invalid_request_does_not_return_a_stale_cached_file(
            self, mock_get_edm, synthetic_map, region_box_atoms):
        # Regression test: an earlier bug checked the on-disk cache before
        # validating box/atoms, so an empty/invalid request for a
        # (pdbid, region, radius) that had already been served correctly
        # could return that unrelated cached file instead of None.
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, atoms = region_box_atoms
        density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=1.6)

        assert density_mask.get_masked_region_map("1cbs", "ligand", box, [], radius=1.6) is None
        assert density_mask.get_masked_region_map("1cbs", "ligand", {}, atoms, radius=1.6) is None

    @patch("core.eds_utils.get_edm")
    def test_different_radii_produce_different_cache_files(
            self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, atoms = region_box_atoms

        out_small = density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=0.5)
        out_large = density_mask.get_masked_region_map("1cbs", "ligand", box, atoms, radius=3.0)

        assert out_small != out_large

    @patch("core.eds_utils.get_edm")
    def test_pdb_and_pdb_redo_sources_are_cached_separately(
            self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, atoms = region_box_atoms

        with patch("core.pdb_redo_utils.get_EDM", return_value=synthetic_map):
            out_pdb = density_mask.get_masked_region_map(
                "1cbs", "ligand", box, atoms, radius=1.6, source="pdb")
            out_redo = density_mask.get_masked_region_map(
                "1cbs", "ligand", box, atoms, radius=1.6, source="pdb_redo")

        assert out_pdb != out_redo


# ---------------------------------------------------------------------------
# get_masked_region_map — error handling
# ---------------------------------------------------------------------------

class TestGetMaskedRegionMapErrorHandling:
    @patch("core.eds_utils.get_edm")
    def test_corrupt_source_map_returns_none_not_an_exception(
            self, mock_get_edm, tmp_path, region_box_atoms):
        bad_map = tmp_path / "corrupt.ccp4"
        bad_map.write_bytes(b"not a real ccp4 file")
        mock_get_edm.return_value = (str(bad_map), 1.0)
        box, atoms = region_box_atoms

        assert density_mask.get_masked_region_map("1cbs", "ligand", box, atoms) is None

    @patch("core.eds_utils.get_edm")
    def test_malformed_atom_center_returns_none_not_an_exception(
            self, mock_get_edm, synthetic_map, region_box_atoms):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        box, _atoms = region_box_atoms
        bad_atoms = [{"center": [1.0, 2.0]}]  # missing z

        assert density_mask.get_masked_region_map("1cbs", "ligand", box, bad_atoms) is None

    @patch("core.eds_utils.get_edm")
    def test_atom_list_is_truncated_to_max_mask_atoms(
            self, mock_get_edm, synthetic_map, region_box_atoms, monkeypatch):
        mock_get_edm.return_value = (synthetic_map, 1.0)
        monkeypatch.setattr(density_mask, "MAX_MASK_ATOMS", 2)
        box, _atoms = region_box_atoms
        many_atoms = [{"center": [12.0 + 0.01 * i, 12.0, 12.0]} for i in range(10)]

        # Must not raise despite far more atoms than MAX_MASK_ATOMS, and
        # must still produce a valid file.
        out = density_mask.get_masked_region_map("1cbs", "ligand", box, many_atoms, use_cache=False)
        assert out is not None and os.path.isfile(out)


# ---------------------------------------------------------------------------
# prefetch_default_masks
# ---------------------------------------------------------------------------

class TestPrefetchDefaultMasks:
    def test_calls_get_masked_region_map_for_every_populated_region(self):
        density_boxes = {
            "ligand": {"min": [0, 0, 0], "max": [1, 1, 1]},
            "binding_site": {"min": [0, 0, 0], "max": [1, 1, 1]},
            "residues_to_examine": None,  # not present for this ligand
        }
        density_atoms = {
            "ligand": [{"center": [0.5, 0.5, 0.5]}],
            "binding_site": [{"center": [0.5, 0.5, 0.5]}],
            "residues_to_examine": [],
        }
        with patch.object(density_mask, "get_masked_region_map") as mock_get:
            density_mask.prefetch_default_masks("1cbs", density_boxes, density_atoms, executor=None)

        called_regions = {c.args[1] for c in mock_get.call_args_list}
        assert called_regions == {"ligand", "binding_site"}

    def test_skips_regions_with_missing_box_or_atoms(self):
        with patch.object(density_mask, "get_masked_region_map") as mock_get:
            density_mask.prefetch_default_masks("1cbs", {}, {}, executor=None)
        mock_get.assert_not_called()

    def test_none_boxes_or_atoms_does_not_raise(self):
        with patch.object(density_mask, "get_masked_region_map") as mock_get:
            density_mask.prefetch_default_masks("1cbs", None, None, executor=None)
        mock_get.assert_not_called()

    def test_uses_default_radius_and_pdb_source_unless_overridden(self):
        density_boxes = {"ligand": {"min": [0, 0, 0], "max": [1, 1, 1]}}
        density_atoms = {"ligand": [{"center": [0.5, 0.5, 0.5]}]}
        with patch.object(density_mask, "get_masked_region_map") as mock_get:
            density_mask.prefetch_default_masks("1cbs", density_boxes, density_atoms, executor=None)

        call = mock_get.call_args
        assert call.kwargs.get("radius", 1.6) == 1.6
        assert call.kwargs.get("source", "pdb") == "pdb"

    def test_submits_to_executor_instead_of_running_inline(self):
        density_boxes = {"ligand": {"min": [0, 0, 0], "max": [1, 1, 1]}}
        density_atoms = {"ligand": [{"center": [0.5, 0.5, 0.5]}]}

        with patch.object(density_mask, "get_masked_region_map") as mock_get:
            with concurrent.futures.ThreadPoolExecutor(max_workers=1) as ex:
                density_mask.prefetch_default_masks(
                    "1cbs", density_boxes, density_atoms, executor=ex)
                ex.shutdown(wait=True)

        mock_get.assert_called_once()

    def test_synchronous_mode_propagates_masking_errors(self):
        # prefetch_default_masks(executor=None) runs inline and does NOT
        # itself add error handling — it relies on get_masked_region_map
        # never raising in practice (see that function's own try/except).
        # The "never breaks the caller" guarantee lives one layer up, in
        # app.routes._prefetch_density_masks, which is what real callers
        # use — see tests/test_routes.py for that guarantee being tested.
        with patch.object(density_mask, "get_masked_region_map", side_effect=RuntimeError("boom")):
            density_boxes = {"ligand": {"min": [0, 0, 0], "max": [1, 1, 1]}}
            density_atoms = {"ligand": [{"center": [0.5, 0.5, 0.5]}]}
            with pytest.raises(RuntimeError):
                density_mask.prefetch_default_masks(
                    "1cbs", density_boxes, density_atoms, executor=None)
