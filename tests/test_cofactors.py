# -*- coding: utf-8 -*-
"""
Unit tests for core.cofactors.

This module only holds two static dicts (`metals`, `ligand_blacklist`) and
three functions to swap/dump/load them. There is no network access; the
only external dependency is the filesystem (CSV read/write), which we
redirect to tmp_path instead of mocking, since these are simple, fast,
real file operations and mocking `open`/`csv` would add noise without
value.
"""
import csv

import pytest

import core.cofactors as cofactors


@pytest.fixture(autouse=True)
def restore_module_state():
    """
    cofactors.update_lists() mutates module-level globals. Snapshot and
    restore them around every test so tests don't leak state into each
    other (the module is a singleton shared across the whole test suite).
    """
    original_metals = dict(cofactors.metals)
    original_blacklist = dict(cofactors.ligand_blacklist)
    yield
    cofactors.metals = original_metals
    cofactors.ligand_blacklist = original_blacklist


class TestStaticData:
    def test_metals_is_a_nonempty_dict(self):
        assert isinstance(cofactors.metals, dict)
        assert len(cofactors.metals) > 0

    def test_ligand_blacklist_is_a_nonempty_dict(self):
        assert isinstance(cofactors.ligand_blacklist, dict)
        assert len(cofactors.ligand_blacklist) > 0

    def test_known_metal_entries(self):
        # Spot-check a couple of well-known entries rather than the whole dict.
        assert cofactors.metals["ZN"] == "Zinc"
        assert cofactors.metals["MG"] == "MAGNESIUM ION"
        assert cofactors.metals["FE"] == "Fe(2+)"

    def test_known_blacklist_entries(self):
        assert cofactors.ligand_blacklist["HOH"] == "Water"
        assert cofactors.ligand_blacklist["SO4"] == "Sulfate Ion"
        # The empty-string key is a real (odd) entry in the current data.
        assert cofactors.ligand_blacklist[""] == ""


class TestUpdateLists:
    def test_update_lists_with_no_args_is_a_noop_on_content(self):
        # Calling with no args re-assigns metals/ligand_blacklist using the
        # function's default arguments (bound at definition time) — current
        # behaviour is that content stays equal to what it was before the
        # call, since nothing new is passed in.
        before_metals = dict(cofactors.metals)
        before_blacklist = dict(cofactors.ligand_blacklist)
        result = cofactors.update_lists()
        assert result is None
        assert cofactors.metals == before_metals
        assert cofactors.ligand_blacklist == before_blacklist

    def test_update_lists_with_explicit_dicts(self):
        new_metals = {"XX": "Test metal"}
        new_blacklist = {"YY": "Test ligand"}
        cofactors.update_lists(new_metals, new_blacklist)
        assert cofactors.metals == {"XX": "Test metal"}
        assert cofactors.ligand_blacklist == {"YY": "Test ligand"}

    def test_update_lists_stores_same_object_reference(self):
        new_metals = {"XX": "Test metal"}
        new_blacklist = {"YY": "Test ligand"}
        cofactors.update_lists(new_metals, new_blacklist)
        assert cofactors.metals is new_metals
        assert cofactors.ligand_blacklist is new_blacklist


class TestDumpLists:
    def test_dump_lists_appends_csv_extension_when_missing(self, tmp_path):
        target = tmp_path / "mylists"
        result = cofactors.dump_lists(str(target))
        assert result == 0
        assert (tmp_path / "mylists.csv").is_file()
        assert not target.is_file()  # no extension-less file was created

    def test_dump_lists_does_not_double_append_extension(self, tmp_path):
        target = tmp_path / "mylists.csv"
        cofactors.dump_lists(str(target))
        assert target.is_file()
        assert not (tmp_path / "mylists.csv.csv").exists()

    def test_dump_lists_default_filename(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        result = cofactors.dump_lists()
        assert result == 0
        assert (tmp_path / "notligands.csv").is_file()

    def test_dump_lists_content_structure(self, tmp_path):
        cofactors.update_lists({"ZZ": "Zinc-like"}, {"AA": "A ligand"})
        target = tmp_path / "out.csv"
        cofactors.dump_lists(str(target))

        with open(target, newline="") as fh:
            rows = list(csv.reader(fh))

        assert rows[0] == ["[Blacklist]"]
        assert ["AA", "A ligand"] in rows
        assert ["[Non-propagating]"] in [rows[i] for i in range(len(rows))] or \
            any(r == ["[Non-propagating]"] for r in rows)
        assert ["ZZ", "Zinc-like"] in rows


class TestLoadLists:
    def test_load_lists_round_trip(self, tmp_path):
        cofactors.update_lists({"ZZ": "Zinc-like"}, {"AA": "A ligand"})
        target = tmp_path / "out.csv"
        cofactors.dump_lists(str(target))

        # Overwrite in-memory state, then reload from disk and check it's restored.
        cofactors.update_lists({}, {})
        result = cofactors.load_lists(str(target))

        assert result == 0
        assert cofactors.metals == {"ZZ": "Zinc-like"}
        assert cofactors.ligand_blacklist == {"AA": "A ligand"}

    def test_load_lists_ignores_rows_before_first_section_header(self, tmp_path):
        target = tmp_path / "custom.csv"
        with open(target, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["ORPHAN", "should be ignored, no section set yet"])
            w.writerow(["[Blacklist]"])
            w.writerow(["BB", "Blacklisted"])
            w.writerow(["[Non-propagating]"])
            w.writerow(["MM", "A metal"])

        cofactors.load_lists(str(target))

        assert cofactors.ligand_blacklist == {"BB": "Blacklisted"}
        assert cofactors.metals == {"MM": "A metal"}

    def test_load_lists_skips_blank_rows(self, tmp_path):
        target = tmp_path / "with_blanks.csv"
        with open(target, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["[Blacklist]"])
            w.writerow([])
            w.writerow(["BB", "Blacklisted"])
            w.writerow(["[Non-propagating]"])

        # Should not raise despite the blank row.
        result = cofactors.load_lists(str(target))
        assert result == 0
        assert cofactors.ligand_blacklist == {"BB": "Blacklisted"}
