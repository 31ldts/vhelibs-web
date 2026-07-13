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


# ---------------------------------------------------------------------------
# get_default_entries
# ---------------------------------------------------------------------------

class TestGetDefaultEntries:
    def test_returns_list_of_dicts_with_expected_keys(self):
        entries = cofactors.get_default_entries()
        assert isinstance(entries, list)
        assert len(entries) > 0
        assert all(set(e) == {"code", "name", "category"} for e in entries)
        assert all(e["category"] in ("blacklist", "metal") for e in entries)

    def test_skips_empty_code_blacklist_entry(self):
        # DEFAULT_LIGAND_BLACKLIST[""] is a real (odd) entry, but it has
        # nothing meaningful to show/toggle in the UI, so it's excluded.
        entries = cofactors.get_default_entries()
        assert not any(e["code"] == "" for e in entries)

    def test_includes_known_blacklist_and_metal_entries(self):
        by_key = {(e["code"], e["category"]): e["name"] for e in cofactors.get_default_entries()}
        assert by_key[("HOH", "blacklist")] == "Water"
        assert by_key[("ZN", "metal")] == "Zinc"

    def test_sorted_by_category_then_code(self):
        entries = cofactors.get_default_entries()
        keys = [(e["category"], e["code"]) for e in entries]
        assert keys == sorted(keys)

    def test_total_count_matches_defaults_minus_empty_key(self):
        entries = cofactors.get_default_entries()
        expected = (len(cofactors.DEFAULT_LIGAND_BLACKLIST) - 1) + len(cofactors.DEFAULT_METALS)
        assert len(entries) == expected

    def test_reflects_builtin_defaults_not_mutated_module_globals(self):
        # get_default_entries is documented to always reflect the built-in
        # DEFAULT_* snapshots, never a single user's runtime customization
        # applied via update_lists().
        cofactors.update_lists({"XX": "Should not appear"}, {"YY": "Should not appear"})
        entries = cofactors.get_default_entries()
        codes = {e["code"] for e in entries}
        assert "XX" not in codes
        assert "YY" not in codes


# ---------------------------------------------------------------------------
# build_effective_lists
# ---------------------------------------------------------------------------

class TestBuildEffectiveLists:
    def test_no_args_returns_copies_of_defaults(self):
        m, lb = cofactors.build_effective_lists()
        assert m == cofactors.DEFAULT_METALS
        assert lb == cofactors.DEFAULT_LIGAND_BLACKLIST
        # Must be fresh instances so callers can't mutate the shared defaults.
        assert m is not cofactors.DEFAULT_METALS
        assert lb is not cofactors.DEFAULT_LIGAND_BLACKLIST

    def test_disabled_codes_are_removed_case_insensitively(self):
        m, lb = cofactors.build_effective_lists(disabled_codes=["zn", " HOH "])
        assert "ZN" not in m
        assert "HOH" not in lb
        # Untouched entries remain.
        assert "MG" in m
        assert "SO4" in lb

    def test_unknown_disabled_code_is_a_noop(self):
        m, lb = cofactors.build_effective_lists(disabled_codes=["NOPE"])
        assert m == cofactors.DEFAULT_METALS
        assert lb == cofactors.DEFAULT_LIGAND_BLACKLIST

    def test_custom_entry_defaults_to_blacklist_category_and_code_as_name(self):
        m, lb = cofactors.build_effective_lists(custom_entries=[{"code": "xyz"}])
        assert lb["XYZ"] == "XYZ"
        assert "XYZ" not in m

    def test_custom_entry_with_metal_category(self):
        m, lb = cofactors.build_effective_lists(
            custom_entries=[{"code": "xx", "name": "Test Metal", "category": "metal"}]
        )
        assert m["XX"] == "Test Metal"

    def test_custom_entries_applied_after_disabled_can_reintroduce_a_code(self):
        m, lb = cofactors.build_effective_lists(
            disabled_codes=["hoh"],
            custom_entries=[{"code": "HOH", "name": "Reintroduced", "category": "blacklist"}],
        )
        assert lb["HOH"] == "Reintroduced"

    def test_custom_entry_missing_code_is_skipped(self):
        m, lb = cofactors.build_effective_lists(custom_entries=[{"name": "no code"}])
        assert m == cofactors.DEFAULT_METALS
        assert lb == cofactors.DEFAULT_LIGAND_BLACKLIST

    def test_replace_fully_replaces_defaults(self):
        m, lb = cofactors.build_effective_lists(
            replace={"metals": {"aa": "Custom metal"}, "ligand_blacklist": {"bb": "Custom ligand"}}
        )
        assert m == {"AA": "Custom metal"}
        assert lb == {"BB": "Custom ligand"}

    def test_replace_makes_disabled_codes_irrelevant(self):
        # Starting from an uploaded file already means "use exactly these
        # entries" — disabled_codes is documented to be ignored in that case.
        m, lb = cofactors.build_effective_lists(
            disabled_codes=["aa"], replace={"metals": {"aa": "Custom metal"}}
        )
        assert m == {"AA": "Custom metal"}

    def test_replace_with_only_metals_key_leaves_blacklist_empty(self):
        m, lb = cofactors.build_effective_lists(replace={"metals": {"aa": "Custom"}})
        assert m == {"AA": "Custom"}
        assert lb == {}

    def test_replace_with_both_keys_empty_falls_back_to_defaults(self):
        m, lb = cofactors.build_effective_lists(replace={"metals": {}, "ligand_blacklist": {}})
        assert m == cofactors.DEFAULT_METALS
        assert lb == cofactors.DEFAULT_LIGAND_BLACKLIST

    def test_custom_entries_still_applied_on_top_of_replace(self):
        m, lb = cofactors.build_effective_lists(
            replace={"metals": {"aa": "Custom"}},
            custom_entries=[{"code": "bb", "name": "Extra", "category": "metal"}],
        )
        assert m == {"AA": "Custom", "BB": "Extra"}

    def test_returned_dicts_are_independent_across_calls(self):
        m1, _ = cofactors.build_effective_lists()
        m1["ZZ"] = "mutated"
        m2, _ = cofactors.build_effective_lists()
        assert "ZZ" not in m2


# ---------------------------------------------------------------------------
# parse_uploaded_list
# ---------------------------------------------------------------------------

class TestParseUploadedList:
    def test_simple_format_code_only_uses_code_as_name(self):
        m, lb = cofactors.parse_uploaded_list("ABC\nDEF")
        assert lb == {"ABC": "ABC", "DEF": "DEF"}
        assert m == {}

    @pytest.mark.parametrize("line,code,name", [
        ("ABC,Some ligand", "ABC", "Some ligand"),
        ("DEF;Another", "DEF", "Another"),
        ("GHI\tTabbed", "GHI", "Tabbed"),
    ])
    def test_simple_format_supports_comma_semicolon_and_tab_separators(self, line, code, name):
        m, lb = cofactors.parse_uploaded_list(line)
        assert lb[code] == name

    def test_simple_format_ignores_blank_and_comment_lines(self):
        text = "# a comment\n\nABC,Ligand\n"
        m, lb = cofactors.parse_uploaded_list(text)
        assert lb == {"ABC": "Ligand"}

    def test_simple_format_entries_always_categorized_as_blacklist(self):
        m, lb = cofactors.parse_uploaded_list("ABC,Something")
        assert "ABC" in lb
        assert "ABC" not in m

    def test_sectioned_format_routes_rows_to_correct_dict(self):
        text = "[Blacklist]\nBB,Blacklisted\n[Non-propagating]\nMM,A metal\n"
        m, lb = cofactors.parse_uploaded_list(text)
        assert lb == {"BB": "Blacklisted"}
        assert m == {"MM": "A metal"}

    def test_sectioned_format_ignores_rows_before_first_header(self):
        text = "ORPHAN,ignored\n[Blacklist]\nBB,Blacklisted\n"
        m, lb = cofactors.parse_uploaded_list(text)
        assert lb == {"BB": "Blacklisted"}
        assert m == {}

    def test_empty_text_returns_empty_dicts(self):
        m, lb = cofactors.parse_uploaded_list("")
        assert m == {}
        assert lb == {}

    def test_code_with_empty_name_falls_back_to_code(self):
        m, lb = cofactors.parse_uploaded_list("ABC,")
        assert lb["ABC"] == "ABC"

    def test_codes_are_uppercased(self):
        m, lb = cofactors.parse_uploaded_list("abc,lower name")
        assert "ABC" in lb
        assert "abc" not in lb

    def test_row_with_no_code_before_separator_is_skipped(self):
        m, lb = cofactors.parse_uploaded_list(",no code here")
        assert lb == {}

    def test_round_trip_with_dump_lists_output(self, tmp_path):
        # parse_uploaded_list's "sectioned" branch is documented to accept
        # exactly what dump_lists() produces (e.g. a previously saved,
        # re-uploaded file).
        cofactors.update_lists({"ZZ": "Zinc-like"}, {"AA": "A ligand"})
        target = tmp_path / "out.csv"
        cofactors.dump_lists(str(target))

        text = target.read_text()
        m, lb = cofactors.parse_uploaded_list(text)
        assert m == {"ZZ": "Zinc-like"}
        assert lb == {"AA": "A ligand"}
