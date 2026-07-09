# -*- coding: utf-8 -*-
"""
Unit tests for core.pdb_utils (post-refactor: download/cache logic now
delegates to core.http_cache).

Network access happens inside core.http_cache (requests.get for
downloads/JSON fetches) except for the direct UniProt search POST, which
pdb_utils still issues itself. Both are mocked; CACHEDIR is redirected to
tmp_path.
"""
import os
from unittest.mock import MagicMock, patch

import pytest
import requests

import core.pdb_utils as pdb_utils


@pytest.fixture(autouse=True)
def isolated_cache_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(pdb_utils, "CACHEDIR", str(tmp_path))
    yield tmp_path


# ---------------------------------------------------------------------------
# set_cache_dir
# ---------------------------------------------------------------------------

class TestSetCacheDir:
    def test_sets_module_global_and_creates_dir(self, tmp_path):
        new_dir = tmp_path / "new_cache"
        pdb_utils.set_cache_dir(str(new_dir))
        assert pdb_utils.CACHEDIR == str(new_dir)
        assert new_dir.is_dir()

    def test_idempotent_on_existing_dir(self, tmp_path):
        existing = tmp_path / "existing"
        existing.mkdir()
        pdb_utils.set_cache_dir(str(existing))  # should not raise
        assert pdb_utils.CACHEDIR == str(existing)


# ---------------------------------------------------------------------------
# _get_nested
# ---------------------------------------------------------------------------

class TestGetNested:
    def test_full_path_resolves(self):
        d = {"a": {"b": {"c": 42}}}
        assert pdb_utils._get_nested(d, "a", "b", "c") == 42

    def test_missing_intermediate_key_returns_default(self):
        d = {"a": {}}
        assert pdb_utils._get_nested(d, "a", "b", "c") is None

    def test_custom_default(self):
        d = {}
        assert pdb_utils._get_nested(d, "a", default="fallback") == "fallback"

    def test_non_dict_intermediate_value_returns_default(self):
        d = {"a": "not-a-dict"}
        assert pdb_utils._get_nested(d, "a", "b") is None

    def test_none_value_short_circuits(self):
        d = {"a": None}
        assert pdb_utils._get_nested(d, "a", "b", default="x") == "x"

    def test_no_keys_returns_input_dict(self):
        d = {"a": 1}
        assert pdb_utils._get_nested(d) == d


# ---------------------------------------------------------------------------
# _extract_report_fields
# ---------------------------------------------------------------------------

RCSB_ENTRY_JSON = {
    "rcsb_entry_info": {"experimental_method": "X-RAY DIFFRACTION"},
    "refine": [{
        "ls_R_factor_R_free": 0.218,
        "ls_R_factor_R_work": 0.175,
        "ls_d_res_high": 1.8,
        "ls_number_reflns_rfree": 1234,
    }],
    "cell": {
        "angle_alpha": 90.0, "angle_beta": 90.0, "angle_gamma": 90.0,
        "length_a": 50.0, "length_b": 60.0, "length_c": 70.0,
    },
}


class TestExtractReportFields:
    def test_extracts_all_fields(self):
        row = pdb_utils._extract_report_fields(RCSB_ENTRY_JSON)
        assert row["rFree"] == 0.218
        assert row["rWork"] == 0.175
        assert row["refinementResolution"] == 1.8
        assert row["nreflections"] == 1234
        assert row["experimentalTechnique"] == "X-RAY DIFFRACTION"
        assert row["unitCellAngleAlpha"] == 90.0
        assert row["lengthOfUnitCellLatticeA"] == 50.0

    def test_falls_back_to_ls_d_res_low_when_high_missing(self):
        data = {
            "refine": [{"ls_d_res_low": 2.5}],
        }
        row = pdb_utils._extract_report_fields(data)
        assert row["refinementResolution"] == 2.5

    def test_missing_refine_block_yields_placeholder_rfree(self):
        data = {"rcsb_entry_info": {"experimental_method": "SOLUTION NMR"}}
        row = pdb_utils._extract_report_fields(data)
        assert row["rFree"] == 9999
        assert row["rWork"] == 9999
        assert row["refinementResolution"] == 0.0
        assert row["experimentalTechnique"] == "SOLUTION NMR"

    def test_malformed_refine_raises_for_caller_to_catch(self):
        with pytest.raises(Exception):
            pdb_utils._extract_report_fields({"refine": "not-a-list"})


# ---------------------------------------------------------------------------
# get_custom_report
# ---------------------------------------------------------------------------

class TestGetCustomReport:
    @patch("core.http_cache.requests.get")
    def test_fetches_and_parses_report(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = RCSB_ENTRY_JSON
        mock_get.return_value = resp

        result = pdb_utils.get_custom_report("1cbs")

        assert "1CBS" in result
        row = result["1CBS"]
        assert row["rFree"] == 0.218
        assert row["refinementResolution"] == 1.8

    @patch("core.http_cache.requests.get")
    def test_uses_disk_cache_on_second_call(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = RCSB_ENTRY_JSON
        mock_get.return_value = resp

        pdb_utils.get_custom_report("1cbs")
        pdb_utils.get_custom_report("1cbs")

        assert mock_get.call_count == 1

    @patch("core.http_cache.requests.get")
    def test_network_failure_returns_empty_dict(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        result = pdb_utils.get_custom_report("1cbs")
        assert result == {}

    @patch("core.http_cache.requests.get")
    def test_malformed_json_is_handled_gracefully(self, mock_get):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {"refine": "not-a-list-of-dicts"}
        mock_get.return_value = resp

        result = pdb_utils.get_custom_report("3xyz")
        assert result == {}

    @patch("core.http_cache.requests.get")
    def test_corrupt_cache_file_triggers_refetch(self, mock_get, tmp_path):
        pdbid = "1cbs"
        cache_dir = os.path.join(str(tmp_path), pdbid)
        os.makedirs(cache_dir, exist_ok=True)
        with open(os.path.join(cache_dir, "pdb_stats.json"), "w") as fh:
            fh.write("{not valid json")

        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = RCSB_ENTRY_JSON
        mock_get.return_value = resp

        result = pdb_utils.get_custom_report(pdbid)
        assert "1CBS" in result
        mock_get.assert_called_once()

    @patch("core.http_cache.requests.get")
    def test_request_uses_lowercased_pdbid_in_url(self, mock_get):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = RCSB_ENTRY_JSON
        mock_get.return_value = resp

        pdb_utils.get_custom_report("1CBS")
        called_url = mock_get.call_args[0][0]
        assert called_url == "https://data.rcsb.org/rest/v1/core/entry/1cbs"


# ---------------------------------------------------------------------------
# _build_uniprot_query
# ---------------------------------------------------------------------------

class TestBuildUniprotQuery:
    def test_builds_expected_shape(self):
        query = pdb_utils._build_uniprot_query("P00734", 50)
        assert query["query"]["parameters"]["value"] == "P00734"
        assert query["query"]["parameters"]["operator"] == "exact_match"
        assert query["return_type"] == "entry"
        assert query["request_options"]["paginate"]["rows"] == 50


# ---------------------------------------------------------------------------
# get_pdbids_for_uniprot
# ---------------------------------------------------------------------------

class TestGetPdbidsForUniprot:
    @patch("core.pdb_utils.requests.post")
    def test_returns_sorted_unique_pdbids(self, mock_post):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {
            "result_set": [
                {"identifier": "3dzu"},
                {"identifier": "1CBS"},
                {"identifier": "1cbs"},
            ]
        }
        mock_post.return_value = resp

        result = pdb_utils.get_pdbids_for_uniprot("p00734")

        assert result == ["1CBS", "3DZU"]
        called_url = mock_post.call_args[0][0]
        assert called_url == pdb_utils.UNIPROT_SEARCH_URL
        sent_query = mock_post.call_args[1]["json"]
        assert sent_query["query"]["parameters"]["value"] == "P00734"

    @patch("core.pdb_utils.requests.post")
    def test_204_means_no_hits(self, mock_post):
        resp = MagicMock()
        resp.status_code = 204
        mock_post.return_value = resp

        result = pdb_utils.get_pdbids_for_uniprot("P99999")
        assert result == []

    @patch("core.pdb_utils.requests.post")
    def test_network_error_returns_empty_list(self, mock_post):
        mock_post.side_effect = requests.exceptions.ConnectionError("down")
        result = pdb_utils.get_pdbids_for_uniprot("P00734")
        assert result == []

    @patch("core.pdb_utils.requests.post")
    def test_uses_cache_on_second_call(self, mock_post):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {"result_set": [{"identifier": "1CBS"}]}
        mock_post.return_value = resp

        first = pdb_utils.get_pdbids_for_uniprot("P00734")
        second = pdb_utils.get_pdbids_for_uniprot("P00734")

        assert first == second == ["1CBS"]
        assert mock_post.call_count == 1

    @patch("core.pdb_utils.requests.post")
    def test_ignores_hits_without_identifier(self, mock_post):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {"result_set": [{"identifier": None}, {}]}
        mock_post.return_value = resp

        result = pdb_utils.get_pdbids_for_uniprot("P00000")
        assert result == []

    @patch("core.pdb_utils.requests.post")
    def test_corrupt_cache_file_triggers_refetch(self, mock_post, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "uniprot")
        os.makedirs(cache_dir, exist_ok=True)
        with open(os.path.join(cache_dir, "P00734.json"), "w") as fh:
            fh.write("{not valid json")

        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {"result_set": [{"identifier": "1CBS"}]}
        mock_post.return_value = resp

        result = pdb_utils.get_pdbids_for_uniprot("P00734")
        assert result == ["1CBS"]
        mock_post.assert_called_once()


# ---------------------------------------------------------------------------
# get_pdb_file
# ---------------------------------------------------------------------------

class TestGetPdbFile:
    @patch("core.http_cache.requests.get")
    def test_builds_rcsb_url_and_filename_by_default(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = b"data"
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        result = pdb_utils.get_pdb_file("1CBS")

        assert result.endswith("1CBS.cif.gz")
        called_url = mock_get.call_args[0][0]
        assert called_url == "https://files.rcsb.org/download/1cbs.cif.gz"

    @patch("core.http_cache.requests.get")
    def test_uses_pdb_redo_url_when_flag_set(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = b"data"
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        result = pdb_utils.get_pdb_file("1CBS", pdb_redo=True)

        called_url = mock_get.call_args[0][0]
        assert called_url == "https://pdb-redo.eu/db/1cbs/1cbs_final.cif"
        assert result.endswith("1cbs_final.cif")

    def test_returns_cached_file_without_downloading(self, tmp_path):
        cached = os.path.join(str(tmp_path), "1CBS.cif.gz")
        with open(cached, "wb") as fh:
            fh.write(b"cached-content")

        with patch("core.http_cache.requests.get") as mock_get:
            result = pdb_utils.get_pdb_file("1CBS")
            mock_get.assert_not_called()
        assert result == cached

    @patch("core.http_cache.requests.get")
    def test_returns_empty_string_on_download_failure(self, mock_get, tmp_path):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        result = pdb_utils.get_pdb_file("9ZZZ")
        assert result == ""
