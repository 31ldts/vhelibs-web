# -*- coding: utf-8 -*-
"""
Unit tests for core.pdb_utils.

All network access (via `requests`) and the module-level CACHEDIR are
mocked/redirected so tests never touch the real network or the real
filesystem outside of pytest's tmp_path.
"""
import json
import os
from unittest.mock import MagicMock, patch

import pytest
import requests

import core.pdb_utils as pdb_utils


@pytest.fixture(autouse=True)
def isolated_cache_dir(tmp_path, monkeypatch):
    """Redirect the module-level CACHEDIR to a throwaway tmp_path per test."""
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
# _download (private, but exercised indirectly and directly since it's the
# shared retry/download primitive used by every public function below)
# ---------------------------------------------------------------------------

class TestDownloadHelper:
    @patch("core.pdb_utils.requests.get")
    def test_download_success_writes_file(self, mock_get, tmp_path):
        mock_resp = MagicMock()
        mock_resp.content = b"filecontent"
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp

        dest = tmp_path / "out.bin"
        result = pdb_utils._download("http://example.com/f", str(dest))

        assert result is True
        assert dest.read_bytes() == b"filecontent"
        mock_get.assert_called_once_with("http://example.com/f", timeout=30, verify=True)

    @patch("core.pdb_utils.requests.get")
    def test_download_retries_then_fails(self, mock_get, tmp_path):
        mock_get.side_effect = requests.exceptions.ConnectionError("boom")
        dest = tmp_path / "out.bin"

        result = pdb_utils._download("http://example.com/f", str(dest), retries=3)

        assert result is False
        assert mock_get.call_count == 3
        assert not dest.exists()

    @patch("core.pdb_utils.requests.get")
    def test_download_succeeds_after_transient_failure(self, mock_get, tmp_path):
        ok_resp = MagicMock()
        ok_resp.content = b"data"
        ok_resp.raise_for_status = MagicMock()
        mock_get.side_effect = [requests.exceptions.Timeout("slow"), ok_resp]

        dest = tmp_path / "out.bin"
        result = pdb_utils._download("http://example.com/f", str(dest), retries=3)

        assert result is True
        assert mock_get.call_count == 2
        assert dest.read_bytes() == b"data"


# ---------------------------------------------------------------------------
# get_custom_report
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


class TestGetCustomReport:
    @patch("core.pdb_utils.requests.get")
    def test_fetches_and_parses_report(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = RCSB_ENTRY_JSON
        mock_get.return_value = resp

        result = pdb_utils.get_custom_report("1cbs")

        assert "1CBS" in result
        row = result["1CBS"]
        assert row["rFree"] == 0.218
        assert row["rWork"] == 0.175
        assert row["refinementResolution"] == 1.8
        assert row["nreflections"] == 1234
        assert row["experimentalTechnique"] == "X-RAY DIFFRACTION"
        assert row["unitCellAngleAlpha"] == 90.0
        assert row["lengthOfUnitCellLatticeA"] == 50.0

    @patch("core.pdb_utils.requests.get")
    def test_uses_disk_cache_on_second_call(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = RCSB_ENTRY_JSON
        mock_get.return_value = resp

        pdb_utils.get_custom_report("1cbs")
        pdb_utils.get_custom_report("1cbs")

        assert mock_get.call_count == 1  # second call served from cache file

    @patch("core.pdb_utils.requests.get")
    def test_network_failure_returns_empty_dict(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        result = pdb_utils.get_custom_report("1cbs")
        assert result == {}

    @patch("core.pdb_utils.requests.get")
    def test_falls_back_to_ls_d_res_low_when_high_missing(self, mock_get):
        data = json.loads(json.dumps(RCSB_ENTRY_JSON))  # deep copy
        del data["refine"][0]["ls_d_res_high"]
        data["refine"][0]["ls_d_res_low"] = 2.5
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = data
        mock_get.return_value = resp

        result = pdb_utils.get_custom_report("1cbs")
        assert result["1CBS"]["refinementResolution"] == 2.5

    @patch("core.pdb_utils.requests.get")
    def test_missing_refine_block_yields_placeholder_rfree(self, mock_get):
        data = {"rcsb_entry_info": {"experimental_method": "SOLUTION NMR"}}
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = data
        mock_get.return_value = resp

        result = pdb_utils.get_custom_report("2abc")
        row = result["2ABC"]
        assert row["rFree"] == 9999
        assert row["rWork"] == 9999
        assert row["refinementResolution"] == 0.0

    @patch("core.pdb_utils.requests.get")
    def test_malformed_json_is_handled_gracefully(self, mock_get):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {"refine": "not-a-list-of-dicts"}
        mock_get.return_value = resp

        # refine[0] on a string indexes the first character -> not a dict ->
        # .get(...) calls inside the try block raise -> caught -> {}
        result = pdb_utils.get_custom_report("3xyz")
        assert result == {}

    @patch("core.pdb_utils.requests.get")
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
                {"identifier": "1cbs"},  # duplicate, different case
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


# ---------------------------------------------------------------------------
# get_pdb_file
# ---------------------------------------------------------------------------

class TestGetPdbFile:
    @patch("core.pdb_utils._download")
    def test_builds_rcsb_url_and_filename_by_default(self, mock_download, tmp_path):
        def fake_download(url, dest_path, retries=3):
            with open(dest_path, "wb") as fh:
                fh.write(b"data")
            return True
        mock_download.side_effect = fake_download

        result = pdb_utils.get_pdb_file("1CBS")

        assert result.endswith("1CBS.cif.gz")
        called_url = mock_download.call_args[0][0]
        assert called_url == "https://files.rcsb.org/download/1cbs.cif.gz"

    @patch("core.pdb_utils._download")
    def test_uses_pdb_redo_url_when_flag_set(self, mock_download, tmp_path):
        def fake_download(url, dest_path, retries=3):
            with open(dest_path, "wb") as fh:
                fh.write(b"data")
            return True
        mock_download.side_effect = fake_download

        result = pdb_utils.get_pdb_file("1CBS", pdb_redo=True)

        called_url = mock_download.call_args[0][0]
        assert called_url == "https://pdb-redo.eu/db/1cbs/1cbs_final.cif"
        assert result.endswith("1cbs_final.cif")

    def test_returns_cached_file_without_downloading(self, tmp_path):
        cached = os.path.join(str(tmp_path), "1CBS.cif.gz")
        with open(cached, "wb") as fh:
            fh.write(b"cached-content")

        with patch("core.pdb_utils._download") as mock_download:
            result = pdb_utils.get_pdb_file("1CBS")
            mock_download.assert_not_called()
        assert result == cached

    @patch("core.pdb_utils._download", return_value=False)
    def test_returns_empty_string_on_download_failure(self, mock_download, tmp_path):
        result = pdb_utils.get_pdb_file("9ZZZ")
        assert result == ""
