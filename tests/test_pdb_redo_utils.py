# -*- coding: utf-8 -*-
"""
Unit tests for core.pdb_redo_utils.

Network calls (requests.get) and the shared CACHEDIR are mocked/redirected;
no real network or persistent filesystem access occurs.
"""
import json
import os
from unittest.mock import MagicMock, patch

import pytest
import requests

import core.pdb_redo_utils as pdb_redo_utils
import core.pdb_utils as pdb_utils


@pytest.fixture(autouse=True)
def isolated_cache_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(pdb_utils, "CACHEDIR", str(tmp_path))
    yield tmp_path


ED_DATA_JSON = [
    {
        "pdb": {"compID": "REA", "strandID": "A", "seqNum": 200},
        "RSR": 0.15,
        "RSCCS": 0.93,
    },
    {
        "pdb": {"compID": "HOH", "strandID": "A", "seqNum": 301},
        # RSR/RSCCS missing -> defaults apply
    },
]

ALLDATA_JSON = {
    "properties": {
        "EXPTYP": "X-RAY DIFFRACTION",
        "RFFIN": 0.21,
        "RFIN": 0.17,
        "RESOLUTION": 1.9,
        "ALPHA": 90.0, "BETA": 90.0, "GAMMA": 90.0,
        "AAXIS": 40.0, "BAXIS": 50.0, "CAXIS": 60.0,
        "NREFCNT": 5000,
    }
}


# ---------------------------------------------------------------------------
# _download (shared private helper, mirrors pdb_utils._download)
# ---------------------------------------------------------------------------

class TestDownloadHelper:
    @patch("core.pdb_redo_utils.requests.get")
    def test_download_success(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = b"content"
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        dest = tmp_path / "f.json"
        assert pdb_redo_utils._download("http://x", str(dest)) is True
        assert dest.read_bytes() == b"content"

    @patch("core.pdb_redo_utils.requests.get")
    def test_download_all_retries_fail(self, mock_get, tmp_path):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        dest = tmp_path / "f.json"
        assert pdb_redo_utils._download("http://x", str(dest), retries=2) is False
        assert mock_get.call_count == 2
        assert not dest.exists()


# ---------------------------------------------------------------------------
# get_ED_data
# ---------------------------------------------------------------------------

class TestGetEDData:
    @patch("core.pdb_redo_utils._download")
    def test_downloads_and_parses(self, mock_download, tmp_path):
        def fake_download(url, dest_path, retries=3):
            with open(dest_path, "w") as fh:
                json.dump(ED_DATA_JSON, fh)
            return True
        mock_download.side_effect = fake_download

        result = pdb_redo_utils.get_ED_data("1cbs")

        assert result["REA A  200"] == {"RSR": 0.15, "RSCC": 0.93}
        # Defaults applied when RSR/RSCCS absent.
        assert result["HOH A  301"] == {"RSR": 100.0, "RSCC": 0.0}
        called_url = mock_download.call_args[0][0]
        assert called_url == "https://pdb-redo.eu/db/1cbs/1cbs_final.json"

    @patch("core.pdb_redo_utils._download", return_value=False)
    def test_download_failure_returns_none(self, mock_download):
        result = pdb_redo_utils.get_ED_data("1cbs")
        assert result is None

    def test_uses_cached_file_without_downloading(self, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        with open(os.path.join(cache_dir, "1cbs_final.json"), "w") as fh:
            json.dump(ED_DATA_JSON, fh)

        with patch("core.pdb_redo_utils._download") as mock_download:
            result = pdb_redo_utils.get_ED_data("1cbs")
            mock_download.assert_not_called()
        assert "REA A  200" in result

    def test_corrupt_cached_file_returns_none(self, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        with open(os.path.join(cache_dir, "1cbs_final.json"), "w") as fh:
            fh.write("not valid json")

        result = pdb_redo_utils.get_ED_data("1cbs")
        assert result is None


# ---------------------------------------------------------------------------
# get_pdbredo_data
# ---------------------------------------------------------------------------

class TestGetPdbredoData:
    @patch("core.pdb_redo_utils.requests.get")
    def test_fetches_and_maps_fields(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = ALLDATA_JSON
        mock_get.return_value = resp

        result = pdb_redo_utils.get_pdbredo_data("1cbs")

        assert result["experimentalTechnique"] == "X-RAY DIFFRACTION"
        assert result["rFree"] == 0.21
        assert result["rWork"] == 0.17
        assert result["refinementResolution"] == 1.9
        assert result["nreflections"] == 5000
        assert result["lengthOfUnitCellLatticeA"] == 40.0

    @patch("core.pdb_redo_utils.requests.get")
    def test_uses_cache_on_second_call(self, mock_get):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = ALLDATA_JSON
        mock_get.return_value = resp

        pdb_redo_utils.get_pdbredo_data("1cbs")
        pdb_redo_utils.get_pdbredo_data("1cbs")

        assert mock_get.call_count == 1

    @patch("core.pdb_redo_utils.requests.get")
    def test_network_failure_returns_none(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        result = pdb_redo_utils.get_pdbredo_data("1cbs")
        assert result is None

    @patch("core.pdb_redo_utils.requests.get")
    def test_missing_properties_key_returns_none(self, mock_get):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {"unexpected": "shape"}
        mock_get.return_value = resp

        result = pdb_redo_utils.get_pdbredo_data("1cbs")
        assert result is None

    def test_corrupt_cache_file_triggers_refetch(self, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        with open(os.path.join(cache_dir, "data.json"), "w") as fh:
            fh.write("{not json")

        with patch("core.pdb_redo_utils.requests.get") as mock_get:
            resp = MagicMock()
            resp.raise_for_status = MagicMock()
            resp.json.return_value = ALLDATA_JSON
            mock_get.return_value = resp

            result = pdb_redo_utils.get_pdbredo_data("1cbs")

            mock_get.assert_called_once()
        assert result["rFree"] == 0.21


# ---------------------------------------------------------------------------
# get_EDM
# ---------------------------------------------------------------------------

class TestGetEDM:
    def test_returns_cached_file_path_without_downloading(self, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        cached = os.path.join(cache_dir, "1cbs_final.mtz")
        with open(cached, "wb") as fh:
            fh.write(b"mtzdata")

        with patch("core.pdb_redo_utils._download") as mock_download:
            result = pdb_redo_utils.get_EDM("1cbs")
            mock_download.assert_not_called()
        assert result == cached

    @patch("core.pdb_redo_utils._download")
    def test_downloads_when_not_cached(self, mock_download, tmp_path):
        def fake_download(url, dest_path, retries=3):
            with open(dest_path, "wb") as fh:
                fh.write(b"mtzdata")
            return True
        mock_download.side_effect = fake_download

        result = pdb_redo_utils.get_EDM("1cbs")

        assert result.endswith("1cbs_final.mtz")
        called_url = mock_download.call_args[0][0]
        assert called_url == "https://pdb-redo.eu/db/1cbs/1cbs_final.mtz"

    @patch("core.pdb_redo_utils._download", return_value=False)
    def test_returns_none_on_failure(self, mock_download):
        result = pdb_redo_utils.get_EDM("1cbs")
        assert result is None
