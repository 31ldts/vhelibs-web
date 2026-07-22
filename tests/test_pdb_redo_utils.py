# -*- coding: utf-8 -*-
"""
Unit tests for core.pdb_redo_utils (post-refactor: download/cache logic
now delegates to core.http_cache; JSON parsing split into
_parse_ed_data / _extract_pdbredo_props).

Network access happens inside core.http_cache (requests.get); mocked
there. CACHEDIR is redirected to tmp_path.
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
# _parse_ed_data
# ---------------------------------------------------------------------------

class TestParseEdData:
    def test_parses_all_entries(self):
        result = pdb_redo_utils._parse_ed_data(ED_DATA_JSON)
        assert result["REA A  200"] == {"RSR": 0.15, "RSCC": 0.93}

    def test_missing_rsr_rsccs_use_defaults(self):
        result = pdb_redo_utils._parse_ed_data(ED_DATA_JSON)
        assert result["HOH A  301"] == {"RSR": 100.0, "RSCC": 0.0}

    def test_empty_list_returns_empty_dict(self):
        assert pdb_redo_utils._parse_ed_data([]) == {}

    def test_malformed_entry_raises_for_caller_to_catch(self):
        with pytest.raises(Exception):
            pdb_redo_utils._parse_ed_data([{"pdb": {"compID": "X"}}])  # missing strandID/seqNum


# ---------------------------------------------------------------------------
# get_ED_data
# ---------------------------------------------------------------------------

class TestGetEDData:
    @patch("core.http_cache.requests.get")
    def test_downloads_and_parses(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = json.dumps(ED_DATA_JSON).encode()
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        result = pdb_redo_utils.get_ED_data("1cbs")

        assert result["REA A  200"] == {"RSR": 0.15, "RSCC": 0.93}
        assert result["HOH A  301"] == {"RSR": 100.0, "RSCC": 0.0}
        called_url = mock_get.call_args[0][0]
        assert called_url == "https://pdb-redo.eu/db/1cbs/1cbs_final.json"

    @patch("core.http_cache.requests.get")
    def test_download_failure_returns_none(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        result = pdb_redo_utils.get_ED_data("1cbs")
        assert result is None

    def test_uses_cached_file_without_downloading(self, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        with open(os.path.join(cache_dir, "1cbs_final.json"), "w") as fh:
            json.dump(ED_DATA_JSON, fh)

        with patch("core.http_cache.requests.get") as mock_get:
            result = pdb_redo_utils.get_ED_data("1cbs")
            mock_get.assert_not_called()
        assert "REA A  200" in result

    def test_empty_cached_file_is_treated_as_missing_and_redownloaded(self, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        # Zero-byte file -> download_if_missing does not consider it cached.
        open(os.path.join(cache_dir, "1cbs_final.json"), "w").close()

        with patch("core.http_cache.requests.get") as mock_get:
            resp = MagicMock()
            resp.content = json.dumps(ED_DATA_JSON).encode()
            resp.raise_for_status = MagicMock()
            mock_get.return_value = resp

            result = pdb_redo_utils.get_ED_data("1cbs")
            mock_get.assert_called_once()
        assert "REA A  200" in result

    def test_corrupt_cached_file_returns_none(self, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        with open(os.path.join(cache_dir, "1cbs_final.json"), "w") as fh:
            fh.write("not valid json")

        # File is non-empty, so download_if_missing considers it "already
        # cached" and never re-fetches -> load_cached_json fails to parse
        # it -> get_ED_data returns None without a network call.
        with patch("core.http_cache.requests.get") as mock_get:
            result = pdb_redo_utils.get_ED_data("1cbs")
            mock_get.assert_not_called()
        assert result is None

    @patch("core.http_cache.requests.get")
    def test_json_shaped_wrong_for_parsing_returns_none(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = json.dumps([{"pdb": {"compID": "X"}}]).encode()  # missing keys
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        result = pdb_redo_utils.get_ED_data("1cbs")
        assert result is None


# ---------------------------------------------------------------------------
# _extract_pdbredo_props
# ---------------------------------------------------------------------------

class TestExtractPdbredoProps:
    def test_extracts_all_fields(self):
        row = pdb_redo_utils._extract_pdbredo_props(ALLDATA_JSON)
        assert row["experimentalTechnique"] == "X-RAY DIFFRACTION"
        assert row["rFree"] == 0.21
        assert row["rWork"] == 0.17
        assert row["refinementResolution"] == 1.9
        assert row["nreflections"] == 5000
        assert row["lengthOfUnitCellLatticeA"] == 40.0

    def test_missing_properties_key_raises_keyerror(self):
        with pytest.raises(KeyError):
            pdb_redo_utils._extract_pdbredo_props({"unexpected": "shape"})

    def test_missing_individual_props_use_defaults(self):
        row = pdb_redo_utils._extract_pdbredo_props({"properties": {}})
        assert row["rFree"] == 9999
        assert row["rWork"] == 9999
        assert row["refinementResolution"] == 0
        assert row["nreflections"] == 0


# ---------------------------------------------------------------------------
# get_pdbredo_data
# ---------------------------------------------------------------------------

class TestGetPdbredoData:
    @patch("core.http_cache.requests.get")
    def test_fetches_and_maps_fields(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = ALLDATA_JSON
        mock_get.return_value = resp

        result = pdb_redo_utils.get_pdbredo_data("1cbs")

        assert result["experimentalTechnique"] == "X-RAY DIFFRACTION"
        assert result["rFree"] == 0.21
        assert result["lengthOfUnitCellLatticeA"] == 40.0

    @patch("core.http_cache.requests.get")
    def test_uses_cache_on_second_call(self, mock_get):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = ALLDATA_JSON
        mock_get.return_value = resp

        pdb_redo_utils.get_pdbredo_data("1cbs")
        pdb_redo_utils.get_pdbredo_data("1cbs")

        assert mock_get.call_count == 1

    @patch("core.http_cache.requests.get")
    def test_network_failure_returns_none(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        result = pdb_redo_utils.get_pdbredo_data("1cbs")
        assert result is None

    @patch("core.http_cache.requests.get")
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

        with patch("core.http_cache.requests.get") as mock_get:
            resp = MagicMock()
            resp.raise_for_status = MagicMock()
            resp.json.return_value = ALLDATA_JSON
            mock_get.return_value = resp

            result = pdb_redo_utils.get_pdbredo_data("1cbs")

            mock_get.assert_called_once()
        assert result["rFree"] == 0.21

    @patch("core.http_cache.requests.get")
    def test_uses_longer_timeout_than_default(self, mock_get):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = ALLDATA_JSON
        mock_get.return_value = resp

        pdb_redo_utils.get_pdbredo_data("1cbs")
        assert mock_get.call_args.kwargs["timeout"] == pdb_redo_utils._DOWNLOAD_TIMEOUT
        assert pdb_redo_utils._DOWNLOAD_TIMEOUT == 60


# ---------------------------------------------------------------------------
# get_EDM
# ---------------------------------------------------------------------------

class TestGetEDM:
    """get_EDM() fetches PDB-REDO's *density map* — not to be confused
    with the *_final.json/data.json fetched by get_ED_data/get_pdbredo_data
    above. It used to (incorrectly) download *_final.mtz, the reflection-
    data file, which no CCP4-map reader can parse as a density grid — see
    the module-level note above PDB_REDO_MAP_MAKER_URL.
    """

    def test_returns_cached_file_path_without_downloading(self, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        cached = os.path.join(cache_dir, "1cbs_final.map")
        with open(cached, "wb") as fh:
            fh.write(b"mapdata")

        with patch("core.http_cache.requests.get") as mock_get:
            result = pdb_redo_utils.get_EDM("1cbs")
            mock_get.assert_not_called()
        assert result == cached

    @patch("core.http_cache.requests.get")
    def test_downloads_when_not_cached(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = b"mapdata"
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        result = pdb_redo_utils.get_EDM("1cbs")

        assert result.endswith("1cbs_final.map")
        called_url = mock_get.call_args[0][0]
        assert called_url == "https://pdb-redo.eu/map-maker/map?id=1cbs&stage=final&type=density"

    @patch("core.http_cache.requests.get")
    def test_pdbid_is_lowercased_in_the_url(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = b"mapdata"
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        pdb_redo_utils.get_EDM("1CBS")

        called_url = mock_get.call_args[0][0]
        assert "id=1cbs" in called_url

    @patch("core.http_cache.requests.get")
    def test_returns_none_on_failure(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        result = pdb_redo_utils.get_EDM("1cbs")
        assert result is None

    @patch("core.http_cache.requests.get")
    def test_use_cache_false_forces_redownload(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = b"mapdata"
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        pdb_redo_utils.get_EDM("1cbs", use_cache=True)
        pdb_redo_utils.get_EDM("1cbs", use_cache=False)

        assert mock_get.call_count == 2

    @patch("core.http_cache.requests.get")
    def test_uses_the_longer_map_download_timeout(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = b"mapdata"
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        pdb_redo_utils.get_EDM("1cbs")

        assert mock_get.call_args.kwargs["timeout"] == pdb_redo_utils._MAP_DOWNLOAD_TIMEOUT
        assert pdb_redo_utils._MAP_DOWNLOAD_TIMEOUT == 120
        # The map-maker computes the map on request, so it gets a longer
        # timeout than the other (static-file) PDB-REDO downloads.
        assert pdb_redo_utils._MAP_DOWNLOAD_TIMEOUT > pdb_redo_utils._DOWNLOAD_TIMEOUT
