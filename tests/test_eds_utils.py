# -*- coding: utf-8 -*-
"""
Unit tests for core.eds_utils.

Network calls (requests.get) and the shared CACHEDIR are mocked/redirected;
no real network or persistent filesystem access occurs.
"""
import os
from unittest.mock import MagicMock, patch

import pytest
import requests

import core.eds_utils as eds_utils
import core.pdb_utils as pdb_utils


@pytest.fixture(autouse=True)
def isolated_cache_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(pdb_utils, "CACHEDIR", str(tmp_path))
    yield tmp_path


VALIDATION_XML = b"""<?xml version="1.0"?>
<wwPDB-validation-information>
  <Entry PDBID="1cbs"/>
  <ModelledSubgroup resname="REA" chain="A" resnum="200" icode=""
                     rsr="0.12" rscc="0.95" owab="18.4" rsrz="0.5" avgoccu="1.0"/>
  <ModelledSubgroup resname="HOH" chain="A" resnum="301" icode=""
                     rsr="0.30" rscc="0.80" avgoccu="0.5"/>
</wwPDB-validation-information>
"""


# ---------------------------------------------------------------------------
# edm_box_url
# ---------------------------------------------------------------------------

class TestEdmBoxUrl:
    def test_builds_url_with_lowercased_pdbid(self):
        box = {"min": [1.0, 2.0, 3.0], "max": [4.0, 5.0, 6.0]}
        url = eds_utils.edm_box_url("1CBS", box)
        assert url.startswith("https://www.ebi.ac.uk/pdbe/densities/x-ray/1cbs/box/")
        assert "1.0,2.0,3.0/4.0,5.0,6.0" in url
        assert url.endswith("?detail=3")

    def test_custom_detail_level(self):
        box = {"min": [0, 0, 0], "max": [1, 1, 1]}
        url = eds_utils.edm_box_url("1cbs", box, detail=6)
        assert url.endswith("?detail=6")


# ---------------------------------------------------------------------------
# get_edm
# ---------------------------------------------------------------------------

class TestGetEdm:
    @patch("core.eds_utils.requests.get")
    def test_downloads_and_caches_map(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.iter_content = MagicMock(return_value=[b"CCP4", b"DATA"])
        mock_get.return_value = resp

        path, sigma = eds_utils.get_edm("1cbs")

        assert sigma == 1.0
        assert os.path.isfile(path)
        with open(path, "rb") as fh:
            assert fh.read() == b"CCP4DATA"

    @patch("core.eds_utils.requests.get")
    def test_uses_cache_when_present(self, mock_get, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        cache_file = os.path.join(cache_dir, "1cbs.ccp4")
        with open(cache_file, "wb") as fh:
            fh.write(b"cached")

        path, sigma = eds_utils.get_edm("1cbs")

        mock_get.assert_not_called()
        assert path == cache_file
        assert sigma == 1.0

    @patch("core.eds_utils.requests.get")
    def test_use_cache_false_forces_redownload(self, mock_get, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        with open(os.path.join(cache_dir, "1cbs.ccp4"), "wb") as fh:
            fh.write(b"old")

        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.iter_content = MagicMock(return_value=[b"new"])
        mock_get.return_value = resp

        path, sigma = eds_utils.get_edm("1cbs", use_cache=False)

        mock_get.assert_called_once()
        with open(path, "rb") as fh:
            assert fh.read() == b"new"

    @patch("core.eds_utils.requests.get")
    def test_404_returns_none_none(self, mock_get):
        resp = MagicMock()
        resp.status_code = 404
        mock_get.return_value = resp

        path, sigma = eds_utils.get_edm("9zzz")
        assert path is None
        assert sigma is None

    @patch("core.eds_utils.requests.get")
    def test_network_error_returns_none_none(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        path, sigma = eds_utils.get_edm("1cbs")
        assert path is None
        assert sigma is None


# ---------------------------------------------------------------------------
# get_EDS / get_validation_data
# ---------------------------------------------------------------------------

class TestGetEDS:
    def test_get_validation_data_is_alias_for_get_EDS(self):
        assert eds_utils.get_validation_data is eds_utils.get_EDS

    @patch("core.eds_utils.requests.get")
    def test_parses_validation_xml(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.content = VALIDATION_XML
        mock_get.return_value = resp

        pdbdict, edd_dict = eds_utils.get_EDS("1cbs")

        assert pdbdict == {"1cbs": True}
        assert "REA A  200" in edd_dict
        rea = edd_dict["REA A  200"]
        assert rea["RSR"] == 0.12
        assert rea["RSCC"] == 0.95
        assert rea["OWAB"] == 18.4
        assert rea["occupancy"] == 1.0

    @patch("core.eds_utils.requests.get")
    def test_missing_optional_attributes_use_defaults(self, mock_get):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.content = VALIDATION_XML
        mock_get.return_value = resp

        _, edd_dict = eds_utils.get_EDS("1cbs")
        hoh = edd_dict["HOH A  301"]
        # OWAB and RSRZ absent from that <ModelledSubgroup> -> defaults apply
        assert hoh["OWAB"] == 1000
        assert hoh["RSRZ"] == 9999
        assert hoh["occupancy"] == 0.5
        assert hoh["RSR"] == 0.30

    @patch("core.eds_utils.requests.get")
    def test_404_sets_pdbdict_false_and_empty_edd(self, mock_get):
        resp = MagicMock()
        resp.status_code = 404
        mock_get.return_value = resp

        pdbdict, edd_dict = eds_utils.get_EDS("9zzz")
        assert pdbdict == {"9zzz": False}
        assert edd_dict == {}

    @patch("core.eds_utils.requests.get")
    def test_network_error_stores_exception_string(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        pdbdict, edd_dict = eds_utils.get_EDS("1cbs")
        assert pdbdict["1cbs"] == "down"
        assert edd_dict == {}

    @patch("core.eds_utils.requests.get")
    def test_uses_cached_xml_file_without_refetching(self, mock_get, tmp_path):
        cache_dir = os.path.join(str(tmp_path), "1cbs")
        os.makedirs(cache_dir, exist_ok=True)
        with open(os.path.join(cache_dir, "1cbs_validation.xml"), "wb") as fh:
            fh.write(VALIDATION_XML)

        pdbdict, edd_dict = eds_utils.get_EDS("1cbs")

        mock_get.assert_not_called()
        assert pdbdict == {"1cbs": True}
        assert "REA A  200" in edd_dict

    @patch("core.eds_utils.requests.get")
    def test_malformed_xml_after_download_is_handled(self, mock_get):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.content = b"not xml at all <<<"
        mock_get.return_value = resp

        pdbdict, edd_dict = eds_utils.get_EDS("1cbs")
        assert edd_dict == {}
        assert isinstance(pdbdict["1cbs"], str)  # exception message stored
