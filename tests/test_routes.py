# -*- coding: utf-8 -*-
"""
Unit tests for app.routes.

Strategy:
  - A real Flask app (via app.create_app) + test_client() for HTTP-level
    behaviour, with CACHE_DIR redirected to tmp_path.
  - core.pdb_utils / core.eds_utils / core.pdb_redo_utils / the imported
    core.rsr_core.parse_binding_site are all mocked — no real network
    access ever happens.
  - Background jobs run on real threads (whatever _run_job actually uses
    internally — a single thread, a ThreadPoolExecutor, etc. — is treated
    as an implementation detail). Tests that need job completion poll
    routes._jobs via wait_for_job() instead of assuming synchronous
    execution, so they don't depend on that internal choice.
"""
import json
import os
import time
from unittest.mock import MagicMock, patch

import pytest

import app.routes as routes
import core.pdb_utils as pdb_utils


def wait_for_job(job_id, timeout=5.0, poll_interval=0.01):
    """Poll routes._jobs until the given job reaches a terminal status.

    Deliberately implementation-agnostic: it doesn't assume _run_job uses a
    single background thread vs. a thread pool internally, just that it
    eventually flips job["status"] to "done" (mirroring the real async
    contract of POST /api/analyse + GET /api/status/<id>).
    """
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = routes._jobs.get(job_id)
        if job and job.get("status") == "done":
            return job
        time.sleep(poll_interval)
    raise AssertionError(f"job {job_id} did not reach 'done' within {timeout}s")


@pytest.fixture
def client(tmp_path, monkeypatch):
    from app import create_app
    application = create_app(cache_dir=str(tmp_path))
    application.config["TESTING"] = True
    monkeypatch.setattr(pdb_utils, "CACHEDIR", str(tmp_path))
    with application.test_client() as c:
        yield c


@pytest.fixture(autouse=True)
def clear_jobs():
    """Job store is a module-level dict shared across the whole test
    session; make sure tests don't see each other's jobs."""
    routes._jobs.clear()
    yield
    routes._jobs.clear()


# ---------------------------------------------------------------------------
# UNIPROT_RE / _expand_ids
# ---------------------------------------------------------------------------

class TestUniprotRegex:
    @pytest.mark.parametrize("accession", ["P00734", "A0A0A0MRZ7", "Q9Y6K9"])
    def test_matches_valid_uniprot_accessions(self, accession):
        assert routes.UNIPROT_RE.match(accession)

    @pytest.mark.parametrize("pdbid", ["1CBS", "3dzu", "9XYZ"])
    def test_does_not_match_plain_pdb_ids(self, pdbid):
        assert not routes.UNIPROT_RE.match(pdbid)


class TestExpandIds:
    def test_plain_pdbids_pass_through_unchanged(self):
        pdbids, origin_map, unresolved = routes._expand_ids(["1cbs", "3DZU"])
        assert pdbids == ["1cbs", "3DZU"]
        assert origin_map == {"1cbs": None, "3dzu": None}
        assert unresolved == []

    @patch("core.pdb_utils.get_pdbids_for_uniprot")
    def test_uniprot_accession_expands_to_pdbids(self, mock_lookup):
        mock_lookup.return_value = ["1CBS", "3DZU"]
        pdbids, origin_map, unresolved = routes._expand_ids(["P00734"])
        assert pdbids == ["1CBS", "3DZU"]
        assert origin_map["1cbs"] == "P00734"
        assert origin_map["3dzu"] == "P00734"
        assert unresolved == []

    @patch("core.pdb_utils.get_pdbids_for_uniprot")
    def test_unresolved_uniprot_accession_is_reported(self, mock_lookup):
        mock_lookup.return_value = []
        pdbids, origin_map, unresolved = routes._expand_ids(["P99999"])
        assert pdbids == []
        assert unresolved == ["P99999"]

    @patch("core.pdb_utils.get_pdbids_for_uniprot")
    def test_deduplicates_preserving_first_appearance_order(self, mock_lookup):
        mock_lookup.return_value = ["1CBS"]
        pdbids, origin_map, unresolved = routes._expand_ids(["1cbs", "P00734", "1CBS"])
        assert pdbids == ["1cbs"]  # second "1CBS" (via uniprot) and third dropped as dup
        assert origin_map["1cbs"] is None  # first occurrence (plain ID) set the origin

    def test_mixed_tokens_preserve_order(self):
        pdbids, origin_map, unresolved = routes._expand_ids(["1cbs", "3dzu"])
        assert pdbids == ["1cbs", "3dzu"]


# ---------------------------------------------------------------------------
# GET /
# ---------------------------------------------------------------------------

class TestIndexRoute:
    def test_returns_200_and_html(self, client):
        resp = client.get("/")
        assert resp.status_code == 200
        assert b"<html" in resp.data.lower() or b"<!doctype" in resp.data.lower()


# ---------------------------------------------------------------------------
# POST /api/analyse
# ---------------------------------------------------------------------------

class TestAnalyseRoute:
    def test_no_pdbids_returns_400(self, client):
        resp = client.post("/api/analyse", json={})
        assert resp.status_code == 400
        assert "No PDB IDs" in resp.get_json()["error"]

    def test_blank_pdbids_string_returns_400(self, client):
        resp = client.post("/api/analyse", json={"pdbids": "   "})
        assert resp.status_code == 400

    @patch("core.rsr_core.parse_binding_site")
    def test_valid_pdbid_starts_job_and_completes(self, mock_parse, client):
        mock_parse.return_value = {"pdbid": "1cbs", "ligands": [], "rejected": {}, "struc_dict": {}}

        resp = client.post("/api/analyse", json={"pdbids": "1cbs"})

        assert resp.status_code == 200
        data = resp.get_json()
        assert "job_id" in data
        assert data["total"] == 1
        job = wait_for_job(data["job_id"])
        assert job["results"][0]["pdbid"] == "1cbs"
        assert job["results"][0]["uniprot"] is None

    @patch("core.rsr_core.parse_binding_site")
    def test_accepts_pdbids_as_a_list(self, mock_parse, client):
        mock_parse.return_value = {"pdbid": "1cbs"}
        resp = client.post("/api/analyse", json={"pdbids": ["1cbs", "3dzu"]})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["total"] == 2
        wait_for_job(data["job_id"])

    @patch("core.rsr_core.parse_binding_site")
    def test_pdbids_split_on_commas_newlines_and_whitespace(self, mock_parse, client):
        mock_parse.return_value = {"pdbid": "x"}
        resp = client.post("/api/analyse", json={"pdbids": "1cbs, 3dzu\n2xyz"})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["total"] == 3
        wait_for_job(data["job_id"])

    @patch("core.pdb_utils.get_pdbids_for_uniprot", return_value=[])
    def test_unresolved_uniprot_with_no_other_ids_returns_400_with_detail(self, mock_lookup, client):
        resp = client.post("/api/analyse", json={"pdbids": "P99999"})
        assert resp.status_code == 400
        assert "P99999" in resp.get_json()["error"]

    @patch("core.rsr_core.parse_binding_site")
    @patch("core.pdb_utils.get_pdbids_for_uniprot", return_value=[])
    def test_unresolved_uniprot_alongside_valid_id_returns_warning(
            self, mock_lookup, mock_parse, client):
        mock_parse.return_value = {"pdbid": "1cbs"}
        resp = client.post("/api/analyse", json={"pdbids": "1cbs P99999"})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["total"] == 1
        assert "warnings" in data
        assert "P99999" in data["warnings"][0]
        wait_for_job(data["job_id"])

    @patch("core.rsr_core.parse_binding_site")
    def test_default_config_values_used_when_not_provided(self, mock_parse, client):
        captured_cfg = {}

        def fake_parse(pdbid, cfg):
            captured_cfg["cfg"] = cfg
            return {"pdbid": pdbid}

        mock_parse.side_effect = fake_parse
        resp = client.post("/api/analyse", json={"pdbids": "1cbs"})
        wait_for_job(resp.get_json()["job_id"])

        cfg = captured_cfg["cfg"]
        assert cfg.rsr_upper == 0.4
        assert cfg.rsr_lower == 0.24
        assert cfg.tolerance == 2
        assert cfg.pdb_redo is False
        assert cfg.check_owab is False

    @patch("core.rsr_core.parse_binding_site")
    def test_custom_config_values_are_forwarded(self, mock_parse, client):
        captured_cfg = {}

        def fake_parse(pdbid, cfg):
            captured_cfg["cfg"] = cfg
            return {"pdbid": pdbid}

        mock_parse.side_effect = fake_parse
        resp = client.post("/api/analyse", json={
            "pdbids": "1cbs",
            "rsr_upper": 0.5,
            "tolerance": 1,
            "use_pdb_redo": True,
            "check_owab": True,
            "owab_max": 30,
        })
        wait_for_job(resp.get_json()["job_id"])

        cfg = captured_cfg["cfg"]
        assert cfg.rsr_upper == 0.5
        assert cfg.tolerance == 1
        assert cfg.pdb_redo is True
        assert cfg.check_owab is True
        assert cfg.owab_max == 30

    @patch("core.rsr_core.parse_binding_site")
    def test_empty_string_numeric_field_falls_back_to_default(self, mock_parse, client):
        captured_cfg = {}

        def fake_parse(pdbid, cfg):
            captured_cfg["cfg"] = cfg
            return {"pdbid": pdbid}

        mock_parse.side_effect = fake_parse
        resp = client.post("/api/analyse", json={"pdbids": "1cbs", "rsr_upper": ""})
        wait_for_job(resp.get_json()["job_id"])

        assert captured_cfg["cfg"].rsr_upper == 0.4  # default, not an error

    @patch("core.rsr_core.parse_binding_site")
    def test_per_entry_exception_does_not_abort_job(self, mock_parse, client):
        def fake_parse(pdbid, cfg):
            if pdbid == "bad1":
                raise RuntimeError("boom")
            return {"pdbid": pdbid}

        mock_parse.side_effect = fake_parse
        resp = client.post("/api/analyse", json={"pdbids": "1cbs bad1"})
        job = wait_for_job(resp.get_json()["job_id"])
        results = {r["pdbid"]: r for r in job["results"]}
        assert results["1cbs"] == {"pdbid": "1cbs", "uniprot": None}
        assert results["bad1"] == {"pdbid": "bad1", "error": "boom", "uniprot": None}


# ---------------------------------------------------------------------------
# GET /api/status/<job_id>
# ---------------------------------------------------------------------------

class TestStatusRoute:
    def test_unknown_job_returns_404(self, client):
        resp = client.get("/api/status/does-not-exist")
        assert resp.status_code == 404

    @patch("core.rsr_core.parse_binding_site")
    def test_known_job_returns_full_status(self, mock_parse, client):
        mock_parse.return_value = {"pdbid": "1cbs"}
        post_resp = client.post("/api/analyse", json={"pdbids": "1cbs"})
        job_id = post_resp.get_json()["job_id"]
        wait_for_job(job_id)

        resp = client.get(f"/api/status/{job_id}")
        assert resp.status_code == 200
        data = resp.get_json()
        assert data["status"] == "done"
        assert data["progress"] == 1
        assert data["total"] == 1
        assert data["results"][0]["pdbid"] == "1cbs"


# ---------------------------------------------------------------------------
# GET /api/edm/<pdbid>
# ---------------------------------------------------------------------------

class TestEdmRoute:
    @patch("core.eds_utils.get_edm")
    def test_default_source_serves_ccp4_file(self, mock_get_edm, client, tmp_path):
        mapfile = tmp_path / "1cbs.ccp4"
        mapfile.write_bytes(b"CCP4DATA")
        mock_get_edm.return_value = (str(mapfile), 1.0)

        resp = client.get("/api/edm/1cbs")
        assert resp.status_code == 200
        assert resp.data == b"CCP4DATA"

    @patch("core.eds_utils.get_edm", return_value=(None, None))
    def test_default_source_404_when_no_map(self, mock_get_edm, client):
        resp = client.get("/api/edm/9zzz")
        assert resp.status_code == 404

    @patch("core.pdb_redo_utils.get_EDM")
    def test_pdb_redo_source_serves_map_file(self, mock_get_edm, client, tmp_path):
        mapfile = tmp_path / "1cbs_final.map"
        mapfile.write_bytes(b"MAPDATA")
        mock_get_edm.return_value = str(mapfile)

        resp = client.get("/api/edm/1cbs?source=pdb_redo")
        assert resp.status_code == 200
        assert resp.data == b"MAPDATA"

    @patch("core.pdb_redo_utils.get_EDM", return_value=None)
    def test_pdb_redo_source_404_when_missing(self, mock_get_edm, client):
        resp = client.get("/api/edm/9zzz?source=pdb_redo")
        assert resp.status_code == 404


# ---------------------------------------------------------------------------
# GET /api/density-box/<pdbid>
# ---------------------------------------------------------------------------

class TestDensityBoxRoute:
    def test_missing_min_max_returns_400(self, client):
        resp = client.get("/api/density-box/1cbs")
        assert resp.status_code == 400

    def test_malformed_min_returns_400(self, client):
        resp = client.get("/api/density-box/1cbs?min=a,b,c&max=1,2,3")
        assert resp.status_code == 400

    @patch("app.routes.requests.get")
    def test_fetches_and_caches_binary_chunk(self, mock_get, client, tmp_path):
        resp_obj = MagicMock()
        resp_obj.raise_for_status = MagicMock()
        resp_obj.content = b"BINARYDATA"
        mock_get.return_value = resp_obj

        resp = client.get("/api/density-box/1cbs?min=0,0,0&max=1,1,1&detail=3")
        assert resp.status_code == 200
        assert resp.data == b"BINARYDATA"
        mock_get.assert_called_once()

        # Second identical request should be served from the on-disk cache.
        resp2 = client.get("/api/density-box/1cbs?min=0,0,0&max=1,1,1&detail=3")
        assert resp2.data == b"BINARYDATA"
        mock_get.assert_called_once()  # still just one call

    @patch("app.routes.requests.get")
    def test_upstream_failure_returns_502(self, mock_get, client):
        import requests
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        resp = client.get("/api/density-box/1cbs?min=0,0,0&max=1,1,1")
        assert resp.status_code == 502

    @patch("app.routes.requests.get")
    def test_default_detail_level_is_3(self, mock_get, client):
        resp_obj = MagicMock()
        resp_obj.raise_for_status = MagicMock()
        resp_obj.content = b"X"
        mock_get.return_value = resp_obj

        client.get("/api/density-box/1cbs?min=0,0,0&max=1,1,1")
        called_url = mock_get.call_args[0][0]
        assert "detail=3" in called_url
