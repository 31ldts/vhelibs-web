# -*- coding: utf-8 -*-
"""
Unit tests for core.http_cache.

This module is the new single home for the "check cache, download with
retries, distinguish 404s, cache JSON" logic shared by pdb_utils.py,
eds_utils.py and pdb_redo_utils.py. All network access (requests.get/post)
is mocked; filesystem access is redirected to pytest's tmp_path.
"""
import json
import os
from unittest.mock import MagicMock, patch

import pytest
import requests

import core.http_cache as http_cache


# ---------------------------------------------------------------------------
# entry_cache_dir
# ---------------------------------------------------------------------------

class TestEntryCacheDir:
    def test_builds_lowercased_subdirectory_and_creates_it(self, tmp_path):
        result = http_cache.entry_cache_dir(str(tmp_path), "1CBS")
        assert result == os.path.join(str(tmp_path), "1cbs")
        assert os.path.isdir(result)

    def test_idempotent_when_dir_already_exists(self, tmp_path):
        first = http_cache.entry_cache_dir(str(tmp_path), "1cbs")
        second = http_cache.entry_cache_dir(str(tmp_path), "1cbs")
        assert first == second
        assert os.path.isdir(second)

    def test_already_lowercase_pdbid(self, tmp_path):
        result = http_cache.entry_cache_dir(str(tmp_path), "3dzu")
        assert result.endswith("3dzu")


# ---------------------------------------------------------------------------
# download_file
# ---------------------------------------------------------------------------

class TestDownloadFile:
    @patch("core.http_cache.requests.get")
    def test_success_writes_response_body(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = b"filedata"
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        dest = tmp_path / "out.bin"
        result = http_cache.download_file("http://example.com/f", str(dest))

        assert result is True
        assert dest.read_bytes() == b"filedata"
        mock_get.assert_called_once_with("http://example.com/f", timeout=30, verify=True)

    @patch("core.http_cache.requests.get")
    def test_all_retries_exhausted_returns_false(self, mock_get, tmp_path):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        dest = tmp_path / "out.bin"

        result = http_cache.download_file("http://example.com/f", str(dest), retries=3)

        assert result is False
        assert mock_get.call_count == 3
        assert not dest.exists()

    @patch("core.http_cache.requests.get")
    def test_succeeds_after_transient_failure(self, mock_get, tmp_path):
        ok_resp = MagicMock()
        ok_resp.content = b"data"
        ok_resp.raise_for_status = MagicMock()
        mock_get.side_effect = [requests.exceptions.Timeout("slow"), ok_resp]

        dest = tmp_path / "out.bin"
        result = http_cache.download_file("http://example.com/f", str(dest), retries=3)

        assert result is True
        assert mock_get.call_count == 2
        assert dest.read_bytes() == b"data"

    @patch("core.http_cache.requests.get")
    def test_http_error_status_counts_as_failed_attempt(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.raise_for_status.side_effect = requests.exceptions.HTTPError("500")
        mock_get.return_value = resp

        dest = tmp_path / "out.bin"
        result = http_cache.download_file("http://example.com/f", str(dest), retries=2)

        assert result is False
        assert mock_get.call_count == 2

    @patch("core.http_cache.requests.get")
    def test_custom_timeout_is_forwarded(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = b"x"
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        http_cache.download_file("http://example.com/f", str(tmp_path / "o"), timeout=99)
        assert mock_get.call_args.kwargs["timeout"] == 99


# ---------------------------------------------------------------------------
# download_if_missing
# ---------------------------------------------------------------------------

class TestDownloadIfMissing:
    def test_existing_nonempty_file_skips_network(self, tmp_path):
        dest = tmp_path / "cached.bin"
        dest.write_bytes(b"already here")

        with patch("core.http_cache.requests.get") as mock_get:
            result = http_cache.download_if_missing("http://example.com/f", str(dest))
            mock_get.assert_not_called()
        assert result is True

    def test_existing_empty_file_triggers_download(self, tmp_path):
        dest = tmp_path / "empty.bin"
        dest.write_bytes(b"")

        with patch("core.http_cache.requests.get") as mock_get:
            resp = MagicMock()
            resp.content = b"real content"
            resp.raise_for_status = MagicMock()
            mock_get.return_value = resp

            result = http_cache.download_if_missing("http://example.com/f", str(dest))
        assert result is True
        assert dest.read_bytes() == b"real content"

    @patch("core.http_cache.requests.get")
    def test_missing_file_downloads_and_returns_true(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.content = b"data"
        resp.raise_for_status = MagicMock()
        mock_get.return_value = resp

        dest = tmp_path / "new.bin"
        result = http_cache.download_if_missing("http://example.com/f", str(dest))
        assert result is True
        assert dest.read_bytes() == b"data"

    @patch("core.http_cache.requests.get")
    def test_download_failure_returns_false(self, mock_get, tmp_path):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        dest = tmp_path / "new.bin"
        result = http_cache.download_if_missing("http://example.com/f", str(dest), retries=2)
        assert result is False
        assert not dest.exists()


# ---------------------------------------------------------------------------
# download_or_404
# ---------------------------------------------------------------------------

class TestDownloadOr404:
    @patch("core.http_cache.requests.get")
    def test_non_streamed_success_returns_true(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.content = b"payload"
        mock_get.return_value = resp

        dest = tmp_path / "f.xml"
        result = http_cache.download_or_404("http://x", str(dest))

        assert result is True
        assert dest.read_bytes() == b"payload"

    @patch("core.http_cache.requests.get")
    def test_404_returns_none_and_writes_nothing(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.status_code = 404
        mock_get.return_value = resp

        dest = tmp_path / "f.xml"
        result = http_cache.download_or_404("http://x", str(dest))

        assert result is None
        assert not dest.exists()

    @patch("core.http_cache.requests.get")
    def test_network_error_returns_false(self, mock_get, tmp_path):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        dest = tmp_path / "f.xml"
        result = http_cache.download_or_404("http://x", str(dest))
        assert result is False

    @patch("core.http_cache.requests.get")
    def test_streamed_download_writes_via_temp_file_and_renames(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.iter_content = MagicMock(return_value=[b"chunk1", b"chunk2"])
        mock_get.return_value = resp

        dest = tmp_path / "f.ccp4"
        result = http_cache.download_or_404("http://x", str(dest), stream=True)

        assert result is True
        assert dest.read_bytes() == b"chunk1chunk2"
        assert not (tmp_path / "f.ccp4.part").exists()  # renamed away

    @patch("core.http_cache.requests.get")
    def test_streamed_download_skips_empty_chunks(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.iter_content = MagicMock(return_value=[b"a", b"", b"b"])
        mock_get.return_value = resp

        dest = tmp_path / "f.ccp4"
        http_cache.download_or_404("http://x", str(dest), stream=True)
        assert dest.read_bytes() == b"ab"

    @patch("core.http_cache.requests.get")
    def test_stream_true_forwards_to_requests(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.iter_content = MagicMock(return_value=[b"x"])
        mock_get.return_value = resp

        http_cache.download_or_404("http://x", str(tmp_path / "f"), stream=True)
        assert mock_get.call_args.kwargs["stream"] is True

    @patch("core.http_cache.requests.get")
    def test_http_error_status_other_than_404_returns_false(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.status_code = 500
        resp.raise_for_status.side_effect = requests.exceptions.HTTPError("500")
        mock_get.return_value = resp

        result = http_cache.download_or_404("http://x", str(tmp_path / "f"))
        assert result is False


# ---------------------------------------------------------------------------
# load_cached_json / save_json
# ---------------------------------------------------------------------------

class TestLoadCachedJson:
    def test_missing_file_returns_none(self, tmp_path):
        assert http_cache.load_cached_json(str(tmp_path / "nope.json")) is None

    def test_empty_file_returns_none(self, tmp_path):
        f = tmp_path / "empty.json"
        f.write_text("")
        assert http_cache.load_cached_json(str(f)) is None

    def test_valid_json_is_parsed(self, tmp_path):
        f = tmp_path / "ok.json"
        f.write_text(json.dumps({"a": 1}))
        assert http_cache.load_cached_json(str(f)) == {"a": 1}

    def test_corrupt_json_returns_none(self, tmp_path):
        f = tmp_path / "bad.json"
        f.write_text("{not valid")
        assert http_cache.load_cached_json(str(f)) is None

    def test_list_payload_round_trips(self, tmp_path):
        f = tmp_path / "list.json"
        f.write_text(json.dumps([1, 2, 3]))
        assert http_cache.load_cached_json(str(f)) == [1, 2, 3]


class TestSaveJson:
    def test_writes_serializable_data(self, tmp_path):
        f = tmp_path / "out.json"
        http_cache.save_json(str(f), {"x": [1, 2]})
        assert json.loads(f.read_text()) == {"x": [1, 2]}

    def test_swallows_write_errors(self, tmp_path):
        # Directory as destination path -> open() raises IsADirectoryError,
        # which save_json must swallow rather than propagate.
        bad_path = tmp_path  # a directory, not a file
        http_cache.save_json(str(bad_path), {"a": 1})  # should not raise

    def test_swallows_non_serializable_data(self, tmp_path):
        f = tmp_path / "out.json"
        http_cache.save_json(str(f), {"bad": object()})  # should not raise
        # File may exist but be incomplete/absent — either is acceptable
        # under the "best effort, never raises" contract.


# ---------------------------------------------------------------------------
# fetch_json
# ---------------------------------------------------------------------------

class TestFetchJson:
    @patch("core.http_cache.requests.get")
    def test_fetches_and_caches_on_first_call(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {"hello": "world"}
        mock_get.return_value = resp

        cache_path = str(tmp_path / "cache.json")
        result = http_cache.fetch_json("http://x", cache_path)

        assert result == {"hello": "world"}
        assert json.loads(open(cache_path).read()) == {"hello": "world"}
        mock_get.assert_called_once()

    @patch("core.http_cache.requests.get")
    def test_second_call_uses_cache_without_network(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {"hello": "world"}
        mock_get.return_value = resp

        cache_path = str(tmp_path / "cache.json")
        http_cache.fetch_json("http://x", cache_path)
        result2 = http_cache.fetch_json("http://x", cache_path)

        assert result2 == {"hello": "world"}
        assert mock_get.call_count == 1

    @patch("core.http_cache.requests.get")
    def test_network_failure_with_no_cache_returns_none(self, mock_get, tmp_path):
        mock_get.side_effect = requests.exceptions.ConnectionError("down")
        result = http_cache.fetch_json("http://x", str(tmp_path / "cache.json"))
        assert result is None

    @patch("core.http_cache.requests.get")
    def test_corrupt_cache_file_triggers_refetch(self, mock_get, tmp_path):
        cache_path = tmp_path / "cache.json"
        cache_path.write_text("{not valid json")

        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {"fresh": True}
        mock_get.return_value = resp

        result = http_cache.fetch_json("http://x", str(cache_path))
        assert result == {"fresh": True}
        mock_get.assert_called_once()

    @patch("core.http_cache.requests.get")
    def test_custom_timeout_is_forwarded(self, mock_get, tmp_path):
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.json.return_value = {}
        mock_get.return_value = resp

        http_cache.fetch_json("http://x", str(tmp_path / "c.json"), timeout=15)
        assert mock_get.call_args.kwargs["timeout"] == 15
