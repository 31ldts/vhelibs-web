# -*- coding: utf-8 -*-
#
#   Copyright 2010-2024 Adrià Cereto Massagué
#   Added during refactor: shared HTTP download + disk-cache helpers.
#
#   The "check local cache dir, download with retries if missing, log and
#   swallow errors" pattern used to be duplicated almost verbatim across
#   pdb_utils.py, eds_utils.py and pdb_redo_utils.py. This module is the
#   single place that pattern now lives; the three call sites keep their
#   own module-specific caching decisions (what to cache, under which
#   filename, how to parse it) but delegate the actual network I/O here.
#
import os
import json
import logging

import requests

logger = logging.getLogger(__name__)


def entry_cache_dir(cachedir_root, pdbid):
    """Return the per-entry cache directory for a PDB id, creating it if needed.

    Args:
        cachedir_root (str): Root cache directory (e.g. ``pdb_utils.CACHEDIR``).
        pdbid (str): PDB identifier. Case is normalized to lower-case, since
            all three callers key their per-entry cache directories by
            lower-cased PDB id.

    Returns:
        str: Path to ``{cachedir_root}/{pdbid.lower()}``, guaranteed to exist.
    """
    path = os.path.join(cachedir_root, pdbid.lower())
    os.makedirs(path, exist_ok=True)
    return path


def download_file(url, dest_path, timeout=30, retries=3):
    """Download ``url`` to ``dest_path``, retrying up to ``retries`` times.

    Args:
        url (str): URL to download.
        dest_path (str): Local filesystem path to write the response body to.
        timeout (int or float, optional): Per-request timeout in seconds.
            Defaults to ``30``.
        retries (int, optional): Maximum number of attempts before giving up.
            Defaults to ``3``.

    Returns:
        bool: ``True`` if the download succeeded and the file was written,
        ``False`` if every attempt failed.

    Raises:
        None: Request and I/O errors are caught internally on each attempt
            and logged as warnings; the function returns ``False`` instead
            of propagating exceptions.
    """
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=timeout, verify=True)
            r.raise_for_status()
            with open(dest_path, "wb") as fh:
                fh.write(r.content)
            return True
        except Exception as exc:
            logger.warning("Download attempt %d/%d failed for %s: %s", attempt, retries, url, exc)
    return False


def download_if_missing(url, dest_path, timeout=30, retries=3):
    """Ensure ``dest_path`` exists locally, downloading ``url`` if it doesn't.

    A file is considered already cached if it exists and is non-empty; in
    that case no network request is made at all.

    Args:
        url (str): URL to download if ``dest_path`` is missing or empty.
        dest_path (str): Local filesystem path to check/write.
        timeout (int or float, optional): Per-request timeout in seconds.
            Defaults to ``30``.
        retries (int, optional): Maximum download attempts. Defaults to ``3``.

    Returns:
        bool: ``True`` if ``dest_path`` is usable afterwards (already
        cached, or freshly downloaded), ``False`` if it had to be
        downloaded and every attempt failed.

    Raises:
        None: See :func:`download_file`.
    """
    if os.path.isfile(dest_path) and os.path.getsize(dest_path) > 0:
        return True
    logger.info("Downloading %s", url)
    return download_file(url, dest_path, timeout=timeout, retries=retries)


def download_or_404(url, dest_path, timeout=60, stream=False, chunk_size=1 << 16):
    """Download ``url`` to ``dest_path`` in a single attempt, distinguishing "not found".

    Unlike :func:`download_file`, this makes exactly one request and treats
    an HTTP 404 as a distinct, non-error outcome (some datasets simply have
    no file for a given entry). Used by callers that want to tell "this
    entry has no such resource" apart from "the request failed".

    Args:
        url (str): URL to download.
        dest_path (str): Local filesystem path to write the response body
            to. When ``stream`` is ``True``, content is first written to
            ``dest_path + ".part"`` and atomically renamed on success, so a
            failed/interrupted download never leaves a partial file at the
            final path.
        timeout (int or float, optional): Request timeout in seconds.
            Defaults to ``60``.
        stream (bool, optional): If ``True``, stream the response body in
            chunks instead of buffering it fully in memory. Use for large
            files (e.g. full electron density maps). Defaults to ``False``.
        chunk_size (int, optional): Chunk size in bytes when ``stream`` is
            ``True``. Defaults to 64 KiB.

    Returns:
        bool or None: ``True`` on success, ``False`` if the request failed
        for a reason other than "not found", or ``None`` if the server
        responded ``404`` (resource does not exist for this entry).

    Raises:
        None: All request and I/O errors are caught internally and logged;
            the function returns ``False`` instead of propagating them.
    """
    try:
        r = requests.get(url, timeout=timeout, verify=True, stream=stream)
        if r.status_code == 404:
            return None
        r.raise_for_status()
        if stream:
            tmp_path = dest_path + ".part"
            with open(tmp_path, "wb") as fh:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        fh.write(chunk)
            os.replace(tmp_path, dest_path)
        else:
            with open(dest_path, "wb") as fh:
                fh.write(r.content)
        return True
    except Exception as exc:
        logger.error("Download error for %s: %s", url, exc)
        return False


def load_cached_json(cache_path):
    """Best-effort load of a JSON cache file.

    Args:
        cache_path (str): Path to a JSON file previously written by
            :func:`save_json` or :func:`fetch_json`.

    Returns:
        The parsed JSON value, or ``None`` if the file is missing, empty,
        or not valid JSON.

    Raises:
        None
    """
    if not (os.path.isfile(cache_path) and os.path.getsize(cache_path) > 0):
        return None
    try:
        with open(cache_path, "rt") as fh:
            return json.load(fh)
    except Exception:
        return None


def save_json(cache_path, data):
    """Best-effort write of ``data`` as JSON to ``cache_path``.

    Failures are swallowed (a cache write is an optimization, not something
    that should ever fail the caller's request).

    Args:
        cache_path (str): Destination path.
        data: JSON-serializable value to write.

    Returns:
        None

    Raises:
        None
    """
    try:
        with open(cache_path, "wt") as fh:
            json.dump(data, fh)
    except Exception:
        pass


def fetch_json(url, cache_path, timeout=30):
    """Return parsed JSON for ``url``, using ``cache_path`` as an on-disk cache.

    If a valid cached copy exists at ``cache_path`` it is loaded and
    returned directly, with no network request. Otherwise ``url`` is
    fetched with a plain GET, the parsed JSON is cached to ``cache_path``
    (best effort) and returned.

    Args:
        url (str): URL to GET and parse as JSON.
        cache_path (str): Local path used both to read an existing cache
            and to write a freshly-fetched one.
        timeout (int or float, optional): Request timeout in seconds.
            Defaults to ``30``.

    Returns:
        The parsed JSON value, or ``None`` if there is no usable cache and
        the fetch/parse failed.

    Raises:
        None: Request and parsing errors are caught internally and logged;
            the function returns ``None`` instead of propagating them.
    """
    cached = load_cached_json(cache_path)
    if cached is not None:
        logger.debug("Loading cached JSON: %s", cache_path)
        return cached

    logger.info("Fetching %s", url)
    try:
        r = requests.get(url, timeout=timeout, verify=True)
        r.raise_for_status()
        data = r.json()
    except Exception as exc:
        logger.error("Could not fetch %s: %s", url, exc)
        return None

    save_json(cache_path, data)
    return data
