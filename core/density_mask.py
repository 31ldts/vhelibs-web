# -*- coding: utf-8 -*-
#
#   New module, added to move per-atom electron-density masking from the
#   browser (Mol* volume-representation .clip(), one isosurface per atom
#   — see the old buildMvsData() in app.js) to the server, where it can
#   be computed once with gemmi and cached to disk like every other
#   downloaded/derived file (see core.eds_utils.get_edm/get_EDS).
#
#   For a given (pdbid, region, radius, source), do the masking
#   once here — crop the source CCP4 map to the region's padded box and
#   zero out every voxel further than `radius` Å from every atom in that
#   region — and write the result out as its own small .ccp4 file. The
#   client then needs exactly ONE isosurface representation per region,
#   built from an already-masked volume, with no clip() at all.
#
import os
import json
import hashlib
import logging

import numpy as np
import gemmi

import core.pdb_utils as pdb_utils
import core.eds_utils as eds_utils
import core.http_cache as http_cache

logger = logging.getLogger(__name__)

# NOTE: this assumes core.pdb_redo_utils exposes a get_EDM(pdbid,
# use_cache=True) -> (filepath_or_None, sigma_or_None) function returning
# the full PDB-REDO CCP4 map, mirroring core.eds_utils.get_edm — that's
# what app.js's pdbRedoMapUrl()/checkPdbRedoMapAvailable() imply on the
# existing /api/edm?source=pdb_redo route. Written without sight of
# pdb_redo_utils.py itself, so please confirm the signature matches
# before wiring in the "source=pdb_redo" path below.
import core.pdb_redo_utils as pdb_redo_utils

# Quantization step (Angstrom) for the atom-mask radius. The 3D viewer's
# slider is continuous (0.1 Å steps), but caching a masked map for every
# float value users drag past would never hit the cache. Snapping to
# 0.25 Å keeps the number of cached files per region small (13 possible
# values across the slider's 0.5-3.5 Å range) while staying visually
# indistinguishable from a truly continuous radius.
RADIUS_STEP = 0.25
RADIUS_MIN = 0.25
RADIUS_MAX = 5.0

# Extra padding (Angstrom), beyond the atom-mask radius itself, kept
# around the region when cropping the map. A small margin avoids cutting
# the isosurface exactly at the crop boundary, which can otherwise leave
# a visible flat face where the voxel grid ends abruptly.
CROP_MARGIN = 1.5

# Sanity cap on how many atoms one masking call will honour. Unlike the
# old per-atom-representation approach, cost here is one native
# set_points_around() call per atom on an already-small cropped grid, so
# this is far cheaper than before — the cap exists to bound worst-case
# gemmi/query-string cost, not GPU draw calls.
MAX_MASK_ATOMS = 500

VALID_REGIONS = ("ligand", "binding_site", "residues_to_examine")


def quantize_radius(radius, step=RADIUS_STEP):
    """Snap a user-supplied atom-mask radius to a cache-friendly step.

    Args:
        radius (float): Requested radius in Angstrom.
        step (float, optional): Quantization step. Defaults to
            :data:`RADIUS_STEP`.

    Returns:
        float: ``radius`` rounded to the nearest multiple of ``step``,
        clamped to ``[RADIUS_MIN, RADIUS_MAX]``.
    """
    try:
        radius = float(radius)
    except (TypeError, ValueError):
        radius = RADIUS_MIN
    radius = max(RADIUS_MIN, min(RADIUS_MAX, radius))
    return round(round(radius / step) * step, 2)


def _region_cache_path(downloaddir, pdbid, region, radius, source, box, atoms):
    tag = "redo" if source == "pdb_redo" else "pdb"
    # A PDB entry can have several ligands/binding sites, each producing
    # its OWN box/atoms for the *same* region name (e.g. two different
    # "ligand" regions for two different bound ligands) — so the cache
    # key can't be just (pdbid, region, radius, source), or the second
    # ligand's request would silently return the first ligand's cached
    # (wrong) masked map. A short hash of the actual box+atom coordinates
    # disambiguates them while still hitting the cache on a genuine
    # repeat request for the very same ligand.
    digest_src = json.dumps(
        {"box": box, "atoms": [[round(c, 3) for c in a["center"]] for a in atoms]},
        sort_keys=True,
    )
    digest = hashlib.sha1(digest_src.encode("utf-8")).hexdigest()[:10]
    return os.path.join(downloaddir, f"{pdbid}_{region}_{tag}_{digest}_r{radius:.2f}.ccp4")


def _full_map_path(pdbid, source, use_cache=True):
    """Resolve the on-disk path of the full source map, downloading/
    caching it first via the same helpers get_edm()/get_EDM() already
    use, so this never duplicates that download logic.

    Returns:
        str or None: Path to the full CCP4 map, or None if unavailable.
    """
    if source == "pdb_redo":
        # core.pdb_redo_utils.get_EDM returns the map path directly, NOT
        # a (path, sigma) tuple like eds_utils.get_edm below — confirmed
        # from its real call site in app/routes.py's /api/edm route:
        #   mapfile = pdb_redo_utils.get_EDM(pdbid)
        # It also isn't called with a use_cache kwarg there, so we don't
        # assume it accepts one either.
        return pdb_redo_utils.get_EDM(pdbid)
    path, _sigma = eds_utils.get_edm(pdbid, use_cache=use_cache)
    return path


def get_masked_region_map(pdbid, region, box, atoms, radius=1.6,
                           source="pdb", use_cache=True):
    """Build (or reuse a cached) masked CCP4 map for one density region.

    Crops the full source map to ``box`` (padded by ``radius`` +
    :data:`CROP_MARGIN`) and zeroes every voxel further than ``radius``
    Angstrom from every atom in ``atoms``, so the returned map already
    shows exactly the segmented density a region should display — no
    further per-atom clipping is needed client-side.

    Args:
        pdbid (str): PDB identifier. Lower-cased before use.
        region (str): One of :data:`VALID_REGIONS`
            (``"ligand"``, ``"binding_site"``, ``"residues_to_examine"``).
            Only used to name the cache file; has no effect on the
            masking itself (that's fully determined by ``box``/``atoms``).
        box (dict): Padded bounding box, ``{"min": [x,y,z], "max": [x,y,z]}``
            — normally ``density_boxes[region]`` from the analysis result
            (see :func:`core.rsr_core.residues_bbox`).
        atoms (list): List of ``{"residue": ..., "center": [x,y,z]}`` dicts
            — normally ``density_atoms[region]`` from the analysis result
            (see :func:`core.rsr_core.residue_atom_centers`). Truncated to
            :data:`MAX_MASK_ATOMS` if longer.
        radius (float, optional): Atom-mask radius in Angstrom, before
            quantization. Defaults to ``1.6``.
        source (str, optional): ``"pdb"`` for the standard RCSB/EBI map,
            or ``"pdb_redo"`` for the PDB-REDO map. Defaults to ``"pdb"``.
        use_cache (bool, optional): Whether an existing cached masked map
            (or the underlying full map) may be reused. Defaults to
            ``True``.

    Returns:
        str or None: Path to the masked CCP4 map file, or ``None`` if no
        map could be produced (missing source map, empty box/atoms, or a
        gemmi/IO error — all logged, never raised).

    Raises:
        None
    """
    pdbid = pdbid.lower()
    if region not in VALID_REGIONS:
        logger.warning("Unknown density region %r", region)
        return None

    # Validate BEFORE touching the cache: otherwise an empty/invalid
    # request for a (pdbid, region) that was previously served correctly
    # for a *different* ligand's box/atoms could — depending on cache-key
    # details — return that unrelated cached file instead of None.
    if not box or not box.get("min") or not box.get("max") or not atoms:
        return None

    atoms = atoms[:MAX_MASK_ATOMS]
    radius = quantize_radius(radius)
    downloaddir = http_cache.entry_cache_dir(pdb_utils.CACHEDIR, pdbid)
    out_path = _region_cache_path(downloaddir, pdbid, region, radius, source, box, atoms)

    if use_cache and os.path.isfile(out_path) and os.path.getsize(out_path) > 0:
        return out_path

    full_map_path = _full_map_path(pdbid, source, use_cache=use_cache)
    if not full_map_path:
        logger.info("No source map available to mask for %s (%s)", pdbid, source)
        return None

    try:
        m = gemmi.read_ccp4_map(full_map_path)
        # Expand to a real-space box covering the region we care about
        # (maps for the asymmetric unit only would otherwise need
        # symmetry expansion first; setup(0.0) is a safe no-op if the
        # map already covers the full cell, which is the case for the
        # EBI entry-files .ccp4 downloaded by get_edm()).
        m.setup(0.0)

        # Build the atom mask at the FULL map's native grid dimensions,
        # BEFORE cropping — and multiply it in before calling set_extent.
        # Both steps matter and must happen in this order:
        #
        #  1. A gemmi.FloatGrid(nu, nv, nw) interprets nu/nv/nw as the
        #     number of sampling points across the *whole* unit cell.
        #     m.grid.nu/nv/nw means that BEFORE cropping too, so building
        #     the mask from those values here reproduces the source map's
        #     true voxel spacing. Building it from the POST-crop nu/nv/nw
        #     instead (nu shrinks after set_extent) would silently
        #     construct a much coarser grid spanning the same full cell,
        #     misaligning every atom position passed to
        #     set_points_around() and masking almost nothing.
        #  2. set_extent() (below) leaves the map unable to recompute its
        #     own header statistics afterwards (update_ccp4_header raises
        #     "run setup() first", and setup() re-expands the grid back
        #     to full size, undoing the crop) — so update_ccp4_header
        #     must be called while the grid is still "full", i.e. before
        #     set_extent, not after.
        mask = gemmi.FloatGrid(m.grid.nu, m.grid.nv, m.grid.nw)
        mask.set_unit_cell(m.grid.unit_cell)
        mask.spacegroup = m.grid.spacegroup
        for atom in atoms:
            x, y, z = atom["center"]
            # Native gemmi primitive for exactly this: mark every voxel
            # within `radius` of a Cartesian point. This is the
            # server-side equivalent of the per-atom "sphere" the old
            # client-side clip() built — the difference is it's applied
            # to a mask grid ONCE, not to N redundant isosurface meshes.
            mask.set_points_around(gemmi.Position(x, y, z), radius=radius, value=1.0)

        grid_arr = np.array(m.grid, copy=False)
        mask_arr = np.array(mask, copy=False)
        grid_arr *= mask_arr  # zero every voxel outside every atom's sphere, in place

        # Recompute header stats (mean/rms/min/max) for the masked data
        # while the map is still in its full, "set up" state — see note
        # above for why this can't be done after set_extent.
        m.update_ccp4_header(2, True)

        pad = radius + CROP_MARGIN
        cell = m.grid.unit_cell
        corners = [
            (box["min"][0] - pad, box["min"][1] - pad, box["min"][2] - pad),
            (box["max"][0] + pad, box["max"][1] + pad, box["max"][2] + pad),
        ]
        fcorners = [cell.fractionalize(gemmi.Position(*c)) for c in corners]
        fbox = gemmi.FractionalBox()
        fbox.minimum = gemmi.Fractional(
            min(fcorners[0].x, fcorners[1].x),
            min(fcorners[0].y, fcorners[1].y),
            min(fcorners[0].z, fcorners[1].z),
        )
        fbox.maximum = gemmi.Fractional(
            max(fcorners[0].x, fcorners[1].x),
            max(fcorners[0].y, fcorners[1].y),
            max(fcorners[0].z, fcorners[1].z),
        )
        # Crops the grid down to the (small) region box. This is what
        # keeps the masked file small — everything past this point is a
        # few thousand voxels, not the whole-cell map. Do NOT call
        # update_ccp4_header() or setup() again after this: both raise/
        # undo the crop (see the comment above the mask construction).
        m.set_extent(fbox)

        os.makedirs(downloaddir, exist_ok=True)
        m.write_ccp4_map(out_path)
        return out_path

    except Exception as exc:
        logger.error("Density masking failed for %s/%s r=%.2f (%s): %s",
                     pdbid, region, radius, source, exc)
        return None


def prefetch_default_masks(pdbid, density_boxes, density_atoms,
                            radius=1.6, source="pdb", executor=None):
    """Fire-and-forget precompute of the default-radius masked maps.

    Intended to be called right after an analysis result is built (see
    core.rsr_core), while the caller already has ``density_boxes``/
    ``density_atoms`` in hand, so that by the time someone actually opens
    the 3D viewer for this ligand the default view needs no gemmi work —
    only the on-demand endpoint has to do anything, and only for radii
    other than the default.

    This intentionally only prefetches the STANDARD (non PDB-REDO) map
    at the default radius: PDB-REDO maps come from a separate, slower
    external map-maker service and aren't always available, so eagerly
    fetching them for every result would slow down/spam that service for
    entries nobody ends up viewing in 3D. The much cheaper thing (an RCSB
    map lookup that's already cached from get_edm()/edm_exists() calls
    elsewhere in the analysis) is prefetched unconditionally.

    Args:
        pdbid (str): PDB identifier.
        density_boxes (dict): The ``density_boxes`` dict for one ligand
            result, as built by core.rsr_core (keys: ``"ligand"``,
            ``"binding_site"``, ``"residues_to_examine"``).
        density_atoms (dict): The matching ``density_atoms`` dict.
        radius (float, optional): Radius to prefetch. Defaults to 1.6,
            matching the 3D viewer's default slider value (see
            index.html's #atomRadiusSlider).
        source (str, optional): Defaults to ``"pdb"`` — see note above
            about why PDB-REDO isn't prefetched here.
        executor (concurrent.futures.Executor, optional): Executor to
            submit the work to. If ``None``, this runs synchronously
            (only useful for tests/CLI use — call sites in a Flask
            request path should always pass a shared background
            executor, or this will block the analysis job).

    Returns:
        None. Errors are caught and logged inside get_masked_region_map;
        this function never raises and never blocks the caller when an
        executor is given.
    """
    def _run():
        for region in VALID_REGIONS:
            box = (density_boxes or {}).get(region)
            atoms = (density_atoms or {}).get(region)
            if box and atoms:
                get_masked_region_map(pdbid, region, box, atoms,
                                       radius=radius, source=source, use_cache=True)

    if executor is None:
        _run()
    else:
        executor.submit(_run)
