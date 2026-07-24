# VHELIBS Web

**Validation Helper for LIgands and Binding Sites** — a web-based rewrite of [VHELIBS](https://doi.org/10.1186/1758-2946-5-36).

VHELIBS Web analyses protein–ligand complexes from PDB structures. For every ligand it finds, it
evaluates the ligand itself and its surrounding binding site against real-space and refinement
quality metrics, then classifies each as **Good**, **Dubious**, or **Bad** based on thresholds you
control. Results can be filtered interactively and inspected residue-by-residue in an integrated
3D viewer, together with the experimental electron density around each region.

Where the original VHELIBS was a Java/Jython + Cython desktop application, this version is a pure
Python (Flask) backend with a lightweight JavaScript frontend — no compilation step, no desktop
dependencies. `pip install`, run, and it opens in your browser.

---

## Table of contents

- [Why VHELIBS](#why-vhelibs)
- [Features](#features)
- [Quick start](#quick-start)
- [Using the app](#using-the-app)
  - [Analysis tab](#analysis-tab)
    - [Editing the ligand blacklist](#editing-the-ligand-blacklist)
    - [Clearing the disk cache](#clearing-the-disk-cache)
  - [Results tab](#results-tab)
    - [Exporting results](#exporting-results)
  - [3D Viewer tab](#3d-viewer-tab)
    - [Manual quality review](#manual-quality-review)
- [How an analysis works](#how-an-analysis-works)
- [Classification criteria](#classification-criteria)
- [Project structure](#project-structure)
- [REST API](#rest-api)
- [Scripting the API / example script](#scripting-the-api--example-script)
- [Data sources](#data-sources)
- [Migration notes](#migration-notes-vs-the-original-vhelibs)
- [Citation](#citation)

---

## Why VHELIBS

X-ray structures deposited in the PDB vary widely in local quality. A ligand can sit in a
beautifully refined pocket — or in electron density that barely supports its pose. Relying on the
PDB ID or the structure's overall resolution alone isn't enough: two ligands in the *same*
structure can have very different real-space support. VHELIBS automates the residue-level checks
a crystallographer would normally do by eye (real-space R-factor, real-space correlation
coefficient, occupancy, refinement statistics) and turns them into a simple, actionable
**Good / Dubious / Bad** verdict per ligand and per binding site, so you can quickly triage which
structures are safe to use for docking, pharmacophore modelling, or structure-based design.

## Features

- **Batch analysis** of PDB IDs and/or UniProt IDs (each UniProt ID is
  auto-expanded, via the RCSB Search API, to every PDB structure that references it).
- **Configurable quality thresholds** with sensible crystallographic defaults (RSR, RSCC,
  occupancy, R-free, tolerance, binding-site distance cutoff).
- **Optional advanced checks**, off by default: OWAB, structure-wide resolution, R-diff
  (|R-free − R-work|), and DPI (diffraction-component precision index).
- **PDB-REDO integration** — analyse the re-refined structure and statistics instead of the
  original RCSB/PDBe deposition.
- **Editable ligand blacklist** — browse and search the built-in non-ligand/metal tables, uncheck
  entries you don't want treated as blacklisted, add custom codes, or replace the whole list from
  an uploaded file, all without touching the shared defaults.
- **Interactive results view**: every ligand and binding site is scored independently, and the
  Results tab lets you toggle Good/Dubious/Bad on each axis (plus a separate toggle for structures
  that couldn't be analysed) to filter down to exactly the combination you care about — e.g. every
  *Bad* ligand sitting in a *Dubious* binding site.
- **Export results to Excel (.xlsx)** — a two-sheet workbook with the parameters a run was
  submitted with and a per-ligand results table.
- **Integrated Mol\* 3D viewer** with independently toggleable protein/ligand/binding-site
  layers, per-region 2Fo-Fc electron density overlay, and live-adjustable contour level
  (isovalue) and atom-mask radius.
- **Manual quality review** — override the computed Good/Dubious/Bad call for any ligand, binding
  site, or individual "component to examine" directly from the 3D Viewer, and write the
  correction back into the Results tab.
- **Disk caching** of every downloaded structure/statistics file, so re-running an analysis (or
  re-opening the 3D viewer) doesn't re-hit external APIs. Nothing is removed from it
  automatically — the Analysis tab's **Clear disk cache** button deletes it on demand.
- **Scriptable REST API** — every action the web UI performs is a plain HTTP call (see
  [REST API](#rest-api)), so batch runs can be driven headlessly without a browser; see
  [Scripting the API / example script](#scripting-the-api--example-script) for a ready-to-run
  example.

## Quick start

> **Tip:** it's good practice to install dependencies inside a virtual environment rather than system-wide. Create and activate one first if you'd like:
>
> **macOS / Linux**
> ```bash
> python3 -m venv venv
> source venv/bin/activate
> ```
>
> **Windows (PowerShell)**
> ```powershell
> python -m venv venv
> venv\Scripts\Activate.ps1
> ```

```bash
pip install -r requirements.txt
python run.py
```

Then open <http://localhost:8000> in your browser.

### Command-line options

```
python run.py [--host HOST] [--port PORT] [--cache-dir PATH] [--no-browser] [--browser NAME] [--debug]
```

| Flag | Default | Description |
|------|---------|-------------|
| `--host` | `127.0.0.1` | Bind address |
| `--port` | `8000` | TCP port |
| `--cache-dir` | `~/.cache/vhelibs` | Directory for downloaded PDB/EDS files |
| `--no-browser` | off | Suppress automatic browser launch |
| `--browser` | system default | Browser to launch (e.g. `firefox`, `google-chrome`, `chromium`, `safari`, `windows-default`). Falls back to the system default if the name isn't recognized — see Python's [`webbrowser`](https://docs.python.org/3/library/webbrowser.html#webbrowser.get) docs for the names registered on your platform |
| `--debug` | off | Enable Flask debug/reloader (dev only) |

### Requirements

- Python 3.8+
- Packages: `flask`, `requests`, `gemmi`, `numpy` (see `requirements.txt`)
- A browser with WebGL support (for the 3D viewer, powered by [Mol\*](https://molstar.org))

## Using the app

The app is a single-page interface with four tabs: **Analysis**, **Results**, **3D Viewer**, and
**About**.

### Analysis tab

1. Enter one or more **PDB IDs and/or UniProt IDs** in the text box, separated by commas,
   whitespace, or newlines (e.g. `1cbs, 3dzu, 4hhb, P00734`) — or load them from a `.txt`/`.csv`
   file. UniProt IDs are automatically resolved to every associated PDB structure.
2. Optionally check **Use PDB-REDO structures** to analyse PDB-REDO's re-refined coordinates and
   statistics instead of the standard RCSB/PDBe deposition.
3. Adjust the **quality thresholds** (RSR, RSCC, R-free, occupancy, tolerance, binding-site
   distance) if the defaults don't suit your use case, and optionally expand **Advanced options**
   to enable OWAB / resolution / R-diff / DPI checks.
4. Optionally customize the **Ligand Blacklist** — see below.
5. Click **Analyse**. A progress bar tracks the job as each structure is downloaded and scored;
   when it finishes you're taken straight to the Results tab.

#### Editing the ligand blacklist

The **Ligand Blacklist** card lists every built-in entry across two columns — non-ligand
blacklist (solvents, buffers, crystallisation additives, …) and metals/ions — that VHELIBS
normally excludes from scoring rather than treating as ligands. You can:

- **Search/filter** the list by code or name, and **uncheck** any entry you don't want excluded
  for this analysis (it will then be scored as a normal ligand instead).
- Use **Select all** / **Select none** to bulk-toggle every visible entry, or **Restore defaults**
  to undo all changes.
- **Add a custom entry** (code, optional description, and category) under *Add a custom entry*.
- **Replace the whole list from a file** — accepts a plain list (one code per line, optionally
  `CODE,Description`) or a file previously exported with `[Blacklist]`/`[Non-propagating]` section
  headers. The app shows a preview of how many entries it parsed before you commit; applying it
  **replaces** the built-in defaults entirely rather than adding to them.

None of this modifies VHELIBS' shared built-in tables — it's scoped to the browser session and
sent along with the next **Analyse** request only.

#### Clearing the disk cache

Every structure, validation report, and density map VHELIBS downloads is cached on the server's
disk, and nothing ever removes it automatically otherwise. The
**Clear disk cache** button under *Use cached downloads* deletes the whole
cache on demand (after a confirmation prompt, since it can't be undone): use it to reclaim disk
space. It
reports how many files were removed and how much space was freed.

### Results tab

Each analysed structure appears as a collapsible card showing its R-free/R-work/resolution and
the number of ligands found; click a card to expand it and see every ligand, its binding-site
residues, low-occupancy warnings, and any residues excluded from scoring due to missing
validation data (listed under "rejected"). Structures that couldn't be analysed (e.g. unavailable
validation data) appear as separate error cards.

At the top of the tab, two rows of **toggle filters** — one for Ligand quality, one for Binding
site quality (Good / Dubious / Bad), each showing how many ligands currently have that
classification — let you narrow the view to exactly the combination you want (e.g. Ligand = Bad
+ Binding site = Dubious). A third toggle shows or hides structures that weren't found /
couldn't be analysed. **Show all** resets every toggle back on.

Each ligand entry has a **View 3D** button that jumps straight to the 3D Viewer tab, loaded with
that ligand, its binding site, and (if available) its electron density.

#### Exporting results

The **Export (.xlsx)** button downloads the current results as a
two-sheet Excel workbook:

- **Parameters** — every threshold/option the analysis was submitted with (RSR, RSCC, R-free,
  occupancy, tolerance, distance, PDB-REDO, advanced checks, blacklist customization, export
  timestamp).
- **Ligands** — one row per ligand: UniProt ID, PDB ID, ligand, ligand
  classification, binding-site classification, R-free, R-work, rejected molecules, and whether an
  electron-density map is available for that structure. Structures that couldn't be analysed are
  still listed, with the fields that have no data left blank.

### 3D Viewer tab

Renders the model with [Mol\*](https://molstar.org). You can:

- Load any PDB ID directly (via the sidebar input), independent of a prior analysis.
- Toggle the **protein**, **ligand**, and **binding site** representations independently.
- Overlay the **2Fo-Fc electron density**, segmented per region (ligand / binding site / residues
  to examine) instead of as one whole-model surface — the backend crops the source map to each
  region and masks it around every atom of that region *before* sending it to the browser, so you only ever see density belonging
  to the residues you're inspecting, and the viewer only has to build one isosurface per region.
- Adjust the **contour level (isovalue, in σ)** and the **per-atom mask radius** live from the
  sidebar sliders.
- Click through a ligand's **components to examine** in the sidebar list to focus the camera on each
  one in turn.

#### Manual quality review

Opening a structure via **View 3D** also populates a **Quality review** panel to the right of the
viewer:

- **Components to examine** — every flagged component gets its own row with Good / Dubious / Bad
  buttons; whichever matches its computed classification starts active. Clicking a different one
  reclassifies that component.
- **Overall classification** — separate Good / Dubious / Bad selectors for the ligand and for the
  binding site as a whole, working the same way.
- **Reset** discards any changes you haven't confirmed yet, reverting to the last confirmed
  classification (or the original computed one, if you haven't confirmed anything for this ligand
  yet).
- **Confirm changes** writes your edits back into the Results tab — the ligand/binding-site
  badges update immediately, and the entry is flagged **✎ edited**.

These overrides live only in the current browser session. They are not sent to the server or
saved to disk, so reloading the page or re-running the analysis discards them.

## How an analysis works

For each PDB ID, VHELIBS downloads the mmCIF model (from RCSB, or from PDB-REDO if that option is
checked) together with per-residue real-space validation statistics (from PDBe's EDS validation
report, or from PDB-REDO's own re-refinement data). Ligands that are known solvents, buffers,
ions, or crystallisation additives — or that are covalently bound to the protein chain — are
filtered out using built-in metal/blacklist tables rather than being scored
as ligands. Remaining ligands are grouped into complexes, and every residue within the binding-site distance cutoff is collected as that
ligand's binding site.

## Classification criteria

Each ligand and binding-site residue is scored against the active criteria; every criterion it
fails adds one "fail point":

| Metric | Default threshold | Effect on score |
|--------|-------------------|------------------|
| RSR (real-space R-factor) | ≤ 0.24 good · 0.24–0.40 → +1 · > 0.40 → +2 | +1 or +2 |
| RSCC (real-space correlation coefficient) | > 0.9 | +1 if below |
| Occupancy | = 1.0 | +1 if below · error if > 1.0 |
| R-free | < 1.0 | +1 if above |

A residue's total score maps to a bucket: `0` → **Good**, `1..tolerance` → **Dubious**,
`> tolerance` → **Bad** (`tolerance` defaults to 2). A ligand or binding site then takes the worst
classification of its own residues (any Bad → Bad; else any Dubious → Dubious; else Good).

**Missing validation data.** A residue with no electron-density validation stats (RSR/RSCC) at all
gets a severe automatic penalty, which almost always
means "Bad" — but this only means *no data was available to judge it*, not necessarily a genuine
density-fit problem. Ligand and binding-site residues are affected differently:

- **Ligand** residues missing data are pruned out of the ligand entirely and reported under
  `rejected` instead — they don't force the *remaining* ligand residues to "Bad" by themselves,
  but if a whole (single-residue) ligand has no data it disappears from the results list
  altogether rather than showing up as "Bad".
- **Binding site** residues missing data are *not* pruned, so a single such residue among many is
  enough to mark the entire binding site "Bad".

Each ligand result carries `ligand_density_data_available`/`binding_site_density_data_available`
booleans so this distinction is visible instead of
being silently folded into an opaque "Bad" — the Results tab shows a ⚠ next to the Ligand/BS
quality badge whenever that data was missing, and its **Not found RSR/RSCC** toggle filters on it.

This is a distinct concept from whether the electron-density **map** itself (the file the 3D
viewer renders, as opposed to the RSR/RSCC stats above) could be found — a structure can have one
without the other.

The Results tab's filter bar has three single on/off toggles alongside the
Ligand/Binding-site quality filters: **Not found complexes** (structures RCSB/PDB-REDO/EDS had no
usable data for at all), **Not found RSR/RSCC** (ligands/binding sites missing validation stats, per above), and
**Not found EDM** (structures missing the density map file). Each is active (visible) by default;
switching one off hides just that subset, independently of the Ligand/Binding-site quality filters
and of each other.

Enabling any of the **advanced checks** below (all off by default) adds further fail points on
the same scale:

| Metric | Default threshold | Notes |
|--------|-------------------|-------|
| OWAB | < 50 | Occupancy-weighted average B-factor, per residue |
| Resolution | ≤ 3.5 Å | Structure-wide, from RCSB or PDB-REDO refinement stats |
| R-diff | ≤ 0.05 | \|R-free − R-work\|, flags possible over-refinement |
| DPI | < 0.42 | Diffraction-component precision index, estimated from cell volume, atom count, reflection count and R-free |

## Project structure

```
vhelibs-web/
├── run.py                   # Entry point
├── requirements.txt
├── app/
│   ├── __init__.py          # Flask app factory
│   ├── routes.py            # REST endpoints + job queue
│   ├── static/
│   │   ├── css/style.css
│   │   └── js/app.js        # SPA frontend + Mol* Viewer integration
│   └── templates/
│       └── index.html
└── core/
    ├── cofactors.py         # Metal/blacklist lookup tables
    ├── pdb_atom.py          # Pure-Python PdbAtom + canonical residue-key formatting
    ├── http_cache.py        # Shared download + on-disk JSON cache helpers
    ├── pdb_utils.py         # PDB/mmCIF download + RCSB stats, UniProt → PDB resolution
    ├── eds_utils.py         # EDS validation XML fetch + electron density map fetch
    ├── pdb_redo_utils.py    # PDB-REDO integration (structures, stats, density map)
    ├── density_mask.py      # Server-side, per-region density masking
    └── rsr_core.py          # Analysis engine — orchestrates the above and returns JSON
```

A few notes on the core modules:

- **`http_cache.py`** is the single place the "check local cache, download with retries, log and
  swallow errors" pattern lives — `pdb_utils.py`, `eds_utils.py`, and `pdb_redo_utils.py` each
  keep their own decisions about *what* to cache and *how* to parse it, but delegate the actual
  network I/O and disk caching here.
- **`pdb_atom.py`**'s `format_reskey()` is the single source of truth for residue-key formatting
  and is used across `rsr_core`, `eds_utils`, and `pdb_redo_utils` — they must agree byte-for-byte
  since these keys are cross-referenced by plain string equality between parsed atoms and
  validation statistics.
- **`pdb_redo_utils.py`**'s `get_EDM()` fetches PDB-REDO's final density map from PDB-REDO's
  [map-maker service](https://pdb-redo.eu/map-maker/map?id=1cbs&stage=final&type=density).
- **`rsr_core.py`** is the analysis engine: it ties structure parsing (via
  [`gemmi`](https://gemmi.readthedocs.io)), validation-statistics lookup, and scoring together,
  and additionally computes the per-region density bounding boxes/atom centers (`density_boxes`,
  `density_atoms`) consumed by the 3D viewer.
- **`density_mask.py`** turns those bounding boxes/atom centers into an actual pre-masked `.ccp4`
  file per region, using `gemmi` to crop and mask the source map.

## REST API

### `POST /api/analyse`

**Request body (JSON):**

```json
{
  "pdbids": "1cbs 3dzu P00734",
  "rsr_upper": 0.4,
  "rsr_lower": 0.24,
  "rscc_min": 0.9,
  "rfree_max": 1.0,
  "occupancy_min": 1.0,
  "tolerance": 2,
  "distance": 4.5,
  "use_pdb_redo": false,
  "check_owab": false,
  "owab_max": 50,
  "check_resolution": false,
  "resolution_max": 3.5,
  "use_rdiff": false,
  "rdiff_max": 0.05,
  "use_dpi": false,
  "dpi_max": 0.42,
  "blacklist": {
    "disabled": ["HOH"],
    "custom": [{ "code": "XYZ", "name": "Custom additive", "category": "blacklist" }],
    "replace": null
  }
}
```

`pdbids` accepts PDB IDs, UniProt IDs, or a mix of both, separated by commas, whitespace,
or newlines. Each UniProt ID is resolved via the RCSB Search API to every PDB entry whose
polymer entities reference it; unresolvable accessions are omitted and reported in an optional
`warnings` array in the response.

`blacklist` is optional and built by the Analysis tab's *Ligand Blacklist* card: `disabled` lists built-in codes to
stop treating as blacklisted for this run, `custom` adds extra `{code, name, category}` entries,
and `replace` — when non-null — substitutes the built-in metal/blacklist tables entirely with the
`{"metals": {...}, "ligand_blacklist": {...}}` shape returned by
[`POST /api/blacklist/parse`](#post-apiblacklistparse).

**Response:**

```json
{ "job_id": "uuid", "total": 2 }
```

`total` reflects the number of PDB entries actually queued for analysis, i.e. after expanding any
UniProt ID in the request.

### `GET /api/blacklist`

Returns the built-in non-ligand-blacklist and metal/ion entries, so the Analysis tab can render
them as toggleable checkboxes rather than leaving them hardcoded/invisible.

```json
{ "entries": [ { "code": "HOH", "name": "Water", "category": "blacklist" }, "..." ] }
```

### `POST /api/blacklist/parse`

Parses an uploaded blacklist file (sent as raw text) into structured entries, so the Analysis tab
can preview it ("this file defines N blacklist + M metal entries") before committing to replacing
the current list with it.

**Request:** `{ "text": "<file contents>" }`

**Response:**

```json
{
  "entries": [ { "code": "XYZ", "name": "...", "category": "blacklist" }, "..." ],
  "metals": { "ZN": "Zinc" },
  "ligand_blacklist": { "HOH": "Water" }
}
```

Returns `400` with `{"error": "..."}` if the text is empty or no valid entries could be parsed.

### `GET /api/status/<job_id>`

**While running:**

```json
{ "status": "running", "progress": 1, "total": 2, "results": null }
```

**When done:**

```json
{
  "status": "done",
  "progress": 2,
  "total": 2,
  "results": [
    {
      "pdbid": "1cbs",
      "uniprot": "P00734",
      "ligands": [
        {
          "ligand_residues": ["REA A  200"],
          "binding_site_residues": ["TYR A   60", "..."],
          "residues_to_examine": ["..."],
          "residue_qualities": { "TYR A   60": "Dubious" },
          "ligand_quality": "Good",
          "binding_site_quality": "Good",
          "source": "PDB",
          "ligand_score": 0,
          "binding_site_score": 0,
          "ligand_density_data_available": true,
          "binding_site_density_data_available": true,
          "low_occupancy": [],
          "other_ligands": [],
          "density_boxes": {
            "ligand": { "min": [0, 0, 0], "max": [0, 0, 0] },
            "binding_site": { "min": [0, 0, 0], "max": [0, 0, 0] },
            "residues_to_examine": { "min": [0, 0, 0], "max": [0, 0, 0] }
          },
          "density_atoms": {
            "ligand": [{ "residue": "REA A  200", "center": [0, 0, 0] }],
            "binding_site": [],
            "residues_to_examine": []
          }
        }
      ],
      "rejected": {},
      "struc_dict": { "rFree": 0.218, "rWork": 0.175 },
      "edm_available": true
    }
  ]
}
```

`density_boxes` gives a padded bounding box per region. `density_atoms` gives the per-atom
coordinates of each region. The 3D viewer sends both, together with the atom-mask radius, to
[`GET /api/density-mask/<pdbid>/<region>`](#get-apidensity-maskpdbidregion), which uses them to
crop and mask the source map server-side. Both fields are only populated for X-ray entries with usable
validation/density data.

`uniprot` is the UniProt ID that produced this PDB entry, or `null` if it was given
directly as a PDB ID. `other_ligands` lists residues belonging to any *other* ligand present in
the same structure (still counted towards this ligand's own scoring, but excluded from
`binding_site_residues`/`residues_to_examine`/`density_boxes`/`density_atoms` so the 3D viewer
never leaks another ligand into the current scene).

`residue_qualities` maps every entry in `residues_to_examine` to the classification it was
computed with (`"Dubious"` or `"Bad"` — by construction, a component never ends up in
`residues_to_examine` if it's `"Good"`). This is what pre-selects the right button for each
component in the 3D Viewer's [manual quality review](#manual-quality-review) panel. Overrides made
there are a browser-only annotation layer. They update the in-memory results the Results tab
renders from, but are never sent back to this API or persisted to disk — re-running the analysis
or reloading the page discards them.

### `GET /api/edm-exists/<pdbid>`

Lightweight existence check for the full 2Fo-Fc electron density map (standard RCSB/EBI source
only), used by the Results tab's [export](#exporting-results) to report map availability without
downloading it. Issues a single `HEAD` request (falling back to a closed-early `GET` if the server
doesn't support `HEAD`); the answer is cached on disk.

```json
{ "pdbid": "1cbs", "exists": true }
```

### `GET /api/density-mask/<pdbid>/<region>`

Returns a single pre-masked `.ccp4` density map for one region (`region` is `ligand`,
`binding_site`, or `residues_to_examine`), cropped and masked server-side.

**Query parameters:**

| Param | Required | Description |
|-------|----------|--------------|
| `min`, `max` | Yes | Region bounding box corners, `"x,y,z"` — from `density_boxes[region]`. |
| `atoms` | Yes | Semicolon-separated atom centers, `"x1,y1,z1;x2,y2,z2;..."` — from `density_atoms[region]`. |
| `radius` | No (default `1.6`) | Atom-mask radius in Å; quantized to 0.25 Å steps and cached. |
| `source` | No (default `pdb`) | `pdb` (standard RCSB/EBI map) or `pdb_redo`. |

Returns the raw `.ccp4` file (`application/octet-stream`) on success, or `404` with
`{"error": "..."}` if the region is invalid or no map could be produced (e.g. no source map
available for this entry/source).

### `POST /api/cache/clear`

Deletes everything under the server's on-disk cache. The cache is otherwise never cleaned up
automatically, so this is the only way to reclaim that disk space. Powers the Analysis tab's
[**Clear disk cache**](#clearing-the-disk-cache) button; takes no body/parameters.

**Response:**

```json
{ "removed_files": 128, "freed_bytes": 47185920, "errors": [] }
```

`errors` lists any cache entries that could not be deleted (e.g. a permissions problem) — the rest
of the cleanup still proceeds.

## Scripting the API / example script

Everything the web UI does goes through the plain HTTP endpoints documented above — you don't
need a browser to run analyses. The basic flow to script is: `POST /api/analyse` to start a job,
then poll `GET /api/status/<job_id>` until `status` is `"done"` and read `results` from there. This makes it straightforward to run VHELIBS Web
headlessly on a server — e.g. as part of a batch pipeline that screens many structures overnight
and exports the results without anyone opening a browser.

`examples/vhelibs_batch_client.py` is a self-contained example of exactly that. It:

1. Starts the Flask server itself (`run.py`) as a subprocess and waits until it responds.
2. Submits the given PDB/UniProt IDs to `/api/analyse` in small batches rather than all at once,
   polling `/api/status/<job_id>` for each batch until it's done.
3. Flattens each batch's results into spreadsheet rows as soon as they arrive — discarding the raw
   analysis JSON (atom coordinates, density boxes, per-residue quality maps, …) for that batch
   right away — and, once every batch is processed, writes a two-sheet `.xlsx` (`Ligands` +
   `Parameters`) in the same layout the Results tab's own
   [Export (.xlsx)](#exporting-results) button produces.

```bash
pip install requests openpyxl
python examples/vhelibs_batch_client.py \
    --repo-root /path/to/vhelibs-web \
    --pdbids 1cbs 3dzu 4hhb P00734 \
    --batch-size 5 \
    --output results.xlsx
```

`--repo-root` is the directory containing `run.py`; the script launches and stops the server for
you, so it doesn't need to be running beforehand. Run
`python examples/vhelibs_batch_client.py --help` for the full list of options, including the
quality thresholds, `--use-pdb-redo`, and the polling interval.

## Data sources

| Purpose | Source |
|---------|--------|
| Structure files | [RCSB PDB](https://www.rcsb.org) (mmCIF format) |
| Validation statistics | [PDBe validation reports](https://www.ebi.ac.uk/pdbe) |
| Re-refined structures & statistics | [PDB-REDO](https://pdb-redo.eu) (optional alternative to RCSB) |
| Electron density maps | EBI's density server (segmented 2Fo-Fc/Fo-Fc volumes) |
| UniProt → PDB resolution | RCSB Search API |
| 3D visualisation | [Mol\* Viewer](https://molstar.org) |

## Citation

Cereto-Massagué A *et al.* **VHELIBS: a validation helper for ligands and binding sites.**
*J Cheminform* 5, 36 (2013). <https://doi.org/10.1186/1758-2946-5-36>
