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
  - [Results tab](#results-tab)
  - [3D Viewer tab](#3d-viewer-tab)
- [How an analysis works](#how-an-analysis-works)
- [Classification criteria](#classification-criteria)
- [Project structure](#project-structure)
- [REST API](#rest-api)
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

- **Batch analysis** of PDB IDs and/or UniProt accessions (each UniProt accession is
  auto-expanded, via the RCSB Search API, to every PDB structure that references it).
- **Configurable quality thresholds** with sensible crystallographic defaults (RSR, RSCC,
  occupancy, R-free, tolerance, binding-site distance cutoff).
- **Optional advanced checks**, off by default: OWAB, structure-wide resolution, R-diff
  (|R-free − R-work|), and DPI (diffraction-component precision index).
- **PDB-REDO integration** — analyse the re-refined structure and statistics instead of the
  original RCSB/PDBe deposition.
- **Interactive results view**: every ligand and binding site is scored independently, and the
  Results tab lets you toggle Good/Dubious/Bad on each axis (plus a separate toggle for structures
  that couldn't be analysed) to filter down to exactly the combination you care about — e.g. every
  *Bad* ligand sitting in a *Dubious* binding site.
- **Integrated Mol\* 3D viewer** with independently toggleable protein/ligand/binding-site
  layers, per-region 2Fo-Fc electron density overlay (streamed on demand and clipped to a sphere
  around each atom rather than the whole structure), and live-adjustable contour level (isovalue)
  and atom-mask radius.
- **Disk caching** of every downloaded structure/statistics file, so re-running an analysis (or
  re-opening the 3D viewer) doesn't re-hit external APIs.

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
python run.py [--host HOST] [--port PORT] [--cache-dir PATH] [--no-browser] [--debug]
```

| Flag | Default | Description |
|------|---------|-------------|
| `--host` | `127.0.0.1` | Bind address |
| `--port` | `8000` | TCP port |
| `--cache-dir` | `~/.cache/vhelibs` | Directory for downloaded PDB/EDS files |
| `--no-browser` | off | Suppress automatic browser launch |
| `--debug` | off | Enable Flask debug/reloader (dev only) |

### Requirements

- Python 3.8+
- Packages: `flask`, `requests`, `gemmi` (see `requirements.txt`)
- A browser with WebGL support (for the 3D viewer, powered by [Mol\*](https://molstar.org))

## Using the app

The app is a single-page interface with four tabs: **Analysis**, **Results**, **3D Viewer**, and
**About**.

### Analysis tab

1. Enter one or more **PDB IDs and/or UniProt accessions** in the text box, separated by commas,
   whitespace, or newlines (e.g. `1cbs, 3dzu, 4hhb, P00734`) — or load them from a `.txt`/`.csv`
   file. UniProt accessions are automatically resolved to every associated PDB structure.
2. Optionally check **Use PDB-REDO structures** to analyse PDB-REDO's re-refined coordinates and
   statistics instead of the standard RCSB/PDBe deposition.
3. Adjust the **quality thresholds** (RSR, RSCC, R-free, occupancy, tolerance, binding-site
   distance) if the defaults don't suit your use case, and optionally expand **Advanced options**
   to enable OWAB / resolution / R-diff / DPI checks.
4. Click **Analyse**. A progress bar tracks the job as each structure is downloaded and scored;
   when it finishes you're taken straight to the Results tab.

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

### 3D Viewer tab

Renders the model with [Mol\*](https://molstar.org). You can:

- Load any PDB ID directly (via the sidebar input), independent of a prior analysis.
- Toggle the **protein**, **ligand**, and **binding site** representations independently.
- Overlay the **2Fo-Fc electron density**, segmented per region (ligand / binding site / residues
  to examine) instead of as one whole-model surface — density is streamed on demand from EBI's
  density server as small boxes around each region, then clipped to a small sphere around every
  atom of that region, so you only ever see density belonging to the residues you're inspecting.
- Adjust the **contour level (isovalue, in σ)** and the **per-atom mask radius** live from the
  sidebar sliders.
- Click through a ligand's **components to examine** in the sidebar list to focus the camera on each
  one in turn.

## How an analysis works

For each PDB ID, VHELIBS downloads the mmCIF model (from RCSB, or from PDB-REDO if that option is
checked) together with per-residue real-space validation statistics (from PDBe's EDS validation
report, or from PDB-REDO's own re-refinement data). Ligands that are known solvents, buffers,
ions, or crystallisation additives — or that are covalently bound to the protein chain — are
filtered out using built-in metal/blacklist tables rather than being scored as ligands. Remaining
ligands are grouped into complexes (covalently linked HETATM groups count as one ligand), and
every residue within the binding-site distance cutoff is collected as that ligand's binding site.

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
Residues with missing validation data are excluded from scoring and reported under `rejected`
instead.

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
- **`rsr_core.py`** is the analysis engine: it ties structure parsing (via
  [`gemmi`](https://gemmi.readthedocs.io)), validation-statistics lookup, and scoring together,
  and additionally computes the per-region density bounding boxes/atom centers (`density_boxes`,
  `density_atoms`) consumed by the 3D viewer.

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
  "dpi_max": 0.42
}
```

`pdbids` accepts PDB IDs, UniProt accessions, or a mix of both, separated by commas, whitespace,
or newlines. Each UniProt accession is resolved via the RCSB Search API to every PDB entry whose
polymer entities reference it; unresolvable accessions are omitted and reported in an optional
`warnings` array in the response.

**Response:**

```json
{ "job_id": "uuid", "total": 2 }
```

`total` reflects the number of PDB entries actually queued for analysis, i.e. after expanding any
UniProt accessions in the request.

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
          "ligand_quality": "Good",
          "binding_site_quality": "Good",
          "source": "PDB",
          "ligand_score": 0,
          "binding_site_score": 0,
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
      "struc_dict": { "rFree": 0.218, "rWork": 0.175 }
    }
  ]
}
```

`density_boxes` gives a padded bounding box per region, used to size the on-demand download
window from EBI's density server. `density_atoms` gives the per-atom coordinates of each region,
used by the 3D viewer to clip displayed density to a small sphere around each atom. Both are only
populated for X-ray entries with usable validation/density data.

`uniprot` is the UniProt accession that produced this PDB entry, or `null` if it was given
directly as a PDB ID. `other_ligands` lists residues belonging to any *other* ligand present in
the same structure (still counted towards this ligand's own scoring, but excluded from
`binding_site_residues`/`residues_to_examine`/`density_boxes`/`density_atoms` so the 3D viewer
never leaks another ligand into the current scene).

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
