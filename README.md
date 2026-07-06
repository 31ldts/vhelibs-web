# VHELIBS Web

**Validation Helper for LIgands and Binding Sites** – web-based version.

## Quick Start

```bash
pip install -r requirements.txt
python run.py
```

Then open <http://localhost:8000> in your browser.

## Requirements

- Python 3.8+
- Packages: `flask`, `requests`, `gemmi`

## Usage

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

## Project Structure

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
    ├── cofactors.py         # Metal/blacklist data (unchanged from original)
    ├── pdb_atom.py          # Pure-Python PdbAtom (replaces Cython cPdbAtom)
    ├── pdb_utils.py         # PDB/mmCIF download + RCSB stats (from PDBfiles.py)
    ├── eds_utils.py         # EDS validation XML fetch (from EDS_parser.py)
    ├── pdb_redo_utils.py    # PDB-REDO integration (from pdb_redo.py)
    └── rsr_core.py          # Analysis engine returning JSON (from rsr_analysis.py)
```

## Migration Notes

| Old file | New location | Key changes |
|----------|-------------|-------------|
| `rsr_analysis.py` | `core/rsr_core.py` + `app/routes.py` | Removed argparse/multiprocessing; returns dicts; mmCIF parsing now via `gemmi` instead of `pdbx`; adds PDB-REDO, OWAB/resolution/R-diff/DPI checks, and per-region density boxes/atoms for the 3D viewer |
| `PDBfiles.py` | `core/pdb_utils.py` | `urllib` → `requests`; proper SSL |
| `EDS_parser.py` | `core/eds_utils.py` | Same; proper SSL; logging |
| `pdb_redo.py` | `core/pdb_redo_utils.py` | Removed Java locale hack |
| `PdbAtom.py` | `core/pdb_atom.py` | Kept pure Python; added `__hash__` |
| `cofactors.py` | `core/cofactors.py` | Essentially unchanged |
| `cPdbAtom.pyx` | *(removed)* | Pure Python is fast enough |
| `visualitzador.py` | `app/templates/index.html` + `app/static/js/app.js` | Jmol/Swing → NGL Viewer → Mol* Viewer |
| `Main.java`, `PdbAtomJava.java` | *(removed)* | No Java required |
| `multithreading.py` | *(removed)* | Flask threaded + Python threading |
| `argparse.py` | *(removed)* | Bundled copy not needed in Python 3 |
| `setup.py` | *(removed)* | No Cython compilation |
| `vhelibs.sh` | *(removed)* | Use `python run.py` |

## API

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
or newlines. Each UniProt accession (e.g. `P00734`) is resolved via the RCSB Search API to every
PDB entry whose polymer entities reference it, and every one of those entries is queued for
analysis alongside any plain PDB IDs given. If a UniProt accession doesn't resolve to any entry,
it's omitted and reported back in an optional `warnings` array in the response (see below).

The last eight fields are the *advanced* (opt-in) checks, off by default. When enabled, each
adds its own pass/fail criterion on top of the core RSR/RSCC/occupancy/R-free scoring:

| Field | Checks | Notes |
|-------|--------|-------|
| `check_owab` / `owab_max` | Occupancy-weighted average B-factor, per residue | fails if OWAB ≥ `owab_max` |
| `check_resolution` / `resolution_max` | Structure-wide resolution | from RCSB or PDB-REDO refinement stats |
| `use_rdiff` / `rdiff_max` | \|R-free − R-work\| | flags possible over-refinement |
| `use_dpi` / `dpi_max` | Diffraction-component precision index | estimated from cell volume, atom count, reflection count and R-free |

When `use_pdb_redo` is true, both the structure model and the per-residue validation statistics
are fetched from [PDB-REDO](https://pdb-redo.eu) instead of RCSB/PDBe.

**Response:**

```json
{ "job_id": "uuid", "total": 2 }
```

If `pdbids` included a UniProt accession that couldn't be resolved to any PDB entry, the
response also includes:

```json
{ "job_id": "uuid", "total": 2, "warnings": ["No PDB entries found for UniProt accession P99999"] }
```

`total` reflects the number of PDB entries actually queued for analysis, i.e. after expanding
any UniProt accessions in the request.

### `GET /api/status/<job_id>`

**Response while running:**

```json
{ "status": "running", "progress": 1, "total": 2, "results": null }
```

**Response when done:**

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
          "binding_site_residues": ["TYR A   60", ...],
          "residues_to_examine": [...],
          "ligand_quality": "Good",
          "binding_site_quality": "Good",
          "source": "PDB",
          "ligand_score": 0,
          "binding_site_score": 0,
          "low_occupancy": [],
          "other_ligands": [],
          "density_boxes": {
            "ligand": { "min": [x, y, z], "max": [x, y, z] },
            "binding_site": { "min": [x, y, z], "max": [x, y, z] },
            "residues_to_examine": { "min": [x, y, z], "max": [x, y, z] }
          },
          "density_atoms": {
            "ligand": [{ "residue": "REA A  200", "center": [x, y, z] }, ...],
            "binding_site": [...],
            "residues_to_examine": [...]
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
window from EBI's density server (see `core/eds_utils.py:edm_box_url`). `density_atoms` gives
the actual per-atom coordinates of each region, used by the 3D viewer to clip the displayed
density to a small sphere around each atom rather than showing everything inside the box.
Both are only populated for X-ray entries with usable validation/density data; a region with no
atoms (e.g. a fully solvent-exposed ligand with `distance: 0`) yields `null`/an empty list.

`uniprot` is the UniProt accession that produced this PDB entry, or `null` if it was given
directly as a PDB ID. `other_ligands` lists residues belonging to any *other* ligand present in
the same structure (e.g. a second ligand close enough to have been pulled into this one's
binding site) — these still count towards this ligand's own scoring, but are deliberately
excluded from `binding_site_residues`, `residues_to_examine`, `density_boxes` and
`density_atoms`, so the 3D viewer for a given ligand never shows another ligand from the same
structure.

## Classification

Each ligand and binding-site residue is scored against the active criteria; every criterion it
fails adds one "fail point":

| Metric | Default threshold | Effect on score |
|--------|-------------------|------------------|
| RSR | ≤ 0.24 good · 0.24–0.40 → +1 · > 0.40 → +2 | +1 or +2 |
| RSCC | > 0.9 | +1 if below |
| Occupancy | = 1.0 | +1 if below · error if > 1.0 |
| R-free | < 1.0 | +1 if above |

A residue's total score maps to a bucket: `0` → **Good**, `1..tolerance` → **Dubious**,
`> tolerance` → **Bad** (`tolerance` defaults to 2). A ligand or binding site then takes the
worst classification of its own residues (any Bad → Bad; else any Dubious → Dubious; else Good).
Residues with missing validation data are excluded from scoring and reported under `rejected`
instead. Enabling any of the advanced checks (OWAB, resolution, R-diff, DPI — see the API section
above) adds further +1 fail points on the same scale.

## 3D Viewer

The Viewer tab renders the model with [Mol*](https://molstar.org), with independent toggleable
layers for protein, ligand, and binding site. When a 2Fo-Fc electron density map is available,
it can be overlaid per region (ligand / binding site / residues to examine) instead of as one
whole-model surface: density is streamed on demand from EBI's density server as small boxes
around each region (`density_boxes`) and clipped to a sphere around every atom of that region
(`density_atoms`), so only density belonging to the residues being inspected is shown. Both the
contour level (isovalue, in σ) and the per-atom mask radius are adjustable live from the sidebar.-->

## Citation

Cereto-Massagué A et al. *VHELIBS: a validation helper for ligands and binding sites.*
J Cheminform 5, 36 (2013). <https://doi.org/10.1186/1758-2946-5-36>

<!--## License

Copyright 2012–2024 Adrià Cereto Massagué.-->
