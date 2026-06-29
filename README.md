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
- Packages: `flask`, `requests`, `pdbx`

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
│   │   └── js/app.js        # SPA frontend + NGL Viewer integration
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
| `rsr_analysis.py` | `core/rsr_core.py` + `app/routes.py` | Removed argparse/multiprocessing; returns dicts |
| `PDBfiles.py` | `core/pdb_utils.py` | `urllib` → `requests`; proper SSL |
| `EDS_parser.py` | `core/eds_utils.py` | Same; proper SSL; logging |
| `pdb_redo.py` | `core/pdb_redo_utils.py` | Removed Java locale hack |
| `PdbAtom.py` | `core/pdb_atom.py` | Kept pure Python; added `__hash__` |
| `cofactors.py` | `core/cofactors.py` | Essentially unchanged |
| `cPdbAtom.pyx` | *(removed)* | Pure Python is fast enough |
| `visualitzador.py` | `app/templates/index.html` + `app/static/js/app.js` | Jmol/Swing → NGL Viewer |
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
  "pdbids": "1cbs 3dzu",
  "rsr_upper": 0.4,
  "rsr_lower": 0.24,
  "rscc_min": 0.9,
  "rfree_max": 1.0,
  "occupancy_min": 1.0,
  "tolerance": 2,
  "distance": 4.5,
  "use_pdb_redo": false
}
```

**Response:**

```json
{ "job_id": "uuid", "total": 2 }
```

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
          "low_occupancy": []
        }
      ],
      "rejected": {},
      "struc_dict": { "rFree": 0.218, "rWork": 0.175 }
    }
  ]
}
```

## Citation

Cereto-Massagué A et al. *VHELIBS: a validation helper for ligands and binding sites.*
J Cheminform 5, 36 (2013). <https://doi.org/10.1186/1758-2946-5-36>

## License

Copyright 2012–2024 Adrià Cereto Massagué.
