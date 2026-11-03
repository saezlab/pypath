# AGENTS.md — omnipath-utils Development Guide

## Project Overview

omnipath-utils is a bioinformatics utility service providing ID translation,
taxonomy resolution, orthology mapping, and reference lists for molecular
biology. It serves both as a local Python library and as an HTTP web service
backed by PostgreSQL.

**Key facts:**
- Python 3.10+, Hatchling build, GPL-3.0
- 7400+ lines of code across 45 modules
- 211 tests (3000+ lines)
- Deployed at https://utils.omnipathdb.org
- PyPI: `pip install omnipath-utils`

## Architecture

```
omnipath_utils/
├── mapping/           # ID translation (core feature)
│   ├── __init__.py    # Public API: map_name, translate, translate_column
│   ├── _translate.py  # translate_core: unified entry point for all APIs
│   ├── _mapper.py     # Mapper: table management, backend dispatch
│   ├── _table.py      # MappingTable: dict wrapper with expiry
│   ├── _reader.py     # MapReader: loads data from backends, pickle cache
│   ├── _id_types.py   # IdTypeRegistry: 110+ ID types from id_types.yaml
│   ├── _cleanup.py    # UniProt cleanup: sec→pri, TrEMBL→SwissProt
│   ├── _special.py    # Fallbacks: gene symbol case, RefSeq, miRNA, chain
│   └── backends/      # 8 data source backends
│       ├── _base.py       # MappingBackend ABC (handles common logic)
│       ├── _uniprot.py    # UniProt REST + pypath fallback
│       ├── _uniprot_ftp.py # UniProt FTP bulk idmapping
│       ├── _biomart.py    # Ensembl BioMart
│       ├── _uploadlists.py # UniProt ID Mapping batch
│       ├── _mirbase.py    # miRBase
│       ├── _unichem.py    # UniChem (small molecules)
│       ├── _ramp.py       # RaMP (metabolites + synonyms)
│       └── _hmdb.py       # HMDB
├── taxonomy/          # Organism name resolution
│   ├── __init__.py    # Public API: ensure_ncbi_tax_id, etc.
│   └── _taxonomy.py   # TaxonomyManager: 28K organisms, synonyms, API loading
├── orthology/         # Cross-species gene translation
│   ├── __init__.py    # Public API: translate, translate_column
│   └── _manager.py    # OrthologyManager: HCOP, OMA, Ensembl, HomoloGene, OrthoDB, Alliance
├── reflists/          # Reference lists (SwissProt, TrEMBL sets)
│   ├── __init__.py    # Public API: all_swissprots, is_swissprot
│   └── _manager.py    # ReferenceListManager
├── db/                # PostgreSQL database
│   ├── _schema.py     # SQLAlchemy 2.0 ORM models
│   ├── _connection.py # Engine + raw psycopg3 connections
│   ├── _build.py      # DatabaseBuilder: presets, populate_*, export_parquet
│   ├── _query.py      # translate_ids, identify_ids, chain_translate
│   ├── _presets.py    # Build presets: minimal, standard, model, full
│   └── _loader.py     # On-demand background table loading
├── server/            # Litestar web service
│   ├── _app.py        # App factory, landing page, health, favicon
│   ├── _routes_mapping.py    # /mapping/* endpoints
│   ├── _routes_taxonomy.py   # /taxonomy/* endpoints
│   ├── _routes_orthology.py  # /orthology/* endpoints
│   └── _routes_reflists.py   # /reflists/* endpoints
├── cli/               # Command-line interface
│   ├── _main.py       # Entry point
│   ├── _build.py      # `omnipath-utils build` with presets
│   └── _serve.py      # `omnipath-utils serve`
├── data/
│   ├── id_types.yaml  # 110+ ID type definitions with backends, CURIE, aliases
│   └── organisms.yaml # 28 seed organisms with synonyms
├── _constants.py      # DEFAULT_ORGANISM, NOT_ORGANISM_SPECIFIC
├── _session.py        # pkg_infra session
└── _metadata.py       # Version
```

## Key design patterns

### Dual-mode: memory + database
All translation functions work in two modes:
- **Memory mode**: downloads data from upstream APIs, caches as pickle files
- **Database mode**: queries pre-populated PostgreSQL

The `translate_core()` function is the single entry point for all APIs
(Python, REST, DataFrame). It uses vectorized batch lookup from tables,
with per-ID fallback through `map_name()` for misses.

### Backend system
Backends inherit from `MappingBackend` (in `_base.py`). The base class
handles ID type column resolution, pypath-vs-direct dispatch, and logging.
Subclasses declare `name` and `yaml_key` and implement `_read_via_pypath()`
and `_read_direct()`.

### pypath integration
Backends try `pypath.inputs` first, fall back to direct HTTP. This means
the package works both with and without pypath installed. Data access
through pypath uses `dlmachine` (not the old `pypath.share.curl`).

### REST fallback chain
The REST API does: DB query → gene symbol fallbacks (via DB) → chain
translation (single SQL JOIN query) → on-demand background loading for
missing tables.

### Build presets
```bash
omnipath-utils build --preset minimal   # human only
omnipath-utils build --preset standard  # human + mouse + metabolites
omnipath-utils build --preset model     # 8 model organisms + everything
omnipath-utils build --preset full      # all organisms via UniProt FTP
```

## Dependencies

### Our packages
- `pkg-infra` — logging, config, session, utility functions
- `cachedir` — file caching with SQLite metadata
- `dlmachine` — HTTP downloads
- `pypath-omnipath` (optional) — data access for backends

### Key external
- `narwhals` — backend-agnostic DataFrame ops (pandas/polars/pyarrow)
- `pydantic` — validation
- `sqlalchemy` + `psycopg` — PostgreSQL (optional, for db mode)
- `litestar` + `uvicorn` — web service (optional)

## Testing
```bash
uv sync --extra tests
uv run pytest tests/ -v
```

Tests use mocks — no live HTTP calls or database required. 211 tests in ~50s.

## Deployment
```bash
# Production at utils.omnipathdb.org
~/deploy/omnipath-utils/restart.sh
~/deploy/omnipath-utils/status.sh

# Build database
omnipath-utils build --preset model \
    --db-url postgresql+psycopg://... \
    --parquet-dir ./parquet

# Serve
omnipath-utils serve --port 8083
```

Dev PostgreSQL: port 5433 (dev), port 5434 (deploy).
NixOS: needs `source ~/dev/.envrc` for libstdc++/libz paths.

## Common tasks

### Adding a new mapping backend
1. Create `mapping/backends/_newbackend.py`
2. Subclass `MappingBackend`, set `name` and `yaml_key`
3. Implement `_read_via_pypath()` and `_read_direct()`
4. Call `register('name', NewBackend)` at module level
5. Add ID types to `data/id_types.yaml` if needed

### Adding a new ID type
1. Add to `data/id_types.yaml` with label, entity_type, curie_prefix, backends, aliases
2. The IdTypeRegistry auto-loads it

### Adding a new REST endpoint
1. Create or update route in `server/_routes_*.py`
2. Register controller in `server/_app.py`
