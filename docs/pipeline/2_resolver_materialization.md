> Historical pre-Parquet pipeline. These commands are retained for context; use the [current architecture](../architecture.md) and [root workflow](../../README.md).

# 2. Resolver Materialization

Resolver materialization builds parquet lookup tables used later by DuckDB to
resolve source identifiers into canonical entity identifiers.

## Command

```bash
make resolver
```

Useful variants:

```bash
make resolver MAX_RECORDS=100000
make resolver PUBCHEM_URL=https://example.org/pubchem.sdf.gz
make resolver RESOLVER_SOURCES=uniprot
```

## What It Produces

By default, resolver files are written under `data/`:

- `data/proteins/protein_identifier_lookup.parquet`
- `data/proteins/identifier_type.parquet`
- `data/chemicals/lookup/<source>.parquet`
- `data/chemicals/ambiguous/<source>.parquet`
- `data/chemicals/identifier_type.parquet`

Protein mappings resolve to primary UniProt accessions and preserve taxonomy for
species-scoped resolution. Chemical mappings resolve identifiers such as ChEBI,
ChEMBL, HMDB, LipidMaps, SwissLipids, and PubChem toward standard InChI keys.

## Main Files

- `Makefile`: `resolver` target.
- `omnipath_build/cli.py`: `build-resolver` command.
- `omnipath_build/resolver/mapping_tables.py`: resolver source orchestration.
- `omnipath_build/resolver/sources/proteins.py`: protein lookup builder.
- `omnipath_build/resolver/sources/chemicals.py`: chemical lookup builder.
- `omnipath_build/resolver/sources/pubchem.py`: PubChem SDF materialization.

## Phase Boundary

This phase does not write PostgreSQL content tables. It only prepares parquet
inputs consumed by the DuckDB load phase.
