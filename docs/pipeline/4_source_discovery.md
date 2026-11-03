> Historical pre-Parquet pipeline. These commands are retained for context; use the [current architecture](../architecture.md) and [root workflow](../../README.md).

# 4. Source Discovery

Source discovery finds runnable `pypath.inputs_v2` datasets without maintaining a
hand-written source list in `omnipath_build`.

## Where It Runs

Discovery happens inside `make load`, through:

```bash
python -m omnipath_build.duckdb_direct_pipeline
```

## What It Does

Discovery imports the configured inputs package, walks its modules, and collects
pypath dataset objects:

- `Resource`
- `Dataset`
- `ArtifactDataset`

The load pipeline selects only entity datasets with raw dataset access.
Id-translation datasets and artifact-only datasets are skipped for evidence
ingest.

## Selection Rules

By default, all discovered sources are eligible except sources in the loader's
default exclusion set, currently `imm1415` and `rampdb`.

Selection can be narrowed with:

```bash
make load SOURCE=bindingdb
make load SOURCES=uniprot,bindingdb,intact
make load SOURCE=bindingdb DATASET=interactions
```

## Main Files

- `omnipath_build/resources.py`: general discovery implementation.
- `omnipath_build/duckdb_direct_pipeline.py`: load-specific filtering.
- `pypath/pypath/inputs_v2/`: source dataset definitions.

## Download Directory

If `PYPATH_DOWNLOAD_DATADIR` is not set, discovery points pypath downloads at
the project-local `pypath-data/` directory.
