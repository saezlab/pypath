> Historical pre-Parquet pipeline. These commands are retained for context; use the [current architecture](../architecture.md) and [root workflow](../../README.md).

# 1. Setup

Setup prepares the local development environment and editable submodules.

## Command

```bash
make setup
```

## What It Does

The `setup` target:

- Adds or updates the `pypath` submodule.
- Adds or updates the `download-manager` submodule.
- Updates submodules recursively.
- Runs `uv sync`.
- Installs local editable packages for `cache-manager`, `download-manager`, and
  `pypath`.

## Main Files

- `Makefile`: `setup` target.
- `pyproject.toml`: Python package and dependency metadata.
- `uv.lock`: locked dependency graph.
- `pypath/`: local pypath checkout used for source datasets.
- `download-manager/`: local download manager checkout.
- `cache-manager/`: local cache manager checkout.

## Notes

The build expects pypath downloads to use a project-local data directory unless
`PYPATH_DOWNLOAD_DATADIR` is already set. Source discovery configures this in
`omnipath_build/resources.py`.

