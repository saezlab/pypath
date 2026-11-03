# omnipath-build

A general database builder on top of pypath.

## Setup

```bash
make setup
```

This project uses `uv` for dependency management.

## DuckDB Direct Pipeline

The `omnipath_build` pipeline streams `inputs_v2` sources into DuckDB evidence
tables, canonicalizes them in DuckDB, and copies the projected evidence and
canonical rows into PostgreSQL.

For the data model, phase boundaries, refresh semantics, and common workflows,
see [docs/pipeline.md](docs/pipeline.md).

For tuning the build's Postgres memory and parallelism (and how to size it for
your deployment, including the lab's docker memory cap), see
[docs/build-tuning.md](docs/build-tuning.md).

The default database URL is:

```bash
postgresql://omnipath:omnipath@localhost:55432/omnipath
```

Override it with `DATABASE_URL=...` when needed.

### Build Resolver Files

Build local resolver parquet files:

```bash
make resolver
```

Limit resolver builds for smoke tests:

```bash
make resolver MAX_RECORDS=100000
```

Use a single PubChem SDF shard during development:

```bash
make resolver PUBCHEM_URL=https://example.org/pubchem.sdf.gz
```

### Prepare Database

Create schema and supporting indexes:

```bash
make db-setup
```

Start from a clean schema. This defers secondary evidence indexes until
canonicalization, so ingest is faster:

```bash
make db-setup DROP_EXISTING=1
```

Use another schema:

```bash
make db-setup SCHEMA=omnipath_test DROP_EXISTING=1
```

Drop and recreate the target schema without loading resolver tables:

```bash
make db-reset
```

### Load And Canonicalize

Load any sources that are not already present. The DuckDB/PostgreSQL loader
projects evidence, canonicalizes it, and copies the result into PostgreSQL:

```bash
make load
```

Load one source if it is not already present:

```bash
make load SOURCE=bindingdb
```

Load multiple missing sources:

```bash
make load SOURCES=uniprot,bindingdb,intact
```

Use a shared staging worker pool:

```bash
make load LOAD_JOBS=5
```

`LOAD_JOBS` is used by the staged loader as one shared pool. Workers are
assigned across sources first; when preparse shards are available, idle workers
can stage shards from the same source while PostgreSQL COPY continues to run
serially.

Large staged loads keep simple preparse parquet shards under
`pypath-data/<source>/preparse/`. These shards are reused by default and are
rebuilt when `FORCE_REFRESH=1` is set.

Refresh existing source content by deleting it first, then loading current
parser output:

```bash
make reload SOURCE=bindingdb
```

Run `make derive` after loading the selected sources to refresh query indexes,
derived count/search tables, and bitmaps.

### Full Scratch Build

Build resolver files, recreate the database schema, load all sources through
the DuckDB/PostgreSQL pipeline, and derive query tables:

```bash
make all DROP_EXISTING=1
```

If resolver files already exist, run:

```bash
make db-setup DROP_EXISTING=1
make load
make derive
```

### Test Runs

Limit source rows per dataset:

```bash
make load SOURCE=bindingdb MAX_RECORDS=200000 SCHEMA=omnipath_test
```

The default load batch size is `BATCH_SIZE=50000`, which is safer
for full loads.

### Maintenance

Reset omnipath_build content tables without dropping resolver tables:

```bash
make reset-content
```

Run the Python tests:

```bash
uv run pytest tests -q
```

Check table sizes when PostgreSQL runs in Docker:

```bash
docker exec -i omnipathv2-main-mmvxvb-omnipathv2-postgres-1 \
  psql -U omnipath -d omnipath -v ON_ERROR_STOP=1 <<'SQL'
SELECT
  table_name,
  pg_size_pretty(pg_total_relation_size(format('%I.%I', table_schema, table_name)::regclass)) AS total_size,
  pg_size_pretty(pg_relation_size(format('%I.%I', table_schema, table_name)::regclass)) AS table_size,
  pg_size_pretty(pg_indexes_size(format('%I.%I', table_schema, table_name)::regclass)) AS indexes_size
FROM information_schema.tables
WHERE table_schema = 'public'
  AND table_type = 'BASE TABLE'
ORDER BY pg_total_relation_size(format('%I.%I', table_schema, table_name)::regclass) DESC;
SQL
```
