# nicesrv serving deployment — 8 October 2026

The API, explorer and HTTPS file service follow branch `parquet-migration`:
`deploy.sh` checks out the branch head, builds images tagged with its short
commit and restarts the stack. The deployed commit is recorded in
`source-commit.txt`.

- Explorer: https://omnipath-metabo-dev.schaul.click/explore
- Files: https://data.omnipath-metabo-dev.schaul.click
- Deployment root: `/root/projects/omnipath-releases/20261007-serving-tables`
  (`deploy.sh`, `serving.env`, `source-commit.txt` and one log per deploy).
- Source checkout used by `deploy.sh`: `/root/projects/omnipath-migration/wt-normalized`
  (detached at `origin/parquet-migration`).
- Artifact root: `/root/data/release-20261007-native` (release `2026.10.9`, all 46
  resources), mounted read-only at `/data`.
- Project: `omnipath-serving-20261007-full`; loopback API/web/files ports:
  `8295` / `8292` / `8290`.
- Maintained checkout for builds and loads: `/root/projects/omnipath-migration/omnipath-build`.

## Deploy, inspect or restart

```sh
ssh -T -o BatchMode=yes -o ForwardAgent=no -o ConnectTimeout=15 nicesrv
deployment_root=/root/projects/omnipath-releases/20261007-serving-tables
# Deploy the current head of parquet-migration (push first).
"$deployment_root/deploy.sh"
# Inspect the running stack.
cd /root/projects/omnipath-migration/wt-normalized
make serving-status COMPOSE_ENV="$deployment_root/serving.env" \
  COMPOSE_FILES='-f compose.serving.yaml -f deploy/compose.traefik.yaml'
make serving-logs COMPOSE_ENV="$deployment_root/serving.env" \
  COMPOSE_FILES='-f compose.serving.yaml -f deploy/compose.traefik.yaml'
```

A deploy builds and starts serving only; it runs no resource build, identity
build or PostgreSQL load. Changing the served data means changing `DATA_DIR` in
`serving.env`.

## Effective API limits

Set in `serving.env`:

- Container: 13 GB RAM and 7 CPUs, leaving 1 CPU and about 2.6 GB of the host's
  8 CPUs and 16 GB for the other apps.
- Query pool: three concurrent DuckDB queries, each with two threads and 3,500 MB.
- Capacity wait: up to 30 seconds, then HTTP 503 with a retry hint.
- Cache warming at startup (`API_WARM_CACHE=true`): examples and relation filter
  counts are prepared before the first request.

`GET /api/status` reports the query slots in use and waiting, recent response
times and coarse host load; the explorer header shows it.

These limits were chosen by load testing on 9 October 2026
(`docs/reports/nicesrv-explorer-capacity-20261009.md`). They serve about 105
explorer interactions a minute with the slowest 5% under about 5 s, or roughly
30–55 people using the explorer at the same time. The previous limits (6 GB,
4 CPUs, two queries of 2,500 MB) are kept on nicesrv as
`serving.env.bak-before-capacity-20261009`. Rerun the test after changing them:

```sh
uv run python scripts/explorer_load_test.py --users 8,12,16,24,32 \
  --abort-step-p95 15 --out load.json
```

## Rollback

Every deployed commit keeps its images (`omnipath-api:serving-<commit>`,
`omnipath-web:serving-<commit>`). To roll back, set `API_IMAGE` and `WEB_IMAGE`
in `serving.env` to an earlier tag and run `make serving-up` with the same
`COMPOSE_ENV` and `COMPOSE_FILES` as above. Running `deploy.sh` again returns to
the branch head.

Backups of the configuration before the switch to the normalized tables are
`deploy.sh.bak-normalized-tables` and `serving.env.bak-tables`. Earlier
deployment roots remain under `/root/projects/omnipath-releases/`. Keep
PostgreSQL, artifacts and other apps.

## HTTPS PostgreSQL input

The data host serves releases and tables with byte ranges, so PostgreSQL can
load a release without copying files:

```sh
make load-postgres \
  RELEASE_MANIFEST=https://data.omnipath-metabo-dev.schaul.click/releases/2026.10.9.json \
  DATA_ROOT=https://data.omnipath-metabo-dev.schaul.click \
  POSTGRES_SCHEMA=an_explicit_new_schema
```

Set `OMNIPATH_DATABASE_URL` for the intended database first. Use the normal user
environment, including `HOME`, so DuckDB can load its HTTP extension. Local
DuckDB staging is temporary; no persistent Parquet download cache is created.
Release `2026.10.9` was loaded this way into dev5 on beauty in 1 h 47 min.

For a large release, copy it first and load the copy. Release `2026.10.14.2` (171
resources) failed twice over HTTPS on beauty: DuckDB reported "Could not resolve
hostname" in the projection staging after thousands of good requests, retries did
not help. The copy (14 GB, checksums checked, files the host does not publish
skipped) took 3.5 minutes:

```sh
uv run python scripts/mirror_release.py 2026.10.14.2 /tmp/release \
  --base-url https://data.omnipath-metabo-dev.schaul.click
make load-postgres \
  RELEASE_MANIFEST=/tmp/release/releases/2026.10.14.2.json \
  DATA_ROOT=/tmp/release \
  POSTGRES_SCHEMA=an_explicit_new_schema
```

The load defers the constraints of all resource partitions in one transaction,
about 60 locks per resource; the loader checks before the projection that
PostgreSQL's lock table can hold them. The default `max_locks_per_transaction` of
64 held 46 resources but not 171: set it higher (the build compose uses 512, dev5
runs with 4096) and restart PostgreSQL. Release `2026.10.14.2` loaded into dev5
from a local copy in 2 h 40 min: 210 GB, 5.2 million entities, 27.9 million
relations.
