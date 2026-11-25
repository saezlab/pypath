# nicesrv serving deployment

The API, explorer and HTTPS file service follow branch `monorepo` of saezlab/pypath:
`deploy.sh` checks out the branch head, builds images tagged with its short
commit (`serving-mono-<commit>`) and restarts the stack. The deployed commit is
recorded in `source-commit.txt`.

- Explorer: https://omnipath-metabo-dev.schaul.click/explore
- Files: https://data.omnipath-metabo-dev.schaul.click
- Deployment root: `/root/projects/omnipath-releases/20261007-serving-tables`
  (`deploy.sh`, `serving.env`, `source-commit.txt` and one log per deploy).
- Source checkout used by `deploy.sh`: `/root/projects/omnipath-monorepo`
  (detached at `origin/monorepo`).
- Artifact root: `/root/data/release-20261007-native` (release `2026.10.25.1`, 173
  resources), mounted read-only at `/data`.
- Project: `omnipath-serving-20261007-full`; loopback API/web/files ports:
  `8295` / `8292` / `8290`.

## Deploy, inspect or restart

```sh
ssh -T -o BatchMode=yes -o ForwardAgent=no -o ConnectTimeout=15 nicesrv
deployment_root=/root/projects/omnipath-releases/20261007-serving-tables
# Deploy the current head of monorepo (push first).
"$deployment_root/deploy.sh"
# Inspect the running stack.
cd /root/projects/omnipath-monorepo
make serving-status COMPOSE_ENV="$deployment_root/serving.env" \
  COMPOSE_FILES='-f packages/api/deploy/compose.serving.yaml -f packages/api/deploy/compose.traefik.yaml'
make serving-logs COMPOSE_ENV="$deployment_root/serving.env" \
  COMPOSE_FILES='-f packages/api/deploy/compose.serving.yaml -f packages/api/deploy/compose.traefik.yaml'
```

A deploy builds and starts serving only; it runs no resource build, identity
build or PostgreSQL load. Changing the served data means changing `DATA_DIR` in
`serving.env`.

## Effective API limits

Set in `serving.env`:

- Container: 13 GB RAM and 7 CPUs (`API_MEMORY_LIMIT`, `API_CPUS`), leaving 1 CPU and
  about 3 GB of the host's 8 CPUs and 16 GB for the other apps.
- Query pool: sized by the API from those limits (no `API_QUERY_*` settings): seven
  query slots, one DuckDB thread for the shared database (each query runs on its
  caller's thread), and half the memory, 6.5 GB, for DuckDB.
- Capacity wait: up to 30 seconds, then HTTP 503 with a retry hint.
- Cache warming at startup (`API_WARM_CACHE=true`): examples and relation filter
  counts are prepared before the first request.

`GET /api/status` reports the query slots in use and waiting, recent response
times and coarse host load; the explorer header shows it. Its `memory` counts the
container's file cache too, so it sits near 1 under load without any pressure.

Measured on 25 November 2026 (commit 7d5361785, database format 13) with
`packages/web/scripts/api_load_test.py` run on the host, cold answer cache, levels of
1, 8, 32 and 64 clients without think time:

| Workload | Requests/s (8–64 clients) | p95 at 8 / 64 clients | API CPU per request |
|---|---|---|---|
| Entity details | 39–40 | 0.38 / 2.4 s | 125 ms |
| First relations page | 23–25 | 0.43 / 3.3 s | 190 ms |
| Mix (entity, relations, resolve, search, relationships) | 23–27 | 0.89 / 3.8 s | 185–225 ms |

No errors at any level. Throughput levels off at 8 clients; more clients only wait.
This is about half of what the same code and database do in the same container
limits on a host with plenty of RAM (51 mixed requests/s at about 100 ms CPU):
the 18.5 GB database does not fit in nicesrv's page cache, so under load the host
reads 140–300 MB/s from disk and swaps. More host memory is the lever, not more
query slots. The settings before this deploy (format 12, six pinned slots of
1,750 MB) are kept as `serving.env.bak-before-mono-7d53617`. Rerun after changing limits:

```sh
python3 packages/web/scripts/api_load_test.py sample 2026.10.25.1.duckdb params.json
python3 packages/web/scripts/api_load_test.py run --port 8295 --params params.json \
  --cgroup /sys/fs/cgroup/system.slice/docker-<api container id>.scope \
  --workload mix --clients 1,8,32,64 --processes 2 --out mix.jsonl
```

`packages/web/scripts/explorer_load_test.py` replays explorer sessions instead; the
October results are in `packages/api/docs/reports/nicesrv-explorer-capacity-20261009.md`.

## Release databases

The API serves each published release from its DuckDB database,
`serving/<version>.duckdb` in the artifact root. The default view (`latest`) is the
newest release that has one. A request for an older release whose database was
pruned answers 410; its files stay downloadable. The API mounts the artifact root
read-only, so the database is written on the host:

The database format must match the deployed code: a deploy that changes the format
needs the database rebuilt first. Build it with the new API image into a scratch
root whose `releases`, `resources` and `references` link to the artifact root, then
hard-link it into `serving/` right before `make serving-up`:

```sh
S=/root/data/scratch/<build>; D=/root/data/release-20261007-native
mkdir -p $S/root/serving; for x in releases resources references; do ln -s $D/$x $S/root/$x; done
docker run --rm --user 0 --cpus 3 --cpu-shares 256 --memory 8g -v $D:$D:ro -v $S:$S \
  --entrypoint python omnipath-api:serving-mono-<commit> -c \
  "from omnipath_api.release_database import build_served_database as b; b('$S/root', '2026.10.25.1', keep=None)"
```

- **Time and size:** a full release (173 resources, 2026.10.25.1) takes about 25 minutes
  at 3 threads and 5 GB (the defaults; less memory runs out in the merges), and the
  format-13 database is 18.5 GB. Format 13 was built this way on 25 November 2026
  (`/root/data/scratch/format13-20261125`); the format-12 database is kept as a hard
  link in `/root/data/scratch/format12-backup-20261125/` for rollback with the
  `serving-mono-01b7dee` images.
- **Disk:** the command refuses to start without the last database's size plus 20 GiB
  free. Three kept databases need about 60 GB.
- **Pickup:** the running API sees a new database as soon as it is in place: no restart
  and no deploy.
- **Query memory:** all slots share one DuckDB database with half the container's
  memory. DuckDB spills to the container's `/tmp`, or to
  `OMNIPATH_DUCKDB_TEMP_DIRECTORY`.

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
uv run python packages/api/scripts/mirror_release.py 2026.10.14.2 /tmp/release \
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
