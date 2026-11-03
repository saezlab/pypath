# Build performance tuning (Postgres memory & parallelism)

The build runs a handful of heavy one-off statements — large sorts and
aggregations for the derived tables and network views, the RDKit structure
substrate parse (omnipath-metabo), and bulk index builds. These benefit a lot
from extra working memory and intra-query parallelism.

## Why it's session-level, not global config

The Postgres instance that the build writes to is often the **same instance
that serves the web API and web app**. Build and serving have opposite needs:

- the build wants large `work_mem` / `maintenance_work_mem` and many parallel
  workers for a few statements;
- serving wants small per-connection footprints across many concurrent
  connections.

So the build sets these as **session GUCs on its own connection** (they vanish
when the build finishes) and never writes them to the global config. Keep
`shared_buffers`, the global `work_mem`, etc. modest for serving; a small bump
to `effective_cache_size` (a planner hint, no allocation) is the only global
change worth making for a shared instance.

## The knobs

Both `omnipath-build` (derive, network-views, index builds) and
`omnipath-metabo` (RDKit substrate build) read the same environment variables.
Set one to an empty string to leave that GUC at the server default.

| Environment variable | GUC | Default | Notes |
|---|---|---|---|
| `OMNIPATH_BUILD_WORK_MEM` | `work_mem` | `512MB` | Per sort/hash node; a heavy query may use a few × this across parallel workers. |
| `OMNIPATH_BUILD_MAINTENANCE_WORK_MEM` | `maintenance_work_mem` | `2GB` | Index/GiST builds; × `max_parallel_maintenance_workers`. |
| `OMNIPATH_BUILD_MAX_PARALLEL_WORKERS_PER_GATHER` | `max_parallel_workers_per_gather` | `6` (build) / `8` (substrate parse) | Capped by the server's `max_parallel_workers` / `max_worker_processes`. |
| `OMNIPATH_BUILD_MAX_PARALLEL_MAINTENANCE_WORKERS` | `max_parallel_maintenance_workers` | `4` | Parallel index builds. |

Server-side prerequisites for the parallel settings to take full effect (these
are global and need a restart, so set them on a build host or a dedicated build
instance, **not** on a small serving instance):

```
max_worker_processes = 16     # must be >= the parallel worker counts you want
max_parallel_workers = 12
```

## The utils / resolver-source Postgres (READ side) — tune it too

There are **two** Postgres containers in a build, and this doc historically only
covered the one the build *writes* to. The build also **reads** its identifier
resolution from a **separate `omnipath-utils` Postgres** over DuckDB `ATTACH`
(`OMNIPATH_BUILD_UTILS_PG_URL`, e.g. `utils2` on `:5102`). The heavy resolver
work — the `resolver_gene` (~52M rows), `resolver_protein` (~300M), and
`resolver_chemical` (~133M) **materialized views**, their rebuilds, and the
full-scan COPYs the build pulls from them — all execute on **that** instance,
using **its** config. The build cannot session-tune it (it doesn't own those
connections), so this PG must be tuned at the **container level**.

Left at the Postgres default (`shared_buffers=128MB`, `work_mem=4MB`), matview
rebuilds crawl single-core and STITCH canonicalization stalls on the chemical
COPY. **Tune it in its compose `command:`, not by ad-hoc `ALTER SYSTEM`** (which
is lost on recreate). The `utils2` compose
(`~/instances/utils2/docker-compose.yml`) now carries the canonical values below.

**`max_connections` scales with `LOAD_JOBS`.** The build's `LOAD_JOBS` (a.k.a.
`--stage-jobs`) parallel staging workers each read the resolver from this PG —
one DuckDB `ATTACH` connection plus a handful of short-lived psycopg2
connections for the keyed lookups (`_fetch_live_utils_rows_for_keys`). Peak
concurrent connections roughly track `LOAD_JOBS × ~4`, so the utils PG's
`max_connections` must comfortably exceed that. The default 100 is fine for
`LOAD_JOBS=16` but overflows at `LOAD_JOBS=32`; keep `max_connections` at
`400` (set below) or raise it alongside any further `LOAD_JOBS` increase.

## Canonical values on the beauty host (503 GB RAM, 64 cores)

Both containers set these via their compose `command:` (build PG:
`docker-compose.postgres18.yml`, env-overridable; utils PG: its own compose):

| GUC | Build PG (writes) | Utils/resolver PG (reads) | Why |
|---|---|---|---|
| `shared_buffers` | 4 GB | **32 GB** | utils holds the multi-hundred-M-row resolver matviews; 128MB default missed the pool entirely |
| `effective_cache_size` | 128 GB | 300 GB | planner hint (no allocation); ~340 GB is in OS page cache |
| `work_mem` | 128 MB (session ↑) | 64 MB | build raises its own per session; utils rebuilds spill less |
| `maintenance_work_mem` | 2 GB | 2 GB | faster index builds on the matviews |
| `max_worker_processes` | 24 | 24 | **the real parallel cap** — restart-only, so must be container-level |
| `max_parallel_workers` | 16 | 24 | |
| `max_parallel_workers_per_gather` | 6 (session) | 8 | |
| `effective_io_concurrency` | 200 | 200 | NVMe |
| `wal_level` | (WAL tuned) | `minimal` | utils build DB is rebuildable → skip WAL |
| `max_connections` | 100 (default) | **400** | each `LOAD_JOBS` stage worker holds a DuckDB `ATTACH` + serial psycopg2 keyed-lookup connections to the utils PG; at `LOAD_JOBS=32` the default 100 overflowed (`FATAL: sorry, too many clients already` → ~50 datasets fail to stage). 400 leaves ~12/worker of headroom (steady-state peak was ~46). Restart-only. |
| `shm_size` (container, not a GUC) | **1 GB** | **16 GB** | `/dev/shm` for parallel workers + large hash/sort temp. Docker's **64 MB default** is far too small (a single parallel segment is 256 MB) and silently applies if omitted → `could not resize shared memory segment … No space left on device`. Set via `shm_size:` on the compose service (restart-only). Web/serving PGs also want ≥1 GB. A heavy *parallel* analysis query can still exceed 1 GB → `SET max_parallel_workers_per_gather=0` in the query instead. |

`shared_buffers` and `max_worker_processes` are **restart-only** — after changing
them, recreate the container (`docker compose up -d`), don't just `ALTER SYSTEM`.
On beauty, all instances share one docker cgroup (`MemoryMax=200GB`), so keep the
sum of `shared_buffers` across concurrently-running build+utils PGs well under that.

## Memory budget / sizing

A single build connection peaks at roughly **`work_mem` × (a few parallel
nodes) + `maintenance_work_mem` × `max_parallel_maintenance_workers`** — with
the defaults, ~8–10 GB at peak. Size the defaults to your host:

- **Lab workstation (beauty):** all Postgres instances + web services run under
  one `docker.service` cgroup with a **hard `MemoryMax = 200 GB`**
  (`MemoryHigh = 150 GB` soft) — see `saez-nixos`
  `modules/services/resource-limits.nix`. One build at ~10 GB is comfortable;
  if several instances build concurrently, keep the sum well under 200 GB.
- **Smaller hosts:** lower `OMNIPATH_BUILD_WORK_MEM` (e.g. `128MB`) and
  `OMNIPATH_BUILD_MAINTENANCE_WORK_MEM` (e.g. `512MB`), and reduce the parallel
  worker counts to match available cores/RAM.

Example (smaller host):

```bash
export OMNIPATH_BUILD_WORK_MEM=128MB
export OMNIPATH_BUILD_MAINTENANCE_WORK_MEM=512MB
export OMNIPATH_BUILD_MAX_PARALLEL_WORKERS_PER_GATHER=2
export OMNIPATH_BUILD_MAX_PARALLEL_MAINTENANCE_WORKERS=2
make all DERIVE=1 DATABASE_URL=... DATA_ROOT=../data PUBCHEM_SHARDS=1
```
