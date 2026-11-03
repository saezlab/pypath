# Serve existing OmniPath artifacts

The serving stack has three services: API, web and an allowlisted HTTPS file
backend. It reads already built data; starting it performs no resource build or
PostgreSQL load. Package Dockerfiles build code images only.

## Local containers

```sh
cp deploy/serving.env.example deploy/serving.env
# Set DATA_DIR to an existing absolute data directory.
make serving-build COMPOSE_ENV=deploy/serving.env
make serving-up COMPOSE_ENV=deploy/serving.env
make serving-status COMPOSE_ENV=deploy/serving.env
make serving-logs COMPOSE_ENV=deploy/serving.env
make serving-stop COMPOSE_ENV=deploy/serving.env
```

API, web and files bind to loopback ports 8085, 8082 and 8080. Port variables can
be changed in the environment file. The API and file service mount the artifact
directory at `/data`, read-only.
Nginx exposes only complete versioned resource files, named releases and taxonomy
Parquets. It hides raw caches, the full reference library and internal indexes.
Public files must be readable by its worker (normally files 0644, directories
0755); do not change permissions recursively over the whole data directory.

## HTTPS and a parallel cutover

On an existing Traefik/Dokploy host, set PUBLIC_DOMAIN, DATA_DOMAIN and
PUBLIC_DATA_URL, then add the routing override:

```sh
make serving-up COMPOSE_ENV=deploy/serving.env \
  COMPOSE_FILES='-f compose.serving.yaml -f deploy/compose.traefik.yaml'
```

The proxy network must already exist. Set PROXY_NETWORK and TLS_RESOLVER for the
host's configuration. Router names use COMPOSE_PROJECT_NAME, so independent
stacks can coexist.

To replace a prototype deployment without rebuilding data:

1. Use a distinct COMPOSE_PROJECT_NAME, distinct loopback ports and versioned
   API_IMAGE/WEB_IMAGE tags. Point DATA_DIR at the existing artifact directory.
2. Build and start the candidate with the base Compose file only. Check health,
   resource inventory, representative queries, explorer pages and range downloads
   through its loopback ports. No public routing changes yet.
3. Add the HTTPS override. ROUTE_PRIORITY must exceed the existing application's
   route priority. Confirm public requests reach the new healthy containers.
4. Stop only the superseded API/web/file containers after the checks pass.
   Retain their images and configuration until rollback is no longer needed.

Rollback: start the retained old containers, remove public routing from the new
stack by starting it with the base Compose file, and confirm the old endpoints.
Never run a broad host-level Docker prune or remove data/volumes during cutover.
Production PostgreSQL is independent and does not participate in this deployment.

Resource publication is separate from a code cutover. Already built new versions
can be copied atomically into `resources/<source>/<version>/` and a pinned release
published afterward. Keep older versions. Nginx forbids symlinks; copying the
entire cache/reference tree is unnecessary. For an existing PostgreSQL snapshot,
retain its original release JSON bytes if publishing the same manifest over HTTPS.

The current nicesrv deployment and exact rollback commands are recorded in
[nicesrv.md](nicesrv.md).

## Optional build worker

The worker is excluded from normal startup. `make worker-up` explicitly builds
and starts it with a writable data mount. Before using it, provision an identity
library (see [the identity layer guide](../packages/build/REFERENCE.md)), point
`OMNIPATH_LIBRARY_DIR` at it, and configure build limits/queues according
to the [API guide](../packages/api/README.md) and
[resource orchestration guide](../packages/build/ORCHESTRATION.md).
Stopping the worker does not stop read-only serving.

The existing PostgreSQL promotion scripts in this directory are a separate
workflow; serving commands never invoke them.

## Query and database budgets

The API pool bounds simultaneous queries and retained DuckDB connections.
`API_QUERY_CONCURRENCY`, `API_QUERY_THREADS` and `API_QUERY_MEMORY` set its
per-process query budget; `API_MEMORY_LIMIT` and `API_CPUS` bound the container.
The default two databases each have a 1 GB DuckDB limit inside a 4 GB container.
Additional Arrow/Python buffers need headroom. These are conservative starting
values. `scripts/explorer_load_test.py` measures how many explorer users a
deployment supports; `docs/reports/nicesrv-explorer-capacity-20261009.md` has the
results for nicesrv.

A development PostgreSQL instance is a separate explicit operation:
copy `postgres-dev.env.example` to a private ignored `.env` file, choose a
password, then use root `make postgres-up POSTGRES_ENV=deploy/postgres.env`.
Its port binds only to loopback. `postgres-server.env.example` contains larger
starting settings for a dedicated server, without changing any live deployment.
