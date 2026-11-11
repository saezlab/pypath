# The API's queries

DuckDB queries over the release databases of a serving root (including a read-only
mounted one), for the REST API and the explorer: entity search and details, grouped
search, relation search, facets, ontology browsing, statistics and Parquet exports.
The API's routes call `ReleaseQueries`; nothing else does.

Each published release the API serves has one DuckDB database,
`serving/<version>.duckdb`: merged entities and relations (one row each, numbered at
build, with the resources that state them), every evidence record with its resource,
lookup tables and indexes. It is the database format the clients
(`packages/query`) write locally too, built by the same code
(`omnipath_client.database`). A request reads one release, by default the newest one
with a database.

```python
from omnipath_api.queries import ReleaseQueries

queries = ReleaseQueries("/path/to/data")
with queries.query_scope(), queries.release_scope("latest"):
    found = queries.search_entities_api("TP53", limit=5)["entities"]
queries.close()
```

## Data directory

- `resources/<resource>/<version>/`: one resource version's published and serving
  tables (`omnipath_core.schema`), as built. Downloads and raw evidence payloads
  (`evidence_payloads.parquet`) are read from here.
- `releases/<version>.json`: published releases, each pinning resource versions.
- `serving/<version>.duckdb` and `serving/<version>.json`: the database of a release and
  how it was built, written by `omnipath-api build-release-db <version>`
  (`omnipath_api.release_database`). The newest few are kept; a request for a release
  without a database gets 410 (its files stay downloadable).
- `references/taxonomy` names taxa.

The databases are snapshots: a resource rebuilt in place, or a new version, is served
once a release pins it and its database is built.

## Queries

Every query is a step file of the clients' package (`packages/query/sql/queries/`,
described in `packages/query/README.md`): the API, the Python client and the R client
run the same SQL, and the answers' JSON is written by DuckDB. `ReleaseQueries.run_query`
runs one on the request's release, with the selected resources as its `resources` parameter, and caches the answers and the steps
marked `cache`. The services turn a request's fields into a query's parameters; the
generic route `POST /query/{name}` (`routers/queries.py`) runs any query by name, which
the clients' API mode uses.

## Modules

- **`release` / `store`**: `ReleaseQueries`, the served databases, request scopes,
  query caches and configured DuckDB connections. `memory_limit` and `threads` set each
  query connection (a full release needs about 3.5 GB per query with its indexes);
  `OMNIPATH_QUERY_CONCURRENCY` sets how many run at once and
  `OMNIPATH_DUCKDB_TEMP_DIRECTORY` where they spill.
- **`services`**: the routes' requests as step-file parameters (entities, relations,
  evidence, facets, ontology, selection, statistics, groups, exports), the resource
  catalog and downloads, and raw evidence payloads (read from the resources' files).
- **`filters`**: the filter and cursor contracts the requests accept.
- **`examples`**: the landing page's examples, kept as a file per release.

## Tests

```bash
uv run pytest packages/api/tests -q
```

`tests/served.py` builds the database of every published release of a test's data root
(publishing its newest resources first when nothing is published).
