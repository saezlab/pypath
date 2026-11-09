# The API's queries

DuckDB queries over the release databases of a serving root (including a read-only
mounted one), for the REST API and the explorer: entity search and details, grouped
search, relation search, facets, ontology browsing, statistics and Parquet exports.
The API's routes call `ReleaseQueries`; nothing else does.

Each published release the API serves has one DuckDB database,
`serving/<version>.duckdb`: every resource's tables in one native table each, with
global ids, a key table and indexes. It is the database format the clients
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
- `query_policy.json` excludes resources from queries; `references/taxonomy` names taxa.

The databases are snapshots: a resource rebuilt in place, or a new version, is served
once a release pins it and its database is built.

## Reading rows

Rows are read by id through the global ids (`_rid`, `_eid`: the resource's index times
2^32 plus its own id), written into the SQL as literal lists, in a MATERIALIZED step of
their own. DuckDB uses an ART index only then: with any other condition in the same
`WHERE`, or with ids from a subquery or join, it scans (6–7 s instead of 0.3 s for a
gene's relations in a full release). Relation keys are looked up in
`relation_key_lookup`, which is sorted by key.

## Modules

- **`release` / `store`**: `ReleaseQueries`, the served databases, request scopes,
  query caches and configured DuckDB connections. `memory_limit` and `threads` set each
  query connection (a full release needs about 3.5 GB per query with its indexes);
  `OMNIPATH_QUERY_CONCURRENCY` sets how many run at once and
  `OMNIPATH_DUCKDB_TEMP_DIRECTORY` where they spill.
- **`services`**: domain query services for entities, relations, ontology, evidence,
  facets, resources, statistics, selection and export.
- **`filters`**: the filter and cursor contracts the queries accept.
- **`annotations` / `shape`**: presentation adapters for compiled identifiers,
  annotations, entities and relations.

## Tests

```bash
uv run pytest packages/api/tests -q
```

`tests/served.py` builds the database of every published release of a test's data root
(publishing its newest resources first when nothing is published).
