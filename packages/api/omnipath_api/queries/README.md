# The API's queries

DuckDB queries over OmniPath's Parquet tables in a serving root (including a read-only
mounted one), for the REST API and the explorer: entity search and details, grouped
search, relation search, facets, ontology browsing, statistics and Parquet exports.
The API's routes call `ParquetQueries`; nothing else does.

The clients (`packages/query`) do not run this code. Their search, relation export and
entity details are SQL macros (`packages/query/sql`) in the database each client
writes; the Python client's tests compare their answers with these queries. The other
queries here can move into SQL the same way, one at a time.

```python
from omnipath_api.queries import ParquetQueries

queries = ParquetQueries("/path/to/data")
with queries.query_scope():
    found = queries.search_entities_api("TP53", limit=5)["entities"]
queries.close()
```

## Data directory

`resources/<resource>/<version>/` holds one resource version: the published tables and
the serving tables (`omnipath_core.schema`). A version is served when all of them are
present; raw evidence payloads (`evidence_payloads.parquet`) are optional. Optional
`releases/<version>.json` manifests pin resource versions, `query_policy.json` excludes
resources from queries, and `references/taxonomy` names taxa.

## Modules

- **`parquet` / `store`**: `ParquetQueries`, inventory generations, request scopes,
  query caches and configured DuckDB connections. `memory_limit` and `threads` set each
  query connection; `OMNIPATH_QUERY_CONCURRENCY` sets how many run at once.
- **`services`**: domain query services for entities, relations, ontology, evidence,
  facets, resources, statistics, selection and export.
- **`filters`**: the filter and cursor contracts the queries accept.
- **`annotations` / `shape`**: presentation adapters for compiled identifiers,
  annotations, entities and relations.

## Tests

```bash
uv run pytest packages/api/tests -q
```
