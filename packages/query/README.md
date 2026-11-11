# The query layer

How OmniPath data is queried outside the hosted API: SQL in a local DuckDB database,
with a Python and an R client that install the data and write that database.

| Folder | What it is |
|---|---|
| [`sql/`](sql) | Every OmniPath query as a step file (`sql/queries/<name>.sql`), the fragments they share (`sql/fragments/`), the shaping macros (`shape.sql`) and `contract.json`: what every client needs to write and query the database. The one place to edit them: the API and both clients run these files. |
| [`python/`](python/README.md) | `omnipath_client` (pip `omnipath-client`): `init`, `connect`, `sql`, `query`, `search`, `relations`, `entity` |
| [`r/`](r/README.md) | `omnipath.client`: `op_init`, `op_connect`, `op_sql`, `op_query`, `op_search`, `op_relations`, `op_entity` |
| [`scripts/`](scripts) | `generate_contract.py` (the contract from `packages/core`, and the copies of `sql/` in both clients) and `compare_with_api.py` (the clients' answers against the API's on a real installation) |

Both clients have two modes: API mode (the default) asks the hosted API; local mode
downloads the same files into the same data directory and writes the same database,
which is also the database the API serves a release from (`omnipath_client.database`,
`build_release_database` on the server): published tables as native DuckDB tables,
merged entities and relations numbered at build (`entity_n`, `relation_n`) with the
resources that state each, every evidence record with its resource, lookup tables
derived with the build's queries, indexes, the query tables of the contract and the
release's taxon names (`op_taxonomy`). There are no per-resource tables or schemas. Either client reads a directory the other installed; neither
needs the API or Biolink at run time in local mode.

The hosted API (`packages/api`) runs the same step files on the same database, so a
query answers the same in the API, in Python and in R. API mode asks
`POST /query/{name}` and gets the local answer's table as Parquet.

## Step files

A query is a sequence of SQL statements, each a step whose result later steps use: a
temporary table (`__name`), or a list or value fetched by the runner and written into
the later SQL as a literal. Literals are what let DuckDB read rows through the ART
indexes on the entity and relation numbers (`entity_n`, `relation_n`) and skip row
groups: an index read is a bare `entity_n IN (<literal numbers>)` with no other
condition in the same `WHERE`; numbers from a
subquery or join make a scan (seconds instead of milliseconds on a full release).

```sql
-- query: entity_core
-- param: entity_key VARCHAR
-- param: resources VARCHAR[] = null
-- output: json answer
-- example: {"entity_key": "..."}
-- step: ids AS list
SELECT entity_n FROM entity WHERE entity_key = {{entity_key}}
-- step: rows AS table
SELECT * FROM entity WHERE entity_n IN ({{ids*}})
-- step: answer AS value
SELECT to_json(...) FROM __rows
```

| Construct | Effect |
|---|---|
| `-- param: name TYPE [= json default]` | Validated and written into the SQL as a literal (never bound). `resources` is set by the API from the request's selection. |
| `-- output: table <step>` / `json <step>` | What `query(name, output=...)` returns: a table step, or a value step holding JSON. |
| `-- step: x AS table\|list\|value\|rows [when <expr>] [cache]` | A temporary table, the first column as a list, one value, or rows for `{{table:x}}`. `when` skips the step (a skipped list is empty); `cache` keeps its result per release and parameters. |
| `-- repeat: v = <init> next <expr> until <value step> max N` … `-- end` | The only loop (search paging, the ontology walk). |
| `-- steps: fragment`, `{{include:fragment table}}` | Shared step sequences and SQL pieces from `sql/fragments/`. |
| `{{p}}`, `{{p*}}`, `{{#p}}`, `{{p[o:l]}}` | A literal; a list's items for `IN (…)` (an empty subquery when empty); its length; a slice. |
| `error('invalid: …')` | A `ValueError` in Python, an R error, HTTP 400 in the API. |

The runners are `python/omnipath_client/steps.py` and `r/R/steps.R`; each connection
gets the contract's settings (`apply_settings`) and the `shape.sql` macros as
temporary macros. Comments are dropped before a step runs.

## Changing the SQL

1. Edit `sql/queries/*.sql`, `sql/fragments/*.sql` or `sql/shape.sql`, or the tables
   or Biolink in `packages/core`. A new query needs at least one `-- example:`.
2. `uv run python packages/query/scripts/generate_contract.py` writes `contract.json`
   and copies `sql/` into `python/omnipath_client/sql/` and `r/inst/sql/`.
3. `uv run pytest packages/query/python/tests -q` (with R installed, this includes the
   R client against the Python client) and
   `Rscript -e 'testthat::test_local("packages/query/r")'`.

The SQL runs on DuckDB 1.5 and 2.0: use `lambda x: ...` lambdas, not `x -> ...`.
