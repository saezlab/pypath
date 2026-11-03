# The query layer

How OmniPath data is queried outside the hosted API: SQL in a local DuckDB database,
with a Python and an R client that install the data and write that database.

| Folder | What it is |
|---|---|
| [`sql/`](sql) | The macros `op_search`, `op_relations` and `op_entity` (the explorer's search, relation export and entity details) and `contract.json`: what every client needs to write and query the database. The one place to edit them. |
| [`python/`](python/README.md) | `omnipath_client` (pip `omnipath-client`): `init`, `connect`, `sql`, `search`, `relations`, `entity` |
| [`r/`](r/README.md) | `omnipath.client`: `op_init`, `op_connect`, `op_sql`, `op_search`, `op_relations`, `op_entity` |
| [`scripts/`](scripts) | `generate_contract.py` (the contract from `packages/core`, and the copies of `sql/` in both clients) and `compare_with_api.py` (the clients' answers against the API's on a real installation) |

Both clients download the same files into the same data directory and write the same
database: published tables as native DuckDB tables (or views), lookup tables derived
with the build's queries, global ids, the installation's lookup macros and the macros
of `sql/`. Either client reads a directory the other installed. Any DuckDB client can
run the SQL; neither client needs the API, its Python queries or Biolink at run time.

The hosted API (`packages/api`) answers the same three questions with its Python
queries (`omnipath_api.queries`), and more besides. The Python client's tests compare
the macros with them, so the two give the same answers.

## Changing the SQL

1. Edit `sql/*.sql`, or the tables or Biolink in `packages/core`.
2. `uv run python packages/query/scripts/generate_contract.py` writes `contract.json`
   and copies `sql/` into `python/omnipath_client/sql/` and `r/inst/sql/`.
3. `uv run pytest packages/query/python/tests -q` (with R installed, this includes the
   R client against the Python client) and
   `Rscript -e 'testthat::test_local("packages/query/r")'`.

The SQL runs on DuckDB 1.5 and 2.0: use `lambda x: ...` lambdas, not `x -> ...`.
