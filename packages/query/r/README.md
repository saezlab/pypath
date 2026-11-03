# omnipath.client (R)

Install OmniPath resources on your machine and get a DuckDB database of them. Query
it with any SQL, or ask the questions the OmniPath website answers and get the same
answers: the website's search, relation export and entity details are SQL macros in
the database. Nothing is queried over the network; the API is contacted only while
installing.

This is the R counterpart of the Python package `omnipath_client`
(`packages/query/python`). Both read and write the same data directory and run the same SQL,
so a directory installed from R can be queried from Python and the other way round.
Neither needs Python, the API or its queries to query.

```r
library(omnipath.client)

op_init(topic = "signaling")                        # once: download and write the database

con <- op_connect()
op_sql(con, "SELECT * FROM signor.relation LIMIT 5")  # your SQL
op_search(con, "TP53")                                 # the website's search
op_relations(con, "TP53", sources = "signor")          # the website's relation export
op_entity(con, op_search(con, "TP53")$entity_key[1])   # one entity's details
op_disconnect(con)
```

## Install

From this checkout:

```r
install.packages(c("DBI", "duckdb", "httr2", "jsonlite"))
install.packages("packages/query/r", repos = NULL, type = "source")
```

Requires R 4.1 or newer and the R package `duckdb` 1.5 or newer. The name
`omnipath.client` keeps it apart from the Bioconductor package OmnipathR, which reads
the legacy web service.

## Initialize

```r
op_init()                              # every resource of the latest release
op_init(topic = "cell_cell_communication")
op_init(c("signor", "intact"))         # by name
op_available()                         # what a release offers: versions, sizes, tags
```

Selection, updates, releases, raw evidence and interruptions work as in the Python
client: resources are chosen by name and by the catalog's tags (`topic`, `content`,
`molecule`, `knowledge_origin`); a later `op_init()` moves every installed resource to
the requested release and downloads only what is missing; `include_evidence = TRUE`
adds raw evidence payloads; an interrupted run keeps the resources it finished.

The data directory is the same as the Python client's: the directory given
(`data_dir`), else `OMNIPATH_DATA_DIR`, else the per-user default
(`~/Library/Application Support/omnipath` on macOS, `~/.local/share/omnipath` on
Linux, `%LOCALAPPDATA%\omnipath\omnipath` on Windows). `op_data_directory()` returns it.

The database is written with **native DuckDB tables** by default (about twice the
files' size in addition); `materialize = FALSE` writes views over the files instead
(small, but the prebuilt queries slow down with many resources). It is written in
DuckDB 1.5's storage format, so clients with an older DuckDB than R's can read it.

## SQL

`op_connect()` returns a read-only DBI connection, so everything DBI and dplyr offer
works on it:

```r
con <- op_connect()
op_sql(con, "SELECT label, _resource FROM entity WHERE identifier = ?", list("P04637"))
DBI::dbGetQuery(con, "SELECT predicate, count(*) FROM relation GROUP BY 1")
dplyr::tbl(con, "relation") |> dplyr::filter(sign == -1) |> dplyr::count(predicate)
```

The tables are those of the Python client's README: `entity`, `relation` and the
others over every installed resource with `_resource` and `_resource_version`, one
resource's tables as `signor.relation`, and `resources`.

## Prebuilt queries

```r
op_search(con, "gluc", entity_types = "chemical_entity", limit = 50)
op_relations(con, "gene:entrez:7157", predicates = "affects", taxonomy_ids = "9606")
op_relations(con, protein_entity_keys = key, isoform_identifiers = "uniprot:P04637-2")
op_entity(con, key, simplify = FALSE)   # the API's JSON as nested lists
```

- `op_search()` and `op_relations()` return data frames; `lazy = TRUE` returns a
  dplyr table instead (needs dplyr and dbplyr). Annotations, evidence and product
  records are list columns.
- Filters are the explorer's, by the same names as in Python; a character vector is a
  list of values. Unknown filters and Biolink names are errors naming them; a source
  that is not installed is an error listing the installed ones.
- The answers are the Python client's and the hosted API's; see the Python client's
  README for the differences from the hosted API.

From SQL, the macros are `op_search(query, ...)`, `op_relations(entities := [...], ...)`
and `op_entity(key)`. A connection that runs them sets
`SET dynamic_or_filter_threshold = 20000` first, as `op_connect()` does.

## Validation

```bash
Rscript -e 'testthat::test_local("packages/query/r")'
uv run pytest packages/query/python/tests/test_r_client.py -q
```

The tests replace the package's two HTTP functions with readers of
`tests/testthat/fixtures` (the Python client tests' resources, written by
`uv run python packages/query/r/tests/make_fixtures.py`). `test_r_client.py` installs
the same resources with both clients and checks that each reads the other's database
without writing it again and gives the same answers. The SQL files in `inst/sql/` are
copied from `packages/query/sql/` by
`packages/query/scripts/generate_contract.py`; do not edit them here.
