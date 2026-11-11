# OmniPath client

Ask the questions the OmniPath website answers and get the same answers, in one of two
modes:

- **API mode** (the default): the hosted API answers the website's search, relation
  export and entity details. Nothing is downloaded but the answers.
- **Local mode**: install resources on your machine and get the DuckDB database the
  API serves a release from. The same three questions are answered from it, and any
  SQL from Python, R or the DuckDB shell runs on it too.

The R package `omnipath.client` (`packages/query/r`) has the same two modes and
installs and reads the same data directory.

```python
import omnipath_client as op

api = op.connect()                                     # API mode
api.search("TP53")                                     # the website's search
api.relations("TP53", sources="signor")                # the website's relation export

op.init(topic="signaling")                             # local mode: download, write the database
db = op.connect(op.data_directory())
db.sql("SELECT * FROM signor.relation LIMIT 5").df()   # your SQL
db.relations("TP53", sources="signor")                 # the same answer as above
```

`connect()` is in API mode unless it gets a data directory (or `OMNIPATH_DATA_DIR` is
set, or `mode="local"`). In API mode `release=` picks the release the API answers for
(default `latest`, the newest one it serves) and `api_url=` the API; `sql()` needs
local mode.

This package lives in `packages/query/python`. It does not implement the retired
PostgreSQL API.

## Install

From this checkout:

```bash
python -m pip install -e './packages/query/python'
# Optional DataFrame conversions:
python -m pip install -e './packages/query/python[pandas]'
```

Requires Python 3.11+. The client depends on DuckDB, httpx and platformdirs only. What
it shares with the build and the API (table names, the lookup tables' queries,
Biolink names) is in `contract.json`, generated from `omnipath-core`; the queries and
the contract are edited in `packages/query/sql/` and copied into this package (see
"Validation").
The distribution and import names are `omnipath-client` and `omnipath_client`, which
are also used by the colleagues'
[saezlab/omnipath-client](https://github.com/saezlab/omnipath-client). Use a separate
environment for this checkout; a bare `pip install omnipath-client` installs their
published package, not this version. The older `omnipath` package uses a different
namespace and can coexist.

## Initialize

```python
op.init(progress=print)                    # every resource, about 5.6 GB
op.init(topic="cell_cell_communication")   # 9 resources, about 20 MB
op.init(["signor", "intact"])              # by name
```

or from a shell:

```bash
omnipath-client list                       # resources, versions, sizes, tags
omnipath-client init --topic signaling
omnipath-client init signor intact
```

- **Selection.** Resources are chosen by name and by the catalog's tags in four
  dimensions: `topic`, `content`, `molecule` and `knowledge_origin`. Several tags of
  one dimension are alternatives; different dimensions must all match. Without a
  selection, `init` installs every resource. `available()` and `omnipath-client list`
  show what a release offers.
- **Location.** Data goes to the per-user data directory (`data_directory()` returns
  it), to `OMNIPATH_DATA_DIR` when that variable is set, or to `data_dir=...`.
  `connect()` reads the same place. The directory can be moved or copied as a whole.
- **Updates.** `init` only downloads versions that are missing. Running it again
  moves every installed resource to the requested release's versions and removes the
  versions they replace, so a directory always holds one release.
- **Releases.** `init(release="2026.9.6.2")` installs the versions that release pins;
  the default `latest` takes the newest version of each resource.
- **Raw evidence** payloads are downloaded only with `include_evidence=True`
  (`--evidence`), and then kept on later updates.
- **Interruptions.** A resource becomes available when all its files are complete;
  an interrupted `init` keeps what it finished, and the next run continues.

Without `api_url`, `init` uses `OMNIPATH_API_URL` if set, otherwise the development
deployment `https://omnipath-metabo-dev.schaul.click/api`. Downloads are checked for
size and Parquet framing. `connect()` checks that every installed file is present
with its recorded size, and fails with a message naming `init()` when the directory
is empty or damaged.

## The database

`init` ends by writing `omnipath.duckdb` into the data directory and prints its path.

| In the database | What it is |
|---|---|
| `entity`, `relation`, ... | merged tables: one row per entity and per relation, with the `resources` that have it |
| `resources` | the installed resources, their versions and whether they have raw evidence |

These are **native DuckDB tables** built from the downloaded Parquet files: the same
database the hosted API serves a release from (`omnipath_client.database`), written
by the same code. There are no per-resource tables or schemas; the per-resource
Parquet files stay as the downloads. A full release (173 resources) takes about 22 minutes at 3 threads and 5 GB, and 18 GB.
Any change to the installed resources writes the database anew; the Parquet files
stay. Entities and relations are numbered when the database is built, `entity_n` and
`relation_n`, in key order; the child tables refer to them. Indexes cover `entity`
(`entity_n`, `entity_key`), `relation` (`relation_n`) and the evidence and annotation
tables' numbers.

To read rows by number quickly, write the numbers into the SQL as a literal list, in
a step of its own: DuckDB uses an index only then. With other conditions in the same
`WHERE`, or numbers from a subquery or join, it scans every row.

The database is written in DuckDB 1.5's storage format, so DuckDB 1.5 and newer
clients, the Python and R packages among them, can share a data directory.

The database is opened read-only, so several sessions and tools can use it at once.
Create your own tables as temporary tables, or attach the database to another one:

```sql
ATTACH '/path/to/omnipath.duckdb' AS omnipath (READ_ONLY);
```


## SQL

```python
with op.connect() as db:
    result = db.sql('''
        SELECT e.label, r.predicate, o.label AS object, r.resources
        FROM entity AS e
        JOIN relation AS r ON r.subject_n = e.entity_n
        JOIN entity AS o ON o.entity_n = r.object_n
        WHERE e.identifier = ?
        LIMIT 20
    ''', ["P04637"])
    print(result.fetchall())
    # result.df()                 # requires the pandas extra
    # result.to_arrow_table()
```

`sql` returns a native `DuckDBPyRelation`; rows are evaluated when you fetch or export
them. `db.connection` is the DuckDB connection itself. It executes local SQL; it is
not a sandbox. Outside Python, open the file with any DuckDB client:

```bash
duckdb -readonly "$(python -c 'import omnipath_client as op; print(op.data_directory())')/omnipath.duckdb"
```

The database has these tables:

| Table | One row per | Refers to |
|---|---|---|
| `entity` | entity, with its label, type, taxon, references and `resources` | |
| `entity_identifier` | identifier of an entity, merged by namespace and id | `entity_n` |
| `entity_annotation` | annotation of an entity (quantities as `quantity_*` columns), merged | `entity_n` |
| `entity_evidence` | source record stating an entity, with its `resource` | `entity_n` |
| `relation` | relation, with its endpoints (`subject_n`, `object_n`), sign, sources, `resources` and qualifiers | |
| `relation_annotation` | annotation of a relation, merged | `relation_n` |
| `relation_evidence` | source record stating a relation, with its `resource` and molecular forms | `relation_n` |
| `evidence_payloads` | raw source record (only where installed) | `relation_key`, `entity_key` |
| `relation_endpoint` | relation under each of its ends, sorted by entity and page order | `entity_n`, `relation_n` |
| `relation_reference` | the same under the ends' reference keys | `relation_n` |
| `entity_group` | gene, reference or connectivity group of an entity | `entity_n` |
| `entity_term` | lowercase label or identifier of an entity, sorted | `entity_n` |
| `ontology_edge`, `term_entity` | ontology edges between terms, and the entities and relations annotated with a term | |
| `relation_facet` | each distinct combination of the relation filter columns, with its relation count | |

`entity_n` and `relation_n` number the rows of `entity` and `relation`; the
64-character `entity_key` and `relation_key` identify the same entity or relation across
databases and releases. An entity or relation that several resources state is **one row**;
filter by resource with `list_contains(resources, 'signor')`, or use the evidence tables'
`resource` column (`SELECT * FROM relation_evidence WHERE resource = 'signor'`; their
`entity_id` / `relation_id` is the resource's own id). Relation counts count merged
relations. The lookup tables are derived by `init`; use them to find entities by name
(`entity_term`) and relations by endpoint (`relation_endpoint`) without scanning, and
which predicates, types, taxa or qualifiers occur together (`relation_facet`).

## Prebuilt queries

Every question the website asks is a query of the package (`omnipath_client/sql/queries/`),
the same SQL the API and the R client run, so the answers are the API's:

```python
with op.connect() as db:
    db.search("TP53")                                  # ranked entities, merged across resources
    db.search("gluc", entity_types="chemical_entity", limit=50)
    tp53 = db.search("TP53").fetchone()[0]             # entity_key

    db.relations("TP53", sources="signor")             # relations with annotations and evidence
    db.relations("gene:entrez:7157", predicates="affects", taxonomy_ids="9606")
    db.relations(protein_entity_keys=[key], isoform_identifiers=["uniprot:P04637-2"])

    db.entity(tp53)                                    # one entity's details, a dict

    # Any query by name, with its parameters: a DuckDB relation, or the API's JSON text.
    db.query("relation_page", entities=[tp53], limit=20)
    db.query("ontology_tree", term_id="GO:0006915", output="json")
```

- `search` and `relations` return DuckDB relations like `sql`, so they convert,
  filter, join and export the same way. `search` returns the entity table's columns:
  `entity_key`, `label`, `entity_type`, `namespace`, `identifier`, `taxon`,
  `relation_count`, `reference_entity_key`, `gene_reference_keys`.
- `relations` takes exact names, identifiers, entity keys or group keys such as
  `gene:entrez:7157`, at either endpoint, and the explorer's filters by name; an
  unknown filter name is an error that lists the valid ones. One value stands for a
  list of one. It returns what the explorer exports: one row per relation with its
  annotations, evidence and the product records its evidence names, at most 100,000
  relations. Larger selections are SQL on `relation`.
- `entity` takes the exact `entity_key`, as `search` returns it, and returns the
  explorer's detail record (the API's JSON as a dict).
- `query` runs any of the queries (`omnipath_client.steps.names()`: entity pages,
  relation pages and facets, grouped search, ontology search and trees, molecular
  context, statistics, ...). Each file's header lists its parameters, examples and
  outputs (some answer only as JSON). In
  API mode the API runs it (`POST /query/{name}`) and the table comes back as Parquet.
- Unknown filters and Biolink names (`predicates`, `entity_types`, ...) are errors
  naming them; a source that is not installed is an error listing the installed ones.

With every resource of a release installed, these give the hosted API's answers;
`init` also downloads the release's taxon names (`references/taxonomy/`), so taxa are
named as in the API. Resources the hosted API excludes from live queries are included
here. With a subset installed, counts and merged entities cover only the installed
resources. Product records in `relations` carry the entity columns of every installed
resource.

### The queries in SQL

A query is a sequence of SQL statements (steps) whose results later steps use as
literals; `packages/query/README.md` describes the format. `connect()` prepares each
connection for them: the contract's settings (`index_scan_max_count`,
`dynamic_or_filter_threshold`, ...) and the shaping macros of `shape.sql` as temporary
macros. Another DuckDB client can run them with `omnipath_client.steps.run(con, name,
params)` on its own connection, or read a step file and run its steps in order.

## Validation

From the repository root:

```bash
uv run pytest packages/query/python/tests -q
uv run python packages/query/python/examples/smoke.py
OP_THREADS=2 nice -n 19 uv run python packages/query/scripts/compare_with_api.py <data_dir>
```

The tests use Parquet files written by the core fixtures, a mocked catalog and file
host and real DuckDB queries; they compare the clients' answers with the API's on the
same files, in local and API mode, and check that the shared contract and its copies
are up to date. After changing the published or serving tables, Biolink or a query or
macro in `packages/query/sql/`, run
`uv run python packages/query/scripts/generate_contract.py`, which copies `sql/` into
this package and the R package. `compare_with_api.py` compares the three queries with
the API's on a real installation; the smoke script contacts
the default live server, installs SIGNOR into a temporary directory and queries it.
No deployment is required for client changes. Client instances are intended for one
thread.
