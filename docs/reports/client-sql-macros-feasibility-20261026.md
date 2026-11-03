# Client queries as DuckDB macros: feasibility (26 October 2026)

> This is the record of the feasibility spike. Its result is in the clients now: the
> macros are `packages/query/sql/` (also shipped by the R package
> `packages/query/r`), the lookups and global ids are in
> `packages/query/python/omnipath_client/database.py`, and
> `packages/query/scripts/compare_with_api.py` replaces the spike's comparison
> scripts. Paths below that name `packages/client_refactored/` refer to the removed
> spike folder; the measurements stand. Since the spike, the macros use DuckDB's
> `lambda x:` syntax (DuckDB 2.0 rejects `x ->`), global ids are the resource's index
> times 2^32 plus its id, and the clients pass macro arguments as literals (a relation
> from `sql()` with bound parameters was ten times slower).

Can the client's prebuilt queries live in the DuckDB database instead of in Python, so
that R, Python, the DuckDB shell and the API all run the same SQL? This folder
answers that for all three, `relations`, `search` and `entity`, with working macros
compared with the engine on real installations (26 October 2026, DuckDB 1.5.6).

**Result: feasible, on the condition that the database uses native tables.**

- **Answers.** `op_relations` gives the engine's export exactly: same columns and types,
  same rows, in every case on every installation (23 cases on the largest). `op_entity`
  gives the engine's detail record exactly, as JSON, in all 14 cases. `op_search` gives
  the client's `search` exactly except for tie orders where the engine's own order
  changes with the page size (below).
- **Speed with native tables** (`materialize=True`). On 72 resources the macros are as
  fast as the engine or faster. On 74 resources with BioGRID and ChEMBL (6.5M relations,
  54M relation annotations), relations are mostly about twice as fast, and up to seven
  times on large selections, entity details are about as fast, and search takes
  0.1–0.5 s against the engine's 0.02–0.5 s. Each connection needs one setting,
  `SET dynamic_or_filter_threshold = 20000`; without it, small results on large tables
  were 2–3 times slower than the engine.
- **Speed with views** (the default). Every query pays a fixed cost that grows with the
  number of installed resources: about 0.3 s with 8 resources and 1.1 s with 72. The full
  release has 173. A view query also opens hundreds of files: the DuckDB shell failed one
  with "Too many open files".
- **No Python needed to query.** The macros run from the DuckDB 1.3 shell on a database
  written by DuckDB 1.5.


## Files

| File | What it is |
|---|---|
| `sql/relations.sql` | `op_relations` (the explorer's relation export), `op_entity_keys` (names, identifiers, CURIEs and group keys to entity keys), helpers |
| `sql/search.sql` | `op_search` (the explorer's entity search, in the client's columns) |
| `sql/entity.sql` | `op_entity` (the explorer's entity details, as the API's JSON) |
| `lookups.py` | the per-installation lookup macros, for views and for native tables; `add_global_ids` for native tables |
| `apply.py` | writes `op_names` (Biolink names), the lookups and `sql/*.sql` into an installed database |
| `compare_relations.py`, `compare_search.py`, `compare_entity.py` | engine against macro: answers and times; `OP_THREADS=2` limits both |

```bash
uv run python packages/client_refactored/apply.py <data_dir>
uv run python packages/client_refactored/compare_relations.py <data_dir>
uv run python packages/client_refactored/compare_search.py <data_dir>
uv run python packages/client_refactored/compare_entity.py <data_dir>
duckdb -readonly <data_dir>/omnipath.duckdb -c "SELECT * FROM op_search('TP53')"
```

Pass `<data_dir>` as a normalised absolute path (see "Found on the way"). `apply.py`
changes the database in place. The installed client still works afterwards, but its next
rewrite of the database removes the macros.

Test installations, in `local-artifacts/client-refactored/` (git-ignored), from the dev
deployment's release on 2026-10-26:

| Directory | Resources | Mode | Disk |
|---|---|---|---|
| `subset` | 8: SIGNOR, STRING, Reactome, ChEBI, CollecTRI, TRRUST, HGNC, CellPhoneDB | views | 442 MB |
| `wide` | those 8 + the 64 resources under 1 MB | views | 460 MB |
| `wide-native` | the same 72 | native tables + global ids | 1.8 GB |
| `big` | those 72 + BioGRID and ChEMBL | native tables + global ids | 11.7 GB database + 2.6 GB files |

## How it is built

Two layers:

1. **Query macros** (`sql/`): the same for every installation and every language. They
   are written against the database's tables and against lookup macros named
   `op_<table>_at(pairs)` (rows by a resource's own id) and `op_<table>_with(values)`
   (rows by key or term).
2. **Lookup macros** (`lookups.py`): written by `init` along with the database, from the
   list of installed resources. This is the only generated SQL. An R `init` would
   produce it with the same loop it uses for the views.

Argument checks stay in SQL: an unknown predicate or entity type is an `error()` naming
it, using the table `op_names` (Biolink classes and predicates under every name the
engine accepts). Per-language wrappers would only turn one value into a list of one.

`op_entity` returns one row holding one JSON value, the record the API and
`Client.entity` return. It rebuilds the explorer's shaping in SQL: identifier
deduplication, the sources list, display names, the ontology block, relation
qualifiers and labels. The main entity and every endpoint entity go through the same
shaping steps, so that SQL exists once. Optional keys, such as an annotation's
`quantity`, are written per item as JSON, so they are left out where the engine leaves
them out.

### Why lookups

The engine reads only the resources that have something to find, through each
resource's file with a list of ids. A SQL query has to name the files it reads when it
is planned, and DuckDB does not push a join's ids through the views over every
resource. So one relation lookup read all 1.3M relation rows of the 8-resource install,
every column, even when nothing matched.

- **Over views**, a lookup is a UNION ALL of the resources' own views, each filtered by
  `IN (SELECT ...)`. DuckDB then passes a Bloom filter of the wanted ids into each file's
  scan and reads only the row groups holding them. STRING read 2,575 of 430,439 rows for
  TP53. But the plan has one branch per resource and lookup, about 10 ms each. Every file
  is also opened and its first chunk read, even when no id is wanted from it.
- **Over native tables**, a lookup is one filter on one table. Each resource numbers its
  relations and entities from 0, so `add_global_ids` adds `_rid` and `_eid` (the
  resource's offset plus its id) and sorts the tables by them. A lookup is then a
  one-column filter that DuckDB checks against each row group's range. TP53's 3,089
  relations read 6,482 rows of 1.4M, a lookup with nothing wanted reads nothing, and
  planning takes about 1 ms, whatever the number of resources.

On large tables two more things matter:

- **`SET dynamic_or_filter_threshold = 20000`** on every connection. A lookup's wanted ids
  reach the scan as an explicit list only when there are at most this many, 50 by
  default; above that, DuckDB passes only their range. Ids scattered over several
  resources span most of a big table, so glucose's 183 relations read 3.6M relation
  rows, 4.8M evidence rows and 30.8M annotation rows. With the setting, glucose went
  from 0.91 s to 0.27 s, ATP from 0.89 s to 0.24 s, and nothing got slower. The setting
  belongs to the connection, not the database: every client's `connect()` sets it, and
  so must a DuckDB shell user.
- **Search conditions as constants.** The prefix range, the CURIE parts and the phase
  are computed from the macro's `query` argument directly, not read from a CTE.
  Constants reach the scan, so the term range prunes by sort order: TP53 read 4,940 of
  16.2M terms instead of all of them. Search on the big install went from 1.1–1.5 s to
  0.1–0.5 s.

Smaller native row groups (20,480 rows instead of 122,880) made no difference at
realistic sizes: 0.55 s against 0.63 s for glucose's evidence.

Tried and dropped:

- One multi-file `read_parquet` per table: it reads 72 rows for one id. But ids from
  different resources select row groups in each other's files. Whether DuckDB passes the
  Bloom filter at all depends on its size estimates: the same lookup read 24,790 rows in
  one form and 1.4M in another.
- Raising `dynamic_or_filter_threshold` for views: no effect through the views. (For
  native tables it is needed; see below.)
- Joining on (file, id) pairs: precise for one pair, but DuckDB uses no Bloom filter for
  two-column joins.

### Two engine-independent fixes

- Lists in ordinal order (annotations, evidence, product records) are collected
  unordered, then sorted with `list_select(list(x), list_grade_up(list(ordinal)))`.
  `list(x ORDER BY ordinal)` over these nested structs is about ten times slower: 3.3 s
  against 0.26 s for SIGNOR's evidence. The engine's export uses the slow form, and this
  is most of why it takes 3 s for all of SIGNOR. The fix applies to the engine as is.
- A macro argument is copied wherever its parameter appears. A lookup therefore
  evaluates its argument once into a materialised CTE, and so do the query CTEs used
  more than once. Without that, one lookup took 1.2 s to plan on 72 resources.

## Results

Warm times in seconds, the median of three runs. "Engine" is
`ParquetServingEngine.export_slice` and `Client.search`, which read the Parquet files in
every mode.

### `op_relations`, 72 resources

| Case | Rows | Engine | Macro, native | Macro, views |
|---|---:|---:|---:|---:|
| TP53 | 3,911 | 1.00 | 0.28 | 1.72 |
| TP53 in SIGNOR | 238 | 0.14 | 0.12 | 0.91 |
| TP53 and MDM2, both ends | 29 | 0.23 | 0.21 | 1.17 |
| gene group of TP53 | 3,381 | 1.13 | 0.30 | 1.72 |
| ATP (CHEBI:30616) | 3,189 | 0.60 | 0.21 | 1.80 |
| glucose | 181 | 0.35 | 0.17 | 1.48 |
| EGFR protein–chemical | 77 | 0.20 | 0.20 | 1.41 |
| EGFR, `biolink:interacts_with` | 566 | 0.46 | 0.26 | 1.60 |
| all of SIGNOR | 35,544 | 3.32 | 0.86 | 1.04 |
| SIGNOR negative, human | 11,521 | 1.33 | 0.37 | 0.68 |
| CollecTRI directed | 43,536 | 3.03 | 0.80 | 1.01 |
| molecular form (protein key) | 1,071 | 0.68 | 0.33 | 1.49 |
| molecular form (isoform, target) | 64 | 0.24 | 0.21 | 1.40 |
| TRRUST, 2+ evidence | 30 | 0.07 | 0.05 | 0.53 |
| no match | 0 | 0.08 | 0.04 | 1.13 |

All 19 cases give the same answer in all three installations. On 8 resources with views,
the floor is 0.3 s, so the cost of views grows with the installation, not the data.

### `op_search`, 72 resources, 2 threads

| Case | Engine | Macro, native |
|---|---:|---:|
| TP53 | 0.08 | 0.08 |
| CHEBI:30616 | 0.06 | 0.10 |
| kinase | 0.08 | 0.07 |
| a (many hits) | 0.25 | 0.16 |
| insulin receptor | 0.08 | 0.08 |
| empty query | 0.13 | 0.12 |
| TP53 in SIGNOR | 0.02 | 0.05 |
| TP53 among MDM2, TP53 | 0.15 | 0.11 |

Views were measured before the search fixes above, at 0.3–0.7 s on all CPUs.

Every case gives the same rows in the same order, except `a` (and `ATP` on the big
install). For `a`, the same 20 entities come out with four equally matching ones in another order. The
engine reads prefix hits in growing batches and orders ties by the relation counts of
the hits read so far, so its own order changes with the page size: `search('a')` puts
"hydrogen acceptor" second at `limit=20` and "adenine" second at `limit=100`. The macro
ranks every hit at once, as the engine would with an unlimited page.

Bound parameters (`op_relations(entities := ?)`) are as fast as literals.

### BioGRID and ChEMBL added: 74 resources, native tables, 2 threads

6.5M relations (ChEMBL 2.7M, BioGRID 2.4M), 9.6M evidence rows, 54.5M relation
annotations, 2.1M entities, 16.2M terms. Engine and macro both limited to 2 threads at
low priority (`OP_THREADS=2`), so the times are higher than on all CPUs.

| `op_relations` | Rows | Engine | Macro |
|---|---:|---:|---:|
| TP53 | 14,488 | 12.36 | 6.38 |
| TP53 in BioGRID | 3,241 | 4.07 | 1.54 |
| TP53 and MDM2, both ends | 69 | 1.34 | 1.39 |
| EGFR | 16,738 | 13.34 | 6.52 |
| EGFR in ChEMBL, human | 10,889 | 7.00 | 3.95 |
| EGFR protein–chemical | 11,047 | 7.52 | 4.70 |
| ATP (CHEBI:30616) | 3,189 | 1.08 | 0.39 |
| glucose | 183 | 0.83 | 0.42 |
| imatinib | 232 | 1.76 | 3.22 |
| all of SIGNOR | 35,544 | 4.93 | 0.76 |
| CollecTRI directed | 43,536 | 5.34 | 0.69 |
| molecular form (protein key) | 9,205 | 8.76 | 7.88 |
| molecular form (isoform, target) | 2,184 | 5.44 | 7.25 |
| no match | 0 | 0.10 | 0.05 |

All 23 cases match. Two are slower than the engine. For imatinib (232 relations), DuckDB
passes no id list into the evidence and annotation lookups, despite the setting above:
they read all 9.7M evidence rows (2.5 s of the 2 s wall time is that scan's CPU) and
12.3M annotation rows. The isoform filter, with no entity given, scans the molecular
forms of all 9.6M evidence rows; the engine does the same scan, faster.

| `op_search` | Engine | Macro |
|---|---:|---:|
| TP53 | 0.11 | 0.20 |
| CHEBI:30616 | 0.08 | 0.23 |
| imatinib | 0.07 | 0.12 |
| a (many hits) | 0.19 | 0.40 |
| empty query | 0.49 | 0.51 |
| xyzzy-nothing (substring phase) | 0.10 | 0.07 |
| TP53 in SIGNOR | 0.02 | 0.14 |

Two lookups in search still get no id list from DuckDB (the hits' entity rows and the
page's other copies) and read 1.6M and 2.0M entity rows; that is most of the remaining
0.1–0.2 s.

| `op_entity` | Relationships | Associations | Engine, first call | Macro |
|---|---:|---:|---:|---:|
| TP53 | 48 | 100 | 0.48 | 0.27 |
| EGFR | 50 | 100 | 0.40 | 0.30 |
| imatinib | 1 | 82 | 0.18 | 0.09 |
| CHEBI:30616 | 50 | 100 | 0.17 | 0.10 |
| in most resources | 50 | 100 | 0.10 | 0.30 |
| unknown key | 0 | 0 | 0.02 | 0.04 |

All 14 cases match. The engine caches details, so its first call is the fair
comparison (its later calls take 0.01–0.09 s).

Install cost: writing the 74 resources' native tables took 74 s. Adding the global ids
afterwards took 65 s and 10 GB of memory (no memory limit set), and doubled the file to
11.7 GB: DuckDB does not shrink a file when a table is replaced. `init` would write the
ids with the tables in the first place.

## Found on the way

- **Client bug:** `data_directory()` makes a path absolute but does not normalise it. A
  path with `..` therefore no longer matches the one stored in the database, and
  `connect()` rewrites the whole database. Fix: `.resolve()`.
- **Opening with older DuckDB:** the DuckDB 1.3 shell opened databases written by DuckDB
  1.5, both views and native. This was an open question in the client notes.
- With 72 resources a view query opens several hundred Parquet files. The DuckDB shell
  hit the system's open-file limit; Python, with a higher limit, did not.

## What this would mean

- **Native tables become the query mode.** The client README measured 13.6 GB extra and
  under 3 minutes at `init` for an older 46-resource release of 5.4 GB. The current
  release, 9.4 GB of files, needs more. 72 resources took 18 s here. The views mode would
  stay for SQL only, or `init` would materialise by default.
- The client's in-place update of native tables (`database._update`) would also have to
  write the global ids and lookups; here `apply.py` rebuilds them.
- Every connection sets `dynamic_or_filter_threshold`: the R and Python `connect()`, and
  the documentation for DuckDB shell users.
- `_rid` and `_eid` are new columns in the native tables. The macros leave them out of
  their answers; plain SQL users see them beside `_resource`.
- For the API to use the same SQL, the API's own database needs the same native tables
  and macros. Otherwise there are two implementations again.

## Not done

- The full release (173 resources, 9.4 GB of files): measured on up to 74 resources with
  the two largest relation tables. The other big resources (MSigDB, FooDB, SwissLipids,
  BindingDB, 2–4M relations each) would roughly double the table sizes again.
- The two slow relation cases (imatinib, the isoform filter) and the two search
  lookups without an id list.
- An R client calling the macros (R is not installed here). The shell call shows the SQL
  needs nothing from Python.
- `scope_entity_ids` (the explorer's selection scope) in `op_relations`, and pagination
  cursors in `op_search` (the client uses neither).
- Biolink normalisation of qualifier values (`annotation_value`): values pass as stored.
