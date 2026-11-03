# OmniPath PostgreSQL

Load an explicit release of already resolved Parquet resources into the current main branch's physical relational model. Entity resolution, resource parsing, Parquet schemas and independent resource updates stay upstream. PostgreSQL keeps its separately selected release and monthly publication schedule.

```json
{"schema_version": 1, "version": "2026.09", "resources": {"signor": "0.1.1"}}
```

```sh
export OMNIPATH_DATABASE_URL='postgresql:///omnipath'
uv run --frozen omnipath-postgres release.json \
  --data-root data --schema release_2026_09
```

The data root contains `resources/<source>/<version>/` and may be a local directory or an HTTPS base URL. The release manifest may independently be a local file or an HTTPS URL:

```sh
uv run --frozen omnipath-postgres https://data.example.org/releases/2026.09.json \
  --data-root https://data.example.org --schema release_2026_09
```

Remote artifact access requires HTTPS and exact HTTP byte ranges (HTTP loopback is supported for development fixtures). Schemas and footer row counts are read by range; SHA256 checks stream complete artifacts both before loading and before the base checkpoint, as for local files. Entity/relation data is then staged once into the local DuckDB working database and shape validation uses those staged tables. Subsequent projection queries make no remote reads. There is no persistent downloaded-Parquet cache; allow local working disk for DuckDB/CSV staging and choose `--temp-directory` when needed. The endpoint must retain the pinned versions for the duration of a build. PostgreSQL-only finish does not require the endpoint.

 Loading requires a new destination schema. The CLI and package-level `load_release` use the main layout: UUID canonical `entity` and `relation` tables, numeric vocabulary/source/dataset dictionaries, global identifier and annotation dictionaries, source-partitioned evidence and links, and main's constraints and indexes. Serving schema 5 loads (the normalized resource tables) add required molecular reference/occurrence context, `annotation_quantity` metadata and small release/resource/phase records. Scientific annotations, qualifiers, measurements, identifiers and source-scoped evidence remain in the normalized tables; quantity metadata preserves units, prefixes, comparators and source fields needed by their readers. Exact published keys, original order and repetition, raw provenance fields and presentation hints remain in the immutable pinned Parquets. Full source bodies stay out of PostgreSQL. Molecular evidence context retains its exact nested occurrence JSON.

For PostgreSQL inspection copies, pass `retain_published_provenance=True` to `load_release` or add `--retain-published-provenance` to a new CLI load. This also creates and loads `parquet_entity`, `parquet_statement`, `parquet_evidence`, `parquet_identifier_occurrence` and `parquet_annotation_occurrence` with their exact published audit rows. These tables have no downstream scientific reader and can be large because they retain every repeated occurrence. Default loads omit their DDL, staging and COPY. Both modes produce the same normalized scientific tables and mandatory context tables; returned loaded counts describe only tables actually copied, while compatibility metadata separately records published input counts and the provenance location.

`entity_reference_context` retains each published product/entity key, its exact reference key and ordered gene reference keys; virtual gene anchors need no corresponding entity row. `statement_reference_context` retains the two endpoint reference keys. `molecular_evidence_context` retains each relation evidence pair in one occurrence row, plus standalone entity evidence, including product keys, isoforms, modifications, variants and source coordinate references. Resource/version/owner/ordinal identify occurrences without deduplicating repeated evidence; relation occurrences also link to their existing main relation-evidence UUID. These tables are loaded even when audit copies are disabled and cleared by content resets. The loader creates no states, gene-to-protein representatives, or product fanout.

DuckDB prepares relational tables once in an on-disk database and writes CSV chunks. Python forwards their bytes through PostgreSQL COPY without converting biological rows into Python objects. Source artifacts are read but never rewritten. Main keys, foreign keys and deferred indexes are restored and validated before the base checkpoint is advertised. Artifact schemas, sizes, row counts and SHA256 checksums and COPY counts remain build checks.

Main's downstream implementations provide search and ontology tables, chemical structure groups, classifications, interactions and resource facts, labels, bitmap filters, overlap summaries, resource metadata, MetSigDB, network presets and COSMOS. Biolink readers adapt predicates, scopes and quantities; the original main table and index contracts remain the reference. The [column map](../../docs/postgres-parquet-column-map.md) and [parity checklist](../../docs/postgres-main-parity-checklist.md) document unavailable original resolver/occurrence facts and deliberate adaptations. Published identities are explicitly marked `published`; they are never falsely marked as originally matched. MetSigDB uses the previously agreed published-chemical policy.

Main's legacy `ontology_terms` base table stays empty on fresh builds, as in the current main projector. The existing `entity_ontology_term` derivation builds serving terms from ontology edges, relational identifiers and annotations. Published ontology statements use the actual `ontology` kind and known input-provider ontology families; resolving a chemical to InChIKey does not change its source ontology.

The verified base commits before downstream work. Each complete product commits together with its phase record; a failed product rolls back its unfinished work. Completed products are retained on a retry:

```sh
# Stop after base tables, constraints and indexes.
uv run --frozen omnipath-postgres release.json --data-root data \
  --schema release_2026_09 --base-only

# Finish derivations and remaining products from PostgreSQL; no repeated COPY.
uv run --frozen omnipath-postgres --finish --schema release_2026_09

# Resume only selected products; completed ones are skipped.
uv run --frozen omnipath-postgres --finish --schema release_2026_09 --products cosmos
```

The Python API accepts `base_only=True` and `products=('cosmos',)`. `finish_release(database_url, schema=...)` uses stored pins and durable phase records. It preserves existing loaded tables and stored counts, including older builds with PostgreSQL audit copies; the new retention option applies only to fresh loads and does not clean up existing schemas. A schema session advisory lock spans commits and prevents overlapping builds. `--checkpoint-base` and `--defer-constraints` remain accepted for CLI compatibility; both behaviors are already enabled in the main layout. An interrupted unvalidated COPY is not a resumable base.

DuckDB defaults to one thread and 1 GB memory. Use `--duckdb-threads`, `--memory-limit` and `--temp-directory` to choose working limits and staging storage. Temporary CSVs and DuckDB files are removed after loading or failure. Returned timings separate artifact checks, projection, CSV staging, COPY, constraints and durable phases.

Use PostgreSQL 18 with main's [PostgreSQL image](../../postgres/Dockerfile) or an equivalent image providing `roaringbitmap` and `pg_trgm`. The development parity tests compare actual catalog definitions and duplicate-sensitive scientific fixtures with the frozen main source. They are bounded development checks; loading has no mandatory exhaustive verification, rebuild-and-rollback run or automatic safety phase.

The earlier resource-record loader (`omnipath_postgres.compatibility.record_layout`) has been removed together with its tests; schemas built by it are not readable by this package. See [the alignment plan](../../docs/postgres-main-parity-plan.md) for the migration boundary.

The canonical loader is `omnipath_postgres.loader`; shared schema and derivations
live in `omnipath_postgres.relational`. MetSigDB, network presets and COSMOS are
owned by [subsets](../subsets/README.md). The dependency runs from PostgreSQL to
subsets to core, and both commands use the same product runner.

## Measuring HTTPS input

Validation reports `artifact_validation_complete` and
`artifact_revalidation_complete` observer events with actual HTTP request counts,
response bytes and elapsed seconds. DuckDB's own projection transfers are separate;
these counters do not claim to measure the full network total. Revalidation time
is separately recorded in `phase_seconds`. Hash and corruption checks are unchanged.

For a bounded single-resource release, use:

```sh
uv run --frozen --package omnipath-postgres python scripts/measure_parquet_input.py \
  single-resource-release.json --data-root https://data.example.org --max-bytes 67108864
```

The command validates, projects in temporary DuckDB storage, and revalidates.
It opens no PostgreSQL connection and publishes nothing. It rejects multi-resource
releases and declared artifacts above the byte budget. Increase the budget explicitly
for a representative larger resource; retain both validation and projection times
when comparing local and HTTPS input. No persistent download cache is created.
