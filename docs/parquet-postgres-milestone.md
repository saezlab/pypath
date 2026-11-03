# Parquet-to-PostgreSQL milestone — 2026-09-30

## Implemented

The new omnipath_postgres workspace package loads an explicit release manifest
into a fresh PostgreSQL schema. It consumes resolved Parquet through DuckDB;
it does not import the builder, pypath or the native resolver.

Before loading, it checks exact resource pins, supported schemas, footer row
counts, file sizes and SHA-256 checksums. It rechecks the selected files before
commit. Schema creation, record loading, indexes, basic derived tables and
release metadata share one transaction. A failed load leaves no partial schema;
an existing destination cannot be overwritten.

The authoritative per-resource tables are entities, relations, identifiers,
evidence, annotations and payloads. They preserve original text keys, nested
records, duplicate occurrences, quantities, qualifiers and literal payload JSON.
Resource and release metadata retain original manifest text and checksums.

Release-wide entity and relation views combine identical published keys while
retaining every resource record. Relations aggregate evidence and sources.
Taxon is a consensus display projection, not an identity constraint. Conflicting
qualifiers under a shared relation key are rejected. Distinct relation keys are
retained even when endpoints and predicates are identical.

Basic identifier lookup, typed annotation quantities and direct relation counts
are available. Endpoint, identifier, predicate, entity-type, taxon, annotation,
evidence and payload indexes are created after bulk loading.

## Validation

- PostgreSQL 16.10 in disposable, Unix-socket-only local clusters.
- 106 PostgreSQL package tests plus 44 shared-core tests: 150 passed.
- Existing SIGNOR 0.1.1 fixture: 2 entities, 2 relations, 20 evidence occurrences
  and 20 payloads, with exact nested-record comparisons.
- Small synthetic fixtures exercise measurements at entity/relation/evidence
  scopes, standalone entity payloads, duplicate occurrences, shared identities,
  distinct qualified relations, conflicting taxa and self-loops.
- Failure cases prove rollback for missing endpoints, bad payloads, evidence-count
  mismatches, identity/qualifier conflicts and changed manifests.
- CLI loading, immutable destinations, manifest recovery, Ruff and formatting pass.
- Built core and PostgreSQL wheels; an independent installation reads the
  20-record sample and shared vocabulary without pypath, build or resolver installed.

No new resource or reference builds ran. Validation reused the existing bounded
sample and directly written tiny Parquet test fixtures. Temporary database
clusters were shut down and removed.

## Reproduce

Run make test-postgres. Set OMNIPATH_POSTGRES_BIN if the PostgreSQL binaries are
outside PATH and the standard Homebrew location. Run make check for lint/format.

For an explicit database load, use omnipath-postgres RELEASE_JSON --data-root ROOT
--schema NEW_SCHEMA with OMNIPATH_DATABASE_URL set. The package README documents
the release manifest and CLI.

## Remaining scope

This is the loading foundation, not a replacement for every legacy PostgreSQL
derived table. Next, port ontology-dependent projections and the remaining tables
required by MetSigDB, network views and COSMOS, then compare representative
protein, chemical, ontology and reaction queries. Full-scale load and query
performance, web backend selection and release activation remain unvalidated.
