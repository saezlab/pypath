# Molecular forms on parquet-migration

5 October 2026

The gene-reference and molecular-form changes are integrated on the existing
`saezlab/omnipath-build:parquet-migration` branch. Pypath uses its existing
`saezlab/pypath:parquet-migration` branch at
`9d61c6e6cf1c9c0355d1287758ade6df2d477455`; the submodule pins that commit.

## Integration

The implementation came from prototype commit
`b3451f3a16ddb12c581ae28c00a3f928faa4057f`, taking only its changes since
`86f978e39bd1f7a40180714acf03747f3426e87c`. The two repositories have independent
Git histories. That feature diff was adapted onto migration commit
`5c55c86b5c140854ce6b433cc703d3024fa0cf0f`, preserving the migration branch's
existing work rather than replacing it with the prototype tree.

The flat Python package layout, resolver/runtime ownership, bounded API query
pool, PostgreSQL loader and subset runner remain. Build compiles reference
artifacts; resolver consumes them without importing build. The new Python client
is a separate package and queries published Parquets without the server.

The three serving Parquets remain. Build manifests retain `schema_version = 1`
and filename-keyed file metadata; their serving contract becomes version `4`.
Optional `.serving/v2` projections belong to individual resource versions.
Cross-resource grouping happens at query time; no cross-resource web build was
introduced. The Group toggle applies gene-reference and chemical-connectivity
grouping together.

## PostgreSQL and subsets

PostgreSQL needed three required context tables:

- `entity_reference_context`: entity/reference keys and supported gene links.
- `statement_reference_context`: both relation endpoint reference keys.
- `molecular_evidence_context`: complete ordered evidence occurrences, keeping
  subject/object forms paired and preserving standalone evidence, exact features,
  specific identifiers, source coordinates, nulls and empty lists.

These are loaded even when optional audit copies are disabled, link to existing
identities/evidence, and participate in content reset. Main graph identities stay
as published. The loader does not resolve again, infer protein participants or
create entities for PTM/variant combinations. The historical record-layout
compatibility reader also retains the new fields in its records.

COSMOS now retains canonical GeneID catalyst labels, including protein-typed gene
references, instead of choosing an arbitrary UniProt alias. Genuine UniProt
products retain their labels. The optional utils path also preserves GeneID and
accepts other catalyst mappings only when one protein target answers. Connector
labels follow the same rule. MetSigDB's chemical memberships and the network
preset recipes need no semantic change.

## Validation

Checks ran in the isolated migration checkout with the pinned pypath commit and
rebuilt native resolver/reference binaries. Broad runs exposed two stale test
expectations after the port; their corrected focused checks passed. No broad
suite is presented as a single all-green rerun after those fixture corrections.

| Check | Result |
| --- | --- |
| Build suite | 382 passed, 2 skipped; the remaining relocated checksum-helper test passed on its focused rerun |
| Core/resolver/API/client suite | 382 passed; the remaining serving-version expectation was corrected and both migration-serving tests passed |
| Final writer/runtime-boundary/publication checks | 30 passed |
| PostgreSQL suite | 375 passed, 204 skipped; opt-in database/external cases are not counted as passed |
| Focused actual PostgreSQL COPY/reset/load/finish | 16 passed against disposable PostgreSQL 16 |
| Full subset suite | 43 passed, including disposable PostgreSQL fixtures |
| Rust reference policy tests | 29 passed with Parquet input enabled |
| Native resolver extension | Built and Rust test/doctest commands passed |
| Web | 32 tests passed; Svelte check and production build passed |
| Isolated wheels | Core, resolver, build, API and client passed installation/import boundary checks |
| Repository contracts | Ruff lint/format, OpenAPI and Biolink generation checks, lockfile check and diff whitespace check passed |

The [earlier full SIGNOR/IntAct pilot](molecular-form-implementation-20261004.md)
and its JSON reports retain their original producer revisions and measurements.
They are not fresh full resource builds from this integration checkout. This
integration performed no production deployment or full monthly PostgreSQL rebuild.
A new selected PostgreSQL load is needed to populate context tables in a release;
existing snapshots and committed product checkpoints are not altered by this code
push. Remaining resource coverage and transcript mapping gaps are recorded in
the [input audit](molecular-form-input-coverage-20261004.md).
