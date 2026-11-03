# Source attributes and PostgreSQL validation

Completed 2026-09-30 on nicesrv. The code and bounded validation cover every
former PostgreSQL product payload reader identified in the
[dependency audit](payload-dependency-audit.md).

## Result

Seven affected input resources now publish the attributes their downstream
products need. PostgreSQL retains evidence, annotations, quantities and source
hashes, while complete original source records remain exclusively in Parquet.
No new generic source-context table or identifier-resolution step was added.

The private server replay consumed **20 original source records per resource,
140 total**, through the actual inputs_v2 mappers, native resolver, writer and
immutable publication. It used the existing prototype reference at
`data/reference/.library-compact-20260912`; it made no source downloads.

| Resource | Original records | Entities | Relations | Raw Parquet rows retained | Reaction evidence with SHA/type |
| --- | ---: | ---: | ---: | ---: | ---: |
| rhea | 20 | 239 | 263 | 310 | 293 |
| recon3d | 20 | 26 | 20 | 30 | 30 |
| metatlas | 20 | 52 | 68 | 113 | 83 |
| kegg | 20 | 1026 | 1197 | 1200 | 64 |
| reactome | 20 | 101 | 396 | 539 | 95 |
| macdb | 20 | 16 | 10 | 20 | 0 |
| connectomedb2025 | 20 | 22 | 20 | 20 | 0 |

The first three resources cover reaction, metabolic and transport datasets.
KEGG includes five pathway records alongside fifteen reactions; MACdb includes
six traits and fourteen linked associations. Multiple membership edges and
merged evidence explain why output counts exceed the original-record cap.

## PostgreSQL checks

Build commit: `7ed2e63`; input commit: `e27d72e09`; PostgreSQL/subset commit:
`4c896eb`. Both migration branches are pushed to GitHub and pulled into the
isolated server checkout.

Private resource/release version: `2026.9.30.3`. Destination:
`omnipath-migration-postgres`, localhost port 5440, database
`omnipath_migration`, new schema `annotation_sample_20260930`.
The existing `pilot_20260930` schema was preserved.

Every projected row matched the published Parquet values: **1,482 entities,
19,781 identifiers, 1,974 relations, 2,226 evidence occurrences and 5,772
annotations**. Import verified and discarded 2,232 raw rows. Neither a payload
table nor a full raw-text column exists in the new schema. Matching reaction
source rows are checked against exact text SHA/type annotations during import;
invalid, stale or conflicting provenance fails atomically.

The Parquet directory was temporarily renamed after loading. With it unavailable:

- Reaction rebuild reproduced the exact same 95 contexts and 565 participants.
  There were no conflict diagnostics; Rhea's unreported direction stayed unknown.
- COSMOS built 1,690 edges over 75 selected source events.
- MetSigDB built 58 Reactome, 77 KEGG and 10 MACdb memberships. MACdb's subtype
  values were retained (`gene abnormality` and `phenotype` in this sample).
- LIANA returned all 20 selected interactions with explicit endpoint roles.
  Reaction views returned a populated page. MetaLinksDB queried successfully but
  had no eligible rows in these samples; populated binary and transport cases
  are covered by the separate integration fixtures.

The directory was restored and its manifests/checksums revalidated. Twenty-one
original prototype Parquet files (656,315,126 bytes) still matched their original
manifests, and the original manifests matched the pre-run selection hashes.
Production PostgreSQL and all prototype containers retained their exact IDs and
start times. Monthly PostgreSQL versus independent Parquet/API scheduling was
unchanged.

## Tests and timings

- Core/build/API regression suite: 394 passed, 1 skipped. The added direct
  reaction provenance case and nonvacuous assertions also passed in the six-case
  focused suite.
- PostgreSQL: 154 passed; the final index-predicate adjustment passed all 39
  reaction tests.
- Subsets: 31 passed against the final loader, including populated products
  after their input Parquets become unavailable.
- New input-context tests: 40 passed, including alignment, conflicting direction,
  symbolic coefficients, canonical endpoint flips and separate source events.
- Ruff and whitespace checks passed for the changed migration code. The input
  review confirmed five existing KEGG E402/F811 findings were unchanged.

Native builds took 3.0–5.7 seconds per source. Database load/commit took
2.823 seconds, followed by a 0.945-second
PG-only product rebuild. Exact phase timings are in the [build log](build-log.md).
These small-sample timings do not predict full-release performance.

## Scope and next release

The capped versions are private validation artifacts and have not replaced the
full published resource versions. WikiPathways and ClassyFire were outside the
seven-resource server sample; their populated paths are tested in fixtures.
These results establish attribute survival and independence from raw records,
not full-resource biological parity or a full 46-resource PostgreSQL import.

For a complete release, rebuild new full versions of the same seven resources,
retaining all previously published datasets in each resource version. Reuse the
other pinned, unaffected resource versions. Legacy reaction artifacts without
source SHA attributes intentionally fail with rebuild guidance.

Private evidence: `/root/projects/omnipath-migration/annotation-migration-20260930/`:
`selection.json`, `data/inspection.json`, `postgres-check.json`,
`original-artifacts-recheck.json`, private scripts and logs. No raw records or
credentials are committed to this repository.
