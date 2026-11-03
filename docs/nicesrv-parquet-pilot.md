# nicesrv resolved-Parquet PostgreSQL pilot — 2026-09-30

The SIGNOR, ChemOnt and Rhea pilot loaded and passed independent verification.
This establishes the migration path on existing full resource artifacts; a full
release import and complete product parity remain separate milestones.

## Inputs and isolation

- All three resource versions: `2026.9.5.17`.
- Private snapshot: `2026.9.30.1`; canonical manifest SHA-256
  `91a853954cf1a5056ebf3bdd219814f258e479f0c22fc77ee97be74a15e1e146`.
- Existing source files: `/root/projects/full_parquet/data`; nine files totaling
  140,970,841 bytes (134.44 MiB).
  Schema, footer counts and all content hashes passed validation before loading
  and were rechecked before commit. Files were read in place; no rebuild,
  resolution, cache access, copying, or release publication ran.
- Loaded code: `fa439b380c692de7531aa9dc067c9a3843bacb3b`. The pypath branch
  remains `fc6c262a6e1127720238b0ebde54f5d5951044af`; PostgreSQL does not import it.
- New container `omnipath-migration-postgres`, new volume
  `omnipath-migration-postgres-data`, network `omnipath-migration-pilot`.
- PostgreSQL listens only on `127.0.0.1:5440`. Database and role
  `omnipath_migration`; schema `pilot_20260930`.
- Cached PostgreSQL 16 image pinned by ID `sha256:52f87503b2fe5ad2e0d209b1a55830db42eb3ecd3eaf1b583accb0091ab30bff`.
  Limits: 2 CPUs, 4 GiB RAM; shared_buffers 256 MB, work_mem 16 MB,
  maintenance_work_mem 256 MB, max_parallel_workers_per_gather 1.
- Separate venv uses dependencies exported from the frozen migration lockfile.
  Only core, PostgreSQL and subset workspace packages are installed. Credentials remain
  in the private pilot directory, outside Git; this report contains no password.

The existing production PostgreSQL was never queried or changed. Its container,
volume and port remain separate. The old PostgreSQL and all three prototype
service IDs, startup times and running states match preflight observations.
Independent Parquet updates and the monthly PostgreSQL schedule are unchanged.

## Loaded rows and runtime

| Table | Rows |
| --- | ---: |
| entities | 276,326 |
| identifiers | 2,507,201 |
| relations | 444,863 |
| evidence | 890,464 |
| annotations | 657,097 |
| payloads | 892,545 |

Counts describe source-scoped rows; shared entities and statements are exposed
through canonical views without losing their original resource records.

| Loader phase | Seconds |
| --- | ---: |
| validate_files | 0.13 |
| copy_and_validate | 313.03 |
| indexes_and_derived | 384.53 |
| recheck_files | 0.13 |

Total successful load including commit: **697.88 s**.
Per-resource projection/copy times include interleaved flat-table inserts:

- chemont: 1.33 s.
- rhea: 174.18 s.
- signor: 31.64 s.

Independent verification: **176.28 s**. It reread every
published entity and relation using Arrow and compared complete nested JSON
and all scalar columns with PostgreSQL. Every raw payload was compared by its
original ordinal and exact text. Base and nested-array row counts were checked
independently, and recorded resource manifest hashes and Parquet hashes matched.

## Derivations and checks

- 22,365 source-scoped ontology paths; direct anchors, nonreflexive
  paths and transitive completeness passed, with source scope and hierarchy kind
  retained. Release-wide paths include every scoped path at an equal or shorter depth.
- 36,686 source-event reaction contexts. Participants match eligible
  original evidence occurrences exactly, including source, dataset, row and role.
- 39 indexes are live, ready and valid. All constraints are validated.
Independent checks returned zero violations:

- flattened_identifiers: zero violations.
- flattened_evidence: zero violations.
- flattened_annotations: zero violations.
- ontology_direct_edges: zero violations.
- ontology_valid_paths: zero violations.
- ontology_release_paths: zero violations.
- ontology_transitive_paths: zero violations.
- reaction_evidence_occurrences: zero violations.
- reaction_payload_text: zero violations.
- reaction_coefficients: zero violations.
- indexes_valid: zero violations.
- constraints_valid: zero violations.

| Reaction role | Context status | Participants | Known compartments | Coefficient annotations |
| --- | --- | ---: | ---: | ---: |
| enzyme | unavailable | 658,434 | 0 | 0 |
| product | verified_member_ordinal | 90,298 | 542 | 0 |
| reactant | verified_member_ordinal | 82,510 | 4,986 | 0 |

All 36,686 Rhea contexts have unknown direction and no structured coefficient
annotations. These three resources have no typed-quantity annotations, so real
numeric-measurement coverage needs another resource. Raw payloads, including
original equations, remain intact; compartment coverage is shown above. This pilot does not validate a
complete MetSigDB/COSMOS build or decide the web backend.

## Storage and next gate

Committed database size: **7.70 GiB**. Compressed Parquet
size is not a database-size estimate: PostgreSQL stores flat indexes and nested
records plus repeated original evidence payloads. Estimate full-release storage
and choose resource/memory limits before expansion.

Next: validate exact pins for the full resource selection, import into a fresh
schema on this isolated instance, then check complete downstream products and
representative queries. Do not replace the production PostgreSQL or publish the
private pilot manifest to the prototype API.

## Performance correction and retained evidence

The first attempt was cancelled at 624.21 s when reaction startup performed
excessive disk work. Its transaction rolled back completely. Commit
`fa439b3` collects planner statistics before integrity joins and
derivation, then collects statistics for derived tables before consumers run.
The PostgreSQL/subset regression suite passed 141 tests, including a real
database test for statistics availability. Separate code review found no issue.
The corrected first reaction fetch took 0.564 s, compared with about 262 s
without returning rows in the cancelled attempt. Complete reaction derivation
took 338.96 s. The captured plan uses indexed lookups and incremental sorting.

Two independent verification passes completed successfully: the original grouped
payload check in 226.72 s total and the revised indexed exact-text check in
176.28 s total. The first completed before the attempted cancellation, and the
second also finished before cancellation reached it. Both logs are retained as
verification-attempt1.log and verification-attempt2.log. The revised check has
the same null/conflict semantics and was separately reviewed.

Server state and scripts are retained at
`/root/projects/omnipath-migration/pilot-20260930`. Files include
`release.json`, `requirements.txt`, `config.json`,
`preflight.json`, `attempt1.log`, `attempt1-result.json`,
`load.log`, `load-result.json`, `reaction-plan.json`,
`verification.log` and `verification.json`.
The schema must be fresh for any load; do not overwrite the committed pilot.
Nonsecret measurements and attestations are also saved alongside this report
in `nicesrv-parquet-pilot-results.json`.
