# Bulk build implementation and measurements — 6 September 2026

The default resource build now uses flat Arrow/DuckDB working tables and
set-based entity resolution. ChEMBL cache preparation and pypath source
mapping were left unchanged, as requested.

The backend described below was subsequently made the only build implementation.
`canonical/bulk.py` is now part of `canonical/match.py`; `bulk_writer.py` is now
`writer.py` (`ParquetWriter`). The old matcher, consolidation path and backend
selector were removed. See the [current system guide](../../packages/omnipath_build/README.md).
The timings below retain their original stage names.

## Implementation

- `canonical/bulk.py` normalizes identifier votes using the existing policy
  functions, then computes consensus, primary fallback, and symbol outcomes
  in SQL. Positive and negative lookup results persist in a disk-backed
  per-resource cache across batches. Conflicting identifiers are not resolved
  by arbitrarily selecting the first candidate.
- `bulk_writer.py` stores flat entity, identifier, entity-annotation, relation,
  and relation-annotation tables. SQL computes canonical entity hashes, remaps
  endpoints and payload links, applies symmetric orientation, and computes
  qualified relation keys. Shared Python policy functions still validate
  qualifiers and labels; ontology hashing uses the existing key semantics.
- Global endpoint labels and hierarchy counts are attached to flat tables before
  nesting, eliminating the wide post-aggregation join and second file rewrite.
- Complete reference aliases are expanded in SQL. Nested lists are constructed
  only during final aggregation, with evidence multiplicity and annotation order
  preserved. Aggregation remains partitioned to bound working memory.
- Default entity/record batch targets are 50,000, with a 256 MiB estimated-byte
  guard and a 50,000-relation guard. The administrative build source code uses
  the same new entity batch default. Explicit smaller batch requests still work.
- Forced collection after every flush/dataset was removed. Python automatic GC
  remains enabled, and one explicit collection remains before finalization.
- DuckDB working stages default to four threads, configurable with
  `OMNIPATH_BUILD_DUCKDB_THREADS`. Existing `--memory-limit` /
  `OMNIPATH_BUILD_DUCKDB_MEMORY` controls remain; the memory default stays 1GB.
- The comparison backend existed during measurement and was removed afterward.
  Builds now always use the flat SQL path; there is no backend selector.

These changes do not drop identifiers, weaken identity rules, reinstate May's
chemical filter, or remove source payloads. The implementation continues to use
the current mapper and Silver extraction semantics.

## Same-sample performance

Each run reads 100,000 molecules, 100,000 activities, and all 7,568 mechanisms.
All runs use nicesrv, warm source data, py-spy at 50 Hz, and the same 8G/10G
cgroup memory high/max guardrails. The resource outputs are isolated from the
live served data. These are individual instrumented runs, not statistical
confidence intervals or full-build forecasts.

| Measurement | Previous deferred implementation | New bulk implementation |
|---|---:|---:|
| Total wall time | 183.10 s | **132.23 s** |
| Entity resolution | 33.97 s | **19.11 s** |
| Forced garbage collection | 19.74 s / 46 calls | **0.19 s / 1 call** |
| Writer finalization | 25.02 s | **20.00 s** |
| Flushes | 45 | **7** |
| Process CPU, user + system | 242.57 s | **169.18 s** |
| Python process peak RSS | 1.67 GiB | **1.87 GiB** |

The total is **27.8% shorter** than the immediately preceding implementation
and **39.5% shorter** than the original warm baseline of 218.68 seconds.
The larger buffers trade some Python peak memory for throughput. Flat working
storage plus final temporary parts occupied about 365 MB at the measured
snapshot; this is not directly comparable to the old writer's chunk-only metric.

The new `append_observations` timer includes resolution and all Arrow/SQL working
table writes, totaling 40.12 seconds. It cannot be compared directly with the
old consolidation-only timer without adding the old separate chunk writes.

## Verification

- All 133 build/core tests passed locally and on nicesrv, plus three subtests.
- Additional resolver tests compare the SQL matcher to the existing matcher
  across conflicting identifiers, taxa, shuffled batches, and repeated calls.
- Complete file comparisons use bidirectional `EXCEPT ALL`, preserving row
  multiplicity and nested values. The 100k sample has zero added or removed
  rows across all 183,431 entities, 98,736 relations, and 206,989 payloads.
- The admin tests pass after the batch-default adjustment.

## Artifacts and reproduction

100k output version: `2026.9.6.105` under
`/root/projects/full_parquet/data/benchmarks/chembl-bulk-100k-20260906/`.
Logs, profile, events, summary, and equivalence comparison are under
`/root/projects/full_parquet/logs/chembl-bulk-100k-20260906/`.

## Larger sample and final memory fix

The 300k-per-dataset sample contains 607,568 source rows. The initial run
(version `2026.9.6.106`) completed ingestion but failed in the old final serving
projection at a 1GB DuckDB limit. A rerun at 2GB (version `2026.9.6.107`)
completed in **410.06 seconds**, producing 478,739 entities, 273,587 relations,
and 606,989 payloads. Its resolution stage took 61.81 seconds; activity mapping
alone took 110.38 seconds. Cache preparation remains outside this work.

The final implementation moves display metadata onto flat tables before
nesting. Replaying finalization from a closed working-database snapshot of the
successful 300k run completed at **1GB in 35.11 seconds**, with 1.84 GiB process
peak RSS. DuckDB's limit is not a total process-memory cap. Reference aliases
were already expanded in this snapshot, so compare this with 36.96 seconds
for the preceding finalization excluding reference expansion, not the entire
42.61-second close timer. The principal gain here is fixing the 1GB failure. Bidirectional `EXCEPT ALL`
found zero differences across all three output files versus version 107,
including nested values and evidence multiplicity.

The 132.23-second end-to-end measurement above predates this final projection
refactor; a complete end-to-end run of the final revision has not been measured.
The final revision passes all 133 build/core tests on both machines.

Larger-run logs: `/root/projects/full_parquet/logs/chembl-bulk-300k-2gb-20260906`.
Finalization replay, timings, and comparisons:
`/root/projects/full_parquet/logs/chembl-bulk-finalization-20260906`.
The API source default changed, but its running container was not rebuilt.

## Reproduction

Run a new isolated sample with:

```bash
TMPDIR=/path/on/disk omnipath-build build chembl \
  --version <new-version> --max-records 100000 \
  --output-dir /path/to/benchmark --memory-limit 1GB
```

Supply the same reference library through `OMNIPATH_LIBRARY_DIR` when the
isolated output root does not contain it. No unrestricted build or new live
resource release was published during this task.
