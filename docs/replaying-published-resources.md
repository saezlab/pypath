# Replaying published resources

`scripts/replay_resource.py` rebuilds resource versions from their existing
`evidence_payloads.parquet`, using the current native inputs_v2 mapper, resolver
and writer. It uses an existing immutable reference and never downloads source
files or builds a reference. Original artifacts stay untouched; output goes to a
new private directory for review before publication.

For the source-attribute migration, replay Rhea, Recon3D, Metatlas/Human-GEM,
KEGG, Reactome, MACdb and ConnectomeDB2025. Supply each exact published resource
version directory. All recorded datasets are included, including datasets with
zero published records.

```bash
uv run --frozen python scripts/replay_resource.py \
  --input-resource-dir /srv/published/resources/rhea/2026.9.5.17 \
  --input-resource-dir /srv/published/resources/recon3d/2026.9.8.1 \
  --output-dir /srv/private/replayed-release \
  --version 2026.9.30.1 \
  --library-dir /srv/reference/immutable-generation \
  --batch-workers 1 \
  --duckdb-threads 1 \
  --resource-ram-gb 3 \
  --memory-limit 512MB
```

Repeat `--input-resource-dir` for the other five resources; they run serially.
Metatlas also requires `--human-gem-xrefs /path/to/Human-GEM-metabolites.tsv`
(or `OMNIPATH_HUMAN_GEM_XREFS`) pointing to the existing cached auxiliary table.
The mapper reads that local file normally. Internet socket connections are
blocked in both the parent and spawned preparation workers. Download manager
and backend entrypoints are also blocked, including native pycurl transfer and
multi-handle paths that bypass Python sockets. Local file reads remain available.

The default has **no record limit**. `--max-records 20` is an optional per-dataset
cap for local tests; a capped result is a sample, not a full resource replacement.
The build reservation and replay deduplication memory limit are separate budgets.
DuckDB deduplicates with disk spill in a private temporary directory. An
ordered Parquet spool is written through a SQL COPY sink, so the raw SQL result
is never materialized for fetchmany. Python reads buffered batches of at most
32 records with column prefetch and memory mapping disabled, converting one
record at a time. Memory still depends on the largest record, page and dictionary.
The spool is removed when its iterator finishes or closes. Preparation uses native record and
byte batch limits. Reserve enough disk space for original artifacts, the private
DuckDB database/spill, observation shards and final Parquets.

Each original `dataset:<nonnegative integer>` row ID is replayed once despite
appearing under multiple payload owners. Gaps and original numeric suffix text
are preserved. Conflicting exact payload strings for one source row abort the
replay, even if their parsed JSON values agree. Unscoped records, mismatched
source/dataset names and mapped-record fallback payloads are rejected. Dataset
scope comes from the old manifest's `datasets` or `payload_origins`. For older
prototype manifests with neither recorded, the utility uses published row IDs
and all currently discovered datasets, recording that inference in its report.
Current datasets with no stored rows are replayed as empty. Published dataset
names without a current native mapper fail clearly; unknown historical empty
datasets cannot be recovered from the payload artifacts.

The adapter changes only the parent native pipeline's input pickle namespace
while a serial build runs. It replaces enumeration indices and reserialized JSON
with the original index and original JSON text. Spawned workers receive ordinary
raw dictionaries and load the real native mapper. Parent bindings restore after
success or failure; input mappers and pipeline files are unchanged.

The utility verifies that every selected original record was consumed and
republished with its exact row ID and UTF-8 payload text. Every output payload
body must remain nonnull, including duplicate owners of one source row. Losing a formerly
published record under the current mapper fails verification. It also checks
reaction evidence's SHA annotation against the original payload text. The native
build manifest records the immutable reference and current build provenance.
`replay_report.json` records original manifest/payload file checksums, original
reference metadata, dataset counts, phase times and the replay script checksum.
`replay_source_records/<source>.parquet` retains each replayed dataset, original
row ID and exact payload SHA-256 without copying raw JSON into the report.

Coverage means **every original record present in published payload artifacts**.
Original source rows that produced no published entity or relation were never
stored there and cannot be recovered by replay. On failure, the report explains
the error; any finished versions remain private and require review. The script
makes no production publication, release selection or PostgreSQL changes.
