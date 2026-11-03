# Build run log

A running record of individual real builds and loads (historically `make load` /
`make reload`; now `make build`, `make load-postgres`, `make finish-postgres`,
`make build-subsets`) against real data or a real database — not build *tuning* (see [`build-tuning.md`](build-tuning.md))
or pipeline *architecture* (see [`pipeline.md`](pipeline.md)), but the history
of actual runs: why each one happened, what parameters it used, and how long
its phases took. The goal is that nobody has to re-derive "why did we run
this" or "how long does this normally take" from scratch by scrolling terminal
history — it's written down once, here, per run.

## Format

One entry per run, newest first. Each entry:

- **Date / who / cycle**: when, and what spec/task motivated it (link the
  task id if there is one).
- **Reason**: why this run happened — what question it answers or what change
  it verifies. Not "ran a build" — the actual motivating question.
- **Parameters**: the exact invocation (`SOURCES`, `MAX_RECORDS`, `LOAD_JOBS`/
  `--stage-jobs`, `THREADS`, `BATCH_SIZE`, target database). Enough to
  reproduce the run verbatim.
- **Phase durations**: wall-clock time per major phase (per-source delete,
  preparse, stage/canonicalize, copy, derive), and any single-dataset outlier
  worth calling out by name. Exact numbers where the log gives them; honest
  ranges/estimates where it doesn't.
- **Outcome**: what changed, what it confirmed or refuted, links to the
  acceptance query results if the run was measured against
  `contracts/coverage-acceptance.sql`.

Keep entries factual and specific — numbers and log excerpts, not vibes.

---



## 2026-10-01 — batched COSMOS writes, COSMOS-only retry attempt 13 (complete; bounded manual review passed)

**Reason**: user-authorized lightweight optimization after bounded backend
samples showed repeated single-row inserts and client waits. COSMOS now batches
edge and label writes through `executemany`, with buffers capped at 1,000 rows
or distinct label keys. Pending labels keep their first occurrence; the existing
database conflict rule handles repeats across batches. The named read cursor
fetches outside each automatic write pipeline. Projection, insertion order,
schemas, indexes, statistics SQL and product transaction boundaries are unchanged.
Full-release performance improvement remains unmeasured.

**Parameters**: runtime `8a163fb5a548ebd8a31c4b80db16b0c0e48da5b2`, branch
`parquet-migration`, pushed to GitHub and fast-forwarded into the clean dedicated
nicesrv checkout. pypath remains `e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`.
Release `2026.9.30.4`, schema `full_20260930`, 46 fixed pins, manifest SHA-256
`06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`.
Invocation: `/root/.local/bin/uv run --frozen --no-sync --directory
/root/projects/omnipath-migration/omnipath-build python
/root/projects/omnipath-migration/full-release-20260930/omnipath-cosmos-retry-attempt13.py
8a163fb5a548ebd8a31c4b80db16b0c0e48da5b2`.
Calls only `build_subsets(..., products=('cosmos',), checkpoint_products=True,
on_product_committed=product_committed)`.

Attempt12's completed MetSigDB and network views remain committed at runtime
`8cf850e72da9cdb55d321fa2c22adf5c4f568428`, with durations **1602.020031 s** and
**0.014368 s**. Their exact metadata/checkpoint are retained separately under
`subset_retry.retained_products`; only newly committed COSMOS belongs to this
retry's `committed_products`. The combined checkpoint is updated after actual
COSMOS commit. Attempt8's base/runtime provenance remains
`05a02a56eaf06980d420a2010d11079dd026f0dd`; cancelled retries9/10/11/12 remain in
history and archives. No resource build, resolution, COPY, shared derivation,
MetSigDB or network-view rebuild occurs.

Private PG5440: work_mem 16 MB, maintenance_work_mem 128 MB, gather 1, JIT off;
container 4 GiB / CPU 2 / WAL 8 GB. Builder/verifier MemoryMax 3 GiB, swap
256 MiB, CPUQuota 200%, Nice 10, stop timeout 90 s, disk floor 64 GiB.
Build unit `omnipath-migration-subsets-attempt13-20261001.service`, invocation
`c9b85a02f07d4126a126a6cb70e8296d`; queued verifier
`omnipath-migration-full-verification-attempt13-20261001.service`, invocation
`3351c3a8caf14604ac301dfdf3bf7ee0`.
Retry helper SHA-256
`9c2eb1f0dc4614b0affe57c31e1aea2df9bf6da9391132a38771cd6f7b83d64f`;
startup helper SHA-256
`16a2b7af503ff6bf229ce9a253932b3314f1da78ce3880f1098280e1578f86a3`;
monitor SHA-256
`b33c1e53039bfd4634b739107c8b9373afecf94ea93a24ce3b3e242d39089f64`;
verifier SHA-256
`9ee43a4fcadf71445693785c0f73460317240f13b35f66ed624e0e7594d323f2`;
final safety helper SHA-256
`0b02c1372e21767c36dbef3dc68788ae0b0db2c188cdc7a354934ebb5811b326`.

**Phase durations**: started **2026-10-01T16:05:40.901768+00:00** (18:05 CEST).
The first monitor snapshot showed COSMOS running at retry elapsed **130.7 s**,
the two retained products committed and verification waiting. Cumulative
executed build time before attempt13 is **30471.239039 s**, excluding gaps,
probes and tests. The previous cancelled COSMOS transaction and its rollback
checks are recorded in the attempt12 entry below; archive report SHA-256 is
`91e6b596d9aa0e4808e8e19c9608d836372dccced3d6105eab3264e34f83b8dc`.

**Validation/outcome**: focused COSMOS/release-guard tests passed (16 tests,
3.29 s). **609 relevant PostgreSQL/subsets/core tests passed in 35.45 s**, with
one existing Starlette deprecation warning. Real PostgreSQL tests compare
complete edges including generated IDs/order, labels and non-timing statistics
against serial transport, exercise named-cursor fetch/write interleaving and
first-label behavior within/across batches, and confirm rollback after an error
in a later batch. Ruff, formatting and diff checks passed. Root source/test
review and independent restart-helper review passed.

**COSMOS committed; scheduled check at 16:25 UTC**: COSMOS completed and
committed at **2026-10-01T16:14:15.095977+00:00** (18:14 CEST), with exact product
duration **514.173461 s** (8 min 34 s). Retry wall-clock time was **514.195215 s**;
the coordinator's nested phase wrapper reports 514.172231 s. A fresh read-only
private PostgreSQL check confirmed all three `subset_build_metadata` records
match the release, manifest and combined report statistics. MetSigDB/network
views retain their original attempt12 timings and provenance. Successful product
durations sum to **2116.207860 s** across the two retries. Cumulative executed
attempt8-through-attempt13 time is **30985.434254 s**, including cancelled
subset work and excluding gaps/probes/tests. The preceding COSMOS attempt was
cancelled after approximately1341.07 s without completing; it does not provide
a completed serial-build baseline or a precise speedup ratio.

The builder is inactive after success; the verifier is active in full-release
verification. Private data remains restored at this snapshot. PostgreSQL-only
rollback verification and final safety remain pending. The 15-minute monitor
continues to follow attempt13, including retained versus new product provenance.
Completion requires both verification reports and the read-only
`omnipath-full-final-safety-attempt13.py` checks to pass, including all three
product metadata/report matches, restored private data and protected-service
baseline matches. This was the pre-review snapshot; final acceptance is recorded below.


**Final outcome — user-requested manual review, 1 October 2026, 19:50 CEST**:
the user removed mandatory post-build verification, the rebuild-and-rollback
check and the historical final-safety gate. Normal build/rebuild commands and
loader validation/constraints are retained. The optional exhaustive verifier
was stopped at **2026-10-01T17:48:07.285832+00:00**, after **5605.098408 s** in
full verification (**5605.846723 s** coordinator wall time). Its final status is
cancelled/KeyboardInterrupt, preserved as a failed optional audit rather than
rewritten as a successful exhaustive result. The queued rollback rebuild never
ran. Stop time was **0.339353 s**. Archive
`attempt13-verification-stopped-20261001` preserves the successful build report
(SHA-256 `f5acf91dc233ceaf98e0ec56d32d18b12216f7ac23b8e8d487cabfa932e6b351`),
cancelled verification report (SHA-256
`11417ca33ad63732ab2a1cc659e8f041c109957b60a7fb89dfd704ba2c6e3757`),
and log (SHA-256
`26c6926bbf63b78f770d0b06f3f87e1c4f8c1d97abef9fa305e9ec98de63b798`).
Fresh private metadata snapshots before/after stopping confirmed every committed
product and the base checkpoint unchanged, both pilots preserved, private data
restored and no remaining verifier backend.

A one-off read-only sample then **passed in29.085445 s**, at unchanged runtime
`8a163fb5a548ebd8a31c4b80db16b0c0e48da5b2`. Unit
`omnipath-migration-manual-sample-20261001.service`, invocation
`164d447d1ed2449980ca78a0440639c2`. Parameters: frozen/no-sync existing UV
environment, literal private PG5440, repeatable-read/read-only transaction,
60 s per-statement timeout, work_mem16 MB/JIT off, MemoryMax2 GiB/swap256 MiB,
CPUQuota100%, Nice10. The sample read at most20 entity and20 relation records
per resource, spread over first/middle/last row-group fronts where available,
without scanning the remaining records. Selected normalized rows were capped
at20 per table/resource. Existing loader integrity results and complete index/
constraint catalog definitions were checked; source parsers, resolver, resource
builds, raw payloads, COPY and derivation were not run.

Coverage: **46 resources;920 entities;840 relations;920 identifiers;840 evidence
occurrences;893 annotation occurrences;130 typed measurement projections;
145 ontology edge orientation checks;15 reaction contexts/84 participants**.
Reaction checks included53 stored coefficient assertions and49 compartments,
direction/type/source/hash links, member ordinals and diagnostics. Product checks
covered **100 MetSigDB memberships across all five sources**, five actual set
sizes,20 COSMOS edges spread through the edge-ID range, and three small network
query pages (LIANA, MetalinksDB, reactions). Product metadata matched release
`2026.9.30.4`, its exact manifest and the original combined build statistics.
This is a bounded value/provenance sample, not an exhaustive value comparison,
closure proof, fresh full-table recount or legacy scientific-parity claim.

The manual service completed in30.382 s, CPU time12.894 s, observed peak memory
approximately1.9 GiB with zero swap. The script SHA-256 is
`168b8905e4365bcbba5124b99d17ff79f99460e83e43f55dbfcd7187d5cdd187`;
sample report SHA-256 is
`87fb2b1edcc653031f0041e5a28d4393db9feadec88380594941b14e348a0ef9`.
The one-off completion snapshot at **2026-10-01T17:50:49.778349+00:00** confirmed
all four protected services match baseline identity/start/running state, private
data is restored, and no build/verifier/sample session remains. It is not a new
mandatory phase. `migration-completion-20261001.json` records acceptance as a
successful committed build with the requested bounded manual review. The canonical
server status helper reports removed gates explicitly; the historical named
attempt13 helpers and reports remain intact. The scheduled monitor is **paused**.

All46 resources, base tables/late constraints/query indexes, shared derivations
and three subset products remain committed. Exact successful product timings
are MetSigDB1602.020031 s, network views0.014368 s and COSMOS514.173461 s.
Recorded attempt8-through-attempt13 build time is30985.434254 s including
cancelled build work, excluding gaps/probes/tests and the optional verification.
The private PostgreSQL build milestone is complete under the user's revised
acceptance criteria; production/prototype services remain unchanged.


Public workflow documentation was committed and pushed on parquet-migration
as f890899 (Document optional post-build audits). Only the PostgreSQL README
was included; private build logs/reports remain local. Server runtime remains
8a163fb5a548ebd8a31c4b80db16b0c0e48da5b2; the documentation commit changes no
loader, derivation or subset implementation.


## 2026-10-01 — ClassyFire final join and product commits, subset-only retry attempt 12 (COSMOS cancelled; MetSigDB and network views retained)

**Reason**: user-authorized ClassyFire final-join optimization and commits after
each complete product (MetSigDB, network views, COSMOS). A SELECT-only probe
capped at 20 direct assignment edges produced 105 memberships with 20 projected
chemicals. The original publication join rescanned the projection 105 times;
the guarded full join used a hash match and one scan, with identical complete
row multisets. This does not measure full publication time. The variant is
ClassyFire-only; the two-sided non-null guard preserves inner-match semantics
for duplicate, unmatched and null keys without reducing the full join.

**Parameters**: revision `8cf850e72da9cdb55d321fa2c22adf5c4f568428`, branch
`parquet-migration`, pushed to GitHub and fast-forwarded into the clean dedicated
server checkout. pypath remains `e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`.
Release `2026.9.30.4`, schema `full_20260930`, 46 fixed pins, canonical manifest
SHA-256 `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`.
Invocation: `/root/.local/bin/uv run --frozen --no-sync --directory
/root/projects/omnipath-migration/omnipath-build python
/root/projects/omnipath-migration/full-release-20260930/omnipath-subset-retry-attempt12.py
8cf850e72da9cdb55d321fa2c22adf5c4f568428`.
Calls only `build_subsets(..., products=('metsigdb','network_views','cosmos'),
checkpoint_products=True, on_product_committed=product_committed)`.

Each complete product's tables and `subset_build_metadata` commit together;
a later failure or cancellation rolls back only the current product. One schema
session advisory lock spans all product transactions and callbacks. The exact
single release ID/digest is checked before and after each product, with a table
lock protecting the release marker inside the transaction. The callback records
`committed_products` and a `subset_checkpoint` only after actual commit.
`product_prepared` and per-source completion events are not publication commits.
Explicitly select remaining products to resume; database metadata is authoritative
if observer reporting fails after a commit. Default API/CLI behavior stays atomic
unless the checkpoint flag is enabled; it is enabled in this run.

Other four publication queries, KEGG-only `source_edges AS NOT MATERIALIZED`,
all `projected AS MATERIALIZED` aggregates and normal planner settings remain
unchanged. No identifier resolution, resource build, COPY or shared derivation
is repeated. Base/runtime provenance remains attempt8
`05a02a56eaf06980d420a2010d11079dd026f0dd`; cancelled attempts9/10/11 are retained
in `subset_retry_history` and their immutable report/log archives.

Private PG5440: work_mem 16 MB, maintenance_work_mem 128 MB, gather 1, JIT off;
container 4 GiB / CPU 2 / WAL 8 GB. Builder/verifier MemoryMax 3 GiB, swap
256 MiB, CPUQuota 200%, Nice 10, stop timeout 90 s, disk floor 64 GiB.
Build unit `omnipath-migration-subsets-attempt12-20261001.service`, invocation
`4a1002cea5fb4a09814513532fdcf64a`; queued verifier
`omnipath-migration-full-verification-attempt12-20261001.service`, invocation
`49aacc10a680461bbf2bff259a9d673a`.
Retry helper SHA-256
`7519837fbfde4915502c84a4280d163f7929574f7ba94b98e2105b26043702ec`;
startup helper SHA-256
`904647f5182e0d7e6ff3c59bf03f4a242849e861522c13cb2284609177bdec2e`;
monitor SHA-256
`ef69166c97c409243b45b5b51a74804c280f909a22933c5b60d30fbe1252d5d0`;
verifier SHA-256
`2876b92cb07b50000eca5b3598af287670a600a5dc8d2ed9c139f1544c2367c5`;
final safety helper SHA-256
`17d267e28a4bf57c1409adbdcec9bb2c158c02f551939dfc57edf53358c4693a`.

**Phase durations**: started **2026-10-01T15:12:44.561620+00:00** (17:12 CEST).
At the startup snapshot, Reactome finished in **24.239005 s** and WikiPathways
in **32.004826 s**; KEGG is running. Source timings run to the next source start;
the last source includes MetSigDB index/statistics work. Exact product durations
will be reported separately. No complete product had committed at this snapshot;
ClassyFire, all product commits and both verification phases remain pending.
Cumulative executed build time retained before attempt12 is **27528.094099 s**,
excluding gaps, probes and tests.

**Latest scheduled check (15:30 UTC)**: KEGG source completed in
151.292255 s and MACdb in 49.254894 s; Reactome/WikiPathways remain
24.239005/32.004826 s. ClassyFire started at retry elapsed 256.833 s and is
still running. A bounded private backend snapshot found an active query at
899.9 s, waiting on DataFileRead, with zero lock blockers. No complete product
has committed; both units are active and verification remains queued. The
committed base is retained; disk reserve was 329.0 GiB at the monitor snapshot.
This is progress/wait-state monitoring, not a full-plan diagnosis or a completed
ClassyFire performance result.

**Durable product progress (1 October, 17:39 CEST)**: MetSigDB committed at
2026-10-01T15:39:26.611509+00:00, exact product duration **1602.020031 s**
(26 min 42 s). Network views committed at 15:39:26.630790+00:00 in
**0.014368 s**. A fresh read-only private PostgreSQL check confirmed both
product metadata records match the release, manifest and reported stats.
COSMOS started at retry elapsed 1602.071 s and is running; both verification
steps are still queued. All five MetSigDB source stages completed:
Reactome24.239005 / WikiPathways32.004826 / KEGG151.292255 / MACdb49.254894 /
ClassyFire1345.211318 s. The ClassyFire coordinator time includes final
MetSigDB index/statistics work; it is not a query-only timing. The private data
root is present and the hidden-verification directory absent. MetSigDB and
network views are now preserved independently if COSMOS later fails.

**COSMOS cancellation**: after the user authorized batching writes and reducing
duplicate label attempts, reviewed helper `omnipath-stop-cosmos-attempt12.py`
SHA-256 `47fb89a70bdf191ead8602770f8c3d2bfb3ef57f511564b637ebb50163835085`
stopped the queued verifier, then only the current COSMOS transaction. Safety
confirmation at **2026-10-01T16:01:48.061377+00:00** (18:01 CEST) showed both
MetSigDB/network views product metadata/stats and partial checkpoint unchanged,
completed base/shared derivations with46 pins, both pilot schemas preserved,
four protected services unchanged, no importer/verifier backend and absent
uncommitted COSMOS tables after rollback. Private data remained restored.
Attempt12 elapsed **2943.144940 s**, including approximately **1341.07 s**
(22 min21 s) of cancelled COSMOS work. Retained product timings remain
MetSigDB1602.020031 s and network views0.014368 s. Cumulative executed build time
is **30471.239039 s**, excluding gaps/probes/tests. Archive
`attempt12-before-subset-retry-20261001` preserves reports/logs; report SHA-256
`91e6b596d9aa0e4808e8e19c9608d836372dccced3d6105eab3264e34f83b8dc`.
Only COSMOS will be selected for the replacement run; no import, resource build,
shared derivation, MetSigDB or network-view rebuild is required.

**Validation/outcome**: 24 focused MetSigDB tests passed in 4.63 s; 25 focused
orchestration tests passed in 4.46 s. **604 relevant PostgreSQL/subsets/core tests
passed in 39.65 s**, with one existing Starlette deprecation warning. Regression
coverage includes complete five-source result multisets, real PostgreSQL projection
scan plans, duplicate/null/unmatched/empty-side semantics, current-product rollback
on errors/interrupts/SQL cancellation, earlier commit durability, release changes,
schema-lock retention/release, callback timing/snapshot isolation and explicitly
resuming only remaining products. Ruff, formatting and diff checks passed.
Root source review and independent stop/retry/startup/final-safety review passed.

The attempt11 cancellation and archive are recorded below. At startup, safety confirmed
all four protected services matched their baseline identities/start times/running
states and private data remained restored. Both attempt12 units were active and
verification was queued. They were subsequently stopped as recorded above;
COSMOS and the existing verifications continue under attempt13. Production
services remain untouched.


## 2026-10-01 — indexed KEGG joins, subset-only retry attempt 11 (cancelled; base and derivations retained)

**Reason**: user-authorized fix after a read-only diagnosis of KEGG's
reaction-to-pathway pairing. PostgreSQL estimated one chemical edge and two
pathway edges; actual counts were 52,622 and 111,787. The nested loop repeatedly
filtered the 641,065-row materialized source-edge result, implying about
5.9 billion pair comparisons. The exact 20-chemical-edge join probe took
12.402 s and read 3,725,268 temporary blocks (about 28.4 GiB). These are logical
temporary reads, not distinct bytes read from physical storage.

The KEGG-only `source_edges AS NOT MATERIALIZED` variant lets predicates and
reaction-key joins reach the existing base-table indexes. A full pairing-only
probe produced all 481,870 candidates in 13.709 s, with zero temporary reads
or writes and parameterized `relations_subject_idx` lookups. An earlier cold
probe hit its 15 s limit; the successful repeat used a 30 s limit. A diagnostic
hash join counted the same candidates in 4.582 s, but discouraging nested loops
for the whole publication query would sort about 26 million global identifiers.
The deployed variant keeps normal planner settings and indexed alias lookups.
These probes do not measure complete MetSigDB publication.

**Parameters**: tested revision `f14e703d2c7ee6ac083bcd7b05ee43a20f665624`,
branch `parquet-migration`; only KEGG selects the new source-edge variant.
The other four source queries and all `projected AS MATERIALIZED` aggregates
are retained. Release `2026.9.30.4`, schema `full_20260930`, 46 exact pins,
manifest SHA-256 `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`.
Base and shared derivations retain `05a02a56eaf06980d420a2010d11079dd026f0dd`;
pypath stays at `e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`.
Invocation: `/root/.local/bin/uv run --frozen --no-sync --directory
/root/projects/omnipath-migration/omnipath-build python
/root/projects/omnipath-migration/full-release-20260930/omnipath-subset-retry-attempt11.py
f14e703d2c7ee6ac083bcd7b05ee43a20f665624`.
Calls only `build_subsets` for MetSigDB, network views and COSMOS in one
transaction. Retains the committed COPY, keys/indexes and shared derivations;
does not rebuild resources or resolve identifiers. Private PG5440,
work_mem 16 MB, maintenance_work_mem 128 MB, gather 1, JIT off;
container 4 GiB / CPU 2 / WAL 8 GB. Builder/verifier MemoryMax 3 GiB,
swap 256 MiB, CPUQuota 200%, Nice 10, stop timeout 90 s, disk floor 64 GiB.

Build unit `omnipath-migration-subsets-attempt11-20261001.service`, invocation
`d4f9044fbd8e428b9afedc8dd5998292`; queued verifier
`omnipath-migration-full-verification-attempt11-20261001.service`, invocation
`3bf00cfa022f496d8b27370b456eeaf3`. Retry helper SHA-256
`7ac5a26fbea6c59c07611f14765c2f8b22a21014cd442fa8f53222524a14776f`;
guarded startup SHA-256
`d78130fa1586c8bc5e6115a0b2cd0d87e82c4d6fe379818c5733def5fe881347`.
Canonical provenance preserves attempts 9 and 10 under `subset_retry_history`,
their archived reports and the original attempt8 base phases/counts.

**Phase durations**: started **2026-10-01T14:00:26.937680+00:00** (16:00 CEST).
Reactome completed in 25.955178 s; WikiPathways in 30.160086 s;
KEGG in 187.541608 s; MACdb in 50.380578 s. The retry was cancelled at
**2026-10-01T15:07:55.557006+00:00** after **4048.619120 s** total retry work
(67 min 29 s), including approximately 62 min 35 s in unfinished ClassyFire.
No complete product committed. Completed source work rolled back with the
single subset transaction; base and shared derivations stayed committed.
Source coordinator times run through the next source start; product timings
are unavailable because MetSigDB did not finish. Cumulative executed build time
is 27528.094099 s, excluding gaps and diagnostic/test runs.

**Cancellation**: user authorized a ClassyFire-only final join optimization
and commits after each complete product for the replacement run. Reviewed helper
`omnipath-stop-subsets-attempt11.py` SHA-256
`f5a91e215667920dedc8ee1681cf07d32dff5e07a2b8c1b400b65e928dfcd35a`
stopped the waiting verifier first, then only the private subset builder.
Safety confirmation at **2026-10-01T15:07:55.819329+00:00** verified both pilot
schemas, the completed base checkpoint, release metadata and 46 pins, four
protected service identities/start times/running states, zero remaining
importer/verifier backends and absent subset metadata. Private Parquet data
remained restored. Archive `attempt11-before-subset-retry-20261001` retains
the report and logs; report SHA-256
`9417f5d768738f03f2aec9779a0dbd3d1f925676a7016ac43bccfaecfec393cb`.
No COPY, shared derivation, resource build or production service change occurred.

**Validation/outcome**: 18 focused tests passed in 3.63 s; 579 relevant
PostgreSQL/subsets/core tests passed in 30.10 s (existing Starlette deprecation
only). Complete SELECT row multisets match across all five source rules.
The additional KEGG fixture has 15 original entity/relation records and covers
both pathway orientations, two reactions, duplicate candidate memberships,
input/output/participant paths, excluded proteins and all relation/via
provenance. Seven exact candidate tuples yield four unchanged published rows;
numeric source row 2 remains the selected evidence. Root source review,
independent cancellation/retry/startup/final-safety review, Ruff and diff checks
passed. The tested commit was pushed and the clean dedicated server checkout
fast-forwarded; pypath stayed fixed.

Attempt10 cancellation was confirmed **2026-10-01T13:56:12.027050+00:00** after
2426.650620 s total retry work (about 39 min 28 s in unfinished KEGG).
Reactome/WikiPathways subset work rolled back with the same subset transaction.
Completed base checkpoint, release metadata, all 46 pins, both pilot schemas and
all four protected service identity/start/running states were retained; no old
importer/verifier backend remained. Archive `attempt10-before-subset-retry-20261001`,
report SHA-256 `46ec0919f335fa1fb9642491867d1e3245c810ef7f733d7c75918d87556fe3cb`.
Cumulative executed build time before retry11 was 23479.474979 s, excluding
gaps and diagnostic/test runs. Attempt11 is stopped. The replacement subset-only run and its product checkpoints
will be recorded separately; full/PG-only verification and final safety are pending.


## 2026-10-01 — materialized MetSigDB projection, subset-only retry attempt 10 (cancelled; base and derivations retained)

**Reason**: user-authorized follow-up after live CPU/plan diagnosis found repeated
identifier aggregation under the final nested loop. The estimated membership
count was1 versus24301 actual Reactome memberships/2178 metabolites. A bounded
SELECT-only probe of20 outer memberships/20 projected chemicals showed aggregate
loops20 while underlying alias retrieval ran once and sorted input was reused.
Forcing `projected AS MATERIALIZED` reduced aggregate executions to1. Probe
outputs were0 because independently sampled chemicals/memberships did not overlap;
the probes establish executor behavior, not result parity or full-scale speed.

**Parameters**: tested revision `79b64a1c6ec8f983047629be1424af52d852b2a1`,
branch `parquet-migration`; production SQL change is only `projected AS MATERIALIZED`.
Same release `2026.9.30.4`, schema `full_20260930`,46 exact pins and canonical
manifest SHA-256 `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`.
Base/derivation revision remains `05a02a56eaf06980d420a2010d11079dd026f0dd`;
pypath remains `e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`.
Invocation: `/root/.local/bin/uv run --frozen --no-sync --directory
/root/projects/omnipath-migration/omnipath-build python
/root/projects/omnipath-migration/full-release-20260930/omnipath-subset-retry-attempt10.py
79b64a1c6ec8f983047629be1424af52d852b2a1`. Calls only `build_subsets` for
MetSigDB, network views and COSMOS, atomically. No Parquet/resource build,
resolution, COPY, key/index recreation, `finish_release` or shared derivation.
Private PG5440; work_mem16MB, maintenance128MB, gather1, JIT off;
container4GiB/CPU2/WAL8GB; build/verifier MemoryMax3GiB, swap256MiB,
CPUQuota200%, Nice10, stop90s, diskfloor64GiB. Same artifact/cache/reference selection.

Build unit `omnipath-migration-subsets-attempt10-20261001.service`, invocation
`cabf4b1221dd4a83b6326939bf8ca265`; queued verifier
`omnipath-migration-full-verification-attempt10-20261001.service`, invocation
`b369dd1c4adf4dc597083a0037713bbe`. Retry helper SHA-256
`43dc6eb71b87c8fc4ff5c6835928532056c70af43aa566612e7f3aedbda4465a`.
Canonical report retains all base8 fields and the cancelled9 metadata under
`subset_retry_history`; new retry10 has its own HEAD/start/timings and archive hash.

**Phase durations**: started **2026-10-01T13:15:45.151003+00:00** (15:15 CEST).
Reactome completed in25.152219 s; WikiPathways33.578872 s. KEGG was cancelled
after about39min28s, with total retry time2426.650620s; only subset work rolled back.
Source timers include publication/counts through the next source start. Final
ClassyFire also includes index/statistics work; product phase_seconds capture
complete builder timings. Remaining sources/products, commit and verification
remain pending. These measured source completions are a full-scale improvement
over the previous unfinished first query, without claiming a completed migration.

**Validation/outcome**: focused16 passed in3.78s; relevant PostgreSQL/subsets/core
577 passed in30.38s (existing Starlette deprecation only). Complete publication
SELECT multisets match across all five fixture rules; every synthetic source has
at most20 original records. Forced nested-loop regression observes old projection
aggregate loops2, materialized producer1, reusable CTE scans2. Root SQL review,
independent coordination review, Ruff/format/diff checks passed. No other SQL
expressions, eligibility, identifiers, labels, taxonomy, set sizes, qualifiers or
provenance were changed.

Attempt9 cancellation was verified **2026-10-01T13:13:23.747845+00:00**, after
4095.859603 s of unfinished Reactome work. Only its subset transaction rolled
back; completed checkpoint/release/all46 pins, copied counts, base indexes,
shared derivations and both prior pilot schemas remain. No old importer/verifier
backend remained; all four protected container identity/start/running states
matched baseline. Archive `attempt9-before-subset-retry-20261001`, report SHA-256
`f79a6ebf133f6cc2d51043da58c3b081dcb096fd6d963aa5228e2e9f41ace5bc`.
The original8 archive/hash is retained and checked recursively. Fresh retry
HEAD/import-path/clean-checkout guards passed. Attempt10 was archived on cancellation;
the15-minute monitor and queued verification now follow attempt11, recorded above.


## 2026-10-01 — MetSigDB optimization and subset-only retry on nicesrv, attempt 9 (cancelled; base and derivations retained)

**Reason**: user-authorized change after SELECT-only diagnosis showed global
canonical-entity sorts and identifier deduplication during MetSigDB enrichment.
Select canonical fields only for chosen entity keys, and look up aliases across
the release pins using the existing resource/version/entity identifier PK.
No additional index, Parquet change, resource build or identifier resolution.

**Parameters**: tested subset revision
`48c6e14325ff831173d95a8eaba1a8d8ba5b920c`, branch `parquet-migration`;
unchanged release `2026.9.30.4`, schema `full_20260930`, exact46 pins and manifest
SHA-256 `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`.
Retained base/derivation runtime is
`05a02a56eaf06980d420a2010d11079dd026f0dd`, pypath remains
`e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`. Invocation:
`/root/.local/bin/uv run --frozen --no-sync --directory
/root/projects/omnipath-migration/omnipath-build python
/root/projects/omnipath-migration/full-release-20260930/omnipath-subset-retry-attempt9.py
48c6e14325ff831173d95a8eaba1a8d8ba5b920c`. Calls only
`build_subsets(dsn, 'full_20260930', products=('metsigdb','network_views','cosmos'))`;
the three products publish in one transaction. No `load_release`, `finish_release`,
COPY, schema/key recreation or shared derivation is invoked. Private PG5440,
work_mem16MB, maintenance_work_mem128MB, gather1, JIT off; private PostgreSQL
4GiB/CPU2, WAL8GB and durability settings unchanged. Both systemd units retain
MemoryMax3GiB, swap256MiB, CPUQuota200%, Nice10, stop90s and64GiB disk reserve.

Build unit `omnipath-migration-subsets-attempt9-20261001.service`, invocation
`33cd80deaf204fb5867a95dd52dbb9c5`; queued verification unit
`omnipath-migration-full-verification-attempt9-20261001.service`, invocation
`f770be27ddd94bf1b4fee698625ca9fe`. Retry coordinator SHA-256
`a6510f1a2bc60116ccfa39768b3ddec5f879c75a5829efbe5d0a95d11fc8806c`.
The canonical report retains original attempt8 timings/code/counts, and records
the new revision and progress separately under `subset_retry` attempt9.

**Phase durations**: started **2026-10-01T12:05:07.590585+00:00** (14:05 CEST).
Cancelled **2026-10-01T13:13:23.747845+00:00**, after4095.859603 s without
completing the first Reactome source. Only the subset transaction rolled back.
The prior attempt8 timings/counts/derivations and its archive remain retained;
cumulative executed build time before retry10 was21052.824359 s, excluding
between-run gaps. This run supplied the repeated-aggregation diagnosis recorded
in the newer attempt10 entry; product publication and verification did not run.

**Validation/outcome**: focused MetSigDB integration14 passed in3.52s; PostgreSQL,
subsets and core suite575 passed in30.61s, with one existing Starlette deprecation.
Root SQL review, independent coordinator review, Ruff/format/diff checks passed.
Regressions cover canonical winners, aliases from another pin, duplicate aliases,
valid canonical InChIKey preference, aliasless/name-only chemicals, scoped set
names and source taxon. An initially wider packages suite was interrupted in an
unrelated reference-builder test after226 passed in270.92s; it is not a complete
suite result. No real upstream resources were rebuilt.

Cancellation verified the completed base checkpoint/release,46 pins, no remaining
importer/verifier backends, absence of committed subset metadata, both prior pilot
schemas and all four protected containers unchanged at
**2026-10-01T12:03:12.509258+00:00**. The original report/log are archived privately
under `attempt8-before-subset-retry-20261001`; original report SHA-256
`3f55e0fda4b231e703486e24c6bf5bdfc01cbddd4d222399009acb95edd9766f`.
Runtime import-path/HEAD and clean-checkout guards passed at retry initialization.
Full-release comparison and PostgreSQL-only rebuild/rollback remain queued, followed
by `omnipath-full-final-safety-attempt9.py`; completion and full-scale speed claims
remain pending. The existing15-minute heartbeat now tracks the newer attempt10 units.


## 2026-10-01 — full resolved-Parquet release on nicesrv, attempt 8 (base/derivations committed; subsets cancelled)

**Reason**: user-authorized restart after attempt 7's confirmed rollback. Move
five copied-table PKs, annotation-owner uniqueness and all seven FKs after COPY,
fully validating before base checkpoint; increase only private PG's WAL budget.
Reuse every selected artifact, cache, reference and exact resource version.

**Parameters**: exact checkout `05a02a56eaf06980d420a2010d11079dd026f0dd`,
late constraints `05a02a5`, checkpoint loader `552e1d0`, frontier ontology
`4fe91d1`, replay `300d06b`, pypath
`e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`. Same46 pins/release `2026.9.30.4`,
schema `full_20260930`, canonical manifest SHA-256
`06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`.
Invocation: `uv run --frozen --no-sync --directory
/root/projects/omnipath-migration/omnipath-build python
/root/projects/omnipath-migration/full-release-20260930/postgres-sql-build.py`.
`load_release(..., batch_size=1024, duckdb_threads=1, memory_limit='1024MB',
validate_source_records=False, checkpoint_base=True, defer_constraints=True,
temp_directory=state)` then all three products. CSV16MB/COPY1MiB unchanged;
work_mem16MB, maintenance_work_mem128MB, max_parallel_workers_per_gather1,
JIT off. Private PostgreSQL16 remains4GiB/CPU2; max_wal_size8192MB loaded via
reload with no durability setting change/restart. Importer/verifier each
MemoryMax3GiB, swap256MiB, CPUQuota200%, Nice10, stop90s, diskfloor64GiB.
Import unit `omnipath-migration-full-postgres-sql-attempt8-20260930.service`,
invocation `5c27deab85d34cba864863a9cd15d82e`; verifier unit
`omnipath-migration-full-verification-sql-attempt8-20260930.service`, invocation
`ed7534f2519a430db493997e366dde4e`. Coordinator SHA-256
`af3cc87924c76da1b14db113151eb31364bab2903d7db394e04fc850b4ae5de4`.

Column NOT NULL/CHECKs and annotation sequence remain active during COPY. The
exact PKs and owner uniqueness are created afterward, narrow join-key statistics
collected, and every FK validated before integrity/query indexes/file recheck and
base commit. Derivation publishes in a separate transaction; a committed base
allows PostgreSQL-only finish without repeated COPY. Requested raw audits remain
available and run after lookup keys in deferred mode. 527 regression tests passed
in28.39s; exact real-PG catalog/data comparisons, all key/FK negative cases,
checkpoint retries, style and independent review passed.

**Phase durations**: started **2026-10-01T07:20:35.203142+00:00**. All46
resources/all230 projections completed: typed validation130.514596 s,
CSV staging1039.734568 s, COPY2325.139185 s. Late keys/FKs and join-key
statistics1284.842577 s; normalized integrity1430.471222 s; ordinary query
indexes840.228613 s; final release/file recheck6.527422 s. Base checkpoint
committed at elapsed7078.557093 s. Source-scoped ontology448.249527 s,
global ancestor351.502475 s, ontology total1346.177070 s (nested in shared
derivation2951.443228 s); final full statistics12.403613 s. Post-checkpoint
finish total3196.274258 s includes checkpoint validation. Shared derivations
and release metadata committed successfully before subsets_start at
elapsed10274.893 s. User-authorized cancellation during first MetSigDB source
was confirmed2026-10-01T12:03:12.509258+00:00; original total16956.964756 s,
including6682.071756 s of uncompleted subset work. Base counts/indexes and
shared derivations were retained; only the subset transaction rolled back.

**Outcome**: base and shared derivations committed; subset-only retry is recorded
in attempt9 above. Fresh private prelaunch at
2026-10-01T07:20:09.881190+00:00 verified exact clean checkout/pypath, WAL8192MB
without pending restart, full schema absent/prior two schemas preserved, no
previous import/verifier backends, data ready and529.051GiB free. Docker identity,
start timestamp and running state matched the baseline for all four protected
services. This is a prelaunch observation; final post-build checks remain pending.
Committed checkpoint counts are unchanged:6383473 entities,80136865 identifiers,
18339902 relations,22887251 evidence,150690040 annotations. No full-release
performance conclusion yet; full comparison and PG-only rebuild/rollback are
queued after build success.

## 2026-10-01 — deferred-constraint local fixture validation (completed)

**Reason**: confirm that post-COPY key/FK construction gives the same valid final
schema and rejects invalid data before restarting the private full import.
No downloads, resolver runs or real resource builds.

**Parameters**: disposable local PostgreSQL16; two synthetic resources, each
20 published raw records; DuckDB512MB/one thread, raw audit false,
checkpoint_base true; three immediate and three deferred loads. Frozen Python
helper `/private/tmp/omnipath-late-constraints-benchmark.py`; aggregate report
`/private/tmp/omnipath-late-constraints-benchmark.json`. Unit tests separately
exercise audit on/off, checkpoint on/off, empty/multiple resources and quoting.

**Phase durations**: six fixture wall times 0.291240, 0.073713, 0.067253,
0.068565, 0.066102 and 0.065055 s. Deferred key/statistics/FK phases were
0.003733, 0.003721 and 0.003710 s. Each retained 6 entities,4 identifiers,
2 relations,40 evidence and90 annotations. Detailed phase placement is in the
local aggregate report; potential overlap with regression tests and disabled
fixture fsync mean these times do not establish comparative or durable throughput.

**Outcome**: all fixture counts matched; disposable files/cluster removed.
136 focused tests passed in6.16s and the broader PostgreSQL/subset suite passed
**527 tests in28.39s** (one existing Starlette deprecation warning).
Real-PG catalog comparisons include every PK/FK/check/default/index definition,
NOT NULL/sequence, FK actions/deferrability/validation and annotation uniqueness.
Negative cases cover all five PKs/all seven FK categories; checkpoint retry works
without files/DuckDB/COPY. Ruff/format, whitespace checks and independent review
passed. Full-release performance remains unmeasured for this change.

## 2026-10-01 — full resolved-Parquet release on nicesrv, attempt 7 (cancelled)

**Reason**: user-authorized restart after attempt 6's confirmed rollback, using
frontier ontology traversal and a persistent validated base/index checkpoint.
Once this checkpoint commits, late derivation failures can use PostgreSQL-only
`finish_release` / `--finish` without repeating COPY. No resource rebuilds,
reference preparation, resolution or version changes.

**Parameters**: exact checkout `552e1d0efdf8e5e969080fada2ff17bc923b6c21`,
checkpoint loader `552e1d0`, ontology frontier `4fe91d1`, replay `300d06b`,
pypath `e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`. Same 46 pins, release
`2026.9.30.4`, schema `full_20260930`, canonical release-manifest SHA-256
`06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`.
Invocation: `uv run --frozen --no-sync --directory
/root/projects/omnipath-migration/omnipath-build python
/root/projects/omnipath-migration/full-release-20260930/postgres-sql-build.py`.
`load_release(..., batch_size=1024, duckdb_threads=1, memory_limit='1024MB',
validate_source_records=False, checkpoint_base=True, temp_directory=state)`
then all three products. CSV16MB rotation / COPY1 MiB byte blocks; work_mem16MB,
maintenance_work_mem128MB, max_parallel_workers_per_gather1, JIT off.
Importer/verifier each MemoryMax3GiB, swap256MiB, CPUQuota200%, Nice10,
stop timeout90s; existing private PostgreSQL16 stays at4GiB/CPU2. Disk floor64GiB.
Import unit `omnipath-migration-full-postgres-sql-attempt7-20260930.service`,
invocation `526ba7385eb4456794247340583c4c66`; verifier unit
`omnipath-migration-full-verification-sql-attempt7-20260930.service`, invocation
`945a61a9ce784be89b7efb15358a31f4`. Coordinator code SHA-256
`ce6c28b559ab1e834f5efb140f34dc0cde012a189145141abfd85ef0322b61f0`.

Base commits only after mandatory file recheck, normalized integrity and indexes;
a session advisory lock spans commits. Derivation publishes `release_metadata`
in its own transaction. A persistent marker captures exact pins/hashes, audit
metadata, counts and base index definitions. Default loads remain atomic.
Network registration/streaming and COSMOS reject unpublished schemas before
mutation/serving. All486 PostgreSQL/subset tests passed in31.81s, Ruff/format
checks passed, and independent review passed after fixing quoted mixed-case
checkpoint lookup. Resume tests remove input files and forbid DuckDB/COPY.

**Phase durations**: started **2026-10-01T06:20:06.778982+00:00** and
cancelled after **1936.394977 s** (32m16s), before the base checkpoint. Release
read completed in 6.325272 s. Partial COPY counters were 2,808,540 entities;
25,497,190 identifiers; 2,172,554 relations; 3,279,000 evidence; 27,245,121
annotations, totaling 61,002,405 rows. These are historical uncommitted counters;
no derivation or verification completed. Completed-table CSV staging totaled
245.901000 s and COPY 1566.713974 s; typed validation, including ChEMBL's completed
validation, totaled 46.993674 s. Six resources completed all five COPY tables.
Their staging/COPY seconds were BindingDB 231.736932/1400.387034 (largest),
BRENDA 1.120558/15.785462, CellChat 0.862258/8.024571,
Cellinker 1.035992/12.741105, CellPhoneDB 0.109843/1.041236 and
ChEBI 11.035417/128.734566. ChEMBL staging was interrupted before a table result
was recorded, so its partial work is excluded from these timing totals.

**Outcome**: user-authorized cancellation to move five base primary keys,
annotation-owner uniqueness and all seven foreign keys after COPY, with complete
validation before the base checkpoint. Both private attempt-7 units stopped.
Fresh `attempt7-rollback.json` at 2026-10-01T06:54:51.539126+00:00 confirms the
full schema absent, both pilot schemas preserved, no importer/verifier backends,
private data present and hidden data absent. Reports/logs/coordinator were archived
server-locally under unique tag `20261001T070135325182Z-e0503e02`.
All 46 resource artifacts, cache, reference and version pins remain reusable.

A separate private-only configuration change at
2026-10-01T06:56:28.286142+00:00 increased `max_wal_size` from 1024MB to 8192MB
via `ALTER SYSTEM` and reload. The applied value has no pending restart;
`fsync`, `synchronous_commit` and `wal_level` were unchanged. No container restart
or production database change was made. Deferred-constraint code and tests are
being prepared for attempt 8; no full-scale speedup has yet been measured.

## 2026-10-01 — PostgreSQL ontology frontier microbenchmark (completed)

**Reason**: replace cumulative closure expansion with breadth-first frontiers,
then confirm unchanged results and measure the contained algorithmic improvement.
No DuckDB, schema, resolver, resource or running-server changes.

**Parameters**: local disposable PostgreSQL 16 cluster, Unix socket only,
shared_buffers 32MB, work_mem 16MB, JIT disabled, statement timeout 60s.
Direct `_closure` calls compare baseline revision `4a643f2` with the frontier
implementation on one synthetic 120-node chain: first 60 links in resource
`first`, remaining 59 in `second`, version `v1`, kind `subclass`, six physical
copies per link (714 rows). Timed calls include adjacency/frontier setup,
traversal, statistics and cleanup; equality is checked outside the timer.
Invocation: `OMNIPATH_TEST_POSTGRES=1 uv run --frozen --no-sync python
/private/tmp/omnipath-frontier-benchmark.py`. This is a function benchmark,
not a resource build. The recorded run occurred after the regression suite
finished; an earlier concurrent run was used only as an initial smoke check.

**Phase durations**: scoped closure: baseline 0.555150 s, frontier 0.049972 s,
3,600 identical pairs (about 11x faster). Global closure: baseline 1.437431 s,
frontier 0.083676 s, 7,140 identical pairs (about 17x faster).

**Outcome**: exact pair/depth equality in both scopes. Seven focused integration
tests passed in 2.60 s; all 424 PostgreSQL/subset tests passed in 27.26 s;
Ruff lint/format, whitespace checks and independent review passed.
The tests include independent Floyd–Warshall equality, cycles, shortcuts,
depth 25, kinds, resource/version isolation, cross-resource paths, duplicate
provenance, repeated/empty rebuilds, search-path decoys and interruption rollback.
One 14-link planner regression verifies frontiers of 14 through 1 pairs:
105 total expansions, with current temporary-table statistics and usable indexes.
This single synthetic graph does not measure full-release time or peak memory.
One deduplicated indexed adjacency table and two frontiers add working storage;
final closure size is unchanged. Attempt 6 continues on its original revision.

## 2026-10-01 — full resolved-Parquet release on nicesrv, attempt 6 (cancelled)

**Reason**: retry the unchanged release after attempt 5's clean rollback, with
ontology traversal indexes and sampled statistics available before closure
expansion. Cancellation followed code review of late index creation and missing
fresh edge/closure statistics. The optional actual-release probe did not obtain
EXPLAIN, so no confirmed stall, production plan or full-scale speedup is claimed.

**Parameters**: exact checkout `4a643f2750bf1ae35e54700c394ad5cb0ae9a054`;
loader `420c407`, ontology `4a643f2`, native replay `300d06b`, pypath
`e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`. Unchanged release `2026.9.30.4`,
all 46 pins and schema `full_20260930`; canonical release-manifest SHA-256
`06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`.
Invocation: `uv run --frozen --no-sync python
/root/projects/omnipath-migration/full-release-20260930/postgres-sql-build.py`.
`load_release(..., batch_size=1024, duckdb_threads=1, memory_limit='1024MB',
validate_source_records=False, temp_directory=state)` then all three products.
CSV/COPY settings, session settings and resource bounds remain those of attempt 5.
The coordinator adds only source-scoped/release-wide ontology timing wrappers;
SHA-256 `f34108ae832110cdfa8946517abe98ee424789bf15a12ffd6962ead28902c554`.
Import unit `omnipath-migration-full-postgres-sql-attempt6-20260930.service`,
invocation `25bfa2dc65804317a4774393b6e98dc0`; queued verifier unit
`omnipath-migration-full-verification-sql-attempt6-20260930.service`, invocation
`76eb4fdaa2764ebba962a0264138e55c`.

Existing edge and ancestor indexes now precede both closures; populated edge,
seeded closure and changed closure statistics are refreshed. Schemas, index
expressions, shortest-pair and cycle semantics remain unchanged. All 420
PostgreSQL/subset tests passed in 23.40 s and independent review passed. Mandatory
artifact, normalized-integrity and reaction-annotation checks remain enabled.

**Phase durations**: started 2026-10-01T01:07:23.133303Z. A single read-only
aggregate snapshot at 2026-10-01T03:44:07.208111Z records all 46 resources / 230
tables copied: typed validation 126.211743 s, CSV staging 1015.314059 s and COPY
8138.151367 s. Completed final phases were release read 6.537909 s and schema
analyze 4.268968 s. A later aggregate snapshot at
2026-10-01T05:17:39.873059Z confirms loaded-data validation completed in
1441.182233 s and indexes in 840.957112 s. Source-scoped ontology closure began
at elapsed 12250.685 s and remains active in the 05:38 UTC monitor snapshot.
The largest source COPY totals were ChEMBL
2172.242730 s (staging 354.214143 s), BindingDB 1424.617466 s (staging 228.086378 s)
and FooDB 793.982437 s (staging 94.331009 s). Scoped/global ontology,
other derivation, product and verification timings remain pending.

**Outcome**: COPY complete for 278,437,531 rows: 6,383,473 entities; 80,136,865
identifiers; 18,339,902 relations; 22,887,251 evidence; 150,690,040 annotations.
CSV output totaled 7,807 chunks / 139,158,600,874 bytes. These are completed
projection/COPY metrics, not committed or independently verified final PostgreSQL
counts. Loaded-data validation and index creation passed; commit/rollback,
derivations, products and full verification remain pending. No source builds or
version changes.

**Final outcome**: user-authorized cancellation to restart with frontier traversal
and a persistent base checkpoint. Attempt 6 ran for **17745.644952 s**; its
source-scoped closure completed in **5340.111027 s**, then global closure began
at elapsed 17590.798 s and was interrupted before completion. No base transaction
committed. Fresh read-only `attempt6-rollback.json` at
2026-10-01T06:05:52.096952+00:00 confirms `full_20260930` absent, both prior pilot
schemas preserved, no importer/verifier backends, private data present and hidden
data absent. Stopped only the two private attempt-6 units. Reports/logs/coordinator
were archived under unique tag `20261001T060641395642Z-0cb8cf39`; free disk was
529.052 GiB after rollback. All resource artifacts and pins remain reusable.
The earlier 05:38 snapshot above is historical; attempt 6 is no longer running.

## 2026-10-01 — bounded actual-release ontology plan probe on nicesrv (incomplete)

**Reason**: inspect fresh-table plans before and after ontology traversal index/statistics
ordering change `4a643f2`, while attempt 5 was running unchanged. This diagnostic did
not reload entities or resources and never connected to the import or production database.

**Parameters**: private release `2026.9.30.4`, all 46 exact pins; canonical manifest
SHA-256 `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`.
Read-only frozen runtime `/root/projects/omnipath-migration/omnipath-build`, HEAD
`c27a37f03d95d99d19f19f468f389d06234cc172` / loader `420c407`.
`uv run --no-cache --frozen --no-sync --directory <frozen-runtime> python
<exclusive-private-directory>/probe.py --workdir=<exclusive-private-directory>
--unit=<exclusive-probe-unit>.service --container=<exclusive-probe-container>`;
the final classification used `--phase=bootstrap` and read no Parquets.
Production DuckDB relation projection, filtered to ontology statements, plus the
production Biolink direction/hierarchy-kind rules; DuckDB 1.5.6, one thread,
512MB working memory, helper cgroup at most 1 GiB / 50% CPU. PostgreSQL was a new
no-network, no-published-port container from existing image
`sha256:52f87503b2fe5ad2e0d209b1a55830db42eb3ecd3eaf1b583accb0091ab30bff`.
PG/controller limits were 512/384 MiB and 40%/10% CPU; shared buffers32MB,
work_mem8MB, maintenance32MB, gather0, jit off. The graph cap was 3,000,000 edges;
combined CSV/PG-data/temp cap2 GiB, disk floor64 GiB, operation limit600 s,
readiness limit60 s, preparation statement limit90 s and candidate limit10 s.
Only one expansion round per case was permitted. The minimal bootstrap also had
an outer systemd limit90 s. No image/package downloads, checkout updates or raw-payload reads.

The intended baseline was empty-edge `ANALYZE`, edge population, then fresh seeded
source/release closures without traversal indexes or refreshed statistics. The new
case would create the unchanged traversal indexes before expansion and analyze
populated edges and each seeded/changed closure. Neither case reached edge COPY.

**Phase durations**: initial wrong-runtime import failed safely in 0.088385 s.
Three exact edge extractions took 6.233286, 6.171601 and 6.121330 s; each produced
987,393 edges from 12 resources: 972,599 subclass and 14,794 part-of edges.
All three CSVs were 361,977,795 bytes with SHA-256
`0d722a8e7aca0d5375c1b5b1df306fe2eb25af7b438bf14e400f610aaa00594d`.
Their isolated bootstrap attempts ended at 29.393892, 69.794769 and 69.691341 s,
including extraction and cleanup. The final no-Parquet bootstrap check took
62.703720 s. Total elapsed across attempts was 231.672107 s. No source build,
entity load, successful edge COPY, closure EXPLAIN or closure timing occurred.

**Outcome**: incomplete; stopped after the minimal bootstrap failed. The helper
corrected an image UID70/socket-owner mismatch and a 110-byte host socket address,
using a held 29-byte descriptor alias. With explicit server/client port5432,
correct PGDATA bind and the expected `.s.PGSQL.5432` present, libpq recorded six
initial missing-file errors followed by 234 connection-refused errors; the
container log had one readiness marker and no permission/FATAL/OOM markers.
OS errno was not exposed separately. This does not establish a closure performance
result or speedup. All created probe containers and private work directories were
removed. Both mounts were binds; no anonymous volumes were created, and final
cleanup also requested removal of owned anonymous volumes. Sampled peak private
files were 402,056,150 bytes; minimum free disk was 359,774,490,624 bytes.
Attempt 5, full/pilot schemas, protected services, source artifacts and checkouts
were untouched. The secret-free aggregate report is retained locally at
`/private/tmp/omnipath-ontology-plan-probe-summary.json`; per-attempt reports remain
in the exclusive server directories, ending with
`/root/projects/omnipath-migration/full-release-20260930/ontology-plan-probe-Oc6omQEk/report.json`.
Aggregate upload was blocked by automatic approval review as private-metadata egress.

## 2026-09-30 — full resolved-Parquet release on nicesrv, attempt 5 (cancelled, rolled back)

**Reason**: retry the fixed release after HMDB entity staging exceeded attempt
4's 512MB DuckDB allowance. All 46 resources / 230 exact production rotated-CSV
projections have now passed with 1024MB and one thread. The preflight staged all
selected records, checked each table's expected count, and removed its temporary
CSV/spill directories. It took 1144.172568 s, process peak RSS 803,992 KiB.
It made no PostgreSQL writes. The separate HMDB entity test staged all 220,400
rows in 8.489402 s and used process peak RSS 768,216 KiB.

**Parameters**: loader `420c407` (DuckDB projection `91317a3`, narrow statistics
`15500f7`, bounded alias index `3c9bab3`); native replay remains `300d06b`,
pypath `e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`.
Unchanged release `2026.9.30.4`, canonical manifest SHA-256
`06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`;
schema `full_20260930`, same 46 exact pins. Invocation:
`uv run --frozen --no-sync python
/root/projects/omnipath-migration/full-release-20260930/postgres-sql-build.py`.
`load_release(..., batch_size=1024, duckdb_threads=1, memory_limit='1024MB',
validate_source_records=False, temp_directory=state)` then all three products.
Database/role `omnipath_migration`, localhost5440, separate PG16 container.
Process cap 3 GiB RAM / 256 MiB swap / two CPUs; PG4 GiB/two CPUs; disk floor
64 GiB. CSV16MB rotation threshold / COPY1 MiB byte blocks. Session work_mem16MB,
maintenance128MB, gather1, jit off. Unit:
`omnipath-migration-full-postgres-sql-attempt5-20260930.service`.

Raw source auditing is now disabled by the approved default; `--validate-source-records`
retains the old deeper audit when requested. File schemas/counts/byte sizes/SHA
before and after loading, normalized integrity, owner constraints, reaction
annotation checks and atomic rollback remain mandatory. All 419 PostgreSQL/subset
tests passed in 23.56 s, including identical stored output in both modes.

**Phase durations**: started 2026-09-30T19:22:54.217010Z at checkout
`c27a37f03d95d99d19f19f468f389d06234cc172`; total 19,498.681013 s before cancellation
and rollback. All 46 resources / 230 tables completed: typed validation 129.240034 s,
CSV staging 1078.313069 s and COPY 8485.646031 s; 7,807 chunks /
139,158,600,874 CSV bytes. Completed final phases: release read 5.939686 s,
schema analyze 3.841620 s, loaded-data validation 1455.021736 s and indexes
832.696578 s. Derivation was unfinished; no subset or full-verification result.

**Outcome**: cancelled and rolled back after all 278,437,531 expected rows copied:
6,383,473 entities; 80,136,865 identifiers; 18,339,902 relations; 22,887,251 evidence;
150,690,040 annotations. Cancellation followed the reviewed late ontology traversal
index/statistics issue. The optional production-edge probe did not reach EXPLAIN;
a stall or production query plan was not confirmed, and no speedup was measured.
At 2026-10-01T00:54:31.617812Z, `attempt5-rollback.json` confirmed `full_20260930`
absent, both old schemas (`pilot_20260930`, `annotation_sample_20260930`) preserved,
no import/verifier backends, private data restored and no hidden-for-verification
directory. No base load committed. Unique archives end
`20261001T005431618899Z-c9f49013`: `postgres-build-attempt5-archive-<suffix>.json` /
`.log` and `verification-wait-attempt4-attempt5-archive-<suffix>.log`.
No resource rebuilds or source-version changes occurred. These checkpoint
observations precede attempt 6.

## 2026-09-30 — full resolved-Parquet release on nicesrv, attempt 4 (failed)

**Reason**: retry the unchanged full release after supplying narrow planner
statistics before every source-owner check. The third attempt copied BindingDB
correctly, then was cancelled during inefficient owner lookups and rolled back.
The new regression shows that unanalysed resource/version index probes filter
98,976 of 100,000 synthetic rows per 1,024-owner batch. After narrow ANALYZE,
both custom and generic plans put all three scoped keys in the index condition.
Cache warming affects timings; predicate placement is the primary evidence.

**Parameters**: loader `3c9bab3` (owner statistics `15500f7`); native replay remains `300d06b`, pypath
`e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`. Unchanged private release
`2026.9.30.4`, manifest SHA-256
`06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`;
schema `full_20260930`; all 46 exact pins in `release.json`.
Invocation: `uv run --frozen --no-sync python
/root/projects/omnipath-migration/full-release-20260930/postgres-sql-build.py`.
Same database/role `omnipath_migration`, localhost 5440, isolated PG16 container.
Same one-thread / 512 MB DuckDB SQL → 16 MB CSV rotation threshold → 1 MiB
PostgreSQL COPY byte blocks; raw batches 1024 dictionary indices / 64 plain rows.
Process limits: 3 GiB RAM / 256 MiB swap / two CPUs; PG: 4 GiB/two CPUs;
64 GiB disk reserve. Session work_mem16MB, maintenance128MB, parallel gather1,
jit off. Unit: `omnipath-migration-full-postgres-sql-attempt4-20260930.service`.

Each resource now gets scalar resource/version/key statistics for entities and
relations before owner queries. Reaction resources sample their additional
scalar join/filter columns in the same call. No preliminary wide-JSON sampling,
query rewrite or change to integrity rules. Base, indexes, derivations and exact
release metadata still share one transaction; all three products follow after
base commit. A queued verifier then checks the complete release read-only and
rebuilds products with private Parquets unavailable, rolling that rebuild back.

**Phase durations**: 7864.979847 s total; release read 5.886399 s; typed SQL
validation 76.155310 s; successful CSV staging 702.586871 s; PostgreSQL COPY
4773.493333 s; raw audits 2286.706690 s for 12,431,857 occurrences; narrow owner
statistics 13.223463 s. Fourteen resources completed all five tables and audits.
The copy total was 161,266,864 rows in 4,726 CSV files / 84,911,231,452 bytes.
ChEMBL's raw audit took 845.336509 s and FooDB's took 925.231693 s. All completed
tables matched their expected counts. No final indexes, derived tables, products
or full verification ran. The changed loader passed all 396 PostgreSQL/subset
tests in 21.86 s; Ruff, diff checks and independent review passed. The separate
39-resource raw validation remains successful.

A read-only screen covered all 46 selected resources: canonical identifiers and
labels had maximum 1240 UTF-8 bytes; 2442 alias index candidates exceeded a
2000-byte composite-field screen. The unbounded alias index failed in a local
disposable PG16 cluster on these actual values (3704-byte entry versus 2704-byte
limit). The index now uses the first 256 lowered characters with text_pattern_ops;
all stored TEXT and nested/view values remain complete. The same 2442 values
passed actual index creation in 0.013856 s and exact full-value comparison.
Candidate-file SHA-256: d510c11aee16ac89e3170cfccac08e941d910e58f2b843ad7ad63e702d84827b.
Shared-prefix and four-byte Unicode regressions preserve complete lookup results.
This preflight prevented a late full-transaction index failure; no source rebuild
or query API change was needed.

**Outcome**: failed while staging HMDB entities, before any HMDB table copied.
HMDB SQL validation passed in 8.789128 s; its CSV projection then needed a
128 MiB allocation with 394.0 MiB already used against the 488.2 MiB effective
DuckDB allowance. `OutOfMemoryException` unwound the outer transaction.
At 2026-09-30T18:45:03.288337Z, a fresh private connection confirmed
`full_20260930` absent, no import/verifier backends, and both `pilot_20260930`
and `annotation_sample_20260930` preserved. Thus no base load committed.
Archived private evidence: `postgres-build-attempt4.json`/`.log`,
`attempt4-rollback.json`, and `verification-wait-attempt4.log`.

The seven uncapped native replays still cover all 300,008 originally published
source records; 39 other pins are reused. No Parquets or production services
were changed by this failed import. Before retrying, test the exact projections
at a larger bounded DuckDB limit. Future loads use `420c407`'s approved default
`validate_source_records=False`; optional raw-audit counts will be empty while
artifact checks, normalized integrity and reaction-annotation checks stay
mandatory. This failed attempt itself retained the original deep audits.

## 2026-09-30 — full resolved-Parquet release on nicesrv, attempt 3 (cancelled)

**Reason**: run the complete fixed release through the restored DuckDB SQL
projection and bulk COPY path after exact equivalence tests and a real-data
comparison. The two Python-loading attempts rolled back without base commits.

**Parameters**: loader `91317a3`; unchanged native replay `300d06b`, pypath
`e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`. Private release `2026.9.30.4`,
manifest SHA-256 `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`,
schema `full_20260930`; exact 46 pins in the private `release.json`.
Invocation: `uv run --frozen --no-sync python
/root/projects/omnipath-migration/full-release-20260930/postgres-sql-build.py`.
The script calls `load_release(..., batch_size=1024, duckdb_threads=1,
memory_limit='512MB', temp_directory=state)` then all three subset adapters.
Database/role `omnipath_migration`, localhost 5440, isolated PG16 container
(two CPUs / 4 GiB RAM). Process: two CPUs / 3 GiB RAM / 256 MiB swap;
64 GiB free-disk reserve. Session: work_mem16MB, maintenance128MB,
max_parallel_workers_per_gather1, jit off. Unit:
`omnipath-migration-full-postgres-sql-20260930.service`.

DuckDB validates typed source shapes/counts in two aggregate scans per resource,
then flattens all five tables and serializes JSON in SQL. One table at a time
is written to CSV chunks (16 MB rotation threshold), forwarded as 1 MiB byte
blocks to PostgreSQL COPY, and deleted. Each chunk is a separate immediate-FK
statement inside the atomic base transaction; complete parent tables precede
children. No Python projected-row expansion. Raw validation remains the separate
dictionary-preserving reader with every pointer/owner/SHA/type/count check.
No source parsing, entity resolution, downloads or resource rebuilds.

**Phase durations**: deliberately cancelled after 2002.583925 s. Release
reading took 5.378590 s; BindingDB SQL shape/count validation took 16.292700 s.
All five BindingDB tables copied: 52,774,666 rows. Their SQL staging total was
238.644504 s and PostgreSQL COPY total 1447.194627 s. Source-owner validation
started at elapsed 1708.965 s but did not finish; indexes, derived tables and
products did not run. A live owner-key query waited on DataFileRead after
36.416660 s. The separate 39-resource raw validation remains successful.

**Outcome**: cancelled to correct missing planner statistics before non-reaction
source-owner checks. The preliminary scalar ANALYZE had covered only reaction
resources; the batch owner lookups also need current entity/relation statistics.
No data or schema change is required. Psycopg cancelled the active query and
rolled back. A fresh private connection confirmed `full_20260930` absent, no
import/verifier backend remaining and both previous private schemas present.
`postgres-build-attempt3.json`/`.log` and `attempt3-rollback.json` retain the
results. The queued verifier stopped before SQL checks. The seven full native
rebuilds, 39 reused versions and fixed release pins remain unchanged. Production
was not published or modified.

## 2026-09-30 — DuckDB bulk projection/COPY comparison on nicesrv (completed)

**Reason**: confirm the restored SQL projection/COPY path is faithful and measure
its cost against the Python implementation before restarting the full release.

**Parameters**: code `91317a3`, private `bulk-benchmark.py` in
`/root/projects/omnipath-migration/full-release-20260930`, invoked with
`uv run --frozen --no-sync python`. The first 50,000 published BindingDB entities
(version `2026.9.8.1`) are written once to a temporary Parquet subset; relations
are empty in both paths. SQL and Python consume identical files. Both load all
five target tables with complete entity JSON, identifier and annotation
occurrences and the same PostgreSQL schema/primary keys/immediate foreign keys.
SQL uses one DuckDB thread, 512 MB working memory, 16 MB CSV rotation threshold
and 1 MiB byte forwarding. Python uses the previous reference projector,
1024-row COPY batches and parent flushing. All writes use one private transaction
on database/role `omnipath_migration`, localhost 5440, with separate new schemas
`bulk_sql_benchmark_20260930` and `bulk_python_benchmark_20260930`; the transaction
is rolled back after exact bidirectional EXCEPT ALL checks on every projected
column. No source parser, resolution, download or product/index derivation.
Process cap: 3 GiB RAM / 256 MiB swap / two CPUs; isolated PG cap: 4 GiB/two CPUs.

**Phase durations**: subset preparation 0.869144 s; SQL projection/COPY
17.944366 s, Python projection/COPY 28.587794 s. SQL staging took 0.997005 s
across five tables; PostgreSQL COPY took 16.913948 s. Total benchmark/check/
rollback script 55.885584 s. SQL ran first once; Python benefited from any warmed
cache. Process peak RSS 532,532 KiB is cumulative and cannot compare the methods'
memory independently. The cgroup peak was 467.8 MiB. This is a projection/
transport comparison, not a full-release performance estimate.

**Outcome**: success. Both paths loaded the same 50,000 entities and 880,917
identifiers. The relation, evidence and annotation outputs were empty; those
populated cases are covered by the integration fixtures. All projected columns,
including full entity JSON, passed exact bidirectional EXCEPT ALL. Both schemas
were rolled back and confirmed absent. SQL reduced measured elapsed time by
about 37%; no full-release speed claim follows from this subset. The corrected
loader passed 389 PostgreSQL/subset tests
(19.57 s) plus a separate actual multi-file CSV rotation regression (1.63 s).
Independent code/benchmark review, Ruff and diff checks passed.

## 2026-09-30 — full resolved-Parquet release on nicesrv, attempt 2 (cancelled)

**Reason**: retry the complete fixed release after correcting repeated raw payload
expansion. The first attempt rolled back before committing; no artifacts changed.

**Parameters**: loader `a71fcf7`, same release `2026.9.30.4`, canonical manifest
SHA-256 `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`,
schema `full_20260930`, `batch_size=1024`, database/role `omnipath_migration` on
localhost 5440. Same `uv run --frozen --no-sync python .../postgres-build.py`,
46 resource pins, 3 GiB process/4 GiB container memory limits, two CPUs,
256 MiB process swap, 64 GiB disk reserve and session settings as attempt 1.
The unit is `omnipath-migration-full-postgres-attempt2-20260930.service`.
Raw reading retains one native dictionary and at most 256 index validation
results; dictionary batches have at most 1024 scalar owners. Plain raw and nested
record batches remain capped at 64. Every occurrence still receives pointer/key,
owner and count checks; SHA/type validation reuses only exact dictionary values.
No raw body is stored in PostgreSQL. All three products are built after base commit.

**Phase durations**: cancelled at the user's efficiency concern after
1040.683043 s; release reading 5.792572 s, completed COPY calls 699.742567 s.
28,084,290 copied rows are enumerated in `postgres-build-attempt2.json`; indexes,
derived tables and subsets did not run. Before restart, the actual 42,536-row
FooDB block passed in 1.010 s with peak RSS 472,936 KiB and one JSON parse/SHA
call for its 42,831,695-byte text. All 297 disposable PostgreSQL/subset tests
passed (14.22 s); independent review, Ruff and diff checks passed.

The separate artifact-only check completed all 39 reused resources, covering
23,561,803 raw rows, in 284.930287 s (completed first-attempt checks reused;
inputs and hashes unchanged). FooDB's 4,211,723 occurrences took 48.381147 s.
The corrected check peaked at 594 MiB inside its 3 GiB/one CPU cap.

**Outcome**: deliberately cancelled to restore the intended DuckDB SQL projection
and bulk COPY pipeline. The Python row-expansion implementation was inefficient
at full scale. Psycopg cancelled active COPY; a fresh private connection confirmed
`full_20260930` absent, no import backend remaining and both earlier private
schemas retained. The queued verifier was stopped before it began SQL checks.
`attempt2-rollback.json` records the rollback. The full artifacts and successful
raw validation remain available; native replay scope is all 300,008 originally
published records across seven resources, with 39 unchanged pins. Nothing was
published to production.

## 2026-09-30 — full resolved-Parquet release on nicesrv, attempt 1 (cancelled)

**Reason**: validate the complete release after moving all former PostgreSQL
payload consumers to published annotations. The seven affected resources were
rebuilt without a record cap; 39 unchanged resources retain their original pins.

**Parameters**: loader code `c2dfc11`, native replay code `300d06b`, pypath
`e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`. Private release `2026.9.30.4`,
canonical manifest SHA-256
`06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`.
Exact resource pins, file checksums and preflight counts are retained in
`/root/projects/omnipath-migration/full-release-20260930/release.json` and
`final-preflight.json`. Data root: that directory's `data/`.

The private `launch-replay.py` ran `scripts/replay_resource.py` against the seven
original published resource directories, with version `2026.9.30.4`, immutable
reference `/root/projects/full_parquet/data/reference/.library-compact-20260912`,
private cached Human-GEM TSV, `--resource-ram-gb 3 --batch-workers 1
--duckdb-threads 1 --memory-limit 512MB --min-free-disk-gb 64`. No `--max-records`
argument. The process group had a 6 GiB memory / 512 MiB swap cap and two CPUs.
All 300,008 published source records, original row IDs and exact raw text passed
replay verification. Unpublished rows discarded by historical builds are outside
this replay's coverage. The 39 reused versions are independent file copies.

Database invocation, from the isolated migration checkout:

```sh
uv run --frozen --no-sync python \
  /root/projects/omnipath-migration/full-release-20260930/postgres-build.py
```

The script calls `load_release(..., schema='full_20260930', batch_size=1024)`
then `build_subsets(..., products=('metsigdb','network_views','cosmos'))` in a
separate transaction after the base release commits. Database/role
`omnipath_migration`, localhost port 5440, isolated container
`omnipath-migration-postgres` (2 CPUs / 4 GiB RAM). Import process: 3 GiB memory,
256 MiB swap, two CPUs, minimum disk reserve 64 GiB on the same filesystem as its
Docker volume. Session: `work_mem=16MB`, `maintenance_work_mem=128MB`,
`max_parallel_workers_per_gather=1`, `jit=off`. Buffered Arrow decoding uses at
most 64 rows; exact-text validation cache at most 256 entries / 64 MiB charged
storage. COPY foreign keys are immediate within the one atomic base transaction.
No source deletion, download, reference preparation, reparse or entity resolution
runs in PostgreSQL. Existing schemas and production services are preserved.

**Phase durations**: full native replay Rhea 329.122207 s, Recon3D 62.496945 s,
MetAtlas 76.459144 s, KEGG 388.689825 s, Reactome 101.476603 s,
MACdb 34.849949 s, ConnectomeDB2025 61.983096 s. Per-dataset coverage, native
prepare/finalize timings and source hashes are in `rebuilt-data/replay_report.json`.
KEGG is the replay outlier. Final fixed-release preflight 10.144114 s.
The first database attempt ran 2527.816089 s before deliberate cancellation.
Release reading took 5.889342 s; completed COPY calls took 1611.889483 s.
BindingDB projection was incomplete; indexes, derived tables and products did
not run. `postgres-build-attempt1.json` and `postgres-build-attempt1.log`
retain its phase/source progress. Source projection/COPY and raw-validation intervals
are approximate (final partial batches and preliminary ANALYZE sit outside those
intervals); loader totals and separate COPY call totals are authoritative.

**Outcome**: cancelled before any base commit after 51,925,536 rows were copied
inside the transaction. A concurrent artifact check identified a 42,831,695-byte
FooDB dictionary value repeated 42,536 times; expanding 64 rows materialized
2,741,228,480 bytes before Python validation. This would make later loading
unsafe/slow. Only the private import and checker units were stopped to replace
the raw reader with dictionary-preserving validation. Psycopg cancelled COPY and
rolled back; a separate connection confirmed `full_20260930` absent, no import
backend remaining, and both existing private schemas retained. Production
container IDs/start times were unchanged. Disk returned from about 491 GiB free
to 530 GiB. No source artifacts or pins changed.

An earlier artifact-only check was OOM-killed while DuckDB materialized FooDB;
it made no database writes. Two bounded Arrow attempts were deliberately stopped
(the second to enable exact-text caching, the third for dictionary preservation).
Completed checks from the first attempt are reusable because input hashes are
unchanged. Retry and full validation remain pending. Preflight expects 6,383,473
entities, 80,136,865 identifiers, 18,339,902 relations, 22,887,251 evidence
occurrences and 150,690,040 annotations; 26,711,866 raw payload rows are validated
and discarded. A failed/interrupted commit boundary must be inspected through
release metadata rather than inferred from a missing completion flag.

## 2026-09-30 — capped annotation migration on nicesrv (completed)

**Reason**: verify that all former PostgreSQL product payload readers can use
published annotations, so raw source bodies can stay exclusively in Parquet.

**Parameters**: build code `7ed2e63`, pypath `e27d72e09`; seven new private
resource versions and private release `2026.9.30.3`. Exactly 20 original records
per resource, 140 total. Sources: Rhea, Recon3D, MetAtlas/Human-GEM, KEGG,
Reactome, MACdb and ConnectomeDB2025. Selected originals and input hashes live in
`/root/projects/omnipath-migration/annotation-migration-20260930/inputs`.

Native rebuild invocation, from the isolated migration checkout:

```sh
uv run --frozen --no-sync python scripts/annotation_migration_smoke.py \
  --inputs /root/projects/omnipath-migration/annotation-migration-20260930/inputs \
  --output-dir /root/projects/omnipath-migration/annotation-migration-20260930/data \
  --library-dir /root/projects/full_parquet/data/reference/.library-compact-20260912 \
  --version 2026.9.30.3
```

One preparation worker and DuckDB thread, 1 GiB resource RAM, batches at most
20 records, 20 relations and 64 MiB. The selected iterator enforces a total
20-record cap across datasets, beyond the pipeline's per-dataset cap. No source
downloads or reference preparation. Native resolution used the existing immutable
prototype reference. An initial launch with an incorrect `data/ref` path failed
before creating output or consuming records; the verified `data/reference` path
above succeeded.

Database invocation is the private `postgres-check.py` in that state directory,
using `uv run --frozen --no-sync python` from the isolated checkout. Loader
`batch_size=64`, DuckDB one thread / 256 MiB, database/role `omnipath_migration`
on localhost port 5440 in `omnipath-migration-postgres`. Fresh schema
`annotation_sample_20260930`; existing pilot schema retained. No source deletion,
preparse or canonicalization runs in the loader. Source bodies are streamed and
validated, including exact SHA/type checks, then discarded.

**Phase durations**: native resource builds: Rhea 5.66 s, Recon3D 3.20 s,
MetAtlas 4.16 s, KEGG 4.63 s, Reactome 4.19 s, MACdb 3.01 s,
ConnectomeDB2025 3.21 s. Database code `4c896eb`; validate_files 0.011331 s,
copy_and_validate 2.239933 s,
indexes_and_derived 0.548047 s,
recheck_files 0.007211 s; load/commit
2.823386 s. Per-resource projection/COPY: connectomedb2025 0.043 s; kegg 0.537 s; macdb 0.045 s; metatlas 0.097 s; reactome 0.148 s; recon3d 0.045 s; rhea 0.260 s.
Index creation and initial derivation were measured together by the loader;
no source delete, preparse or canonicalization phase ran in PostgreSQL.
PG-only product rebuild 0.944530 s; full load/check script
5.039105 s. KEGG was the sample's copy outlier (largest fan-out).

**Outcome**: all seven capped native builds passed strict manifest/schema/hash
inspection. The new schema committed successfully; every projected row matched Parquet
(1,482 entities, 19,781 identifiers, 1,974 relations, 2,226 evidence occurrences,
5,772 annotations). All 2,232 raw rows were verified then discarded. Reaction,
COSMOS, MetSigDB and network checks passed with the Parquet directory hidden;
source artifacts were subsequently restored and revalidated. MetaLinksDB had no
eligible server sample rows; populated cases passed separate integration fixtures.
Original pilot metadata, all 21 original resource artifacts and production
container identities/start times were preserved. These private capped versions
must not replace full production artifacts. See the
[validation report](source-attributes-validation.md) for counts and scope limits.


## 2026-09-30 — nicesrv resolved-Parquet pilot, attempt 2

**Reason**: repeat the identical pinned pilot after fixing statistics timing
identified by the cancelled first attempt.

**Parameters**: code `fa439b380c692de7531aa9dc067c9a3843bacb3b`; private release
`2026.9.30.1`, SIGNOR/ChemOnt/Rhea each `2026.9.5.17`;
`batch_size=1024`, DuckDB one thread / 256 MB. Same isolated new container,
volume and localhost port `5440`, database `omnipath_migration`,
fresh schema `pilot_20260930`; 2 CPUs / 4 GiB RAM. Exact configuration,
artifact hashes and scripts are retained in the private server pilot directory.

**Phase durations**: validate_files 0.13 s, copy_and_validate 313.03 s, indexes_and_derived 384.53 s, recheck_files 0.13 s;
total successful load including commit 697.88 s.
Per-resource projection/copy: chemont 1.33 s, rhea 174.18 s, signor 31.64 s.
No delete, extraction/preparse or canonicalization ran. Independent verification
took 176.28 s. Rhea is the source copy outlier.

**Outcome**: committed pilot; every base record/raw payload matched its original
Parquet, all flat-array, ontology, reaction, index and constraint checks passed.
Database size 7.70 GiB from 134.44 MiB compressed input.
Production PostgreSQL and prototype services were unchanged. Complete evidence
and scope limits are in [the pilot report](nicesrv-parquet-pilot.md).

## 2026-09-30 — nicesrv resolved-Parquet pilot, attempt 1 (rolled back)

**Reason**: test the migration loader against existing prototype artifacts before
expanding to a complete release. No resource extraction or resolution ran.

**Parameters**: migration commit `0417420`, private release `2026.9.30.1`
pinning SIGNOR, ChemOnt and Rhea at `2026.9.5.17`; source data
`/root/projects/full_parquet/data`. Loader `batch_size=1024`, DuckDB one thread
and 256 MB limit. A new PostgreSQL 16 container `omnipath-migration-postgres`
uses only its new volume `omnipath-migration-postgres-data` and localhost port
`5440`, database/role `omnipath_migration`, schema `pilot_20260930`; 2 CPUs,
4 GiB memory, 256 MB shared buffers, 16 MB work_mem. Private scripts, manifest,
credentials and logs are in `/root/projects/omnipath-migration/pilot-20260930`.

**Phase durations**: preflight schema/footer/hash validation 0.14 s;
ChemOnt copy/projection 1.31 s, Rhea 174.87 s, SIGNOR 30.28 s; integrity
validation 115.67 s; ordinary indexes 15.85 s; ontology derivation 12.96 s.
Reaction derivation's first FETCH remained blocked for about 262 s and the
attempt was cancelled at 624.21 s total. Rhea is the copy outlier (842,046
payload rows). At 140 s into the reaction FETCH, sampled container block I/O
was 70.4 GB read / 21.7 GB written; those are cumulative container counters,
not reaction-only I/O. Nine input files total 140,970,841 bytes.

**Outcome**: base integrity checks passed, but fresh-table statistics were only
collected after derivation. The reaction join therefore ran without column
statistics and did excessive disk work. Cancelled only the pilot backend;
the transaction rolled back and the pilot schema is absent. The production
PostgreSQL and prototype containers were not accessed for SQL or changed.
Attempt logs are retained as `attempt1.log` and `attempt1-result.json`.
Moving ANALYZE before integrity joins/derivation and including derived tables
in final statistics collection is the first correction; a retry will determine
whether additional query changes are needed.

## 2026-09-10 — cycle 011's final full uncapped rebuild

**Reason**: land the two fixes the cycle had deliberately deferred to its
final rebuild (`annotation_object_entity`'s unpromoted duplicate chemical
entities, `chemical_fallback.py`'s unnormalized name-tier canonical key),
then verify the whole cycle (WP1-WP8) against `contracts/coverage-
acceptance.sql` on one clean, complete build — the gating step before cycle
011 can be considered closed.

**Parameters**: `make reset-content` then `make reload` (`RELOAD_EXISTING=1`,
all sources, uncapped, default `LOAD_JOBS=1`/`THREADS=4`/`BATCH_SIZE=50000`)
against `DATABASE_URL=postgresql://omnipath:omnipath@chemres1-omnipath-build-postgres:5432/omnipath`
and `OMNIPATH_BUILD_UTILS_PG_URL=postgresql://omnipath:omnipath_utils_chemres1-utils_dev@omnipath-utils-chemres1-utils-db:5432/omnipath_utils`,
then `make derive`, then (`omnipath-metabo`) a one-time
`TRUNCATE metabo_lipid_name_resolution` (the ~1.9M-row incremental cache still
carrying pre-T121 chain counts) followed by `post_build_metabo(force=True,
conflicts=False)`.

**Four real bugs found and fixed along the way** (all committed, all
verified live before proceeding to the next step):

1. `84a3ff1` — `annotation_object_entity` self-typed a chemical annotation
   object's raw `(object_id_type, object_id)` instead of promoting it through
   the resolver, minting permanent unpromoted duplicates (the deferred fix
   this rebuild existed to land). Fixed via a new `annotation_object_resolution`
   table (mirrors `ontology_relation_endpoint_key`'s existing pattern for "a
   bare external identifier that is not a full mention"); three other read
   sites needed the identical fix to stay consistent with the entity this now
   mints (`sources` attribution, `annotation_projected`,
   `pq_annotation_relation_evidence_resolved`).
2. `d5c2cda` — WP5's T092 (widening the resolver's chemical namespace list to
   drugbank/reactome) registered both in `RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE`
   but never in `resolver/identifier_types.py`'s separate, must-stay-in-sync
   stable-id registry every resolver setup call validates against. Crashed
   **every single dataset** identically on the first reload attempt
   (`source_rows=0` across all 45 sources, `ValueError: Unknown resolver
   identifier type: 'Drugbank:MI:2002'`) — caught immediately since a full
   reload is exactly the kind of run T092's own narrower verification never
   exercised.
3. `ef6a75b` — `structure_consistency_finding`/`_summary` (WP6, T098) were
   never added to `CONTENT_TABLES`; `reset-content` truncates the whole list
   in one statement, so the unlisted table (with a `data_source` FK) made
   Postgres refuse the entire truncate.
4. `a141963` + `1667fc3` — two separate `vocab_identifier_type` seeding bugs,
   both only reachable on a database with real build history (this sandbox's
   table accumulates across every build ever run here — it is not a content
   table, `reset-content` never touches it): `_ensure_static_identifier_types`
   unconditionally upserted by id, which either crashed or would have silently
   renamed an unrelated row wherever a newly-statically-listed name (T092's
   drugbank/reactome, T116's lipid name) already held a different id from an
   earlier partial run; separately, WP7's own lipid-name seed's
   `... FROM vocab_identifier_type WHERE NOT EXISTS (...)` guard didn't
   actually skip the insert when false (an aggregate without GROUP BY always
   returns one row), so it still attempted an insert with a bogus id and
   crashed on the same already-populated row. Both made purely additive/
   idempotent; verified live as true no-ops against the populated table
   before rerunning.

One further **non-code** issue: the first reload attempt ran clean for
1h48m (through chembl/foodb/go/intact/kegg/pfocr/stitch) then hard-crashed on
a `UnicodeDecodeError` copying a SwissLipids shard to Postgres. Root cause: a
stale preparse cache (`pypath-data/swisslipids/preparse/lipids/batch_50000_all`,
dated 2026-09-08, predating this cycle's SwissLipids `InChIKey=` prefix fix)
had corrupted bytes baked in from before the fix landed — not a code bug.
Deleting the cache directory (forcing a fresh preparse, correctly applying
the `encoding='latin-1'` the loader already declares) resolved it; verified
with a small scoped reload before restarting the full one. This also
suggests the earlier-documented "SwissLipids reload blocked by a `cachedir`
bug" note may no longer reproduce — a fresh reload succeeded cleanly this
time with no sign of that error.

**Phase durations** (second, successful full attempt): reload
**6536.2s (~1h49m)** — `sources=45 skipped_sources=0 datasets=100
failed_sources=3 failed_datasets=5 source_rows=9,920,814
identifiers=63,614,139 annotations=92,029,431` (matches an earlier cycle
baseline almost exactly); the 5 failed datasets are the already-documented
pre-existing bugs (`bindingdb.interactions`, `metatlas.metabolites`,
`mirbase.matures`, `mirbase.precursors`, `ptfi.foods`) — none new. `derive`
**~1h15m** (chem-resolution-level 232s: `chemical_entities=2,401,031
members=7,203,093 groups=6,877,860 relations=15,421,707`; ambiguous-name
candidates 71s; bitmaps 589.6s). `post_build_metabo` **4816.0s (~80m)** —
structure substrate 2,196,732 molecules; Goslin lipid labels on 1,027,283
entities (1,124,040 names resolved) with **`partially_specified=6,258`
appearing as its own real level for the first time** (direct, real-data
confirmation of WP7's T109 fix — the removed species-downgrade used to
collapse every one of these); QC structure-consistency: internal=503,958
cross_reference=85,974 cross_reference_pair=21,125.

**Outcome — the full `coverage-acceptance.sql` gate**:

| SC | Check | Baseline → target | Result | Verdict |
|---|---|---|---|---|
| SC-001 | ChEBI-canonical entities | 54,323 → ≤23,000 | 50,301 | fail (moved from 51,832; not expected to close alone, per `deferred-items.md`) |
| SC-002/003 | Reactome structure reach | ~0.1% → ≥80% | 63.6% | fail (known, legitimate ChEBI class-node structurelessness, not a bug) |
| SC-002/003 | zero-reach group ≥60% | — | rhea 81.9%, pfocr 81.9%, intact 80.6% | pass (tcdb has <100 chemical entities, below the query's own threshold) |
| SC-004 | Reactome↔HMDB shared | 4 → ≥1,300 | 484 | fail (expected; not this fix's target) |
| SC-005 | unresolved chemical entities | 13,631 → ≤5,500, accounted | 365, 365 accounted | **pass** |
| SC-006 | KEGG structure reach ≥70% | — | 81.9% | **pass** |
| SC-006 | KEGG pair-check agreement ≥90% | — | no rows (kegg is itself a structure authority, so it never appears as a *structureless citing* party under WP6's redesigned pair check — a scoping mismatch with this criterion's original, pre-redesign assumption, not a new failure) | n/a |
| SC-007 | every identity syntactically valid | — | 0 invalid | **pass** |
| SC-007 | no identity exceeds the hub threshold | — | largest = 2,650 | **pass** |
| SC-011 | no lipid name collision across chain assignment | — | 0 | **pass** |
| SC-012 | name-canonical entities | 171,852 → ≤40,000 | 47,282 | fail, but close (a ~3.6x reduction from baseline; the normalize_name fix worked, remaining gap is other near-duplicate forms a simple case/whitespace fold doesn't catch) |
| SC-013 | manifest carries coverage | — | true | **pass** |
| SC-014 | db size ≤82GB × 1.10 | — | 91 GB (budget 90.2 GB) | fail, narrowly (~1% over; WP6/WP7 added more than R11's ~3% estimate) |

SC-008/009/010/015/016 are verified elsewhere per the acceptance file's own
header (contract tests, not this SQL gate).

**Net**: 6 of 13 directly-gated criteria pass outright (SC-005, SC-006a,
SC-007 ×2, SC-011, SC-013); SC-001/004/002-003/012/014 remain open, all
either already-documented-as-out-of-scope for these two fixes or improved
substantially without fully closing. No new, unexplained acceptance
regressions.

---

## 2026-09-09 — omnipath-utils WP5: RefMet/Reactome loads, HMDB extension, SwissLipids fix, exemptions (T084-T094)

**Reason**: verify spec 011 WP5 (namespace coverage) against the live
`chemres1` utils DB — a different database
(`omnipath-utils-chemres1-utils-db:5432/omnipath_utils`) than
`omnipath-build`'s own, logged here per the standing rule of recording
every real run against a project database.

**Parameters**: not a `make` invocation — direct Python calls (via
`DatabaseBuilder`) against
`postgresql://omnipath:omnipath_utils_chemres1-utils_dev@omnipath-utils-chemres1-utils-db:5432/omnipath_utils`:
`_long_refmet()`, `_populate_reactome()`, `_long_hmdb()` (re-run with the
CAS/IUPAC/traditional-IUPAC extension), `_clear_wrong_content_rows()`,
`populate_namespace_exemptions()`, plus a one-time SQL cleanup of
already-stored SwissLipids rows and two runs of
`rebuild_resolver_chemical_duckdb` (the first with a since-fixed Reactome
bug, so its output was discarded and rebuilt).

**Phase durations**: RefMet load (208,170 raw rows fetched from a 22.9 MB
CSV, cached fast) a few seconds. Reactome load (`reactome_chebis()`, a
streamed TSV) well under a minute. HMDB re-run (2.2M+ metabolite records,
the largest single source) a few minutes. `resolver_chemical` rebuild via
DuckDB: ~65s to compute 132.8M rows, ~90-100s to write the Postgres staging
table, then a native `CREATE MATERIALIZED VIEW` + index over that — the
first rebuild attempt's own `DROP TABLE IF EXISTS` (lacking `CASCADE`)
failed because the old materialized view still depended on the staging
table; fixed by dropping the old view `CASCADE` first (a pre-existing
script limitation, not something this cycle introduced).

**Outcome**: RefMet — 208,170 `refmet_id<->name` rows plus 98,502
cross-reference rows to chebi/hmdb/kegg/lipidmaps/pubchem (a raw-CSV
quirk cost one round trip: RefMet's first CSV header column carries a
leading space, `" refmet_id"` not `"refmet_id"`, silently producing zero
rows for the primary id until caught and fixed). Reactome — a new
`reactome` id_type, 38,329 `chebi<->reactome` `id_mapping` rows (a second
bug caught and fixed: Reactome's raw `chebi_id` is a bare number, but
`id_mapping`'s own chebi convention is `CHEBI:NNNN` — verified against
existing `bigg->chebi` and `chebi->inchikey` rows before fixing), 31,406 of
which now bridge to a real InChIKey through the rebuilt
`resolver_chemical`. DrugBank — no new load needed (already had 15,293
`id_mapping` rows from an earlier run), now reads 14,138 rows through
`resolver_chemical` after `duckdb_load.py`'s namespace-slug list was
widened (T092). SwissLipids — 777,811 stored `swisslipids -> inchikey`
rows cleaned to 593,217 well-formed keys (184,509 `InChIKey=none`
placeholders deleted, 592,254 prefixed values stripped in place, 85 other
garbage values deleted); the loader fix itself is verified by a unit test
since the live reload is still blocked by the pre-existing `cachedir` bug.
HMDB — extended to also reach CAS (15,672-216,844 rows depending on the
pairing), IUPAC name and traditional IUPAC name, fields its raw records
already carried but the loader never read. Rhea/PDBe — 3,219 and 242,934
wrong-content `id_mapping` rows respectively cleared (UniChem's Rhea
"compound id" is really the participant's ChEBI id; its PDBe "id" is a
display label plus conformer type, not the ligand code), and both excluded
from future auto-discovery. 21 namespaces recorded as reviewed exemptions
(`namespace_exemption`, new table) — id_types.yaml declares each reachable
via hmdb/ramp/chalmers_gem, but a live check of each source's actual raw
data confirms none of the fields exist. Full namespace-by-namespace
reasoning in `saezverse` `deferred-items.md`.

**Verification**: every declared small-molecule namespace (53 total) now
has either real `id_mapping`/`id_mapping_long` rows or a recorded, reviewed
exemption, and no exemption is stale (both checked directly via SQL against
the live database — the equivalent `pytest` suite,
`tests/test_namespace_coverage.py`, is correct but impractically slow over
this DB's network round-trip time, ~5-14s per parametrized case across 60
cases; the two fast, non-DB unit tests in `test_swisslipids_backend.py`
pass in under half a second).

---

## 2026-09-09 — omnipath-metabo WP6: structure consistency diagnostics (T098-T105)

**Reason**: verify WP6's three structure-consistency checks
(`omnipath-metabo/omnipath_metabo/postbuild/_qc_layer.py`) against the live
build — a different repo/service (`omnipath-metabo`, freshly cloned this
session) than `omnipath-build`'s own `make load`/`make reload`, but the same
target database, so logged here per the standing rule of recording every real
run against it.

**Parameters**: not a `make` invocation — direct Python calls against
`postgresql://omnipath:omnipath@chemres1-omnipath-build-postgres:5432/omnipath`,
schema `public`: `build_structure_consistency_findings(conn, schema='public')`
(three times, iterating on a performance fix) and `post_build_metabo(conn,
schema='public', force=True, conflicts=False)` (the full post-build pipeline,
including the new QC step).

**Phase durations**: first version of the cross-reference/cross-reference-pair
queries — driven from `entity_evidence` (all chemical mentions, ~tens of
millions across 43 partitions) before filtering down to structure-bearing
ones — did not finish in 15 minutes (`timeout 900` killed it). Rewritten to
start from the ~2.9M structure-authority-typed `identifier_evidence` rows
(and a precomputed authority-identifier-to-entity map, built once, reused by
both checks): full three-check pass **671s (~11m12s)** —
`internal=1,199 cross_reference=44,303 cross_reference_pair=16,979
summary_rows=194`. The full `post_build_metabo(force=True)` pipeline
(substrate rebuild + specificity + facet + lipid labels + QC, no RaMP
conflicts): structure substrate **2,223,612 molecules**; specificity
`cis_trans_only=102,451 constitution_only=752,411 no_structure=288,224
stereospecific=1,173,513 unknown_constitution=387,377
variable_constitution=195,237`; Goslin lipid labels on 1,072,117 entities
(1,125,099 names resolved, 770,019 unresolved).

**Outcome**: the checks are correctly derived from `identifier_authority`
(no hand-listed sources), classify sensibly (`different_structure` well
under half of every summary row, e.g. a real `refmet`/`kegg` namespace
disagreement at 3.3% agreement against markedly-higher agreement for other
pairs — the same class of defect research.md R8's original KEGG
investigation found), and are served with no request-time computation
(`test_qc_layer.py`, 13/13 passing). The internal check is narrower than
`contracts/quality-control-api.md` describes — this build's rdkit cartridge
has no InChI support (`mol_from_inchi` does not exist; `mol_inchikey`
returns the literal string `InChI not available`), and `omnipath-metabo`
deliberately runs chemistry only through the cartridge, never Python rdkit —
so it compares a record's own duplicate SMILES assertions against each
other, not SMILES vs. InChI vs. InChIKey. Full writeup in `saezverse`
`specs/011-chemical-resolution/deferred-items.md`.

**Self-inflicted incident, noted for the record**: an intermediate
`test_postbuild.py::test_post_build_refuses_on_build_id_mismatch` run was
killed by its own 30-minute `timeout` wrapper mid-way through a `force=True`
rebuild (the QC step roughly tripled that test's total cost, since it calls
`post_build_metabo` three times). This left `metabo_build_state` stuck at
the test's injected sentinel `deadbeefdead` and an orphaned
`CREATE INDEX ... metabo_entity_structure_mol_idx` still finishing
server-side (harmless — any subsequent `build_structure_substrate` call
drops and rebuilds it anyway). Fixed with one more `force=True` run once the
orphaned index finished; not a code defect.

---

## 2026-09-09 — derive run for entity_name/T072-T074, and the shared-database rename

**Reason**: verify T072-T074 (entity_name table, its population, the new
label-cascade pass) against real data by running `derive` on the live build
DB, per T065's test. Mid-run, a colleague reported that `netds1` and this
sandbox (`chemres1`) appeared to be sharing the same build database
container. Confirmed on `beauty`: both sandboxes' Docker Compose projects
defaulted to the same project name (`omnipath-build`, from the directory
basename), so both had been writing to the exact same
`omnipath-build-postgres-1` container and volume the whole cycle, unknown
to either side. Full story in the `shared-db-incident` memory.

**Parameters**: `make derive` (no `MAX_RECORDS`, full uncapped) against
`DATABASE_URL=postgresql://omnipath:omnipath@omnipath-build-postgres-1:5432/omnipath`
(the hostname *at the time this ran* — see below, it no longer resolves).
Launched ~17:07 UTC.

**Phase durations**: total `event=done seconds=3073.780` (~51m14s),
`failed_steps=` empty. Major phases: derive-tables (schema/PK/index setup
plus `entity_identifier_lookup`/`entity_relation_counts`/
`entity_ontology_term`/`entity_source_count`) 351s; chemical resolution
levels 185s; chemical ambiguous-name candidates 62s; `classify_chemical_class`
308s; `classify_metabolic_domain` 439s; interactions (the largest single
phase) 679s; `entity_labels` 145s; **`entity_name` 77s, 4,157,246 rows
written**; **`chemical_labels` 158s**, `chemical_name=1,194,082
chemical_iupac_name=4,605 chemical_preferred_name=0 chemical_identifier=1,699,773
without_real_label=753`; bitmaps (the second-largest phase) 452s; metsigdb
158s.

**Outcome**: `entity_name`/T072-T073 confirmed working against real,
full-scale data. `chemical_preferred_name=0` is expected, not a bug — the
new T074 pass is a safety net for names `entity_name` reaches but the
pre-existing tiered cascade doesn't; today both read identical underlying
data, so pass 1 already covers everything pass 1b would have caught.
`tests/test_preferred_label.py` (T065) 4/4 GREEN against this run's
resulting state.

**Database hostname changed as a result of the shared-DB fix**: the
container was renamed `omnipath-build-postgres-1` →
`chemres1-omnipath-build-postgres` and its old network aliases dropped
(disconnect+reconnect on the running container — Postgres itself was never
restarted, no data touched). Every command above used the *old* name,
still correct history for what actually ran; **any new command from this
point on must use the new name** — `db-access.env` and the `cycle011-status`
memory are both already updated. `netds1`'s own isolation (a distinct
`COMPOSE_PROJECT_NAME`) is still pending on their side as of this writing.

---

## 2026-09-09 — WP2 full uncapped rebuild (spec 011, cycle 011)

**Reason**: WP2 (T048-T060, candidate arbitration — skeleton collapse,
neutral-parent/least-specified-stereo tie-break, authority arbitration,
`resolution_conflict` recording) is implemented and unit-tested
(`main@8d9b0f1`), but not yet checked against real, full-scale data. This
run exists to answer: does WP2 actually move SC-001/002/003/004/005 past
their targets now, and does T061/T062's threshold hold (unresolved <5,500,
KEGG structure reach >70%)?

First had to fix a database schema issue found while preparing this run:
`entity`, `entity_evidence`, `identifier_evidence`, `relation`,
`relation_evidence`, and `annotation` all predated their declared
`PRIMARY KEY` on this database (`CREATE TABLE IF NOT EXISTS` is a no-op on
an existing table, so a table created before the PK was added to its DDL
never retrofits one) — blocked every `make db-setup`. Separately,
`gene_protein_representative`/`state`/`state_component`/`evidence_state`
were never in `CONTENT_TABLES`, so `reset-content` never truncated them;
stale rows survived every prior wipe this cycle and eventually collided
with `evidence_state`'s own PK add. Both fixed in `main@b7e16b0`, verified
with a clean `make reset-content` (`entity` count confirmed 0) followed by
two consecutive clean `make db-setup` runs.

**Parameters**: identical to the two T047 runs below —

```
make reload \
  DATABASE_URL="postgresql://omnipath:omnipath@omnipath-build-postgres-1:5432/omnipath" \
  OMNIPATH_BUILD_UTILS_PG_URL="postgresql://omnipath:omnipath_utils_chemres1-utils_dev@omnipath-utils-chemres1-utils-db:5432/omnipath_utils"
```

All ~47 sources, uncapped (`--max-records 0`), `--stage-jobs 1`,
`BATCH_SIZE=50000`, `THREADS=4`. Only Postgres content was wiped
(`reset-content`), not the on-disk preparse shard cache, so preparse reused
the cache from the 2026-09-08 runs and was fast, same as the second T047
run below.

**Finished** ~13:58 UTC. Pipeline's own reported total: `sources=45
skipped_sources=0 datasets=100 failed_sources=3 failed_datasets=4
source_rows=9920814 identifiers=63614139 annotations=92029431
total=6455.728s` (~1h47m36s) — same 4 pre-existing unrelated failures as
both 2026-09-08 runs (bindingdb.interactions, metatlas.metabolites,
mirbase.precursors, ptfi.foods), plus the already-documented
`hormone2cell` discovery failure. `source_rows`/`identifiers`/`annotations`
are byte-for-byte identical to the second 2026-09-08 run, as expected —
same source data, cache reused.

**Outcome — first four acceptance blocks, before vs. after WP2**:

| Criterion | Target | Before WP2 (both 2026-09-08 runs) | After WP2 (this run) |
|---|---|---|---|
| SC-001 ChEBI-canonical | ≤23,000 | 51,832 (fail) | 51,832 (fail, unchanged — tied to the deferred `annotation_object_entity` bug, not WP2's scope) |
| SC-002/003 Reactome structure reach | ≥80% | 63.1% (fail) | 63.7% (fail, ~unchanged — legitimate ChEBI generic/class-node structurelessness, confirmed non-bug) |
| SC-002/003 intact/rhea/pfocr | ≥60% | pass | pass (80.6% / 81.9% / 81.9%) |
| SC-004 Reactome↔HMDB shared | ≥1,300 | 493 (fail) | 494 (fail, unchanged — same root cause as SC-001) |
| **SC-005 unresolved chemical entities** | **≤5,500** | **11,424 (fail)** | **365 (PASS)** |
| **SC-006 KEGG structure reach** | **≥70%** | not yet passing | **81.9% (PASS)** |

T061/T062's two thresholds (unresolved <5,500, KEGG structure reach >70%)
both cleared. SC-001/002-003(reactome)/004 remain exactly where they were —
expected, since none of them are in WP2's scope (they trace to the deferred
`annotation_object_entity` bug and legitimate ChEBI class-node
structurelessness, both already root-caused in the 2026-09-08 entry below).

**Bug found via this run's own diagnostic column**: SC-005's
`with_conflict_record` came back 0 despite 11,424→365 unresolved entities
and WP2's own unit tests proving the `resolution_conflict` CTE logic is
correct. Root cause: `resolution_conflict` was created correctly inside
each shard's DuckDB session, but was never added to
`STAGED_LOAD_TABLES` (`duckdb_direct_pipeline.py`), the `pq_*` view mapping,
or `_bulk_copy_canonical`'s COPY step (`duckdb_load.py`) — so it was
silently discarded before ever reaching Postgres. `resolution_conflict` was
0 rows database-wide after this run despite genuine candidate disagreements
existing. Fixed in `main@4e6f5c2`; verified against real data with a
standalone `kegg`-only reload (35,307 correctly-translated conflict rows
written for kegg alone). Full rebuild re-run below to backfill this for
every source.

---

## 2026-09-09 — WP2 rebuild #2, after fixing the resolution_conflict copy bug

**Reason**: the run above proved WP2's resolution logic works
(SC-005/SC-006 both clear target), but `resolution_conflict` was silently
empty database-wide due to the COPY bug fixed in `main@4e6f5c2` (see above).
A `kegg`-only standalone reload confirmed the fix works (35,307 rows for
kegg alone), but left the database in a mixed state — only `kegg` carries
`resolution_conflict` rows, the other 44 sources still reflect the buggy
pre-fix run. Re-running the full rebuild for a coherent, fully-correct
state and an accurate `with_conflict_record` figure.

**Parameters**: identical to the run above, preceded by a verified-clean
`make reset-content` (`entity` count confirmed 0).

**Finished**: pipeline's own reported total: `sources=45 skipped_sources=0
datasets=100 failed_sources=3 failed_datasets=4 source_rows=9920814
identifiers=63614139 annotations=92029431 total=7213.428s` (~2h0m13s — a
little slower than run #1's 1h47m36s, ordinary run-to-run variance on the
shared host, not a regression). `source_rows`/`identifiers`/`annotations`
identical to run #1 and to both 2026-09-08 runs, as expected — same source
data, same resolution logic, only the `resolution_conflict` COPY path
changed. Same 4 pre-existing unrelated failures again
(bindingdb.interactions, metatlas.metabolites, mirbase.precursors,
ptfi.foods).

**`resolution_conflict` now populated for every source with real
disagreements**, not just `kegg`:

| source | conflict rows |
|---|---|
| kegg | 35,307 |
| recon3d | 23,514 |
| metatlas | 3,472 |
| wikipathways | 2,117 |
| pmidb | 416 |
| phenol_explorer | 268 |
| intact | 168 |
| mrclinksdb | 68 |
| chebi | 12 |
| lipidmaps | 6 |
| mebocost | 4 |

Total 65,352 rows across 11 sources.

**Outcome — full re-run of the acceptance gate**: identical to run #1 on
every criterion (SC-001 51,832 fail, SC-002/003 Reactome 63.7% fail /
intact+rhea+pfocr pass their own bar, SC-004 494 fail, SC-005 365 unresolved
pass, SC-006 KEGG 81.9% pass) — confirming the copy-bug fix only affected
the `resolution_conflict` diagnostic table, not any resolution outcome.
**SC-005's `with_conflict_record` is now 365, exactly matching `unresolved`
(365)** — every unresolved chemical entity is accounted for by a recorded
conflict, closing the gap run #1 exposed. This is the fully-correct,
coherent state of the database for cycle 011's WP2 work: T061/T062 both
clear target, and the "recorded, never silently dropped" guarantee the
spec's User Story 2 asks for is now actually true in Postgres, not just in
the DuckDB session that computes it.

---

## 2026-09-08 — T047 full uncapped rebuild (spec 011, cycle 011)

**Reason**: T046's capped reload (`reactome,hmdb,chebi,kegg`, `MAX_RECORDS=5000`)
surfaced and fixed two real bugs beyond WP1's own T027-T045 —
`omnipath-build main@d879ba9` (resolver_candidate joined evidence's raw
identifier instead of the T040/T041-normalized one) and `main@656ad24`
(`multigene_split._explode_one` silently dropped `identifier_normalized` when
rebuilding `entity_identifier_raw` for shards with a multi-gene UniProt
mention). Both verified fixed on the capped run (Reactome structure reach
~0.1% → ~86.6%, Reactome↔HMDB shared metabolites 0 → 263 in the 5,000-row-
capped sample). T047 is the uncapped, all-sources run needed to actually check
`contracts/coverage-acceptance.sql`'s first four blocks (SC-001, SC-002/003,
SC-004, SC-005) — the capped numbers are indicative, not the acceptance gate
itself.

**Parameters**:

```
make reload \
  DATABASE_URL="postgresql://omnipath:omnipath@omnipath-build-postgres-1:5432/omnipath" \
  OMNIPATH_BUILD_UTILS_PG_URL="postgresql://omnipath:omnipath_utils_chemres1-utils_dev@omnipath-utils-chemres1-utils-db:5432/omnipath_utils"
```

i.e. Makefile defaults: `SOURCES=` (all ~47 discovered source modules, not a
named subset), `MAX_RECORDS=` (uncapped, `--max-records 0`), `BATCH_SIZE=50000`,
`THREADS=4`, `LOAD_JOBS=1` (`--stage-jobs 1`, sequential staging — not tuned up
for this run; `build-tuning.md` documents the utils PG is provisioned for
`LOAD_JOBS` up to 32, worth trying on a future uncapped run if wall-clock time
matters), `--reload-existing --append`. Started ~10:27 UTC.

**Phase durations** (this run is a long alphabetical sweep through sources;
timings below are what was actually observed, not estimates, unless marked):

- Per-source delete (existing content wipe before reload): fast, not a
  bottleneck — ran through all ~47 sources in the first few minutes.
- Preparse pass (raw pypath records → cached shard files on disk, one source
  at a time, `--stage-jobs 1`): the dominant cost so far.
  - `chembl.molecules`: **3211s (~53.5 min)** end to end — this includes
    downloading `chembl_36_sqlite.tar.gz` (5.23 GB) from the upstream EBI
    mirror at a bursty, throttled ~700 kB/s–6.5 MB/s (not a steady rate),
    producing 1,422,531 rows across 29 shards.
  - `chembl.activities` (same download, separate dataset): 3,459,458 rows
    across 70 shards, exported as 50k-row Parquet chunks — the single largest
    dataset preparsed so far.
  - `chembl.mechanisms`: 7,568 rows, fast.
  - `intact.interactions`: `human.zip` (1.09 GB) download alone took ~9 min at
    a similarly throttled rate; 1,240,727 rows.
  - `stitch.mouse_interactions`: still building as of the last check (350,000+
    rows across 6+ shards, 330s+ elapsed for this one dataset) — noticeably
    slower than `stitch.human_interactions` (338,521 rows, fast).
  - By contrast, the acceptance-relevant chemical sources preparsed quickly
    once reached: `chebi.molecules` 216,787 rows, `kegg.reactions` 87,381 rows,
    `reactome` (reactions 14,953 / pathways 2,883 / controls 10,080 /
    control_groups 1,130), `rhea` and `pfocr` — all sub-minute.
  - Preparse pass finished entirely at 12:26 UTC (~1h59m elapsed): 103
    (source, dataset) pairs preparsed across all ~47 sources, alphabetically
    `bindingdb` → `wikipathways`.
- Stage/canonicalize/copy pass (the real chemical-resolution pipeline —
  project → canonicalize → COPY to Postgres, one dataset/shard at a time,
  50k-row shards): began 12:26 UTC, much faster throughput than preparse.
  `task_done` reached 40 within 5 minutes of starting (vs. 95 datasets over
  ~2h for preparse) — most individual shard canonicalize times are single-digit
  seconds. `chebi.molecules` (216,787 rows, the acceptance-critical source for
  T040/T041/d879ba9) staged and canonicalized cleanly across 5 shards, 2–16s
  canonicalize each, **no BinderException** — the resolver-candidate fix is
  holding under real, uncapped data. All 7 acceptance-critical sources
  (chebi/hmdb/reactome/kegg/rhea/pfocr/tcdb/intact) finished staging clean,
  no failures, no `BinderException` anywhere in the run.

**Finished** 13:55 UTC, **total 3h27m39s** (12,464.7s). Final summary:
`sources=45 datasets=99 failed_sources=4 failed_datasets=5 source_rows=9,901,436
identifiers=63,468,308 annotations=91,865,478`. 5 failed datasets:
`drugcentral.interactions` (disk — see the follow-up entry below),
`bindingdb.interactions` (`BadZipFile`), `metatlas.metabolites`
(`ArrowTypeError`), `mirbase.precursors` (`IndexError`), `ptfi.foods` (no JSON
files found) — the last four are pre-existing, unrelated to this cycle.

**Outcome**: ran `contracts/coverage-acceptance.sql`'s first four blocks —
**none pass yet**: SC-001 51,832 (target ≤23,000), SC-002/003 Reactome 63.1%
(target ≥80%; `intact`/`pfocr`/`rhea` did clear their own ≥60% bar), SC-004
shared 493 (target ≥1,300), SC-005 unresolved 11,424 (target ≤5,500). All four
show large real improvement over baseline (SC-004 alone is 123x baseline), just
short of target. See the follow-up entry below for the investigation into why,
and two more real bugs found along the way.

---

## 2026-09-08 — T047 clean rebuild, after a full `reset-content` wipe

**Reason**: suspected the 51,832/63.1%/493/11,424 numbers above were inflated
by orphaned entities left over from this sandbox's many prior test reloads
across this whole cycle (24,187 of the 51,832 "ChEBI-canonical" entities have
zero `entity_evidence_resolution` rows referencing them). Wiped all content
tables (`make reset-content`, not just per-source deletes) and reran the exact
same full uncapped reload to test whether a genuinely clean database gives
different numbers.

**What actually happened**: `reset-content` itself was broken —
`CONTENT_TABLES` (the list it truncates) was missing 6 tables carrying foreign
keys into it, so Postgres refused the TRUNCATE. Fixed all 6 in two commits
(`main@bf34258` `interaction_fact_resource`; `main@36544a9` five more found by
querying `information_schema` systematically rather than fixing one crash at a
time: `data_source_license`, `entity_annotation_relation_default`,
`identifier_authority`, `identifier_role`, `resource_overlap_summary`). Also
had to `pg_terminate_backend` a zombie diagnostic query that had outlived my
own `pkill` by over an hour and was holding a lock blocking the truncate.

**Parameters**: identical to the run above —
`make reload DATABASE_URL=... OMNIPATH_BUILD_UTILS_PG_URL=...`, all sources,
uncapped, `--stage-jobs 1` — preceded by a verified-clean `reset-content`
(`entity` count confirmed 0 before launch). Started ~15:45 UTC.

**Phase durations**: dramatically faster — preparse was **instant** (reused
the shard cache from the run above, since only Postgres content was wiped,
not the on-disk preparse cache). `task_done` reached 97 within 23 minutes
(vs. ~2h for the first run to even start real staging). `foodb.foods` was
again the standout slow dataset (16.8min projection this time, vs. 6min in
the first run — 992 rows exploding into 4.66M identifiers, a huge fan-out,
not a bug). All 7 acceptance-critical sources confirmed clean again,
including `tcdb`.

**Finished** 17:38 UTC, **total 1h54m41s** (6,881.5s) — 45% faster than the
first run thanks to cache reuse. `sources=45 datasets=100 failed_sources=3
failed_datasets=4 source_rows=9,920,814 identifiers=63,614,139
annotations=92,029,431`. Only 4 failed datasets this time — `drugcentral`
succeeded (the first run's failure was transient shared-host disk pressure,
not a real bug); the same 4 pre-existing unrelated bugs (bindingdb, metatlas,
mirbase, ptfi) reproduced identically.

**Outcome — the wipe theory was wrong, but led to the real bug.** Re-ran the
acceptance gate: **identical numbers** to the contaminated run (51,832 /
63.1% / 493 / 11,424), down to the exact same orphaned entity UUID
(`1f9f0a48-1292-d05d-6c0b-894283589c84`, CHEBI:107644) reappearing byte for
byte after a complete wipe. Entity ids are content-addressed hashes, so this
proved the orphans are a **reproducible pipeline bug, not stale data**.

Root-caused it: `duckdb_load.py`'s `annotation_object_entity` CTE (inside
`batch_entity_candidate`) self-types an annotation-relation's object by its
raw `(object_id_type, object_id)` verbatim — no resolver join, no InChIKey
promotion, and only CV-term objects are excluded. A chemical annotation
object (e.g. a reaction annotated with a bare ChEBI id, not a full
`entity_evidence` mention) mints its own dead-end entity, permanently
separate from whatever properly-resolved InChIKey entity the same real
molecule gets via the normal mention pipeline. Documented in code
(`main@95ee298`, comment only — **not yet pushed**, GitHub token expired
mid-session, needs a sandbox credential refresh) and in memory as a follow-up,
not fixed this cycle per an explicit decision to defer it.

Separately, investigated the Reactome 63.1% shortfall directly: legitimate,
not a bug — 269 of a 270-id sample of stuck-at-ChEBI entities are either
R-group-wildcard generic templates or have no structural data in ChEBI at
all. Matches the spec's own control that ChEBI class-nodes must never gain a
fabricated structure.

**Infra note**: this session's `/workspace/instances/chemres1` (a CephFS
mount, separate from the container's local `/scratch`-backed overlay) hit
genuine zero-bytes-available for about an hour (14:13-15:25 UTC) mid-wipe —
corrupted `uv.lock`, blocked all writes including `git checkout`. Self-cleared;
nothing fixable from inside the sandbox if it recurs, just poll and wait.

T048+ (WP2, candidate arbitration) now looks clearly necessary to actually
clear SC-001/002/003/004/005's targets, not just a nice-to-have — the
unresolved counts are in the right range to match WP2's already-documented
candidate-ambiguity problem.


## 2026-10-01: main alignment setup and authorized old-schema cleanup

The user authorized using nicesrv and removing only the previous approximately 190 GiB full migration build to make room for alignment with main. No Parquet outputs, resolution reference or cache were changed. A guarded drop removed `full_20260930` from private PostgreSQL16 at loopback5440; both prior pilot schemas and all four production/prototype containers remained intact. The schema occupied206,661,582,848bytes (192.47GiB). Free filesystem space rose to562,186,121,216bytes. A schema-only DDL backup and copied release/build/sample reports are retained in `/root/projects/omnipath-migration/main-alignment-20261001`; cleanup receipt records four matching protected container IDs/start times.

A separate PostgreSQL18 image was built from main9f9bb709c764's Dockerfile, imageSHA256e261f1bed70e7d1001653fa69875e8ce2ea79ed2ef2b99834fffd77d26c32201. Isolated container`omnipath-migration-main-postgres` binds only127.0.0.1:5441 with separate volume`omnipath-migration-main-postgres-data`, memory8GiB, CPU3, shared_buffers1GB, work_mem32MB, maintenance_work_mem512MB and two parallel maintenance workers. Existing privatePG16 and protected services were retained.

Read-only scalar audit of the unchanged pinned46-resource release2026.9.30.4 scanned6,383,473published entity rows. It found no missing canonical type/namespace/identifier fields;2,188keys differed only taxonomy, of which2,181contain one known taxon plus the prototype's empty-string unknown marker and only7contain multiple actual known taxa.39natural-key conflict groups are distinct source-scoped name fallback identities. The adapter retains published identity and original occurrences, selects sole known taxonomy where unambiguous and NULL on actual conflict, and uses explicit published fallback identity only for colliding source-scoped names. No upstream recanonicalization occurred.

Development checks are ongoing. The first private main-oracle test run passed20non-DB checks but two setup tests failed because the migration pypath dependency lacks main's LIPID_NAME CV enum. The oracle now uses the exact CV package from main's pinned pypath33f37fbaab59993d24f5c32bb4b3e7587085bcf4; the migration dependency is untouched. This is a bounded development fixture run, not a full resource build or mandatory post-build gate. Full load has not started.


## 2026-10-01: bounded main-compatible PostgreSQL pilots

All loads below use the isolated PostgreSQL18 instance at loopback5441, unchanged synthetic serving Parquets, DuckDB1thread/512MB, source fixture with3entities/1statement/2evidence entries, and a new destination schema. No resource parser, resolver or real resource rebuild ran. Each base used all20 canonical/evidence/provenance COPY tables with late main keys, FKs and indexes.

First base pilot `aligned_pilot_20261001` succeeded: artifact checks0.003178s, DuckDB projection0.380812s, CSV staging0.019629s, COPY0.011752s, constraints0.125269s and base post-copy0.236685s. Its first shared finish failed at build-manifest cost serialization after successful earlier shared steps; main's exact cost-shaping helpers corrected the type mismatch. PostgreSQL-only finish succeeded without repeated COPY: shared derivations0.241692s, MetSigDB0.012112s, network presets0.005879s, COSMOS0.019649s. The tiny fixture has no eligible MetSigDB/reaction product rows; it tests orchestration rather than product science.

Second pilot `aligned_pilot2_20261001` succeeded through base: artifact0.002926s, projection0.363939s, CSV0.015982s, COPY0.009198s, constraints0.099499s, base post-copy0.194888s. Its shared derive failed when the bitmap extension was found installed under the first pilot schema. The isolated extension was explicitly relocated to public with dependent objects preserved; code now installs extensions in public and reports mismatched existing namespaces without implicit relocation. PostgreSQL-only finish retained the copied base: shared derivations0.244538s, MetSigDB0.011740s, network presets0.005744s, COSMOS0.023541s. All three complete products, subset_build_metadata and parquet_phase committed with the pinned pilot release. A separately selected network product rebuild through the public subsets entrypoint succeeded in0.005278s, with callback after durable commit and other products retained.

Focused local regression suite:393passed/326skipped in5.06s; server-only tests are skipped locally by explicit opt-in. Guarded isolated PostgreSQL fixtures:110passed in7.78s, including exact main catalogue definitions and late indexes/statistics, source-specific nullable/sign facts, merged reactions/transport/stoichiometry, quantities/publications, ontology shortest closure/labels, source-local miRBase maturation, full COSMOS row/statistics parity, protein labels, ChEBI display identity, two-schema bitmap types and actual product checkpoint durability. The server fixture override targets only loopback5441/databaseomnipath_migration and creates disposable random schemas. This is bounded development verification; no mandatory full-release verification/rebuild/rollback/safety phase was reinstated.

Reports are retained under `/root/projects/omnipath-migration/main-alignment-20261001`:pilot-main-base.json,pilot-main-finished.json,pilot2-main-finished.json,pilot2-selected-product.json. Exact main reference9f9bb709c764 and pypath CV dependency33f37fbaab59 are pinned independently from the migration pypath. Genuine unavailable original inputs are declared in build_manifest/build_capability and the column map.

Final frozen main-oracle rerun: 15 passed in 3.27s, including BIGINT source IDs and executed view/function definitions. Ruff and patch whitespace checks passed.


## 2026-10-01: full main-compatible release load started

Runtime commit `d0a1c80c33699353a589e06b957ad49960f01ef2` is pushed on `parquet-migration`. A615-file source snapshot was verified byte-for-byte on nicesrv; its manifest SHA256 is66c69d10b3d4113fd9e00ba1e102443788af970b08d0886eba362dce0fbdd4b2. Copied main SQL retains two inherited whitespace warnings; the implementation's Ruff checks passed.

The full build started at2026-10-01T20:25:13.932388+00:00, using unchanged release2026.9.30.4 (46resources, canonical SHA25606f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e) in new schema`aligned_full_20261001`, isolated mainPG18 at loopback5441. Parameters: DuckDB2threads/2GB, on-disk projection and CSV staging, late main keys/FKs/indexes, durable base and complete-product commits. Systemd unit`omnipath-migration-main-aligned-20261001.service` caps the builder at5GiB/2CPU and runs atNice10. No mandatory verification, rollback-rebuild or safety job is queued.

Preflight confirmed all615runtime files, mainPG18 image/loopback binding, public extension namespaces, no competing active private backend, both prior pilot schemas in privatePG16 and all four running protected services. Free disk560,904,417,280bytes. Initial monitor at20:27UTC showed healthy DuckDBvalidate_and_stage, matching release pins and no committed base yet. Memory was approximately2.34GBanonymous plus reclaimable file cache; zero OOM events. The existing15-minute check is active for this new run and reports meaningful changes only. Full timings, product results and final storage are pending.


## 2026-10-01: full attempt1 diagnosis and taxonomy validation correction

The first aligned full run stopped beforeCOPY at20:34:50UTC after576.513767s. No base, product or COPY checkpoint was committed; the fresh schema contains only setup/seed metadata (approximately1.5MiB). Its reports/events/log/error and source manifest are archived privately under main-alignment-20261001/attempt1-d0a1c80.

Read-only scalar audit of18,339,902published statement rows confirmed exactly35shared statement keys differ only in scalar taxonomy, all betweenChEMBL andDrugCentral. Endpoints, predicate, statement kind, direction, sign, category and interaction class agree. The publisher computes taxonomy as a resource-specific consensus and excludes it from the statement key. The projection check now allows only that legitimate taxonomy variation while preserving both source-owned statements/evidence/taxon annotations; true identity and all other existing consistency checks remain. No Parquet was changed and no arbitrary biological value was selected. The bounded aggregate diagnostic completed in31.611060s after restricting its nine-field detail comparison to the35already identified keys; earlier512MiB/1GiB diagnostic attempts exceeded working memory.

Five new projection regressions verify source-owned taxonomy/evidence and reject changed subject/object/predicate/kind. Focused projection/loader/MetSigDB tests:69passed,6opt-in server tests skipped in3.08s. The new exact-main MetSigDB oracle passed all3tests on private5441 in2.59s:20synthetic source records cover allfive complete memberships/provenance, productDDL/indexes, idempotence and source-local stale-row removal. The failed empty schema will be preserved under an attempt-specific name; restart uses a fresh target and cannot reuseCOPY because none ran.


## 2026-10-01: corrected full attempt2 running overnight

The failed attempt1 schema was confirmed empty of copied scientific tables and renamed to `aligned_full_20261001_attempt1_empty`; its setup metadata and private archived diagnostics are retained. No COPY or committed base existed to resume.

Runtime `bbe8ea40f90f60d5a5bfaad0c99de2baac856fc5` was pushed and deployed as a verified 616-file snapshot, manifest SHA256 `9ccb45f39eb83dd577ef55b508a2ff300a398d685733a2001611a2edec841060`. The corrected full load started at 2026-10-01T20:53:47.904529+00:00 in `omnipath-migration-main-aligned-attempt2-20261001.service`, using a fresh `aligned_full_20261001` schema. The unchanged 46-resource release, private PostgreSQL18 at loopback5441, DuckDB2threads/2GB, builder5GiB/2CPU/Nice10, late main keys/FKs/indexes and durable base/product checkpoints remain as recorded above. Preflight free disk was557,452,517,376bytes.

At21:03:07UTC the run had passed shared-statement validation and reached DuckDB dimensions; no base had committed yet. Free disk538,289,545,216bytes. The user requested overnight progress; the existing15-minute monitor was updated to this single attempt, with authorization to diagnose clear failures, test narrow corrections and preserve completed checkpoints. No mandatory exhaustive verification, rebuild/rollback or safety job is queued. Completion, final phase durations, representative samples and storage remain pending.


## 2026-10-01: attempt2 memory diagnosis and bounded dimension preparation

Attempt2 was deliberately stopped at2026-10-01T21:13:08.295055+00:00 after1160.390526s, still in DuckDB dimensions and beforeCOPY. It had passed the corrected statement identity validation. No base/product checkpoint or scientific COPY rows existed; the setup-only schema is retained as `aligned_full_20261001_attempt2_empty`. Reports/log/source receipt are archived privately under `main-alignment-20261001/attempt2-bbe8ea4-dimension-memory`. Fresh checks confirmed zero remaining builder backends, both original pilots and all four protected services unchanged.

The dimension collector fetched every entity type, predicate, category and annotation scope occurrence before deduplicating in Python. At diagnosis the process used approximately4.9GBanonymous memory under its5GiB service limit, with no OOM kills but heavy memory and disk pressure. The SQL now selects distinct non-null names before fetching. The existing sorted vocabulary allocation, seeded IDs, dataset pairs, scientific projection and published Parquets are unchanged.

A bounded10,000-row-per-input fixture compares exact prior/new dimension and dataset outputs and reduces transferred names from35,006to14. Subagent focused tests used a disposable local PostgreSQL16 cluster: `OMNIPATH_TEST_POSTGRES=1 uv run pytest -q packages/omnipath_postgres/tests/test_aligned_projection.py packages/omnipath_postgres/tests/test_aligned_loader.py packages/omnipath_postgres/tests/test_duckdb_projection.py packages/omnipath_subsets/tests/test_main_adapter.py`;129passed in4.70s, Ruff clean. Root reviewed the patch and its seeded-ID/full-value parity comparator. The next full retry keeps the same pinned release and checkpoint parameters; no mandatory verification or rebuild/rollback job is introduced.


## 2026-10-01: full attempt3 started with bounded dimensions

Corrected runtime `2942d38385c1319f41a4ee1d5559c65604fb0893` is pushed on `parquet-migration`. Its616committed source files were verified on nicesrv; runtime manifest SHA256 `be1269009cf961a20a1bd9b11dfe628e6420ef85edc024d53cba5e94d408fedc`. The single `omnipath-migration-main-aligned-attempt3-20261001.service` started at2026-10-01T21:19:23.787442+00:00 in fresh `aligned_full_20261001`, same unchanged release2026.9.30.4/46resources/canonical SHA25606f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e, privatePG18loopback5441, DuckDB2threads/2GB, builder5GiB/2CPU/Nice10 and durable base/complete-product commits. Preflight free disk538,289,389,568bytes; both prior pilots, public extension placement, private image/loopback and protected services were checked.

Initial progress at21:19:46UTC: validate_and_stage, healthy active unit, matching release pins and no base yet. The existing15-minute heartbeat now follows attempt3 and retains prior failure/cancellation history. A separate bounded read-only manual sample helper is prepared for completion; no mandatory verifier/rebuild/rollback/safety phase was added. Full phase times, actual product results, final sample and storage remain pending.


## 2026-10-01: attempt3 memory limit and equivalent entity aggregate

Attempt3 failed beforeCOPY at2026-10-01T21:29:32.914704+00:00 after609.127313s with DuckDBOutOfMemoryException in the entities phase (1.8GiB/1.8GiB pool used). Dimension preparation completed in3.421316s, confirming the SQL deduplication fix. The exact failing SELECT was not captured, so no more specific attribution is claimed. No base/product checkpoint or scientific COPY rows existed. The setup-only schema is retained as `aligned_full_20261001_attempt3_empty`; diagnostics/reports are archived under `main-alignment-20261001/attempt3-2942d38-duckdb-memory`. All four protected service IDs/start times matched their baseline.

The sole-known taxonomy aggregate now usesMIN=MAX rather thanCOUNT(DISTINCT)=1; both ignore NULL and yield NULL for zero or multiple distinct knownBIGINT values. A20-occurrence comparator executes the entire old/new entity projection and compares complete rows in eight entity/identifier tables/views, including empty, invalid and overlapping source cases; upstream invalid-taxonomy validation is unchanged.25focused tests passed in3.01s. This eliminates unnecessary distinct state in that aggregate without changing any Parquet, identity, main schema or scientific result.

Independent coverage was completed for structured measurement identities, actual Parquet-to-COPY ontology terms/axioms and frozen-main term/shortest-closure outputs, plus exact main preset registry descriptors/DDL/defaults/upserts. Main presets register metadata only; resource/organism/license filtering and folds remain the separate API consumer's responsibility, outside this build oracle. The pinned release already usesconnectomedb2025 andmetatlas. Root ran all seven new checks with `OMNIPATH_TEST_POSTGRES=1 uv run --frozen --no-sync pytest -q packages/omnipath_postgres/tests/test_main_network_presets_parity.py packages/omnipath_postgres/tests/test_main_parquet_semantics.py packages/omnipath_postgres/tests/test_aligned_entity_memory.py`:7passed in2.38s on disposable local PostgreSQL16; Ruff clean.

Next retry is configured for DuckDB1thread/4GB and a7GiB builder cap, retaining the same2CPU/Nice10, exact pinned release and durable base/product checkpoints. The standalone helper now records its own hash and private traceback filenames/functions/line numbers for precise future failure diagnosis. No mandatory verification/rebuild/rollback phase was added; full completion and actual full-scale memory/storage/timings remain pending.


## 2026-10-01: full attempt4 started with revised memory settings

Full attempt4 started at2026-10-01T21:43:49.757124+00:00 as `omnipath-migration-main-aligned-attempt4-20261001.service`. Runtime commit `5703acf50630e8c1716b0b183c0510087200b66f` is pushed;742committed source files (including the completed independent fixtures and workspace package files) were verified, manifest SHA256 `74abecfd60313583ff2414cac4ddfc201cee26e22feb70ed769c3807de15bd4c`. Standalone run helper SHA256 `44e773253b3bc7db11e2a2e56408edc2398396c95dd594bc301ac583862a8146`.

Parameters: unchanged release2026.9.30.4/46resources/canonical SHA25606f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e; fresh `aligned_full_20261001` in isolated mainPG18loopback5441; DuckDB1thread/4GB, on-disk staging, builderMemoryMax7GiB/CPU2/Nice10; exact-main late keys/FKs/indexes and complete base/product checkpoints. Admission confirmed previous unit/backend terminal, fresh target, source hashes, private image/loopback, public extensions, both original pilots and protected services. Free disk538,278,174,720bytes; hostMemAvailable12,998,672KiB. Initial21:45UTCprogress wasvalidate_and_stage, active unit, matching pins and no base yet.

The existing15-minute overnight heartbeat now follows only attempt4. No mandatory verification/rebuild/rollback/safety phase is queued. Prior failed/cancelled empty schemas and diagnostics are preserved. Full outcome, timings, representative samples and final storage remain pending.


## 2026-10-01 — Main alignment full attempt4: annotation CTAS memory failure

- Runtime: `5703acf50630e8c1716b0b183c0510087200b66f`; unit `omnipath-migration-main-aligned-attempt4-20261001.service`; private PostgreSQL18 at loopback5441; release `2026.9.30.4`, manifest `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`. DuckDB1thread/4GB, builder7GiB, base/product checkpoints enabled, no mandatory verifier.
- Started `2026-10-01T21:43:49.757124+00:00`; failed `2026-10-01T22:39:50.912637+00:00`; wall **3361.156325s**. Canonical/entity/statement/dictionary preparation passed. Prepared staged links included42,991,511 entity identifiers,42,607,573 source evidence identifiers,8,562,132 entity annotation links and22,887,251 relation evidence rows. These were not PostgreSQL COPY results.
- Exact failed CTAS: `relation_evidence_annotation`, `prepare_aligned_release` line705. DuckDB OOM requested1GiB at3.4/3.7GiB used. No COPY began and no checkpoint was committed, so a new load is required; there is no copied base to repeat or lose.
- Archived reports/logs/runtime manifest under `/root/projects/omnipath-migration/main-alignment-20261001/attempt4-5703acf-annotation-memory`; run report SHA256 `bb02c438281d7bcff0b20ef366132ff0f6559ab8ad6a1f637dc756fd3312c27a`. Proven-empty schema retained as `aligned_full_20261001_attempt4_empty`. Both original pilots, both completed aligned pilots and all four protected service IDs/start times were confirmed intact. Temporary projection/spill cleanup restored free disk to538,317,664,256bytes.
- Next implementation: narrow primitive ownership inputs clustered by complete-owner shard; preserve existing source/dataset/scope and evidence-versus-statement rules. Independent full-row old-query regression and review must finish before launching the next attempt. No Parquet, reference, cache, resolver or main scientific algorithm changes.


## 2026-10-01 — Bounded owner-aware annotation preparation fix

- Two narrow DuckDB link inputs carry only primitive resource/version/owner, evidence ordinal/UUID, annotation UUID, raw source/dataset, normalized scope and typed/synthetic flags. They are physically ordered by a stored 128-group owner hash; each group stays complete and filter predicates target the stored group column. Values, quantities, arrays and upstream record fields stay out of the ownership joins.
- Separate evidence-owned and surviving statement-owned paths preserve the frozen query's exact null-safe suppression, raw wildcard matching, scope normalization, synthetic evidence and final four-column DISTINCT. The staged result is reused for counts/COPY. No main schema/index/scientific algorithm or source artifact changed.
- Independent full-row/decoded-expectation fixture: 13 passed in3.77s; combined projection/quantity suite:38 passed in8.58s over nine synthetic records. Loader/checkpoint/product adapter suite:60 passed in1.48s, including real DuckDB CSV and disposable PostgreSQL checkpoint visibility. Ruff passed; independent review found no remaining row-semantic issue. Progress field `event` initially collided with real loader `_emit`; renamed to `state` and tested through that exact callback.
- Exact full-release memory behavior is pending the next load. An exceptionally large complete owner cannot be split without a separate reviewed semantic approach; the fix does not claim a fixed bound independent of owner skew.


## 2026-10-01 — Main alignment full attempt5 launched

- Runtime `afd4ccb9b8f722d5dee885d567a7f378a9a1e88e`, reviewed 128-group narrow annotation-link fix; unit `omnipath-migration-main-aligned-attempt5-20261001.service`, invocation `279b6708d46c450d98ef3d105f4c5ac9`. Started **2026-10-01T22:58:29.704298+00:00**.
- Fresh schema `aligned_full_20261001` on private PostgreSQL18 loopback5441; unchanged46-resource release2026.9.30.4, manifest `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`. No artifact, cache, reference or resolver changes.
- Verified743committed runtime files, source manifest SHA256 `796cfb701cc383eecb42fc0cd5ab688a01849d510b0f3e23067c923890311c19`; archive SHA256 `31f88ba48b6adb2b21b73a2a3a72d0332fb4c21bd5e4fb6b7f3520815be61354`. Helper SHA256 `44e773253b3bc7db11e2a2e56408edc2398396c95dd594bc301ac583862a8146`.
- DuckDB1thread/4GB; builderMemoryMax7GiB/CPUQuota200%/Nice10. Base and complete-product checkpoints enabled; no mandatory verifier/rebuild/rollback job. Preflight free538,317,246,464bytes, hostMemAvailable12,989,360KiB; exact protected service IDs/start times and both original/completed aligned pilot pairs preserved.
- At22:58:44UTC authoritative service state was active/running in validate_and_stage. Actual completion/phase times/storage and representative manual samples remain pending. Monitor target and15-minute heartbeat updated to this unit/revision.


## 2026-10-02 — Main alignment attempt5: annotation batches passed, ontology staging failed

- Runtime `afd4ccb9b8f722d5dee885d567a7f378a9a1e88e`; unit `omnipath-migration-main-aligned-attempt5-20261001.service`; private PostgreSQL18 loopback5441; unchanged release2026.9.30.4/46resources/manifest `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`. DuckDB1thread/4GB, builder7GiB/2CPU/Nice10; complete base/product checkpoints enabled; no mandatory verification/rebuild/rollback.
- Started `2026-10-01T22:58:29.704298+00:00`; failed `2026-10-01T23:55:16.280352+00:00`; wall **3406.576844s**. All128 annotation groups finished, staging **78,521,399** relation-evidence annotation links. Annotation preparation ran from23:48:10.225314 to23:54:44.863942UTC (**394.638628s**): annotation input228.658486s; evidence input36.269404s; grouped publication approximately129.70s including coordinator overhead. These are DuckDB staging results, not copied PostgreSQL data.
- The later `ontology_terms` CTAS failed at `prepare_aligned_release` line802 with DuckDB OutOfMemoryException (4KiB allocation at3.7/3.7GiB used). Its query aggregates synonyms for all owners and joins all annotations before filtering descriptions. No COPY began and no base/product checkpoint exists.
- Archived private reports/log/runtime receipts in `/root/projects/omnipath-migration/main-alignment-20261001/attempt5-afd4ccb-ontology-memory`; run report SHA256 `6fd717f0ce2c8639a79f18948fcee7fa02af11b04ffa7832d014ed3e434cf8cf`. The terminal unit hasMainPID0; actual empty biological tables and no builder backend were checked before retaining setup as `aligned_full_20261001_attempt5_empty`. Both original pilots, both completed aligned pilots and all four protected service IDs/start times were preserved. Temporary staging cleanup restored free disk to538,314,575,872bytes.
- A review also found that the publisher emits `statement_kind='ontology'`, while the adapter's ontology COPY/capability gates expect only `ontology_axiom`; this explains the staged ontology-edge count0 and requires a downstream adapter correction. Published files remain unchanged. Narrow owner-scoped ontology preparation, complete-row old-query/independent fixtures and actual publisher-kind coverage are in progress before any retry.


## 2026-10-02 — Main-aligned ontology staging and serving compatibility

Attempt5 failed before COPY at the invented wide ontology_terms enrichment. Tracing frozen main 9f9bb709c764 showed that its eight-writer projector never populates ontology_terms_raw and records zero terms on a fresh build. Keep the exact empty legacy base and use main's existing entity_ontology_term derivation from relational ontology edges, identifiers and annotations. No Parquet, parser, reference, resolver, physical main table/index contract or scientific projection changed.

A read-only aggregate audit of all46 pinned relation Parquets found987,393 published ontology statements and zero ontology_axiom. Accept the actual publisher kind (plus legacy compatibility) and recover the twelve known input ontology families independently of post-resolution namespaces. Unknown provider scope remains an explicit capability gap. Narrow primitive endpoint lookups add exact-string controlled-vocabulary aliases only for eligible source/type/namespace owners; genuine scalar labels can add Name lookups only when no meaningful global Name/Synonym exists and the label is not an identifier fallback. Generated aliases have no fabricated occurrence evidence or authority. The serving reader uses one materialized entity/value proof and retains main's canonical-ID priority, original definition evidence and chemical/protein selection rules.

Final exact frozen development suite:131passed,11explicit private-DSN skips in17.54s across14files; the ten ontology checks passed in8.01s with disposable local PostgreSQL16. Coverage includes actual SilverExtractor/ParquetWriter output, empty base COPY, frozen-main complete serving row/definition parity, multiple roots without row multiplication, negative source/type/namespace and original CV-alias guards, global meaningful-name eligibility and unchanged raw-taxonomy counts. Ruff passed all five changed Python files. Independent source, plan and full-row review found no remaining mismatch. Earlier private PostgreSQL18 catalogue/oracle tests cover the unchanged physical DDL. These bounded fixtures do not establish full-release completion or full-scale memory behavior; attempt6 will test that. Attempt5's empty schema and diagnostics are retained.


## 2026-10-02 — Main-aligned full release attempt6 launched

Started at2026-10-02T00:51:32.984886+00:00, unit omnipath-migration-main-aligned-attempt6-20261001.service, invocation685ae471d47a4cac984b55454cb5e5ea. Exact runtime commit8ae7e8ba12113232e2484cde596847cf0bf8f7b9 was pushed on parquet-migration;744committed files were checked against manifestSHA88b8d4bdaae2b16550ef0854c9f856ba7a85e9fb8233d5f27d424e0465c45f4f and archiveSHA65152f7f8919de2b481c0097e6b7dced4f7ec21d6aa0805058f0dd5fee5dc71b. HelperSHA44e773253b3bc7db11e2a2e56408edc2398396c95dd594bc301ac583862a8146. Attempt5 was authoritatively terminal and its actual empty schema retained before launch; no imported or committed work was discarded.

Destination is private PostgreSQL18 loopback5441/schema aligned_full_20261001, unchanged46-resource release2026.9.30.4/canonicalSHA06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e. DuckDB1thread/4GB; builder7GiB/2CPU, base/complete-product checkpoints. Preflight free disk538313969664bytes, hostMemAvailable13037680KiB. Both original pilots, both completed aligned pilots and all four protected service IDs/start times matched. At00:51:41UTC unit active in validate_and_stage, release preparing, no committed phases. The15-minute automation now names attempt6 and exact runtime provenance. No mandatory verifier, rebuild/rollback or safety job exists. This is a running launch, not successful completion; actual staging/COPY/index/shared/product durations and optional bounded final samples remain to be recorded.


## 2026-10-02 — Attempt6 staged all20 tables and entered COPY

All20 table queries were prepared by2026-10-02T01:59:52.702000+00:00 at frozen runtime8ae7e8ba12113232e2484cde596847cf0bf8f7b9. All128 annotation groups completed:78,521,399links in381.136341s. The ontology lookup stage took18.287307s; ontology edges987,393 exactly cover the published axioms, and the legacy ontology_terms base remains empty as in fresh main. Prepared quantities1,577,215; identifier/annotation occurrences80,136,865/150,690,040. Dictionary/link lookups include the declared computed CV/name aliases; original evidence links and canonical/source entity counts are unchanged. This passed the prior staging failure point on the full release without Parquet or runtime changes.

At02:07:02UTC eleven COPY tables had completed: entity, identifier_evidence, annotation, entity_evidence, entity_evidence_resolution, entity_identifier, entity_evidence_identifier, entity_evidence_annotation, relation, relation_evidence and relation_evidence_relation. Actual COPY alone for the two large identifier-link tables took60.542292/55.133320s; their CSV staging took12.389521/11.657077s. COPY is still one unfinished transaction until restored keys/FKs; no durable base or product was claimed. Remaining COPY, keys/FKs, main query indexes, base checkpoint, main shared derivations and products are pending.

Independent read-only status audit: no active main_compat reader unintentionally excludes published5; shared status_id<>2 admits it, and all-five-source MetSigDB parity explicitly uses5 under the agreed policy. Retired matched-only MetaLinks templates are bypassed by current presets, matching frozen main; external API consumers were outside this audit. Existing immutable BRENDA audit covers7,337endpoint owners, all with genuine CV/CURIE candidates, so no additional scan or alias expansion was needed. No additional tests, deployments or server writes resulted from these audits.


## 2026-10-02 — Optional published inspection copies and exact deferred indexes

- Reviewed default fresh-load policy: 14 main COPY tables plus required `annotation_quantity`; five `parquet_*` inspection copies are opt-in. No Parquet, reference, cache, resource build or resolution changes. Scientific ownership, annotations, scopes and measurement fields remain normalized; original inspection keys/order/repetition stay in pinned Parquets. Existing twenty-table checkpoints and counts remain compatible; no existing tables were removed.
- The first real disposable PostgreSQL16 retention test exposed a pre-existing late-index defect: captured `pg_get_indexdef` includes `ON ONLY`, recreating an invalid molecular-type parent without source/default child indexes. Initial focused suite: 99 passed, 1 failed in16.67s; data/serving/finish assertions passed, exact catalogue assertion failed.
- Fixed capture/replay only for catalog-confirmed partitioned index roots, preserving original saved triples. A narrow additive molecular-type repair attaches valid equivalent existing child indexes or creates missing ones, preserves parent/child OIDs and rows, and verifies final parent validity. No copied main science or DDL definitions changed. Independent review found no remaining issue on the corrected actual path.
- Final command: `OMNIPATH_TEST_POSTGRES=1 uv run pytest -q --tb=short packages/omnipath_postgres/tests/test_aligned_index_replay.py packages/omnipath_postgres/tests/test_aligned_provenance_retention.py packages/omnipath_postgres/tests/test_aligned_projection.py packages/omnipath_postgres/tests/test_aligned_loader.py packages/omnipath_postgres/tests/test_aligned_ontology_memory.py packages/omnipath_postgres/tests/test_main_parquet_semantics.py`. Result: **104 passed in16.54s**; Ruff clean. Eighteen synthetic entity/statement records exercise full normalized DuckDB/PostgreSQL rows, nonempty serving derivation, quantities/scopes, exact frozen-main parent/source/default catalogue in both modes and legacy checkpoint finish. Four additional replay/repair cases retain rows and completed index OIDs. No server tests or data writes for this change.
- Existing full attempt6 remains frozen at8ae7e8b, PID396212. All20 COPY operations completed2026-10-02T02:17:53.749151UTC: summed COPY **803.379790s**, CSV staging **264.891105s**, 565,375,453 rows and81,215,151,078 CSV bytes. At02:55:13UTC the same unit/PID was active in constraints, release preparing, no advertised base; free disk335,349,739,520bytes. A bounded activity diagnosis found foreign-key validation with zero blocking sessions. Optional storage observation times out under constraint locks and is tolerated; this is not a terminal build error.
- Pending: let constraints/indexes and main derivations/products finish; apply the narrow missing-index repair to the current private target only after a durable checkpoint/completion and schema-lock/release guards, without COPY or science repetition. Inspect reasonable final samples and actual timings/storage. New retention defaults do not change this twenty-table run or delete its inspection copies. Fifteen-minute monitor remains active; no mandatory verifier/rebuild/rollback/safety job.


## 2026-10-02 — Full attempt6 identifier-index failure and grounded diagnosis

- Full attempt6 stopped at2026-10-02T03:01:05.206952UTC after **7772.222781s**, SQLSTATE54000: `identifier_evidence_type_value_idx` row3640bytes exceeds main B-tree limit2704bytes. All20 COPYs had completed (803.379790s; CSV264.891105s), but captured index restoration ran in the same uncommitted transaction. The error rolled back COPY and constraints; no base/phase/product was committed. At03:04UTC exact EXISTS checks found no visible entity, identifier_evidence or annotation_quantity rows. Schema setup/dimensions and saved plan remain. Both completed aligned pilots and all four protected service IDs/start times matched baseline. No source/artifact/service changes.
- Immutable failure archive: `/root/projects/omnipath-migration/main-alignment-20261001/attempt6-8ae7e8b-identifier-index-size`, created03:17:39UTC. Run report SHA2562923b73bce06d8e0b754e4f6481530d0cdda1a0e8da59b0648ac040ac4fd5944. Private error frames were from the older attempt5; the sanitized current run report is authoritative. The next helper records safe current frames for PostgreSQL exceptions too. No schema rename/drop has yet occurred.
- Read-only pinned-release identifier-size diagnosis: **169.39s**,46resources. Raw occurrences>1800UTF8bytes: InChI4096(max4038), SMILES300(max3723), synonyms82(max2837). No canonical identifiers or scalar labels>1800bytes. Main `ingest.common.include_identifier` explicitly excludes Standard InChI from generic IDs; the adapter had included it. Main has no general length cap, and its raw unique plus lower-value B-trees must both remain unchanged.
- Narrow actual-value probe: **86.34s**, private PostgreSQL18 temporary session tables only, all2114distinct long values. InChI1923values caused96raw-index and99lower-index failures54000. All150distinct SMILES and41synonyms passed both exact main index definitions, including UUID INCLUDE on lower(value). No source values were printed or persisted in reports; no existing schema/product/metadata changed. Probe report SHA256dc0bb2004a8b5c988e87425e2c7b50fc9de9a64e5bdfd744973c0cb9b55e8313.
- Narrow reviewed correction in progress: main's exact Standard InChI exclusion from generic dictionary/links, with complete raw occurrences retained in opt-in inspection copies and pinned Parquets; raw SMILES gets exact `Smiles:MI:0239` spelling required by main MetSigDB. Canonical identities, all accepted long identifiers, scientific projection and exact main indexes remain unchanged. Independent code review found no remaining issue; focused tests and one guarded default15-table restart are pending. No general chemical validity filter is added: old raw normalizers reject published outer CURIE spellings that the documented lookup aliases legitimately support. Future changes require an explicit occurrence/canonical policy, rather than silently rejecting authoritative published identities.


## 2026-10-02 — Exact-main InChI/SMILES correction ready for guarded retry

- Frozen correction excludes only raw `inchi` and exact `Standard Inchi:MI:2010` from generic ID dictionaries/normalized links, matching main. Opt-in audit rows retain full namespace/value/order/source/duplicate/null occurrences with NULL generic dictionary refs. Canonical UUIDs/values and scientific annotations/quantities are unchanged. Raw SMILES maps to exact `Smiles:MI:0239`; no length cap, truncation, broad chemical filter, schema or index exception.
- Independent reviewer found no remaining issue. New five-record synthetic fixture executes SHA-verified frozen inclusion policy and compares full default/opt-in dictionary/occurrence/quantity outputs; preserves >4KB accepted non-InChI values. Actual local PostgreSQL16 COPY plus complete frozen/current MetSigDB publication rows/attributes retain typed SMILES, using the existing approved published-status eligibility adaptation.
- New-file command `OMNIPATH_TEST_POSTGRES=1 uv run pytest -q --tb=short packages/omnipath_postgres/tests/test_aligned_identifier_policy.py`: **3 passed in3.25s**. Final command adds that file to the previous six-file replay/retention/projection/loader/ontology/semantics suite: **107 passed in17.60s**, disposable local PostgreSQL16 only; Ruff clean. Early test-only fixtures were corrected for an incompressible unknown-namespace value outside the real retained data and the already approved status1/status5 eligibility comparison; no production cap or status change was introduced.
- Next build will use the new default15-table policy. Reviewed standalone build helper explicitly disables PostgreSQL inspection copies and records safe current PostgreSQL exception frames (SHA25655628f1ed6dfec94fc3c46022dd8a52aae6c0722ce39a5d82e2a61da88757ab0). The optional bounded manual sample now checks normalized relation-evidence endpoints, so it supports default loads (SHA2569640c55fa7b4971a4ee3309ca44c2999cbb706ac778d4e1e90c07b150e2097c8). Both syntax checked; neither was launched yet. Old failure/schema/diagnostics and all artifacts/pilots/services remain preserved. Guarded source packaging/deployment and full attempt7 launch are pending.


## 2026-10-02 — Guarded full attempt7 launched with exact main identifier policy

- Reviewed runtime **7bb2ce7dc9c85bd0ebf6ab26a841470181887ca4** was committed/pushed and staged only after attempt6 stopped. All **747** committed files matched manifestSHA256723d92acad97e372e1d01bdb51eb74f792e703b016dcfd612fdda13fdd34c5ce; archiveSHA256399e49d7624705267fa2a0083c6ed7ba6df1f66fe5c247aad29d904b867908ca. Build helperSHA25655628f1ed6dfec94fc3c46022dd8a52aae6c0722ce39a5d82e2a61da88757ab0; manual sampleSHA2569640c55fa7b4971a4ee3309ca44c2999cbb706ac778d4e1e90c07b150e2097c8.
- Preflight confirmed the old importer stopped, no other private database work, exact public extension placement, both original and both completed aligned pilots, four protected service IDs/start times,46unchanged2026.9.30.4pins, private PostgreSQL18 image/loopback5441 and sufficient resources: free **487,605,088,256bytes**, MemAvailable **13,176,272KiB**. Under the schema advisory lock, renamed only the proven rolled-back failed target to **aligned_full_20261001_attempt6_empty**; preserved setup/dimensions/plan/pages and failure archive. No schema was dropped or source artifact changed.
- Started **omnipath-migration-main-aligned-attempt7-20261001.service**, invocation **a742b50adbe4466f8ec42058918e4744**, at **2026-10-02T03:29:57.210086UTC**. Destination **aligned_full_20261001**, same release; DuckDB1thread/4GB, builder7GiB/2CPU, complete base/product checkpoints enabled, optional audit copies explicitly disabled (15COPYtables). No mandatory verifier/rebuild/rollback/safety job.
- Authoritative initial check at03:30:21UTC: **PID509418 active/running** in validate_and_stage, expected7bb2runtime/helper/pins, releasepreparing, no committed phases. PostgreSQL idle while DuckDB stages is normal. Free disk487,602,835,456bytes. ObserverSHA256284abdc779ce9247a263152f54366dcc4be72c6707ac67d4277aa5610f873cf1 tolerates optional storage lock timeouts; the15-minute heartbeat now follows only attempt7. Full staging/COPY/indexes/base/shared derivations/products and reasonable final samples remain pending; no completion claim.


## 2026-10-02 — Main-aligned attempt7 COPY completed; constraints pending

The authorized isolated nicesrv run at runtime 7bb2ce7dc9c85bd0ebf6ab26a841470181887ca4 completed all 15 COPY operations at 04:30:21.884362UTC. Bounded event-log extraction measured 277,233,095 rows and 24,522,976,857 CSV bytes; summed COPY was 309.928782 seconds and CSV staging 83.141180 seconds. Initial DuckDB projection is separate. Core counts remain 4,794,341 entities, 6,383,473 source entity-evidence rows, 78,521,399 relation annotation links, 987,393 ontology axioms and 1,577,215 structured quantities. Generic identifiers are 26,213,189 after the exact main Standard InChI exclusion; fresh-main ontology_terms remains empty. The five optional inspection copies were not loaded, and source artifacts/canonical identities are unchanged.

At 04:35:24UTC, unit omnipath-migration-main-aligned-attempt7-20261001.service and original PID509418 were active restoring constraints/indexes, with no reported error or durable base checkpoint. COPY plus captured constraints/indexes is still uncommitted. No completion claim, mandatory verifier, rebuild or rollback job was added. Shared derivations, durable products, bounded final samples and final storage/timings remain pending.


## 2026-10-02 — Attempt7 committed the main-aligned base

The private base checkpoint completed_at is 05:14:09.308667UTC; shared derivations started at 05:14:15.514200UTC. A separate read-only connection confirmed committed imported rows at 05:10, followed by status base, all 15 loaded counts and the durable phase at 05:23. COPY is retained; any later retry must use PostgreSQL-only finish and preserve completed phases/products.

Recorded timings: artifact validation 6.563978s; initial DuckDB projection 3220.882356s; COPY 309.928782s; CSV staging 83.141180s; captured constraint/index restoration 2263.296688s. Base phase 3022.195104s starts after projection and includes COPY/CSV/restoration plus remaining query indexes, analysis and checkpoint work. Do not sum its components into it again. Bounded catalogue reads found 1,870 indexes, all valid/ready, and 7,424 constraints, all validated. Exact comparisons of all 114 captured constraint definitions and six captured index definitions found zero missing/mismatched objects; the molecular-type parent had 47 attached valid/ready children. Shared derivations and all complete products remain pending; no full-release completion or final storage claim. No mandatory verifier/rebuild/rollback/safety job was introduced.


## 2026-10-02 — Attempt7 shared derivations and first two products committed

- Frozen runtime remains `7bb2ce7dc9c85bd0ebf6ab26a841470181887ca4`; original unit/PID `omnipath-migration-main-aligned-attempt7-20261001.service`/509418, private PostgreSQL18 loopback5441, schema `aligned_full_20261001`, unchanged 46-resource release `2026.9.30.4` and manifest `06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e`. No new deployment or repeat COPY.
- Shared `derived` phase took **7222.042703s**, with actual after-commit observer event **2026-10-02T07:14:37.561622UTC**. Its component times are: identifier_authority0.008251s; derived_tables1277.045308s; chemical_resolution_levels445.912745s; chemical_ambiguous_name_candidates118.588081s; chemical_class755.962660s; metabolic_domain894.426454s; interaction_class0.019227s; interactions1668.977273s; entity_labels413.442918s; entity_name198.410665s; chemical_labels499.640786s; bitmaps947.903007s; resource_overlap0.307190s; resources0.587652s; licenses0.046916s; build_manifest0.704750s. Components are included in the derived total.
- MetSigDB took **350.714549s** and committed at **07:20:28.333959UTC**. PostgreSQL records **3,598,709 memberships** and five nonempty sources: Reactome24,301rows/2,272sets/2,178metabolites/20.3s; WikiPathways11,513/963/4,101/13.8s; KEGG92,179/806/6,177/67.9s; MACdb20,291/269/5,389/22.8s; ClassyFire3,450,425/2,863/145,793/224.5s. These are rounded per-source measurements; the complete product time is authoritative.
- Main network preset registry took **0.023042s** and committed at **07:20:28.364609UTC**. MetSigDB and network_views subset metadata independently match the exact release and manifest. This is registry/build parity, not a claim about the separate API's serving filters/folding.
- COSMOS started **07:20:28.366038UTC**. At **07:31:40UTC** the same PID was active without a terminal error; release status was network_views and all prior phases/products were durably committed. One backend was reading a temporary spill file. Storage snapshot **97,908,023,296bytes**, free disk **439,253,340,160bytes**; neither is final. Only COSMOS and the bounded final samples/report remain. A later failure must resume only COSMOS through PostgreSQL-only finish, preserving every completed phase/product.
- Phase existence proves durability. `parquet_phase.completed_at` defaults to PostgreSQL `now()` and may represent transaction start; actual commit chronology above comes from after-commit `phase_committed` observer events. Prepared/source-done events are not commits. The 15-minute monitor now follows COSMOS only; no mandatory verification, rebuild/rollback or safety job is queued.


## 2026-10-02 — Main-aligned full attempt7 completed and sampled

- The unchanged runtime `7bb2ce7dc9c85bd0ebf6ab26a841470181887ca4` completed at **07:48:10.874505UTC**; the build_complete observer event followed at **07:48:10.900235UTC**. Unit `omnipath-migration-main-aligned-attempt7-20261001.service` exited0/success, MainPID0, no builder backend. Wall **15493.664462s (4h18m14s)**. Release2026.9.30.4/manifest06f55112b371055b9b539c27eb7886d403b0131098393d1bafe4723ce006c03e is complete in aligned_full_20261001/private PostgreSQL18 loopback5441. All base/derived/metsigdb/network_views/cosmos phases and product metadata are committed; no restart or server code deployment was needed.
- Actual COSMOS phase **1662.485741s**; after-commit event **07:48:10.859545UTC**. It published2,432,572reaction edges +530,195connectors;47,292reactions,20,033reversible,31,158orphan. Earlier measured projection/base/shared/MetSigDB/registry phases remain as recorded; COPY/CSV/restoration are included in base and must not be double-counted. Result report SHA2568faad46672607d9c133ed00b87e54396a7cfbe604f45994dd10ecd7f703e460e.
- Optional manual sample at07:52:08UTC:26aggregate results, no query notes; final index count2,317 valid/ready, constraints7,548 all validated. Protein/chemical/activity, graph/source endpoints, ontology/interaction facts,100reaction parties,100quantities, source-owned KEGG/metatlas/Reactome annotations,100memberships from each of five MetSigDB sources,100COSMOS references and exact preset registry checked. Source-only type-column scan0.471194s confirmed6,383,473published entity rows usechemical_entity rather than small_molecule, explaining the empty category probe. Additional100-row reverse/orphan/transport samples passed at08:03:44UTC (0.007176/0.495890/0.001348s), including forward counterpart/reversibility, pseudo-enzyme/canonical references and compartments. Sample file SHA2562f44ddc6758cdae6cc57513b329ca2c32632c5caabb1ce8dc8bd8c30a4cf968a.
- Bounded final preservation read07:57:18UTC confirmed46exact version/build-manifest pins, canonical release hash, all four protected service identities/start times/running state, original pilot_20260930 and annotation_sample_20260930, both completed aligned pilots and all six retained failed setup schemas. At08:08:32UTC schema725physical tables/98,933,350,400bytes (92.14GiB), free438,284,460,032bytes. No further cleanup, source build, resolver run, cache/reference change or mandatory verifier/rebuild/rollback/safety phase.
- Read-only COSMOS diagnosis while the original run remained active found the exact adapter's OR endpoint join planned as Nested Loop over materialized claims. The run completed before any intervention. Separate future SQL optimization materializes the identical owner-qualified three-column membership lookup once and uses disjoint subject/object equality branches with IS DISTINCT FROM for NULL/self-loop overlap. Local PG16 focused fixture:3passed1.31s, complete claim multiset/independent expectations/full rows/statistics and bounded insertion order; main-science regression file2passed4.70s; Ruff and independent review passed. Final planner-only branch estimated41rows/cost4,522,288.38, oneproducer/twoscans. No full-release speed claim or global order guarantee. Production schema/indexes/science are unchanged; completed server runtime remains7bb. Diagnostic scaffolding errors (Decimal JSON serialization, local prefix extraction syntax and missing resources path component) were corrected without affecting the build.
- Detailed aggregate evidence saved to docs/reports/nicesrv-main-aligned-completion-20261002.json; report updated with actual successful timings and retained attempt history. Prior attempts1–6 ended before any durable base, including attempt6's rolled-back COPY/index transaction. All are retained; successful attempt7 base/products are preserved. Completion uses reasonable samples, not an exhaustive verification/rebuild.


## 2026-10-03: package consolidation and developer workflow

Implemented short project folders and flat Python packages, while preserving public distribution/import names. Runtime matching belongs to resolver; PostgreSQL owns projection/shared derivation; subsets owns scientific products and one checkpoint runner. Production entry points use the current normalized layout, with historical record-layout APIs explicitly isolated. Shared publication contracts, installed-file provenance, reproducible vocabulary generation, bounded API queries, role-specific packaging and source-only developer checks are in place.

Validation: combined Python run871 passed/32 skipped/323 deselected; complete disposable PostgreSQL/subsets suite747 passed/11 skipped; final API147 passed, plus2 worker/publisher environment checks; native29 and web26 tests passed. Clean source-oracle and isolated wheel checks pass; all12 product SQL files and scientific algorithms remain unchanged. Full details and intermediate findings are in [the implementation record](reviews/refactoring-progress-20261003.md). No live service, published data or server database changed.

## 2026-10-03 — Remove historical subset compatibility

Removed `omnipath_subsets.compatibility`, its optional driver extra, the obsolete
full-release verification/rebuild script, and their dedicated tests. Current
MetSigDB, network views, COSMOS, SQL and checkpoint runner are unchanged. Retained
current scientific parity and transaction tests; moved publication guard coverage
to the current runner tests and detached the older PostgreSQL identifier-index
regression from the removed audit script. Updated package and architecture docs.

Validation: 601 passed, 11 existing opt-in skips in 43.59 seconds against disposable
PostgreSQL plus the API migration-serving fixture. All 1,083 remaining tests collect;
isolated subsets wheel, Ruff, formatting and lock checks pass. No deployment,
resource build, server access or persistent database changes.

## 2026-10-03 — Refresh nicesrv serving with API resource limits

Deployed API/web from commit `2b04935098e256880e83eed33eb27c6e5f82dfab` in the
archived source `/root/projects/omnipath-releases/20261003-serving-2b04935/source`.
Project `omnipath-serving-20261003-2b04935` uses loopback ports 8285/8282/8280 and
reuses `/root/projects/full_parquet/data` read-only. Container budget: API 4 GiB,
4 CPUs; query pool two databases, each 1 GB/two threads, 30-second admission wait.
The actual installed query settings and Docker limits were checked.

Candidate tests compared complete small entity/relation results from SIGNOR,
Rhea and ChEBI, relation evidence and one filtered CSV export with the preceding
API. All matched. Eight requests from four concurrent clients passed. All 46
Latest pins and release selections were preserved; explorer, web proxy, HTTPS
health and byte-range file reads passed. The first routing probe ran while the
new file container was still becoming healthy; it correctly did not confirm the
new file route. After all services were healthy, fresh marked public requests
were confirmed in the new API/file container logs before stopping the old stack.

The three previous serving containers were stopped and retained for rollback.
All 30 other pre-existing containers retained their IDs, images, start times and
states. All new services were healthy with no restarts or OOM events; API memory
was 232.5 MiB after smoke checks. Public checks passed again after old services
stopped. No resource/reference/PostgreSQL build, data publication or worker start.
This is bounded functional validation, not a full-load performance benchmark.
See `docs/reports/nicesrv-serving-limits-20261003.json` and `deploy/nicesrv.md`.
