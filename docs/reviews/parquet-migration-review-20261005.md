# `parquet-migration` branch review

Date: 5 October 2026
Branches:
- omnipath-build `parquet-migration` at `8dbe0d3` (68 commits / 890 files ahead of `main`)
- pypath submodule `parquet-migration` at `dd23d862b` (20 commits / 97 files ahead of `master`)

Status: review complete; fixes applied 5 October 2026 (see Status column and "Follow-up" at the end)

## Scope and approach

Correctness-focused code review of both branches. The architecture/package
layout was already reviewed in
[architecture-review-20261003.md](architecture-review-20261003.md) and
implemented per [refactoring-progress-20261003.md](refactoring-progress-20261003.md);
this review does not repeat those findings. Generated files (OpenAPI, Biolink
vocab, shadcn UI components, lockfiles) and run logs/reports are not reviewed
line by line.

Severity: **High** = wrong data/results, security exposure or broken workflow;
**Medium** = latent bug or fragile behaviour likely to bite; **Low** = cleanup,
docs or minor inconsistency.

## Summary

| ID | Repo | Severity | Finding | Status |
| --- | --- | --- | --- | --- |
| P1 | pypath | High | `omnipath-core` dependency unresolvable outside omnipath-build; stale `uv.lock` path | Fixed: inputs_v2-only deps moved to an `inputs-v2` extra; uv source `../packages/core`; both locks regenerated |
| B1 | build | Medium | Native resolver `.so` not rebuilt after Rust changes → 41 test-collection errors | Fixed: `[tool.uv] cache-keys` for Rust sources |
| B2 | build | Medium | Capped/sample builds can be pinned into releases and loaded into PostgreSQL | Fixed: publisher and loader refuse sample builds unless `"partial_resources": true` |
| P3 | pypath | Medium | Mapper exceptions now abort whole builds (previously swallowed) | Mitigated: failures now name source/dataset/row; fail-fast kept, full-source sweep still advisable |
| P4 | pypath | Medium | 12 new pypath test files (183 tests) are not run by any target | Fixed: `make test-pypath` runs all 18 branch test files (327 tests) |
| P5 | pypath | Medium | Obsolete ontology terms silently dropped | Documented as intended policy in `BIOLINK_MODELING.md` (rule 7) |
| B9 | build | Medium | Job-queue first-use race (`no such table: jobs`) fails `/admin/status` | Fixed: readers treat a missing `jobs` table as an empty queue; regression test |
| B3 | build | Low | Release identity pins names/versions, not content digests | Fixed: releases record `resource_manifests` digests; loader verifies them |
| B4 | build | Low | Admin secret kept in `localStorage`; unused cookie auth path | Fixed: `sessionStorage` only (legacy copy removed); cookie auth removed in API and proxy |
| B5 | build | Low | Failed base load leaves a schema that neither `load` nor `finish` can recover | Fixed: self-created schema dropped on pre-`base` failure; writable spool fallback; integration test. Per-step derivation checkpoints not added |
| B6 | build | Low | Late cancel marks an already-published build as cancelled | Fixed: a normal return records `done` and its result |
| B7 | build | Low | Client defaults to a personal dev API host | Partly: `OMNIPATH_API_URL` override; dev-host fallback kept until an official endpoint exists |
| B8 | build | Low | `Claude.md`/docs reference removed make targets and old paths | Fixed |
| B10 | build | Low | ESLint error in `MolecularContext.svelte` (bare expression); `check-web` did not run lint | Fixed: `void entityKey`; lint added to `make check-web` |
| P2 | pypath | Low | `make test-pypath` red: stale BindingDB release in test | Fixed: test reads the pinned release default |
| P6 | pypath | Info | Many former annotation fields now raw-payload only | No change (documented policy) |
| P7 | pypath | Low | Formatting churn mixed with semantic changes | No change (needs history rewrite) |

## Automated checks

| Check | Result |
| --- | --- |
| `ruff check` / `ruff format --check` | pass |
| OpenAPI / Biolink generated-file checks | pass |
| `pnpm check` / web unit tests | pass (0 errors; 32/32 tests) |
| `make test` (before rebuilding resolver) | **41 collection errors** (B1) |
| Per-package runs after rebuild | core 82, client 43, subsets 34, resolver 1, postgres 373 (+31 skipped) passed; api 256 passed / **1 failed** (B9); build not completed (≈50 min, see test-suite section) |
| `make test-pypath` | **1 failed**, 143 passed (P2) |
| New pypath test files (run explicitly) | 183 passed |
| `uv lock --check` in `pypath/` | **fails** (P1) |

## Findings

### omnipath-build

#### B1 — Native resolver extension is not rebuilt after Rust changes (Medium, reproduced)

`packages/resolver` is an editable maturin package, but no
`[tool.uv] cache-keys` covers `src/**/*.rs` / `Cargo.*`. uv only rebuilds an
editable package when `pyproject.toml` changes, so after pulling `8ddc4c0`
(which added `resolve_molecular_batch` in `src/precomputed.rs`) a normal
`uv sync` / `uv run` keeps the Oct 3 `.so`. Result: `make test` fails at
collection with 41 `ImportError: cannot import name 'resolve_molecular_batch'`
errors in build/postgres tests.

Fix: add to `packages/resolver/pyproject.toml`

```toml
[tool.uv]
cache-keys = [
  { file = "pyproject.toml" }, { file = "Cargo.toml" }, { file = "Cargo.lock" },
  { file = "src/**/*.rs" }, { file = "rust/reference/Cargo.toml" }, { file = "rust/reference/src/**/*.rs" },
]
```

(workaround: `uv sync --reinstall-package omnipath-resolver`).

#### B2 — Releases can pin capped sample builds (Medium, confirmed by reading)

`build_resource` records `max_records` (and a `datasets` subset) in
`build_manifest.json`, and the API surfaces it as `sample_build`
(`packages/api/omnipath_api/queries/resources.py:193`). Nothing downstream
acts on it: `ReleaseStore.publish` (`packages/api/omnipath_api/store/inventory.py:96`)
only checks that files exist and are readable Parquet, and
`omnipath_postgres.releases._resource` validates checksums/schemas but never
looks at `max_records`/`datasets`. A `make build` result (default
`MAX_RECORDS=20`) can therefore be pinned into an OmniPath release and loaded
into PostgreSQL as if it were the complete resource, silently producing a
partial database.

Fix: reject (or require an explicit `allow_sample=True`) resources whose
manifest has a non-null `max_records` or a `datasets` restriction, in both
`ReleaseStore.publish` and `releases._resource`.

#### B3 — Release identity is name/version only, not content (Low, design)

The release JSON pins `{resource: version}`; checksums come from whatever
`build_manifest.json` is present at load time. PostgreSQL records the per-resource
manifest digest in `parquet_resource`, which is good, but `parquet_release.manifest_sha256`
(the identity used for checkpoint resume) is the digest of the name/version
JSON only. If a resource version is deleted and rebuilt under the same version
(possible on any data root where the release file isn't present, e.g. a remote
release manifest loaded from outside `data_root`), the same release identity
binds to different data and the checkpoint check in `_identity` won't notice.
Admin delete does refuse versions referenced by *local* releases, so this is
only a hazard across hosts. Consider pinning each resource's
`build_manifest.json` sha256 inside the release manifest.

#### B4 — Admin password persisted in `localStorage`; unused cookie auth path (Low)

`packages/web/src/lib/features/admin/auth.ts` stores the admin secret in
`localStorage` indefinitely (and actively migrates any `sessionStorage` value
into it), so it survives browser restarts and is readable by any script on the
origin. `require_admin` also accepts an `omnipath_admin` cookie containing the
raw secret, and the web proxy forwards it, but nothing sets that cookie. CORS
is `allow_origins=["*"]` without credentials, so the cookie path isn't
exploitable cross-site today, but it is dead code that would become a CSRF
surface if credentials were ever enabled. Prefer `sessionStorage` (or an
expiring token) and drop the cookie branch.

#### B5 — Failed base load leaves an unrecoverable schema (Low)

`load_release` creates the target schema, then COPYs. If anything fails before
the `base` checkpoint, the schema exists without a `base` phase: `load_release`
refuses ("requires a new schema") and `finish_release` refuses ("No committed
base checkpoint"). The operator must know to `DROP SCHEMA … CASCADE` manually.
Either document this in the CLI error message or drop the schema on failure
before `base` is committed. (Also: the default spool is created next to the
release manifest, i.e. inside `data_root/releases/main-staging`, which fails if
`releases/` is a read-only mount.)

Related: `run_main_derivations` is checkpointed as a single `derived` phase
(~17 sub-steps). A failure in a late step (e.g. `bitmaps`) re-runs every
earlier rebuild on `finish`. The steps are idempotent rebuilds, so this is a
cost rather than a correctness issue; per-step checkpoints would make long
full-release retries much cheaper.

#### B6 — A late cancel marks a successfully published build as "cancelled" (Low)

`AdminService.run_job` (`packages/api/omnipath_api/admin.py:516`) calls
`should_cancel()` *after* `_execute` returns. If a cancel arrives after
`build_resource` has already renamed the staged directory into
`resources/<source>/<version>`, the job is recorded as `cancelled` with no
result, while the resource version is in fact published (and blocks a re-run
with the same version: 409 "already exists"). Record `done` when `_execute`
returned normally, or at least keep `result` and note that publication
completed.

#### B7 — Client defaults to a personal dev host (Low)

`packages/client/omnipath_client/client.py:28`:
`DEFAULT_API_URL = "https://omnipath-metabo-dev.schaul.click/api"`. A
distributable client should not silently target a personal dev domain; require
`api_url` or default to the official endpoint once one exists. (Remote lazy
reads via httpfs are not checksum-verified, only offline snapshots are; that is
a reasonable trade-off but worth stating in the client README.)

#### B9 — Job queue first-use race: `no such table: jobs` (Medium, reproduced as a test failure)

`api/omnipath_api/jobs/store.py`: readers (`list`/`get`) only check that
`jobs/queue.sqlite3` exists and then query it read-only, but a writer creates
that file at `sqlite3.connect()` and only afterwards runs `CREATE TABLE IF NOT
EXISTS jobs`. The API process and the build worker share this file, so the
first enqueue can make a concurrent `/admin/status` fail with HTTP 500. It
surfaced in this review's run as
`test_admin.py::TestAdminAPI::test_status_inspects_resolver_files` failing
(`sqlite3.OperationalError: no such table: jobs`) while the test's worker thread
started. Fix: create the schema once in `JobStore.__init__` (or create the DB in
a temp file and `rename` it into place), and/or treat "no such table" as an
empty queue in read paths.

#### B8 — Project instructions and some docs reference removed commands/paths (Low)

- `Claude.md` (project instructions) mandates a build-log entry for
  `make load` / `make reload` / `make all`, none of which exist on this branch.
  The equivalent real runs are now `make build`/`omnipath-build build`,
  `make load-postgres`, `make finish-postgres` and `make build-subsets`;
  update the instruction so the logging rule keeps applying.
- `docker-compose.postgres18.yml:4` points to `omnipath_build/db/bitmaps.py`
  (now `packages/postgres/omnipath_postgres/relational/db/bitmaps.py`).
- `docs/payload-dependency-audit.md:22,74,84` use the pre-consolidation
  `packages/omnipath_*/src/...` paths.

### pypath

#### P1 — pypath is not installable on its own: `omnipath-core` dependency is unresolvable (High for merging to master)

`pyproject.toml` now requires `omnipath-core>=0.1.0,<0.2` and
`biolink-model==4.4.4`. `omnipath-core` lives in omnipath-build
(`packages/core`) and is **not on PyPI**; `uv lock --check` in the submodule
fails with "omnipath-core was not found in the package registry". The
committed `uv.lock` still points to `../../full_parquet/packages/omnipath_core`,
a path that no longer exists (left over after fc6c262a6 removed the
`[tool.uv.sources]` entry). Consequences:

- pypath standalone (`pip install`, its own CI, `uv sync` inside the submodule) is broken;
- there is a circular repository dependency (omnipath-build → pypath → omnipath-core in omnipath-build).

Options: publish `omnipath-core` (or the tiny subset pypath uses: naming,
keys, biolink, molecular_forms, measurements, source_attributes) to PyPI;
or make it an optional extra (`pypath-omnipath[inputs-v2]`) so legacy pypath
users are unaffected. Regenerate `uv.lock` either way.

#### P2 — `make test-pypath` fails: stale BindingDB release in test (Low, reproduced)

`pypath/inputs_v2/bindingdb.py:39,43` bumped the default release to `202610`
(9d61c6e6c / dd23d862b) but
`test/test_measurement_review_fixes.py:224` still asserts
`BindingDB_All_202605_tsv.zip`. Result: 1 failed, 143 passed.

#### P3 — Callback/extraction exceptions are no longer swallowed (Medium, behaviour change)

`tabular_builder._CallableSource.extract` and `_BaseCvBuilder._safe_extract`
previously caught every exception and returned `[]`; they now propagate. This
is the right default for data quality, but any existing inputs_v2 mapper that
relied on the silent fallback (e.g. a lambda indexing a missing column on a
minority of rows) will now abort the whole resource build. Worth one full
`make build` sweep over all sources (uncapped or with a large cap) before
merging, and/or a per-row error budget in `two_phase.prepare_worker` that
records and skips failing rows instead of failing the shard.

#### P4 — Most new pypath tests are not run by any target (Medium)

`make test-pypath` lists six files, and the root `pytest` `testpaths` exclude
`pypath/`. The 12 test files added on this branch —
`test_inputs_v2_{communication,pathway,ontology}_migration.py`,
`test_input_source_context.py`, `test_molecular_{catalogue_coverage,forms,sequence_coverage}.py`,
`test_reactome_molecular_coverage.py`, `test_chembl_molecular_queries.py`,
`test_resource_representation_review.py`, `test_download_timeouts.py`,
`test_cv_terms.py` — are never executed by `make test`/`make test-pypath`.
They currently pass (**183 passed** in 8.6 s when run explicitly), so adding them
is cheap. Prefer a glob/marker over a hand-maintained file list.

#### P5 — Obsolete ontology terms are now dropped (Medium, behaviour change)

`inputs_v2/base.py::ontology_term_to_entity` returns `None` for
`term.is_obsolete` (master emitted them with an `is_obsolete` annotation).
Annotations in other resources that still cite an obsolete GO/HPO/MONDO term
now have no ontology entity to attach to, and `replaced_by`/`consider` hints are
not carried anywhere. If this is intended, it should be documented in
`BIOLINK_MODELING.md`; otherwise keep obsolete terms (flagged) or emit
`replaced_by` edges so references can be redirected.

#### P6 — Source fields that were annotations on master are now raw-payload only (Info)

Per `docs/inputs-v2-chemical-migration.md`, fields without a reviewed BioLink
representation intentionally stay only in the evidence payload. For the
record, the BioLink migration (2b7253c4b and follow-ups) stops mapping, among
others: HMDB `cellular/tissue/biospecimen_locations` and `diseases`; ChEMBL
molecule flags (`max_phase`, `withdrawn_flag`, `therapeutic_flag`,
`natural_product`, …) and assay context; Guide to Pharmacology ligand flags
(`Approved`, `Endogenous`, `Withdrawn`, …); LIPID MAPS / RefMet / FooDB /
RaMP class levels; MACdb study context; SIGNOR complex/phenotype/stimulus
names; UniProt `Length`/`Mass`/`Protein families`; SLC tables
`Substrates`/expression. These are no longer filterable or visible as typed
attributes in the API/PostgreSQL. Worth confirming none of the downstream
products or the explorer relied on them (e.g. class levels for chemical
classification, `Approved` for drug filters).

#### P7 — Large formatting-only churn (Low)

Several modules (e.g. `inputs_v2/base.py`) are reflowed line-by-line alongside
semantic changes, which makes the pypath diff (14k+/7k−) much harder to review
upstream. If this is going to pypath `master`, consider splitting the
formatting into its own commit.

## Reviewed without findings

These areas were read with the specific risk in mind and looked sound:

- **Atomic resource publication** (`build/pipeline.py`, `two_phase.py`): per-version
  lock, staging under the data root, manifest written before a same-filesystem
  `rename`, worker errors/exit codes and the free-disk reserve all abort before
  publication; worker shards get disjoint `event_id` ranges on import.
- **Writer aggregation** (`build/writer.py`): undirected-edge flipping also swaps
  scoped annotations and molecular forms; per-event de-duplication of expanded
  resolutions; empty-resource schemas; bucketed nested aggregation.
- **Release validation for PostgreSQL** (`postgres/releases.py`): duplicate JSON
  keys/NaN rejected, symlink and path-escape checks, streamed checksums with
  inode/mtime stability checks, and re-verification before the base checkpoint.
- **Loader locking/checkpoints** (`postgres/loader.py`): session-level advisory
  lock, new-schema requirement, identity check before each checkpoint, COPY row
  counts compared with staged counts.
- **Admin API**: secret compared with `hmac.compare_digest`; mutations refused in
  read-only mode; log/resource/version names validated before building paths;
  deletes refuse versions pinned by a release and run under the engine lock.
- **Static data server** (`deploy/data.nginx.conf`): allow-list of exactly the
  published artifacts, symlinks disabled, everything else 404.
- **Gene/connectivity grouping query** (`api/queries/gene_reference.py`):
  placeholder/parameter ordering in both branches, keyset cursor matches the
  `ORDER BY`.
- **API query pool** (`api/store/query_pool.py`): bounded slots, generation-based
  connection retirement.
- **Molecular forms** (`core/molecular_forms.py`, pypath `_molecular_forms.py`,
  Reactome/MITAB feature parsing, Rust `resolve_molecular_batch`): no invented
  positions for fuzzy ranges, conflicting substitutions collapse to a
  description, entry accessions don't assert a canonical isoform, gene IDs never
  select a protein.

## Not covered

Line-by-line review of the ported main derivation SQL
(`postgres/relational/db/derived_tables.py`, `schema.py`, `bitmaps.py`) and the
subset product SQL; these are guarded by the frozen-main parity tests rather
than by reading. Per-source mapping semantics of each of the ~60 rewritten
inputs_v2 modules were only spot-checked (HMDB, ChEMBL types, Reactome,
MITAB). The Svelte UI was covered by `pnpm check`/tests only.

## Test-suite reduction analysis

Static inventory (AST scan of `packages/*/tests`, plus collection counts):
**797 test functions → 1,391 collected items** (1,207 run by `make test`,
184 `integration`-marked for `make test-postgres`). Timings from one
per-package run (packages in parallel, so contended): core 1.4 s, client 2.1 s,
subsets 1.0 s, postgres 31 s, api 149 s, **build ≈ 50 min** (not completed).

### What does *not* help much

- **Copy-paste duplication is essentially nil.** No exact duplicate test
  bodies; only 2 groups (5 tests, all in `api/tests/test_server.py`) differ
  only in literals → 3 removable by `parametrize`.
- **Argument-validation fan-out** (~150 items such as
  `test_missing_and_whitespace_resolved_keys_are_rejected[16 cases]`,
  `test_nonboolean_modes_before_io[12]`, `test_rejects_*_before_connecting`)
  inflates the count but runs in well under a second in total. Folding them
  into table-driven single tests would shrink the item count without saving
  time; not recommended.
- **Review-regression files** (`test_api_review_regressions.py`,
  `test_build_review_regressions.py`, `test_mapping_review_fixes.py`,
  `test_regression_contract.py`) each test distinct behaviours. Move them next
  to their owning modules for discoverability; no reduction.

### 1. Remove tests of the historical `record_layout` loader (≈ 300+ items) — done

No production module, entry point, Makefile target or deployment imports
`omnipath_postgres.compatibility.record_layout` (2,686 lines); the README keeps
it only "for audited older-schema callers". These test files exercise only that
code:

| File | Items (of which integration) |
| --- | --- |
| `test_loader.py` | 47 |
| `test_duckdb_projection.py` | 45 |
| `test_reactions.py` | 45 (37) |
| `test_duckdb_loader.py` | 34 (34) |
| `test_checkpoint.py` | 25 (25) |
| `test_projection.py` | 24 (1) |
| `test_postgres.py` | 21 (21) |
| `test_payload_validation_cache.py` | 19 |
| `test_payload_dictionary.py` | 14 |
| `test_source_record_audit.py` | 9 (7) |
| `test_sql_contract.py` | 7 |
| `subsets/tests/test_ontology.py` | 7 (7) |
| `test_loader_memory.py` | 5 (5) |
| `test_projection_streaming.py`, `test_payload_owner_statistics.py`, `test_identifier_index.py` | 3 each |
| `test_cli_source_audit.py`, `test_serving_independence.py` | 2, 1 |

Total ≈ **315 items**, including **146 of the 184 integration tests** (only
3 integration items target the current pipeline exclusively; 35 more in
`test_late_constraints.py`/`test_aligned_molecular_context.py` mix both and
need splitting). Deleting the package and these tests together removes ~23 %
of the suite and most of `make test-postgres`'s runtime. **Decision needed:**
is any live PostgreSQL schema still on the historical layout? If yes, keep them
behind an explicit `historical` marker excluded from default runs.

Note the flip side: current-pipeline PostgreSQL coverage lives mostly in
`postgres_dsn`-fixture tests (`test_aligned_*`, `test_main_*`) that are
*skipped* rather than marked `integration`. Mark them `integration` so
`-m integration` actually selects the current pipeline's database tests.

### 2. Share reference-library builds in build/API tests (largest time saving)

A single `build_fixture_library` template build took **91 s** (it is
`lru_cache`d once per process). Five places build their own full library
instead of reusing it: `test_gene_targets.py:210`,
`test_chemical_targets.py:31` (module fixture), `test_mapping_review_fixes.py:214`,
`test_global_reference.py:17`, `test_build_review_regressions.py:110-126`.
At ~90 s each that is roughly **6–8 min** of the build suite. Where their extra
hub rows are additive (Ensembl products, gene targets, chemical targets), add
them to the shared template; for deliberately conflicting fixtures
(`test_global_reference`'s external CID claim) use one second session-scoped
template. Also pre-build the template once per session (session fixture)
rather than inside the first test, so its cost is visible in `--durations`.

### 3. Consolidate full spawned builds

About 16 tests run a complete `build_resource`/`build_all` with `spawn`
worker processes (`test_two_phase.py` 7, `test_versioning.py` 5,
`test_partition_contract.py` 2, `test_pipeline.py`, `test_publication_reader_contract.py`),
each paying interpreter + pypath/resolver import per worker. Most assert
different invariants on the *same* kind of small successful build; share one
module-scoped published build for the success-path assertions and keep
separate runs only for failure/cancel/lock scenarios. Expected: 16 → ~6 full
builds. Exact savings need a `--durations=30` run of `packages/build/tests`.

### 4. Small items

- Parametrize the 3 literal-only clones in `api/tests/test_server.py`.
- Add the 12 unrun pypath test files to `make test-pypath` (+183 items, +9 s) — see P4.
- Fix B9 so `test_admin` is not timing-dependent.

### Expected outcome

| Step | Items | Default-run time |
| --- | --- | --- |
| Today | 1,207 default / 184 integration | ≈ 55 min |
| 1. drop historical loader tests | −~170 default, −146 integration | ~unchanged (they are fast) |
| 2. shared libraries | 0 | −6–8 min |
| 3. shared spawned builds | −0 (assertions merge) | −? (measure) |
| 4. small items | −3 | — |

Step 1 is the big *count* reduction; steps 2–3 are where the *time* is.

## Follow-up: changes made (5 October 2026)

### Historical loader and its tests removed

- Deleted `omnipath_postgres.compatibility` (record-layout loader, 2,686 lines),
  18 PostgreSQL test files, `subsets/tests/test_ontology.py` and one historical
  test in `test_aligned_molecular_context.py`. The one current CLI test in
  `test_late_constraints.py` moved to `test_cli_checkpoint.py`.
- Shared Parquet builders used by current tests moved from `test_projection.py`
  to `packages/postgres/tests/published_fixture.py`; `release_fixture.py` moved
  to `packages/api/tests/` (its only user) and now describes complete builds.
- Dropped the unused psycopg3 dev dependency and the `record-layout` extra.
- Root `conftest.py` marks every test using `postgres_dsn` as `integration`, so
  `-m integration` selects the current database tests (previously they were
  only skipped).
- Collected tests: **1,391 → 1,056** (default 1,207 → 1,031; integration
  184 → 25, all current-pipeline). Integration suite on a disposable
  PostgreSQL 16: **25 passed in ~13 s**.

### New regression tests

`test_queue_readers_tolerate_file_created_before_its_table` (B9),
`test_cancel_after_completed_work_still_records_result` (B6),
`test_publish_pins_build_manifest_bytes`, `test_publish_refuses_unmarked_sample_builds`
(B2/B3, API), `test_rejects_malformed_content_pins`,
`test_content_pin_must_match_build_manifest_bytes`,
`test_sample_builds_require_an_explicit_partial_release` (B2/B3, PostgreSQL),
`test_failed_base_load_removes_its_schema_so_a_rerun_succeeds` (B5, integration),
`test_api_url_falls_back_to_environment` (B7). The B5, B6 and B9 tests were
confirmed to fail without their fixes.

### Test-suite speed

`make test` (1,018 passed, 13 skipped): **≈ 11 min serial → 3 min** (179 s on
4 xdist workers, quiet machine). Before the polling fix it was roughly 50 min.

- `build_reference.Build.stage` polled each stage subprocess with
  `time.sleep(1)`, so each of 62 stages cost ≥ 1 s; it now uses
  `p.wait(timeout=1)`. Stages run in-process when the total hub input is
  ≤ 1 MiB (`OMNIPATH_REFERENCE_INLINE_STAGE_BYTES`, 0 forces subprocesses); real
  builds keep one subprocess per stage. The direct compact index reuses one
  compiler and one process pool per build and lists shard directories directly;
  small Goslin batches parse in-process. Fixture library build: **88 s → ~19 s**.
  The library's outputs (all Parquet rows, every LMDB key/value, normalised JSON)
  and all 292 intermediate stage Parquets were compared before and after and
  are identical apart from compressed byte sizes.
- Tests share reference libraries: three variants (`base`,
  `external-cid-claim`, `exact-structures`) are built once per run, across
  xdist workers via a file lock, and hard-linked read-only into each test
  (previously six builds per process, plus a 45 MB copy per call).
- pytest-xdist added to the dev group; `make test`/`test-api` use
  `-n $(PYTEST_WORKERS)` (default 4), PostgreSQL targets `-n $(PYTEST_PG_WORKERS)`
  (default 2); `PYTEST_WORKERS=0` runs serially. The PostgreSQL fixture uses a
  per-worker port.
- Other checks at the end: `make test-postgres` 234 passed; integration
  (`-n 2`) 25 passed; `make test-pypath` 327 passed; `make check`,
  `check-generated`, `check-web` (incl. lint) and `test-wheels` (7 wheels) pass;
  both `uv.lock` files are current.

Considered but not done: skipping `env.sync(True)` for empty LMDB partitions
(~7 s CPU per fixture library; it is the durability barrier before checkpoint
JSONs) and merging the two chemistry fixture variants (would require renaming
identifiers across chemistry tests).

### Not done

- Consolidating spawned-worker build tests (each is 3–10 s; not worth it).
- Per-step checkpoints for `run_main_derivations` (B5 note).
- A full uncapped `make build` sweep over all sources to flush out mapper errors
  that used to be swallowed (P3).
- Publishing `omnipath-core` to PyPI (P1): until then `pypath-omnipath[inputs-v2]`
  only resolves inside the omnipath-build checkout.
