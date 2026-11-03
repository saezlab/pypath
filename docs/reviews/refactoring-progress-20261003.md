# Consolidation implementation log

Baseline: `5eb26de5aae1998335c2f7e2a35afb9f4d687639`, branch `parquet-migration`.
Status: implementation and local validation complete; no deployment performed.

## Agreed destination

Short outer package names (`core`, `resolver`, `build`, `postgres`, `subsets`, `api`, `web`), flat Python packages with unchanged public import/distribution names. Rust and Svelte retain their language-specific `src` directories. PostgreSQL consumes resolved Parquet; subsets owns the three scientific products and a shared runner. Dependency direction: postgres -> subsets -> core. The extra main-layout adapter is removed. Finish resumes completed work; the subsets command explicitly rebuilds requested products. Preserve schema/index/scientific output and publication schedules.

## Scope from the architecture review

- R1/R4: tracked immutable test oracles, one current PostgreSQL pipeline, isolated historical compatibility, actual subset ownership and one durable product runner.
- R3/R7: shared models describe existing publication files; exact software provenance in installed distributions.
- R5: independent resolver runtime; construction remains in build; no biological policy change.
- R2/R6/R8/R9/R10/R11/R14: safe local PostgreSQL defaults, honest pinned dependencies, green developer checks/CI, reproducible vocabulary generation, visible web content, bounded API query resources, current docs and role-specific packaging.
- R12: preserve provenance of tracked exploratory data before removing generated outputs from the source checkout; no Git history rewrite or server cleanup.
- R13: account for remote validation traffic and provide a bounded measurement path; preserve immutable pins and corruption detection.

No resource/reference rebuild, full PostgreSQL build, production deployment or mandatory post-build verification is part of this refactor.

## Progress

### Layout and parallel implementation

- Renamed all seven project directories and flattened the six Python import packages, retaining names such as `omnipath_core`.
- Updated workspace, package discovery, Docker/Makefile/source paths and ignore patterns. Explicit package discovery protects flat-layout wheels from accidental tests/tooling inclusion.
- Assigned independent PostgreSQL/subsets, resolver, and core/provenance work. Parent owns integration, API/web/tooling/docs and reviews all changes.
- Refreshing the workspace environment for the new paths; validation remains pending.

### Implemented boundaries and independent review

- Resolver runtime moved to resolver with no build dependency; construction and replay remain in build. Native extension no longer enables build-only Parquet support. Biological algorithms were checked by AST comparison and frozen output fixtures.
- Current PostgreSQL modules are canonical loader/projection/relational modules. Historical record-layout code is explicitly isolated. Actual scientific products and all 12 SQL assets moved to subsets without changing algorithms or SQL bytes. One runner owns publication, release locks and product checkpoints; no reverse subsets -> postgres dependency.
- Shared core manifests now describe existing files; callers retain their storage/integrity policies. Installed-file provenance covers core/build/resolver/pypath/Biolink, native binary and vocabulary even without Git.
- API request pool bounds simultaneous and retained DuckDB databases, with configurable budgets and retryable overload responses. Independent review caught an admin refresh bypass; reload now uses the same pool and has a regression test. A second finding corrected role commands to retain their selected package environment rather than resyncing the full workspace.
- Web content is package-owned; editor skills link to it. Restored vocabulary generation preserves every existing entry. Worker images select only API/build dependencies. Local PostgreSQL binds to loopback and requires an explicit password; large-server settings are separate examples.
- Tracked exploratory Parquet/notebook bytes were preserved locally under ignored local-artifacts; a stripped notebook and hash/provenance manifest remain tracked. No server cleanup or Git history rewrite occurred.
- HTTP validation reports actual response bytes/request counts and separate timings, while projection traffic remains explicitly separate. Added a bounded single-resource measurement command.

### Validation underway

- Resolver: 131 Python tests passed, one optional skip; 29 native tests passed; isolated resolver wheel with real tiny-reference matching passed without build, DuckDB or pypath.
- API: all 145 tests passed before the final admin-refresh fix; its six focused pool tests then passed, including exhaustion/reuse/generation cleanup/concurrent result parity/admin refresh.
- Every package passed an isolated wheel installation with only declared locked dependencies, including SQL/YAML assets and real wheel-only software provenance.
- Web checks, 26 unit tests, lint, formatting and production build passed. Generated OpenAPI/Biolink checks, frozen lock and Compose syntax passed.
- Combined PostgreSQL fixtures initially exposed historical test callers still pointing at new canonical APIs. Those test callers are being routed to explicit compatibility entry points; current CLI/rebuild coverage is also being added. This is not a production data error; the tests use disposable local databases.
- Full combined Python and PostgreSQL fixture runs remain in progress. Completion has not been declared.

## Final validation and review

- Combined Python fixture run: **871 passed, 32 skipped, 323 integration tests deselected**, 515.88 seconds. Final API changes: **147 passed**, plus the new worker-environment and publisher-environment tests (**2 passed**). Two obsolete proxy-rewriting tests were removed with their unused helper. The current tree collects **1,231 tests**.
- Complete PostgreSQL/subsets suite against a disposable local cluster: **747 passed, 11 existing opt-in skips**, 65.10 seconds. This includes physical/scientific parity, local/HTTP projection, constraints, shared runner locking, durable callbacks, atomic rollback, resume, and explicit rebuild through current CLIs. No live database was used.
- Resolver matching/reference tests: **131 passed, one optional skip**; **29 native policy/component tests** passed. Independent runtime wheel performs synthetic LMDB lookup, native matching and enrichment without build, pypath or DuckDB.
- Fresh source snapshot: frozen oracle tests pass with no `.git` and no ignored `main-reference`; all 76 imported frozen source files are hash-checked.
- Every current package passes an isolated installation with only declared locked dependencies. Wheels are now built from clean source distributions and checked for exact Python/SQL/YAML/dictionary membership and bytes, catching stale modules left in build/lib after moves. Final build/API wheels were rerun after the last worker-environment fix. Real installed builder provenance includes core/native/pypath/vocabulary identities without Git.
- Web: **26 tests**, Svelte check, ESLint, formatting and production build passed. Python Ruff and format checks pass. OpenAPI and Biolink generation are current; Biolink values are identical to the original generated file. Frozen dependency lock, Compose syntax and current documentation links pass.
- Independent review found and fixed two API admission edges (admin refresh and delete-before-refresh) and role commands accidentally restoring the full workspace. PostgreSQL/subsets review confirmed unchanged scientific algorithms and all 12 product SQL files byte-identical, with preserved identity/lock/checkpoint guarantees.

## Remaining operational choices

The implementation is local. No resource/reference rebuild, PostgreSQL rebuild, deployment, remote service change, Git history rewrite or mandatory production verification phase was introduced. API budget defaults are configurable starting values, with synthetic concurrent-query parity coverage; full-release throughput/load benchmarks remain an operational sizing task. The bounded HTTPS measurement tool and actual validation traffic counters are ready for a chosen representative resource; verification policy remains unchanged.

Published cachedir/dlmachine pins are retained. Their optional source checkouts remain available to developers but are not selected or copied into the runtime image. The known cachedir compatibility patch stays until a selected upstream release contains the fix; no dependency release was silently replaced. Optional large-file splitting, a new Python client and removal of pypath's metadata dependency are outside this cleanup.

## Follow-up: remove historical subsets (2026-10-03)

At the user's request, removed `omnipath_subsets.compatibility` and its optional
psycopg3 dependency, historical-only product tests, and the obsolete full-release
verification/rebuild script and its tests. Current product implementations, SQL,
shared runner, checkpoint behavior and independent frozen parity references remain
unchanged. Renamed current runner tests to match their actual owner. The retained
PostgreSQL identifier-index regression queries the catalog directly instead of
importing the deleted audit script. Updated the package guide to describe only the
current products. The separate historical PostgreSQL loader is outside this cleanup.

Validation: **601 passed, 11 existing opt-in skips** in 43.59 seconds across the
PostgreSQL/subsets suite and API migration-serving fixture, using a disposable local
PostgreSQL cluster. Full test collection: **1,083 tests**, no dangling imports.
The isolated subsets wheel passes with neither the compatibility namespace nor
psycopg3 installed. Ruff, formatting, dependency lock and whitespace checks pass.
Current scientific implementations, all product SQL and the PostgreSQL runtime
have no diff. No deployment or persistent database changes.
