# OmniPath repository architecture review

Date: 3 October 2026
Branch: `parquet-migration`
Baseline commit: `5eb26de5aae1998335c2f7e2a35afb9f4d687639`
Status: review complete; implementation tracked in [refactoring-progress-20261003.md](refactoring-progress-20261003.md)

## Scope and approach

Review the root workspace, every active package (core, resolver, build, PostgreSQL, subsets, API and web), shared dependencies, packaging, tests, documentation and deployment. Assess the `packages/name/src/name` layout against actual dependency and release boundaries. Preserve the existing Parquet contract and independent Parquet/PostgreSQL publication schedules. This review makes no code or server changes.

The intermediate review log is preserved at the end. Findings below distinguish confirmed defects, maintenance debt and unmeasured operational risks.

## Subsequent decisions

The user chose shorter project folders and removal of the Python `src` layer, with explicit wheel discovery and installed-wheel checks. Subset implementations will live in `omnipath_subsets`, with dependency direction `postgres -> subsets -> core`; one runner replaces the main-layout adapter. These decisions supersede the alternatives evaluated below.

## Recommendation

Keep one repository, a small Python workspace, the separate web application and the Parquet boundary. The main design is sound: resources publish independently, serving can run without the builder, and PostgreSQL consumes explicitly pinned, resolved artifacts. There is no reason to change Parquet schemas, entity resolution policy or PostgreSQL's release schedule for this cleanup.

The most useful work is to finish consolidating migration-era implementations and make the checked-in contracts, commands and documentation agree. Directory depth is a secondary concern. The recommendations below do not require a new scientific pipeline or another full data build.

## Do we need `packagename/src/packagename`?

No. These directories represent three different things:

| Level | Purpose |
| --- | --- |
| `packages/omnipath_core/` | Independently buildable project: metadata, tests and documentation |
| `src/` | Keeps importable code separate from project tooling and test files |
| `omnipath_core/` | The namespace used by `import omnipath_core` |

`src` makes it harder for tests to accidentally import code from the checkout that would be missing from an installed wheel. That is useful here because the API and worker are installed as wheels and the resolver is a native extension. This is the conventional tradeoff described by [PyPA](https://packaging.python.org/en/latest/discussions/src-layout-vs-flat-layout/); [uv workspaces](https://docs.astral.sh/uv/concepts/projects/workspaces/) support this multi-project layout with one lockfile.

**My preference:** retain `src`, and optionally shorten outer project directories to `packages/core/src/omnipath_core`, `packages/build/src/omnipath_build`, etc. Public imports and distribution names can remain unchanged. Do this only after fixing boundaries: it requires updating workspace paths, Dockerfiles, tests, docs and the source-checkout lookup for the Rust helper, but has no inherent runtime benefit.

A flat layout per project is valid, but mainly removes one path segment. A single root distribution with every package under one `src/` is simpler to configure, but makes the independent serving installation and native build dependencies harder to separate. I would not combine all Python packages into one distribution.

## What already works well

- The API's required dependencies are core/serving libraries; build integration is optional. The API image does not install pypath or compile the resolver.
- Immutable resource directories, explicit releases, resource checksums and a local DuckDB projection stage form a useful boundary between servers.
- PostgreSQL keeps its own selected release. Base and complete-product checkpoints preserve completed work, and the implementation records information that cannot be recovered from aggregated Parquets instead of inventing it.
- The builder uses bounded worker shards, one final aggregation and atomic publication. Resource parsing, resolution and serving remain distinguishable stages.
- The web app already separates generated transport types, domain adapters, selection state, query cancellation and the server proxy. Keep these improvements.
- Package-local tests and Dockerfiles plus a thin root Makefile are a good arrangement. The PostgreSQL wheel correctly contains all 15 expected SQL/YAML assets.

## Prioritized findings

Priority 1 means fix before presenting the repository as a clean, reusable setup. Priority 2 is the next consolidation work. Priority 3 is optional cleanup or an improvement to make when touching that code.

### R1 — Make the tests reproducible from tracked files (priority 1, reproduced)

**Evidence:** `packages/omnipath_postgres/tests/test_aligned_ontology_memory.py:189` directly reads ignored `main-reference/9f9bb709c764/omnipath_build/duckdb_load.py`; another direct read is at line 326. `test_main_parity_contract.py` already has a tracked `legacy/postgres` fallback. `test_main_cosmos_direction_lookup.py` reads an old revision through `git show`.

**Impact:** the current checkout passes, but a colleague's fresh checkout can fail without any code regression. I reproduced the ontology test's `FileNotFoundError` in a temporary snapshot containing the tracked PostgreSQL tests and legacy reference. The history-based test also needs a full enough Git history.

**Change:** use one immutable, hash-checked test oracle from tracked fixtures or the existing tracked legacy source. Route all parity tests through one fixture loader. Keep this oracle independent from candidate code. Do not solve the problem by silently skipping the test.

### R2 — Correct PostgreSQL's checked-in startup defaults (priority 1, confirmed configuration issue)

**Evidence:** `docker-compose.postgres18.yml` publishes `55432:5432` without a loopback host and includes a fixed default database password. Its memory/parallelism defaults are sized for beauty rather than a typical developer machine.

**Impact:** following this Compose recipe on another host can expose a database more widely than intended and apply unsuitable resource settings. This review did not inspect or change the live PostgreSQL deployments; it is not a finding that those services are exposed.

**Change:** bind to `127.0.0.1` by default, require explicit credentials, and provide named developer/server settings. Include a short PostgreSQL setup section in the root workflow or a `postgres-up` wrapper so users can find the intended recipe. Keep deployment-specific production settings separate.

### R3 — Make core describe the actual artifact contract (priority 2, confirmed)

**Evidence:** `omnipath_core/versioning.py:26` publicly exports `ManifestFile` with `bytes`/`num_rows`; `BuildManifest` uses `schema_version="2.0"` and a list of files. `omnipath_build/pipeline.py:68` publishes integer schema version 1, serving schema version 3 and a filename-keyed mapping containing `size_bytes`/`rows`. `omnipath_postgres/releases.py` reads that latter shape. Core and build also duplicate resource/version validation.

**Impact:** a developer using the advertised shared model produces a different format from the functioning pipeline. Current builds avoid the problem by bypassing the model.

**Change:** define accurate resource/release metadata models and version constants in a small core contract module, matching existing files byte-for-byte where identity matters. Remove or deprecate the obsolete model. Reuse shared structural validation, while keeping API discovery and PostgreSQL's stricter artifact verification as separate policies. Add a publisher-to-reader contract fixture. **Do not change Parquet outputs.**

### R4 — Consolidate PostgreSQL and product orchestration (priority 2, confirmed duplication)

**Evidence:** the package-level PostgreSQL API uses `aligned_loader.py`, while `loader.py`, the earlier schema/projection/derivation modules and historical product readers remain active importable code. Current scientific implementations live under `omnipath_postgres/main_compat`. `omnipath_subsets/main_adapter.py` repeats release checks, locks, product publication and checkpoint metadata logic. Even the new loader imports `validate_schema` from the old loader; the subset adapter imports private `_widen_temp_buffers`.

**Impact:** there are several places to maintain transaction/checkpoint behavior, two PostgreSQL drivers, and a package tree that suggests a separation the implementation does not have. Changes risk reaching only one path.

**Change:** establish one current loader, one product runner and shared small database helpers. Preserve the distinct **resume** and **explicit rebuild** operations. Move historical record-layout code behind an explicitly named compatibility namespace; retain it until its consumers are audited. Keep the frozen main oracle as tests, not as a second apparent production implementation.

For the current scope, I would keep the PostgreSQL-specific products under `omnipath_postgres.products` and reduce `omnipath_subsets` to a compatibility facade, retiring that distribution only after callers migrate. If separate product ownership is important, moving the scientific implementations into `omnipath_subsets` is possible, but the top-level orchestrator must then avoid a postgres↔subsets dependency cycle. I would not introduce that larger change solely to preserve the original proposed folder tree.

### R5 — Give the resolver boundary a clear meaning (priority 2, design debt)

**Evidence:** `omnipath_resolver` exposes one native batch function. `LibraryMatcher`, `EntityResolver`, normalization, lookup and enrichment live in build's `canonical/`, `resolver.py` and `reference/full_index_runtime.py`. Runtime lookup imports namespace codes and observation helpers from the replay utility module. The extension also enables the Rust crate's `parquet-input` feature even though the runtime policy does not read Parquet.

**Impact:** the nominal resolver cannot be used as the complete matcher without installing/importing build internals. Runtime policy is harder to find and shares module ownership with conversion/replay tooling.

**Change:** keep reference construction under build, but put the stable matching facade, index reader and identifier/key policy under resolver, with explicit observation/result contracts. Move shared namespace encoding out of the replay script. Alternatively, explicitly name/document the existing package as a native kernel if independent matching is not intended. Separate the build-only Parquet feature from the kernel where feasible. Preserve current biological behavior and parity tests.

### R6 — Make dependency and submodule policy honest (priority 2, confirmed)

**Evidence:** `.gitmodules` contains pypath, cache-manager, download-manager and omnipath-utils. The frozen lock resolves pypath from its editable checkout, cachedir 0.1.3 and dlmachine 0.0.3 from PyPI, and no omnipath-utils distribution. Setup initializes all submodules and the worker Dockerfile copies them. `docs/pipeline/1_setup.md` claims local editable cache/download packages are used.

**Impact:** local changes in three apparent dependency checkouts do not control the locked runtime. This is especially confusing because build carries a runtime monkeypatch for cachedir.

**Change:** deliberately choose published, pinned dependencies or local development sources for each repository. Remove unused submodule/container coupling only after that choice; do not silently swap dependency versions. Fix the cachedir defect in its owning package and replace the monkeypatch when the selected release contains the fix. Keep environment-specific proxy handling in a launcher/configuration layer rather than the scientific build function.

### R7 — Preserve exact software provenance in wheels and containers (priority 2, confirmed limitation)

**Evidence:** `omnipath_build/provenance.py:24` records dependency version strings and a Git revision when available; its fallback hash covers builder Python files only. Installed worker wheels have no checkout `.git`. Core, resolver and build currently all use version 0.1.0, and pypath branch changes retain the same distribution version.

**Impact:** two artifacts can report the same version information while using different native/core/pypath code. The serving deployment report records image revisions, but that is external to an individual resource manifest.

**Change:** generate a small build-info record when constructing wheels/images, containing relevant source revisions or wheel digests and the vocabulary version. Reuse it in resource provenance. Keep package versions, resource versions, manifest schema versions and monthly release identifiers distinct.

### R8 — Make quality commands green and automate the intended coverage (priority 2, confirmed)

**Evidence:** root Ruff reports 38 unused-import errors, all in the immutable CV oracle's `__init__.py`. Formatting reports 63 files. No tracked CI configuration was found. Root `make test` does not itself execute native Rust, web/browser, pypath-specific or isolated-wheel checks; those have separate commands.

**Change:** exclude immutable vendored fixtures from style transformations, format maintained code separately, then add a small CI workflow or documented equivalent. Include Python tests, native policy tests, web checks and OpenAPI/generated vocabulary freshness. Add a bounded installed-wheel smoke test so workspace-wide dependencies cannot mask missing declarations. Run relevant disposable database fixtures in a dedicated opt-in CI job. These are developer checks, not a mandatory post-build production verification phase.

### R9 — Restore generated-file ownership and remove hidden web dependencies (priority 2, confirmed)

**Evidence:** `omnipath_web/src/lib/domain/biolink.generated.ts:1` names absent `scripts/sync_biolink.py` as its generator. `src/lib/skills/index.ts:1` imports production content from `.cursor/skills`; the web Dockerfile consequently copies `.cursor`. OpenAPI generation is present and its freshness check passes.

**Change:** restore a generator/check for web vocabulary from the pinned Biolink source. Put shipped skill content under a normal content directory and have editor integrations reference/copy it. Expose generation commands through the root Makefile. Keep generated transport/domain artifacts separate from handwritten adapters.

### R10 — Set an API resource budget across concurrent queries (priority 2, capacity risk)

**Evidence:** `omnipath_api/store/connection.py:12` creates an independent DuckDB database with a default 4 GB limit. `engine.py:165` creates one per worker thread. The query adapter has no shared query-admission limit, DuckDB threads are not explicitly bounded there, and serving Compose has no API CPU/memory budget.

**Impact:** the per-connection limit is not a service-wide memory limit; concurrent heavy queries can multiply memory and CPU demand. I did not run a load test or observe a live failure.

**Change:** expose query thread/memory/concurrency settings together and use a small bounded query pool/admission limit, with a container budget as a backstop. Benchmark overlapping search, facets and exports against representative artifacts before choosing defaults. Clean closed connections out of the tracking collection during generation changes. Keep immutable request inventory snapshots and bounded result caches.

### R11 — Establish one current documentation path (priority 2, confirmed)

**Evidence:** `docs/pipeline/README.md` calls the old pre-Parquet flow current and recommends removed targets such as `make all`, `make load` and `make derive`. Root README says the subsets CLI is historical-only, although it supports aligned explicit rebuilds. API README still shows the prototype image for index preparation and describes cloud input, while its inventory root is a local filesystem path. Two reference-guide benchmark links are broken.

**Change:** keep root README as the launchpad; maintain current architecture, reference methods and operations guides; move historical plans/run logs into clearly labelled history/reports. Document local-mounted API artifacts separately from PostgreSQL's new HTTPS input. Document `finish` versus `rebuild`. Keep deployment commands next to their package/Compose implementation and link them from root.

### R12 — Remove generated data from the source checkout (priority 2, confirmed)

**Evidence:** a tracked root filename consisting of one space is a 49,700,475-byte Parquet file. `test.ipynb` is 37,658,590 bytes with substantial saved outputs. There is also an unused-looking root `package-lock.json` with no root Node package, and a Nix shell containing Meilisearch but lacking the Node/pnpm prerequisites now documented for setup.

**Change:** establish the artifact/notebook provenance, then move data to artifact storage, strip notebook outputs or archive the notebook, and remove obsolete root tooling through a separate review. The two large files account for 87.4 MB of current tracked content. Deleting them from the current tree does not shrink historical Git objects; no history rewrite is needed for the initial cleanup. Do not delete `legacy/` indiscriminately: it is still used as a parity oracle.

### R13 — Treat remote verification traffic as a measured cost (priority 2, performance follow-up)

**Evidence:** `omnipath_postgres/releases.py:200` streams every artifact for SHA256 at initial validation and again in `verify_release`; this includes `evidence_payloads.parquet`, which is not loaded into normalized PostgreSQL. DuckDB separately reads/stages entities and relations. The small successful HTTPS test is a functional check only.

**Change:** report validation bytes/time separately from projection, and benchmark a representative large pinned resource over HTTPS. Then choose a deliberate verification policy or verified staging strategy if transfer dominates. Retain immutable pins and corruption detection. Avoid promising that direct reads transfer every byte only once, and do not change the Parquet format to address this.

### R14 — Keep role-specific builds and avoid unnecessary rewrites (priority 3)

- Select API/build dependencies for the worker image instead of `--all-packages`, which also installs PostgreSQL and historical subset packages. Keep the already slim serving image.
- Offer clear setup recipes for serving, build and PostgreSQL roles; the convenience `make setup` can still install everything for maintainers.
- Split large files when changing their responsibility: builder writer/reference orchestration, PostgreSQL schema/derivations and web details/admin components. File length alone is not a reason to rewrite tested SQL or create more packages.
- API query mixins rely on a substantial implicit shared engine interface. Small typed protocols/context objects would make this visible; a wholesale service framework rewrite is unnecessary.
- Keep pypath input modules in pypath. PostgreSQL's pypath use is CV/resource/license metadata, not a second resolver. A small packaged metadata snapshot could eventually remove that dependency while preserving frozen main identifiers and physical contracts.

## Package-by-package conclusion

| Package | Keep | Improve | Boundary recommendation |
| --- | --- | --- | --- |
| `omnipath_core` | Shared Parquet schemas, keys, quantities, Biolink names and lightweight metadata | Correct manifest models; deduplicate validators; distinguish input models from published contracts; clarify bundled/generated vocabulary ownership | Keep a separate shared package. Avoid turning it into a store for unrelated pipeline behavior. |
| `omnipath_resolver` | Rust policy kernel, Python adapter, native policy tests | Complete runtime ownership or name it explicitly as a kernel; separate build-only Parquet dependency; expose a public batch API | Keep separate because native packaging is a real boundary. |
| `omnipath_build` | Input discovery, reference construction, bounded preparation, finalization and atomic publication | Move matching internals to their chosen owner; reduce replay/compiler/runtime cross-imports; upstream compatibility patches; embed provenance | Own everything before immutable resource publication. Release publication should be usable without installing a web server. |
| `omnipath_postgres` | Local/HTTPS pins, one DuckDB stage, bulk COPY, main-compatible schema, checkpoints and scientific products | Consolidate active loader/product state machine, isolate historical layout, reduce pypath metadata coupling and dual-driver leakage | Own the current PostgreSQL model and database-specific products; keep the main scientific contract. |
| `omnipath_subsets` | Existing explicit-rebuild entry point and compatibility promises | Remove duplicated orchestration; make rebuild/resume distinction obvious; isolate historical readers | Prefer a thin compatibility facade now; retire it after consumers move, or extract real product ownership only when needed. |
| `omnipath_api` | Slim install, read-only default, optional worker integration, release scopes, query caches, serving indexes | Bound concurrency/resources; keep publication contracts in shared/build code; reduce implicit mixin interface | Keep separate from builder and PostgreSQL. Current artifact inventory uses a local mount. |
| `omnipath_web` | Independent Node package, HTTP boundary, generated contracts, domain adapters and feature state | Restore Biolink generation, move shipped content out of `.cursor`, split the largest UI controllers as they evolve | Keep separate; it does not need Python packaging or a nested Python-style namespace. |

The external repositories were reviewed at their integration boundary (gitlinks, dependency metadata and usage), not as a full audit of every pypath/cache/download utility module. The Python client is still a later milestone; no empty client package is needed now.

## Suggested target structure

This is an optional destination after consolidation, not a prerequisite for fixes. The shorter outer folder names do not change public Python imports.

```text
omnipath-build/
├── pyproject.toml, uv.lock, Makefile, README.md
├── packages/
│   ├── core/       # pyproject.toml + src/omnipath_core + tests
│   ├── resolver/   # Python facade, Rust policy kernel and tests
│   ├── build/      # inputs, reference construction, resource publication
│   ├── postgres/   # release input, projection, schema, derived, products
│   ├── api/        # query service, HTTP routes, optional job adapters
│   ├── web/        # SvelteKit application and shipped content
│   └── subsets/    # temporary compatibility facade, if still needed
├── scripts/        # small developer/generation entry points
├── deploy/         # Compose, environment examples and operating instructions
├── docs/           # current architecture/methods; separate history/reports
└── tests/oracles/  # explicitly frozen external compatibility fixtures
```

Do not add a new distribution for every internal concern. `projection`, `schema`, `derived` and `products` can be modules inside PostgreSQL. Reference compilation and matching can be distinct owners without requiring separate repositories. Keep server separation at the artifact and HTTP boundaries.

## Proposed sequence

1. **Make a clean checkout trustworthy.** Repair frozen fixtures and quality commands; correct current docs and database startup defaults; restore missing generation tooling. No data changes.
2. **Consolidate contracts and provenance.** Make core match the existing manifest format; centralize small validators; embed software build information; clarify dependency/submodule selection.
3. **Consolidate PostgreSQL/product APIs.** One current loader and product runner, explicit resume/rebuild modes, compatibility namespace for the old layout. Preserve schema/index/output and checkpoint tests.
4. **Clarify the resolver runtime.** Move the facade/index reader/shared key encoding together behind one public API. Keep reference construction and resource orchestration in build.
5. **Tune operational packaging and budgets.** Narrow worker dependencies, add role-specific setup examples and API query budgets, then measure a representative HTTPS transfer.
6. **Finish cosmetic cleanup.** Move/archive tracked data, simplify outer folder names if desired, and split large modules incrementally. Keep this separate from scientific or SQL changes.

Each step should be a small reviewable change. Existing fixture parity tests protect the scientific behavior; none of these steps requires automatically rebuilding all resources or restoring the removed mandatory full-release verification phase.

## Validation performed

| Check | Result |
| --- | --- |
| Default Python test collection | 1,177 tests collected |
| Non-integration suite on current checkout | 822 passed; 32 skipped; 323 deselected; 499.54 s |
| Clean tracked snapshot: one frozen-main ontology test | Reproduced missing ignored-reference file failure |
| Native reference policy tests | 28 passed |
| Web unit tests | 26 passed |
| Svelte check | 0 errors, 0 warnings |
| Web ESLint / Prettier | Passed |
| OpenAPI freshness | Passed |
| Python Ruff | 38 F401 errors in one frozen oracle file |
| Python format check | 63 files would be reformatted |
| PostgreSQL wheel contents | All 15 expected SQL/YAML assets included |
| Relative Markdown link scan | Two missing benchmark-report links |

No production services, live databases or published artifacts were changed. Tests used temporary fixtures; no live resource download/build or full PostgreSQL run was requested. Full database integration, browser end-to-end tests, live concurrency testing and representative large HTTPS benchmarks were not rerun for this structural review. Passing the local suite does not remove the reproduced fresh-checkout issue.

## Intermediate review log

### 1. Baseline and inventory

- The consolidated repository is clean at the baseline above.
- Six Python packages and one Svelte web package are present under `packages/`; four Python source repositories are included as submodules.
- Root setup, serving deployment and HTTPS PostgreSQL support were added in the most recent commits.
- Detailed package and cross-package review is underway.

### 2. Workspace and package boundaries — initial findings

- `packages/<project>/src/<import_name>` is a conventional layout for separately installable Python distributions. The outer directory owns configuration/tests/docs; `src` isolates importable code; the inner name is Python's import namespace. It is optional, but deleting `src` would mostly shorten paths rather than simplify ownership. A lower-risk cosmetic option is `packages/core/src/omnipath_core`, etc.
- The workspace and single lockfile suit the current tightly coordinated development. Keep independent resource versions, release manifests and deployment schedules; Python distribution versioning is a separate concern.
- The resolver package is a small Rust/Python kernel. Reference lookup, normalization, enrichment and the runtime facade live under `omnipath_build/canonical` and `omnipath_build/reference`. Its name promises a larger boundary than it currently owns. Investigate either making this explicit or moving the complete runtime resolver behind its public API.
- PostgreSQL contains both the earlier resource-record loader and the newer main-layout pipeline. `main_compat` also owns current MetSigDB, network and COSMOS implementations, while `omnipath_subsets` still contains the earlier implementations plus aligned adapters. This is migration debt worth resolving before cosmetic directory changes.
- The root README describes the subsets CLI as historical-only, but the package README says it detects aligned schemas and can explicitly rebuild them. Trace the implementation to establish the correct contract.
- No tracked GitHub Actions configuration was found. Check existing test/lint commands and whether package isolation, native Rust policy tests and browser tests are exercised elsewhere.

Packaging references: [PyPA src vs flat layout](https://packaging.python.org/en/latest/discussions/src-layout-vs-flat-layout/), [uv workspaces](https://docs.astral.sh/uv/concepts/projects/workspaces/).

### 3. Contracts, PostgreSQL and developer checks

- Confirmed: the subsets CLI detects `parquet_release` and delegates aligned schemas to `main_compat` builders. Its semantics are **explicit rebuild**, whereas `omnipath-postgres --finish` **skips committed products**. Root documentation currently conflates these; retain both operations but name them clearly.
- Confirmed contract drift: publicly exported `omnipath_core.BuildManifest` uses `schema_version="2.0"`, a list of files and `bytes`/`num_rows`; the actual publisher and PostgreSQL reader use integer version 1, a filename-keyed mapping and `size_bytes`/`rows`. The live path works because it bypasses the exported model. Fix the Python contract to describe existing Parquets; do not change published artifacts.
- Confirmed duplicated resource/version validators in core and build; release validation also lives independently in API and PostgreSQL. A small shared contract module would reduce divergence without merging publication schedules or package responsibilities.
- Default test collection succeeds: 1,177 tests collected. The non-integration suite is running.
- Root lint currently fails: 38 F401 diagnostics in the frozen CV oracle; formatting check reports 63 files. Preserve immutable oracle content and exclude vendored fixtures from style checks, then format maintained code separately.
- Some PostgreSQL tests directly read ignored `main-reference/9f9bb709c764` files; other parity tests have a tracked `legacy/postgres` fallback. Two direct reads need the same reproducible fixture strategy. A separate COSMOS fixture uses `git show` on an older commit and therefore also depends on available history.
- PostgreSQL still depends on pypath for CV/resource/license metadata, despite not parsing or resolving inputs. This is a dependency weight/ownership issue, not evidence that resolution is being repeated. Both psycopg2 and psycopg3 remain because current and historical paths use different drivers.
- Lockfile inspection: pypath is editable, but cachedir/dlmachine are installed from PyPI and omnipath-utils is absent. Their submodules are still initialized/copied, and old setup documentation claims they are editable. Make source policy explicit; local edits to those submodules currently do not control the locked runtime.

### 4. Documentation and deployment review — initial findings

- The numbered `docs/pipeline` guide still presents the pre-Parquet build as current and recommends removed Make targets. Archive/label historical migration material and route users to one current architecture guide.
- Two links in `docs/reference-resolver.md` point to missing benchmark reports.
- The web app imports shipped skill Markdown from repository `.cursor/skills`, so an editor-specific hidden directory is a production build dependency. Move canonical shipped content under web/content or docs and treat editor copies as adapters.
- The PostgreSQL development Compose file publishes its port on all host interfaces and supplies a fixed default password. Prefer loopback binding and explicit credentials; this finding is about the checked-in startup default, not a claim about the protected deployed PostgreSQL services.
- HTTPS loading correctly stages entity/relation inputs once for projection, but hash validation streams all artifacts before loading and again before the base checkpoint, including evidence payloads not loaded into PostgreSQL. This is real extra network transfer, not a Parquet cache. Keep integrity requirements explicit and benchmark representative transfer size before promising performance.
- Tracked root clutter includes a roughly 50 MB file whose name is a single space and a roughly 38 MB notebook. Their purpose/contents are being checked before recommending disposition; nothing is deleted.

### 5. Completed developer checks and deeper package findings

- Non-integration suite: **822 passed, 32 skipped, 323 deselected**, 499.54 seconds. The current checkout has the ignored frozen reference available.
- Clean-snapshot reproduction: the frozen-main ontology writer test fails with `FileNotFoundError` when run from tracked tests plus tracked legacy oracle, because it hardcodes ignored `main-reference`. This is now reproduced, not just inferred.
- Native reference policy tests: **28 passed**. Web: **26 tests passed**, Svelte **0 errors/0 warnings**, ESLint and Prettier passed. OpenAPI freshness check passed.
- PostgreSQL wheel built successfully and contains **all 15 expected SQL/YAML assets**. Its basic `src` packaging is working; flattening directories would not fix the more important issues.
- The single-space root file is an **Apache Parquet artifact, 49,700,475 bytes**. The notebook is **37,658,590 bytes**, with 28 code cells and large saved outputs. Together they add 87.4 MB of current tracked content. Move data out of the source tree and strip/archive notebook output through a separate reviewed cleanup.
- The generated web Biolink file names `scripts/sync_biolink.py` as its generator, but that script is absent. Keep generated outputs, restore their reproducible generation/check commands, and use the pinned Biolink version as the source.
- The current API correctly makes build dependencies optional. However, its worker image builds every workspace package, including PostgreSQL and historical subsets. Select only API/build dependencies for that image.
- API query connections are separate per thread, default to 4 GB each, and have no explicit DuckDB thread setting or shared query admission limit. Compose supplies no API memory/CPU budget. This is a capacity risk under concurrent work, not a measured live failure. Budget concurrency and per-query resources together.
- Wheel/container build provenance cannot rely on `git rev-parse`: the worker wheels omit `.git`. The fallback hashes only builder Python files; core/native/pypath remain identified mainly by unchanged package version numbers. Embed source revisions or wheel digests when building the image so standalone artifacts retain exact software provenance.
- Current and historical PostgreSQL paths are coupled even at small helpers: the aligned loader imports `validate_schema` from the historical loader, and subsets imports a private MetSigDB tuning helper. Isolate shared helpers and consolidate product checkpoint orchestration before moving large SQL modules.
