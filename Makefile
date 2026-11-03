.PHONY: help setup setup-python setup-web test test-pypath test-postgres test-subsets test-api check check-web build sample hubs reference publish-release api web dev load-postgres finish-postgres build-subsets serving-build serving-up serving-status serving-logs serving-stop worker-up docs

export PKG_INFRA_CONFIG ?= $(CURDIR)/config/pkg_infra_quiet.yaml
SOURCE ?= signor
VERSION ?= 0.1.0
MAX_RECORDS ?= 20
DATA_ROOT ?= data
POSTGRES_SCHEMA ?= omnipath
POSTGRES_ARGS ?=
# Identity layer (packages/build/REFERENCE.md): hub exports -> per-hub indexes -> identity library.
HUBS_DIR ?= $(DATA_ROOT)/reference/hubs
HUB_INDEX_ROOT ?= $(DATA_ROOT)/reference/hub-index
IDENTITY_DIR ?= $(DATA_ROOT)/reference/identity
# Every hub the identity layer indexes; ontology only supplies term labels and has no kv store.
IDENTITY_HUBS ?= $(shell uv run --frozen --package omnipath-build python -c 'from omnipath_build.identity.common import HUBS, ONTOLOGY; print(*HUBS, *ONTOLOGY)')
HUB_INDEX_ARGS ?=
HUB_KV_ARGS ?=
IDENTITY_ARGS ?=
# pytest-xdist: PYTEST_WORKERS=0 runs serially; PYTEST_WORKERS=auto uses every core. 4 is the
# default because the build tests already fork worker processes and DuckDB threads of their own.
# Integration targets get their own (smaller) pool: each worker starts a private PostgreSQL.
PYTEST_WORKERS ?= 4
PYTEST_PG_WORKERS ?= 2
PYTEST_ARGS ?=
SUBSET_PRODUCTS ?= exclusion_metabolites metsigdb network_views cosmos
API_PORT ?= 8085
WEB_PORT ?= 5173
COMPOSE_FILES ?= -f compose.serving.yaml
COMPOSE_ENV ?=
COMPOSE = docker compose $(if $(COMPOSE_ENV),--env-file "$(COMPOSE_ENV)") $(COMPOSE_FILES)

help:
	@echo 'setup | sample | build | hubs | reference | publish-release'
	@echo 'dev | api | web | serving-build | serving-up | serving-status | serving-logs | serving-stop'
	@echo 'load-postgres | finish-postgres | build-subsets | test | check | check-web | docs'
	@echo 'See README.md for variables and deploy/README.md for container deployment.'

setup: setup-python setup-web
setup-python:
	git submodule update --init pypath
	uv sync --frozen --all-packages
setup-web:
	pnpm --dir packages/web install --frozen-lockfile

test:
	uv run --frozen pytest -m 'not integration' -n $(PYTEST_WORKERS) $(PYTEST_ARGS)
# Offline inputs_v2 tests maintained with this workspace (no downloads).
PYPATH_TESTS = pypath/test/test_chembl_molecular_queries.py \
	pypath/test/test_cv_terms.py \
	pypath/test/test_download_timeouts.py \
	pypath/test/test_input_source_context.py \
	pypath/test/test_inputs_v2_chemical_migration.py \
	pypath/test/test_inputs_v2_communication_migration.py \
	pypath/test/test_inputs_v2_ontology_migration.py \
	pypath/test/test_inputs_v2_pathway_migration.py \
	pypath/test/test_interaction_profiles.py \
	pypath/test/test_measurement_review_fixes.py \
	pypath/test/test_model_rule_review_fixes.py \
	pypath/test/test_molecular_catalogue_coverage.py \
	pypath/test/test_molecular_forms.py \
	pypath/test/test_molecular_sequence_coverage.py \
	pypath/test/test_reactome_molecular_coverage.py \
	pypath/test/test_resource_representation_review.py \
	pypath/test/test_signor_identifiers.py \
	pypath/test/test_tabular_execution.py
test-pypath:
	uv run --frozen pytest $(PYPATH_TESTS) -q
test-postgres:
	OMNIPATH_TEST_POSTGRES=1 uv run --frozen pytest packages/postgres/tests -q -n $(PYTEST_PG_WORKERS) $(PYTEST_ARGS)
test-subsets:
	OMNIPATH_TEST_POSTGRES=1 uv run --frozen pytest packages/postgres/tests packages/subsets/tests -q -n $(PYTEST_PG_WORKERS) $(PYTEST_ARGS)
test-api:
	uv run --frozen pytest packages/api/tests -q -n $(PYTEST_WORKERS) $(PYTEST_ARGS)
check:
	uv run --frozen ruff check packages scripts conftest.py
	uv run --frozen ruff format --check packages scripts conftest.py
check-web:
	pnpm --dir packages/web check
	pnpm --dir packages/web test
	pnpm --dir packages/web lint

# Resource processing is bounded unless explicitly overridden.
build:
	uv run --frozen --package omnipath-build omnipath-build build "$(SOURCE)" --version "$(VERSION)" --max-records $(MAX_RECORDS) --output-dir "$(DATA_ROOT)"
sample:
	uv run --frozen --package omnipath-build python scripts/migration_smoke.py --output-dir "$(DATA_ROOT)/migration-smoke" --version "$(VERSION)" --max-records $(MAX_RECORDS)
hubs:
	uv run --frozen --package omnipath-build omnipath-build export-hubs --output-dir "$(HUBS_DIR)" $(HUB_ARGS)
# Hub indexes are reused while their hub export is unchanged; build-identity also writes its kv store.
reference:
	@test -n "$(IDENTITY_HUBS)" || (echo 'Cannot list the identity hubs'; exit 2)
	for hub in $(IDENTITY_HUBS); do \
		uv run --frozen --package omnipath-build omnipath-build build-hub-index --hub $$hub --hubs-dir "$(HUBS_DIR)" --output-root "$(HUB_INDEX_ROOT)" $(HUB_INDEX_ARGS) || exit 1; \
	done
	for hub in $(filter-out ontology,$(IDENTITY_HUBS)); do \
		uv run --frozen --package omnipath-build omnipath-build build-hub-kv --hub $$hub --hub-index-root "$(HUB_INDEX_ROOT)" $(HUB_KV_ARGS) || exit 1; \
	done
	uv run --frozen --package omnipath-build omnipath-build build-identity --hub-index-root "$(HUB_INDEX_ROOT)" --output-dir "$(IDENTITY_DIR)" $(IDENTITY_ARGS)
publish-release:
	@test -n "$(RELEASE_MANIFEST)" || (echo 'Set RELEASE_MANIFEST to an explicit pinned JSON file'; exit 2)
	uv run --frozen --package omnipath-api python -m omnipath_api.publish_release "$(RELEASE_MANIFEST)" --data-root "$(DATA_ROOT)"

api:
	uv run --frozen --package omnipath-api omnipath-api --host 127.0.0.1 --port $(API_PORT) --data-root "$(DATA_ROOT)"
web:
	API_SERVICE_URL=http://127.0.0.1:$(API_PORT) pnpm --dir packages/web dev --host 127.0.0.1 --port $(WEB_PORT)
dev:
	OMNIPATH_DATA_ROOT="$(DATA_ROOT)" ./run_all.sh $(API_PORT) $(WEB_PORT)

load-postgres:
	@test -n "$(RELEASE_MANIFEST)" || (echo 'Set RELEASE_MANIFEST to a pinned JSON file or HTTPS URL'; exit 2)
	uv run --frozen --package omnipath-postgres omnipath-postgres "$(RELEASE_MANIFEST)" --data-root "$(DATA_ROOT)" --schema "$(POSTGRES_SCHEMA)" $(POSTGRES_ARGS)
finish-postgres:
	uv run --frozen --package omnipath-postgres omnipath-postgres --finish --schema "$(POSTGRES_SCHEMA)" $(POSTGRES_ARGS)
# Use the aligned PostgreSQL pipeline; durable completed phases are skipped.
build-subsets:
	uv run --frozen --package omnipath-postgres omnipath-postgres --finish --schema "$(POSTGRES_SCHEMA)" --products $(SUBSET_PRODUCTS)

serving-build:
	$(COMPOSE) build api web
serving-up:
	$(COMPOSE) up -d --wait api web data
serving-status:
	$(COMPOSE) ps
serving-logs:
	$(COMPOSE) logs --tail 50 api web data
serving-stop:
	$(COMPOSE) stop api web data
# Only this explicit target starts resource-job processing.
worker-up:
	$(COMPOSE) --profile build up -d --build worker

.PHONY: setup-serving setup-build setup-postgres test-native check-generated generate test-wheels postgres-up rebuild-subsets
setup-serving:
	uv sync --frozen --package omnipath-api
	$(MAKE) setup-web
setup-build:
	git submodule update --init pypath
	uv sync --frozen --package omnipath-build
setup-postgres:
	git submodule update --init pypath
	uv sync --frozen --package omnipath-postgres
test-native:
	cargo test --locked --manifest-path packages/resolver/rust/reference/Cargo.toml
check-generated:
	uv run --frozen python scripts/generate_openapi.py --check
	uv run --frozen python scripts/sync_biolink.py --check
	uv run --frozen python scripts/generate_resource_metadata.py --check
	uv run --frozen python packages/query/scripts/generate_contract.py --check
generate:
	uv run --frozen python scripts/generate_openapi.py
	uv run --frozen python scripts/sync_biolink.py
	uv run --frozen python scripts/generate_resource_metadata.py
	uv run --frozen python packages/query/scripts/generate_contract.py
	pnpm --dir packages/web generate:types
test-wheels:
	uv run --frozen python scripts/check_wheels.py
POSTGRES_ENV ?= deploy/postgres.env
postgres-up:
	docker compose --env-file "$(POSTGRES_ENV)" -f docker-compose.postgres18.yml up -d --build --wait
rebuild-subsets:
	uv run --frozen --package omnipath-subsets omnipath-subsets build --schema "$(POSTGRES_SCHEMA)" --products $(SUBSET_PRODUCTS) --checkpoint-products
# Core documentation viewer (Markdown sources in core_documentation/).
DOCS_PORT ?= 8090
docs:
	uv run --frozen python -m http.server $(DOCS_PORT) --bind 127.0.0.1 --directory core_documentation
