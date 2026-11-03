.PHONY: setup resolver db-setup db-reset reset-content drop-source derive load reload pipeline all test

INSTANCE_ENV_FILE ?= ../.env
ifneq ($(wildcard $(INSTANCE_ENV_FILE)),)
include $(INSTANCE_ENV_FILE)
endif

export PKG_INFRA_CONFIG ?= $(CURDIR)/config/pkg_infra_quiet.yaml

DATA_ROOT ?= $(if $(DATA_DIR),$(DATA_DIR),data)
DATABASE ?= omnipath
DEFAULT_DATABASE_URL ?= postgresql://omnipath:omnipath@localhost:55432/omnipath
ifneq ($(POSTGRES_PORT),)
DEFAULT_DATABASE_URL := postgresql://omnipath:omnipath@localhost:$(POSTGRES_PORT)/omnipath
else ifneq ($(DATABASE_URL),)
DEFAULT_DATABASE_URL := $(DATABASE_URL)
endif

ifeq ($(origin DATABASE_URL),command line)
BUILD_DATABASE_URL ?= $(DATABASE_URL)
else
BUILD_DATABASE_URL ?= $(DEFAULT_DATABASE_URL)
DATABASE_URL := $(BUILD_DATABASE_URL)
endif
SCHEMA ?= public
INPUTS_PACKAGE ?= pypath.inputs_v2
SOURCE ?=
SOURCES ?=
DATASET ?=
source ?=
sources ?=
SELECTED_SOURCES = $(strip $(if $(SOURCES),$(SOURCES),$(if $(SOURCE),$(SOURCE),$(if $(sources),$(sources),$(source)))))
LOAD_SOURCES = $(strip $(SELECTED_SOURCES))

MAPPING_DIR ?= $(DATA_ROOT)
RESOLVER_SOURCES ?=
BATCH_SIZE ?= 50000
THREADS ?= 4
LOAD_JOBS ?= 1
JOBS ?= 1
STAGING_DIR ?=
RESOLVER_BATCH_SIZE ?= 100000
DROP_EXISTING ?=
RESET_DROP_INDEXES ?=
DERIVE ?=
PROFILE ?=
MAX_RECORDS ?=
FORCE_REFRESH ?=
RELOAD_EXISTING ?=
PUBCHEM_URL ?=
PUBCHEM_SHARDS ?=
DUCKDB_MAX_RECORDS_ARG = $(if $(MAX_RECORDS),--max-records "$(MAX_RECORDS)",--max-records 0)

setup:
	git submodule add -b main https://github.com/saezlab/pypath.git pypath || true
	git submodule add -b main https://github.com/saezlab/download-manager.git download-manager || true
	git submodule update --init --recursive --remote
	uv sync
	-uv pip uninstall pypath-omnipath dlmachine cachedir
	uv pip install -e ./cache-manager -e ./download-manager -e ./pypath

resolver:
	@echo "[omnipath_build] build-resolver output_dir=$(MAPPING_DIR) sources=$(if $(RESOLVER_SOURCES),$(RESOLVER_SOURCES),ALL)"
	@PYTHONUNBUFFERED=1 uv run python -m omnipath_build.cli build-resolver \
		--output-dir "$(MAPPING_DIR)" \
		$(if $(MAX_RECORDS),--max-records $(MAX_RECORDS)) \
		$(if $(PUBCHEM_URL),--pubchem-url "$(PUBCHEM_URL)") \
		$(if $(PUBCHEM_SHARDS),--pubchem-shards $(PUBCHEM_SHARDS)) \
		$(if $(JOBS),--jobs $(JOBS)) \
		$(if $(FORCE_REFRESH),--no-skip-existing) \
		$(RESOLVER_SOURCES)

db-setup:
	@if [ -z "$(DATABASE_URL)" ]; then \
		echo "DATABASE_URL is required, e.g. make db-setup DATABASE_URL=postgresql://user:pass@host:5432/dbname"; \
		exit 1; \
	fi
	@set -e; \
	echo "[omnipath_build] init-db schema=$(SCHEMA)"; \
	PYTHONUNBUFFERED=1 uv run python -m omnipath_build.cli \
		--database-url "$(DATABASE_URL)" \
		--schema "$(SCHEMA)" \
		init-db \
		$(if $(DROP_EXISTING),--drop-existing --no-indexes); \
	if [ -z "$(DROP_EXISTING)" ]; then \
		echo "[omnipath_build] create secondary indexes schema=$(SCHEMA)"; \
		PYTHONUNBUFFERED=1 uv run python -m omnipath_build.cli \
			--database-url "$(DATABASE_URL)" \
			--schema "$(SCHEMA)" \
			derive \
			--no-tables \
			--no-bitmaps; \
	else \
		echo "[omnipath_build] deferring secondary indexes until after load"; \
	fi

db-reset:
	@if [ -z "$(DATABASE_URL)" ]; then \
		echo "DATABASE_URL is required, e.g. make db-reset DATABASE_URL=postgresql://user:pass@host:5432/dbname"; \
		exit 1; \
	fi
	@echo "[omnipath_build] reset database schema=$(SCHEMA)"
	@PYTHONUNBUFFERED=1 uv run python -m omnipath_build.cli \
		--database-url "$(DATABASE_URL)" \
		--schema "$(SCHEMA)" \
		init-db \
		--drop-existing

reset-content:
	@if [ -z "$(DATABASE_URL)" ]; then \
		echo "DATABASE_URL is required, e.g. make reset-content DATABASE_URL=postgresql://user:pass@host:5432/dbname"; \
		exit 1; \
	fi
	@echo "[omnipath_build] reset content tables schema=$(SCHEMA)"
	@PYTHONUNBUFFERED=1 uv run python -m omnipath_build.cli \
		--database-url "$(DATABASE_URL)" \
		--schema "$(SCHEMA)" \
		reset-content \
		$(if $(RESET_DROP_INDEXES),--drop-indexes)

drop-source:
	@if [ -z "$(DATABASE_URL)" ]; then \
		echo "DATABASE_URL is required, e.g. make drop-source SOURCE=uniprot DATABASE_URL=postgresql://user:pass@host:5432/dbname"; \
		exit 1; \
	fi
	@if [ -z "$(SELECTED_SOURCES)" ]; then \
		echo "SOURCE or SOURCES is required, e.g. make drop-source SOURCE=uniprot"; \
		exit 1; \
	fi
	@set -e; \
	SOURCE_LIST=$$(printf '%s' "$(SELECTED_SOURCES)" | tr ',' ' '); \
	for source in $$SOURCE_LIST; do \
		echo "[omnipath_build] drop-source source=$$source schema=$(SCHEMA)"; \
		PYTHONUNBUFFERED=1 uv run python -m omnipath_build.cli \
			--database-url "$(DATABASE_URL)" \
			--schema "$(SCHEMA)" \
			drop-source \
			--source "$$source"; \
	done

derive:
	@if [ -z "$(DATABASE_URL)" ]; then \
		echo "DATABASE_URL is required, e.g. make derive DATABASE_URL=postgresql://user:pass@host:5432/dbname"; \
		exit 1; \
	fi
	@echo "[omnipath_build] refresh derived tables and bitmaps schema=$(SCHEMA)"
	@PYTHONUNBUFFERED=1 uv run python -m omnipath_build.cli \
		--database-url "$(DATABASE_URL)" \
		--schema "$(SCHEMA)" \
		derive \
		$(if $(MAX_RECORDS),--max-records "$(MAX_RECORDS)")

network-views:
	@if [ -z "$(DATABASE_URL)" ]; then \
		echo "DATABASE_URL is required, e.g. make network-views DATABASE_URL=postgresql://user:pass@host:5432/dbname"; \
		exit 1; \
	fi
	@echo "[omnipath_build] apply specialized network views schema=$(SCHEMA)"
	@PYTHONUNBUFFERED=1 uv run python -m omnipath_build.cli \
		--database-url "$(DATABASE_URL)" \
		--schema "$(SCHEMA)" \
		network-views \
		$(if $(REFRESH),--refresh)

load:
	@if [ -z "$(DATABASE_URL)" ]; then \
		echo "DATABASE_URL is required, e.g. make load SOURCE=uniprot DATABASE_URL=postgresql://user:pass@host:5432/dbname"; \
		exit 1; \
	fi
	@echo "[omnipath_build] load sources=$(LOAD_SOURCES) schema=$(SCHEMA) batch_size=$(BATCH_SIZE) threads=$(THREADS) load_jobs=$(LOAD_JOBS)"
	@PYTHONUNBUFFERED=1 uv run python -m omnipath_build.duckdb_direct_pipeline \
		--database-url "$(DATABASE_URL)" \
		--schema "$(SCHEMA)" \
		--database "$(DATABASE)" \
		--inputs-package "$(INPUTS_PACKAGE)" \
		--sources "$(LOAD_SOURCES)" \
		--resolver-dir "$(MAPPING_DIR)" \
		$(if $(DATASET),--dataset "$(DATASET)") \
		$(DUCKDB_MAX_RECORDS_ARG) \
		--batch-size "$(BATCH_SIZE)" \
		--threads "$(THREADS)" \
		--stage-jobs "$(LOAD_JOBS)" \
		$(if $(STAGING_DIR),--staging-dir "$(STAGING_DIR)") \
		$(if $(FORCE_REFRESH),--force-refresh) \
		$(if $(RELOAD_EXISTING),--reload-existing) \
		--append

reload: RELOAD_EXISTING=1
reload: load

pipeline: load
	@if [ "$(DERIVE)" != "" ]; then \
		$(MAKE) derive; \
	fi

all:
	@$(MAKE) resolver
	@$(MAKE) db-setup
	@$(MAKE) pipeline

test:
	@uv run python -m compileall -q omnipath_build
