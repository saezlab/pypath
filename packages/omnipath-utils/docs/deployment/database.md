# Database Build

## Prerequisites

- PostgreSQL 14+
- Python with `omnipath-utils[db]` installed

## Setup

```bash
# Start a PostgreSQL instance (example with Docker)
docker run -d --name omnipath-db -p 5433:5432 \
  -e POSTGRES_PASSWORD=dev -e POSTGRES_DB=omnipath_utils \
  postgres:16

# Set the connection URL
export OMNIPATH_UTILS_DB_URL="postgresql+psycopg://postgres:dev@localhost:5433/omnipath_utils"
```

## Build

```bash
# Build reference tables + human mappings
omnipath-utils build --organisms 9606

# Build for multiple organisms
omnipath-utils build --organisms 9606 10090 10116

# Build only reference tables (fast)
omnipath-utils build --ref-only
```

## What gets built

| Table | Content | Human rows |
|-------|---------|-----------|
| id_type | 97 identifier types with metadata | 97 |
| backend | 11 data source backends | 11 |
| organism | 22 organisms with all name forms | 22 |
| id_mapping | ID translation pairs | ~165K per mapping type |
| build_info | Build audit log | 1 per build target |

## Schema

All tables are created in the `omnipath_utils` PostgreSQL schema.
The build process:

1. Creates the schema if it does not exist
2. Creates all tables via SQLAlchemy ORM
3. Populates reference tables (id_type, backend, organism) from YAML files
4. For each organism, loads ID mapping pairs from upstream APIs
   and bulk-inserts via PostgreSQL COPY

## Programmatic usage

```python
from omnipath_utils.db._build import DatabaseBuilder

builder = DatabaseBuilder(db_url='postgresql+psycopg://postgres:dev@localhost:5433/omnipath_utils')

# Build everything for human
builder.build_all(organisms=[9606])

# Or just reference tables
builder.build_reference_tables()

# Or a single mapping
builder.populate_mapping('genesymbol', 'uniprot', 9606, 'uniprot')
```
