# Web Service

omnipath-utils includes a Litestar HTTP API for serving ID translation,
taxonomy, and reference list queries.

## Starting the server

```bash
# Ensure the database is built first
omnipath-utils build --organisms 9606

# Start the web service
omnipath-utils serve --port 8082
```

Options:

- `--host` -- bind address (default: `0.0.0.0`)
- `--port` -- listen port (default: `8082`)
- `--db-url` -- PostgreSQL connection URL (or set `OMNIPATH_UTILS_DB_URL`)
- `-v` -- verbose logging

## Endpoints

### Health check

```
GET /health
```

Returns `{"status": "ok", "service": "omnipath-utils"}`.

### Taxonomy

```
GET /taxonomy/resolve?organism=human
```

Resolves an organism name to all name forms (NCBI tax ID, common name,
latin name, KEGG code, Ensembl name, etc.).

```
GET /taxonomy/organisms
```

Lists all known organisms.

### ID Translation

```
GET /mapping/translate?identifiers=TP53,EGFR&id_type=genesymbol&target_id_type=uniprot
```

Translates comma-separated identifiers. Optional `ncbi_tax_id` parameter
(default: 9606).

```
POST /mapping/translate
Content-Type: application/json

{
  "identifiers": ["TP53", "EGFR"],
  "id_type": "genesymbol",
  "target_id_type": "uniprot",
  "ncbi_tax_id": 9606
}
```

Response:

```json
{
  "results": {"TP53": ["P04637"], "EGFR": ["P00533"]},
  "unmapped": [],
  "meta": {
    "id_type": "genesymbol",
    "target_id_type": "uniprot",
    "ncbi_tax_id": 9606,
    "total_input": 2,
    "total_mapped": 2
  }
}
```

```
GET /mapping/id-types
```

Lists all 97 supported identifier types with metadata.

### Reference Lists

```
GET /reflists/list-names
```

Returns available reference list names.

```
GET /reflists/{list_name}?ncbi_tax_id=9606
```

Returns identifiers in a reference list.

## OpenAPI schema

The server auto-generates an OpenAPI schema at `/schema/openapi.json`
(or explore interactively via the Swagger UI if enabled).

## Programmatic usage

```python
from omnipath_utils.server._app import create_app

app = create_app(db_url='postgresql+psycopg://postgres:dev@localhost:5433/omnipath_utils')
```
