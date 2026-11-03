# omnipath-utils

!!! tip "For most users: use omnipath-client"

    If you just want to translate IDs, resolve organisms, or access reference
    lists from Python, install the lightweight client instead:

    ```bash
    pip install omnipath-client
    ```

    ```python
    from omnipath_client.utils import map_name, translate_column
    map_name('TP53', 'genesymbol', 'uniprot')  # {'P04637'}
    ```

    The client queries [utils.omnipathdb.org](https://utils.omnipathdb.org) --
    no database setup required, same API as the local package.

ID translation, taxonomy, orthology, and reference lists for molecular biology.

omnipath-utils provides a unified API for translating between biological
identifiers (UniProt, gene symbols, Ensembl, HGNC, etc.), resolving organism
names, and accessing reference lists. It works as a standalone Python library
or as an HTTP web service backed by PostgreSQL.

## Features

- **ID Translation** -- Translate between 97 identifier types across 7+ backends (UniProt, Ensembl BioMart, HMDB, RaMP, etc.)
- **Taxonomy** -- Resolve organism names across naming systems (NCBI, Ensembl, KEGG, OMA, miRBase)
- **Reference Lists** -- Complete sets of identifiers (all human SwissProt IDs, etc.)
- **Web Service** -- Litestar HTTP API with auto-generated OpenAPI schema
- **Database** -- PostgreSQL-backed for fast lookups at scale
