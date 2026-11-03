# Installation

## omnipath-client (recommended for most users)

If you only need to translate IDs, resolve organisms, or check reference
lists, install the client. It queries the
[utils.omnipathdb.org](https://utils.omnipathdb.org) web service and requires
no local database:

```bash
pip install omnipath-client
```

```python
from omnipath_client.utils import map_name, ensure_ncbi_tax_id
map_name('TP53', 'genesymbol', 'uniprot')  # {'P04637'}
ensure_ncbi_tax_id('human')  # 9606
```

See the [omnipath-client docs](https://saezlab.github.io/omnipath-client/)
for details.

## omnipath-utils (local library / server)

Install `omnipath-utils` directly if you need to run the service locally,
build the database yourself, or develop backends.

### Basic (Python API only)

```bash
pip install omnipath-utils
```

### With database support

```bash
pip install "omnipath-utils[db]"
```

### With web service

```bash
pip install "omnipath-utils[server]"
```

### With pypath backends (recommended for building)

```bash
pip install "omnipath-utils[pypath]"
```

### Development

```bash
git clone https://github.com/saezlab/omnipath-utils
cd omnipath-utils
uv sync --all-extras
uv run pytest
```
