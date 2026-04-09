# Installation

## From PyPI

```bash
pip install pypath-omnipath
```

Or with [uv](https://docs.astral.sh/uv/):

```bash
uv pip install pypath-omnipath
```

## Optional dependencies

### curl support

For HTTP downloads via pycurl:

```bash
pip install "pypath-omnipath[curl]"
```

### Visualization

For plotting capabilities:

```bash
pip install "pypath-omnipath[vis]"
```

### Graph analysis

For igraph-based network analysis:

```bash
pip install "pypath-omnipath[graph]"
```

### Metabolomics

For metabolite structure handling:

```bash
pip install "pypath-omnipath[metabo]"
```

## Development installation

```bash
git clone https://github.com/saezlab/pypath.git
cd pypath
uv sync --extra dev --extra tests
```
