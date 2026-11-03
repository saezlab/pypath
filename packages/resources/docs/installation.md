# Installation

## From PyPI

```bash
pip install omnipath-resources
```

Or with [uv](https://docs.astral.sh/uv/):

```bash
uv pip install omnipath-resources
```

## Optional dependencies

### curl support

For HTTP downloads via pycurl:

```bash
pip install "omnipath-resources[curl]"
```

### Visualization

For plotting capabilities:

```bash
pip install "omnipath-resources[vis]"
```

### Graph analysis

For igraph-based network analysis:

```bash
pip install "omnipath-resources[graph]"
```

### Metabolomics

For metabolite structure handling:

```bash
pip install "omnipath-resources[metabo]"
```

## Development installation

```bash
git clone https://github.com/saezlab/pypath.git
cd pypath
uv sync --extra dev --extra tests
```
