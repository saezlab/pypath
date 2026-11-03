# CLAUDE.md - PyPath Development Guide

## Project Overview

PyPath is a Python library for molecular biology and biomedical prior knowledge processing. It aggregates data from ~200 biological databases into unified structures for network biology, annotations, protein complexes, and enzyme-substrate relationships.

**Key facts:**
- 12+ years of development history with heterogeneous code quality
- Python 3.9+ (dropped Python 2 in 2021)
- Poetry for dependency management
- GPLv3 license
- Ongoing modularization into separate packages

## Project Structure

```
pypath/
├── core/          # Database classes (network, complex, annot, enz_sub, intercell)
├── inputs/        # Data acquisition: 189 modules for ~150 data sources
├── share/         # Shared utilities (curl, session/logging, settings, cache)
├── utils/         # Standalone utilities (mapping, taxonomy, orthology, GO, sequences)
├── resources/     # Data resources (urls.py with 1000+ endpoints, descriptions)
├── formats/       # Format parsers (sqlite, xml, sdf, obo, sqldump)
├── internals/     # Internal infrastructure (input_formats, resource definitions)
├── omnipath/      # High-level database management, web server, export
├── legacy/        # Deprecated code (do not use)
└── obsolete/      # Old code (do not use)
```

**External companion packages:**
- `pypath-common` - Shared utilities across pypath projects
- `download-manager` - Generic download infrastructure
- `cachemanager` - Cache management

## Code Style and Conventions

### Naming

- **Classes**: `PascalCase`. Resource names as single words: `ProteinatlasAnnotation`
- **Functions/methods/variables**: `snake_case`
- **Resource names in functions**: Single word, no underscores: `proteinatlas_annotations()` not `protein_atlas_annotations()`
  - Underscores separate primary/secondary sources: `PhosphoSite_ProtMapper`

### Input Function Naming Pattern

```
{resource}_{suffix}
```

Standard suffixes:
- `_interactions` - Protein-protein interactions
- `_enz_sub` - Enzyme-substrate relationships
- `_complexes` - Protein complex membership
- `_annotations` - Functional annotations
- `_raw` - Near-raw data dumps
- `_mapping` - ID cross-references

### Formatting

- **Two blank lines** between methods and classes
- **One blank line** before/after logic blocks, between logical segments
- **Single quotes** preferred unless reason for double quotes
- **Operators** surrounded by spaces: `a = a * 4 + 3`

**Argument lists:**
```python
def function(
        arg1: str,
        arg2: int,
        **kwargs,
    ) -> list:
```
- New line after opening parenthesis
- Each element on own line, indented one level
- Trailing comma after each element
- Closing parenthesis on own line at original indentation

**Long comprehensions:**
```python
result = [
    transform(item)
    for item in items
    if condition(item)
]
```

### Docstrings

Google (Napoleon) style with type hints:
```python
def function(a: int) -> list[int]:
    """
    Brief description of function.

    Args:
        a: A number. Description continues
            on next line with indent.

    Returns:
        A list of integers.
    """
```

### Imports

Order (with blank lines between sections):
1. `from __future__ import annotations`
2. Standard library
3. Typing imports
4. Third-party packages
5. PyPath internal imports

```python
from __future__ import annotations

import collections
import re

from typing import Iterable, Optional

import pandas as pd

import pypath.share.curl as curl
import pypath.utils.mapping as mapping
```

**Avoid in new code** (legacy compatibility):
```python
# DON'T use these in new code:
from future.utils import iteritems  # Use .items() directly
from past.builtins import xrange    # Use range() directly
```

### Logging

Classes inherit from `pypath.share.session.Logger`:
```python
class MyClass(session_mod.Logger):
    def __init__(self):
        session_mod.Logger.__init__(self, name='mymodule')
        self._log('Initialized')
```

Or module-level logger:
```python
import pypath.share.session as session_mod
_logger = session_mod.Logger(name='inputs.mymodule')
_logger._log('Processing data...')
```

## The `inputs` Module

### Structure

- **189 Python files** for ~150 data sources
- **16 sub-packages** for complex resources (chembl, brenda, bindingdb, ramp, reactome, rhea, disgenet, hmdb, lipidmaps, etc.)
- Simple resources: single file `inputs/{resource}.py`
- Complex resources: sub-package `inputs/{resource}/`

### Standard Data Flow

1. **Download** via `pypath.share.curl.Curl`
2. **Parse** to Python structures (lists, dicts, namedtuples)
3. **Map IDs** to UniProt using `pypath.utils.mapping`
4. **Return** standardized format

### Basic Input Module Pattern

```python
"""
Client and parser for ResourceName database.
"""

from __future__ import annotations

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


ResourceRecord = collections.namedtuple(
    'ResourceRecord',
    ['source', 'target', 'effect'],
)


def resource_interactions() -> list[ResourceRecord]:
    """
    Retrieves interactions from ResourceName.

    Returns:
        List of interaction records.
    """

    url = urls.urls['resource']['url']
    c = curl.Curl(url, silent=False, large=True)

    result = []

    for line in c.result.split('\n')[1:]:  # Skip header

        fields = line.strip().split('\t')
        if len(fields) < 3:
            continue

        result.append(
            ResourceRecord(
                source=fields[0],
                target=fields[1],
                effect=fields[2],
            )
        )

    return result
```

### Complex Resource Sub-Package Pattern

```
inputs/resource/
├── __init__.py      # Public API exports
├── _common.py       # Shared utilities and logger
├── _sqlite.py       # SQLite database access (if applicable)
├── _molecules.py    # Specific data extractors
├── _targets.py
└── _raw.py          # Raw data access
```

### SQLite Database Pattern

```python
import pypath.formats.sqlite as _sqlite
import pypath.share.curl as curl
import pypath.resources.urls as urls


def resource_sqlite(version: int = 1, connect: bool = True):
    """
    Downloads and connects to ResourceName SQLite database.
    """

    url = urls.urls['resource']['sqlite'] % version

    def _download() -> curl.Curl:
        return curl.Curl(url, large=True, silent=False)

    def _extractor(curl_obj: curl.Curl) -> str:
        return curl_obj.result

    return _sqlite.download_sqlite(
        download_callback=_download,
        extractor=_extractor,
        database='ResourceName',
        version=str(version),
        connect=connect,
    )
```

### Common Utilities

**ID Mapping:**
```python
import pypath.utils.mapping as mapping

# Map multiple IDs
uniprot_ids = mapping.map_names(gene_symbols, 'genesymbol', 'uniprot')

# Map single ID
uniprot = mapping.map_name('EGFR', 'genesymbol', 'uniprot')
```

**Taxonomy:**
```python
import pypath.utils.taxonomy as taxonomy

name = taxonomy.taxids[9606]  # 'Homo sapiens'
```

**pypath_common utilities:**
```python
from pypath_common import _misc, _constants as _const

first_item = _misc.first(iterable)
```

**Glom for JSON processing:**
```python
import glom
from pypath_common import _constants as _const

spec = {'field': 'nested.path.to.field'}
result = glom.glom(record, spec, default=_const.GLOM_ERROR)
```

## Cache and Downloads

**Cache location:** `~/.pypath/cache/` (hash-named files)

**Cache control context managers:**
```python
import pypath.share.curl as curl

# Force fresh download
with curl.cache_delete_on():
    data = fetch_data()

# Enable caching
with curl.cache_on():
    data = fetch_data()

# Debug mode
with curl.debug_on():
    data = fetch_data()
```

## Testing Input Modules

```bash
# Test specific module
python input_module_maintenance/test_input_modules.py --module signor

# Test specific function
python input_module_maintenance/test_input_modules.py --function signor.signor_interactions

# List all modules
python input_module_maintenance/test_input_modules.py --list-modules
```

**Quick smoke test:**
```python
from pypath.inputs import signor
ints = signor.signor_interactions()
print(len(ints), ints[0])
```

## Adding a New Resource

1. Create module in `pypath/inputs/` (file or sub-package)
2. Add URL to `pypath/resources/urls.py`
3. Add metadata to `pypath/resources/descriptions.py`
4. Update `resources.json` with license information
5. Create core database input definitions if needed:
   - `network`: Format definition in `resources/data_formats.py`
   - `annot`: Resource class in `core/annot.py`
   - `complex`: Resource class in `core/complex.py`
   - `enz_sub`: Definition in `resources.json`
6. Write tests
7. Update documentation

## Common Gotchas

### Network Issues
- SSL certificate errors common - use `curl.debug_on()` context
- Clear corrupted cache with `curl.cache_delete_on()`
- Some resources have rate limiting or require retries

### ID Mapping
- Species parameter is critical: `ncbi_tax_id=9606` for human
- Old/deleted UniProt IDs may fail to map
- Always validate mapping results

### Format Drift
- Data sources change formats without warning
- Compare cached vs fresh downloads when debugging
- Manual browser verification needed for failed tests

### Code Quality Variance
- ~30% of codebase is older "messy" style from pre-2019
- Follow modern conventions in new code
- Don't fix unrelated style issues in the same PR

### Large Datasets
- Some resources timeout during download
- May need extended timeout configuration
- Use `large=True` parameter in Curl for big files

## Key Files Reference

- `pypath/resources/urls.py` - All 1000+ data source URLs
- `pypath/resources/descriptions.py` - Resource metadata
- `pypath/resources/data_formats.py` - Input format definitions
- `pypath/share/curl.py` - Download infrastructure
- `pypath/utils/mapping.py` - ID translation
- `CONTRIBUTING.md` - Detailed contribution guide
- `input_module_maintenance/` - Testing and maintenance tools
