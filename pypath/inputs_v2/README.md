# inputs_v2 - Declarative Data Input Modules

This module provides a declarative framework for downloading, parsing, and normalizing biological data from external resources into a unified `Entity` schema.

## Overview

The `inputs_v2` system replaces the legacy input modules with a composable, declarative approach. Each data source is defined as a **Resource** containing one or more **Datasets**, which combine:

1. **Download** - Configuration for fetching remote data
2. **Raw Parser** - Function to convert raw files into dictionaries
3. **Schema (EntityBuilder)** - Declarative mapping from dictionaries to `Entity` records

## Directory Structure

```
inputs_v2/
├── __init__.py           # Module initialization and get_method() helper
├── base.py               # Core building blocks (Resource, Dataset, Download, ResourceConfig)
├── parsers/              # Raw data parsers for various file formats
│   ├── base.py           # Common parsers: iter_csv, iter_tsv, iter_json, iter_jsonl
│   ├── hmdb.py           # HMDB XML parser
│   ├── lipidmaps.py      # LIPID MAPS SDF parser
│   ├── reactome.py       # Reactome BioPAX parser
│   └── bindingdb.py      # BindingDB parser
├── ontologies/           # Ontology-specific modules
│   ├── shared.py         # Common OBO parsing utilities
│   ├── gene_ontology.py  # Gene Ontology
│   ├── psi_mi.py         # PSI-MI ontology
│   ├── omnipath.py       # OmniPath internal ontology
│   └── uniprot_keywords.py  # UniProt Keywords
├── bindingdb.py          # BindingDB interactions & compounds
├── guidetopharma.py      # Guide to Pharmacology
├── hmdb.py               # Human Metabolome Database
├── intact.py             # IntAct interactions
├── lipidmaps.py          # LIPID MAPS Structure Database
├── reactome.py           # Reactome pathways & reactions
├── signor.py             # SIGNOR signaling network
├── swisslipids.py        # SwissLipids
└── uniprot.py            # UniProt proteins
```

## Core Concepts

### 1. ResourceConfig

Defines metadata about a data source:

```python
from pypath.inputs_v2.base import ResourceConfig
from pypath.internals.cv_terms import ResourceCv, LicenseCV, UpdateCategoryCV

config = ResourceConfig(
    id=ResourceCv.UNIPROT,           # CV term identifying the resource
    name='UniProt',                   # Human-readable name
    url='https://www.uniprot.org/',   # Resource URL
    license=LicenseCV.CC_BY_4_0,      # License type
    update_category=UpdateCategoryCV.REGULAR,  # How often data updates
    pubmed='33237286',                # Optional PubMed reference
    description='UniProt is a comprehensive resource for protein...',
)
```

### 2. Download

Configures how to fetch remote data:

```python
from pypath.inputs_v2.base import Download

download = Download(
    url='https://example.com/data.tsv',  # URL (can be a callable for dynamic URLs)
    filename='data.tsv',                  # Local filename
    subfolder='example',                  # Cache subfolder
    large=True,                           # Whether to show progress
    encoding='utf-8',                     # File encoding
    default_mode='r',                     # File open mode
    ext='gz',                             # Optional: file extension hint
    download_kwargs={'post': True},       # Optional: additional download options
)
```

### 3. Raw Parsers

Functions that convert raw file data into dictionaries. Standard parsers are available in `parsers/base.py`:

```python
from pypath.inputs_v2.parsers.base import iter_tsv, iter_csv, iter_json, iter_jsonl

# Signature: (opener, **kwargs) -> Generator[dict[str, Any], None, None]
def custom_parser(opener, **_kwargs):
    """Parse a custom file format."""
    if not opener or not opener.result:
        return
    for line in opener.result:
        yield {'field': line.strip()}
```

### 4. Schema Definition (EntityBuilder)

Declaratively map dictionary fields to Entity structure using `FieldConfig` and `CV`:

```python
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import EntityTypeCv, IdentifierNamespaceCv

# Create a field configuration (reusable across multiple CVs)
f = FieldConfig(
    extract={
        'chebi': r'^(?:CHEBI:)?(\d+)$',  # Named regex patterns
    },
    transform={
        'chebi': lambda v: f'CHEBI:{v}' if v else None,
    },
)

# Define the schema
schema = EntityBuilder(
    entity_type=EntityTypeCv.SMALL_MOLECULE,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.CHEBI, value=f('chebi_id', extract='chebi')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('synonyms', delimiter=';')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.DESCRIPTION, value=f('description')),
    ),
)
```

### 5. Dataset

Combines download, parser, and schema:

```python
from pypath.inputs_v2.base import Dataset

dataset = Dataset(
    download=download,
    mapper=schema,        # EntityBuilder schema
    raw_parser=iter_tsv,  # Parser function
)
```

### 6. Resource

Container for multiple related datasets:

```python
from pypath.inputs_v2.base import Resource

resource = Resource(
    config,
    proteins=proteins_dataset,    # Datasets are accessed as attributes
    interactions=interactions_dataset,
)

# Access datasets
for entity in resource.proteins():
    print(entity)

# Get resource metadata
metadata = resource.metadata()
```

## Creating a New Input Module

Follow these steps to add a new data source:

### Step 1: Define ResourceConfig

Create a new file (e.g., `myresource.py`) and define the resource metadata:

```python
"""
Parse MyResource data and emit Entity records.

This module converts MyResource data into Entity records using the
declarative schema pattern.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig

config = ResourceConfig(
    id=ResourceCv.MY_RESOURCE,  # Add to cv_terms.py if needed
    name='My Resource',
    url='https://myresource.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='12345678',
    description='Description of the resource...',
)
```

### Step 2: Define the Schema

```python
f = FieldConfig(
    extract={
        # Named extraction patterns (regex or callables)
        'id_pattern': r'^ID:(\d+)$',
    },
    transform={
        # Named transformations
        'normalize': str.lower,
    },
)

schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,  # Or other entity type
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('uniprot_id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('organism_id')),
    ),
)
```

### Step 3: Create or Select a Parser

For standard formats, use existing parsers:

```python
from pypath.inputs_v2.parsers.base import iter_tsv
```

For custom formats, create a parser in `parsers/`:

```python
# parsers/myresource.py
def _raw(opener, **_kwargs):
    """Parse MyResource format."""
    if not opener or not opener.result:
        return
    # Parse and yield dictionaries
    for record in parse_data(opener.result):
        yield {
            'uniprot_id': record.id,
            'name': record.name,
            'organism_id': record.taxon,
        }
```

### Step 4: Create the Resource

```python
resource = Resource(
    config,
    entities=Dataset(
        download=Download(
            url='https://myresource.org/data.tsv',
            filename='myresource.tsv',
            subfolder='myresource',
        ),
        mapper=schema,
        raw_parser=iter_tsv,  # or custom parser
    ),
)
```

## FieldConfig Reference

`FieldConfig` provides a fluent interface for field extraction:

```python
f = FieldConfig(
    extract={
        'pattern_name': r'regex_pattern',  # Named regex
    },
    map={
        'enum_name': {'value': CV_TERM},   # Named mapping dict
    },
    transform={
        'func_name': lambda x: x.upper(),  # Named transform
    },
    delimiter=',',  # Default delimiter for multi-value fields
)

# Usage:
f('column_name')                           # Simple column access
f('column_name', delimiter=';')            # Split by semicolon
f('column_name', extract='pattern_name')   # Apply named regex
f('column_name', map='enum_name')          # Apply named mapping
f('column_name', transform='func_name')    # Apply named transform
```

## CV Reference

`CV` specifies a controlled vocabulary field:

```python
CV(term=IdentifierNamespaceCv.UNIPROT, value=f('accession'))
CV(term=f('type_column', map='type_mapping'))  # Dynamic term
CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed_ids', delimiter=';'))
```

## Advanced Patterns

### Dynamic URLs

```python
def get_url(**kwargs):
    taxon = kwargs.get('taxon', 9606)
    return f'https://example.com/data?taxon={taxon}'

download = Download(
    url=get_url,  # Callable
    filename=lambda **kw: f'data_{kw.get("taxon", 9606)}.tsv',
    subfolder='example',
)
```

### Complex Entities with Membership

For entities with members (complexes, families, interactions):

```python
from pypath.internals.tabular_builder import (
    Member,
    MembershipBuilder,
    MembersFromList,
)

schema = EntityBuilder(
    entity_type=EntityTypeCv.COMPLEX,
    identifiers=IdentifiersBuilder(...),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.PROTEIN,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.UNIPROT,
                   value=f('members', delimiter=',')),
            ),
        )
    ),
)
```

### Interactions with Participants

```python
schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(...),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('protein_a')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.SOURCE),  # Role annotation
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('protein_b')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.TARGET),
            ),
        ),
    ),
)
```

## Best Practices

1. **Use CV terms**: Always use `IdentifierNamespaceCv`, `EntityTypeCv`, and other CV enums from `cv_terms.py` rather than raw strings.

2. **Document parsers**: Include docstrings explaining the expected input format and output fields.

3. **Handle missing data**: Use `extract` patterns and `transform` functions to normalize and validate data.

4. **Reuse parsers**: Place common parsing logic in `parsers/` and share between modules.

5. **Test incrementally**: Use `max_records` parameter (if supported by parser) during development.

6. **Export resource**: Each module should export a `resource` variable at module level for automatic discovery.

## Usage Example

```python
from pypath.inputs_v2.uniprot import resource

# Get resource metadata
for meta in resource.metadata():
    print(meta)

# Iterate over proteins
for entity in resource.proteins():
    print(entity.identifiers)

# Access raw data before mapping
for record in resource.proteins.raw():
    print(record)  # dict from parser
```

## Discovery

The `omnipath_build` loaders automatically discover `Resource` objects exported at module level. Each resource's datasets are then processed to build the silver layer.
