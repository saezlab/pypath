# inputs_v2 - Declarative Data Input Modules

This module provides a declarative framework for downloading, parsing, and normalizing biological data from external resources into a unified `Entity` schema.

Read [Biolink modeling conventions](BIOLINK_MODELING.md) before changing an input module.

`inputs_v2` needs the `inputs-v2` extra (`pypath-omnipath[inputs-v2]`), which adds
`biolink-model`, `polars` and `omnipath-core`. `omnipath-core` is not on PyPI yet:
it lives in the [omnipath-build](https://github.com/saezlab/omnipath-build)
superproject (`packages/core`), which contains pypath as its `pypath` submodule, and
uv resolves it from there (`../packages/core`). The rest of pypath installs without it.

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
├── bindingdb.py          # BindingDB interactions & compounds
├── cellinker.py          # Cellinker ligand-receptor & metabolite-protein data
├── drugcentral.py        # DrugCentral drug-target interactions
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
from omnipath_core.resource_metadata import Content
from pypath.inputs_v2.base import ResourceConfig
from pypath.internals.cv_terms import ResourceCv, LicenseCV, UpdateCategoryCV

config = ResourceConfig(
    id=ResourceCv.UNIPROT,           # CV term identifying the resource
    name='UniProt',                   # Display name
    url='https://www.uniprot.org/',   # Resource URL
    license=LicenseCV.CC_BY_4_0,      # License term, with its academic/commercial use
    update_category=UpdateCategoryCV.REGULAR,  # How often data updates
    pubmed='33237286',                # Optional PubMed reference
    # Tags: what this module imports (see below)
    content=(Content.MOLECULAR_ANNOTATIONS, Content.IDENTIFIER_MAPPINGS),
    # One line on what this module imports
    summary='Protein identities, functional annotations and reference identifier translations.',
    description='UniProt is a comprehensive resource for protein...',  # The resource itself
)
```

The website's resources page is built from this object: `name` is the name it shows,
`summary` the line under it, the three tag fields its filters, and `license` its
license and permitted use.

- `name` is the resource's own short spelling (`HMDB`, `Mondo`). Put a longer name in
  `full` (`Human Metabolome Database`).
- `summary` says what this module imports, in one line. `description` describes the
  resource itself and can be longer.
- Academic and commercial use are not set per resource. Each `LicenseCV` term carries
  them (`LicenseCV.CC_BY_NC_4_0.commercial` is `LicenseUse.NOT_ALLOWED`); to classify a
  new license, add a term with both uses in `pypath/internals/cv_terms/resource_metadata.py`.

#### Tags

Three fields say what the module imports, one per group of filters on the website's
resources page. Each takes members of its own enum from `omnipath_core.resource_metadata`:

| Field | Enum | Says | Examples |
| --- | --- | --- | --- |
| `topics` | `Topic` | The area of biology | `SIGNALING`, `METABOLISM`, `PHARMACOLOGY` |
| `content` | `Content` | The kind of records | `PATHWAYS`, `LIGAND_RECEPTOR_INTERACTIONS`, `ONTOLOGIES` |
| `knowledge_origin` | `KnowledgeOrigin` | Where the knowledge comes from | `LITERATURE_CURATED`, `EXPERIMENTAL`, `TEXT_MINED` |

They are a controlled vocabulary: the enums have one member per entry of
`omnipath_core/vocab/resource_tags.yaml` (`packages/core` of the omnipath-build
superproject), grouped by the entry's `dimension`. A value outside its field's enum
fails when the module is imported. The YAML file has every tag with its dimension and
description. To list one group's IDs:

```bash
uv run python -c "from omnipath_core.resource_metadata import Content; print(*Content, sep='\n')"
```

`config.tags` is the three fields joined, in that order.

The page's fourth group, Entities, is not declared: it lists the entity types the
build found in the data (`Protein`, `ChemicalEntity`, `Pathway`, ...), so the mappers'
`entity_type` choices are what it shows.

Rules for assigning them:

- Tag what the module imports, not everything the provider offers.
- A knowledge-origin tag means the resource contains that kind of knowledge, not that
  every record has that evidence.
- Leave out a tag you are unsure of; a field may stay empty, but every resource needs
  at least one tag.

To add a tag, add it to `resource_tags.yaml` with a label, dimension and description;
it becomes a member of that dimension's enum.
After changing any `ResourceConfig`, regenerate the copy the API reads, from the
superproject root:

```bash
uv run python scripts/generate_resource_metadata.py
```

`make check-generated` fails while that copy is stale.

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
from biolink_model.datamodel.model import ChemicalEntity, slots
from omnipath_core.naming import Namespace

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
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.CHEBI, value=f('chebi_id', extract='chebi')),
        CV(term=Namespace.NAME, value=f('name')),
        CV(term=Namespace.SYNONYM, value=f('synonyms', delimiter=';')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.description, value=f('description')),
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

# Optional special kinds:
# Dataset(..., kind='id_translation')  # excluded from regular silver parquet output
```

A resource may also expose a narrow tabular dataset for a specific downstream use case.
For example, some small-molecule resources export an `id_translation` dataset with rows like:

```python
{
    'source': 'chebi',
    'key_type': Namespace.CHEBI,
    'key_value': '15377',
    'standard_inchi': 'InChI=...'
}
```

Use this when a source can provide a stable, source-owned translation view directly,
rather than forcing downstream code to reinterpret broader entity datasets.

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

## Newly Migrated Inputs

The following legacy inputs have been migrated into `inputs_v2` with known
download entry points and module-level `resource` definitions.

| Resource | Download entry point | `inputs_v2` datasets | Legacy module |
| --- | --- | --- | --- |
| Cellinker | https://www.cellknowledge.com.cn/cellinker/Download.html | Human and mouse protein ligand-receptor, metabolite receptor, enzyme, and transporter datasets | `pypath.inputs.cellinker` |
| DrugCentral | https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz | Drug-target interaction dataset | `pypath.inputs.drugcentral` |

### Step 1: Define ResourceConfig

Create a new file (e.g., `myresource.py`) and define the resource metadata:

```python
"""
Parse MyResource data and emit Entity records.

This module converts MyResource data into Entity records using the
declarative schema pattern.
"""

from __future__ import annotations

from biolink_model.datamodel.model import Protein, MacromolecularComplex, slots
from omnipath_core.naming import Namespace

from pypath.internals.cv_terms import (
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
from omnipath_core.resource_metadata import Content, Topic
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig

config = ResourceConfig(
    id=ResourceCv.MY_RESOURCE,  # Add to cv_terms.py if needed
    name='My Resource',
    url='https://myresource.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='12345678',
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    summary='One line on what this module imports.',
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
    entity_type=Protein,  # Or other entity type
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=f('uniprot_id')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=lambda row: f"NCBITaxon:{row['organism_id']}"),
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
CV(term=Namespace.UNIPROT, value=f('accession'))
CV(term=f('type_column', map='type_mapping'))  # Dynamic term
CV(term=slots.publications, value=lambda row: [f'PMID:{v}' for v in row['pubmed_ids'].split(';')])
```

## Advanced Patterns

### Efficient repeated mappings

Caching is automatic in the shared mapping executor. Source modules describe
what to extract and transform; they do not need cache decorators or duplicate
lists of cache-key columns.

```python
participant = EntityBuilder(
    entity_type=Protein,
    identifiers=IdentifiersBuilder(
        CV.from_pairs(f('cross_references', transform=parse_identifier_pairs)),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.name, value=f('label'))),
)
# parse_identifier_pairs returns ((Namespace.UNIPROT, 'P04637'), ...)
```

`CV.from_pairs(Column(...))` projects a cell's sequence of pairs into aligned
namespace/value columns, evaluating the source only once per row. The
`CV.from_pairs(row_callback)` form remains supported.

Builders evaluate every row; nothing is memoized across rows. Within a row,
fields that share a source extract it once. Cross-row caches were removed
because, on measured inputs, building their keys and copying results cost more
than the repeated work they saved. Callbacks should still be pure
transformations, and mapping specifications/constants must remain stable.
Declare reusable fields when defining a mapper rather than constructing them
inside per-row callbacks.

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
    entity_type=MacromolecularComplex,
    identifiers=IdentifiersBuilder(...),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=Protein,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.UNIPROT,
                   value=f('members', delimiter=',')),
            ),
        )
    ),
)
```

### Interactions (Relations)

```python
schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=f('protein_a')),
        ),
    ),
    predicate=slots.interacts_with,  # or callable / Column mapping
    object=EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=f('protein_b')),
        ),
    ),
    identifiers=IdentifiersBuilder(...),
    annotations=AnnotationsBuilder(...),
)
```

## Best Practices

1. **Use Biolink primitives**: Follow [the modeling conventions](BIOLINK_MODELING.md). Use native classes and slots, Biolink enum values, and `Namespace` for identifiers.

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

## Command Line Smoke Test

To print the first mapped Entity from the human MEBOCOST dataset:

```bash
poetry run python -c "from pypath.inputs_v2.mebocost import resource; print(next(resource.human()))"
```

## Discovery

The `omnipath_build` loaders automatically discover `Resource` objects exported at module level. Each resource's datasets are then processed to build the silver layer.

## Quantitative attributes

Keep assay measurements structured. The annotation term describes the measured
endpoint; its value can be a native Biolink `QuantityValue`. A `Measurement` wrapper
preserves the source column and original comparison operator alongside it:

```python
from biolink_model.datamodel.model import QuantityValue
from omnipath_core.measurements import Measurement

# BAO:0000190 is IC50. A source-specific parser has read "<=10" nM.
attribute = CV(
    term='BAO:0000190',
    value=lambda row: Measurement(
        quantity=QuantityValue(has_numeric_value=10, has_unit='nM'),
        source_field='IC50 (nM)',
        comparator='<=',
    ),
)
```

A real mapper parses each source value, keeps its units and handles malformed cells
explicitly. Raw IC50 is not `pIC50`: native log-affinity slots are only used for
source log values. Source-field provenance distinguishes reported minimum, maximum
and median without inventing a statistical interpretation. Quantities are retained
in the serving annotations' optional typed `quantity` struct and in each evidence
occurrence. Effects, mechanisms and other assertion qualifiers remain separate.

Use `description` for actual prose. Prefer a faithful native slot or a reviewed
external ontology attribute for structured fields; document fields that cannot yet
be exposed faithfully. Do not replace useful measurements with descriptions.

Nested participants may specify an explicit native `predicate` on `Member` or
`MembersFromList` (for example, reaction `has_input`/`has_output`/`enabled_by`). The
default remains `has_member`. A dynamic predicate returning `None` omits that member;
use `preserve_indices=True` for aligned columns so missing roles cannot shift the
predicates onto different participants.

Source names, synonyms and abbreviations are identifier values: use `Namespace.NAME`
and `Namespace.SYNONYM` in the identifier builder, after source accessions. They are
not emitted as name/synonym annotations; the serving label is derived from these
identifiers. This convention also applies to nested members and ontology terms.
