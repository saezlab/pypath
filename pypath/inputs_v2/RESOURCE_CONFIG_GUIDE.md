# How to Write Resource Configurations (inputs_v2 Style)

## Overview

Resource configurations in `pypath/inputs_v2/` use a declarative schema pattern to convert tabular data into `SilverEntity` records.

## Basic Structure

```python
from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity as SilverEntity, Resource
from pypath.internals.cv_terms import EntityTypeCv, IdentifierNamespaceCv, LicenseCV, UpdateCategoryCV
from .tabular_builder import Annotations, Column, Entity, Identifiers
import csv

def get_resource() -> Resource:
    """Define resource metadata."""
    return Resource(
        id='resource_id',
        name='Resource Name',
        license=LicenseCV.CC_BY_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        publication='PMID:12345678',
        url='https://example.org/',
        description='Brief description of the resource.',
    )

def resource_data() -> Generator[SilverEntity]:
    """Download and parse resource data."""
    # 1. Download
    opener = download_and_open(
        'https://example.org/data.tsv',
        filename='data.tsv',
        subfolder='resource_name',
    )

    # 2. Define schema
    map = Entity(
        entity_type=EntityTypeCv.PROTEIN,
        identifiers=Identifiers(...),
        annotations=Annotations(...),
    )

    # 3. Parse and yield
    for row in csv.DictReader(opener.result, delimiter='\t'):
        yield map(row)
```

## Column Basics

A `Column` extracts and processes values from a single field:

```python
Column(
    'Column Name',              # Field selector (string for dict, int for tuple)
    delimiter=';',              # Split multi-value fields
    cv=IdentifierNamespaceCv.X, # Controlled vocabulary term
)
```

### Processing Options

```python
Column(
    'Field',
    processing={
        'extract_value': r'^([^(]+)',      # Regex to extract value (group 1)
        'extract_prefix': r'^(\w+):',      # Extract prefix
        'extract_unit': r'(\w+)$',         # Extract units
    },
    cv=IdentifierNamespaceCv.NAME,
)
```

## Identifiers

Map columns to entity identifiers:

```python
identifiers=Identifiers(
    Column('ID', cv=IdentifierNamespaceCv.UNIPROT),
    Column('Name', cv=IdentifierNamespaceCv.NAME),
    Column('Synonyms', delimiter=';', cv=IdentifierNamespaceCv.SYNONYM),

    # Extract primary name before parentheses
    Column('Full Name',
           processing={'extract_value': r'^([^(]+)'},
           cv=IdentifierNamespaceCv.NAME),

    # Extract synonyms from parentheses
    Column('Full Name',
           processing={'extract_value': r'\(([^)]+)\)'},
           cv=IdentifierNamespaceCv.SYNONYM),
)
```

**Note**: Multiple columns can reference the same field with different processing.

## Annotations

Map columns to metadata annotations:

```python
annotations=Annotations(
    # Simple annotations
    Column('Length', cv=AnnotationTypeCv.LENGTH),
    Column('Function', cv=AnnotationTypeCv.FUNCTION),

    # Multi-value annotations (split by delimiter)
    Column('GO IDs', delimiter=';', cv=AnnotationTypeCv.CV_TERM_ACCESSION),
    Column('Keywords', delimiter=';', cv=AnnotationTypeCv.CV_TERM_ACCESSION),

    # Special annotations
    Column('Organism ID', cv=AnnotationTypeCv.ORGANISM),
    Column('PubMed IDs', delimiter=';', cv=AnnotationTypeCv.PUBMED),
)
```

**Key Points**:
- GO terms and keywords use `AnnotationTypeCv.CV_TERM_ACCESSION`
- Organism and references are also annotations
- The value extracted becomes the annotation value

## Members (for Complexes/Families)

Define entity membership relationships:

```python
from .tabular_builder import Members, MembersFromList

members=Members(
    MembersFromList(
        entity_type=EntityTypeCv.PROTEIN,
        identifiers=Identifiers(
            Column('Member IDs', delimiter=',', cv=IdentifierNamespaceCv.UNIPROT),
        ),
    )
)
```

## Common Patterns

### Multi-value Fields
```python
# Semicolon-separated cross-references
Column('Ensembl IDs', delimiter=';', cv=IdentifierNamespaceCv.ENSEMBL)

# Space-separated gene synonyms
Column('Gene Synonyms', delimiter=' ', cv=IdentifierNamespaceCv.GENESYMBOL)
```

### Name Parsing
```python
# Primary name (before parentheses)
Column('Protein Name',
       processing={'extract_value': r'^([^(]+)'},
       cv=IdentifierNamespaceCv.NAME)

# Synonyms (inside parentheses) - extracts first match
Column('Protein Name',
       processing={'extract_value': r'\(([^)]+)\)'},
       cv=IdentifierNamespaceCv.SYNONYM)
```

### CV Term Annotations
```python
# GO terms, keywords, etc.
Column('GO IDs', delimiter=';', cv=AnnotationTypeCv.CV_TERM_ACCESSION)
```

## Examples

### Simple Resource (like SIGNOR complexes)
```python
def signor_complexes() -> Generator[SilverEntity]:
    opener = download_and_open(url, ...)

    map = Entity(
        entity_type=EntityTypeCv.PROTEIN_COMPLEX,
        identifiers=Identifiers(
            Column('SIGNOR ID', cv=IdentifierNamespaceCv.SIGNOR),
            Column('COMPLEX NAME', cv=IdentifierNamespaceCv.NAME),
        ),
        members=Members(
            MembersFromList(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=Identifiers(
                    Column('LIST OF ENTITIES', delimiter=',', cv=IdentifierNamespaceCv.UNIPROT),
                ),
            )
        ),
    )

    for row in csv.DictReader(opener.result, delimiter=';'):
        yield map(row)
```

### Complex Resource (like UniProt)
See [pypath/inputs_v2/uniprot.py](pypath/inputs_v2/uniprot.py) for a comprehensive example with:
- Multiple identifier types (primary, synonyms, cross-references)
- Regex-based name parsing
- Multi-value fields (GO terms, keywords, PubMed IDs)
- Extensive annotations (function, location, PTMs, etc.)

## Key Principles

1. **Declarative over imperative**: Define schema, don't iterate manually
2. **Reuse columns**: Same field can be processed multiple ways
3. **Use delimiters**: Let `Column` handle splitting multi-value fields
4. **Annotations for everything**: Organism, references, GO terms, keywords all use annotations
5. **CV terms**: Always use appropriate controlled vocabulary enums
