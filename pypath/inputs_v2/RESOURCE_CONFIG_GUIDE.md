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

### Prefix-Based Identifier Mapping

When dealing with data that contains multiple identifier types with prefixes (e.g., `uniprotkb:P12345`, `chebi:CHEBI:12345`), use a mapping dictionary instead of hardcoding a single identifier type:

```python
# Define mapping from prefixes to CV terms
identifier_cv_mapping = {
    'uniprotkb': IdentifierNamespaceCv.UNIPROT,
    'chebi': IdentifierNamespaceCv.CHEBI,
    'ensembl': IdentifierNamespaceCv.ENSEMBL,
    'refseq': IdentifierNamespaceCv.REFSEQ,
}

# Processing to extract prefix and value
general_id_processing = {
    'extract_prefix': r'^([^:]+):',
    'extract_value': r'^[^:]+:([^|"]+)'
}

# Use in schema
identifiers=Identifiers(
    Column('ID Column', delimiter='|',
           processing=general_id_processing,
           cv=identifier_cv_mapping),
)
```

**Best Practice**: Programmatically analyze the data to find all unique prefixes:

```python
# Analysis script to find all prefixes
import csv
from collections import defaultdict

prefixes = set()
with open('data_file.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        ids = row.get('ID Column', '')
        for id_str in ids.split('|'):
            if ':' in id_str:
                prefix = id_str.split(':', 1)[0].lower()
                prefixes.add(prefix)

print("All prefixes found:", sorted(prefixes))
```

Then create CV terms for any missing prefixes and build the complete mapping.

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

### Interaction Resource with Prefix-Based IDs (like SIGNOR or IntAct)

```python
def signor_interactions() -> Generator[SilverEntity]:
    opener = download_and_open(url, ...)

    # Mapping of identifier prefixes to CV terms
    identifier_cv_mapping = {
        'uniprotkb': IdentifierNamespaceCv.UNIPROT,
        'chebi': IdentifierNamespaceCv.CHEBI,
        'complexportal': IdentifierNamespaceCv.COMPLEXPORTAL,
        'pubchem': IdentifierNamespaceCv.PUBCHEM,
        'signor': IdentifierNamespaceCv.SIGNOR,
    }

    # Processing patterns
    general_id_processing = {'extract_prefix': r'^([^:]+):', 'extract_value': r'^[^:]+:([^|"]+)'}
    tax_processing = {'extract_value': r'taxid:([-\d]+)'}
    mi_term_processing = {'extract_term': r'(MI:\d+)'}

    schema = Entity(
        entity_type=EntityTypeCv.INTERACTION,
        identifiers=Identifiers(
            Column('Interaction identifier(s)', delimiter='|',
                   processing=general_id_processing,
                   cv=identifier_cv_mapping),
        ),
        annotations=Annotations(
            Column('Interaction type(s)', delimiter='|', processing=mi_term_processing),
        ),
        members=Members(
            Member(
                entity=Entity(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=Identifiers(
                        Column('#ID(s) interactor A', delimiter='|',
                               processing=general_id_processing,
                               cv=identifier_cv_mapping),
                        Column('Alt. ID(s) interactor A', delimiter='|',
                               processing=general_id_processing,
                               cv=identifier_cv_mapping),
                    ),
                    annotations=Annotations(
                        Column('Taxid interactor A', delimiter='|',
                               processing=tax_processing,
                               cv=IdentifierNamespaceCv.NCBI_TAX_ID),
                    ),
                ),
            ),
            Member(
                entity=Entity(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=Identifiers(
                        Column('ID(s) interactor B', delimiter='|',
                               processing=general_id_processing,
                               cv=identifier_cv_mapping),
                        Column('Alt. ID(s) interactor B', delimiter='|',
                               processing=general_id_processing,
                               cv=identifier_cv_mapping),
                    ),
                    annotations=Annotations(
                        Column('Taxid interactor B', delimiter='|',
                               processing=tax_processing,
                               cv=IdentifierNamespaceCv.NCBI_TAX_ID),
                    ),
                ),
            ),
        ),
    )

    for row in csv.DictReader(opener.result, delimiter='\t'):
        yield schema(row)
```

**Key Points**:
- Uses `identifier_cv_mapping` dict to handle multiple ID types
- Same mapping works for both ID and Alt ID columns
- Different processing patterns for different field types (IDs, taxonomy, MI terms)

### Complex Resource (like UniProt)
See [pypath/inputs_v2/uniprot.py](pypath/inputs_v2/uniprot.py) for a comprehensive example with:
- Multiple identifier types (primary, synonyms, cross-references)
- Regex-based name parsing
- Multi-value fields (GO terms, keywords, PubMed IDs)
- Extensive annotations (function, location, PTMs, etc.)

## Working with Controlled Vocabulary (CV) Terms

### Finding Existing CV Terms

Before adding new CV terms, always search for existing ones. CV terms are organized in `pypath/internals/cv_terms/`:

```
cv_terms/
├── __init__.py           # Main exports
├── core.py              # Base CvEnum class
├── entity_types.py      # Entity type terms (proteins, complexes, etc.)
├── identifiers.py       # Identifier namespaces (databases, IDs)
├── annotations.py       # Annotation terms (roles, methods, parameters)
└── resource_metadata.py # Licenses, update frequencies
```

**Search Strategy**:

1. **Grep for existing terms** in the codebase:
   ```bash
   # Search for a specific term (e.g., "ki" or "patent")
   grep -ri "patent\|ki value" pypath/internals/cv_terms/

   # List all available CV classes
   grep "^class.*Cv(CvEnum)" pypath/internals/cv_terms/*.py
   ```

2. **Check PSI-MI ontology** for standard terms (most interaction-related terms):
   ```bash
   # Look for binding affinity terms
   grep -i "ki value\|kd\|ic50\|ec50" pypath-data/obo/psi-mi.obo

   # Search by MI accession
   grep -A 5 "^id: MI:0640" pypath-data/obo/psi-mi.obo
   ```

3. **Examine similar resources**: Check existing resource configurations for similar data types:
   ```bash
   # Find resources using similar CV terms
   grep -l "InteractionParameterCv\|DetectionMethodCv" pypath/inputs_v2/*.py
   ```

### When to Add New CV Terms

**Use existing PSI-MI terms when available** (preferred):
- Check `psi-mi.obo` first for interaction-related terms
- PSI-MI terms are standardized and widely recognized
- Most interaction parameters, detection methods, and roles are already in PSI-MI

**Add OmniPath-specific terms (OM:XXXX) when**:
1. No suitable PSI-MI term exists
2. The data is OmniPath-specific (e.g., custom annotations)
3. You need a different unit or representation (e.g., Celsius instead of Kelvin)

**Accession Ranges**:
- `MI:XXXX` - PSI-MI standard terms (don't create new ones)
- `OM:0001-0009` - Specialized identifiers
- `OM:0010-0099` - Entity types
- `OM:0100-0209` - Database identifiers and names
- `OM:0300-0399` - Membership roles
- `OM:0400-0499` - Curation and update metadata
- `OM:0500-0599` - Licenses
- `OM:0600-0699` - Molecule annotations
- `OM:0700-0799` - Interaction parameters (custom)

### Adding New CV Terms

**Example: Adding binding affinity parameters**

1. **Search PSI-MI ontology**:
   ```bash
   grep -B 2 -A 5 "dissociation constant\|ki value" pypath-data/obo/psi-mi.obo
   ```

2. **Create CV class** in appropriate module (e.g., `annotations.py`):
   ```python
   class InteractionParameterCv(CvEnum):
       """Interaction parameter terms from PSI-MI.

       Describes kinetic and thermodynamic parameters for enzymatic or binding
       studies, including affinity measurements, rate constants, and experimental
       conditions.
       """

       parent_cv_term = "MI:0640"  # parameter type

       # PSI-MI standard terms
       KI = "MI:0643"  # Equilibrium constant for dissociation of inhibitor
       KD = "MI:0646"  # Equilibrium dissociation constant
       IC50 = "MI:0641"  # 50% maximum inhibitory response

       # OmniPath-specific term (when no PSI-MI equivalent exists)
       TEMPERATURE_CELSIUS = ("OM:0701", "Temperature in Celsius (C)")
   ```

3. **Export the CV class** in `cv_terms/__init__.py`:
   ```python
   from .annotations import (
       ...,
       InteractionParameterCv
   )
   ```

4. **Use in resource configuration**:
   ```python
   from pypath.internals.cv_terms import InteractionParameterCv

   annotations=Annotations(
       Column('Ki (nM)', cv=InteractionParameterCv.KI),
       Column('Temp (C)', cv=InteractionParameterCv.TEMPERATURE_CELSIUS),
   )
   ```

### CV Term Guidelines

1. **Prefer PSI-MI terms**: They're standardized and interoperable
2. **Document parent terms**: Always specify `parent_cv_term` for hierarchy
3. **Include definitions**: For OmniPath terms, provide clear descriptions
4. **Use tuple format for OmniPath terms**: `("OM:XXXX", "Description", "URL")`
5. **Check accession ranges**: Use correct OM range for the term type
6. **Add comments**: Include units and clarifications in comments

**Example with best practices**:
```python
class InteractionParameterCv(CvEnum):
    """Interaction parameter terms from PSI-MI.

    Describes kinetic and thermodynamic parameters for enzymatic or binding
    studies, including affinity measurements, rate constants, and experimental
    conditions (e.g., pH, temperature).
    """

    parent_cv_term = "MI:0640"  # parameter type - Parameter for enzymatic or binding kinetic studies

    # Binding affinity measurements (equilibrium constants)
    KI = "MI:0643"  # Equilibrium constant for dissociation of an inhibitor. Unit Molar.
    KD = "MI:0646"  # The equilibrium dissociation constant. Unit Molar.

    # Rate constants
    KON = "MI:0834"  # Association rate constant. Unit M-1 s-1
    KOFF = "MI:0835"  # Dissociation rate constant. Unit s-1

    # Experimental conditions
    PH = "MI:0837"  # pH at which interaction was determined
    TEMPERATURE = "MI:0836"  # Temperature at which interaction was determined. Unit KELVIN (K)
    TEMPERATURE_CELSIUS = ("OM:0701", "Temperature at which interaction was determined. Unit CELSIUS (C)")
```

## Key Principles

1. **Declarative over imperative**: Define schema, don't iterate manually
2. **Reuse columns**: Same field can be processed multiple ways
3. **Use delimiters**: Let `Column` handle splitting multi-value fields
4. **Annotations for everything**: Organism, references, GO terms, keywords all use annotations
5. **CV terms**: Always use appropriate controlled vocabulary enums
6. **Search before creating**: Check PSI-MI ontology and existing CV terms before adding new ones
