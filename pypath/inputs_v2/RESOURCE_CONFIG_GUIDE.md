Writing inputs_v2 Resource Configurations — Concise Guide

Overview

Resource configurations in pypath/inputs_v2/ convert tabular data into SilverEntity objects using a declarative schema. You define:
	1.	Resource metadata
	2.	A schema (identifiers, annotations, members)
	3.	Iterate rows → yield entities

⸻

Minimal Structure
```python
from __future__ import annotations

from collections.abc import Generator
import csv

from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    AnnotationTypeCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceAnnotationCv,
    ResourceCv,
    # Import other CV terms as needed
)
from ..internals.tabular_builder import (
    EntityBuilder,
    IdentifiersBuilder,
    AnnotationsBuilder,
    FieldConfig,
    CV,
    Member,
    MembershipBuilder,
)


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing resource metadata.
    """
    yield Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.RESOURCE_ID),
            Identifier(type=IdentifierNamespaceCv.NAME, value='Resource Name'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='12345678'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://example.org/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'Brief description of the resource, what it contains, '
                'and what kind of data it provides.'
            )),
        ],
    )


def resource_data() -> Generator[Entity]:
    """
    Download and parse resource data as Entity records.

    Yields:
        Entity records with the appropriate type (e.g., PROTEIN, INTERACTION)
    """
    opener = download_and_open(
        'https://example.org/data.tsv',
        filename='data.tsv',
        subfolder='resource_name',
    )

    f = FieldConfig()
    schema = EntityBuilder(
        entity_type=EntityTypeCv.PROTEIN,
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.UNIPROT, value=f('Entry')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=AnnotationTypeCv.LENGTH, value=f('Length')),
        ),
    )

    if opener and opener.result:
        for row in csv.DictReader(opener.result, delimiter='\t'):
            yield schema(row)

```
⸻

Core Components: FieldConfig and CV

The configuration DSL has two main building blocks:

1. **FieldConfig**: Reusable column factory
   - Centralizes extraction rules, mappings, and default delimiters
   - Produces `Column` instances via the shorthand `f(...)`
   - Keeps transformation logic close to schemas without repeating patterns

```python
f = FieldConfig(
    extract={
        'prefix_lower': [r'^([^:]+):', str.lower],
        'value': r'^[^:]+:([^|"]+)',
    },
    map={
        'identifier_cv': identifier_cv_mapping,
    },
    delimiter='|',
)

f('Field')                          # Extract single value
f('Field', delimiter=';')           # Split multi-value field
f('ID Column', extract='value')     # Apply named extract steps
f('ID Column', extract='prefix_lower', map='identifier_cv')  # Convert to CV terms
```

2. **CV**: Semantic specification
   - Links a CV term (what it means) to a value source (where it comes from)
   - Separates term, value, and optional unit
   - Used by both IdentifiersBuilder and AnnotationsBuilder

```python
CV(
    term=IdentifierNamespaceCv.UNIPROT,  # What type of identifier
    value=f('ID'),                       # Where to get the value
)

CV(
    term=AnnotationTypeCv.MOLECULAR_WEIGHT,
    value=f('MW'),
    unit=f('MW_Unit'),                   # Optional unit specification
)
```

⸻

Identifiers

Define primary names, synonyms, and cross-references using CV objects:

```python
f = FieldConfig(
    extract={
        'name': [r'^([^(]+)', str.strip],
        'synonym': r'\(([^)]+)\)',
    },
)

identifiers = IdentifiersBuilder(
    # Simple identifier
    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('Entry')),

    # Simple name
    CV(term=IdentifierNamespaceCv.NAME, value=f('Name')),

    # Multi-value synonyms
    CV(
        term=IdentifierNamespaceCv.SYNONYM,
        value=f('Synonyms', delimiter=';')
    ),

    # Extract name before parentheses
    CV(
        term=IdentifierNamespaceCv.NAME,
        value=f('Full Name', extract='name')
    ),

    # Extract parenthesized synonyms
    CV(
        term=IdentifierNamespaceCv.SYNONYM,
        value=f('Full Name', extract='synonym')
    ),
)
```

**Key principle**: Multiple CV definitions can reference the same field with different extraction logic.

⸻

Annotations

Annotations capture extra metadata, including references:

```python
annotations = AnnotationsBuilder(
    # Simple value annotation
    CV(term=AnnotationTypeCv.LENGTH, value=f('Length')),

    # Text annotation
    CV(term=AnnotationTypeCv.FUNCTION, value=f('Function')),

    # Multi-value annotations
    CV(
        term=AnnotationTypeCv.CV_TERM_ACCESSION,
        value=f('GO IDs', delimiter=';')
    ),
    CV(
        term=AnnotationTypeCv.PUBMED,
        value=f('PubMed IDs', delimiter=';')
    ),

    # With unit specification
    CV(
        term=AnnotationTypeCv.MOLECULAR_WEIGHT,
        value=f('MW'),
        unit=f('MW_Unit'),  # e.g., "kDa"
    ),

    # Organism taxonomy ID
    CV(term=AnnotationTypeCv.ORGANISM, value=f('Organism ID')),
)
```


⸻

Members (complexes/families)

For group relationships, use MembershipBuilder with either individual Members or MembersFromList:

**Single Member Pattern**:
```python
membership = MembershipBuilder(
    Member(
        entity=EntityBuilder(
            entity_type=EntityTypeCv.PROTEIN,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.UNIPROT, value=f('Interactor_A')),
            ),
            annotations=AnnotationsBuilder(
                CV(term=AnnotationTypeCv.ORGANISM, value=f('TaxID_A')),
            ),
        )
    ),
    Member(
        entity=EntityBuilder(
            entity_type=EntityTypeCv.PROTEIN,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.UNIPROT, value=f('Interactor_B')),
            ),
        )
    ),
)
```

**List-Based Member Pattern** (for multi-value fields):
```python
membership = MembershipBuilder(
    MembersFromList(
        entity_type=EntityTypeCv.PROTEIN,
        identifiers=IdentifiersBuilder(
            CV(
                term=IdentifierNamespaceCv.UNIPROT,
                value=f('Member IDs', delimiter=',')
            ),
        ),
        annotations=AnnotationsBuilder(
            # Annotations applied to each member
            CV(term=AnnotationTypeCv.ORGANISM, value=f('Organism')),
        ),
    )
)
```


⸻

Boolean Fields

For boolean fields (e.g., "yes"/"no" or "true"/"false"), use a mapping where the true value maps to the CV term itself. When the field value matches the key, the CV term is added to annotations; otherwise, it's omitted.

```python
annotations = AnnotationsBuilder(
    CV(term={'yes': MoleculeAnnotationsCv.APPROVED}, value=f('Approved')),
    CV(term={'yes': MoleculeAnnotationsCv.WITHDRAWN}, value=f('Withdrawn')),
    CV(term={'yes': MoleculeAnnotationsCv.LABELLED}, value=f('Labelled')),
)
```

**Result**:
• If Approved = "yes" → annotation includes APPROVED term
• If Approved = "" or "no" → no APPROVED term added

This avoids storing redundant values like "APPROVED=yes" and instead produces cleaner output where presence of the term indicates true.


⸻

Multi-value and Prefix-Based Fields

**Simple multi-value fields**:
```python
CV(
    term=IdentifierNamespaceCv.GENESYMBOL,
    value=f('Gene Synonyms', delimiter=' ')
)
```

**Prefix-based ID mapping**:

Use FieldConfig with dynamic CV term extraction for heterogeneous identifier fields like `uniprotkb:P12345|chebi:15377`:

```python
# Define mapping from prefix to CV term
identifier_cv_mapping = {
    'uniprotkb': IdentifierNamespaceCv.UNIPROT,
    'chebi': IdentifierNamespaceCv.CHEBI,
    'ensembl': IdentifierNamespaceCv.ENSEMBL,
    'refseq': IdentifierNamespaceCv.REFSEQ,
}

# Configure a reusable FieldConfig helper
f = FieldConfig(
    extract={
        'prefix_lower': [r'^([^:]+):', str.lower],
        'value': r'^[^:]+:([^|"]+)',
    },
    map={
        'identifier_cv': identifier_cv_mapping,
    },
    delimiter='|',
)

# Use f for both term and value extraction
CV(
    term=f('ID Column', extract='prefix_lower', map='identifier_cv'),
    value=f('ID Column', extract='value'),
)
```

**Finding all prefixes in a file**:

Before creating your mapping, discover what prefixes exist:

```python
prefixes = set()
for row in csv.DictReader(open('file.tsv'), delimiter='\t'):
    for id_str in row['ID Column'].split('|'):
        if ':' in id_str:
            prefixes.add(id_str.split(':', 1)[0].lower())
print(sorted(prefixes))
```


⸻

Example: Interaction Resource (MITAB Format)

```python
# Mapping for heterogeneous identifier prefixes
identifier_cv_mapping = {
    'uniprotkb': IdentifierNamespaceCv.UNIPROT,
    'chebi': IdentifierNamespaceCv.CHEBI,
    'intact': IdentifierNamespaceCv.INTACT,
}

f = FieldConfig(
    extract={
        'prefix_lower': [r'^([^:]+):', str.lower],
        'value': r'^[^:]+:(.+)',
        'mi': r'(MI:\d+)',
        'tax': r'taxid:([-\d]+)',
    },
    map={
        'identifier_cv': identifier_cv_mapping,
    },
    delimiter='|',
)

schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,

    identifiers=IdentifiersBuilder(
        # Dynamic identifier type based on prefix
        CV(
            term=f('Interaction identifier(s)', extract='prefix_lower', map='identifier_cv'),
            value=f('Interaction identifier(s)', extract='value'),
        ),
    ),

    annotations=AnnotationsBuilder(
        # Extract MI term accessions
        CV(
            term=AnnotationTypeCv.CV_TERM_ACCESSION,
            value=f('Interaction type(s)', extract='mi'),
        ),
    ),

    membership=MembershipBuilder(
        # First interactor
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=f('#ID(s) interactor A', extract='prefix_lower', map='identifier_cv'),
                        value=f('#ID(s) interactor A', extract='value'),
                    ),
                ),
                annotations=AnnotationsBuilder(
                    CV(
                        term=AnnotationTypeCv.ORGANISM,
                        value=f('Taxid interactor A', extract='tax'),
                    ),
                ),
            )
        ),
        # Second interactor (similar pattern)
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=f('ID(s) interactor B', extract='prefix_lower', map='identifier_cv'),
                        value=f('ID(s) interactor B', extract='value'),
                    ),
                ),
            )
        ),
    ),
)
```


⸻

Working with CV Terms

Where CV terms live

pypath/internals/cv_terms/ contains:
	•	entity_types.py
	•	identifiers.py
	•	annotations.py
	•	resource_metadata.py

Search first before adding new terms:

grep -ri "KI\|patent\|affinity" pypath/internals/cv_terms/

Also search PSI-MI:

grep -i "Kd\|Ki\|IC50" pypath-data/obo/psi-mi.obo


⸻

Adding CV Terms (only if no PSI-MI term exists)

Prefer PSI-MI (MI:XXXX). Add OmniPath (OM:XXXX) only when necessary.

Example:

class InteractionParameterCv(CvEnum):
    parent_cv_term = "MI:0640"  # parameter type

    KI = "MI:0643"
    KD = "MI:0646"
    KON = "MI:0834"
    KOFF = "MI:0835"
    PH = "MI:0837"
    TEMPERATURE = "MI:0836"  # Kelvin

    # Custom OmniPath extension
    TEMPERATURE_CELSIUS = ("OM:0701", "Temperature in Celsius (C)")


⸻

Advanced Patterns

**Dynamic Entity Types**

Entity types can be determined from the data:

```python
type_mapping = {
    'protein': EntityTypeCv.PROTEIN,
    'complex': EntityTypeCv.COMPLEX,
    'smallmolecule': EntityTypeCv.SMALL_MOLECULE,
}

f = FieldConfig()

schema = EntityBuilder(
    entity_type=f('Type', map=type_mapping),
    identifiers=IdentifiersBuilder(...),
)
```

**Callable Functions for Complex Parsing**

For complex parsing logic that can't be expressed with regex:

```python
def parse_ligand_name(row):
    """Custom function to extract ligand name."""
    raw = row.get('Ligand Name', '')
    # Complex parsing logic here
    return cleaned_name

f = FieldConfig()

identifiers = IdentifiersBuilder(
    CV(
        term=IdentifierNamespaceCv.NAME,
        value=f(parse_ligand_name)  # Function reference
    ),
)
```

**Chained Extractions with FieldConfig**

Apply multiple transformation steps in sequence:

```python
f = FieldConfig()

CV(
    term=AnnotationTypeCv.AFFINITY,
    value=f(
        'Affinity',
        extract=[
            r'([0-9.]+)',        # Extract number
            float,               # Convert to float
        ],
    ),
    unit=f('Affinity', extract=[r'([a-zA-Z]+)$']),  # Extract unit at end
)
```

⸻

Core Principles
	1.	**Declarative schema**: Define mappings, do not manually build entities
	2.	**Separation of concerns**: FieldConfig extracts/transforms, CV assigns meaning
	3.	**Reusability**: Same column can be used with different transformations
	4.	**CV term priority**: Always prefer existing CV terms (PSI-MI first, then OmniPath)
	5.	**Multi-value fields**: Use `delimiter` in FieldConfig for split operations
	6.	**Dynamic extraction**: Use FieldConfig with dictionaries for prefix-based or conditional CV terms
	7.	**Search first**: Always search before adding new CV terms
	8.	**Type annotations**: Use Generator[Entity] for return types

⸻
