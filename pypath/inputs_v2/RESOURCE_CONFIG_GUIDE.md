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
from pypath.internals.silver_schema import Entity, Identifier, Annotation, Resource
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    AnnotationTypeCv,
    # Import other CV terms as needed
)
from ..internals.tabular_builder import (
    EntityBuilder,
    IdentifiersBuilder,
    AnnotationsBuilder,
    Column,
    CV,
    Map,
    Member,
    MembershipBuilder,
)

def get_resource() -> Resource:
    return Resource(
        id='resource_id',
        name='Resource Name',
        license=LicenseCV.CC_BY_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        publication='PMID:12345678',
        url='https://example.org/',
        description='Brief description.',
    )

def resource_data() -> Generator[Entity]:
    """Download and parse resource data as Entity records."""

    opener = download_and_open(
        'https://example.org/data.tsv',
        filename='data.tsv',
        subfolder='resource_name',
    )

    schema = EntityBuilder(
        entity_type=EntityTypeCv.PROTEIN,
        identifiers=IdentifiersBuilder(...),
        annotations=AnnotationsBuilder(...),
    )

    if opener and opener.result:
        for row in csv.DictReader(opener.result, delimiter='\t'):
            yield schema(row)

```
⸻

Core Components: Column, Map, and CV

The configuration DSL has three main building blocks:

1. **Column**: Pure data extraction from a field
   - Extracts values by column name (str), index (int), or callable
   - Optionally splits multi-value fields with `delimiter`
   - Does NOT include CV term or transformation logic

```python
Column('Field')                    # Extract single value
Column('Field', delimiter=';')     # Split multi-value field
Column(0)                          # Extract by index
```

2. **Map**: Value transformation and extraction
   - Applies regex patterns or callable functions to transform values
   - Maps extracted values to CV terms using dictionaries
   - Chains multiple extraction steps

```python
Map(
    col=Column('Field'),
    extract=[r'^([^(]+)', str.strip],  # Apply regex, then strip
    map={'value': 'mapped_value'},      # Optional mapping dict
    default='fallback',                 # Optional default value
)
```

3. **CV**: Semantic specification
   - Links a CV term (what it means) to a value source (where it comes from)
   - Separates term, value, and optional unit
   - Used by both IdentifiersBuilder and AnnotationsBuilder

```python
CV(
    term=IdentifierNamespaceCv.UNIPROT,  # What type of identifier
    value=Column('ID'),                  # Where to get the value
)

CV(
    term=AnnotationTypeCv.MOLECULAR_WEIGHT,
    value=Column('MW'),
    unit=Column('MW_Unit'),              # Optional unit specification
)
```

⸻

Identifiers

Define primary names, synonyms, and cross-references using CV objects:

```python
identifiers = IdentifiersBuilder(
    # Simple identifier
    CV(term=IdentifierNamespaceCv.UNIPROT, value=Column('Entry')),

    # Simple name
    CV(term=IdentifierNamespaceCv.NAME, value=Column('Name')),

    # Multi-value synonyms
    CV(
        term=IdentifierNamespaceCv.SYNONYM,
        value=Column('Synonyms', delimiter=';')
    ),

    # Extract name before parentheses
    CV(
        term=IdentifierNamespaceCv.NAME,
        value=Map(
            col=Column('Full Name'),
            extract=[r'^([^(]+)', str.strip]
        )
    ),

    # Extract parenthesized synonyms
    CV(
        term=IdentifierNamespaceCv.SYNONYM,
        value=Map(
            col=Column('Full Name'),
            extract=[r'\(([^)]+)\)']
        )
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
    CV(term=AnnotationTypeCv.LENGTH, value=Column('Length')),

    # Text annotation
    CV(term=AnnotationTypeCv.FUNCTION, value=Column('Function')),

    # Multi-value annotations
    CV(
        term=AnnotationTypeCv.CV_TERM_ACCESSION,
        value=Column('GO IDs', delimiter=';')
    ),
    CV(
        term=AnnotationTypeCv.PUBMED,
        value=Column('PubMed IDs', delimiter=';')
    ),

    # With unit specification
    CV(
        term=AnnotationTypeCv.MOLECULAR_WEIGHT,
        value=Column('MW'),
        unit=Column('MW_Unit'),  # e.g., "kDa"
    ),

    # Organism taxonomy ID
    CV(term=AnnotationTypeCv.ORGANISM, value=Column('Organism ID')),
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
                CV(term=IdentifierNamespaceCv.UNIPROT, value=Column('Interactor_A')),
            ),
            annotations=AnnotationsBuilder(
                CV(term=AnnotationTypeCv.ORGANISM, value=Column('TaxID_A')),
            ),
        )
    ),
    Member(
        entity=EntityBuilder(
            entity_type=EntityTypeCv.PROTEIN,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.UNIPROT, value=Column('Interactor_B')),
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
                value=Column('Member IDs', delimiter=',')
            ),
        ),
        annotations=AnnotationsBuilder(
            # Annotations applied to each member
            CV(term=AnnotationTypeCv.ORGANISM, value=Column('Organism')),
        ),
    )
)
```


⸻

Boolean Fields

For boolean fields (e.g., "yes"/"no" or "true"/"false"), use a mapping where the true value maps to the CV term itself. When the field value matches the key, the CV term is added to annotations; otherwise, it's omitted.

```python
annotations = AnnotationsBuilder(
    CV(term={'yes': MoleculeAnnotationsCv.APPROVED}, value=Column('Approved')),
    CV(term={'yes': MoleculeAnnotationsCv.WITHDRAWN}, value=Column('Withdrawn')),
    CV(term={'yes': MoleculeAnnotationsCv.LABELLED}, value=Column('Labelled')),
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
    value=Column('Gene Synonyms', delimiter=' ')
)
```

**Prefix-based ID mapping**:

Use Map with dynamic CV term extraction for heterogeneous identifier fields like `uniprotkb:P12345|chebi:15377`:

```python
# Define mapping from prefix to CV term
identifier_cv_mapping = {
    'uniprotkb': IdentifierNamespaceCv.UNIPROT,
    'chebi': IdentifierNamespaceCv.CHEBI,
    'ensembl': IdentifierNamespaceCv.ENSEMBL,
    'refseq': IdentifierNamespaceCv.REFSEQ,
}

# Use Map for both term and value extraction
CV(
    term=Map(
        col=Column('ID Column', delimiter='|'),
        extract=[r'^([^:]+):', str.lower],  # Extract prefix, lowercase
        map=identifier_cv_mapping,           # Map to CV term
    ),
    value=Map(
        col=Column('ID Column', delimiter='|'),
        extract=[r'^[^:]+:([^|"]+)'],        # Extract value after colon
    ),
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

schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,

    identifiers=IdentifiersBuilder(
        # Dynamic identifier type based on prefix
        CV(
            term=Map(
                col=Column('Interaction identifier(s)', delimiter='|'),
                extract=[r'^([^:]+):', str.lower],
                map=identifier_cv_mapping,
            ),
            value=Map(
                col=Column('Interaction identifier(s)', delimiter='|'),
                extract=[r'^[^:]+:(.+)'],
            ),
        ),
    ),

    annotations=AnnotationsBuilder(
        # Extract MI term accessions
        CV(
            term=AnnotationTypeCv.CV_TERM_ACCESSION,
            value=Map(
                col=Column('Interaction type(s)', delimiter='|'),
                extract=[r'(MI:\d+)'],
            ),
        ),
    ),

    membership=MembershipBuilder(
        # First interactor
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=Map(
                            col=Column('#ID(s) interactor A', delimiter='|'),
                            extract=[r'^([^:]+):', str.lower],
                            map=identifier_cv_mapping,
                        ),
                        value=Map(
                            col=Column('#ID(s) interactor A', delimiter='|'),
                            extract=[r'^[^:]+:(.+)'],
                        ),
                    ),
                ),
                annotations=AnnotationsBuilder(
                    CV(
                        term=AnnotationTypeCv.ORGANISM,
                        value=Map(
                            col=Column('Taxid interactor A'),
                            extract=[r'taxid:([-\d]+)'],
                        ),
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
                        term=Map(
                            col=Column('ID(s) interactor B', delimiter='|'),
                            extract=[r'^([^:]+):', str.lower],
                            map=identifier_cv_mapping,
                        ),
                        value=Map(
                            col=Column('ID(s) interactor B', delimiter='|'),
                            extract=[r'^[^:]+:(.+)'],
                        ),
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

schema = EntityBuilder(
    entity_type=Map(
        col=Column('Type'),
        map=type_mapping,
    ),
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

identifiers = IdentifiersBuilder(
    CV(
        term=IdentifierNamespaceCv.NAME,
        value=parse_ligand_name  # Function reference
    ),
)
```

**Chained Extractions with Map**

Apply multiple transformation steps in sequence:

```python
CV(
    term=AnnotationTypeCv.AFFINITY,
    value=Map(
        col=Column('Affinity'),
        extract=[
            r'([0-9.]+)',        # Extract number
            float,               # Convert to float
        ],
    ),
    unit=Map(
        col=Column('Affinity'),
        extract=[r'([a-zA-Z]+)$'],  # Extract unit at end
    ),
)
```

⸻

Core Principles
	1.	**Declarative schema**: Define mappings, do not manually build entities
	2.	**Separation of concerns**: Column extracts, Map transforms, CV assigns meaning
	3.	**Reusability**: Same column can be used with different transformations
	4.	**CV term priority**: Always prefer existing CV terms (PSI-MI first, then OmniPath)
	5.	**Multi-value fields**: Use `delimiter` in Column for split operations
	6.	**Dynamic extraction**: Use Map with dictionaries for prefix-based or conditional CV terms
	7.	**Search first**: Always search before adding new CV terms
	8.	**Type annotations**: Use Generator[Entity] for return types

⸻
