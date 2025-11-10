Writing inputs_v2 Resource Configurations — Concise Guide

Overview

Resource configurations in pypath/inputs_v2/ convert tabular data into SilverEntity objects using a declarative schema. You define:
	1.	Resource metadata
	2.	A schema (identifiers, annotations, members)
	3.	Iterate rows → yield entities

⸻

Minimal Structure

from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity as SilverEntity, Resource
from pypath.internals.cv_terms import *
from .tabular_builder import Entity, Identifiers, Annotations, Column
import csv

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

def resource_data():
    opener = download_and_open(
        'https://example.org/data.tsv',
        filename='data.tsv',
        subfolder='resource_name',
    )

    schema = Entity(
        entity_type=EntityTypeCv.PROTEIN,
        identifiers=Identifiers(...),
        annotations=Annotations(...),
    )

    for row in csv.DictReader(opener.result, delimiter='\t'):
        yield schema(row)


⸻

Column Basics

Column extracts values from a field and optionally splits/filters them:

Column(
    'Field',                    # key or index
    delimiter=';',              # split multi-value fields
    cv=IdentifierNamespaceCv.X, # or AnnotationTypeCv.X
)

Regex-based processing

Column(
    'Field',
    processing={
        'extract_value': r'^([^(]+)',
        'extract_prefix': r'^(\w+):',
        'extract_unit': r'(\w+)$',
    },
    cv=IdentifierNamespaceCv.NAME,
)


⸻

Identifiers

Define primary names, synonyms, and cross-references:

identifiers = Identifiers(
    Column('ID', cv=IdentifierNamespaceCv.UNIPROT),
    Column('Name', cv=IdentifierNamespaceCv.NAME),
    Column('Synonyms', delimiter=';', cv=IdentifierNamespaceCv.SYNONYM),

    # Parse name before parentheses
    Column('Full Name',
           processing={'extract_value': r'^([^(]+)'},
           cv=IdentifierNamespaceCv.NAME),

    # Parenthesized synonyms
    Column('Full Name',
           processing={'extract_value': r'\(([^)]+)\)'},
           cv=IdentifierNamespaceCv.SYNONYM),
)

Multiple Column definitions for a single field allow different processing logic.

⸻

Annotations

Annotations capture extra metadata, including references:

annotations = Annotations(
    Column('Length', cv=AnnotationTypeCv.LENGTH),
    Column('Function', cv=AnnotationTypeCv.FUNCTION),
    Column('GO IDs', delimiter=';', cv=AnnotationTypeCv.CV_TERM_ACCESSION),
    Column('PubMed IDs', delimiter=';', cv=AnnotationTypeCv.PUBMED),
    Column('Organism ID', cv=AnnotationTypeCv.ORGANISM),
)


⸻

Members (complexes/families)

For group relationships:

from .tabular_builder import Members, MembersFromList

members = Members(
    MembersFromList(
        entity_type=EntityTypeCv.PROTEIN,
        identifiers=Identifiers(
            Column('Member IDs', delimiter=',', cv=IdentifierNamespaceCv.UNIPROT),
        ),
    )
)


⸻

Multi-value and Prefix-Based Fields

Simple multi-value fields

Column('Gene Synonyms', delimiter=' ', cv=IdentifierNamespaceCv.GENESYMBOL)

Prefix-based ID mapping

Use a dictionary to map prefixes → CV terms:

identifier_cv_mapping = {
    'uniprotkb': IdentifierNamespaceCv.UNIPROT,
    'chebi': IdentifierNamespaceCv.CHEBI,
    'ensembl': IdentifierNamespaceCv.ENSEMBL,
    'refseq': IdentifierNamespaceCv.REFSEQ,
}

processing = {
    'extract_prefix': r'^([^:]+):',
    'extract_value': r'^[^:]+:([^|"]+)',
}

Identifiers(
    Column('ID Column', delimiter='|',
           processing=processing,
           cv=identifier_cv_mapping),
)

Finding all prefixes in a file

prefixes = set()
for row in csv.DictReader(open('file.tsv'), delimiter='\t'):
    for id_str in row['ID Column'].split('|'):
        if ':' in id_str:
            prefixes.add(id_str.split(':', 1)[0].lower())
print(sorted(prefixes))


⸻

Example: Interaction Resource

schema = Entity(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=Identifiers(
        Column('Interaction identifier(s)', delimiter='|',
               processing=general_id_processing,
               cv=identifier_cv_mapping),
    ),
    annotations=Annotations(
        Column('Interaction type(s)', delimiter='|',
               processing={'extract_term': r'(MI:\d+)'}),
    ),
    members=Members(
        Member(
            entity=Entity(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=Identifiers(
                    Column('#ID(s) interactor A', delimiter='|',
                           processing=general_id_processing,
                           cv=identifier_cv_mapping),
                ),
                annotations=Annotations(
                    Column('Taxid interactor A', delimiter='|',
                           processing={'extract_value': r'taxid:([-\d]+)'},
                           cv=IdentifierNamespaceCv.NCBI_TAX_ID),
                ),
            )
        ),
        ...
    ),
)


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

Core Principles
	1.	Declarative schema: define mappings, do not manually build entities.
	2.	Columns are reusable with different processing.
	3.	Always prefer existing CV terms (PSI-MI first).
	4.	Use delimiters for multi-value fields.
	5.	Annotations for GO, keywords, organism IDs, references.
	6.	Prefix-based processing for heterogeneous IDs.
	7.	Search before adding new CV terms.

⸻
