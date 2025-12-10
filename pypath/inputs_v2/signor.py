"""
Parse SIGNOR data and emit Entity records.

This module converts SIGNOR data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import EntityTypeCv, IdentifierNamespaceCv, LicenseCV, UpdateCategoryCV, ResourceAnnotationCv, ResourceCv
from ..internals.tabular_builder import (
    AnnotationsBuilder,
    Column,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    Map,
    Member,
    MembershipBuilder,
    MembersFromList,
)
import csv


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing SIGNOR metadata.
    """
    yield Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.SIGNOR),
            Identifier(type=IdentifierNamespaceCv.NAME, value='SIGNOR'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='31665520'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://signor.uniroma2.it/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'SIGNOR (SIGnaling Network Open Resource) is a comprehensive '
                'resource of causal relationships between biological entities '
                'with a focus on signaling pathways. It provides manually curated '
                'interactions with mechanistic details including protein-protein '
                'interactions, post-translational modifications, transcriptional '
                'regulation, and small molecule effects.'
            )),
        ],
    )


def signor_complexes() -> Generator[Entity]:
    """
    Download and parse SIGNOR complex data as Entity records.

    Yields:
        Entity records with type PROTEIN_COMPLEX, containing member proteins
    """
    # Download and open the file
    url = 'https://signor.uniroma2.it/download_complexes.php'
    opener = download_and_open(
        url,
        filename='signor_complexes.txt',
        subfolder='signor',
        query={'submit': 'Download complex data'},
        post=True,
    )

    # Define the schema mapping
    map = EntityBuilder(
        entity_type=EntityTypeCv.COMPLEX,
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.SIGNOR, value=Column('SIGNOR ID')),
            CV(term=IdentifierNamespaceCv.NAME, value=Column('COMPLEX NAME')),
        ),
        membership=MembershipBuilder(
            MembersFromList(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.UNIPROT,
                        value=Column('LIST OF ENTITIES', delimiter=','),
                    ),
                ),
            )
        ),
    )

    # Parse and yield entities
    for entity in csv.DictReader(opener.result, delimiter=';'):
        yield map(entity)


def signor_protein_families() -> Generator[Entity]:
    """
    Download and parse SIGNOR protein family data as Entity records.

    Yields:
        Entity records with type PROTEIN_FAMILY, containing member proteins
    """
    # Download and open the file
    url = 'https://signor.uniroma2.it/download_complexes.php'
    opener = download_and_open(
        url,
        filename='signor_protein_families.txt',
        subfolder='signor',
        query={'submit': 'Download protein family data'},
        post=True,
    )

    # Define the schema mapping
    map = EntityBuilder(
        entity_type=EntityTypeCv.PROTEIN_FAMILY,
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.SIGNOR, value=Column('SIGNOR ID')),
            CV(term=IdentifierNamespaceCv.NAME, value=Column('PROT. FAMILY NAME')),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
        ),
        membership=MembershipBuilder(
            MembersFromList(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.UNIPROT,
                        value=Column('LIST OF ENTITIES', delimiter=','),
                    ),
                ),
            )
        ),
    )

    # Parse and yield entities
    for entity in csv.DictReader(opener.result, delimiter=';'):
        yield map(entity)


def signor_phenotypes() -> Generator[Entity]:
    """
    Download and parse SIGNOR phenotype data as Entity records.

    Yields:
        Entity records with type PHENOTYPE
    """
    # Download and open the file
    url = 'https://signor.uniroma2.it/download_complexes.php'
    opener = download_and_open(
        url,
        filename='signor_phenotypes.txt',
        subfolder='signor',
        query={'submit': 'Download phenotype data'},
        post=True,
    )

    # Define the schema mapping
    map = EntityBuilder(
        entity_type=EntityTypeCv.PHENOTYPE,
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.SIGNOR, value=Column('SIGNOR ID')),
            CV(term=IdentifierNamespaceCv.NAME, value=Column('PHENOTYPE NAME')),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            CV(term=IdentifierNamespaceCv.SYNONYM, value=Column('PHENOTYPE DESCRIPTION')),
        ),
    )

    # Parse and yield entities
    for entity in csv.DictReader(opener.result, delimiter=';'):
        yield map(entity)


def signor_stimuli() -> Generator[Entity]:
    """
    Download and parse SIGNOR stimulus data as Entity records.

    Yields:
        Entity records with type STIMULUS
    """
    # Download and open the file
    url = 'https://signor.uniroma2.it/download_complexes.php'
    opener = download_and_open(
        url,
        filename='signor_stimuli.txt',
        subfolder='signor',
        query={'submit': 'Download stimulus data'},
        post=True,
    )

    # Define the schema mapping
    map = EntityBuilder(
        entity_type=EntityTypeCv.STIMULUS,
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.SIGNOR, value=Column('SIGNOR ID')),
            CV(term=IdentifierNamespaceCv.NAME, value=Column('STIMULUS NAME')),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            CV(term=IdentifierNamespaceCv.SYNONYM, value=Column('STIMULUS DESCRIPTION')),
        ),
    )

    # Parse and yield entities
    for entity in csv.DictReader(opener.result, delimiter=';'):
        yield map(entity)


def signor_interactions() -> Generator[Entity, None, None]:
    """
    Download SIGNOR causalTab interactions and yield `Entity` objects.
    """

    opener = download_and_open(
        url='https://signor.uniroma2.it/download_entity.php',
        filename='signor_all_causalTab.txt',
        subfolder='signor',
        query={
            'format': 'causalTab',
            'submit': 'Download',
        },
        post=True,
    )

    # Mapping of identifier prefixes (from SIGNOR data) to CV terms
    # All prefixes found in ID and Alt ID columns: chebi, complexportal, pubchem, signor, uniprotkb
    identifier_cv_mapping = {
        'chebi': IdentifierNamespaceCv.CHEBI,
        'complexportal': IdentifierNamespaceCv.COMPLEXPORTAL,
        'pubchem': IdentifierNamespaceCv.PUBCHEM,
        'signor': IdentifierNamespaceCv.SIGNOR,
        'uniprotkb': IdentifierNamespaceCv.UNIPROT,
    }

    # Processing helpers
    prefix_regex = r'^([^:]+):'
    general_value_regex = r'^[^:]+:([^|"]+)'
    interaction_value_regex = r'^[^:]+:(.*)'
    tax_regex = r'taxid:([-\d]+)'
    pubmed_regex = r'(?i)pubmed:(\d+)'
    mi_regex = r'(MI:\d+)'

    term_mapping = {
        'signor': IdentifierNamespaceCv.SIGNOR,
        'signor-interaction': IdentifierNamespaceCv.SIGNOR,
    }

    def interaction_identifier_cv() -> CV:
        column = Column('Interaction identifier(s)', delimiter='|')
        return CV(
            term=Map(
                col=column,
                extract=[prefix_regex, str.lower],
                map=term_mapping,
            ),
            value=Map(col=column, extract=[interaction_value_regex]),
        )

    def general_identifier_cv(column_name: str) -> CV:
        column = Column(column_name, delimiter='|')
        return CV(
            term=Map(
                col=column,
                extract=[prefix_regex, str.lower],
                map=identifier_cv_mapping,
            ),
            value=Map(col=column, extract=[general_value_regex]),
        )

    def mi_term_cv(column_name: str) -> CV:
        column = Column(column_name, delimiter='|')
        return CV(term=Map(col=column, extract=[mi_regex]))

    def mi_term_string(column_name: str) -> Map:
        """Extract MI term as a plain string for entity_type field."""
        column = Column(column_name, delimiter='|')
        return Map(col=column, extract=[mi_regex])

    def pubmed_annotation(column_name: str) -> CV:
        column = Column(column_name, delimiter='|')
        return CV(term=IdentifierNamespaceCv.PUBMED, value=Map(col=column, extract=[pubmed_regex]))

    def tax_cv(column_name: str) -> CV:
        column = Column(column_name, delimiter='|')
        return CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=Map(col=column, extract=[tax_regex]))

    schema = EntityBuilder(
        entity_type=EntityTypeCv.INTERACTION,
        identifiers=IdentifiersBuilder(
            interaction_identifier_cv(),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            mi_term_cv('Interaction type(s)'),
            mi_term_cv('Interaction detection method(s)'),
            mi_term_cv('Causal statement'),
            mi_term_cv('Causal Regulatory Mechanism'),
            pubmed_annotation('Publication Identifier(s)'),
        ),
        membership=MembershipBuilder(
            Member(
                entity=EntityBuilder(
                    entity_type=mi_term_string('Type(s) interactor A'),
                    identifiers=IdentifiersBuilder(
                        general_identifier_cv('\ufeff#ID(s) interactor A'),
                        general_identifier_cv('Alt. ID(s) interactor A'),
                    ),
                    annotations=AnnotationsBuilder(tax_cv('Taxid interactor A')),
                ),
                annotations=AnnotationsBuilder(
                    mi_term_cv('Biological role(s) interactor A'),
                    mi_term_cv('Experimental role(s) interactor A'),
                ),
            ),
            Member(
                entity=EntityBuilder(
                    entity_type=mi_term_string('Type(s) interactor B'),
                    identifiers=IdentifiersBuilder(
                        general_identifier_cv('ID(s) interactor B'),
                        general_identifier_cv('Alt. ID(s) interactor B'),
                    ),
                    annotations=AnnotationsBuilder(tax_cv('Taxid interactor B')),
                ),
                annotations=AnnotationsBuilder(
                    mi_term_cv('Biological role(s) interactor B'),
                    mi_term_cv('Experimental role(s) interactor B'),
                ),
            ),
        ),
    )
    reader = csv.DictReader(opener.result, delimiter='\t')
    for entity in reader:
        #print(entity)
        yield schema(entity)
