"""
Parse SIGNOR data and emit Entity records.

This module converts SIGNOR data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv

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
    Member,
    MembershipBuilder,
    MembersFromList,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig


def _iter_semicolon(opener, **_kwargs: object):
    if not opener or not opener.result:
        return
    yield from csv.DictReader(opener.result, delimiter=';')


def _iter_tsv(opener, **_kwargs: object):
    if not opener or not opener.result:
        return
    yield from csv.DictReader(opener.result, delimiter='\t')


_IDENTIFIER_CV_MAPPING = {
    'chebi': IdentifierNamespaceCv.CHEBI,
    'complexportal': IdentifierNamespaceCv.COMPLEXPORTAL,
    'pubchem': IdentifierNamespaceCv.PUBCHEM,
    'signor': IdentifierNamespaceCv.SIGNOR,
    'uniprotkb': IdentifierNamespaceCv.UNIPROT,
}

_PREFIX_REGEX = r'^([^:]+):'
_GENERAL_VALUE_REGEX = r'^[^:]+:([^|"]+)'
_INTERACTION_VALUE_REGEX = r'^[^:]+:(.*)'
_TAX_REGEX = r'taxid:([-\d]+)'
_PUBMED_REGEX = r'(?i)pubmed:(\d+)'
_MI_REGEX = r'(MI:\d+)'

_TERM_MAPPING = {
    'signor': IdentifierNamespaceCv.SIGNOR,
    'signor-interaction': IdentifierNamespaceCv.SIGNOR,
}

f = FieldConfig(
    extract={
        'prefix_lower': [_PREFIX_REGEX, str.lower],
        'general_value': _GENERAL_VALUE_REGEX,
        'interaction_value': _INTERACTION_VALUE_REGEX,
        'tax': _TAX_REGEX,
        'pubmed': _PUBMED_REGEX,
        'mi': _MI_REGEX,
    },
    map={
        'identifier_cv': _IDENTIFIER_CV_MAPPING,
        'term_cv': _TERM_MAPPING,
    },
    delimiter='|',
)


def interaction_identifier_cv() -> CV:
    return CV(
        term=f('Interaction identifier(s)', extract='prefix_lower', map='term_cv'),
        value=f('Interaction identifier(s)', extract='interaction_value'),
    )


def general_identifier_cv(column_name: str) -> CV:
    return CV(
        term=f(column_name, extract='prefix_lower', map='identifier_cv'),
        value=f(column_name, extract='general_value'),
    )


def mi_term_cv(column_name: str) -> CV:
    return CV(term=f(column_name, extract='mi'))


def mi_term_string(column_name: str):
    return f(column_name, extract='mi')


def pubmed_annotation(column_name: str) -> CV:
    return CV(
        term=IdentifierNamespaceCv.PUBMED,
        value=f(column_name, extract='pubmed'),
    )


def tax_cv(column_name: str) -> CV:
    return CV(
        term=IdentifierNamespaceCv.NCBI_TAX_ID,
        value=f(column_name, extract='tax'),
    )


config = ResourceConfig(
    id=ResourceCv.SIGNOR,
    name='SIGNOR',
    url='https://signor.uniroma2.it/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='31665520',
    description=(
        'SIGNOR (SIGnaling Network Open Resource) is a comprehensive '
        'resource of causal relationships between biological entities '
        'with a focus on signaling pathways. It provides manually curated '
        'interactions with mechanistic details including protein-protein '
        'interactions, post-translational modifications, transcriptional '
        'regulation, and small molecule effects.'
    ),
)

complexes_schema = EntityBuilder(
    entity_type=EntityTypeCv.COMPLEX,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.SIGNOR, value=f('SIGNOR ID')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('COMPLEX NAME')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.PROTEIN,
            identifiers=IdentifiersBuilder(
                CV(
                    term=IdentifierNamespaceCv.UNIPROT,
                    value=f('LIST OF ENTITIES', delimiter=','),
                ),
            ),
        )
    ),
)

protein_families_schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN_FAMILY,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.SIGNOR, value=f('SIGNOR ID')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('PROT. FAMILY NAME')),
    ),
    annotations=AnnotationsBuilder(),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.PROTEIN,
            identifiers=IdentifiersBuilder(
                CV(
                    term=IdentifierNamespaceCv.UNIPROT,
                    value=f('LIST OF ENTITIES', delimiter=','),
                ),
            ),
        )
    ),
)

phenotypes_schema = EntityBuilder(
    entity_type=EntityTypeCv.PHENOTYPE,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.SIGNOR, value=f('SIGNOR ID')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('PHENOTYPE NAME')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('PHENOTYPE DESCRIPTION')),
    ),
)

stimuli_schema = EntityBuilder(
    entity_type=EntityTypeCv.STIMULUS,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.SIGNOR, value=f('SIGNOR ID')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('STIMULUS NAME')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('STIMULUS DESCRIPTION')),
    ),
)

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        interaction_identifier_cv(),
    ),
    annotations=AnnotationsBuilder(
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

resource = Resource(
    config,
    complexes=Dataset(
        download=Download(
            url='https://signor.uniroma2.it/download_complexes.php',
            filename='signor_complexes.txt',
            subfolder='signor',
            download_kwargs={'query': {'submit': 'Download complex data'}, 'post': True},
        ),
        mapper=complexes_schema,
        raw_parser=_iter_semicolon,
    ),
    protein_families=Dataset(
        download=Download(
            url='https://signor.uniroma2.it/download_complexes.php',
            filename='signor_protein_families.txt',
            subfolder='signor',
            download_kwargs={'query': {'submit': 'Download protein family data'}, 'post': True},
        ),
        mapper=protein_families_schema,
        raw_parser=_iter_semicolon,
    ),
    phenotypes=Dataset(
        download=Download(
            url='https://signor.uniroma2.it/download_complexes.php',
            filename='signor_phenotypes.txt',
            subfolder='signor',
            download_kwargs={'query': {'submit': 'Download phenotype data'}, 'post': True},
        ),
        mapper=phenotypes_schema,
        raw_parser=_iter_semicolon,
    ),
    stimuli=Dataset(
        download=Download(
            url='https://signor.uniroma2.it/download_complexes.php',
            filename='signor_stimuli.txt',
            subfolder='signor',
            download_kwargs={'query': {'submit': 'Download stimulus data'}, 'post': True},
        ),
        mapper=stimuli_schema,
        raw_parser=_iter_semicolon,
    ),
    interactions=Dataset(
        download=Download(
            url='https://signor.uniroma2.it/download_entity.php',
            filename='signor_all_causalTab.txt',
            subfolder='signor',
            download_kwargs={'query': {'format': 'causalTab', 'submit': 'Download'}, 'post': True},
        ),
        mapper=interactions_schema,
        raw_parser=_iter_tsv,
    ),
)
