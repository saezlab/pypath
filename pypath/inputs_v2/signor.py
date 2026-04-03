"""
Parse SIGNOR data and emit Entity records.

This module converts SIGNOR data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv
import re

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
_GENERAL_VALUE_REGEX = r'^[^:]+:\"?([^|\"]+)\"?'
_INTERACTION_VALUE_REGEX = r'^[^:]+:(.*)'
# Tax fields can look like:
# taxid:9606(human)|taxid:9606(Homo sapiens)
# We only need the first taxon number appearing in the field.
_TAX_REGEX = r'(-?\d+)'
_PUBMED_REGEX = r'(?i)pubmed:(\d+)'
_MI_REGEX = r'(MI:\d+)'

_TERM_MAPPING = {
    'signor': IdentifierNamespaceCv.SIGNOR,
    'signor-interaction': IdentifierNamespaceCv.SIGNOR,
}

_INTERACTOR_TYPE_MAPPING = {
    'MI:0326': EntityTypeCv.PROTEIN,
    'MI:0314': EntityTypeCv.COMPLEX,
    'MI:0328': EntityTypeCv.SMALL_MOLECULE,
    'MI:2261': EntityTypeCv.PHENOTYPE,
    'MI:2260': EntityTypeCv.STIMULUS,
    'MI:2258': EntityTypeCv.SMALL_MOLECULE,  # xenobiotic
    'MI:1304': EntityTypeCv.PROTEIN_FAMILY,  # molecule set
}

# The SIGNOR complex and protein-family exports currently used here do not
# expose an organism column. The downloaded tables contain human UniProt
# accessions, so we tax-scope these member proteins to human.
SIGNOR_DEFAULT_TAX_ID = '9606'

f = FieldConfig(
    extract={
        'prefix_lower': [_PREFIX_REGEX, str.lower],
        'general_value': _GENERAL_VALUE_REGEX,
        'interaction_value': _INTERACTION_VALUE_REGEX,
        'tax': _TAX_REGEX,
        'pubmed': _PUBMED_REGEX,
        'mi': _MI_REGEX,
        'signor_member': r'^(SIGNOR-[A-Z0-9-]+)$',
        'uniprot_member': r'^((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?)$',
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


def _normalize_signor_identifier(prefix: str, value: str) -> tuple[object | None, str | None]:
    prefix = prefix.lower()
    value = value.strip().strip('"')
    lower_value = value.lower()

    if value.startswith('URS'):
        return IdentifierNamespaceCv.RNACENTRAL, value
    if value.startswith('SIGNOR-'):
        return IdentifierNamespaceCv.SIGNOR, value
    if value.startswith('CHEBI:'):
        return IdentifierNamespaceCv.CHEBI, value
    if value.startswith('DB') and value[2:].isdigit():
        return IdentifierNamespaceCv.DRUGBANK, value
    if value.startswith('PUBCHEM:'):
        return IdentifierNamespaceCv.PUBCHEM, value.split(':', 1)[1]
    if prefix == 'pubchem' and value.upper().startswith('CID:'):
        return IdentifierNamespaceCv.PUBCHEM, value.split(':', 1)[1]
    if prefix == 'pubchem' and value.upper().startswith('SID:'):
        return None, None
    if value.startswith('uniprotkb:'):
        return _normalize_signor_identifier('uniprotkb', value.split(':', 1)[1])
    if prefix == 'uniprotkb' and '-PRO_' in value:
        value = value.split('-PRO_', 1)[0]
    if prefix == 'uniprotkb' and value.startswith('SIGNOR-'):
        return IdentifierNamespaceCv.SIGNOR, value
    if prefix == 'uniprotkb' and value.startswith('URS'):
        return IdentifierNamespaceCv.RNACENTRAL, value

    return _IDENTIFIER_CV_MAPPING.get(prefix), value or None


def _parse_signor_identifier_pairs(raw: object) -> list[tuple[object, str]]:
    pairs: list[tuple[object, str]] = []
    for item in _split_signor_field(raw):
        if ':' not in item:
            continue
        prefix, value = item.split(':', 1)
        mapped, normalized_value = _normalize_signor_identifier(prefix, value)
        if mapped is not None and normalized_value:
            pairs.append((mapped, normalized_value))
    return pairs


def general_identifier_cv(column_name: str) -> CV:
    return CV(
        term=lambda row: [term for term, _ in _parse_signor_identifier_pairs(row.get(column_name))],
        value=lambda row: [value for _, value in _parse_signor_identifier_pairs(row.get(column_name))],
    )


def mi_term_cv(column_name: str) -> CV:
    return CV(term=f(column_name, extract='mi'))


def _infer_signor_interactor_type(row: dict[str, object], suffix: str) -> EntityTypeCv:
    raw_type = str(row.get(f'Type(s) interactor {suffix}') or '').strip()
    if raw_type and raw_type != '-':
        mi_match = re.search(_MI_REGEX, raw_type)
        if mi_match:
            mapped = _INTERACTOR_TYPE_MAPPING.get(mi_match.group(1))
            if mapped is not None:
                return mapped

    id_fields = [
        f'\ufeff#ID(s) interactor {suffix}' if suffix == 'A' else f'ID(s) interactor {suffix}',
        f'Alt. ID(s) interactor {suffix}',
        f'Alias(es) interactor {suffix}',
    ]
    values: list[str] = []
    for field in id_fields:
        raw = row.get(field)
        values.extend(_split_signor_field(raw))

    lower_values = [value.lower() for value in values]

    if any(value.startswith('chebi:') or value.startswith('pubchem:') for value in lower_values):
        return EntityTypeCv.SMALL_MOLECULE
    if any('signor-c' in value for value in lower_values):
        return EntityTypeCv.COMPLEX
    if any('signor-pf' in value or 'signor-fp' in value for value in lower_values):
        return EntityTypeCv.PROTEIN_FAMILY
    if any('mir' in value or 'urs' in value for value in lower_values):
        return EntityTypeCv.RNA
    if any(value.startswith('uniprotkb:') for value in lower_values):
        return EntityTypeCv.PROTEIN
    if any(value.startswith('signor:') for value in lower_values):
        return EntityTypeCv.PHYSICAL_ENTITY

    return EntityTypeCv.PHYSICAL_ENTITY


def interactor_entity_type(suffix: str):
    return lambda row: _infer_signor_interactor_type(row, suffix)


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


def _split_signor_field(raw: object) -> list[str]:
    if raw is None:
        return []
    text = str(raw).strip()
    if not text or text == '-':
        return []
    return [part.strip() for part in text.split('|') if part and part.strip() and part.strip() != '-']


def _clean_signor_value(value: str) -> str | None:
    cleaned = value.strip().strip('"').strip()
    return cleaned or None


def _parse_term_value_annotations(
    raw: object,
    *,
    default_term: str | None = None,
) -> list[tuple[str, str | None]]:
    annotations: list[tuple[str, str | None]] = []
    for item in _split_signor_field(raw):
        if ':' in item:
            term, value = item.split(':', 1)
            parsed_term = _clean_signor_value(term)
            parsed_value = _clean_signor_value(value)
            if parsed_term:
                annotations.append((parsed_term, parsed_value))
        elif default_term:
            parsed_value = _clean_signor_value(item)
            if parsed_value:
                annotations.append((default_term, parsed_value))
    return annotations


def parsed_annotation_terms(column_name: str, *, default_term: str | None = None):
    return lambda row: [
        term
        for term, _ in _parse_term_value_annotations(row.get(column_name), default_term=default_term)
    ]


def parsed_annotation_values(column_name: str, *, default_term: str | None = None):
    return lambda row: [
        value
        for _, value in _parse_term_value_annotations(row.get(column_name), default_term=default_term)
    ]


config = ResourceConfig(
    id=ResourceCv.SIGNOR,
    name='SIGNOR',
    url='https://signor.uniroma2.it/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='31665520',
    primary_category='interactions',
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
                    value=f('LIST OF ENTITIES', delimiter=',', extract='uniprot_member'),
                ),
                CV(
                    term=IdentifierNamespaceCv.SIGNOR,
                    value=f('LIST OF ENTITIES', delimiter=',', extract='signor_member'),
                ),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=SIGNOR_DEFAULT_TAX_ID),
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
                    value=f('LIST OF ENTITIES', delimiter=',', extract='uniprot_member'),
                ),
                CV(
                    term=IdentifierNamespaceCv.SIGNOR,
                    value=f('LIST OF ENTITIES', delimiter=',', extract='signor_member'),
                ),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=SIGNOR_DEFAULT_TAX_ID),
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
        CV(
            term=parsed_annotation_terms(
                'Interaction annotation(s)',
                default_term='signor:interaction_annotation',
            ),
            value=parsed_annotation_values(
                'Interaction annotation(s)',
                default_term='signor:interaction_annotation',
            ),
        ),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=interactor_entity_type('A'),
                identifiers=IdentifiersBuilder(
                    general_identifier_cv('\ufeff#ID(s) interactor A'),
                    general_identifier_cv('Alt. ID(s) interactor A'),
                ),
                annotations=AnnotationsBuilder(tax_cv('Taxid interactor A')),
            ),
            annotations=AnnotationsBuilder(
                mi_term_cv('Biological role(s) interactor A'),
                mi_term_cv('Experimental role(s) interactor A'),
                CV(
                    term=parsed_annotation_terms('Feature(s) interactor A'),
                    value=parsed_annotation_values('Feature(s) interactor A'),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=interactor_entity_type('B'),
                identifiers=IdentifiersBuilder(
                    general_identifier_cv('ID(s) interactor B'),
                    general_identifier_cv('Alt. ID(s) interactor B'),
                ),
                annotations=AnnotationsBuilder(tax_cv('Taxid interactor B')),
            ),
            annotations=AnnotationsBuilder(
                mi_term_cv('Biological role(s) interactor B'),
                mi_term_cv('Experimental role(s) interactor B'),
                CV(
                    term=parsed_annotation_terms('Feature(s) interactor B'),
                    value=parsed_annotation_values('Feature(s) interactor B'),
                ),
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
