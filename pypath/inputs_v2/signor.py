"""
Parse SIGNOR data and emit Entity records.

This module converts SIGNOR data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

import csv
import re

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
    RelationBuilder,
)


def _iter_semicolon(opener, **_kwargs: object):
    if not opener or not opener.result:
        return
    yield from csv.DictReader(opener.result, delimiter=';')


def _iter_tsv(opener, **_kwargs: object):
    if not opener or not opener.result:
        return
    yield from csv.DictReader(opener.result, delimiter='\t')


_IDENTIFIER_CV_MAPPING = {
    'chebi': Namespace.CHEBI,
    'complexportal': Namespace.COMPLEXPORTAL,
    'pubchem': Namespace.PUBCHEM,
    'signor': Namespace.SIGNOR,
    'uniprotkb': Namespace.UNIPROT,
    'uniprot': Namespace.UNIPROT,
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
    'signor': Namespace.SIGNOR,
    'signor-interaction': Namespace.SIGNOR,
}

_INTERACTOR_TYPE_MAPPING = {
    'MI:0326': model.Protein,
    'MI:0314': model.MacromolecularComplex,
    'MI:0328': model.ChemicalEntity,
    'MI:2261': model.PhenotypicFeature,
    'MI:2260': model.ExposureEvent,
    'MI:2258': model.ChemicalEntity,  # xenobiotic does not imply molecular size
    'MI:1304': model.ProteinFamily,  # molecule set
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
        term=f(
            'Interaction identifier(s)', extract='prefix_lower', map='term_cv'
        ),
        value=f('Interaction identifier(s)', extract='interaction_value'),
    )


def _normalize_signor_identifier(
    prefix: str, value: str
) -> tuple[object | None, str | None]:
    prefix = prefix.lower()
    value = value.strip().strip('"')
    lower_value = value.lower()

    if value.startswith('URS'):
        return Namespace.RNACENTRAL, value
    if value.startswith('SIGNOR-'):
        return Namespace.SIGNOR, value
    if lower_value.startswith('chebi:'):
        match = re.fullmatch(r'CHEBI:(\d+)', value, flags=re.IGNORECASE)
        return Namespace.CHEBI, match.group(1) if match else None
    if value.startswith('DB') and value[2:].isdigit():
        return Namespace.DRUGBANK, value
    if lower_value.startswith('pubchem:sid:'):
        return Namespace.PUBCHEM_SUBSTANCE, value.rsplit(':', 1)[1]
    if lower_value.startswith('pubchem:cid:'):
        return Namespace.PUBCHEM, value.rsplit(':', 1)[1]
    if lower_value.startswith('pubchem:'):
        return Namespace.PUBCHEM, value.split(':', 1)[1]
    if prefix == 'pubchem' and value.upper().startswith('CID:'):
        return Namespace.PUBCHEM, value.split(':', 1)[1]
    if prefix == 'pubchem' and value.upper().startswith('SID:'):
        return Namespace.PUBCHEM_SUBSTANCE, value.split(':', 1)[1]
    if lower_value.startswith('uniprotkb:'):
        return _normalize_signor_identifier('uniprotkb', value.split(':', 1)[1])
    if prefix == 'uniprotkb' and '-PRO_' in value:
        value = value.split('-PRO_', 1)[0]
    if prefix == 'uniprotkb' and value.startswith('SIGNOR-'):
        return Namespace.SIGNOR, value
    if prefix == 'uniprotkb' and value.startswith('URS'):
        return Namespace.RNACENTRAL, value

    return _IDENTIFIER_CV_MAPPING.get(prefix), value or None


def _parse_signor_identifier_pairs(raw: object) -> list[tuple[object, str]]:
    return list(_identifier_pairs_from_text('' if raw is None else str(raw)))


def _identifier_pairs_from_text(text: str) -> tuple[tuple[object, str], ...]:
    pairs: list[tuple[object, str]] = []
    for item in _split_signor_field(text):
        if ':' not in item:
            continue
        prefix, value = item.split(':', 1)
        mapped, normalized_value = _normalize_signor_identifier(prefix, value)
        if mapped is not None and normalized_value:
            pairs.append((mapped, normalized_value))
    return tuple(pairs)


def general_identifier_cv(column_name: str) -> CV:
    return CV.from_pairs(
        f(column_name, transform=_identifier_pairs_from_text),
    )


def _infer_signor_interactor_type(
    row: dict[str, object], suffix: str
) -> type[model.NamedThing]:
    """Read the source type; repair only documented missing-type cases.

    The cached export has 533 untyped endpoints: RNAcentral identifiers mislabeled
    UniProt, SIGNOR FP families, DrugBank entities and a few small molecules.
    See docs/biolink-boundaries.md for the evidence and modeling decisions.
    """
    raw_type = str(row.get(f'Type(s) interactor {suffix}') or '').strip()
    if raw_type and raw_type != '-':
        match = re.search(_MI_REGEX, raw_type)
        if match and match.group(1) in _INTERACTOR_TYPE_MAPPING:
            return _INTERACTOR_TYPE_MAPPING[match.group(1)]
        raise ValueError(f'Unknown SIGNOR interactor type: {raw_type!r}')

    primary = next(
        (
            row.get(key)
            for key in (
                f'\ufeff#ID(s) interactor {suffix}',
                f'#ID(s) interactor {suffix}',
                f'ID(s) interactor {suffix}',
            )
            if row.get(key)
        ),
        None,
    )
    pairs = _parse_signor_identifier_pairs(primary)
    types = set()
    for namespace, identifier in pairs:
        if namespace == Namespace.RNACENTRAL:
            types.add(model.RNAProduct)
        elif namespace in (
            Namespace.DRUGBANK,
            Namespace.CHEBI,
            Namespace.PUBCHEM,
        ):
            types.add(model.ChemicalEntity)
        elif namespace == Namespace.SIGNOR and re.fullmatch(
            r'SIGNOR-(?:PF|FP)\d+', identifier
        ):
            types.add(model.ProteinFamily)
    if len(types) == 1:
        return types.pop()
    raise ValueError(
        f'SIGNOR endpoint {suffix} has no supported type: {primary!r}'
    )


def interactor_entity_type(suffix: str):
    return lambda row: _infer_signor_interactor_type(row, suffix)


def interactor_tax_cv(suffix: str) -> CV:
    return CV(
        term=slots.in_taxon,
        value=f(
            f'Taxid interactor {suffix}',
            extract='tax',
            transform=lambda taxon: (
                f'NCBITaxon:{taxon}'
                if str(taxon).isdigit() and int(taxon) > 0
                else None
            ),
        ),
    )


def pubmed_annotation(column_name: str) -> CV:
    return CV(
        term=slots.publications,
        value=f(
            column_name, extract='pubmed', transform=lambda pmid: f'PMID:{pmid}'
        ),
    )


def _split_signor_field(raw: object) -> list[str]:
    if raw is None:
        return []
    text = str(raw).strip()
    if not text or text == '-':
        return []
    parts: list[str] = []
    current: list[str] = []
    in_quotes = False
    for char in text:
        if char == '"':
            in_quotes = not in_quotes
        if char == '|' and not in_quotes:
            part = ''.join(current).strip()
            if part and part != '-':
                parts.append(part)
            current = []
            continue
        current.append(char)
    part = ''.join(current).strip()
    if part and part != '-':
        parts.append(part)
    return parts


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

def _group_member_type(value):
    if re.fullmatch(r'CHEBI:\d+', value):
        return model.ChemicalEntity
    if re.fullmatch(r'SIGNOR-C\d+', value):
        return model.MacromolecularComplex
    if re.fullmatch(r'SIGNOR-(?:PF|FP)\d+', value):
        return model.ProteinFamily
    if f('LIST OF ENTITIES', extract='uniprot_member').extract(
        {'LIST OF ENTITIES': value}
    ):
        return model.Protein
    return model.NamedThing if value.startswith('SIGNOR-') else None


def _group_members():
    return MembershipBuilder(
        MembersFromList(
            entity_type=f(
                'LIST OF ENTITIES', delimiter=',',
                transform=_group_member_type, preserve_indices=True,
            ),
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.UNIPROT, value=f(
                    'LIST OF ENTITIES', delimiter=',',
                    extract='uniprot_member', preserve_indices=True,
                )),
                CV(term=Namespace.SIGNOR, value=f(
                    'LIST OF ENTITIES', delimiter=',',
                    extract='signor_member', preserve_indices=True,
                )),
                CV(term=Namespace.CHEBI, value=f(
                    'LIST OF ENTITIES', delimiter=',',
                    extract=r'^CHEBI:(\d+)$', preserve_indices=True,
                )),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=slots.in_taxon,
                   value=f(
                       'LIST OF ENTITIES', delimiter=',', preserve_indices=True,
                       transform=lambda value: (
                           None if _group_member_type(value) == model.ChemicalEntity
                           else f'NCBITaxon:{SIGNOR_DEFAULT_TAX_ID}'
                       ),
                   )),
            ),
        ),
    )


complexes_schema = EntityBuilder(
    entity_type=model.MacromolecularComplex,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.SIGNOR, value=f('SIGNOR ID')),
        CV(term=Namespace.NAME, value=f('COMPLEX NAME')),
    ),
    annotations=AnnotationsBuilder(),
    membership=_group_members(),
)

protein_families_schema = EntityBuilder(
    entity_type=model.ProteinFamily,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.SIGNOR, value=f('SIGNOR ID')),
        CV(term=Namespace.NAME, value=f('PROT. FAMILY NAME')),
    ),
    annotations=AnnotationsBuilder(),
    membership=_group_members(),
)

phenotypes_schema = EntityBuilder(
    entity_type=model.PhenotypicFeature,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.SIGNOR, value=f('SIGNOR ID')),
        CV(term=Namespace.NAME, value=f('PHENOTYPE NAME')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.description, value=f('PHENOTYPE DESCRIPTION')),
    ),
)

stimuli_schema = EntityBuilder(
    entity_type=model.ExposureEvent,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.SIGNOR, value=f('SIGNOR ID')),
        CV(term=Namespace.NAME, value=f('STIMULUS NAME')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.description, value=f('STIMULUS DESCRIPTION')),
    ),
)

# Source ontology crosswalk, not a new vocabulary. The source's causal
# accessions encode direction and aspect; all targets are native Biolink enums.
_DIRECTION = model.DirectionQualifierEnum
_ASPECT = model.GeneOrGeneProductOrChemicalEntityAspectEnum
_CAUSAL_QUALIFIERS = {
    'MI:2235': (_DIRECTION.increased, None),
    'MI:2236': (_DIRECTION.increased, _ASPECT.activity),
    'MI:2237': (_DIRECTION.increased, _ASPECT.abundance),
    'MI:2238': (_DIRECTION.increased, _ASPECT.expression),
    'MI:2239': (_DIRECTION.increased, _ASPECT.stability),
    'MI:2240': (_DIRECTION.decreased, None),
    'MI:2241': (_DIRECTION.decreased, _ASPECT.activity),
    'MI:2242': (_DIRECTION.decreased, _ASPECT.abundance),
    'MI:2243': (_DIRECTION.decreased, _ASPECT.stability),
    'MI:2244': (_DIRECTION.decreased, _ASPECT.expression),
}


def _causal_accession(row):
    raw = str(row.get('Causal statement') or '')
    match = re.search(_MI_REGEX, raw)
    if not match or match.group(1) not in _CAUSAL_QUALIFIERS:
        raise ValueError(f'Unsupported SIGNOR causal statement: {raw!r}')
    return match.group(1)


def signor_predicate(row):
    _causal_accession(row)
    # SIGNOR includes exogenous chemical effects: do not assert evolved regulation
    # for every record. Direction/aspect refine this model-native affects edge.
    return slots.affects


def _participant_builder(suffix):
    return EntityBuilder(
        entity_type=interactor_entity_type(suffix),
        identifiers=IdentifiersBuilder(
            general_identifier_cv(
                '\ufeff#ID(s) interactor A'
                if suffix == 'A'
                else 'ID(s) interactor B'
            ),
            general_identifier_cv(f'Alt. ID(s) interactor {suffix}'),
        ),
        annotations=AnnotationsBuilder(
            interactor_tax_cv(suffix),
            CV(
                term=slots.description,
                value=f(f'Feature(s) interactor {suffix}'),
            ),
        ),
    )


interactor_a_builder = _participant_builder('A')
interactor_b_builder = _participant_builder('B')
interactions_schema = RelationBuilder(
    subject=interactor_a_builder,
    predicate=signor_predicate,
    object=interactor_b_builder,
    identifiers=IdentifiersBuilder(interaction_identifier_cv()),
    annotations=AnnotationsBuilder(
        # The source causal accession remains in raw evidence, not a mapping-only Biolink slot.
        CV(
            term=slots.object_direction_qualifier,
            value=lambda row: _CAUSAL_QUALIFIERS[_causal_accession(row)][0],
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=lambda row: _CAUSAL_QUALIFIERS[_causal_accession(row)][1],
        ),
        CV(
            term=slots.has_evidence_of_type,
            value=f('Interaction detection method(s)', extract='mi'),
        ),
        pubmed_annotation('Publication Identifier(s)'),
        CV(term=slots.description, value=f('Interaction annotation(s)')),
    ),
)

resource = Resource(
    config,
    complexes=Dataset(
        download=Download(
            url='https://signor.uniroma2.it/download_complexes.php',
            filename='signor_complexes.txt',
            subfolder='signor',
            download_kwargs={
                'query': {'submit': 'Download complex data'},
                'post': True,
            },
        ),
        mapper=complexes_schema,
        raw_parser=_iter_semicolon,
    ),
    protein_families=Dataset(
        download=Download(
            url='https://signor.uniroma2.it/download_complexes.php',
            filename='signor_protein_families.txt',
            subfolder='signor',
            download_kwargs={
                'query': {'submit': 'Download protein family data'},
                'post': True,
            },
        ),
        mapper=protein_families_schema,
        raw_parser=_iter_semicolon,
    ),
    phenotypes=Dataset(
        download=Download(
            url='https://signor.uniroma2.it/download_complexes.php',
            filename='signor_phenotypes.txt',
            subfolder='signor',
            download_kwargs={
                'query': {'submit': 'Download phenotype data'},
                'post': True,
            },
        ),
        mapper=phenotypes_schema,
        raw_parser=_iter_semicolon,
    ),
    stimuli=Dataset(
        download=Download(
            url='https://signor.uniroma2.it/download_complexes.php',
            filename='signor_stimuli.txt',
            subfolder='signor',
            download_kwargs={
                'query': {'submit': 'Download stimulus data'},
                'post': True,
            },
        ),
        mapper=stimuli_schema,
        raw_parser=_iter_semicolon,
    ),
    interactions=Dataset(
        download=Download(
            url='https://signor.uniroma2.it/download_entity.php',
            filename='signor_all_causalTab.txt',
            subfolder='signor',
            download_kwargs={
                'query': {'format': 'causalTab', 'submit': 'Download'},
                'post': True,
            },
        ),
        mapper=interactions_schema,
        raw_parser=_iter_tsv,
    ),
)
