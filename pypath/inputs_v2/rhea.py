"""Rhea master reactions as molecular activities with explicit input, output and enzyme edges.

Transport sides (in/out) and reported conversion direction attach to membership
evidence as ordinary attributes, without gene-effect direction qualifiers.
Reaction identifiers, equations, EC topics, publications and cross-references
are preserved.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Mapping
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import TRANSPORT_SIDE

from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembershipBuilder,
    MembersFromList,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2._source_context import conversion_direction_cv
from pypath.inputs_v2.parsers.rhea import _raw


config = ResourceConfig(
    id=ResourceCv.RHEA,
    name='Rhea',
    url='https://www.rhea-db.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='34755880',
    primary_category='reactions',
    description=(
        'Rhea is an expert-curated knowledgebase of chemical and transport '
        'reactions of biological interest.'
    ),
)


_direction_map = {
    'LR': 'LEFT-TO-RIGHT',
    'RL': 'RIGHT-TO-LEFT',
    'BI': 'REVERSIBLE',
    'UN': None,
}

_role_map = {
    'reactant': slots.has_input,
    'product': slots.has_output,
}

_CHEBI_RE = re.compile(r'^(?:CHEBI:)?(\d+)$')
_RHEA_ID_RE = re.compile(r'(?:RHEA:)?(\d+)')


_PARTICIPANT_FIELDS = (
    'participant_role', 'participant_chebi',
    'participant_display_name', 'participant_compartment',
)


def _participant_tokens(row: Mapping[str, Any], field: str) -> list[str | None]:
    value = row.get(field)
    if value is None:
        return []
    values = value if isinstance(value, (list, tuple)) else str(value).split('||')
    return [str(item).strip().strip('"') if item is not None else None for item in values]


def _participant_field(
    field: str,
    transform: Callable[[str], Any] | None = None,
) -> Callable[[Mapping[str, Any]], list[Any]]:
    """Pad source arrays so a one-item field cannot broadcast to other members."""

    def extract(row: Mapping[str, Any]) -> list[Any]:
        arrays = {name: _participant_tokens(row, name) for name in _PARTICIPANT_FIELDS}
        size = max(map(len, arrays.values()), default=0)
        values = arrays[field] + [None] * (size - len(arrays[field]))
        return [
            transform(value) if value and transform else value or None
            for value in values
        ]

    return extract


def _participant_chebi(value: str) -> str | None:
    match = _CHEBI_RE.fullmatch(value)
    return match[1] if match else None


def _reaction_associations() -> AssociationsBuilder:
    return AssociationsBuilder(
        AssociationBuilder(
            object_entity_type=model.OntologyClass,
            object_identifier_type=Namespace.GO,
            object_identifier=f('go'),
        ),
        AssociationBuilder(
            object_entity_type=model.MolecularActivity,
            object_identifier_type=Namespace.REACTOME,
            object_identifier=f('reactome'),
        ),
    )


# ── reactions ─────────────────────────────────────────────────────────────────

f = FieldConfig(
    delimiter=';',
    map={
        'role': lambda value: _role_map.get(value),
    },
)


reactions_download = Download(
    url=(
        'https://www.rhea-db.org/rhea/?query=&columns=rhea-id,equation,chebi,'
        'chebi-id,ec,uniprot,go,pubmed,reaction-xref(EcoCyc),reaction-xref(MetaCyc),'
        'reaction-xref(KEGG),reaction-xref(Reactome),reaction-xref(M-CSA)'
        '&format=tsv&limit=1000000'
    ),
    filename='rhea_reactions.tsv',
    subfolder='rhea',
    ext='tsv',
    default_mode='r',
)


reactions_schema = EntityBuilder(
    entity_type=model.MolecularActivity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.RHEA, value=f('rhea_id')),
        CV(term=Namespace.NAME, value=f('equation')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.has_topic,
            value=lambda row, source=f('ec'): [
                'EC:' + str(v).removeprefix('EC:') for v in source.extract(row)
            ],
        ),
        CV(
            term=slots.publications,
            value=f(
                'pubmed',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:')
                if str(v).removeprefix('PMID:').isdigit()
                and int(str(v).removeprefix('PMID:')) > 0
                else None,
            ),
        ),
        CV(
            term=slots.has_topic,
            value=f(
                'ecocyc',
                transform=lambda v: 'EcoCyc:' + str(v).removeprefix('EcoCyc:'),
            ),
        ),
        CV(
            term=slots.has_topic,
            value=f(
                'metacyc',
                transform=lambda v: 'MetaCyc:'
                + str(v).removeprefix('MetaCyc:'),
            ),
        ),
        CV(
            term=slots.has_topic,
            value=f(
                'kegg',
                transform=lambda v: 'KEGG.REACTION:'
                + str(v).removeprefix('KEGG.REACTION:').removeprefix('KEGG:'),
            ),
        ),
    ),
    associations=_reaction_associations(),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=model.ChemicalEntity,
            predicate=_participant_field('participant_role', _role_map.get),
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.CHEBI,
                    value=_participant_field('participant_chebi', _participant_chebi),
                ),
                CV(
                    term=Namespace.NAME,
                    value=_participant_field('participant_display_name'),
                ),
            ),
            annotations=AnnotationsBuilder(
                conversion_direction_cv(),
                CV(
                    term=TRANSPORT_SIDE,
                    value=_participant_field('participant_compartment'),
                ),
            ),
            entity_annotations=AnnotationsBuilder(),
        ),
        MembersFromList(
            entity_type=model.Protein,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.UNIPROT, value=f('uniprot'))
            ),
            annotations=AnnotationsBuilder(conversion_direction_cv()),
            predicate=slots.enabled_by,
        ),
    ),
)


# ── transport_reactions ───────────────────────────────────────────────────────

transport_reactions_schema = EntityBuilder(
    entity_type=model.MolecularActivity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.RHEA, value=f('rhea_id')),
        CV(term=Namespace.NAME, value=f('equation')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.has_topic,
            value=lambda row, source=f('ec'): [
                'EC:' + str(v).removeprefix('EC:') for v in source.extract(row)
            ],
        ),
        CV(
            term=slots.publications,
            value=f(
                'pubmed',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:')
                if str(v).removeprefix('PMID:').isdigit()
                and int(str(v).removeprefix('PMID:')) > 0
                else None,
            ),
        ),
        CV(
            term=slots.has_topic,
            value=f(
                'ecocyc',
                transform=lambda v: 'EcoCyc:' + str(v).removeprefix('EcoCyc:'),
            ),
        ),
        CV(
            term=slots.has_topic,
            value=f(
                'metacyc',
                transform=lambda v: 'MetaCyc:'
                + str(v).removeprefix('MetaCyc:'),
            ),
        ),
        CV(
            term=slots.has_topic,
            value=f(
                'kegg',
                transform=lambda v: 'KEGG.REACTION:'
                + str(v).removeprefix('KEGG.REACTION:').removeprefix('KEGG:'),
            ),
        ),
    ),
    associations=_reaction_associations(),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=model.ChemicalEntity,
            predicate=_participant_field('participant_role', _role_map.get),
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.CHEBI,
                    value=_participant_field('participant_chebi', _participant_chebi),
                ),
                CV(
                    term=Namespace.NAME,
                    value=_participant_field('participant_display_name'),
                ),
            ),
            annotations=AnnotationsBuilder(
                conversion_direction_cv(),
                CV(
                    term=TRANSPORT_SIDE,
                    value=_participant_field('participant_compartment'),
                ),
            ),
            entity_annotations=AnnotationsBuilder(),
        ),
        MembersFromList(
            entity_type=model.Protein,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.UNIPROT, value=f('uniprot'))
            ),
            annotations=AnnotationsBuilder(conversion_direction_cv()),
            predicate=slots.enabled_by,
        ),
    ),
)


# ── catalysis ─────────────────────────────────────────────────────────────────

g = FieldConfig(
    map={
        'direction': lambda value: _direction_map.get(value),
    },
)


catalysis_download = Download(
    url='https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv',
    filename='rhea2uniprot.tsv',
    subfolder='rhea',
    ext='tsv',
    default_mode='r',
)


# ── resource ─────────────────────────────────────────────────────────────────

resource = Resource(
    config,
    reactions=Dataset(
        download=reactions_download,
        mapper=reactions_schema,
        raw_parser=lambda opener, force_refresh=False, **kwargs: _raw(
            opener,
            uniprot_opener=catalysis_download.open(force_refresh=force_refresh),
            force_refresh=force_refresh,
            **kwargs,
        ),
    ),
    metabolic_reactions=Dataset(
        download=reactions_download,
        mapper=reactions_schema,
        raw_parser=lambda opener, force_refresh=False, **kwargs: _raw(
            opener,
            data_type='metabolic_reactions',
            uniprot_opener=catalysis_download.open(force_refresh=force_refresh),
            force_refresh=force_refresh,
            **kwargs,
        ),
    ),
    transport_reactions=Dataset(
        download=reactions_download,
        mapper=transport_reactions_schema,
        raw_parser=lambda opener, force_refresh=False, **kwargs: _raw(
            opener,
            data_type='transport_reactions',
            uniprot_opener=catalysis_download.open(force_refresh=force_refresh),
            force_refresh=force_refresh,
            **kwargs,
        ),
    ),
)
