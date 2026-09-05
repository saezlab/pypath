"""Rhea master reactions as molecular activities with explicit input, output and enzyme edges.

Transport compartment labels and conversion direction remain in parsed evidence;
these are not gene-effect direction qualifiers. Reaction identifiers and equations,
EC topics, publications and source cross-references are preserved.
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace

import re

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

_CHEBI_RE = re.compile(r'(?:CHEBI:)?(\d+)')
_RHEA_ID_RE = re.compile(r'(?:RHEA:)?(\d+)')


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
                + str(v).removeprefix('KEGG.REACTION:'),
            ),
        ),
    ),
    associations=_reaction_associations(),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=model.ChemicalEntity,
            predicate=f(
                'participant_role',
                delimiter='||',
                map='role',
                preserve_indices=True,
            ),
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'participant_chebi', delimiter='||', extract=_CHEBI_RE
                    ),
                ),
                CV(
                    term=Namespace.NAME,
                    value=f('participant_display_name', delimiter='||'),
                ),
            ),
            annotations=AnnotationsBuilder(),
            entity_annotations=AnnotationsBuilder(),
        ),
        MembersFromList(
            entity_type=model.Protein,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.UNIPROT, value=f('uniprot'))
            ),
            annotations=AnnotationsBuilder(),
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
                + str(v).removeprefix('KEGG.REACTION:'),
            ),
        ),
    ),
    associations=_reaction_associations(),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=model.ChemicalEntity,
            predicate=f(
                'participant_role',
                delimiter='||',
                map='role',
                preserve_indices=True,
            ),
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'participant_chebi', delimiter='||', extract=_CHEBI_RE
                    ),
                ),
                CV(
                    term=Namespace.NAME,
                    value=f('participant_display_name', delimiter='||'),
                ),
            ),
            annotations=AnnotationsBuilder(),
            entity_annotations=AnnotationsBuilder(),
        ),
        MembersFromList(
            entity_type=model.Protein,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.UNIPROT, value=f('uniprot'))
            ),
            annotations=AnnotationsBuilder(),
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
