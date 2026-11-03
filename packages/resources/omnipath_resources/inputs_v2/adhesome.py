"""
Parse the integrin adhesome network and emit Relation records.

Two sources are read.

The literature curated network of Zaidel-Bar et al. 2007 is in the two
supplementary tables of the publication. Table S2 has one interaction per row,
from a source component to a target component, with its kind and one PubMed ID:

    Source  Interaction (1=activating; 2=inhibiting; 3=binding)  target  PMID
    ABL     1                                                    PAXILLIN 9603926

Table S1 describes the components the interactions name: a protein with its
gene name and Swiss-Prot accession, a group of proteins with several
accessions (``Integrin-A``, ``PKC``), or one of four lipids and calcium.

The adhesome.org website publishes the same network together with interactions
taken from other collections (``interactions.csv``):

    Source,Target,Effect,Type,PMID,data_source
    "LPXN","GIT1","0","Binding","12674328","adhesome.org"

The website lists every interaction twice, once in each orientation (its
``process/clean_edges.py`` appends all edges reversed), so ``Source`` and
``Target`` of that file do not say which component acts on which. The rows of
the 2007 network (``data_source`` ``adhesome.org``) are therefore taken from the
publication, which has their direction, and the remaining rows are emitted as
interactions without a direction.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel.model import (
    ChemicalEntity,
    DirectionQualifierEnum,
    Protein,
    ProteinFamily,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import is_pmid, pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.adhesome import (
    iter_interactions,
    iter_website_pairs,
    uniprot_accessions,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
    Relation,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.ADHESOME,
    name='Adhesome',
    url='https://www.adhesome.org/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='17671451',
    topics=(Topic.SIGNALING,),
    content=(Content.CAUSAL_REGULATION, Content.PHYSICAL_INTERACTIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Activating, inhibiting and binding interactions between the '
        'components of integrin adhesions, with publications.'
    ),
    description=(
        'The adhesome is the network of proteins that build and regulate '
        'integrin-mediated cell-matrix adhesions. Zaidel-Bar et al. 2007 '
        'collected its components from the literature, 90 that reside in the '
        'adhesion and 66 that associate with it, and 690 interactions between '
        'them, each classified as binding, activating or inhibiting and '
        'supported by a publication. The adhesome.org website presents the '
        'network with an updated list of 232 components and further '
        'protein-protein interactions between them.'
    ),
)

# Supplementary Tables S2 (interactions) and S1 (components) of Zaidel-Bar et
# al. 2007, from the publisher. The website offered this network with its
# direction until 2016 (``download/fa.sig``); that file is no longer there.
interactions_download = Download(
    url='https://static-content.springer.com/esm/art%3A10.1038%2Fncb0807-858/MediaObjects/41556_2007_BFncb0807858_MOESM3_ESM.xls',
    filename='zaidelbar2007_table_s2_interactions.xls',
    subfolder='adhesome',
    default_mode='rb',
    encoding=None,
)
components_download = Download(
    url='https://static-content.springer.com/esm/art%3A10.1038%2Fncb0807-858/MediaObjects/41556_2007_BFncb0807858_MOESM2_ESM.xls',
    filename='zaidelbar2007_table_s1_components.xls',
    subfolder='adhesome',
    default_mode='rb',
    encoding=None,
)
website_download = Download(
    url='https://www.adhesome.org/interactions.csv',
    filename='adhesome_org_interactions.csv',
    subfolder='adhesome',
)

# The accessions of Table S1 are human Swiss-Prot entries and the website's
# symbols are human gene symbols.
TAXON = 'NCBITaxon:9606'
# The components of Table S1 that are not proteins: four lipids and calcium
# ions. They have no gene name and no accession; ``Synonyms`` holds their
# full name.
CHEMICALS = frozenset({'PIP2', 'PIP3', 'IP3', 'DAG', '[Ca+]'})
# Interaction code -> predicate and direction of the effect on the target, by
# the key in the column name. Activation and inhibition are effects of the
# source on the target; the table does not say what is changed (activity,
# abundance), so no aspect is set. Binding has no direction.
EFFECTS = {
    '1': (slots.affects, DirectionQualifierEnum.increased),
    '2': (slots.affects, DirectionQualifierEnum.decreased),
    '3': (slots.interacts_with, None),
}

# ``Type`` of the website's other rows -> predicate. As the file does not say
# which component is the source, an enzymatic type does not tell the enzyme
# from its substrate, and a signed ``Effect`` does not tell the regulator from
# its target: these rows state only that the two components interact.
WEBSITE_TYPES = {
    'Binding': slots.interacts_with,
    'direct_interaction': slots.interacts_with,
    'Phosphorylation': slots.interacts_with,
    'phosphorylation': slots.interacts_with,
    'phosphorylation_reaction': slots.interacts_with,
    'Dephosphorylation': slots.interacts_with,
    'dephosphorylation': slots.interacts_with,
    'dephosphorylation_reaction': slots.interacts_with,
    'enzymatic_reaction': slots.interacts_with,
    'GEF': slots.interacts_with,
    'inhibition': slots.interacts_with,
    # An effect through an intermediate compound.
    'compound': slots.interacts_with,
    'colocalization': slots.colocalizes_with,
    # Not emitted: the kind is not known and the reference of each such row is
    # the placeholder below.
    'Unknown': None,
}
# Written in place of a PubMed ID by the three ``Unknown`` interactions.
PLACEHOLDER_PMID = '12345678'
# The website's file was last changed in October 2022, when PubMed IDs were
# below this. A few of its rows have larger numbers, which are not PubMed IDs.
PMID_LIMIT = 37_000_000


def _protein(accession: str, row: dict[str, Any] | None = None, role: str = '') -> Entity:
    identifiers = [Identifier(type=Namespace.UNIPROT, value=accession)]
    if row:
        # The gene name is the symbol of the time (``C20orf42``, ``DDEF1``).
        for namespace, field in (
            (Namespace.GENESYMBOL, 'gene'),
            (Namespace.NAME, 'name'),
        ):
            if value := row.get(f'{role}_{field}'):
                identifiers.append(Identifier(type=namespace, value=value))
    return Entity(
        type=Protein,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


def _component(row: dict[str, Any], role: str) -> Entity | None:
    """The component as a protein, a group of proteins or a chemical."""
    name = row.get(f'{role}_name')
    if not name:
        return None
    if name in CHEMICALS:
        identifiers = [Identifier(type=Namespace.NAME, value=name)]
        if full_name := row.get(f'{role}_synonyms'):
            identifiers.append(Identifier(type=Namespace.SYNONYM, value=full_name))
        return Entity(type=ChemicalEntity, identifiers=identifiers)
    accessions = uniprot_accessions(row.get(f'{role}_uniprot'))
    if not accessions:
        return None
    if len(accessions) == 1:
        return _protein(accessions[0], row, role)
    # Several accessions: the component stands for any of these proteins. The
    # gene names of such a row are not aligned with its accessions.
    return Entity(
        type=ProteinFamily,
        identifiers=[Identifier(type=Namespace.NAME, value=name)],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
        membership=[
            Membership(member=_protein(accession)) for accession in accessions
        ],
    )


def map_interaction(row: dict[str, Any]) -> Relation | None:
    """An interaction of the 2007 network, from its source to its target."""
    code = row.get('interaction')
    if code not in EFFECTS:
        raise ValueError(f'Unknown Adhesome interaction code: {code!r}')
    predicate, direction = EFFECTS[code]
    source = _component(row, 'source')
    target = _component(row, 'target')
    if source is None or target is None:
        return None
    annotations = []
    if direction is not None:
        annotations.append(
            Annotation(term=slots.object_direction_qualifier, value=direction)
        )
    if curie := pmid_curie(row.get('pmid')):
        annotations.append(Annotation(term=slots.publications, value=curie))
    return Relation(
        subject=source,
        predicate=predicate,
        object=target,
        annotations=annotations,
    )


def _website_pubmed(value: str) -> str | None:
    if (
        not is_pmid(value)
        or int(value) >= PMID_LIMIT
        or value == PLACEHOLDER_PMID
    ):
        return None
    return pmid_curie(value)


f = FieldConfig(transform={'pmid_curie': _website_pubmed})


def _website_protein(column: str) -> EntityBuilder:
    # The website names its components by official gene symbol.
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(column)),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=TAXON),
        ),
    )


def _website_predicate(row: dict[str, Any]) -> Any:
    """The predicate of a website interaction, None if it is not emitted."""
    types = row.get('types') or []
    unknown = [value for value in types if value not in WEBSITE_TYPES]
    if unknown or not types:
        raise ValueError(f'Unknown Adhesome interaction type: {unknown or types!r}')
    predicates = [WEBSITE_TYPES[value] for value in types]
    # Mirrored rows that differ: an interaction includes being in one place.
    if any(predicate is slots.interacts_with for predicate in predicates):
        return slots.interacts_with
    if any(predicate is slots.colocalizes_with for predicate in predicates):
        return slots.colocalizes_with
    return None


# The components are in alphabetical order, not source and target.
website_interactions_schema = RelationBuilder(
    subject=_website_protein('component_a'),
    predicate=_website_predicate,
    object=_website_protein('component_b'),
    annotations=AnnotationsBuilder(
        CV(term=slots.publications, value=f('PMID', transform='pmid_curie')),
    ),
)


def map_website_interaction(row: dict[str, Any]) -> Relation | None:
    """An interaction the website added, without a direction."""
    if _website_predicate(row) is None:
        return None
    return website_interactions_schema(row)

resource = Resource(
    config,
    interactions=Dataset(
        download=interactions_download,
        mapper=map_interaction,
        raw_parser=functools.partial(
            iter_interactions, components_download=components_download
        ),
    ),
    website_interactions=Dataset(
        download=website_download,
        mapper=map_website_interaction,
        raw_parser=iter_website_pairs,
    ),
)
