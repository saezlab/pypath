"""
Parse SignaLink 3 and emit pathway annotations and signed, directed interactions.

SignaLink is a signaling pathway resource of the Korcsmaros group. Its pathway layer
(layer 0) lists the proteins of the signaling pathways with the pathway and the
function (ligand, receptor, mediator, scaffold, co-factor, transcription factor,
endocytosis related) of each, and the directed interactions between them, merged
from SignaLink's own curation from pathway reviews and from ACSN, InnateDB,
Reactome and SIGNOR.

The module has three datasets:

- ``terms``: the pathways and functions of the nodes file, as terms of SignaLink;
- ``annotations``: each protein that has a pathway or a function, associated with
  the term of each;
- ``interactions``: the directed pairs, from the first node to the second.

A pair with an effect term of PSI-MI (stimulation, inhibition, up- and
down-regulation, the enzymatic reactions) is an ``affects`` relation with the
direction and aspect of the effect; a pair with terms of both directions is one
relation for each. A pair that says only that the nodes bind or associate is an
interaction, a direct one where SignaLink says it is direct. The layers of
transcriptional, post-transcriptional and lncRNA regulators of the website are not
in the files read.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2 import _mitab
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.signalink import (
    iter_interactions,
    iter_nodes,
    iter_terms,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Association,
    Entity,
    EntityRef,
    Identifier,
    Relation,
)

config = ResourceConfig(
    id=ResourceCv.SIGNALINK,
    name='SignaLink',
    full='SignaLink 3',
    url='https://signalink.org/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_SA_3_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='34634810',
    topics=(Topic.SIGNALING,),
    content=(
        Content.CAUSAL_REGULATION,
        Content.PHYSICAL_INTERACTIONS,
        Content.MOLECULAR_ANNOTATIONS,
    ),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Signaling pathway members with their pathways and functions, and '
        'directed, signed interactions between them with publications.'
    ),
    description=(
        'SignaLink is an integrated resource of signaling pathways of human, '
        'fruit fly, nematode and zebrafish. The proteins of the major pathways '
        'are assigned from the full texts of pathway reviews, classified by '
        'their function in the pathway, and joined by directed interactions, '
        'signed where the literature gives the effect, which come from the '
        'curation of SignaLink and from ACSN, InnateDB, Reactome and SIGNOR.'
    ),
)

PATHWAY = 'signalink_pathway'
FUNCTION = 'signalink_function'
NAMESPACES = {'pathway': PATHWAY, 'function': FUNCTION}
_ANNOTATION_NAMESPACES = {'pathways': PATHWAY, 'functions': FUNCTION}

IS_DIRECT = 'signalink:is_direct'

# The databases SignaLink merges, as ``source_db`` labels them. SLKv2.0, SLKv2.1 and
# SLKv3.0 are the releases of the curation of SignaLink itself, SignaFish its
# zebrafish resource and TCRcuration the curation of T cell receptor signaling.
SOURCES = {
    label: label
    for label in (
        'ACSN',
        'InnateDB',
        'Reactome',
        'SLKv2.0',
        'SLKv2.1',
        'SLKv3.0',
        'SignaFish',
        'Signor',
        'TCRcuration',
    )
}
_DIRECTIONS = {
    'increased': model.DirectionQualifierEnum.increased,
    'decreased': model.DirectionQualifierEnum.decreased,
}
_ASPECTS = {
    aspect: getattr(model.GeneOrGeneProductOrChemicalEntityAspectEnum, key)
    for aspect, key in {
        'activity': 'activity',
        'expression': 'expression',
        'stability': 'stability',
        'abundance': 'abundance',
        'phosphorylation': 'phosphorylation',
        'ubiquitination': 'ubiquitination',
        'acetylation': 'acetylation',
        'methylation': 'methylation',
        'sumoylation': 'sumoylation',
        'ADP-ribosylation': 'ADP-ribosylation',
        'cleavage': 'cleavage',
        'localization': 'localization',
        'guanyl_nucleotide_exchange': 'guanyl_nucleotide_exchange',
    }.items()
}


def map_term(row: dict[str, Any]) -> Entity:
    """A pathway or a function of the nodes file."""
    return Entity(
        type=model.OntologyClass,
        identifiers=[
            Identifier(type=NAMESPACES[row['kind']], value=row['term']),
            Identifier(type=Namespace.NAME, value=row['term']),
        ],
    )


def _node(node: dict[str, Any], annotated: bool = False) -> Entity:
    """The protein or RNA of a node, with the associations of its annotations."""
    if node['database'] == 'Uniprot':
        entity_type, namespace = model.Protein, Namespace.UNIPROT
    elif node['database'] == 'RNACentral':
        entity_type, namespace = model.RNAProduct, Namespace.RNACENTRAL
    else:
        raise ValueError(f'Unknown SignaLink node database: {node["database"]!r}')
    associations = None
    if annotated:
        associations = [
            Association(
                object=EntityRef(model.OntologyClass, _ANNOTATION_NAMESPACES[field], term),
                predicate=slots.associated_with,
            )
            for field in ('pathways', 'functions')
            for term in node[field]
        ]
    return Entity(
        type=entity_type,
        identifiers=[Identifier(type=namespace, value=node['accession'])],
        annotations=[
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{node["taxon"]}')
        ],
        associations=associations,
    )


def map_annotation(row: dict[str, Any]) -> Entity:
    """A node associated with its pathways and functions."""
    return _node(row, annotated=True)


def _predicate(row: dict[str, Any]):
    if row['affects']:
        return slots.affects
    predicate = _mitab.physical_predicate(row['interaction_types'])
    if predicate is None and row['interaction_types']:
        # Only "interaction type", "bind" or "biochemical": the pair interacts.
        predicate = slots.interacts_with
    if predicate is None:
        # The file gives the pair no type: it is joined in the pathway, and
        # directly where the row says so.
        predicate = slots.associated_with
        if row['is_direct']:
            return slots.directly_physically_interacts_with
        return predicate
    if row['is_direct'] and predicate == slots.interacts_with:
        return slots.directly_physically_interacts_with
    return predicate


def map_interaction(row: dict[str, Any]) -> Relation:
    """A directed pair of nodes, with the effect of the first on the second."""
    annotations = [
        Annotation(term=IS_DIRECT, value='true' if row['is_direct'] else 'false')
    ]
    if row['direction']:
        annotations.append(
            Annotation(
                term=slots.object_direction_qualifier,
                value=_DIRECTIONS[row['direction']],
            )
        )
    if row['aspect']:
        annotations.append(
            Annotation(
                term=slots.object_aspect_qualifier, value=_ASPECTS[row['aspect']]
            )
        )
    annotations += [
        Annotation(term=slots.has_evidence_of_type, value=term)
        for term in row['methods']
    ]
    annotations += [
        Annotation(term=slots.original_predicate, value=term)
        for term in row['interaction_types']
    ]
    annotations += [
        Annotation(term=slots.publications, value=pmid) for pmid in row['pmids']
    ]
    for source in row['sources']:
        if source not in SOURCES:
            raise ValueError(f'Unknown SignaLink source database: {source!r}')
        annotations.append(
            Annotation(term=slots.supporting_data_source, value=SOURCES[source])
        )
    return Relation(
        subject=_node(row['a']),
        predicate=_predicate(row),
        object=_node(row['b']),
        annotations=annotations,
    )


# The files of the 2020-03-23 release of SignaLink 3 that OmniPath keeps. The
# website (signalink.org/download) now offers v3.1 only through a form that selects
# species, pathways and layers, and the full SignaLink 2.0 as one archive.
nodes_download = Download(
    url='https://rescued.omnipathdb.org/SLK3_nodes_2020_03_23.csv',
    filename='SLK3_nodes_2020_03_23.csv',
    subfolder='signalink',
    large=True,
)
edges_download = Download(
    url='https://rescued.omnipathdb.org/SLK3_layer0_2020_03_23.csv',
    filename='SLK3_layer0_2020_03_23.csv',
    subfolder='signalink',
    large=True,
)

resource = Resource(
    config,
    terms=Dataset(download=nodes_download, mapper=map_term, raw_parser=iter_terms),
    annotations=Dataset(
        download=nodes_download, mapper=map_annotation, raw_parser=iter_nodes
    ),
    interactions=Dataset(
        download=edges_download,
        mapper=map_interaction,
        raw_parser=functools.partial(iter_interactions, nodes_download=nodes_download),
    ),
)
