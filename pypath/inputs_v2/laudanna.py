"""
Parse the signaling network edge attributes of the Laudanna Lab and emit Relation records.

The Laboratory of Cell Trafficking and Signal Transduction (University of
Verona) offers data sets for network analysis in Cytoscape formats. Two of its
edge attribute files describe a directed human protein signaling network, the
edges named by HGNC gene symbols:

    Signaling_flow
    F2R (pp) PAWR = state_change

    SignalingDirection
    EDNRA (pp) NBN = activation

- ``SignalingFlow.EA`` ("directed PPI signaling networks with edge attributes
  of type state_change"): the interaction type of each edge. The lab's page
  names the Pathway Commons binary interaction rules as the meaning of these
  attributes: STATE_CHANGE, "the first entity controls a reaction that changes
  the state of the second ... This interaction is directed". The file also
  has ``complex_formation`` edges, which the page does not define.
- ``SignalingDirection.EA`` ("information of type activation-inhibition-docking
  activity"): ``activation``, ``inhibition`` or ``docking`` for each edge.

The lab's page says the data sets "are compiled from public data-bases and from
literature and manually curated"; it does not say how these two files were
made. The flow file is made of the Pathway Commons rules, which infer binary
interactions from pathway records, and the legacy description calls both
"directionality and effect information inferred from multiple sources". The
files are therefore tagged as inferred, not as curated.

The first node is the subject and the second the object. The two files give
separate statements and are separate datasets; an edge can have several lines,
and each line is one relation.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel.model import DirectionQualifierEnum, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.laudanna import (
    DIRECTION_VALUES,
    FLOW_VALUES,
    SIGNALING_DIRECTION,
    SIGNALING_FLOW,
    iter_edges,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.LAUDANNA,
    name='Laudanna',
    full='Compiled Datasets for Network Analysis from Laudanna Lab',
    url='http://dp.univr.it/~laudanna/LCTST/downloads/',
    # The page states no license.
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    topics=(Topic.SIGNALING,),
    content=(Content.CAUSAL_REGULATION, Content.PHYSICAL_INTERACTIONS),
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary=(
        'Inferred directed links between human signaling proteins: the type '
        'of each link, and its activating, inhibiting or docking effect.'
    ),
    description=(
        'The Laboratory of Cell Trafficking and Signal Transduction at the '
        'University of Verona offers data sets for network analysis in '
        'Cytoscape formats, compiled from public databases and the '
        'literature. Two edge attribute files describe a directed human '
        'signaling network of about 77,500 and 63,400 links between HGNC '
        'gene symbols: the type of each link (state change or complex '
        'formation, after the Pathway Commons binary interaction rules) and '
        'its effect (activation, inhibition or docking). The lab does not '
        'describe how the two files were derived.'
    ),
)

# The lab's pages went through several moves; the legacy addresses are
# http://dp.univr.it/~laudanna/LCTST/downloads/files/SignalingFlow.EA and
# .../SignalingDirection.EA, which still answer. This is OmniPath's copy of the
# 2014 files: the flow file is byte for byte the original, and the direction
# file is the original with its old Mac line endings changed to newlines.
download_flow = Download(
    url='https://rescued.omnipathdb.org/Laudanna_SignalingFlow.EA',
    filename='Laudanna_SignalingFlow.EA',
    subfolder='laudanna',
    large=True,
)
download_direction = Download(
    url='https://rescued.omnipathdb.org/Laudanna_SignalingDirection.EA',
    filename='Laudanna_SignalingDirection.EA',
    subfolder='laudanna',
    large=True,
)

# The lab's page: "Taxid is 9606 Homo sapiens. Gene/protein ID is provided by
# HGNC."
TAXON = 'NCBITaxon:9606'
_UP = DirectionQualifierEnum.increased
_DOWN = DirectionQualifierEnum.decreased
# ``SignalingFlow`` value -> predicate and direction of the effect on the target.
# STATE_CHANGE is defined by the Pathway Commons rule on the lab's page: the
# first entity controls a change of state of the second (a modification or a
# change of location), directed. What changes and in which direction is not
# stated, so neither an aspect nor a direction is set. The page does not define
# ``complex_formation``; the name says that the two form a complex, which is an
# interaction without a direction.
FLOW = {
    'state_change': (slots.affects, None),
    'complex_formation': (slots.interacts_with, None),
}
# ``SignalingDirection`` value -> predicate and direction of the effect on the
# target. Activation and inhibition are effects of the first protein on the
# second; the file does not say what is changed, so no aspect is set (as in
# ``ca1.py``). The page gives no definition of ``docking``; it is taken as
# binding, an interaction without an effect or a direction.
DIRECTION = {
    'activation': (slots.affects, _UP),
    'inhibition': (slots.affects, _DOWN),
    'docking': (slots.interacts_with, None),
}


def _crosswalk(table: dict[str, Any], name: str):
    def lookup(row: dict[str, Any]) -> Any:
        value = row.get('value')
        if value not in table:
            raise ValueError(f'Unknown Laudanna {name} value: {value!r}')
        return table[value]

    return lookup


def _participant(column: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=lambda row: row.get(column)),
        ),
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=TAXON)),
    )


def _schema(table: dict[str, Any], name: str) -> RelationBuilder:
    effect = _crosswalk(table, name)
    # A line without a second node yields no relation: the participant has no
    # identifier, so the builder returns None.
    return RelationBuilder(
        subject=_participant('source'),
        predicate=lambda row: effect(row)[0],
        object=_participant('target'),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.object_direction_qualifier,
                value=lambda row: effect(row)[1],
            ),
        ),
    )


resource = Resource(
    config,
    signaling_flow=Dataset(
        download=download_flow,
        mapper=_schema(FLOW, 'signaling flow'),
        raw_parser=functools.partial(
            iter_edges, title=SIGNALING_FLOW, values=FLOW_VALUES
        ),
    ),
    signaling_direction=Dataset(
        download=download_direction,
        mapper=_schema(DIRECTION, 'signaling direction'),
        raw_parser=functools.partial(
            iter_edges, title=SIGNALING_DIRECTION, values=DIRECTION_VALUES
        ),
    ),
)
