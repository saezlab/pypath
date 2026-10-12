"""
Parse the curated interactions of the NRF2 regulatory network and emit Relation records.

Source file (comma separated, no header; one interaction per line, with a
comma ending every line):

    Q14145,Q16236,1,1,-1,9887101,

The columns are not named in the file. Their order is that of the legacy
reader, and is that of the curated file of the Autophagy Regulatory Network
(``arn.py``), made by the same group with the same export, minus the two
autophagy flags:

    source, target     UniProt accessions of the two proteins
    is_direct          1 for a direct interaction, 0 for an indirect one
    is_directed        1 for a directed interaction, 0 for an undirected one
    effect             1 stimulation, -1 inhibition, 0 not known
    references         PubMed IDs joined by colons

The meaning of the columns is verified for ARN against AutophagyNet, not for
this file: NRF2ome offers no named table any more.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel.model import DirectionQualifierEnum, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import PMID_EXTRACT, pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_csv
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.NRF2OME,
    name='NRF2ome',
    full='NRF2 regulatory network',
    url='https://nrf2ome.org/',
    # The site states no license.
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_SA_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='23710289',
    topics=(Topic.SIGNALING,),
    content=(Content.CAUSAL_REGULATION, Content.PHYSICAL_INTERACTIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Curated inhibiting and undirected interactions of the NRF2 '
        'pathway proteins, with PubMed references.'
    ),
    description=(
        'NRF2-ome is a web resource on the transcription factor NRF2 '
        '(NFE2L2), the master regulator of the response to oxidative stress, '
        'and on its interactions and regulators: protein-protein and '
        'signaling interactions curated from the literature, and '
        'transcriptional and post-transcriptional regulators, with '
        'references. It was made by the Korcsmaros and Csermely groups. '
        'This module reads the curated protein interactions.'
    ),
)

# The legacy address http://nrf2.elte.hu/ does not answer. nrf2ome.org still
# serves the site's front page, but no download. This is OmniPath's copy of the
# curated interactions the legacy reader used.
download = Download(
    url='https://rescued.omnipathdb.org/nrf2ome.csv',
    filename='nrf2ome.csv',
    subfolder='nrf2ome',
)
COLUMNS = (
    'source',
    'target',
    'is_direct',
    'is_directed',
    'effect',
    'references',
    # Every line ends with a comma.
    'trailing',
)

# The file has no organism column. The UniProt entries of all its accessions
# are human proteins (checked against UniProt, not stated by the file).
TAXON = 'NCBITaxon:9606'
_UP = DirectionQualifierEnum.increased
_DOWN = DirectionQualifierEnum.decreased
# ``is_directed`` -> whether the interaction has a direction.
DIRECTED = {'1': True, '0': False}
# ``effect`` -> direction of the effect of the source on the target; the file
# does not say what is changed, so no aspect is set (as in ``ca1.py``). 0 is
# "not known". The file has no stimulation (1), but the column is the one of
# ARN, which has.
EFFECTS = {'1': _UP, '-1': _DOWN, '0': None}


def _statement(row: dict[str, Any]) -> tuple[Any, Any]:
    """Predicate and direction. An undirected interaction has no source or
    target and no direction of effect."""
    directed, effect = row.get('is_directed'), row.get('effect')
    if directed not in DIRECTED:
        raise ValueError(f'Unknown NRF2ome is_directed value: {directed!r}')
    if effect not in EFFECTS:
        raise ValueError(f'Unknown NRF2ome effect: {effect!r}')
    if not DIRECTED[directed]:
        return slots.interacts_with, None
    return slots.affects, EFFECTS[effect]


f = FieldConfig()


def _protein(column: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=f(column))),
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=TAXON)),
    )


interactions_schema = RelationBuilder(
    subject=_protein('source'),
    predicate=lambda row: _statement(row)[0],
    object=_protein('target'),
    annotations=AnnotationsBuilder(
        CV(term=slots.object_direction_qualifier, value=lambda row: _statement(row)[1]),
        CV(
            term=slots.publications,
            value=f('references', delimiter=':', extract=PMID_EXTRACT, transform=pmid_curie),
        ),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=interactions_schema,
        raw_parser=functools.partial(
            iter_csv,
            fieldnames=COLUMNS,
            strict=True,
            quotechar=None,
            required=COLUMNS,
        ),
    ),
)
