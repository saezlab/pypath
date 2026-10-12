"""
Parse the curated interactions of the Autophagy Regulatory Network and emit Relation records.

Source file (comma separated, no header; one interaction per line, with a
comma ending every line):

    O75385,O75143,1,1,1,1,1,20562859:19225151:19287211:19211835,

The columns are not named in the file. Their order is that of the legacy
reader and the names are checked against AutophagyNet, the successor of the
site (its table of the whole ARN 1.0 database lists the same interactions
with named columns):

    source, target     UniProt accessions of the two proteins
    is_direct          1 for a direct interaction, 0 for an indirect one
                       (the AutophagyNet column ``directness``)
    is_directed        1 for a directed interaction, 0 for an undirected one
                       (``PPI directed`` or ``PPI undirected`` of the column
                       ``interaction_type``)
    effect             1 stimulation, -1 inhibition, 0 not known
                       (``effect`` of AutophagyNet)
    source_autophagy,  1 if the protein is a core autophagy protein
    target_autophagy
    references         PubMed IDs joined by colons

Of the 91 pairs that AutophagyNet also lists, all agree with this reading.
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
    id=ResourceCv.ARN,
    name='ARN',
    full='Autophagy Regulatory Network',
    url='https://www.autophagynet.org/',
    # Neither the site nor its download page states a license.
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_SA_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='25635527',
    topics=(Topic.SIGNALING,),
    content=(Content.CAUSAL_REGULATION, Content.PHYSICAL_INTERACTIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Curated stimulating, inhibiting and undirected interactions of '
        'human autophagy proteins, with PubMed references.'
    ),
    description=(
        'The Autophagy Regulatory Network (ARN) is a systems-level resource '
        'on the components of autophagy and their regulation, made by the '
        'Korcsmaros, Vellai and Lenti groups. Its core is a manually curated '
        'set of interactions between autophagy proteins and their '
        'regulators, with their effect and PubMed references; it was '
        'extended with pathway, post-translational, transcriptional and '
        'post-transcriptional regulators. The site is now AutophagyNet. This '
        'module reads the curated interactions of ARN 1.0.'
    ),
)

# The legacy address http://autophagyregulation.org/ now leads to AutophagyNet,
# which offers the whole ARN 1.0 database (190 MB, with its tables of
# regulators) but not this file. This is OmniPath's copy of the curated
# interactions the legacy reader used.
download = Download(
    url='https://rescued.omnipathdb.org/arn_curated.csv',
    filename='arn_curated.csv',
    subfolder='arn',
)
COLUMNS = (
    'source',
    'target',
    'is_direct',
    'is_directed',
    'effect',
    'source_autophagy',
    'target_autophagy',
    'references',
    # Every line ends with a comma.
    'trailing',
)

# The file has no organism column; the interactions of the curated set are
# human (the 93 lines of AutophagyNet that match them have taxon 9606).
TAXON = 'NCBITaxon:9606'
_UP = DirectionQualifierEnum.increased
_DOWN = DirectionQualifierEnum.decreased
# ``is_directed`` -> whether the interaction has a direction.
DIRECTED = {'1': True, '0': False}
# ``effect`` -> direction of the effect of the source on the target. The
# stimulation or inhibition is the effect on the target protein; the file does
# not say what is changed, so no aspect is set (as in ``ca1.py``). 0 is "not
# known" (``unknown`` in AutophagyNet).
EFFECTS = {'1': _UP, '-1': _DOWN, '0': None}


def _statement(row: dict[str, Any]) -> tuple[Any, Any]:
    """Predicate and direction. An undirected interaction has no source or
    target and no direction of effect, so its sign, which the file gives to a
    few of them, is not carried."""
    directed, effect = row.get('is_directed'), row.get('effect')
    if directed not in DIRECTED:
        raise ValueError(f'Unknown ARN is_directed value: {directed!r}')
    if effect not in EFFECTS:
        raise ValueError(f'Unknown ARN effect: {effect!r}')
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
