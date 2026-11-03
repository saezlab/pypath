"""
Parse the STRING protein networks and emit Relation records.

STRING publishes three networks for each organism, each as a table of protein
pairs with a combined score and the scores of the evidence channels behind
it. They assert different things and are separate datasets:

- the functional network: the two proteins contribute to a shared function.
  Its evidence includes coexpression, genomic context and co-mention in the
  literature, so a link is an association, not an interaction;
- the physical subnetwork: the two proteins bind or are part of the same
  complex;
- the regulatory network (since version 12.5): the first protein regulates
  the second, with separate scores for positive and negative regulation.

Scores are probabilities multiplied by 1000. Each dataset names the
confidence level it reads; nothing is filtered silently.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel.model import (
    DirectionQualifierEnum,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2._measurements import measurement as _measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.string import iter_links
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.STRING,
    name='STRING',
    url='https://string-db.org/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='39558183',
    content=(Content.PHYSICAL_INTERACTIONS, Content.CAUSAL_REGULATION),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
        KnowledgeOrigin.TEXT_MINED,
    ),
    summary=(
        'High-confidence human protein links of the STRING functional, '
        'physical and regulatory networks, with the combined score and the '
        'score of each evidence channel.'
    ),
    description=(
        'STRING is a database of known and predicted protein-protein '
        'associations. It integrates experimental interaction data, curated '
        'pathway and complex databases, automated text mining, coexpression '
        'and genomic context, transfers evidence between organisms by '
        'homology, and summarises the evidence for each protein pair in a '
        'combined confidence score. Besides the functional association '
        'network it provides a physical subnetwork of proteins that bind or '
        'share a complex and, since version 12.5, a network of directed '
        'regulatory relationships.'
    ),
)

# Each version has its own files; this is the version the module reads.
VERSION = '12.5'
_URL = (
    'https://stringdb-downloads.org/download/'
    'protein.{network}links.detailed.v{version}/'
    '{taxon}.protein.{network}links.detailed.v{version}.txt.gz'
)
# STRING's names for its combined score cutoffs. The functional and physical
# files contain no pair below 150.
CONFIDENCE = {'low': 150, 'medium': 400, 'high': 700, 'highest': 900}
# Evidence channels of each network, as the ``detailed`` files name them.
FUNCTIONAL_CHANNELS = (
    'neighborhood',
    'fusion',
    'cooccurence',
    'coexpression',
    'experimental',
    'database',
    'textmining',
)
PHYSICAL_CHANNELS = ('experimental', 'database', 'textmining')
REGULATORY_CHANNELS = ('experimental', 'database', 'textmining')
# Scores of the two signs of a regulatory link.
REGULATORY_SIGNS = ('positive', 'negative')


def regulatory_direction(row: dict[str, Any]) -> DirectionQualifierEnum | None:
    """The sign of a regulatory link that has a score for one sign only.

    About one link in nine has a score for both signs and one in four for
    neither; these get no direction of effect.
    """
    positive = int(row.get('positive') or 0) > 0
    negative = int(row.get('negative') or 0) > 0
    if positive == negative:
        return None
    return (
        DirectionQualifierEnum.increased
        if positive
        else DirectionQualifierEnum.decreased
    )


f = FieldConfig(
    extract={
        # A STRING identifier is the NCBI taxon and the protein's identifier
        # in its source genome, joined by a dot: ``9606.ENSP00000000233``.
        'string': r'^(\d+\.\S+)$',
        'taxon': r'^(\d+)\.\S+$',
        'ensp': r'^\d+\.(ENS[A-Z]*P\d{11})$',
    },
    transform={'taxon': lambda v: f'NCBITaxon:{v}'},
)


def _protein(column: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.STRING, value=f(column, extract='string')),
            CV(term=Namespace.ENSP, value=f(column, extract='ensp')),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.in_taxon,
                value=f(column, extract='taxon', transform='taxon'),
            ),
        ),
    )


def _score(column: str) -> CV:
    """A score on the source's scale of 0 to 1000; 0 is a channel without evidence."""
    return CV(
        term=slots.has_confidence_score,
        value=f(
            column,
            transform=lambda v: (
                None if v == '0' else _measurement(v, source_field=column)
            ),
        ),
    )


def _schema(predicate, channels: tuple[str, ...], *qualifiers: CV) -> RelationBuilder:
    return RelationBuilder(
        subject=_protein('protein1'),
        predicate=predicate,
        object=_protein('protein2'),
        annotations=AnnotationsBuilder(
            *qualifiers,
            _score('combined_score'),
            *map(_score, channels),
        ),
    )


# A functional link is an association: its evidence may be coexpression,
# genomic context or co-mention alone.
functional_schema = _schema(slots.associated_with, FUNCTIONAL_CHANNELS)
# Proteins that bind or share a complex; their contact may be indirect.
physical_schema = _schema(slots.interacts_with, PHYSICAL_CHANNELS)
# The first protein regulates the second. The files do not say what is
# regulated (activity, abundance, a modification), so no aspect is given.
regulatory_schema = _schema(
    slots.affects,
    (*REGULATORY_CHANNELS, *REGULATORY_SIGNS),
    CV(term=slots.object_direction_qualifier, value=regulatory_direction),
)


def _dataset(
    network: str,
    schema: RelationBuilder,
    columns: tuple[str, ...],
    confidence: str,
    taxon: int = 9606,
) -> Dataset:
    filename = _URL.format(
        network=network, version=VERSION, taxon=taxon
    ).rsplit('/', 1)[1]
    return Dataset(
        download=Download(
            url=_URL.format(network=network, version=VERSION, taxon=taxon),
            filename=filename,
            subfolder='string',
            large=True,
            ext='gz',
        ),
        mapper=schema,
        raw_parser=functools.partial(
            iter_links,
            columns=columns,
            minimum_score=CONFIDENCE[confidence],
            directed=schema is regulatory_schema,
        ),
    )


resource = Resource(
    config,
    human_physical_high_confidence=_dataset(
        'physical.', physical_schema, PHYSICAL_CHANNELS, 'high'
    ),
    human_functional_high_confidence=_dataset(
        '', functional_schema, FUNCTIONAL_CHANNELS, 'high'
    ),
    human_regulatory_high_confidence=_dataset(
        'regulatory.',
        regulatory_schema,
        (*REGULATORY_CHANNELS, *REGULATORY_SIGNS),
        'high',
    ),
)
