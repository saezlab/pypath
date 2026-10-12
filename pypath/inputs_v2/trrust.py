"""
Parse the TRRUST regulatory relationships and emit Relation records.

Source files (tab separated, no header; one file per species):

    trrust_rawdata.human.tsv
    trrust_rawdata.mouse.tsv

    AATF    BAX     Repression  22909821
    ABL1    BAX     Activation  11753601
    ABL1    BCL2    Repression  11753601
    AATF    CDKN1A  Unknown     17157788

The columns are the transcription factor, the target gene (both official gene
symbols of the file's species), the mode of regulation and the PubMed IDs of
the articles the relationship was curated from, joined by ``;``.

A pair has one line per mode of regulation the curators found for it, each
with its own articles, so the same pair can be listed as activated, repressed
and of unknown mode. Every line is one statement here; the conflict is the
source's and is not resolved.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.trrust import iter_relationships
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
    RelationBuilder,
)

BASE_URL = 'https://www.grnpedia.org/trrust/data/'
# The files have no organism column: the species is the file.
SPECIES = {'human': '9606', 'mouse': '10090'}

_DIRECTION = model.DirectionQualifierEnum
_ASPECT = model.GeneOrGeneProductOrChemicalEntityAspectEnum
# Mode of regulation -> effect on the target gene's expression. "Unknown" is
# the curators' value for a relationship whose articles do not say which way.
_MODE_OF_REGULATION = {
    'Activation': _DIRECTION.increased,
    'Repression': _DIRECTION.decreased,
    'Unknown': None,
}


def _direction(row: dict[str, Any]):
    """The curated mode of regulation, None where the source has "Unknown"."""
    mode = row.get('mode')
    if mode not in _MODE_OF_REGULATION:
        raise ValueError(f'Unknown TRRUST mode of regulation: {mode!r}')
    return _MODE_OF_REGULATION[mode]


def _aspect(row: dict[str, Any]):
    return _ASPECT.expression if _direction(row) is not None else None


config = ResourceConfig(
    id=ResourceCv.TRRUST,
    name='TRRUST',
    full=(
        'Transcriptional Regulatory Relationships Unraveled by Sentence-based '
        'Text mining'
    ),
    url='https://www.grnpedia.org/trrust/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_SA_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='29087512',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.CAUSAL_REGULATION,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Human and mouse transcription factor–target gene regulation with '
        'literature references and the mode of regulation where the articles '
        'state it.'
    ),
    description=(
        'TRRUST (version 2) is a manually curated database of human and mouse '
        'transcriptional regulatory relationships. Sentences that may '
        'describe the regulation of a gene by a transcription factor were '
        'found in PubMed abstracts by text mining and then read by curators, '
        'who kept the relationships reported by small-scale experimental '
        'studies. Each relationship has the PubMed IDs of its articles and a '
        'mode of regulation: activation, repression or unknown. A pair is '
        'listed once for each mode the literature reports.'
    ),
)

f = FieldConfig(
    extract={
        # One token of the human file, "1724219e", is not a PubMed ID; it is
        # left out, not repaired.
        'pmid': PMID_EXTRACT,
    },
    transform={
        'pmid_curie': lambda v: f'PMID:{v}',
        'taxon_curie': lambda v: f'NCBITaxon:{v}',
    },
)


def _participant(entity_type, column: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=entity_type,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(column)),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.in_taxon,
                value=f('taxon_id', transform='taxon_curie'),
            ),
        ),
    )


# A row without a transcription factor or a target yields no relation: the
# participant has no identifier, so the builder returns None.
interactions_schema = RelationBuilder(
    subject=_participant(model.Protein, 'tf'),
    predicate=slots.affects,
    object=_participant(model.Gene, 'target'),
    annotations=AnnotationsBuilder(
        # Activation or repression of the target gene's expression, as curated
        # from the row's articles. "Unknown" gets neither qualifier.
        CV(term=slots.object_direction_qualifier, value=_direction),
        CV(term=slots.object_aspect_qualifier, value=_aspect),
        CV(
            term=slots.publications,
            value=f(
                'pmids',
                delimiter=';',
                extract='pmid',
                transform='pmid_curie',
            ),
        ),
    ),
)


def _dataset(species: str) -> Dataset:
    return Dataset(
        download=Download(
            url=f'{BASE_URL}trrust_rawdata.{species}.tsv',
            filename=f'trrust_rawdata.{species}.tsv',
            subfolder='trrust',
        ),
        mapper=interactions_schema,
        raw_parser=functools.partial(
            iter_relationships,
            taxon_id=SPECIES[species],
        ),
    )


resource = Resource(
    config,
    human_interactions=_dataset('human'),
    mouse_interactions=_dataset('mouse'),
)
