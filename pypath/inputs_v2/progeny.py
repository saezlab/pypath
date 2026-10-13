"""
Pathway responsive genes of PROGENy.

PROGENy takes a compendium of signaling perturbation experiments and yields,
for each of 14 pathways, the genes whose expression responds to the pathway's
activity, each with a weight and a p-value. The module imports the full human
model (``model_human_full``, the legacy module's choice) as two datasets:

* ``terms``: the 14 pathways of the model as entities with names and the
  descriptions of the package vignette;
* ``responses``: a statement for each row of the model, a gene associated with
  a pathway, with the weight and the p-value. No row is left out: the legacy
  module's top-N and p-value filters are not applied.

The pathways are PROGENy's own vocabulary and are source terms. The mouse
model of the package (``model_mouse_full``, MGI symbols written as HGNC
symbols in its documentation) is not imported.
"""

from __future__ import annotations

from biolink_model.datamodel.model import Gene, OntologyClass, slots
from omnipath_core.measurements import Measurement
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.progeny import iter_model, iter_pathways
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Relation
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.PROGENY,
    name='PROGENy',
    full='Pathway RespOnsive GENes',
    url='https://saezlab.github.io/progeny/',
    # The repository's LICENSE file is the Apache License 2.0 and the DESCRIPTION
    # of the package says "License: Apache License (== 2.0) | file LICENSE"
    # (https://github.com/saezlab/progeny, checked 2026-10-12), as the legacy
    # registry records.
    license=LicenseCV.APACHE_2_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='29295995',
    topics=(Topic.SIGNALING, Topic.GENE_REGULATION),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL, KnowledgeOrigin.COMPUTATIONALLY_PREDICTED),
    summary='Genes whose expression responds to the activity of 14 signaling pathways, with weights and p-values.',
    description=(
        'PROGENy leverages a large compendium of publicly available signaling '
        'perturbation experiments to yield a common core of pathway responsive '
        'genes. For each of 14 pathways it gives a collection of genes with '
        'their contribution (a weight) and significance (a p-value). Coupled '
        'with a statistical method, the weights infer pathway activities from '
        'bulk or single cell transcriptomics.'
    ),
)

# ================================== DOWNLOAD ==================================

# The legacy URL (github.com/.../raw/master/data/model_human_full.rda)
# redirects here.
download = Download(
    url='https://raw.githubusercontent.com/saezlab/progeny/master/data/model_human_full.rda',
    filename='progeny_model_human_full.rda',
    subfolder='progeny',
    large=True,
    default_mode='rb',
    encoding=None,
)

# =================================== SCHEMA ===================================

PATHWAY = 'progeny_pathway'
# The model is the human one and its genes are HGNC symbols (the package
# documentation of ``model_human_full``).
TAXON = 'NCBITaxon:9606'

# Scale of the two values, from the package documentation of
# ``model_human_full``: the weight is the "z-scores for a given gene" of a
# pathway, positive where the gene is up-regulated when the pathway signals and
# negative where it is down-regulated (the vignette: "the gene CBX6 has a
# negative response to EGFR, meaning that when there is EGFR signaling this gene
# is down-regulated"); the p-value is the "significance of gene in pathway".
# The weights span -83.9 to 45.4 and the p-values 1.9e-71 to 1.0 in the file.
WEIGHT = 'progeny:weight'
P_VALUE = 'progeny:p_value'

# What the pathways are, from the vignette of the package (the table of
# the model has only their names).
DESCRIPTIONS = {
    'Androgen': 'involved in the growth and development of the male reproductive organs',
    'EGFR': (
        'regulates growth, survival, migration, apoptosis, proliferation, and '
        'differentiation in mammalian cells'
    ),
    'Estrogen': 'promotes the growth and development of the female reproductive organs',
    'Hypoxia': 'promotes angiogenesis and metabolic reprogramming when O2 levels are low',
    'JAK-STAT': 'involved in immunity, cell division, cell death, and tumor formation',
    'MAPK': 'integrates external signals and promotes cell growth and proliferation',
    'NFkB': 'regulates immune response, cytokine production and cell survival',
    'p53': 'regulates cell cycle, apoptosis, DNA repair and tumor suppression',
    'PI3K': 'promotes growth and proliferation',
    'TGFb': 'involved in development, homeostasis, and repair of most tissues',
    'TNFa': (
        'mediates haematopoiesis, immune surveillance, tumour regression and '
        'protection from infection'
    ),
    'Trail': 'induces apoptosis',
    'VEGF': 'mediates angiogenesis, vascular permeability, and cell migration',
    'WNT': 'regulates organ morphogenesis during development and tissue repair',
}

f = FieldConfig()


def _measure(field: str):
    def value(row) -> Measurement | None:
        return measurement(row.get(field), source_field=field)

    return value


# A gene whose expression responds to a pathway is associated with it: the
# model fits the response of the gene to the perturbation of the pathway, and
# the sign of the weight is the direction of that response, which is carried as
# the weight. No effect of the gene on the pathway, and no direction between the
# two, is stated, so the predicate is ``associated_with``.
response_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=Gene,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('gene')),
        ),
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=TAXON)),
    ),
    predicate=slots.associated_with,
    object=EntityBuilder(
        entity_type=OntologyClass,
        identifiers=IdentifiersBuilder(CV(term=PATHWAY, value=f('pathway'))),
    ),
    annotations=AnnotationsBuilder(
        CV(term=WEIGHT, value=_measure('weight')),
        CV(term=P_VALUE, value=_measure('p_value')),
    ),
)


def map_response(row: dict) -> Relation | None:
    """A gene-pathway statement; a pathway the module does not know raises."""
    if row.get('pathway') not in DESCRIPTIONS:
        raise ValueError(f'Unknown PROGENy pathway: {row.get("pathway")!r}')
    return response_schema(row)


pathway_schema = EntityBuilder(
    entity_type=OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=PATHWAY, value=f('pathway')),
        CV(term=Namespace.NAME, value=f('pathway')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.description, value=lambda row: DESCRIPTIONS.get(row['pathway'])),
    ),
)


def map_pathway(row: dict):
    """The entity of a pathway; a pathway the module does not know raises."""
    if row.get('pathway') not in DESCRIPTIONS:
        raise ValueError(f'Unknown PROGENy pathway: {row.get("pathway")!r}')
    return pathway_schema(row)


# ================================= RESOURCE ===================================

resource = Resource(
    config,
    terms=Dataset(
        download=download,
        mapper=map_pathway,
        raw_parser=iter_pathways,
    ),
    responses=Dataset(
        download=download,
        mapper=map_response,
        raw_parser=iter_model,
    ),
)

# ================================= REFERENCE ==================================

# gene  pathway  weight  p.value
# e.g. RFC2  EGFR  1.470647  0.001655
