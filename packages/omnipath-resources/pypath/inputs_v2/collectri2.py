"""
Parse CollecTRI2 and emit Relation records of transcription factor - target gene regulation.

CollecTRI2 is the successor of CollecTRI: a collection of transcriptional
regulatory interactions between human transcription factors and their target
genes, from the text mining of ExTRI2 and from curated resources. It is a
separate resource from CollecTRI (``collectri.py``). Two files, read together
by ``parsers/collectri2.py``, give for each pair the sources that support it,
with the PubMed IDs and the signs of their evidence, and the sign (``weight``)
of the signed regulome that the authors made from them.

The relation follows ``collectri.py``: the factor ``affects`` the target gene,
and the effect on the expression of the gene (increased or decreased) is set
only where the source states the sign for the pair. The weight of the
regulome is such a statement only where the signed evidence of the pair has a
majority for it; a tie of the evidence, a pair without any signed evidence
(whose weight is a default, or the sign of most of the factor's other pairs)
and a weight against the evidence give a bare ``affects``. The rule that makes
the weight is in the docstring of the parser.

A pair is supported by one or more of twelve sources, which the relation names
as ``supporting_data_source``. ExTRI2 is the text mining; the others are
curated databases and curation, so the source list tells text-mined support from
curated support. The publications are the PubMed IDs of all of them.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.collectri2 import PMID_MAJORITY, iter_pairs
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
    id=ResourceCv.COLLECTRI2,
    name='CollecTRI2',
    full='Collection of Transcriptional Regulatory Interactions 2',
    url='https://zenodo.org/records/19773408',
    # Zenodo record 19773408: license CC BY 4.0.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    topics=(Topic.GENE_REGULATION,),
    content=(Content.CAUSAL_REGULATION,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.TEXT_MINED,
    ),
    summary=(
        'Transcription factor-target gene regulation from text mining and '
        'curated resources, with the supporting sources, literature '
        'references and the mode of regulation where the evidence states it.'
    ),
    description=(
        'CollecTRI2 aggregates the transcriptional regulatory interactions '
        'between human transcription factors and target genes of ExTRI2, '
        'which extracts them from PubMed abstracts with transformer models '
        'and labels the mode of regulation of each sentence, and of '
        'curated resources: HTRI, TRRUST, TFactS, GOA, IntAct, SIGNOR, '
        'CytReg, GEREDB, Pavlidis2021, the curated level A of DoRothEA and '
        'manual curation at NTNU. The file of the aggregated interactions '
        'gives for each pair the supporting sources with their publications '
        'and signs, and a derived regulome gives a consensus sign of every '
        'pair for the inference of transcription factor activity; where a '
        'pair has no sign in its evidence, activation is the default.'
    ),
)

# Zenodo record 19773408 ("CollecTRI2 dataset and regulome", 2026-04-25) holds
# the two files. Zenodo answers repeated requests with gateway errors and
# refusals for a while.
_RECORD = 'https://zenodo.org/records/19773408/files'
download_pairs = Download(
    url=f'{_RECORD}/CollecTRI2.tsv.gz',
    filename='CollecTRI2.tsv.gz',
    subfolder='collectri2',
    large=True,
    ext='gz',
)
download_regulome = Download(
    url=f'{_RECORD}/CollecTRI2_regulome.tsv.gz',
    filename='CollecTRI2_regulome.tsv.gz',
    subfolder='collectri2',
    large=True,
    ext='gz',
)

# The factors and targets are human gene symbols; the file has no organism
# column (ExTRI2 maps the genes of human, mouse and rat to human orthologs).
TAXON = 'NCBITaxon:9606'
# As in CollecTRI: these two names are the transcription factors that are
# dimers, not gene symbols. The file does not list their subunits.
COMPLEXES = frozenset({'AP1', 'NFKB'})

_DIRECTION = model.DirectionQualifierEnum
_ASPECT = model.GeneOrGeneProductOrChemicalEntityAspectEnum
# ``weight`` of the regulome -> effect on the expression of the target gene.
_MODE_OF_REGULATION = {
    '1': _DIRECTION.increased,
    '-1': _DIRECTION.decreased,
}
# ``sign_decision`` of the parser -> whether the weight is a statement of the
# evidence of the pair (see the parser).
_SIGN_DECISIONS = {
    PMID_MAJORITY: True,
    'pmid_tie': False,
    'no_signed_evidence': False,
    'contradicted': False,
}
# The label of a source as the file writes it -> the label carried. The
# spelling is that of ``collectri.py`` where it names the same database
# (``NTNU.Curated``); ExTRI2 is the text mining of this resource.
SOURCES = {
    'ExTRI2': 'ExTRI2',
    'HTRI': 'HTRI',
    'TRRUST': 'TRRUST',
    'TFactS': 'TFactS',
    'GOA': 'GOA',
    'IntAct': 'IntAct',
    'SIGNOR': 'SIGNOR',
    'CytReg': 'CytReg',
    'GEREDB': 'GEREDB',
    'NTNU Curated': 'NTNU.Curated',
    'Pavlidis2021': 'Pavlidis2021',
    'DoRothEA_A': 'DoRothEA_A',
}


def _is_complex(row: dict[str, Any]) -> bool:
    return row.get('tf') in COMPLEXES


def _direction(row: dict[str, Any]):
    """The mode of regulation, where the evidence of the pair states it."""
    weight, decision = row.get('weight'), row.get('sign_decision')
    if weight not in _MODE_OF_REGULATION:
        raise ValueError(f'Unknown CollecTRI2 weight: {weight!r}')
    if decision not in _SIGN_DECISIONS:
        raise ValueError(f'Unknown CollecTRI2 sign decision: {decision!r}')
    return _MODE_OF_REGULATION[weight] if _SIGN_DECISIONS[decision] else None


def _aspect(row: dict[str, Any]):
    return _ASPECT.expression if _direction(row) is not None else None


def _sources(row: dict[str, Any]) -> list[str]:
    unknown = [label for label in row['sources'] if label not in SOURCES]
    if unknown:
        raise ValueError(f'Unknown CollecTRI2 source: {unknown[0]!r}')
    return [SOURCES[label] for label in row['sources']]


f = FieldConfig()

transcription_factor = EntityBuilder(
    entity_type=lambda row: (
        model.MacromolecularComplex if _is_complex(row) else model.Protein
    ),
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.GENESYMBOL,
            value=lambda row: None if _is_complex(row) else row.get('tf'),
        ),
        CV(
            term=Namespace.NAME,
            value=lambda row: row.get('tf') if _is_complex(row) else None,
        ),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=TAXON)),
)

# Protein coding, microRNA and other RNA genes alike are gene symbols here.
target_gene = EntityBuilder(
    entity_type=model.Gene,
    identifiers=IdentifiersBuilder(CV(term=Namespace.GENESYMBOL, value=f('target'))),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=TAXON)),
)

interactions_schema = RelationBuilder(
    subject=transcription_factor,
    predicate=slots.affects,
    object=target_gene,
    annotations=AnnotationsBuilder(
        # Activation or repression of the target gene's expression, only where
        # the evidence of the pair decides it. A tie, a default or an inferred
        # sign is not carried over, so those pairs have neither qualifier.
        CV(term=slots.object_direction_qualifier, value=_direction),
        CV(term=slots.object_aspect_qualifier, value=_aspect),
        CV(
            term=slots.publications,
            value=lambda row: [f'PMID:{pmid}' for pmid in row['pmids']],
        ),
        # The sources that support the pair: ExTRI2 is text mined, the others
        # curated.
        CV(term=slots.supporting_data_source, value=_sources),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=download_pairs,
        mapper=interactions_schema,
        raw_parser=functools.partial(iter_pairs, regulome=download_regulome),
    ),
)
