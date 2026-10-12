"""
Parse the DoRothEA database and emit Relation records.

Source: the ``entire_database`` data frame of the DoRothEA R package, one row
per human transcription factor - target gene pair with everything the package
knows about it:

    tf    target  mor  confidence  is_evidence_curated  which_curated   pubmed_id
    AHR   CYP1A1  1    A           TRUE                 htri_db,trrust  9698073,...

- ``mor``: mode of regulation, 1 activation, -1 repression, 0 not known;
- ``confidence``: level A (highest) to E, from the number and kind of
  supporting evidence;
- ``is_evidence_curated``, ``_chip_seq``, ``_tfbs``, ``_inferred``: the four
  lines of evidence (curated databases, ChIP-seq peaks, binding motifs in the
  promoter, and inference from GTEx expression), with the data sources of each
  in ``which_curated``, ``which_chip_seq``, ``which_tfbs``, ``which_inferred``.

The package's regulon tables (``dorothea_hs``, ``dorothea_mm``) are made from
this one: they have a filled in mode of regulation for every pair and no
evidence columns, so they are not read here.

What a pair states depends on the lines of evidence it has (Garcia-Alonso et
al. 2019: "(1) literature-curated resources, (2) ChIP-seq binding data, (3)
TFBS predictions on human promoters, and (4) reverse-engineered regulons
inferred from normal tissue gene expression profiles from [GTEx]"). A row with
several takes the strongest:

- curated: a regulatory statement of a curated resource, so ``affects``, with
  the expression aspect and the direction where ``mor`` gives the sign;
- ChIP-seq peaks or a motif match in the promoter, without a curated
  statement: the factor binds at or near the gene, which does not show an
  effect on its expression, so ``interacts_with``;
- inference from GTEx expression alone: a statistical association, so
  ``associated_with``.

The lines of evidence of the row are on the relation (``dorothea:evidence``)
so that this rule can be applied again. The table gives a sign only to
curated pairs; a sign on a pair without curated evidence is an error.
"""

from __future__ import annotations

import functools
import re
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import PMID_PATTERN
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.dorothea import CONFIDENCE_LEVELS, iter_level
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

# The table lists human interactions, with HGNC symbols; it has no organism
# column.
DOROTHEA_TAX_ID = '9606'

_DIRECTION = model.DirectionQualifierEnum
_ASPECT = model.GeneOrGeneProductOrChemicalEntityAspectEnum
# mor -> effect on the target gene's expression; 0 is "not known".
_MODE_OF_REGULATION = {
    1: _DIRECTION.increased,
    -1: _DIRECTION.decreased,
    0: None,
}
# Lines of evidence with an exact Evidence & Conclusion Ontology term. Curated
# databases and expression based inference have none; their data sources are
# still listed with the others.
_EVIDENCE_TYPES = {
    # chromatin immunoprecipitation-seq evidence
    'is_evidence_chip_seq': 'ECO:0000229',
    # position-specific scoring matrix motif search evidence
    'is_evidence_tfbs': 'ECO:0005536',
}
# The four lines of evidence, by the name of the flag column without its
# ``is_evidence_`` prefix.
_EVIDENCE_KINDS = ('curated', 'chip_seq', 'tfbs', 'inferred')
# Source-prefixed attribute: the lines of evidence of the row, with the
# source's names of the flags.
_EVIDENCE = 'dorothea:evidence'
_SOURCE_COLUMNS = (
    'which_curated',
    'which_chip_seq',
    'which_tfbs',
    'which_inferred',
)
# "none" stands for no data source of that line of evidence.
_NO_SOURCE = 'none'
# A reference is a PubMed ID, or IntAct's "imex:IM-14110|pubmed:20478527".
_PUBMED_REGEX = re.compile(rf'^(?:pubmed:)?({PMID_PATTERN})$')


def _confidence(row: dict[str, Any]) -> str:
    level = row.get('confidence')
    if level not in CONFIDENCE_LEVELS:
        raise ValueError(f'Unknown DoRothEA confidence level: {level!r}')
    return level


def _flag(row: dict[str, Any], kind: str) -> bool:
    column = f'is_evidence_{kind}'
    flag = row.get(column)
    if not isinstance(flag, bool):
        raise ValueError(f'DoRothEA {column} is not a boolean: {flag!r}')
    return flag


def _evidence_kinds(row: dict[str, Any]) -> list[str]:
    """The lines of evidence of the row, by the names of the source's flags."""
    return [kind for kind in _EVIDENCE_KINDS if _flag(row, kind)]


def _predicate(row: dict[str, Any]):
    """What the strongest line of evidence of the row states.

    Curated: a regulatory statement. ChIP-seq or a motif match: binding, not
    an effect on expression. Only inference from expression: an association.
    """
    kinds = _evidence_kinds(row)
    if 'curated' in kinds:
        return slots.affects
    if 'chip_seq' in kinds or 'tfbs' in kinds:
        return slots.interacts_with
    if kinds:
        return slots.associated_with
    raise ValueError('DoRothEA row without a line of evidence')


def _direction(row: dict[str, Any]):
    """The mode of regulation the curated sources give, None if not known."""
    mor = row.get('mor')
    if isinstance(mor, bool) or mor not in _MODE_OF_REGULATION:
        raise ValueError(f'Unknown DoRothEA mode of regulation: {mor!r}')
    direction = _MODE_OF_REGULATION[mor]
    # Only a curated statement has a sign; the table gives none to the rest.
    if direction is not None and not _flag(row, 'curated'):
        raise ValueError(
            f'DoRothEA mode of regulation {mor!r} without curated evidence'
        )
    return direction


def _aspect(row: dict[str, Any]):
    return _ASPECT.expression if _direction(row) is not None else None


def _evidence_types(row: dict[str, Any]) -> list[str]:
    return [
        term
        for column, term in _EVIDENCE_TYPES.items()
        if _flag(row, column.removeprefix('is_evidence_'))
    ]


def _data_sources(row: dict[str, Any]) -> list[str]:
    """The data sources of all four lines of evidence, as the table writes them."""
    return [
        source.strip()
        for column in _SOURCE_COLUMNS
        for source in str(row.get(column) or '').split(',')
        if source.strip() and source.strip() != _NO_SOURCE
    ]


def _publications(row: dict[str, Any]) -> list[str]:
    """PubMed IDs of the reference field; "-", "N/A" and IMEx IDs are none."""
    pmids = []
    for token in re.split(r'[,|]', str(row.get('pubmed_id') or '')):
        match = _PUBMED_REGEX.match(token.strip())
        if match:
            pmids.append(f'PMID:{match.group(1)}')
    return pmids


config = ResourceConfig(
    id=ResourceCv.DOROTHEA,
    name='DoRothEA',
    url='https://saezlab.github.io/dorothea/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.COMPOSITE,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='31340985',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.CAUSAL_REGULATION,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Human transcription factor–target gene pairs with confidence '
        'levels A to E, the lines of evidence and their sources; curated '
        'pairs are regulation (with the mode where known), pairs with only '
        'ChIP-seq or motif evidence are binding, and pairs with only '
        'expression inference are associations.'
    ),
    description=(
        'DoRothEA is a gene regulatory network of human transcription '
        'factors and their target genes, built from four lines of evidence: '
        'literature curated databases, ChIP-seq binding peaks, transcription '
        'factor binding motifs in promoters and regulation inferred from '
        'gene expression across GTEx tissues. Each interaction has a '
        'confidence level from A (the best supported curated interactions) '
        'to E (computational predictions only). The mode of regulation, '
        'activation or repression, is known only for part of the curated '
        'interactions. A pair with a curated line of evidence is a '
        'regulation of the target by the factor; a pair with only ChIP-seq '
        'peaks or motif matches is a binding of the factor at the gene, '
        'which does not show an effect on expression; a pair with only '
        'inference from expression is an association.'
    ),
)

f = FieldConfig()


def _participant(entity_type, column: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=entity_type,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(column)),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=f'NCBITaxon:{DOROTHEA_TAX_ID}'),
        ),
    )


# A row without a transcription factor or a target yields no relation: the
# participant has no identifier, so the builder returns None.
interactions_schema = RelationBuilder(
    subject=_participant(model.Protein, 'tf'),
    predicate=_predicate,
    object=_participant(model.Gene, 'target'),
    annotations=AnnotationsBuilder(
        # Activation or repression of the target gene's expression where the
        # table states it. Most pairs have mor 0 and get neither qualifier:
        # a binding site or a correlation does not say which way.
        CV(term=_EVIDENCE, value=_evidence_kinds),
        CV(term=slots.object_direction_qualifier, value=_direction),
        CV(term=slots.object_aspect_qualifier, value=_aspect),
        CV(term=slots.has_confidence_level, value=_confidence),
        CV(term=slots.has_evidence_of_type, value=_evidence_types),
        CV(term=slots.supporting_data_source, value=_data_sources),
        CV(term=slots.publications, value=_publications),
    ),
)

_download = Download(
    url='https://github.com/saezlab/dorothea/raw/master/data/entire_database.rda',
    filename='dorothea_entire_database.rda',
    subfolder='dorothea',
    large=True,
    default_mode='rb',
    encoding=None,
)


def _level(level: str) -> Dataset:
    return Dataset(
        download=_download,
        mapper=interactions_schema,
        raw_parser=functools.partial(iter_level, level=level),
    )


# One dataset per confidence level, so that leaving out the low levels is a
# choice of datasets. Together they are the whole table.
resource = Resource(
    config,
    level_a=_level('A'),
    level_b=_level('B'),
    level_c=_level('C'),
    level_d=_level('D'),
    level_e=_level('E'),
)
