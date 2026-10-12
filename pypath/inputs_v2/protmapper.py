"""
Parse the ProtMapper corpus of phosphosite annotations and emit Relation records.

The corpus is the supplementary data of Bachman et al. (bioRxiv
10.1101/822668, version 4). It has two tables, joined here on their ``ID``
column:

- ``export.csv``: one row per site annotation, a controller that
  phosphorylates a residue of a human target protein;
- ``evidences.csv``: one row per piece of evidence for an annotation, from a
  curated database or a text mining reader.

Each evidence row becomes one record, so that the source, the publication
and the site of an observation stay together. ProtMapper is an aggregator:
the source of each record is the database or the reader it names.

An evidence row whose own site (``MAPPED_*``, else ``ORIG_*``) differs from
the annotation's is not attached to it, because its publication and source
do not support that site; ``parsers/protmapper.py`` explains and counts them.
An annotation left without evidence is still a record, without source and
publication. The controller's ``CTRL_IS_KINASE`` flag and the ``DESCRIPTION``
of the evidence's site mapping are attributes of the relation.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.protmapper import (
    ANNOTATIONS_FILE,
    EVIDENCES_FILE,
    iter_evidences,
)
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

# The values of ``CTRL_NS``, the identifier system of the controller, with the
# class its identifiers name. A FamPlex entry is a protein family or a
# complex; the table does not say which, so it gets no more specific class.
_CONTROLLER_TYPES = {
    'UP': model.Protein,
    'HGNC': model.Gene,
    'FPLX': model.NamedThing,
}
# FamPlex has no member in ``Namespace``.
FAMPLEX_NAMESPACE = 'fplx'
# The values of ``SOURCE``, kept as the source writes them. The first six are
# curated databases (PhosphoSitePlus, SIGNOR, HPRD, NCI-PID, Reactome and the
# BEL Large Corpus), the last three text mining readers (REACH, Sparser and
# RLIMS-P).
SOURCES = frozenset(
    {
        'psp',
        'signor',
        'hprd',
        'pid',
        'reactome',
        'bel',
        'reach',
        'sparser',
        'rlimsp',
    }
)
# The authors kept annotations whose substrate is a human protein and gave
# every site its position in the human UniProt sequence. The controller is
# not restricted to human: its UniProt accessions include mouse, rat, viral
# and bacterial proteins, and the tables have no organism column.
PROTMAPPER_TARGET_TAX_ID = '9606'
# ``CTRL_IS_KINASE`` of an annotation: whether the controller is a known human
# kinase. Every value of version 4 of the corpus is listed.
_IS_KINASE = {'True': 'true', 'False': 'false'}
# ``DESCRIPTION`` of an evidence row: "VALID" if its site was valid as it
# was, otherwise a code or a note for the mapping that gave the annotation's
# site. The values of version 4 of the corpus, kept as the source writes them: the
# codes of ProtMapper (INFERRED_*, REMAPPED_FROM_PSP_SEQUENCE,
# SIGNAL_PEPTIDE_REMOVED) and notes of its curated site corrections.
SITE_MAPPINGS = frozenset(
    {
        'VALID',
        'INFERRED_ALTERNATIVE_ISOFORM',
        'INFERRED_MOUSE_SITE',
        'INFERRED_METHIONINE_CLEAVAGE',
        'INFERRED_RAT_SITE',
        'REMAPPED_FROM_PSP_SEQUENCE',
        'SIGNAL_PEPTIDE_REMOVED',
        'numbering from mouse sequence',
        'numbering after 24-residue signaling peptide removed',
        'numbering after 30-residue signaling peptide removed',
        'numbering from isoform alpha 2',
        'numbering from isoform p52',
        'numbering from isoform p46',
        'numbering from isoform zeta',
        'numbering from chicken sequence',
        'numbering from deprecated sequence GenBank A39248',
        'numbering from off-by-one mouse sequence',
        'numbering from isoform 3, supported by UniProt, Phosphosite',
        'numbering from isoform alpha 2; wrong residue',
        'position from isoform 8',
        'mouse sequence',
        'off by one',
        'off by 1',
        'off by 2',
        'off by 3',
        'off by four',
        'off by 39 (cleavage?)',
        'off by 1 (methionine removed)',
        'wrong residue',
        'probable wrong residue; numbering from isoform alpha 2',
        'immunogen using misnumbered SGpSVS sequence',
        'Listed as T172/183 in some online resources, not clear where T172 comes from',
    }
)
_IS_KINASE_TERM = 'protmapper:ctrl_is_kinase'
_SITE_MAPPING_TERM = 'protmapper:site_mapping'
# Residue letter and position of a site, as in ``TARGET_RES`` and ``TARGET_POS``.
_SITE_REGEX = re.compile(r'^[A-Z][1-9]\d*$')


def _controller_type(row: dict[str, Any]) -> Any:
    namespace = row.get('CTRL_NS')
    if namespace not in _CONTROLLER_TYPES:
        raise ValueError(f'Unknown ProtMapper controller namespace: {namespace!r}')
    return _CONTROLLER_TYPES[namespace]


def _controller_id(namespace: str):
    """``CTRL_ID`` of the rows whose controller is of this identifier system."""
    return lambda row: row.get('CTRL_ID') if row.get('CTRL_NS') == namespace else None


def _name(id_column: str, name_column: str):
    """The name of a participant that has an identifier; a name alone names none."""

    def name(row: dict[str, Any]) -> str | None:
        if not (row.get(id_column) or '').strip():
            return None
        return (row.get(name_column) or '').strip() or None

    return name


def _controller_taxon(row: dict[str, Any]) -> str | None:
    # An HGNC identifier names a human gene. A UniProt controller can be of
    # any species, and the source does not give it.
    return 'NCBITaxon:9606' if row.get('CTRL_NS') == 'HGNC' else None


def _source(row: dict[str, Any]) -> str | None:
    value = row.get('SOURCE')
    if value is None:  # an annotation without evidence rows
        return None
    if value not in SOURCES:
        raise ValueError(f'Unknown ProtMapper evidence source: {value!r}')
    return value


def _is_kinase(row: dict[str, Any]) -> str:
    value = row.get('CTRL_IS_KINASE')
    if value not in _IS_KINASE:
        raise ValueError(f'Unknown ProtMapper CTRL_IS_KINASE value: {value!r}')
    return _IS_KINASE[value]


def _site_mapping(row: dict[str, Any]) -> str | None:
    value = row.get('DESCRIPTION')
    if value is None:  # an annotation without evidence rows
        return None
    if value not in SITE_MAPPINGS:
        raise ValueError(f'Unknown ProtMapper site mapping description: {value!r}')
    return value


def _target_form(row: dict[str, Any]) -> dict | None:
    """The phosphorylated residue, a feature of the target in this record.

    ``TARGET_RES`` and ``TARGET_POS`` are the site after ProtMapper's
    normalization: a position in the sequence of ``TARGET_UP_ID``, which is
    an isoform accession where the source gives one. The site an evidence
    first reported (``ORIG_*``) is not used.
    """
    if not row.get('TARGET_UP_ID'):
        return None
    residue, position = row.get('TARGET_RES') or '', row.get('TARGET_POS') or ''
    if not _SITE_REGEX.match(f'{residue}{position}'):
        raise ValueError(f'Malformed ProtMapper site: {residue!r} {position!r}')
    return {
        'modifications': [
            {
                'term': 'phosphorylated residue',
                'residue': residue,
                'position': int(position),
                'end_position': int(position),
                'coordinate_reference': {
                    'identifier': {'ns': 'uniprot', 'id': row['TARGET_UP_ID']},
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
            }
        ],
    }


config = ResourceConfig(
    id=ResourceCv.PROTMAPPER,
    name='ProtMapper',
    url='https://www.biorxiv.org/content/10.1101/822668',
    # License as recorded in the legacy registry (resources.json); the previous value (CC BY-NC-SA 4.0) had no source statement recorded.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    topics=(Topic.SIGNALING,),
    content=(Content.ENZYME_PTM_RELATIONSHIPS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.TEXT_MINED,
    ),
    summary=(
        'Regulators of human phosphorylation sites from curated databases '
        'and text mining, with the source and publication of each evidence.'
    ),
    description=(
        'ProtMapper is a tool of the INDRA project that maps the positions '
        'of phosphorylation sites to human reference sequences. Its authors '
        'used it to assemble a corpus of annotations, each a regulator that '
        'phosphorylates one residue of a human protein, from six curated '
        'databases (PhosphoSitePlus, SIGNOR, HPRD, NCI-PID, Reactome and the '
        'BEL Large Corpus) and three text mining systems (REACH, Sparser '
        'and RLIMS-P). The regulators are proteins, genes and FamPlex '
        'protein families and complexes. Every annotation lists its '
        'evidence by source, with publications, and has a belief score '
        'calculated by INDRA.'
    ),
)

f = FieldConfig()

controller = EntityBuilder(
    entity_type=_controller_type,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=_controller_id('UP')),
        CV(term=Namespace.HGNC, value=_controller_id('HGNC')),
        CV(term=FAMPLEX_NAMESPACE, value=_controller_id('FPLX')),
        # The name INDRA gave the controller: a gene symbol for a human
        # protein or gene, another gene or protein name otherwise.
        CV(term=Namespace.NAME, value=_name('CTRL_ID', 'CTRL_GENE_NAME')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=_controller_taxon),
    ),
)

target = EntityBuilder(
    entity_type=model.Protein,
    # The phosphorylated residue is in the target's molecular form.
    molecular_form=_target_form,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=f('TARGET_UP_ID')),
        CV(term=Namespace.NAME, value=_name('TARGET_UP_ID', 'TARGET_GENE_NAME')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=f'NCBITaxon:{PROTMAPPER_TARGET_TAX_ID}'),
    ),
)

# Every annotation of the corpus is a phosphorylation: the authors assembled
# it from INDRA Phosphorylation statements with a residue and a position.
enzyme_substrate_schema = RelationBuilder(
    subject=controller,
    predicate=slots.affects,
    object=target,
    annotations=AnnotationsBuilder(
        # The controller increases the target's phosphorylation; the corpus
        # says nothing about the target's activity, and not every controller
        # is a kinase.
        CV(
            term=slots.object_direction_qualifier,
            value=model.DirectionQualifierEnum.increased,
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=model.GeneOrGeneProductOrChemicalEntityAspectEnum.phosphorylation,
        ),
        CV(
            term=slots.publications,
            value=f('PMID', extract=PMID_EXTRACT, transform=lambda pmid: f'PMID:{pmid}'),
        ),
        # The database or reader this evidence comes from.
        CV(term=slots.supporting_data_source, value=_source),
        # Whether the controller is a known human kinase (annotation level),
        # and how the site of this evidence was mapped to the annotation's.
        CV(term=_IS_KINASE_TERM, value=_is_kinase),
        CV(term=_SITE_MAPPING_TERM, value=_site_mapping),
        # INDRA's belief in the annotation, between 0 and 1; it is the same
        # on every evidence row of an annotation.
        CV(
            term=slots.has_confidence_score,
            value=lambda row: measurement(row.get('BELIEF'), source_field='BELIEF'),
        ),
    ),
)

resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=Download(
            url=(
                'https://www.biorxiv.org/content/biorxiv/early/2022/10/29/'
                '822668/DC1/embed/media-1.zip?download=true'
            ),
            filename='protmapper_822668v4_media-1.zip',
            subfolder='protmapper',
            large=True,
            ext='zip',
            needed=[ANNOTATIONS_FILE, EVIDENCES_FILE],
        ),
        mapper=enzyme_substrate_schema,
        raw_parser=iter_evidences,
    ),
)
