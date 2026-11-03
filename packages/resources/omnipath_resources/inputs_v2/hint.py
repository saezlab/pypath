"""
Parse the HINT high-quality interactomes and emit Relation records.

HINT publishes, for each organism, one table of binary interactions and one
of co-complex associations. A row is a protein pair with the list of its
evidence, each item written as ``pmid:method:quality:type``. The two tables
are separate datasets, because a co-complex association does not establish
that the two proteins touch.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.molecular_forms import molecular_form_from_identifiers
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2._publications import PMID_PATTERN
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv


config = ResourceConfig(
    id=ResourceCv.HINT,
    name='HINT',
    full='High-quality INTeractomes',
    url='https://hint.yulab.org/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='22846459',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    summary=(
        'High-quality human binary protein interactions and co-complex '
        'associations compiled by HINT, with their publications, detection '
        'methods and the kind of study (high-throughput or literature '
        'curated) of the evidence.'
    ),
    description=(
        'HINT (High-quality INTeractomes) is a compilation of protein-protein '
        'interactions from eight interactome resources (BioGRID, MINT, '
        'iRefWeb, DIP, IntAct, HPRD, MIPS and the PDB). Interactions are '
        'filtered to keep those from high-throughput studies that validated '
        'their screens and those from small-scale studies that at least two '
        'publications report. Binary interactions and co-complex '
        'associations are distributed separately for ten organisms.'
    ),
)

# Each release has its own address; this is the release the module reads.
_RELEASE = '2026-07'
_URL = 'https://hint.yulab.org/download-raw/{release}/{filename}'
_EVIDENCE = 'pmid:method:quality:type'
# The quality flag of an evidence item: where the evidence comes from. HINT's
# download page names the two kinds of each table "High-throughput" and
# "Literature curated". The flags of the evidence items of a pair are an
# attribute of its relation, since the union of the items would lose them.
_QUALITY = 'hint:quality'
_QUALITY_VALUES = {
    'HT': 'High-throughput',
    'LC': 'Literature curated',
}
_EVIDENCE_ITEM_RE = re.compile(r'[^:]*:\d{4}:([^:]*):[^:]*')
# A UniProt accession, with the isoform when the source names one. A chain
# (``P05067-PRO_0000000092``) is the entry's accession and a feature identifier.
_ACCESSION_RE = re.compile(
    r'^((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(?:-\d+)?)(?:-PRO_\d+)?$'
)


def _molecular_form(column: str):
    """The isoform or chain the source names; the builder sees the accession."""

    def _form(row: dict[str, Any]) -> dict | None:
        accession = str(row.get(column) or '').strip()
        return molecular_form_from_identifiers(
            [{'ns': 'uniprot', 'id': accession}]
        )

    return _form


def _qualities(row: dict[str, Any]) -> list[str]:
    """The distinct quality flags of the evidence items, named as above.

    As for the publications and methods, an item that is not laid out as
    ``pmid:method:quality:type`` contributes nothing; a flag that is not in
    ``_QUALITY_VALUES`` raises ``ValueError``.
    """
    values = []
    for item in str(row.get(_EVIDENCE) or '').split('|'):
        match = _EVIDENCE_ITEM_RE.fullmatch(item)
        if match is None:
            continue
        if match.group(1) not in _QUALITY_VALUES:
            raise ValueError(f'Unknown HINT quality flag: {match.group(1)!r}')
        values.append(_QUALITY_VALUES[match.group(1)])
    return list(dict.fromkeys(values))


f = FieldConfig(
    extract={
        # An evidence item is ``pmid:method:quality:type``; any other layout
        # contributes nothing. In the reference position the source also puts
        # PDB entries (``PDB_1JTI``) and the placeholder numbers of its own
        # sources (``888800001068``), which are not PubMed IDs.
        'pmid': rf'^({PMID_PATTERN}):\d{{4}}:[^:]*:[^:]*$',
        # The PSI-MI accession number of the detection method.
        'method': r'^[^:]*:(\d{4}):[^:]*:[^:]*$',
        'accession': _ACCESSION_RE,
        'taxon': r'^(\d+)$',
    },
    transform={
        'pmid': lambda v: f'PMID:{v}',
        'method': lambda v: f'MI:{v}',
        'taxon': lambda v: f'NCBITaxon:{v}',
    },
    delimiter='|',
)


def _protein(suffix: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=Protein,
        molecular_form=_molecular_form(f'Uniprot_{suffix}'),
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.UNIPROT,
                value=f(f'Uniprot_{suffix}', extract='accession'),
            ),
            CV(term=Namespace.GENESYMBOL, value=f(f'Gene_{suffix}')),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.in_taxon,
                value=f('taxid', extract='taxon', transform='taxon'),
            ),
        ),
    )


def _schema(predicate) -> RelationBuilder:
    return RelationBuilder(
        subject=_protein('A'),
        predicate=predicate,
        object=_protein('B'),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.has_evidence_of_type,
                value=f(_EVIDENCE, extract='method', transform='method'),
            ),
            CV(
                term=slots.publications,
                value=f(_EVIDENCE, extract='pmid', transform='pmid'),
            ),
            CV(term=_QUALITY, value=_qualities),
        ),
    )


def _interaction(schema: RelationBuilder):
    def _map(row: dict[str, Any]):
        # A gene name alone does not identify a participant.
        if not all(
            _ACCESSION_RE.match(str(row.get(column) or '').strip())
            for column in ('Uniprot_A', 'Uniprot_B')
        ):
            return None
        return schema(row)

    return _map


# Pairs HINT classifies as binary: a method that tests two proteins found them.
binary_schema = _schema(slots.directly_physically_interacts_with)
# Pairs found in the same purified complex. Their contact may be indirect.
cocomplex_schema = _schema(slots.interacts_with)


def _dataset(filename: str, mapper: RelationBuilder) -> Dataset:
    return Dataset(
        download=Download(
            url=_URL.format(release=_RELEASE, filename=filename),
            filename=f'{_RELEASE}_{filename}',
            subfolder='hint',
            large=True,
        ),
        mapper=_interaction(mapper),
        raw_parser=iter_tsv,
    )


resource = Resource(
    config,
    human_binary=_dataset('HomoSapiens_binary_hq.txt', binary_schema),
    human_cocomplex=_dataset('HomoSapiens_cocomp_hq.txt', cocomplex_schema),
)
