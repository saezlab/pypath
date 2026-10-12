"""
Parse the DIP human interactions and emit Relation records.

DIP (Database of Interacting Proteins) catalogues experimentally determined
protein-protein interactions. Its PSI-MI TAB file of the interactions that
involve a human protein (read by ``parsers/dip.py``) has one row for each DIP
interaction: two DIP nodes, the organism of each, and the evidence items that
support the interaction (detection method, interaction type, publication and
processing status), listed in parallel in four columns.

This module emits every row as one relation between the two proteins. The
quality status ``core`` that DIP gives to its most reliable interactions is
kept as the confidence level of the relation. The core interactions are a
subset of this file, which is the only one read, so they are not read twice.

DIP's own download pages need a login. The file is the copy of the release of
5 February 2017 (``Hsapi20170205``) on the OmniPath rescued-data server.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.molecular_forms import molecular_form_from_identifiers
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2 import _mitab
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.dip import iter_dip
from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.DIP,
    name='DIP',
    url='https://dip.doe-mbi.ucla.edu/',
    # The home page of DIP links the Creative Commons Attribution-NoDerivs 3.0
    # license (CC BY-ND 3.0). The footer of the site adds that redistribution
    # of the data, or encapsulation of it in another database, needs the
    # permission of DIP.
    license=LicenseCV.CC_BY_ND_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='11752321',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    summary=(
        'Experimentally determined protein-protein interactions involving '
        'human proteins, with detection methods, interaction types, '
        'publications and the DIP quality status, February 2017 release.'
    ),
    description=(
        'The Database of Interacting Proteins (DIP) at UCLA catalogues '
        'protein-protein interactions determined by experiment, curated '
        'from the literature and checked by a second curator. Its "core" '
        'set is the subset of interactions that DIP assesses as most '
        'reliable. The file read here holds the interactions in which at '
        'least one protein is human, including interactions with proteins '
        'of other species.'
    ),
)

# The interaction type column of DIP holds some terms outside the direct
# interaction, association and physical association types that the shared
# MITAB rules map:
# - MI:0218, "physical interaction", is obsolete ("can be direct or
#   indirect"); PSI-MI split it into association and physical association.
# - MI:0254, "genetic interference", is a detection method, not an
#   interaction type. DIP lists the record as an interaction between the two
#   proteins; only the broad predicate holds.
_BROAD_INTERACTIONS = {'MI:0218', 'MI:0254'}
# Processing status of an evidence item: code and label as the file writes
# them (``dip:0002(small scale)``).
_PROCESSING_STATUSES = {
    ('0002', 'small scale'),
    ('0004', 'small scale'),
    ('0005', 'high throughput'),
}
_PROCESSING_STATUS = 'dip:processing_status'

_DIP_NODE_RE = re.compile(r'DIP-\d+N')
_DIP_INTERACTION_RE = re.compile(r'DIP-\d+E')
# A UniProt accession, with the isoform number or the PRO identifier of a
# chain where the record is on one.
_UNIPROT_RE = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(?:-\d+|-PRO_\d+)?'
)


def _identifier_pairs(side: str):
    """The participant's UniProt and RefSeq accessions, then its DIP node.

    An accession of another form is left out: the file has a UniProt
    accession written ``Q13936,159``, and RefSeq entries that hold the word
    ``refseq``, a GenBank protein (``AAH30133``) or a UniProt accession.
    """

    def _pairs(row: dict[str, Any]) -> list[tuple[Any, str]]:
        items = _mitab.parse_mitab_cell(row.get(f'id_{side}'))
        uniprot = [
            (Namespace.UNIPROT, item.accession)
            for item in items
            if item.database == 'uniprotkb' and _UNIPROT_RE.fullmatch(item.accession)
        ]
        refseq = [
            pair
            for item in items
            if item.database == 'refseq'
            and (pair := _mitab.identifier_pair(item)) is not None
        ]
        node = [
            (Namespace.DIP, item.accession)
            for item in items
            if item.database is None and _DIP_NODE_RE.fullmatch(item.accession)
        ]
        return list(dict.fromkeys(uniprot + refseq + node))

    return _pairs


def _participant_form(side: str):
    """The isoform or chain the UniProt accession names, if it names one.

    The RefSeq accessions of the file are cross-references of the DIP node
    and are not taken as the sequence of the participant.
    """
    pairs = _identifier_pairs(side)

    def _form(row: dict[str, Any]) -> dict | None:
        return molecular_form_from_identifiers(
            [
                {'ns': str(ns), 'id': accession}
                for ns, accession in pairs(row)
                if ns == Namespace.UNIPROT
            ]
        )

    return _form


def _interaction_records(row: dict[str, Any]) -> list[tuple[Any, str]]:
    return [
        (Namespace.DIP, item.accession)
        for item in _mitab.parse_mitab_cell(row.get('interaction_ids'))
        if item.database is None and _DIP_INTERACTION_RE.fullmatch(item.accession)
    ]


def _processing_statuses(row: dict[str, Any]) -> list[str]:
    """The processing statuses of the evidence items, with DIP's code."""
    statuses = []
    for item in _mitab.parse_mitab_cell(row.get('processing_status')):
        if item.database != 'dip' or (
            item.accession,
            item.description,
        ) not in _PROCESSING_STATUSES:
            raise ValueError(
                f'Unknown DIP processing status: {item.database}:{item.accession}'
                f'({item.description})'
            )
        statuses.append(f'{item.description} (dip:{item.accession})')
    return list(dict.fromkeys(statuses))


def _predicate(row: dict[str, Any]):
    """
    Predicate for the PSI-MI interaction types of the row.

    The types are those of the evidence items. A direct interaction or one
    of its kinds (the enzymatic reactions, covalent binding) among them is a
    direct physical interaction; an association, a physical association, the
    obsolete "physical interaction" or genetic interference alone says only
    that the proteins interact. Any other type gives None.
    """
    types = _mitab.mi_terms(row.get('interaction_types'))
    predicate = _mitab.physical_predicate(types)
    if predicate is None and types and set(types) <= _BROAD_INTERACTIONS:
        return slots.interacts_with
    return predicate


def _participant(side: str) -> EntityBuilder:
    return EntityBuilder(
        # DIP is a database of proteins; the file has no interactor type.
        entity_type=model.Protein,
        molecular_form=_participant_form(side),
        identifiers=IdentifiersBuilder(CV.from_pairs(_identifier_pairs(side))),
        annotations=AnnotationsBuilder(
            *_mitab.participant_annotations(taxon=f'taxid_{side}'),
        ),
    )


interactions_schema = RelationBuilder(
    subject=_participant('a'),
    predicate=_predicate,
    object=_participant('b'),
    identifiers=IdentifiersBuilder(CV.from_pairs(_interaction_records)),
    annotations=AnnotationsBuilder(
        *_mitab.interaction_annotations(
            detection_method='detection_methods',
            interaction_type='interaction_types',
            publications='publications',
        ),
        # "core" for the interactions DIP marks as most reliable; the others
        # have no status.
        CV(term=slots.has_confidence_level, value=lambda row: row.get('quality_status')),
        CV(term=_PROCESSING_STATUS, value=_processing_statuses),
    ),
)


def exclusion_reason(row: dict[str, Any]) -> str | None:
    """Why a row gives no relation; None for a row that is mapped."""
    for side in 'ab':
        if not _identifier_pairs(side)(row):
            return 'participant without an identifier'
    if not _interaction_records(row):
        return 'no DIP interaction identifier'
    if _predicate(row) is None:
        return 'unknown interaction type'
    return None


def map_interaction(row: dict[str, Any]):
    return None if exclusion_reason(row) else interactions_schema(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/dip/Hsapi20170205.txt.gz',
            filename='Hsapi20170205.txt.gz',
            subfolder='dip',
            large=True,
            ext='gz',
        ),
        mapper=map_interaction,
        raw_parser=iter_dip,
    ),
)
