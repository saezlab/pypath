"""
Parse the InnateDB curated interactions and emit Relation records.

InnateDB distributes the interactions its own curators annotated from the
literature as a PSI-MI TAB 2.7 file with a header line. One row is one
experimental observation of an interaction between two participants, each
identified by an InnateDB gene record. The file has protein-protein,
protein-DNA and protein-RNA interactions.

A second download, "all interactions", adds the content InnateDB imported
from BioGRID, IntAct, MINT, DIP, BIND and other databases. It is not read
here: those records are the evidence of the databases they come from.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

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
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2 import _mitab
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.innatedb import iter_interactions


config = ResourceConfig(
    id=ResourceCv.INNATEDB,
    name='InnateDB',
    url='https://www.innatedb.com/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.DSL,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='23180781',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    summary=(
        'Human and mouse molecular interactions of the innate immune '
        'response curated by InnateDB from the literature, with their '
        'detection methods, interaction types and publications.'
    ),
    description=(
        'InnateDB is a database of the genes, proteins, experimentally '
        'verified interactions and signaling pathways involved in the innate '
        'immune response of humans, mice and bovines to microbial infection. '
        'Its curators annotate interactions of innate immunity genes from the '
        'literature in detail, including protein-protein, protein-DNA and '
        'protein-RNA interactions, and a later curation project added '
        'interactions relevant to allergy and asthma. InnateDB also '
        'integrates interactions and pathways from other public databases.'
    ),
)

# Identifier systems of InnateDB without a registered namespace: its gene
# records (``IDBG-25842``) and its interaction records (``IDB-113240``).
_INNATEDB_GENE = 'innatedb_gene'
_INNATEDB_INTERACTION = 'innatedb'
# PSI-MI name of InnateDB as a source database.
_INNATEDB_SOURCE = 'MI:0974'

_GENE_RECORD_RE = re.compile(r'IDBG-\d+')
_INTERACTION_RECORD_RE = re.compile(r'IDB-\d+')
_UNIPROT_ACCESSION_RE = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
)
_UNIPROT_ENTRY_RE = re.compile(r'[A-Z0-9]{1,6}_[A-Z0-9]{1,5}')

# InnateDB names a DNA participant by the gene whose DNA took part, and all
# its participants by gene records.
_INTERACTOR_TYPES = {
    'MI:0326': model.Protein,  # protein
    'MI:0319': model.Gene,  # deoxyribonucleic acid
    'MI:0320': model.RNAProduct,  # ribonucleic acid
}
# Accessions the file has in the interaction type column that are not PSI-MI
# interaction types: an obsolete feature type and a detection method. The
# record still reports an interaction, so only the broad predicate holds.
_UNTYPED_INTERACTIONS = {
    'MI:0179',  # other modification (obsolete)
    'MI:0401',  # biochemical
}


def _id_column(suffix: str) -> str:
    # The header line starts with a hash.
    prefix = '#' if suffix == 'A' else ''
    return f'{prefix}unique_identifier_{suffix}'


def _alias_pair(item: _mitab.MitabItem) -> tuple[Any, str] | None:
    """
    One entry of the alias column.

    The column lists the UniProt accessions and entry names and the RefSeq
    proteins of the gene, and its symbol with the nomenclature it comes from
    (``hgnc:TLR4(display_short)``, ``mgi:Tlr4(display_short)``; from Ensembl
    where neither has one). A gene without a symbol gets its own record
    number as display name, which is not a name.
    """
    if item.database == 'uniprotkb':
        if _UNIPROT_ACCESSION_RE.fullmatch(item.accession):
            return Namespace.UNIPROT, item.accession
        if _UNIPROT_ENTRY_RE.fullmatch(item.accession):
            return Namespace.UNIPROT_ENTRY, item.accession
        return None
    if item.description == 'display_short':
        if item.database in {'hgnc', 'mgi', 'ensembl'}:
            return Namespace.GENESYMBOL, item.accession
        return None
    return _mitab.identifier_pair(item)


def _identifier_pairs(suffix: str):
    """The InnateDB gene record, then the cross-references the row lists."""

    def _pairs(row: dict[str, Any]) -> list[tuple[Any, str]]:
        pairs = [
            (_INNATEDB_GENE, item.accession)
            for item in _mitab.parse_mitab_cell(row.get(_id_column(suffix)))
            if item.database == 'innatedb'
            and _GENE_RECORD_RE.fullmatch(item.accession)
        ]
        if not pairs:
            return []
        pairs += _mitab.identifier_pairs(row.get(f'alt_identifier_{suffix}'))
        aliases = _mitab.parse_mitab_cell(row.get(f'alias_{suffix}'))
        pairs += filter(None, map(_alias_pair, aliases))
        return list(dict.fromkeys(pairs))

    return _pairs


def _entity_type(suffix: str):
    column = f'interactor_type_{suffix}'
    return lambda row: _mitab.interactor_type(row.get(column), _INTERACTOR_TYPES)


def _interaction_records(row: dict[str, Any]) -> list[str]:
    # The database name is written "innatedb", "InnateDB" or, for the allergy
    # and asthma curation project, "InnateDB Allergy".
    return [
        item.accession
        for item in _mitab.parse_mitab_cell(row.get('idinteraction_in_source_db'))
        if item.database in {'innatedb', 'innatedb allergy'}
        and _INTERACTION_RECORD_RE.fullmatch(item.accession)
    ]


def _predicate(row: dict[str, Any]):
    """
    Predicate for the PSI-MI interaction type of the row.

    A direct interaction or one of its kinds (the enzymatic reactions,
    covalent binding, self interaction) is a direct physical interaction;
    an association or physical association says only that the participants
    interact. Any other interaction type gives None.
    """
    types = _mitab.mi_terms(row.get('interaction_type'))
    predicate = _mitab.physical_predicate(types)
    if predicate is None and types and set(types) <= _UNTYPED_INTERACTIONS:
        return slots.interacts_with
    return predicate


def _participant(suffix: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=_entity_type(suffix),
        identifiers=IdentifiersBuilder(
            CV.from_pairs(_identifier_pairs(suffix)),
        ),
        annotations=AnnotationsBuilder(
            *_mitab.participant_annotations(
                taxon=f'ncbi_taxid_{suffix}',
                biological_role=f'biological_role_{suffix}',
                experimental_role=f'exp_role_{suffix}',
                identification_method=(
                    f'participant_identification_method_{suffix}'
                ),
            ),
        ),
    )


interactions_schema = RelationBuilder(
    subject=_participant('A'),
    predicate=_predicate,
    object=_participant('B'),
    identifiers=IdentifiersBuilder(
        CV(term=_INNATEDB_INTERACTION, value=_interaction_records),
    ),
    annotations=AnnotationsBuilder(
        *_mitab.interaction_annotations(
            detection_method='interaction_detection_method',
            interaction_type='interaction_type',
            publications='pmid',
            host_organism='ncbi_taxid_host_organism',
        ),
    ),
)


def exclusion_reason(row: dict[str, Any]) -> str | None:
    """Why a row gives no relation; None for a row that is mapped."""
    # A record imported from another database is not InnateDB evidence.
    if _mitab.mi_terms(row.get('source_database')) != [_INNATEDB_SOURCE]:
        return 'source database is not InnateDB'
    # A negative result must never become an interaction.
    if _mitab.is_negative(row.get('negative')):
        return 'negative result'
    for suffix in 'AB':
        if _entity_type(suffix)(row) is None:
            return 'unknown interactor type'
        if not _identifier_pairs(suffix)(row):
            return 'participant without an InnateDB gene record'
    if _predicate(row) is None:
        return 'unknown interaction type'
    return None


def map_interaction(row: dict[str, Any]):
    return None if exclusion_reason(row) else interactions_schema(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url=(
                'https://www.innatedb.com/download/interactions/'
                'innatedb_all.mitab.gz'
            ),
            filename='innatedb_all.mitab.gz',
            subfolder='innatedb',
            large=True,
            ext='gz',
        ),
        mapper=map_interaction,
        raw_parser=iter_interactions,
    ),
)
