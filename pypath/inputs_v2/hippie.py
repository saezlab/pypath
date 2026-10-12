"""
Parse the HIPPIE human protein-protein interactions and emit Relation records.

HIPPIE (Human Integrated Protein-Protein Interaction rEference) integrates the
human interactions of several interaction databases and studies and gives each
pair a confidence score, computed from the number of studies, the quality of
the detection techniques and the number of other organisms in which the
interaction was found. Its table (read by ``parsers/hippie.py``) has one row
for a pair of UniProt entries, with the Entrez Gene ID of each, the score, and
the techniques, publications and source databases of the pair.

This module emits every row as one relation between the two proteins. It does
not apply the score threshold that the legacy reader applied: the score is the
confidence score of the relation and a user can filter on it.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.hippie import iter_hippie
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.HIPPIE,
    name='HIPPIE',
    full='Human Integrated Protein-Protein Interaction rEference',
    url='https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/',
    # Neither the download page nor the information page of HIPPIE states
    # terms of use (checked 2026-10-12).
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='27794551',
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    summary=(
        'Human protein-protein interactions from several interaction '
        'databases and studies, with the HIPPIE confidence score, '
        'publications and source databases of each pair (release 2.4).'
    ),
    description=(
        'HIPPIE integrates the human protein-protein interactions of ten '
        'source databases and eleven studies and scores each pair by the '
        'amount and reliability of its evidence: the number of studies that '
        'detected it, the number and quality of the experimental techniques, '
        'and the number of non-human organisms in which it was reproduced. '
        'The quality of each technique was assigned by a group of experts.'
    ),
)

_HUMAN = 'NCBITaxon:9606'
# The scale of the confidence score, which the source gives without a unit.
_SCORE_FIELD = 'HIPPIE confidence score (0 to 1)'
# Other organisms in which the interaction was found, named as the source does.
_CONSERVED_IN = 'hippie:conserved_in_species'

# The source databases and studies the table names, as it writes them. They
# are carried under these labels (as ``lrdb.py`` and ``talklr.py`` carry
# theirs); a name that is not listed raises. "MIPS" is the table's own label:
# it does not say which MIPS database is meant.
SOURCES = frozenset(
    {
        'BioGRID', 'DIP', 'HPRD', 'IntAct', 'STRING', 'Albers05', 'Bell09',
        'BIND', 'BioCarta', 'Colland04', 'Goehler04', 'HomoMINT', 'I2D',
        'Kaltenbach07', 'KEGG', 'Lehner04', 'Lim06', 'MINT', 'MIPS',
        'Nakayama02', 'Ortho', 'Rual05', 'Stelzl05', 'Venkatesan09',
    }
)


def _protein(suffix: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT_ENTRY, value=lambda row: row[f'entry_names_{suffix}']),
            CV(term=Namespace.ENTREZ, value=lambda row: row[f'entrez_{suffix}']),
        ),
        # The table is the human interactome: every entry name ends in _HUMAN.
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_HUMAN)),
    )


def _publications(row: dict[str, Any]) -> list[str]:
    """PubMed IDs only: the table also holds 2147483647, a placeholder."""
    return [curie for curie in map(pmid_curie, row['pmids']) if curie]


def _sources(row: dict[str, Any]) -> list[str]:
    """The source databases and studies of the row, under the table's labels."""
    for source in row['sources']:
        if source not in SOURCES:
            raise ValueError(f'Unknown HIPPIE source database: {source!r}')
    return list(row['sources'])


# The table mixes physical and functional evidence (the STRING, KEGG and
# BioCarta pairs, interactions inferred through orthologs, co-localization
# techniques), so no row is a direct physical interaction by itself.
interactions_schema = RelationBuilder(
    subject=_protein('a'),
    predicate=slots.interacts_with,
    object=_protein('b'),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.has_confidence_score,
            value=lambda row: measurement(row['score'], source_field=_SCORE_FIELD),
        ),
        CV(term=slots.publications, value=_publications),
        CV(term=slots.supporting_data_source, value=_sources),
        CV(term=_CONSERVED_IN, value=lambda row: row['species']),
    ),
)


resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt',
            filename='hippie_current.txt',
            subfolder='hippie',
            large=True,
        ),
        mapper=interactions_schema,
        raw_parser=iter_hippie,
    ),
)
