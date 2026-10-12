"""
Parse the INstruct human interactions and emit Relation records.

INstruct combines high-quality binary protein-protein interactions with
structural information from co-crystal structures: for an interaction it gives
the pair of Pfam domains through which the two proteins interact, inferred
from the interfaces of co-crystal structures, and the PDB entries that support
the domain pair. Its table of human interactions (read by
``parsers/instruct.py``) has a row for each protein pair and domain pair, with
the publications of the interaction, the residues where each domain lies in
its protein, and the PDB entries.

This module emits every row as one relation between the two proteins, with the
domain of each protein as a region of it.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.instruct import iter_interactions
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.INSTRUCT,
    name='INstruct',
    url='http://instruct.yulab.org/',
    # The site of INstruct is gone (checked 2026-10-12) and states no terms.
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='23599502',
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Human protein-protein interactions with the Pfam domain pair that '
        'mediates them, where each domain lies in its protein, the '
        'supporting PDB entries and the publications of the interaction.'
    ),
    description=(
        'INstruct is a database of high-quality, three-dimensional, '
        'structurally resolved protein interactome networks of human and six '
        'model organisms. It combines the scale of available high-quality '
        'binary protein interaction data with atomic-resolution structural '
        'information derived from co-crystal evidence, using a tested '
        'interaction interface inference method.'
    ),
)

_HUMAN = 'NCBITaxon:9606'
# Pfam has no member of Namespace yet.
_PFAM = 'pfam'
_INTERACTING_DOMAIN = 'interacting domain'
_PDB_ID = 'instruct:pdb_id'


def _protein_form(suffix: str):
    """The domain of the protein, at each place where it lies in the sequence.

    The residues are counted from 1 on the UniProt sequence; the file does not
    say which isoform, so the sequence is not named. A domain that occurs
    several times has a region for each place: the source does not say which
    of them binds.
    """

    def _form(row: dict[str, Any]) -> dict[str, Any]:
        return {
            'regions': [
                {
                    'type': _INTERACTING_DOMAIN,
                    'identifier': {'ns': _PFAM, 'id': row[f'pfam_{suffix}']},
                    'position': start,
                    'end_position': end,
                    'coordinate_reference': {
                        'identifier': None,
                        'coordinate_system': 'protein',
                        'position_base': 1,
                    },
                    'description': row[f'domain_{suffix}'],
                }
                for start, end in row[f'locations_{suffix}']
            ],
        }

    return _form


def _protein(suffix: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=Protein,
        molecular_form=_protein_form(suffix),
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=lambda row: row[f'uniprot_{suffix}']),
            CV(term=Namespace.ENTREZ, value=lambda row: row[f'entrez_{suffix}']),
            CV(term=Namespace.GENESYMBOL, value=lambda row: row[f'symbol_{suffix}']),
        ),
        # The file is the one for Homo sapiens (sapiens.sin).
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_HUMAN)),
    )


def _publications(row: dict[str, Any]) -> list[str]:
    """PubMed IDs only: the column also holds text such as a consortium name."""
    return [curie for curie in map(pmid_curie, row['publications']) if curie]


# The file has no interaction type. The publication describes the pairs as
# high-quality binary interactions and the domain pair as inferred from
# co-crystal evidence; the file does not say whether a listed structure holds
# these two proteins, so only the broad interaction is stated.
interactions_schema = RelationBuilder(
    subject=_protein('a'),
    predicate=slots.interacts_with,
    object=_protein('b'),
    annotations=AnnotationsBuilder(
        CV(term=slots.publications, value=_publications),
        CV(term=_PDB_ID, value=lambda row: row['pdb_ids']),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            # The original, http://instruct.yulab.org/download/sapiens.sin,
            # no longer exists. This is the copy on the OmniPath rescued-data
            # server, which the legacy reader used as well.
            url='https://rescued.omnipathdb.org/instruct-sapiens.sin',
            filename='instruct-sapiens.sin',
            subfolder='instruct',
            large=True,
        ),
        mapper=interactions_schema,
        raw_parser=iter_interactions,
    ),
)
