"""
Parse the Interactome3D human interactions and emit Relation records.

Interactome3D annotates protein-protein interactions with three-dimensional
structural data. Its table of the human representative set (read by
``parsers/i3d.py``) has one row for each interaction: the best ranked
experimental structure of the two proteins together, or, where there is none,
a homology model of the complex on a template of the PDB, or, where there is no
such template either, a model on a domain-domain template of 3did. A row also
gives, for each protein, the chain, the percentage of sequence identity to the
structure or template, the percentage of the protein it covers and the residues
of the UniProt sequence it covers.

This module emits every row as one relation between the two proteins. The
three row types do not have the same standing and get a different predicate
and confidence level (see ``TYPES``). The residues a structure or model covers
are a region of each protein.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.i3d import iter_interactions
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.I3D,
    name='Interactome3D',
    url='https://interactome3d.irbbarcelona.org/',
    # The site states no terms of use (home, download, help, references and
    # about pages, checked 2026-10-12).
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='23399932',
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Human protein-protein interactions with the experimental structure '
        'or homology model of the complex, its PDB entry and the residues of '
        'each protein it covers (representative set, release 2024_12).'
    ),
    description=(
        'Interactome3D is a web service for the structural annotation of '
        'protein-protein interaction networks. Its pre-built datasets are '
        'made of binary interactions from public interaction databases and '
        'the structural data for them: experimental structures of the '
        'complex from the PDB, homology models of the complex built with '
        'Modeller on structural templates of the PDB, and models on '
        'domain-domain templates of 3did. For an interaction only one of '
        'these kinds is given, the first available in that order.'
    ),
)

# Source TYPE -> (predicate, confidence level). The level repeats the words of
# the help of Interactome3D.
# - Structure: an experimental structure of the two proteins together, so the
#   chains touch: a direct physical interaction.
# - Model, Dom_dom_model: a model built on the structure of other proteins. It
#   shows how the two proteins could bind, not that they do, so it supports
#   only the broad interaction. The help says that the pairs of the pre-built
#   datasets come from public interaction databases; the table does not name
#   them.
TYPES = {
    'Structure': (
        slots.directly_physically_interacts_with,
        'experimental structure',
    ),
    'Model': (
        slots.interacts_with,
        'homology model from a global template',
    ),
    'Dom_dom_model': (
        slots.interacts_with,
        'model from a domain-domain template',
    ),
}
_COVERED_REGION = {
    'Structure': 'region covered by the structure',
    'Model': 'region covered by the model',
    'Dom_dom_model': 'region covered by the model',
}
_HUMAN = 'NCBITaxon:9606'
_PDB_ID = 'i3d:pdb_id'
_CHAIN = 'i3d:chain'
_SEQUENCE_IDENTITY = 'i3d:sequence_identity'
_COVERAGE = 'i3d:coverage'


def _type(row: dict[str, Any]) -> str:
    if row['type'] not in TYPES:
        raise ValueError(f'Unknown Interactome3D structure type: {row["type"]!r}')
    return row['type']


def _predicate(row: dict[str, Any]):
    return TYPES[_type(row)][0]


def _region(suffix: str):
    """The residues of the UniProt sequence that the structure or model covers.

    The residues are counted from 1 on the UniProt sequence, which the help
    says is the canonical isoform; the sequence is not named. A domain-domain
    model also gives the Pfam name of the modeled domain.
    """

    def _form(row: dict[str, Any]) -> dict[str, Any]:
        return {
            'regions': [
                {
                    'type': _COVERED_REGION[_type(row)],
                    'position': row[f'seq_begin_{suffix}'],
                    'end_position': row[f'seq_end_{suffix}'],
                    'coordinate_reference': {
                        'identifier': None,
                        'coordinate_system': 'protein',
                        'position_base': 1,
                    },
                    'description': row[f'domain_{suffix}'],
                }
            ],
        }

    return _form


def _protein(suffix: str) -> EntityBuilder:
    number = '1' if suffix == 'a' else '2'
    return EntityBuilder(
        entity_type=Protein,
        molecular_form=_region(suffix),
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=lambda row: row[f'uniprot_{suffix}']),
        ),
        # The dataset is the human one.
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=_HUMAN),
            CV(term=_CHAIN, value=lambda row: row[f'chain_{suffix}']),
            CV(
                term=_SEQUENCE_IDENTITY,
                value=lambda row: measurement(
                    row[f'seq_identity_{suffix}'],
                    unit='%',
                    source_field=f'SEQ_IDENT{number}',
                ),
            ),
            CV(
                term=_COVERAGE,
                value=lambda row: measurement(
                    row[f'coverage_{suffix}'],
                    unit='%',
                    source_field=f'COVERAGE{number}',
                ),
            ),
        ),
    )


interactions_schema = RelationBuilder(
    subject=_protein('a'),
    predicate=_predicate,
    object=_protein('b'),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_confidence_level, value=lambda row: TYPES[_type(row)][1]),
        # The PDB entry of the structure, or of the template of a model.
        CV(term=_PDB_ID, value=lambda row: row['pdb_id']),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url=(
                'https://interactome3d.irbbarcelona.org/user_data/human/'
                'download/representative/interactions.dat'
            ),
            filename='interactome3d_human_representative_interactions.dat',
            subfolder='i3d',
            large=True,
        ),
        mapper=interactions_schema,
        raw_parser=iter_interactions,
    ),
)
