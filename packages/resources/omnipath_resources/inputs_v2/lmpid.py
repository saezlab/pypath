"""
Parse LMPID linear motif mediated interactions and emit Relation records.

Each LMPID record is a linear motif instance in one protein (the bait) and a
protein whose domain binds it (the prey), with the publications that report
the interaction. The module emits the protein-protein interaction, with the
motif instance (its residues and its position in the bait) as a region of the
bait and the binding domain, where the source names it, as a region of the
prey.

The LMPID website offers XML and spreadsheet exports without literature
references. This module reads the export with references that the LMPID
authors provided to OmniPath, from OmniPath's archive of such files.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.lmpid import iter_lmpid


config = ResourceConfig(
    id=ResourceCv.LMPID,
    name='LMPID',
    full='Linear Motif mediated Protein Interaction Database',
    url='https://bicresources.jcbose.ac.in/ssaha4/lmpid/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_NC_2_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='25776024',
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Protein-protein interactions mediated by a linear motif binding '
        'a domain, with the motif instance and its position, the domain '
        'and the publications.'
    ),
    description=(
        'LMPID (Linear Motif mediated Protein Interaction Database) is a '
        'manually curated collection of experimentally validated linear '
        'motifs that mediate protein-protein interactions, in proteins of '
        'any organism. Each entry gives the motif instance and its position '
        'in the protein that carries it, the interacting protein and, where '
        'known, the domain that binds the motif.'
    ),
)

_ACCESSION_RE = re.compile(
    r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$'
)

f = FieldConfig(
    extract={'accession': _ACCESSION_RE},
    transform={'pmid': lambda v: f'PMID:{v}'},
)

# The source's words for the two regions: the motif instance of the bait and
# the domain of the prey that binds it (``Interacting_Domain``).
_MOTIF_REGION = 'linear motif'
_DOMAIN_REGION = 'interacting domain'


def _motif_form(row: dict[str, Any]) -> dict[str, Any] | None:
    """The motif instance: its residues, and where they lie in the bait.

    The positions are the source's ``start-end`` text; of the 2,203 records
    all but ten have an end that makes the range as long as the residues, so
    they count residues from 1 and include both ends. The export does not say
    which sequence of the bait they count, so none is named.
    """
    start, end = row.get('motif_start'), row.get('motif_end')
    if start is None or end is None or end < start:
        start = end = None
    if start is None and not row.get('motif_instance'):
        return None
    return {
        'regions': [
            {
                'type': _MOTIF_REGION,
                'position': start,
                'end_position': end,
                'coordinate_reference': {
                    'identifier': None,
                    'coordinate_system': 'protein',
                    'position_base': 1,
                }
                if start is not None
                else None,
                'description': row.get('motif_instance'),
            }
        ],
    }


def _domain_form(row: dict[str, Any]) -> dict[str, Any] | None:
    """The domain of the prey that binds the motif, by the name the source gives.

    The source gives no accession and no position for it.
    """
    domain = row.get('interacting_domain')
    if not domain:
        return None
    return {'regions': [{'type': _DOMAIN_REGION, 'description': domain}]}


def _protein(column: str, form) -> EntityBuilder:
    # The export names no organism, so none is asserted.
    return EntityBuilder(
        entity_type=Protein,
        molecular_form=form,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=f(column, extract='accession')),
        ),
    )


# The export states that the motif of the bait mediates an interaction with the
# prey (with the domain that binds it); it has no interaction type or method
# per record, and the LMPID help lists co-immunoprecipitation and affinity
# purification among the experiments of the database. A record is therefore an
# interaction, not a stated direct one.
interactions_schema = RelationBuilder(
    subject=_protein('bait_uniprot', _motif_form),
    predicate=slots.interacts_with,
    object=_protein('prey_uniprot', _domain_form),
    annotations=AnnotationsBuilder(
        CV(term=slots.publications, value=f('pmids', transform='pmid')),
    ),
)


resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/LMPID_DATA_pubmed_ref.xml',
            filename='LMPID_DATA_pubmed_ref.xml',
            subfolder='lmpid',
            large=True,
            encoding=None,
            default_mode='rb',
        ),
        mapper=interactions_schema,
        raw_parser=iter_lmpid,
    ),
)
