"""
Parse PDZBase interactions and emit Relation records.

PDZBase lists protein-protein interactions in which a PDZ domain of one
protein binds the other protein, its ligand. Each record names the two
proteins by UniProt accession, says which PDZ domain of the first protein
binds, and cites one publication. The module emits the interaction between
the two proteins. The number of the PDZ domain and the organism label stay
in the parsed record.

The PDZBase website is offline. This module reads the copy of its "all
interactions" page that OmniPath keeps in its archive of rescued files.
"""

from __future__ import annotations

from collections.abc import Callable
import re
from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.pdzbase import iter_pdzbase


config = ResourceConfig(
    id=ResourceCv.PDZBASE,
    name='PDZBase',
    url='http://abc.med.cornell.edu/pdzbase',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='15513994',
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Protein-protein interactions mediated by a PDZ domain, with their '
        'publications.'
    ),
    description=(
        'PDZBase is a manually curated database of protein-protein '
        'interactions mediated by PDZ domains. Its curators extracted about '
        '340 interactions from the literature, taking only interactions shown '
        'by in vivo or in vitro experiments such as co-immunoprecipitation '
        'and pull-down, mediated directly by the PDZ domain and with the '
        'binding sites known on both proteins. Each record names the protein '
        'that has the PDZ domain, the domain that binds and the ligand '
        'protein. The database dates from 2004 and its website is no longer '
        'online.'
    ),
)

# A UniProt accession, with the isoform number some records give.
_ACCESSION_RE = re.compile(
    r'^((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(?:-\d+)?)$'
)

# PDZBase's numbers of its interaction records, an identifier system without
# a registered namespace.
_PDZBASE_INTERACTION = 'pdzbase'

f = FieldConfig(
    extract={
        'accession': _ACCESSION_RE,
        'pmid': PMID_EXTRACT,
        'record': r'^(\d+)$',
    },
    transform={'pmid': lambda v: f'PMID:{v}'},
)


def _protein(role: str) -> Callable[[dict[str, Any]], Any]:
    # No taxon: the organism label of a record is not the species of both
    # proteins (a rat PDZ protein with a human or a viral ligand is labelled
    # ``human``).
    builder = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.UNIPROT,
                value=f(f'{role}_accession', extract='accession'),
            ),
            CV(term=Namespace.NAME, value=f(f'{role}_name')),
        ),
    )

    def build(row: dict[str, Any]) -> Any:
        # The accession identifies the protein; PDZBase's name alone does not.
        accession = str(row.get(f'{role}_accession') or '').strip()
        return builder.build(row) if _ACCESSION_RE.match(accession) else None

    return build


# The record says that the PDZ domain of the first protein mediates the
# interaction with the ligand. The page has no interaction type or method per
# record; that the curators took only interactions mediated directly by the
# domain is their description of the database, not a statement of the record,
# so the relation is an interaction and not a direct one. The PDZ protein is
# the subject and the ligand the object.
interactions_schema = RelationBuilder(
    subject=_protein('pdz'),
    predicate=slots.interacts_with,
    object=_protein('ligand'),
    identifiers=IdentifiersBuilder(
        CV(
            term=_PDZBASE_INTERACTION,
            value=f('interaction_id', extract='record'),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.publications,
            value=f('reference', extract='pmid', transform='pmid'),
        ),
    ),
)


resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/PDZBase.html',
            filename='PDZBase.html',
            subfolder='pdzbase',
            large=True,
        ),
        mapper=interactions_schema,
        raw_parser=iter_pdzbase,
    ),
)
