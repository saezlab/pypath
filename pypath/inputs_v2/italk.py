"""
Parse the iTALK ligand-receptor database and emit Relation records.

The iTALK R package ships one data frame of human ligand-receptor pairs. Each
row names the pair and gives the approved gene symbol and a name for the
ligand and for the receptor.
"""

from __future__ import annotations
from collections.abc import Callable
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, Topic
from omnipath_core.source_attributes import PARTICIPANT_ROLE
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.italk import iter_pairs

config = ResourceConfig(
    id=ResourceCv.ITALK,
    name='iTALK',
    url='https://github.com/Coolgenome/iTALK',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_NC_SA_3_0,
    update_category=UpdateCategoryCV.STATIC,
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS,),
    summary='Human ligand–receptor pairs from the iTALK package.',
    description='iTALK is an R package for characterizing and illustrating intercellular communication from bulk or single-cell RNA sequencing data. Its built-in database lists known human ligand-receptor pairs, each classified as a cytokine, growth factor, immune checkpoint or other pair.',
)
download_interactions = Download(
    url='https://raw.githubusercontent.com/Coolgenome/iTALK/master/data/LR_database.rda',
    filename='italk_LR_database.rda',
    subfolder='italk',
    default_mode='rb',
    encoding=None,
)
# The table has no species column. The symbols are HGNC approved symbols.
TAXON = 'NCBITaxon:9606'

f = FieldConfig()


def _participant(role: str) -> Callable[[dict[str, Any]], Any]:
    symbol = f'{role}.ApprovedSymbol'
    builder = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(symbol)),
            CV(term=Namespace.NAME, value=f(f'{role}.Name')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=TAXON),
            # The columns are named for the role (`Ligand.*`, `Receptor.*`).
            CV(term=PARTICIPANT_ROLE, value=role.lower()),
        ),
    )

    def build(row: dict[str, Any]) -> Any:
        # Seventeen receptors are viral proteins named only by an alias
        # (``US28``, ``KSHV``). Without an approved symbol the participant and
        # its species are not identified, so the pair is not emitted.
        return builder.build(row) if (row.get(symbol) or '').strip() else None

    return build


# The ligand is the subject and the receptor the object.
interactions_schema = RelationBuilder(
    subject=_participant('Ligand'),
    predicate=slots.interacts_with,
    object=_participant('Receptor'),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('Pair.Name')),
    ),
)
resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=iter_pairs,
    ),
)
