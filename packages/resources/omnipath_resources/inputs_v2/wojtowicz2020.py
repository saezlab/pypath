"""
Parse the cell-surface interactome of Wojtowicz et al. 2020 and emit Relation records.

Data S4 of the paper lists the prey-bait pairs that bound in the screen, one
row per tested orientation. The pairs are binding partners without ligand or
receptor roles; the prey is the subject and the bait the object, as the table
orders them.
"""

from __future__ import annotations
import functools
import re
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Relation,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_excel

config = ResourceConfig(
    id=ResourceCv.WOJTOWICZ2020,
    name='Wojtowicz2020',
    full='Wojtowicz 2020 IgSF cell-surface interactome',
    url='https://doi.org/10.1016/j.cell.2020.07.025',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='32822567',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL,),
    summary='Binding pairs of human cell-surface and secreted proteins from an extracellular interaction screen.',
    description='Wojtowicz et al. 2020 screened the extracellular domains of 564 human cell-surface and secreted proteins, most of them members of the immunoglobulin superfamily, for binding in all pairwise combinations with an automated ELISA-based assay. Data S4 lists the prey-bait pairs that bound, with their absorbance readings and whether the interaction had been reported before.',
)
# The article page's attachment link is behind bot protection; this is the same
# file on the publisher's content server.
download = Download(
    url='https://ars.els-cdn.com/content/image/1-s2.0-S0092867420309338-mmc4.xlsx',
    filename='wojtowicz2020_data_s4.xlsx',
    subfolder='wojtowicz2020',
    default_mode='rb',
    encoding=None,
)
TAXON = 'NCBITaxon:9606'
# The columns of Data S4 the mapper reads, for the prey and for the bait.
COLUMNS = tuple(
    f'{side}_{field}'
    for side in ('Prey', 'Bait')
    for field in ('Gene_Name', 'UniProt_Name', 'Common_Name')
)
# The `UniProt_Name` columns hold entry names without their `_HUMAN` ending,
# which for an unreviewed entry is its accession.
_accession_re = re.compile(
    '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}'
)


def _protein(row: dict[str, Any], side: str) -> Entity | None:
    gene = row[f'{side}_Gene_Name']
    if not gene:
        return None
    uniprot_name = row.get(f'{side}_UniProt_Name')
    common_name = row.get(f'{side}_Common_Name')
    identifiers = []
    if uniprot_name and _accession_re.fullmatch(uniprot_name):
        identifiers.append(
            Identifier(type=Namespace.UNIPROT, value=uniprot_name)
        )
    identifiers.append(Identifier(type=Namespace.GENESYMBOL, value=gene))
    # The common name repeats the gene name or the accession where there is none.
    if (
        common_name
        and common_name != gene
        and not _accession_re.fullmatch(common_name)
    ):
        identifiers.append(Identifier(type=Namespace.NAME, value=common_name))
    return Entity(
        type=Protein,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


def map_wojtowicz_interaction(row: dict[str, Any]) -> Relation | None:
    prey = _protein(row, 'Prey')
    bait = _protein(row, 'Bait')
    if prey is None or bait is None:
        return None
    return Relation(
        subject=prey,
        predicate=slots.interacts_with,
        object=bait,
        identifiers=[
            Identifier(
                type=Namespace.NAME,
                value=f'{row["Prey_Gene_Name"]} - {row["Bait_Gene_Name"]}',
            )
        ],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=map_wojtowicz_interaction,
        raw_parser=functools.partial(iter_excel, sheet=0, required=COLUMNS),
    ),
)
