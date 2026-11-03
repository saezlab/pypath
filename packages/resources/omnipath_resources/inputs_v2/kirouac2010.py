"""
Parse the ligand-receptor mapping of Kirouac et al. 2010 and emit Relation records.

Supplementary Table S12 of the paper has one row per ligand: a free-text
ligand cell, usually `description (SYMBOL / alias)`, and a receptor cell in
the authors' shorthand, where `BMPR 1A/1B/2` stands for BMPR1A, BMPR1B and
BMPR2 and `FZD 1-10` for FZD1 to FZD10. The parser expands the shorthand into
one ligand-receptor pair per listed receptor; symbols keep the table's spelling.
A receptor cell with a second designation in brackets (`IL8RB (CXCR2)`) gives it as
a synonym of the receptor.
"""

from __future__ import annotations
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from omnipath_core.source_attributes import PARTICIPANT_ROLE
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Relation,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.kirouac2010 import iter_pairs

config = ResourceConfig(
    id=ResourceCv.KIROUAC2010,
    name='Kirouac2010',
    full='Kirouac 2010 ligand-receptor mapping',
    url='https://doi.org/10.1038/msb.2010.71',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_NC_ND_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='20924352',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Manually curated pairs of secreted ligands and their receptors.',
    description='Kirouac et al. 2010 reconstructed the intercellular signalling networks of human blood stem and progenitor cell cultures from gene expression. For this the authors manually curated a mapping of secreted ligands, mostly cytokines, chemokines and growth factors, to their receptors, largely based on the COPE cytokine database, published as Supplementary Table S12.',
)
# The publisher's download page is behind bot protection; Europe PMC serves the
# article's supplementary files as one archive.
download = Download(
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2990637/supplementaryFiles?includeInlineImage=false',
    filename='kirouac2010_supplementary_files.zip',
    subfolder='kirouac2010',
    default_mode='rb',
    encoding=None,
    ext='zip',
    needed=['msb201071-s3.xls'],
)
TAXON = 'NCBITaxon:9606'


def _protein(symbol: str, role: str, synonym: str | None = None) -> Entity:
    identifiers = [Identifier(type=Namespace.GENESYMBOL, value=symbol)]
    if synonym:
        identifiers.append(Identifier(type=Namespace.SYNONYM, value=synonym))
    return Entity(
        type=Protein,
        identifiers=identifiers,
        annotations=[
            Annotation(term=slots.in_taxon, value=TAXON),
            # The table has a ligand column and a receptor column.
            Annotation(term=PARTICIPANT_ROLE, value=role),
        ],
    )


def map_kirouac_interaction(row: dict[str, Any]) -> Relation | None:
    ligand = row.get('ligand')
    receptor = row.get('receptor')
    if not ligand or not receptor:
        return None
    return Relation(
        subject=_protein(ligand, 'ligand'),
        predicate=slots.interacts_with,
        object=_protein(receptor, 'receptor', row.get('receptor_synonym')),
        identifiers=[
            Identifier(type=Namespace.NAME, value=f'{ligand} - {receptor}')
        ],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=map_kirouac_interaction,
        raw_parser=iter_pairs,
    ),
)
