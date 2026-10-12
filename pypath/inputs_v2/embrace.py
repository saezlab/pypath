"""
Parse the EMBRACE ligand-receptor pairs and emit Relation records.

Table S11 of Sheikh et al. 2019 is a one-sheet workbook. Each row is a
ligand-receptor pair named by mouse gene symbols, with flags saying in which
of four cell types of the embryonic mouse brain the ligand and the receptor
are expressed.
"""

from __future__ import annotations
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
from pypath.inputs_v2.parsers.embrace import iter_pairs

config = ResourceConfig(
    id=ResourceCv.EMBRACE,
    name='EMBRACE',
    url='https://doi.org/10.1016/j.isci.2019.10.026',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='31677479',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS,),
    summary='Mouse ligand–receptor pairs mapped in the embryonic brain by Sheikh et al. 2019.',
    description='EMBRACE (embryonic brain cell extraction using FACS) is a method to isolate neural, mural, endothelial and microglial cells from the embryonic mouse brain. Sheikh et al. 2019 sequenced the RNA of the four populations at embryonic day 14.5 and mapped the ligand-receptor pairs of Ramilowski et al. 2015 onto them, giving 1,710 mouse ligand-receptor pairs with the cell types that express each ligand and each receptor.',
)
# The journal's own address (www.cell.com/cms/...) refuses scripted
# downloads; this is the same file on the publisher's content server.
download_interactions = Download(
    url='https://ars.els-cdn.com/content/image/1-s2.0-S2589004219304055-mmc2.xlsx',
    filename='embrace_table_s11.xlsx',
    subfolder='embrace',
    default_mode='rb',
    encoding=None,
)
# The table has no species column. The cells are from mouse embryos and the
# symbols are mouse gene symbols.
TAXON = 'NCBITaxon:10090'

f = FieldConfig()


def _participant(role: str) -> EntityBuilder:
    # The symbol is the only identifier the source gives; a row without it
    # has no participant and yields no relation.
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(f'{role}_symbol')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=TAXON),
            # The columns are named for the role (`ligand_*`, `receptor_*`).
            CV(term=PARTICIPANT_ROLE, value=role),
        ),
    )


# The ligand is the subject and the receptor the object. The cell type flags
# are not part of the relation.
interactions_schema = RelationBuilder(
    subject=_participant('ligand'),
    predicate=slots.interacts_with,
    object=_participant('receptor'),
)
resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=iter_pairs,
    ),
)
