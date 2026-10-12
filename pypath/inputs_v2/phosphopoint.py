"""
Parse the PhosphoPOINT human kinase interactome and emit Relation records.

The table pairs each kinase with a protein and puts the pair into one of four
categories (semicolon separated, with a header):

    Kinase;Entrez GeneID;Gene;Entrez GeneID;Category
    AKT1;207;ESR2;2100;Category 2

The categories are the "four kinds of links" of PhosphoPOINT between a kinase
and another protein (PMID 18689816): physically interacting proteins
(Category 1), interacting phospho-proteins (Category 2), substrates (Category
3) and substrates as well as interacting phospho-proteins (Category 4). The
table has no phosphorylation sites and no references.
"""

from __future__ import annotations
from collections.abc import Callable
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    DirectionQualifierEnum,
    GeneOrGeneProductOrChemicalEntityAspectEnum,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phosphopoint import iter_pairs

config = ResourceConfig(
    id=ResourceCv.PHOSPHOPOINT,
    name='PhosphoPOINT',
    url='https://doi.org/10.1093/bioinformatics/btn297',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_SA_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='18689816',
    topics=(Topic.SIGNALING,),
    content=(Content.ENZYME_PTM_RELATIONSHIPS, Content.PHYSICAL_INTERACTIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Human kinases with their substrates and their interacting proteins, '
        'with the category of each link.'
    ),
    description='PhosphoPOINT is a human kinase interactome and phospho-protein database. It combined the phosphorylation sites of Phospho.ELM, HPRD and Swiss-Prot and 400 kinase-substrate pairs curated by its authors with the protein-protein interactions of the POINT database, and sorted the proteins linked to each of 441 kinases into interacting proteins, interacting phospho-proteins, substrates, and substrates that also interact with the kinase.',
)
# The database's site (kinase.bioinformatics.tw) no longer exists and the
# Internet Archive did not keep its download; this is the copy OmniPath kept.
download = Download(
    url='https://rescued.omnipathdb.org/phosphopoint.csv',
    filename='phosphopoint.csv',
    subfolder='phosphopoint',
    default_mode='rb',
    encoding=None,
)
# PhosphoPOINT is a database of human kinases; the Gene IDs are human.
TAXON = 'NCBITaxon:9606'
# Category -> whether the protein is a substrate of the kinase and whether the
# two interact physically. A phospho-protein that interacts with a kinase
# (Category 2) is not thereby its substrate.
CATEGORIES = {
    'Category 1': {'substrate': False, 'interaction': True},
    'Category 2': {'substrate': False, 'interaction': True},
    'Category 3': {'substrate': True, 'interaction': False},
    'Category 4': {'substrate': True, 'interaction': True},
}


# The category of the pair, as the table writes it ("Category 1"), is an
# attribute of the relation. The categories are kinds of links, not confidence
# classes. PhosphoPOINT's introduction (web.archive.org capture of
# kinase.bioinformatics.tw/introduction.html, April 2012) says it aims "to
# delineate four kinds of links among kinases": their interacting proteins
# (Category 1), substrates (Category 3), and substrates as well as interacting
# phospho-proteins (Category 4); "some of these interacting proteins for
# kinases are phospho-proteins (...Category 2), which might have the
# potential to serve as substrates". Category 1 and 2 are both interactions;
# 2 marks a phospho-protein partner.
_CATEGORY = 'phosphopoint:category'


def _category(row: dict[str, Any]) -> dict[str, bool]:
    category = row.get('Category')
    if category not in CATEGORIES:
        raise ValueError(f'Unknown PhosphoPOINT category: {category!r}')
    return CATEGORIES[category]


f = FieldConfig(
    extract={
        'gene_id': r'^(\d+)$',
        # A spreadsheet turned three symbols into dates (``01-Sep``); those
        # proteins are identified by their Gene ID only.
        'symbol': r'^(?!\d{1,2}-[A-Z][a-z]{2}$)(\S+)$',
    },
)


def _protein(role: str) -> EntityBuilder:
    # A row without a Gene ID or symbol yields no participant.
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.ENTREZ,
                value=f(f'{role} Entrez GeneID', extract='gene_id'),
            ),
            CV(term=Namespace.GENESYMBOL, value=f(role, extract='symbol')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=TAXON),
        ),
    )


# The kinase increases the phosphorylation of its substrate; the table says
# nothing about the substrate's activity and names no site.
enzyme_substrate_schema = RelationBuilder(
    subject=_protein('Kinase'),
    predicate=slots.affects,
    object=_protein('Gene'),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.object_direction_qualifier,
            value=DirectionQualifierEnum.increased,
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=GeneOrGeneProductOrChemicalEntityAspectEnum.phosphorylation,
        ),
        CV(term=_CATEGORY, value=f('Category')),
    ),
)
interactions_schema = RelationBuilder(
    subject=_protein('Kinase'),
    predicate=slots.interacts_with,
    object=_protein('Gene'),
    annotations=AnnotationsBuilder(
        CV(term=_CATEGORY, value=f('Category')),
    ),
)


def _only(kind: str, schema: RelationBuilder) -> Callable[[dict[str, Any]], Any]:
    def mapper(row: dict[str, Any]) -> Any:
        return schema.build(row) if _category(row)[kind] else None

    return mapper


# A Category 4 pair is in both datasets: the source states both that the
# protein is a substrate of the kinase and that the two interact.
resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=download,
        mapper=_only('substrate', enzyme_substrate_schema),
        raw_parser=iter_pairs,
    ),
    interactions=Dataset(
        download=download,
        mapper=_only('interaction', interactions_schema),
        raw_parser=iter_pairs,
    ),
)
