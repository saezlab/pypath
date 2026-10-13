"""
Marker genes of the cell types of the CellTypist immune cell encyclopedia.

CellTypist is a tool for automated cell type annotation of single cell RNA
sequencing data. Its encyclopedia describes the 98 cell types of the Pan-immune
CellTypist atlas (version 2) in one table, each with its Cell Ontology term, the
tissues and datasets it was found in, curated marker genes and the ten genes
with the largest weight in the CellTypist model. The module imports one
dataset, ``markers``: a statement for each gene of the two gene columns, with
the Cell Ontology term of the cell type as its object.

Every cell type of the table has a Cell Ontology ID, so the object of a
statement is the Cell Ontology's own term entity (``Namespace.CL``, as ``cl.py``
writes it) and the module has no vocabulary of its own. The labels of the table
(``Follicular B cells``, and the group ``B cells``) are carried on the
statement, because several cell types of the table share a Cell Ontology term
(for example four cell types are ``CL:0001056``, dendritic cell).
"""

from __future__ import annotations

from biolink_model.datamodel.model import Gene, OntologyClass, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.celltypist import GENE_COLUMNS, iter_markers
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Relation
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.CELLTYPIST,
    name='CellTypist',
    full='CellTypist cell type encyclopedia',
    url='https://www.celltypist.org/',
    # Every page of celltypist.org, the encyclopedia page included, links "Read our Terms &
    # Conditions" to https://www.sanger.ac.uk/legal/ (checked 2026-10-12), which states:
    # "Unless otherwise stated, all online content created and hosted by the Wellcome Sanger
    # Institute is licensed under a Creative Commons Attribution-NonCommercial NoDerivatives
    # 4.0 International License (CC BY-NC-ND 4.0)". The table itself is a file of the GitHub
    # repository Teichlab/celltypist_wiki, which has no license file and states none; the
    # encyclopedia page says its information "comes from our encyclopedia table". The legacy
    # registry (resources.json) records CC BY 3.0 for CellTypist; the site's statement of
    # today is used instead.
    license=LicenseCV.CC_BY_NC_ND_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='35549406',
    topics=(Topic.CELL_TYPES,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED, KnowledgeOrigin.COMPUTATIONALLY_PREDICTED),
    summary='Curated marker genes and model top genes of immune cell types, with Cell Ontology terms.',
    description=(
        'CellTypist is a tool for the automated annotation of cell types in '
        'single cell RNA sequencing data. Its cell type encyclopedia lists the '
        'immune cell types of the Pan-immune CellTypist atlas with a '
        'description, the Cell Ontology term, the tissues and datasets the cell '
        'type was found in, curated marker genes and the ten genes with the '
        'largest weight for the cell type in the CellTypist model.'
    ),
)

# ================================== DOWNLOAD ==================================

# The legacy URL is unchanged. The file holds one worksheet, ``Main``.
download = Download(
    url=(
        'https://raw.githubusercontent.com/Teichlab/celltypist_wiki/main/atlases/'
        'Pan_Immune_CellTypist/v2/encyclopedia/encyclopedia_table.xlsx'
    ),
    filename='celltypist_encyclopedia_table_Pan_Immune_v2.xlsx',
    subfolder='celltypist',
    large=False,
    default_mode='rb',
    encoding=None,
)

# =================================== SCHEMA ===================================

# The table has no species column. Its atlas is the human cross-tissue immune
# atlas (celltypist.org cites "Cross-tissue immune cell analysis reveals
# tissue-specific features in humans"), and its genes are HGNC symbols.
TAXON = 'NCBITaxon:9606'

# What the table's gene columns say. A curated marker is a gene chosen by the
# curators as marking the cell type; ``expressed_in`` is used for it as for the
# markers of other resources. The model genes are the ten genes with the
# largest weight for the cell type in a regression model: they are the genes
# that tell the cell type from the others in that model, and many are not
# markers of the cell type (``MT-RNR2``, ``HBG2`` for B cells), so they are
# only associated with it.
PREDICATES = {
    'curated_marker': slots.expressed_in,
    'model_top_gene': slots.associated_with,
}
assert set(PREDICATES) == set(GENE_COLUMNS)

SOURCE = 'celltypist:marker_source'
CELL_TYPE = 'celltypist:cell_type'
CELL_TYPE_GROUP = 'celltypist:cell_type_group'
CELL_TYPE_DESCRIPTION = 'celltypist:cell_type_description'
TISSUE = 'celltypist:tissue'
DATASET = 'celltypist:dataset'

f = FieldConfig()

marker_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=Gene,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('gene')),
        ),
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=TAXON)),
    ),
    predicate=lambda row: PREDICATES[row['marker_source']],
    object=EntityBuilder(
        entity_type=OntologyClass,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.CL, value=f('cell_ontology_id')),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(term=SOURCE, value=f('marker_source')),
        # The labels of the table's cell type: the Cell Ontology term can be
        # shared by several of them.
        CV(term=CELL_TYPE, value=f('cell_type')),
        CV(term=CELL_TYPE_GROUP, value=f('cell_type_group')),
        CV(term=CELL_TYPE_DESCRIPTION, value=f('description')),
        # Where the cell type was found, as the table names them: tissue
        # names and "author year" labels of datasets, not ontology terms.
        CV(term=TISSUE, value=f('tissues')),
        CV(term=DATASET, value=f('datasets')),
    ),
)


def map_marker(row: dict) -> Relation | None:
    """A marker statement; an unknown marker source raises."""
    if row.get('marker_source') not in PREDICATES:
        raise ValueError(f'Unknown CellTypist marker source: {row.get("marker_source")!r}')
    return marker_schema(row)


# ================================= RESOURCE ===================================

resource = Resource(
    config,
    markers=Dataset(
        download=download,
        mapper=map_marker,
        raw_parser=iter_markers,
    ),
)

# ================================= REFERENCE ==================================

# High-hierarchy cell types  Low-hierarchy cell types  Description  Cell Ontology ID
# Tissues  Datasets  Curated markers  Top 10 important genes from the CellTypist model
# e.g. B cells  Follicular B cells  resting mature B lymphocytes ...  CL:0000843
#      Intestine  Smillie et al. 2019  CXCR5, TNFRSF13B, CD22  MT-RNR2, IGHA1, JCHAIN, ...
