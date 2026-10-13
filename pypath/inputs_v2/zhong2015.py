"""
Parse the human cell adhesion molecules of Zhong et al. 2015 and emit gene-to-type records.

Zhong et al. updated the set of human genes likely to encode cell adhesion
molecules (CAMs) and divided them into functional subclasses: CAMs with a
primary function in information exchange (iCAMs), CAMs of matrix adhesion and
focal adhesion, of cell-cell adhesion, of adherens junctions and tight
junctions, of axonal guidance and of myelin interactions. Supplementary Table S1
lists the genes with their subclass code, and whether the searches of InterPro,
Gene Ontology and NCBI found them.

The subclasses are terms of their own (namespace ``zhong2015_cam_type``, the
code as the table writes it, with the name of the subclass; the ``terms``
dataset); each gene is related to the type the table gives it.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Gene, OntologyClass, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.zhong2015 import iter_cam_types, iter_genes
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Relation

config = ResourceConfig(
    id=ResourceCv.ZHONG2015,
    name='Zhong2015',
    url='https://doi.org/10.1111/nyas.12776',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    # Europe PMC lists the article as not open access (an author manuscript, PMC4564344;
    # checked 2026-10-13), so nothing from the source supports the registry's entry.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='25988664',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Human cell adhesion molecule genes with their functional subtype.',
    description=(
        'Zhong et al. updated the set of human genes likely to encode cell '
        'adhesion molecules (CAMs) with searches of databases, literature '
        'reviews and annotations, and described functional subclasses: CAMs '
        'with a primary function in information exchange (iCAMs), CAMs of '
        'focal adhesions, of adherens junctions and other connections between '
        'cells, and smaller classes.'
    ),
)

# Supplementary Table S1 of the paper: OmniPath's copy.
download = Download(
    url='https://rescued.omnipathdb.org/zhong2015_S1.tsv',
    filename='zhong2015_S1.tsv',
    subfolder='zhong2015',
    default_mode='r',
)

TYPE_NAMESPACE = 'zhong2015_cam_type'
TAXON = 'NCBITaxon:9606'

# Code of the CAM type -> its name. The table gives the codes only; the names are
# those the legacy OmniPath reader used. The paper's abstract names the iCAMs, focal
# adhesions and adherens junctions; the other names could not be checked against
# the table's legend, which is not in the file.
CAM_TYPES = {
    'i': 'iCAM',
    'm': 'matrix adhesion',
    'ag': 'axonal guidance',
    'aj': 'adherens junction',
    'c': 'cell-cell adhesion',
    'fa': 'focal adhesion',
    'tj': 'tight junction',
    'my': 'myelin interactions',
}

# Source attributes of the statement: whether the search of a database found the
# gene (the table's 1/0 columns, written 'true'/'false'), and in how many of the
# three (the table's Sum; a Measurement, a count without a unit).
_FOUND_BY = {
    'interpro': 'zhong2015:found_by_interpro',
    'go': 'zhong2015:found_by_go',
    'ncbi': 'zhong2015:found_by_ncbi',
}
_SOURCES = 'zhong2015:n_searches'


def _type_name(code: str) -> str:
    if code not in CAM_TYPES:
        raise ValueError(f'Unknown Zhong2015 CAM type: {code!r}')
    return CAM_TYPES[code]


def _flag(value: bool) -> str:
    """A 1/0 column of the table, read by the parser into a boolean, as 'true'/'false'."""
    if value is True:
        return 'true'
    if value is False:
        return 'false'
    raise ValueError(f'Not a flag: {value!r}')


def _n_searches(row: dict[str, Any]):
    value = measurement(row['sources'], source_field='Sum')
    if value is None:
        raise ValueError(f'Zhong2015 {row["symbol"]}: Sum is {row["sources"]!r}')
    return value


def map_gene(row: dict[str, Any]) -> Relation:
    """A gene with its CAM type."""
    _type_name(row['cam_type'])
    identifiers = [Identifier(type=Namespace.GENESYMBOL, value=row['symbol'])]
    if row['entrez']:
        identifiers.append(Identifier(type=Namespace.ENTREZ, value=row['entrez']))
    if row['name']:
        identifiers.append(Identifier(type=Namespace.NAME, value=row['name']))
    return Relation(
        subject=Entity(
            type=Gene,
            identifiers=identifiers,
            annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
        ),
        # The type of a cell adhesion molecule is a functional class of the gene.
        predicate=slots.associated_with,
        object=Entity(
            type=OntologyClass,
            identifiers=[Identifier(type=TYPE_NAMESPACE, value=row['cam_type'])],
        ),
        annotations=[
            *(
                Annotation(term=term, value=_flag(row[field]))
                for field, term in _FOUND_BY.items()
            ),
            Annotation(term=_SOURCES, value=_n_searches(row)),
        ],
    )


def map_cam_type(row: dict[str, Any]) -> Entity:
    """A CAM type."""
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=TYPE_NAMESPACE, value=row['cam_type']),
            Identifier(type=Namespace.NAME, value=_type_name(row['cam_type'])),
        ],
    )


resource = Resource(
    config,
    terms=Dataset(
        download=download,
        mapper=map_cam_type,
        raw_parser=iter_cam_types,
    ),
    genes=Dataset(
        download=download,
        mapper=map_gene,
        raw_parser=iter_genes,
    ),
)
