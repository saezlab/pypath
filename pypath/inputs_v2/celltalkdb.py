"""
Parse CellTalkDB ligand-receptor pairs and emit Relation records.

CellTalkDB distributes one R data frame per species. Each row is a
ligand-receptor pair named by gene symbols, with NCBI Gene and Ensembl
cross-references for both partners and the PubMed IDs supporting the pair.
"""

from __future__ import annotations
import functools
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
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
from pypath.inputs_v2.parsers.celltalkdb import iter_lr_pairs

config = ResourceConfig(
    id=ResourceCv.CELLTALKDB,
    name='CellTalkDB',
    url='https://github.com/ZJUFanLab/CellTalkDB',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.GPL_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='33147626',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Literature-supported ligand–receptor pairs for human and mouse.',
    description='CellTalkDB is a manually curated database of ligand-receptor interactions in human and mouse. Candidate pairs were taken from protein-protein interactions in STRING between annotated ligands and receptors and kept after manual verification in the literature, so every pair carries the PubMed IDs that support it.',
)
BASE_URL = (
    'https://raw.githubusercontent.com/ZJUFanLab/CellTalkDB/master/database/'
)
SPECIES = {'human': '9606', 'mouse': '10090'}


def _download(species: str) -> Download:
    return Download(
        url=f'{BASE_URL}{species}_lr_pair.rds',
        filename=f'celltalkdb_{species}_lr_pair.rds',
        subfolder='celltalkdb',
        default_mode='rb',
        encoding=None,
    )


def _parser(species: str) -> Any:
    return functools.partial(iter_lr_pairs, taxon_id=SPECIES[species])


f = FieldConfig(
    extract={
        # The evidence column lists PubMed IDs. Stricter than the shared
        # PubMed identifier rule (1 to 8 digits), because the files have cut
        # numbers: every intact ID has seven or eight digits, and the four
        # others are fragments or typos that are left out, not repaired. At
        # NCBI, 91 is a 1975 paper on a membrane ATPase (it follows 9486650 in
        # the cell of HLA-E and KLRC2), 93 a 1975 paper on lysozyme adsorption
        # (it follows 9405670, HLA-G and KLRD1), 943050 a 1976 paper on
        # calcification of the basal ganglia (IL12A and IL12RB2); 321961154
        # has nine digits.
        'pmid': r'^([1-9]\d{6,7})$',
    },
    transform={
        'pmid_curie': lambda v: f'PMID:{v}',
        'taxon_curie': lambda v: f'NCBITaxon:{v}',
    },
)


def _participant(role: str) -> EntityBuilder:
    # The source names a pair by its gene symbols, so the symbol leads; the
    # Gene and Ensembl identifiers are the source's cross-references for it.
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(f'{role}_gene_symbol')),
            CV(term=Namespace.ENTREZ, value=f(f'{role}_gene_id')),
            CV(term=Namespace.ENSG, value=f(f'{role}_ensembl_gene_id')),
            CV(term=Namespace.ENSP, value=f(f'{role}_ensembl_protein_id')),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.in_taxon,
                value=f('taxon_id', transform='taxon_curie'),
            ),
            # The columns are named for the role (`ligand_*`, `receptor_*`).
            CV(term=PARTICIPANT_ROLE, value=role),
        ),
    )


# The ligand is the subject and the receptor the object.
interactions_schema = RelationBuilder(
    subject=_participant('ligand'),
    predicate=slots.interacts_with,
    object=_participant('receptor'),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('lr_pair')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.publications,
            value=f(
                'evidence',
                delimiter=',',
                extract='pmid',
                transform='pmid_curie',
            ),
        ),
    ),
)
resource = Resource(
    config,
    human_interactions=Dataset(
        download=_download('human'),
        mapper=interactions_schema,
        raw_parser=_parser('human'),
    ),
    mouse_interactions=Dataset(
        download=_download('mouse'),
        mapper=interactions_schema,
        raw_parser=_parser('mouse'),
    ),
)
