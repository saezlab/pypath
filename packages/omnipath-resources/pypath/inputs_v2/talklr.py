"""
Parse the talklr ligand-receptor database and emit Relation records.

The talklr R package ships one data frame of human ligand-receptor pairs, the
pair list of Ramilowski et al. 2015. Each row names the pair, gives the
approved gene symbol and a name for the ligand and for the receptor, flags the
databases the pair was found in, and carries the authors' verdict: literature
supported or putative. The two verdicts are separate datasets.
"""

from __future__ import annotations
from collections.abc import Callable
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
from pypath.inputs_v2.parsers.talklr import DATABASE_COLUMNS, iter_pairs

config = ResourceConfig(
    id=ResourceCv.TALKLR,
    name='talklr',
    url='https://github.com/yuliangwang/talklr',
    # The DESCRIPTION file of the package states "GPL (>= 2)" (checked when
    # the module was written, October 2026); the legacy registry recorded GPLv3.
    license=LicenseCV.GPL_2_0,
    update_category=UpdateCategoryCV.STATIC,
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Human ligand–receptor pairs from the talklr package, literature supported and putative.',
    description='talklr is an R package that uses the Kullback-Leibler divergence to find ligand-receptor pairs with cell type specific or perturbed expression in single-cell RNA sequencing data. Its built-in database is the human ligand-receptor list of Ramilowski et al. 2015: pairs from DLRP, IUPHAR and HPMR and interactions between secreted and plasma membrane proteins from HPRD and STRING, each marked as supported by primary literature or as putative.',
)
download_interactions = Download(
    url='https://raw.githubusercontent.com/yuliangwang/talklr/master/data/receptor_ligand.rda',
    filename='talklr_receptor_ligand.rda',
    subfolder='talklr',
    default_mode='rb',
    encoding=None,
)
# The table has no species column. The symbols are HGNC approved symbols.
TAXON = 'NCBITaxon:9606'
# One column per database (`DATABASE_COLUMNS`). A cell repeats the column's name
# when the pair was found there and is empty otherwise; the reader raises on
# any other value. Both STRING columns name STRING.
DATABASES = {
    'DLRP': 'DLRP',
    'HPMR': 'HPMR',
    'IUPHAR': 'IUPHAR',
    'HPRD': 'HPRD',
    'STRING.binding': 'STRING',
    'STRING.experiment': 'STRING',
}

f = FieldConfig(
    extract={
        # Stricter than the shared PubMed identifier rule (1 to 8 digits): every
        # intact PMID in the table has seven or eight digits. The one other
        # token, 943050 for IL12A_IL12RB2, is a fragment and is left out, not
        # repaired; at NCBI it is a 1976 paper on calcification of the basal
        # ganglia.
        'pmid': r'^([1-9]\d{6,7})$',
    },
    map={'database': DATABASES},
    transform={
        'pmid_curie': lambda v: f'PMID:{v}',
    },
)


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
        # The approved symbol identifies the participant; a name alone does not.
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
    annotations=AnnotationsBuilder(
        *(
            CV(
                term=slots.supporting_data_source,
                value=f(column, map='database'),
            )
            for column in DATABASE_COLUMNS
        ),
        CV(
            term=slots.publications,
            value=f(
                'PMID.Manual',
                delimiter=',',
                extract='pmid',
                transform='pmid_curie',
            ),
        ),
    ),
)
resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=functools.partial(iter_pairs, evidence='literature supported'),
    ),
    putative_interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=functools.partial(iter_pairs, evidence='putative'),
    ),
)
