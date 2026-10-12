"""
Parse the ligand-receptor pairs of Ramilowski et al. 2015 and emit Relation records.

Supplementary Data 2 of the paper is a workbook. Its `All.Pairs` sheet lists
every candidate pair with the databases it was found in and the authors'
verdict; the `HPMR`, `IUPHAR` and `HPRD` sheets carry the PubMed IDs those
databases give for a pair. Only pairs the authors mark as literature supported
are emitted, with the databases the table flags for them.
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
from pypath.inputs_v2.parsers.ramilowski2015 import (
    PMID_SHEETS,
    iter_pairs,
    split_pubmed_cell,
)

config = ResourceConfig(
    id=ResourceCv.RAMILOWSKI2015,
    name='Ramilowski2015',
    full='Ramilowski 2015 ligand-receptor pairs',
    url='https://doi.org/10.1038/ncomms8866',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='26198319',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Literature-supported human ligand–receptor pairs with their source databases and publications.',
    description='Ramilowski et al. 2015 drew a draft map of cell-to-cell communication between human primary cell types for the FANTOM5 project. Known ligand-receptor pairs from DLRP, IUPHAR and HPMR were merged with interactions between secreted and plasma membrane proteins found in HPRD and STRING, and curated by hand into 1,894 pairs with primary literature support, the reference set, and 528 putative pairs without it.',
)
download = Download(
    url='https://static-content.springer.com/esm/art%3A10.1038%2Fncomms8866/MediaObjects/41467_2015_BFncomms8866_MOESM611_ESM.xlsx',
    filename='ramilowski2015_supplementary_data_2.xlsx',
    subfolder='ramilowski2015',
    default_mode='rb',
    encoding=None,
)
TAXON = 'NCBITaxon:9606'
PMID_COLUMNS = ('PMID.Manual', *(f'{sheet}.PMID' for sheet in PMID_SHEETS))
# `Pair.Evidence`, the authors' verdict on a candidate pair: is it emitted?
# Putative pairs have no primary literature in any sheet, and the pairs marked
# EXCLUDED were rejected because a partner is not a ligand or not a receptor.
EVIDENCE = {
    'literature supported': True,
    'putative': False,
    'EXCLUDED': False,
    'EXCLUDED not ligand': False,
    'EXCLUDED not receptor': False,
}
# One column per database the pair was found in. A cell repeats the column's
# name when the pair was found there and is empty otherwise (the parser raises
# on any other value); both STRING columns name STRING (the same identifiers
# as `lrdb.py` and `talklr.py`).
DATABASE_COLUMNS = {
    'DLRP': 'DLRP',
    'HPMR': 'HPMR',
    'IUPHAR': 'IUPHAR',
    'HPRD': 'HPRD',
    'STRING.binding': 'STRING',
    'STRING.experiment': 'STRING',
}


def _databases(row: dict[str, Any]) -> list[str]:
    """Databases that list the pair, each once, in the order of the columns."""
    out: list[str] = []
    for column, database in DATABASE_COLUMNS.items():
        if row.get(column) is not None and database not in out:
            out.append(database)
    return out


def _pubmeds(row: dict[str, Any]) -> list[str]:
    """PubMed IDs of a pair in source order, each once."""
    out: list[str] = []
    for column in PMID_COLUMNS:
        for pmid in split_pubmed_cell(row.get(column), column):
            if pmid not in out:
                out.append(pmid)
    return out


def _protein(symbol: str, name: Any, role: str) -> Entity:
    identifiers = [Identifier(type=Namespace.GENESYMBOL, value=symbol)]
    if name:
        identifiers.append(Identifier(type=Namespace.NAME, value=str(name)))
    return Entity(
        type=Protein,
        identifiers=identifiers,
        annotations=[
            Annotation(term=slots.in_taxon, value=TAXON),
            # The columns are named for the role (`Ligand.*`, `Receptor.*`).
            Annotation(term=PARTICIPANT_ROLE, value=role),
        ],
    )


def map_ramilowski_interaction(row: dict[str, Any]) -> Relation | None:
    evidence = row.get('Pair.Evidence')
    if evidence not in EVIDENCE:
        raise ValueError(f'Unknown Ramilowski 2015 pair evidence {evidence!r}')
    ligand = row['Ligand.ApprovedSymbol']
    receptor = row['Receptor.ApprovedSymbol']
    if not EVIDENCE[evidence] or not ligand or not receptor:
        return None
    return Relation(
        subject=_protein(ligand, row.get('Ligand.Name'), 'ligand'),
        predicate=slots.interacts_with,
        object=_protein(receptor, row.get('Receptor.Name'), 'receptor'),
        identifiers=[
            Identifier(type=Namespace.NAME, value=name)
            for name in [row.get('Pair.Name')]
            if name
        ],
        annotations=[
            Annotation(term=slots.in_taxon, value=TAXON),
            *[
                Annotation(term=slots.supporting_data_source, value=database)
                for database in _databases(row)
            ],
            *[
                Annotation(term=slots.publications, value=f'PMID:{pmid}')
                for pmid in _pubmeds(row)
            ],
        ],
    )


resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=map_ramilowski_interaction,
        raw_parser=iter_pairs,
    ),
)
