"""
Parse MEBOCOST DB data and emit Entity records.

MEBOCOST DB is a curated resource of metabolite-sensor interactions collected
through computational text-mining and manual curation from PubMed abstracts
and databases like HMDB, Recon2, and GPCRdb.
"""

from __future__ import annotations

import re
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    ChemicalEntity,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.interaction_profiles import TRANSPORT_QUALIFIERS
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv

config = ResourceConfig(
    id=ResourceCv.MEBOCOST,
    name='MEBOCOST DB',
    url='https://github.com/kaifuchenlab/MEBOCOST',
    license=LicenseCV.BSD_3,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='40568942',
    primary_category='interactions',
    description='MEBOCOST DB is a curated resource of metabolite-sensor interactions collected through computational text-mining and manual curation from PubMed abstracts and databases like HMDB, Recon2, and GPCRdb.',
)
EVIDENCE_SOURCES = {'hmdb': 'HMDB', 'recon2': 'Recon2', 'cellinker': 'Cellinker',
                    'cellphonedb': 'CellPhoneDB', 'cellchat': 'CellChat'}
f = FieldConfig(delimiter='; ')


def _evidence(row):
    """The Evidence column's PubMed IDs, source databases, URLs and remaining notes."""
    parsed = {'pubmed': [], 'source': [], 'url': [], 'comment': []}
    for item in re.split(r'[;,]\s*', str(row.get('Evidence') or '')):
        item = item.strip()
        if not item:
            continue
        if item.isdigit():
            parsed['pubmed'].append(f'PMID:{item}')
        elif item.lower() in EVIDENCE_SOURCES:
            parsed['source'].append(EVIDENCE_SOURCES[item.lower()])
        elif item.startswith(('http://', 'https://')):
            parsed['url'].append(item)
        else:
            parsed['comment'].append(item)
    return parsed


def _is_transporter(row):
    return 'transporter' in {v.strip().lower() for v in str(row.get('Annotation') or '').split(';')}


def get_interactions_schema(taxon_id: str) -> RelationBuilder:
    """
    Generate the interaction schema for a specific taxon.

    Args:
        taxon_id: NCBI taxonomy ID.

    Returns:
        RelationBuilder for MEBOCOST interactions.
    """
    metabolite_builder = EntityBuilder(
        entity_type=ChemicalEntity,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.HMDB, value=f('HMDB_ID', extract='(HMDB\\d+)')),
            CV(term=Namespace.NAME, value=f('standard_metName')),
            CV(term=Namespace.SYNONYM, value=f('metName', delimiter='; ')),
        ),
        annotations=AnnotationsBuilder(),
    )
    sensor_builder = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('Gene_name')),
            CV(term=Namespace.NAME, value=f('Protein_name')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}')
        ),
    )
    return RelationBuilder(
        subject=lambda row: (sensor_builder if _is_transporter(row) else metabolite_builder).build(row),
        predicate=lambda row: slots.affects if _is_transporter(row) else slots.interacts_with,
        object=lambda row: (metabolite_builder if _is_transporter(row) else sensor_builder).build(row),
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.MEBOCOST, value=f('ID'))
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.original_predicate, value=f('Annotation')),
            *(CV(term=term, value=lambda row, value=value: value if _is_transporter(row) else None) for term, value in TRANSPORT_QUALIFIERS),
            CV(term=slots.publications, value=lambda row: _evidence(row)['pubmed']),
            CV(term=slots.supporting_data_source, value=lambda row: _evidence(row)['source']),
            CV(term=slots.source_record_urls, value=lambda row: _evidence(row)['url']),
            CV(term=slots.description, value=lambda row: _evidence(row)['comment']),
        ),
    )


resource = Resource(
    config,
    human=Dataset(
        download=Download(
            url='https://raw.githubusercontent.com/kaifuchenlab/MEBOCOST/main/data/mebocost_db/human/human_met_sensor_update_Oct21_2025.tsv',
            filename='mebocost_human.tsv',
            subfolder='mebocost',
        ),
        mapper=get_interactions_schema('9606'),
        raw_parser=iter_tsv,
    ),
    mouse=Dataset(
        download=Download(
            url='https://raw.githubusercontent.com/kaifuchenlab/MEBOCOST/main/data/mebocost_db/mouse/mouse_met_sensor_update_Oct21_2025.tsv',
            filename='mebocost_mouse.tsv',
            subfolder='mebocost',
        ),
        mapper=get_interactions_schema('10090'),
        raw_parser=iter_tsv,
    ),
)
