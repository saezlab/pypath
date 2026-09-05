"""
Parse MEBOCOST DB data and emit Entity records.

MEBOCOST DB is a curated resource of metabolite-sensor interactions collected
through computational text-mining and manual curation from PubMed abstracts
and databases like HMDB, Recon2, and GPCRdb.
"""

from __future__ import annotations
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    ChemicalEntity,
    slots,
)
from omnipath_core.naming import Namespace
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
evidence_source_map = {
    name: name
    for name in ('HMDB', 'Recon2', 'Celllinker', 'CellPhoneDB', 'CellChat')
}
f = FieldConfig(
    extract={
        'pubmed': '^(\\d+)$',
        'source': lambda v: evidence_source_map.get(v),
        'comment': lambda v: v
        if v.startswith('http')
        or (not v.isdigit() and v not in evidence_source_map)
        else None,
    },
    delimiter='; ',
)


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
        subject=metabolite_builder,
        predicate=slots.interacts_with,
        object=sensor_builder,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.MEBOCOST, value=f('ID'))
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.publications,
                value=f(
                    'Evidence',
                    extract='pubmed',
                    transform=lambda v: f'PMID:{v}',
                ),
            ),
            CV(
                term=slots.supporting_data_source,
                value=f('Evidence', extract='source'),
            ),
            CV(term=slots.description, value=f('Evidence', extract='comment')),
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
