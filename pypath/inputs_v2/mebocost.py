"""
Parse MEBOCOST DB data and emit Entity records.

MEBOCOST DB is a curated resource of metabolite-sensor interactions collected
through computational text-mining and manual curation from PubMed abstracts
and databases like HMDB, Recon2, and GPCRdb.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    CurationCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv


config = ResourceConfig(
    id=ResourceCv.MEBOCOST,
    name='MEBOCOST DB',
    url='https://github.com/kaifuchenlab/MEBOCOST',
    license=LicenseCV.GPL_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='36573906',
    description=(
        'MEBOCOST DB is a curated resource of metabolite-sensor interactions '
        'collected through computational text-mining and manual curation '
        'from PubMed abstracts and databases like HMDB, Recon2, and GPCRdb.'
    ),
)

# Evidence source mapping
evidence_source_map = {
    'HMDB': IdentifierNamespaceCv.HMDB,
    'Recon2': IdentifierNamespaceCv.RECON2,
    'Celllinker': IdentifierNamespaceCv.CELLINKER,
    'CellPhoneDB': IdentifierNamespaceCv.CELLPHONEDB,
    'CellChat': IdentifierNamespaceCv.CELLCHAT,
}

f = FieldConfig(
    extract={
        'pubmed': r'^(\d+)$',
        'source': lambda v: evidence_source_map.get(v),
        'comment': lambda v: v if v.startswith('http') 
        or (not v.isdigit() and v not in evidence_source_map) else None,
    },
    delimiter='; ',
)


def get_interactions_schema(taxon_id: str) -> EntityBuilder:
    """
    Generate the interaction schema for a specific taxon.

    Args:
        taxon_id: NCBI taxonomy ID.

    Returns:
        EntityBuilder for MEBOCOST interactions.
    """
    return EntityBuilder(
        entity_type=EntityTypeCv.INTERACTION,
        # Generic MEBOCOST ID which is a number, but prefix with "MEBOCOST:" to avoid namespace collisions
        # Not a Stable Identifier since it's not guaranteed to be stable across releases, 
        # but still useful for tracing back to the source record
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.MEBOCOST, value=f('ID')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=IdentifierNamespaceCv.PUBMED, value=f('Evidence', extract='pubmed')),
            CV(term=f('Evidence', extract='source'), value=f('Evidence')),
            CV(term=CurationCv.COMMENT, value=f('Evidence', extract='comment')),
        ),
        membership=MembershipBuilder(
            # Member 1: Metabolite (Small Molecule)
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.SMALL_MOLECULE,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.HMDB, value=f('HMDB_ID')),
                        CV(term=IdentifierNamespaceCv.NAME, value=f('standard_metName')),
                        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('metName', delimiter='; ')),
                    ),
                ),
            ),
            # Member 2: Sensor (Protein)
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('Gene_name')),
                        CV(term=IdentifierNamespaceCv.NAME, value=f('Protein_name')),
                    ),
                    annotations=AnnotationsBuilder(
                        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
                    ),
                ),
            ),
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
