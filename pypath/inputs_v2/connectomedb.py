"""
Parse ConnectomeDB2025 data and emit Entity records.

This module converts ConnectomeDB2025 interactions and complexes into Entity
records using the declarative schema pattern.
"""

from __future__ import annotations
import re
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
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
from pypath.inputs_v2.parsers.base import iter_csv

config = ResourceConfig(
    id=ResourceCv.CONNECTOMEDB,
    name='ConnectomeDB2025',
    url='https://connectomedb.org/',
    license=LicenseCV.MIT,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='41171146',
    primary_category='interactions',
    description='ConnectomeDB is a comprehensive and ongoing project that provides ahigh-quality manually curated database of interacting ligand-receptorpairs for use in cell-cell communication analysis. First released in2015 (Ramilowski, et al.), and subsequently updated in 2020 (Hou, et al.)and 2025 (Liu, Maezono et al.), it aims to enhance the understanding ofcell-cell communication in humans and other mammals, supporting biologicaland medical research.',
)
BASE_URL = 'https://connectomedb.org/downloads/Current-Release/CSV/'
download_interactions = Download(
    url=BASE_URL + 'all_species.csv',
    filename='connectomedb_all_species_interactions.csv',
    subfolder='connectomedb2025',
)
_symbols_pat = re.compile('^([^,(]+)(?:\\s*\\((.+)\\))?')


def _extract_primary_gene(token: str):
    match = _symbols_pat.search(token or '')
    return match.group(1).strip() if match else None


def _extract_gene_alias(token: str):
    match = _symbols_pat.search(token or '')
    return (
        match.group(2).replace(', ', ';') if match and match.group(2) else None
    )


_species_taxon = {
    'human': '9606',
    'mouse': '10090',
    'chimp': '9598',
    'macaque': '9544',
    'marmoset': '9483',
    'rat': '10116',
    'pig': '9823',
    'cow': '9913',
    'dog': '9615',
    'horse': '9796',
    'sheep': '9940',
    'chicken': '9031',
    'frog': '8364',
    'zebrafish': '7955',
}


def _species_to_taxon(species: str) -> str | None:
    taxon = _species_taxon.get((species or '').strip().lower())
    return f'NCBITaxon:{taxon}' if taxon else None


_cdb_pat = re.compile('^CDB\\d{2}:(\\d+)', re.IGNORECASE)
_hgnc_pat = re.compile('^HGNC:(\\d+)$', re.IGNORECASE)


def _extract_cdb(val: str) -> str | None:
    return val if _cdb_pat.match(val) else None


def _extract_hgnc_id(val: str) -> str | None:
    match = _hgnc_pat.match((val or '').strip())
    return match.group(1) if match else None


f = FieldConfig(
    extract={
        'primary_gene': _extract_primary_gene,
        'gene_alias': _extract_gene_alias,
        'cdb': _extract_cdb,
        'hgnc_id': _extract_hgnc_id,
    },
    map={'species_taxon': _species_to_taxon},
)
ligand_builder = EntityBuilder(
    entity_type=Protein,
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.GENESYMBOL,
            value=f('Ligand Symbols', extract='primary_gene'),
        ),
        CV(
            term=Namespace.HGNC, value=f('Ligand Species ID', extract='hgnc_id')
        ),
        CV(term=Namespace.ENSG, value=f('Ligand ENSEMBL ID')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=f('Species', map='species_taxon'))
    ),
)
receptor_builder = EntityBuilder(
    entity_type=Protein,
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.GENESYMBOL,
            value=f('Receptor Symbols', extract='primary_gene'),
        ),
        CV(
            term=Namespace.HGNC,
            value=f('Receptor Species ID', extract='hgnc_id'),
        ),
        CV(term=Namespace.ENSG, value=f('Receptor ENSEMBL ID')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=f('Species', map='species_taxon'))
    ),
)
interactions_schema = RelationBuilder(
    subject=ligand_builder,
    predicate=slots.interacts_with,
    object=receptor_builder,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.CDB, value=f('Interaction ID', extract='cdb')),
        CV(term=Namespace.NAME, value=f('LR Pair')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=f('Species', map='species_taxon'))
    ),
)
resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=iter_csv,
    ),
)
