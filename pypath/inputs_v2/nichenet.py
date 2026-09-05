"""
Parse NicheNet v2 ligand-receptor networks and emit Entity records.
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
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Relation,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.nichenet import iter_lr_network

config = ResourceConfig(
    id=ResourceCv.NICHENET,
    name='NicheNet v2 ligand-receptor network',
    url='https://zenodo.org/records/7074291',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='31819264',
    primary_category='interactions',
    description='NicheNet v2 ligand-receptor networks provide prior knowledge for modeling intercellular communication by linking sender-cell ligands to receiver-cell receptors and downstream target gene regulation.',
)
BASE_URL = 'https://zenodo.org/records/7074291/files/'
SPECIES = {'human': '9606', 'mouse': '10090'}


def _download(species: str) -> Download:
    return Download(
        url=f'{BASE_URL}lr_network_{species}_21122021.rds?download=1',
        filename=f'nichenet_lr_network_{species}_21122021.rds',
        subfolder='nichenet',
        default_mode='rb',
        encoding=None,
    )


def _parser(species: str) -> Any:
    return functools.partial(iter_lr_network, taxon_id=SPECIES[species])


def _protein(gene: str, taxon_id: str, role: str) -> Entity:
    return Entity(
        type=Protein,
        identifiers=[Identifier(type=Namespace.GENESYMBOL, value=gene)],
        annotations=[
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}')
        ],
    )


def map_nichenet_interaction(row: dict[str, Any]) -> Relation | None:
    ligand = str(row.get('ligand') or '').strip()
    receptor = str(row.get('receptor') or '').strip()
    taxon_id = str(row.get('taxon_id') or '').strip()
    if not ligand or not receptor:
        return None
    return Relation(
        subject=_protein(ligand, taxon_id, 'ligand'),
        predicate=slots.interacts_with,
        object=_protein(receptor, taxon_id, 'receptor'),
        annotations=[],
        identifiers=[
            Identifier(type=Namespace.NAME, value=f'{ligand} - {receptor}')
        ],
    )


resource = Resource(
    config,
    human_interactions=Dataset(
        download=_download('human'),
        mapper=map_nichenet_interaction,
        raw_parser=_parser('human'),
    ),
    mouse_interactions=Dataset(
        download=_download('mouse'),
        mapper=map_nichenet_interaction,
        raw_parser=_parser('mouse'),
    ),
)
