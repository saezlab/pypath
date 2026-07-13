"""
Parse Hormone2Cell data and emit Entity records.

This module converts hormone-cell interaction information from Hormone2Cell into
Entity records using the declarative schema pattern.
"""

from __future__ import annotations

import pandas as pd

from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers import brenda as _parsers
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    OntologyCv,
    UpdateCategoryCV,
    ResourceCv,
)

# =================================== SET-UP ===================================

URL = 'https://hormonecellatlas.org.uk/accounts/downloads/file/table-s2b/'

config = ResourceConfig(
    id=ResourceCv.HORMONE2CELL,
    name='Hormone2Cell',
    url='https://hormonecellatlas.org.uk/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='42207862',
    primary_category='intercell',
    description=(
        'The Hormone Cell Atlas uses single cell transcriptomics to map '
        'hormone production and action across 47 disease-free human tissues at '
        'cellular resolution.'
    ),
)


download = Download(
    url=URL,
    filename='hormone2cell.xlsx',
    subfolder='hormone2cell',
    large=True,
    ext='xlsx',
    default_mode='rb',
)

def parser(opener):

    pd.read_excel(opener.path)

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={},
    map={},
    transform={},
)

schema = EntityBuilder()

# ================================= RESOURCE ===================================

#resource = Resource(
#    config=config,
#    data=Dataset(
#        download=download,
#        mapper=schema,
#        raw_parser=parser,
#    ),
#)

# ================================= REFERENCE ==================================
