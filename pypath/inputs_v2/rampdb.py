"""
Parse RaMP-DB data and emit Entity records.

This module converts annotations of lipids and metabolites into Entity records
using the declarative schema pattern.
"""

from pypath.inputs_v2.base import ResourceConfig
from pypath.internals.cv_terms import ResourceCv, LicenseCV, UpdateCategoryCV

config = ResourceConfig(
    id=ResourceCv.RAMPDB,
    name='RaMP-DB',
    url='https://rampdb.nh.gov/',
    license=LicenseCV.GPL_2_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='36373969',
    description=(
        'RaMP-DB (Relational database of Metabolomic Pathways) is a '
        'multi-sourced integrated database with comprehensive annotations on '
        'biological pathways, structure/chemistry, disease and ontology '
        'annotations for genes, proteins, and metabolites.'
    )
)