"""
Parse Human Serum Metabolome data and emit Entity records.

This module converts small molecule metabolite information from the Human Serum
Metabolite Database into Entity records using the declarative schema pattern.
"""

from __future__ import annotations
import re
from functools import partial

import pandas as pd
import lxml.etree as etree

from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
    ontology_entity_mapper,
)
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    MembershipBuilder,
    MembersFromList,
    Member,
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    InterCellAnnotations,
    MoleculeAnnotationsCv,
    AssayAnnotationsCv,
    BiologicalRoleCv,
    MoleculeSubtypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    ParticipantMetadataCv,
    LicenseCV,
    ProteinFunctionalClassCv,
    OntologyCv,
    UpdateCategoryCV,
    ResourceCv,
)

# =================================== SET-UP ===================================

URL = 'https://www.serummetabolome.ca/downloads/%s'

config = ResourceConfig(
    id=ResourceCv.HSM,
    name='Human Serum Metabolome',
    url='https://www.serummetabolome.ca/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='21359215',
    primary_category='small_molecules',
    description=(
        'The Serum Metabolome database is a freely available electronic '
        'database containing detailed information about 4651 small molecule '
        'metabolites found in human serum along with 10895 concentration '
        'values.'
    ),
)

download_metabolites = Download(
    url=URL,
    filename='',
    subfolder='',
    large=True,
    ext='',
    default_mode='r',
)


def parser(opener):

    return


# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={},
    map={},
    transform={},
)

#schema = EntityBuilder()

# ================================= RESOURCE ===================================

#resource = Resource(
#    config=config,
#    data=Dataset(
#        download=download,
#        mapper=schema,
#        raw_parser=parser
#    )
#)

# ================================= REFERENCE ==================================
