"""
Parse Human Serum Metabolome data and emit Entity records.

This module converts small molecule metabolite information from the Human Serum
Metabolite Database into Entity records using the declarative schema pattern.
"""

from __future__ import annotations
import re
from functools import partial
import xml.etree.ElementTree as ET

import pandas as pd

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

URL = 'https://www.serummetabolome.ca/system/downloads/current/%s'

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

download = Download(
    url=URL % 'serum_metabolites.zip',
    filename='serum_metabolites.zip',
    subfolder='serum_metabolome',
    large=True,
    ext='.zip',
    default_mode='rb',
    needed=['serum_metabolites.xml']
)


def _get(element, tag, xmlns='{http://www.hmdb.ca}'):

    child = element.find(f'{xmlns}{tag}')

    return '' if child is None else child.text.strip()


def parser(opener):

    tree = ET.parse(opener.result['serum_metabolites.xml'])
    root = tree.getroot()

    yield from [
        {i: _get(elem, i) for i in [
            '' # Fields - see reference below
        ]}
        for elem in root
    ]


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

# All fields in entries of serum_metabolites.xml

# abnormal_concentrations
# accession
# average_molecular_weight
# bigg_id
# biocyc_id
# biological_properties
# cas_registry_number
# chebi_id
# chemical_formula
# chemspider_id
# creation_date
# description
# diseases
# drugbank_id
# experimental_properties
# fbonto_id
# foodb_id
# general_references
# inchi
# inchikey
# iupac_name
# kegg_id
# knapsack_id
# metlin_id
# monisotopic_molecular_weight
# name
# normal_concentrations
# ontology
# pdb_id
# phenol_explorer_compound_id
# predicted_properties
# protein_associations
# pubchem_compound_id
# secondary_accessions
# smiles
# spectra
# state
# status
# synonyms
# synthesis_reference
# taxonomy
# traditional_iupac
# update_date
# version
# vmh_id
# wikipedia_id