"""
Parse BRENDA data and emit Entity records.

This module converts enzyme and ligand information from BRENDA into Entity
records using the declarative schema pattern.
"""

from __future__ import annotations

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

URL = 'https://www.brenda-enzymes.org/download.php'
_ENZYME_CLASSIFICATION_ONTOLOGY_ID = 'enzyme_classification'

config = ResourceConfig(
    id=ResourceCv.BRENDA,
    name='BRENDA',
    url='https://www.brenda-enzymes.org/index.php',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='41206471',
    primary_category='enzsub',
    annotation_ontologies=(OntologyCv.ENZYME_CLASSIFICATION,),
    description=(
        'BRENDA is the main collection of enzyme functional data available to '
        'the scientific community.'
    ),
)


download = Download(
    url=URL,
    filename='brenda.txt.tar.gz',
    subfolder='brenda',
    large=True,
    ext='tar.gz',
    default_mode='rb',
    download_kwargs={
        'post': True,
        'query': {'dlfile': 'dl-textfile', 'accept-license': '1'},
    }
)

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={},
    map={},
    transform={},
)

enzyme_ontology_schema = ontology_entity_mapper(
    _parsers.term_record_to_term,
    ontology_id=_ENZYME_CLASSIFICATION_ONTOLOGY_ID,
)

schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('UniProt')),
    ),
    associations=AssociationsBuilder(
        AssociationBuilder(
            object_entity_type=EntityTypeCv.CV_TERM,
            object_identifier_type=IdentifierNamespaceCv.CV_TERM_ACCESSION,
            object_identifier=f(_parsers.ec_term_id),
        ),
    ),
)

# ================================= RESOURCE ===================================

resource = Resource(
    config=config,
    enzyme_ontology=Dataset(
        download=download,
        mapper=enzyme_ontology_schema,
        raw_parser=_parsers.iter_enzyme_ontology,
    ),
    protein_enzyme_class_annotations=Dataset(
        download=download,
        mapper=schema,
        raw_parser=_parsers.iter_protein_enzyme_class_annotations,
    ),
    data=Dataset(
        download=download,
        mapper=schema,
        raw_parser=_parsers.parser,
    ),
)

# ================================= REFERENCE ==================================
#{ Example entry
# X   'EC': '1.1.1.1', # Enzyme classification
# X   'UniProt': set(), # UniProt IDs as set
#    '#': '6', # Internal BRENDA identifier, not used
# X   'Organism': {'Mus musculus'}, # Species name
# X   'Refs': [ # References as a list of PMIDS
#        '14756569',
#        ...
#    ],
# X    'Activator': { # Activating compounds as set
#        'S-nitrosoglutathione',
#        'Valeramide',
#        'butyramide',
#        'capronamide',
#        'tert-butanol'
#    },
# X   'Inhibitor': { # Inhibiting compounds as set
#        '1,10-phenanthroline',
#        '4-Methylpyrazole',
#        'Vanillin',
#        'caffeic acid',
#        'dodecanoic acid',
#        'ellagic acid',
#        'pyrazole',
#        'syringaldehyde'
#    },
# X    'Cofactor': {'NAD+', 'NADH'}, # Cofactors as set
# X   'Ki': { # Ki constant as "concentration {compound}"
#        '0.00008 {caffeic acid}',
#        '0.0051 {pyrazole}',
#        '0.0079 {Vanillin}',
#        '0.0156 {syringaldehyde}',
#        '0.022 {ellagic acid}'
#    },
# X   'Km': { # KM constant as "concentration {compound}"
#        '-999 {more}',
#        '0.006 {Hexanol}',
#        '0.085 {Hexanol}',
#        '0.48 {ethanol}',
#        '0.63 {Hexanol}',
#        '0.83 {ethanol}',
#        '1.9 {Hexanol}',
#        '1625 {ethanol}',
#        '255 {ethanol}'
#    }
#}
