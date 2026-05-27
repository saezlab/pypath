"""
Parse BRENDA data and emit Entity records.

This module converts enzyme and ligand information from BRENDA into Entity
records using the declarative schema pattern.
"""

from functools import partial

from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.parsers.macdb import iter_associations
from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
    OntologyDataset
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
from pypath.internals.ontology_builder import (
    OntologyBuilder,
    RelationshipBuilder
)
from pypath.internals.ontology_schema import OntologyTypedef
from pypath.internals.cv_terms import (
    EntityTypeCv,
    InteractionTypeCv,
    OntologyCv,
    DiseaseAnnotationCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    MoleculeSubtypeCv,
    MoleculeAnnotationsCv,
    AssayAnnotationsCv,
)

# =================================== SET-UP ===================================

URL = 'https://www.brenda-enzymes.org/download.php'
# XXX: As specified in inputs v1
RECORDS_ENABLED = {
    'ID', # EC-class
    'PR', # protein 	sequence identifier and name of source in {...}
    'AC', # activating compound
    'IN', # inhibitors
    'CF', # cofactor
    'KI', # Ki-value	inhibitor in {...}
    'KM', # KM-value	substrate in {...}
    'RF', # references
}

config = ResourceConfig(
    id=ResourceCv.BRENDA,
    name='BRENDA',
    url='https://www.brenda-enzymes.org/index.php',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='41206471',
    primary_category='enzsub',
    description=(
        'BRENDA is the main collection of enzyme functional data available to '
        'the scientific community.'
    ),
)

# ================================== DOWNLOAD ==================================

def parser(opener, **kwargs):

    keys = sorted(opener.result)
    # Serialized raw string
    result = opener.result[keys[0]].read().decode('utf-8')

    # NOTE: Continuation lines start with \t character -> merged from serialized
    # string by replacing \n\t with empty string before splitting the lines
    entries = [
        {   # Record key: value(s)
            split[0]: split[1:] if len(split) > 2 else split[1]
            for line in entry.split('\n')
            if (split := line.split('\t'))[0] in RECORDS_ENABLED
        }
        # Each entry spans multiple lines, separated by "///"
        for entry in result.replace('\n\t', '').split('///')[1:] # Ignore first
    ]
    # Controlling for empty, deleted, transferred or dummy records
    entries = [e for e in entries if (e and '(' not in e.get('ID', ''))]



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

association_schema = EntityBuilder(
)

# ================================= RESOURCE ===================================

#resource = Resource(
#    config=config,
#    associations=Dataset(
#        download=download,
#        mapper=schema,
#        raw_parser=None,
#    ),
#)

# ================================= REFERENCE ==================================
