"""
Parse RefMet data and emit Entity records.

This module converts metabolites reference nomenclature into Entity
records using the declarative schema pattern.
"""

from pypath.inputs_v2.parsers.base import iter_csv
from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    MoleculeAnnotationsCv,
)

# =================================== SET-UP ===================================

BASE_URL = 'https://www.metabolomicsworkbench.org/databases/refmet/refmet_download.php'

config = ResourceConfig(
    id=ResourceCv.REFMET,
    name='RefMet',
    url='https://www.metabolomicsworkbench.org/databases/refmet/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='33199890',
    primary_category='metabolism',
    description=(
        'The main objective of RefMet is to provide a standardized reference '
        'nomenclature for both discrete metabolite structures and metabolite '
        'species identified by spectroscopic techniques in metabolomics '
        'experiments'
    ),
)

# ================================== DOWNLOAD ==================================

download = Download(
        url=BASE_URL,
        filename='refmet.csv',
        subfolder='refmet',
        large=True,
        ext='.csv',
        default_mode='r',
)

# =================================== SCHEMA ===================================

def is_lipid(x):

    return 'lipid' if 'lipid' in x.lower() else 'small_molecule'

CLASS_ENTITY = {
    'lipid': EntityTypeCv.LIPID,
    'small_molecule': EntityTypeCv.SMALL_MOLECULE,
}
SUPERCLASS = {
    'lipid': MoleculeAnnotationsCv.LIPID_CATEGORY,
    'small_molecule': MoleculeAnnotationsCv.COMPOUND_SUPERCLASS,
}
CLASS={
    'lipid': MoleculeAnnotationsCv.LIPID_MAIN_CLASS,
    'small_molecule': MoleculeAnnotationsCv.COMPOUND_CLASS,
}
SUBCLASS = {
    'lipid': MoleculeAnnotationsCv.LIPID_SUB_CLASS,
    'small_molecule': MoleculeAnnotationsCv.COMPOUND_SUBCLASS,
}

f = FieldConfig(
    extract={},
    map={
        'class_to_entity': CLASS_ENTITY,
        'super_class': SUPERCLASS,
        'class': CLASS,
        'sub_class': SUBCLASS,
    },
    transform={
        'is_lipid': lambda x: is_lipid(x),
    },
)

schema = EntityBuilder(
    entity_type=f('super_class', map='class_to_entity', transform='is_lipid'),
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.REFMET, value='refmet_id'),
        CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA, value='formula'),
        CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value='pubchem_id'),
        CV(term=IdentifierNamespaceCv.CHEBI, value='chebi_id'),
        CV(term=IdentifierNamespaceCv.HMDB, value='hmdb_id'),
        CV(term=IdentifierNamespaceCv.LIPIDMAPS, value='lipidmaps_id'),
        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value='kegg_id'),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value='inchi_key'),
        CV(term=IdentifierNamespaceCv.NAME, value='refmet_name'),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.MASS_DALTON, value='exactmass'),
        CV(
            term=f('super_class', transform='is_lipid', map='super_class'),
            value='super_class'
        ),
        CV(
            term=f('super_class', transform='is_lipid', map='class'),
            value='main_class'
        ),
        CV(
            term=f('super_class', transform='is_lipid', map='sub_class'),
            value='sub_class'
        ),
    ),
)

# ================================= RESOURCE ===================================

resource = Resource(
    config=config,
    metabolites=Dataset(
        download=download,
        mapper=schema,
        raw_parser=iter_csv,
    )
)

# ================================= REFERENCE ==================================
# refmet_id     RM0202329
# refmet_name   7,11-Dimethyl-3-methylene-1,6E,10-dodecatriene
# super_class   Fatty Acyls
# main_class    Hydrocarbons
# sub_class     Hydrocarbons
# formula       C15H24
# exactmass     204.1878
# pubchem_cid   5281517
# chebi_id      10418
# hmdb_id       HMDB0062763
# lipidmaps_id  LMFA11000040
# kegg_id       C09666
# inchi_key     JSNRRGGBADWTMC-NTCAYCPXSA-N
