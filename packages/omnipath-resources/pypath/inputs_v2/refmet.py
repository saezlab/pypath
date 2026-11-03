"""
Parse RefMet data and emit Entity records.

This module converts metabolites reference nomenclature into Entity
records using the declarative schema pattern.
"""

from biolink_model.datamodel.model import ChemicalEntity, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, Topic

from pypath.inputs_v2._class_terms import (
    Level,
    class_associations,
    iter_class_terms,
    map_class_term,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_csv
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    AssociationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)

BASE_URL = (
    'https://www.metabolomicsworkbench.org/databases/refmet/refmet_download.php'
)
config = ResourceConfig(
    id=ResourceCv.REFMET,
    name='RefMet',
    url='https://www.metabolomicsworkbench.org/databases/refmet/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='33199890',
    topics=(Topic.METABOLISM,),
    content=(Content.MOLECULAR_ANNOTATIONS, Content.IDENTIFIER_MAPPINGS),
    summary='Reference metabolite names, identifiers and chemical classifications.',
    description='The main objective of RefMet is to provide a standardized reference nomenclature for both discrete metabolite structures and metabolite species identified by spectroscopic techniques in metabolomics experiments',
)
download = Download(
    url=BASE_URL,
    filename='refmet.csv',
    subfolder='refmet',
    large=True,
    ext='.csv',
    default_mode='r',
)
# RefMet's classification, by name; a level per namespace, since a name can
# be a class at two levels ("Fatty Acyls").
CLASS_LEVELS = [
    Level('super_class', 'refmet_super_class'),
    Level('main_class', 'refmet_main_class'),
    Level('sub_class', 'refmet_sub_class'),
]
f = FieldConfig(extract={'chebi': '^(?:CHEBI:)?(\\d+)$'})
schema = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.REFMET, value=f('refmet_id')),
        CV(term=Namespace.PUBCHEM, value=f('pubchem_cid')),
        CV(term=Namespace.CHEBI, value=f('chebi_id', extract='chebi')),
        CV(term=Namespace.HMDB, value=f('hmdb_id')),
        CV(term=Namespace.LIPIDMAPS, value=f('lipidmaps_id')),
        CV(term=Namespace.KEGG, value=f('kegg_id')),
        CV(term=Namespace.INCHIKEY, value=f('inchi_key')),
        CV(term=Namespace.NAME, value=f('refmet_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_chemical_formula, value=f('formula')),
        CV(term='chemrof:monoisotopic_mass', value=f('exactmass')),
    ),
    associations=AssociationsBuilder(*class_associations(CLASS_LEVELS)),
)


def _raw(opener, **kwargs):
    for row in iter_csv(opener, **kwargs):
        yield {
            key.strip(): value for key, value in row.items() if key is not None
        }


resource = Resource(
    config=config,
    metabolites=Dataset(download=download, mapper=schema, raw_parser=_raw),
    terms=Dataset(
        download=download,
        mapper=map_class_term,
        raw_parser=lambda opener, max_records=None, **kwargs: iter_class_terms(
            _raw(opener, **kwargs), CLASS_LEVELS, max_records=max_records
        ),
    ),
)
