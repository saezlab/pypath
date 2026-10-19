"""
Parse LIPID MAPS Structure Database (LMSD) data and emit Entity records.

This module downloads the LIPID MAPS Structure Database in SDF format and
converts lipid structures into Entity records using the schema defined in
pypath.internals.silver_schema.
"""

from __future__ import annotations

from biolink_model.datamodel.model import ChemicalEntity, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, Topic

from pypath.inputs_v2._class_terms import (
    Level,
    class_associations,
    coded,
    iter_class_terms,
    map_class_term,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.lipidmaps import _raw
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    AssociationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)

config = ResourceConfig(
    id=ResourceCv.LIPIDMAPS,
    name='LIPID MAPS',
    full='LIPID MAPS Structure Database',
    url='https://lipidmaps.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='37855672',
    topics=(Topic.METABOLISM,),
    content=(Content.MOLECULAR_ANNOTATIONS, Content.IDENTIFIER_MAPPINGS),
    summary='Lipid structures, identifiers and chemical annotations from LIPID MAPS.',
    description='The LIPID MAPS Structure Database (LMSD) is a comprehensive database of lipid structures, annotations, and cross-references. It contains over 47,000 unique lipid structures, classified according to a comprehensive lipid classification system. The database includes structure-based identifiers (InChI, InChIKey, SMILES), chemical properties (formula, exact mass), lipid classification (category, main class, sub class), and cross-references to other databases (ChEBI, PubChem, SwissLipids, HMDB).',
)
f = FieldConfig(
    extract={'chebi': '^(?:CHEBI:)?(\\d+)$', 'hmdb': '^(HMDB\\d{5,8})$'},
    transform={'hmdb': lambda v: v.upper() if v else None},
)
# The LIPID MAPS classification, written "Fatty Acyls [FA]", "Fatty Acids and
# Conjugates [FA01]", "Straight chain fatty acids [FA0101]": the codes are
# unique across the levels, so the levels share one namespace.
CLASS_NAMESPACE = 'lipidmaps_class'
CLASS_LEVELS = [
    Level('CATEGORY', CLASS_NAMESPACE, coded),
    Level('MAIN_CLASS', CLASS_NAMESPACE, coded),
    Level('SUB_CLASS', CLASS_NAMESPACE, coded),
]
lipids_schema = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.LIPIDMAPS, value=f('LM_ID')),
        CV(term=Namespace.INCHIKEY, value=f('INCHI_KEY')),
        CV(term=Namespace.INCHI, value=f('INCHI')),
        CV(term=Namespace.SMILES, value=f('SMILES')),
        CV(term=Namespace.CHEBI, value=f('CHEBI_ID', extract='chebi')),
        CV(term=Namespace.PUBCHEM, value=f('PUBCHEM_CID')),
        CV(
            term=Namespace.HMDB,
            value=f('HMDB_ID', extract='hmdb', transform='hmdb'),
        ),
        CV(term=Namespace.SWISSLIPIDS, value=f('SWISSLIPIDS_ID')),
        CV(term=Namespace.NAME, value=f('COMMON_NAME')),
        CV(term=Namespace.SYNONYM, value=f('SYSTEMATIC_NAME')),
        CV(term=Namespace.SYNONYM, value=f('ABBREVIATION')),
        CV(term=Namespace.SYNONYM, value=f('SYNONYMS', delimiter=';')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_chemical_formula, value=f('FORMULA')),
        CV(term='chemrof:monoisotopic_mass', value=f('EXACT_MASS')),
    ),
    associations=AssociationsBuilder(*class_associations(CLASS_LEVELS)),
)
download = Download(
    url='https://lipidmaps.org/files/?file=LMSD&ext=sdf.zip',
    filename='structures.zip',
    subfolder='lipidmaps',
    large=True,
    ext='zip',
    default_mode='rb',
)


def _id_translation_row(row: dict) -> dict | None:
    lipidmaps_id = row.get('LM_ID')
    standard_inchi = row.get('INCHI')
    if not lipidmaps_id or not standard_inchi:
        return None
    return {
        'source': 'lipidmaps',
        'key_type': Namespace.LIPIDMAPS,
        'key_value': lipidmaps_id,
        'standard_inchi': standard_inchi,
    }


resource = Resource(
    config,
    lipids=Dataset(download=download, mapper=lipids_schema, raw_parser=_raw),
    terms=Dataset(
        download=download,
        mapper=map_class_term,
        raw_parser=lambda opener, max_records=None, **kwargs: iter_class_terms(
            _raw(opener, **kwargs), CLASS_LEVELS, max_records=max_records
        ),
    ),
    id_translation=Dataset(
        download=download,
        mapper=lambda row: row,
        raw_parser=lambda opener, **kwargs: (
            row
            for raw_row in _raw(opener, **kwargs)
            if (row := _id_translation_row(raw_row)) is not None
        ),
        kind='id_translation',
    ),
)
