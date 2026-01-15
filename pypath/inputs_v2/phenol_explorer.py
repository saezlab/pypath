"""
Parse Phenol-Explorer data and emit Entity records.

Phenol-Explorer is a comprehensive database on polyphenol content in foods.
This module creates Food entities with polyphenol compounds as members,
where each compound member includes full identifiers (SMILES, CHEBI, etc.)
and concentration/measurement annotations.

Data sources:
- http://phenol-explorer.eu/
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    MoleculeAnnotationsCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembershipBuilder,
    MembersFromList,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phenol_explorer import _raw, MEMBER_DELIMITER


config = ResourceConfig(
    id=ResourceCv.PHENOL_EXPLORER,
    name='Phenol-Explorer',
    url='http://phenol-explorer.eu/',
    license=LicenseCV.CC_BY_NC_4_0,  # TODO: Verify license
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='20428313',  # Original Phenol-Explorer publication
    description=(
        'Phenol-Explorer is the first comprehensive database on polyphenol '
        'content in foods. It contains data on the content and composition of '
        'polyphenols and other bioactive compounds in foods, with detailed '
        'classification of both compounds and food sources.'
    ),
)

# =============================================================================
# Download configurations
# =============================================================================

download_compounds = Download(
    url='http://phenol-explorer.eu/system/downloads/current/compounds.csv.zip',
    filename='compounds.csv.zip',
    subfolder='phenol_explorer',
    large=True,
    ext='zip',
)

download_compounds_structures = Download(
    url='http://phenol-explorer.eu/system/downloads/current/compounds-structures.csv.zip',
    filename='compounds-structures.csv.zip',
    subfolder='phenol_explorer',
    large=True,
    ext='zip',
)

download_foods = Download(
    url='http://phenol-explorer.eu/system/downloads/current/foods.csv.zip',
    filename='foods.csv.zip',
    subfolder='phenol_explorer',
    large=True,
    ext='zip',
)

download_composition = Download(
    url='http://phenol-explorer.eu/system/downloads/current/composition-data.xlsx.zip',
    filename='composition-data.xlsx.zip',
    subfolder='phenol_explorer',
    large=True,
    ext='zip',
    default_mode='rb',
)


# =============================================================================
# Field configuration
# =============================================================================

f = FieldConfig(delimiter=MEMBER_DELIMITER, preserve_indices=True)


# =============================================================================
# Food Schema with Compound Membership (declarative)
# =============================================================================

foods_schema = EntityBuilder(
    entity_type=EntityTypeCv.FOOD,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.PHENOL_EXPLORER, value=f('id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.FOOD_CLASS, value=f('food_group')),
        CV(term=MoleculeAnnotationsCv.FOOD_SUBCLASS, value=f('food_subgroup')),
        CV(term=MoleculeAnnotationsCv.SCIENTIFIC_NAME, value=f('scientific_name')),
        CV(term=MoleculeAnnotationsCv.BOTANICAL_FAMILY, value=f('botanical_family')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.SMALL_MOLECULE,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.PHENOL_EXPLORER, value=f('member_compound_id')),
                CV(term=IdentifierNamespaceCv.NAME, value=f('member_compound_name')),
                CV(term=IdentifierNamespaceCv.CHEBI, value=f('member_chebi')),
                CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('member_pubchem')),
                CV(term=IdentifierNamespaceCv.CAS, value=f('member_cas')),
                CV(term=IdentifierNamespaceCv.SMILES, value=f('member_smiles')),
                CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA, value=f('member_formula')),
                CV(term=IdentifierNamespaceCv.SYNONYM, value=f('member_synonyms')),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=MoleculeAnnotationsCv.COMPOUND_CLASS, value=f('member_compound_class')),
                CV(term=MoleculeAnnotationsCv.COMPOUND_SUBCLASS, value=f('member_compound_subclass')),
                CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('member_molecular_weight')),
                CV(term=MoleculeAnnotationsCv.AGLYCONE, value=f('member_aglycones')),
            ),
            annotations=AnnotationsBuilder(
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_MEAN, value=f('member_mean')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_MIN, value=f('member_min')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_MAX, value=f('member_max')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_SD, value=f('member_sd')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_UNIT, value=f('member_units')),
                CV(term=MoleculeAnnotationsCv.SAMPLE_COUNT, value=f('member_n')),
                CV(term=MoleculeAnnotationsCv.DATA_POINT_COUNT, value=f('member_N')),
                CV(term=MoleculeAnnotationsCv.EXPERIMENTAL_METHOD, value=f('member_experimental_method')),
                CV(term=IdentifierNamespaceCv.PUBMED, value=f('member_pubmed')),
            ),
        ),
    ),
)


# =============================================================================
# Resource definition
# =============================================================================

resource = Resource(
    config,
    foods=Dataset(
        download=download_foods,
        mapper=foods_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener,
            download_compounds=download_compounds,
            download_compounds_structures=download_compounds_structures,
            download_composition=download_composition,
            **kwargs,
        ),
    ),
)


__all__ = ['config', 'resource']
