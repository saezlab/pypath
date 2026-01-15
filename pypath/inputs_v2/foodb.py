"""
Parse FooDB data and emit Entity records.

FooDB is the world's largest and most comprehensive resource on food 
constituents, chemistry and biology. This module creates Food entities 
with compound members, where each compound member includes full identifiers 
(SMILES, InChIKey, ChEBI, KEGG, CAS, etc.) and concentration annotations.

Data sources:
- https://foodb.ca/
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
from pypath.inputs_v2.parsers.foodb import _raw, MEMBER_DELIMITER


config = ResourceConfig(
    id=ResourceCv.FOODB,
    name='FooDB',
    url='https://foodb.ca/',
    license=LicenseCV.CC_BY_NC_4_0,  # From website terms
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='29059383',  # FooDB publication
    description=(
        'FooDB is the world\'s largest and most comprehensive resource on food '
        'constituents, chemistry and biology. It provides information on both '
        'macronutrients and micronutrients, including compounds that give foods '
        'their taste, colour, texture and aroma.'
    ),
)

# =============================================================================
# Download configurations
# =============================================================================

download_csv = Download(
    url='https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz',
    filename='foodb_2020_4_7_csv.tar.gz',
    subfolder='foodb',
    large=True,
    ext='tar',  # Actually a plain tar despite .tar.gz extension
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
        CV(term=IdentifierNamespaceCv.FOODB, value=f('public_id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('ncbi_taxonomy_id')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.SCIENTIFIC_NAME, value=f('name_scientific')),
        CV(term=MoleculeAnnotationsCv.FOOD_CLASS, value=f('food_group')),
        CV(term=MoleculeAnnotationsCv.FOOD_SUBCLASS, value=f('food_subgroup')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.SMALL_MOLECULE,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.FOODB, value=f('member_compound_public_id')),
                CV(term=IdentifierNamespaceCv.NAME, value=f('member_compound_name')),
                CV(term=IdentifierNamespaceCv.CAS, value=f('member_cas')),
                CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=f('member_inchikey')),
                CV(term=IdentifierNamespaceCv.SMILES, value=f('member_smiles')),
                CV(term=IdentifierNamespaceCv.CHEBI, value=f('member_chebi')),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('member_kegg')),
                CV(term=IdentifierNamespaceCv.IUPAC_NAME, value=f('member_iupac')),
                CV(term=IdentifierNamespaceCv.SYNONYM, value=f('member_synonyms')),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('member_mass')),
                CV(term=MoleculeAnnotationsCv.COMPOUND_CLASS, value=f('member_klass')),
                CV(term=MoleculeAnnotationsCv.COMPOUND_SUBCLASS, value=f('member_subklass')),
            ),
            annotations=AnnotationsBuilder(
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_MEAN, value=f('member_content')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_MIN, value=f('member_min')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_MAX, value=f('member_max')),
                CV(term=MoleculeAnnotationsCv.CONCENTRATION_UNIT, value=f('member_unit')),
                CV(term=MoleculeAnnotationsCv.EXPERIMENTAL_METHOD, value=f('member_method')),
            ),
        ),
    ),
)


# =============================================================================
# Resource definition
# =============================================================================

# Path to extracted FooDB data (extracted from tar archive)
FOODB_DATA_DIR = '/Users/jschaul/Code/omnipath_build/pypath-data/foodb/foodb_2020_04_07_csv'

resource = Resource(
    config,
    foods=Dataset(
        download=download_csv,
        mapper=foods_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            data_dir=FOODB_DATA_DIR,
            **kwargs,
        ),
    ),
)


__all__ = ['config', 'resource', 'download_csv']
