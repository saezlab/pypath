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
    primary_category='foods',
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
    needed=[
        'foodb_2020_04_07_csv/Food.csv',
        'foodb_2020_04_07_csv/Content.csv',
        'foodb_2020_04_07_csv/Compound.csv',
        'foodb_2020_04_07_csv/CompoundSynonym.csv',
        'foodb_2020_04_07_csv/CompoundExternalDescriptor.csv',
    ],
)


# =============================================================================
# Field configuration
# =============================================================================

f = FieldConfig(
    delimiter=MEMBER_DELIMITER,
    preserve_indices=True,
    extract={
        'taxid': r'^(-?\d+)',
        'cas': r'(\d{1,7}-\d{2}-\d)',
        'chebi': r'^(?:CHEBI:)?(\d+)$',
    },
)


# =============================================================================
# Food Schema with Compound Membership (declarative)
# =============================================================================

foods_schema = EntityBuilder(
    entity_type=EntityTypeCv.FOOD,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.FOODB, value=f('public_id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.SCIENTIFIC_NAME, value=f('name_scientific')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('ncbi_taxonomy_id', extract='taxid')),
        CV(term=MoleculeAnnotationsCv.DESCRIPTION, value=f('description')),
        CV(term=MoleculeAnnotationsCv.FOOD_CLASS, value=f('food_group')),
        CV(term=MoleculeAnnotationsCv.FOOD_SUBCLASS, value=f('food_subgroup')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.CHEMICAL,
            # Only the FOODB public_id linking key. The compound's full
            # identifier set (incl. InChIKey, SMILES, ChEBI, KEGG), synonyms and
            # class annotations are projected ONCE by the `compounds` dataset
            # below; canonicalize merges this membership compound into that
            # shared entity by FOODB id. Emitting any further identifier here
            # would re-project it for each of the ~4M food↔compound content rows
            # (the FooDB projection bottleneck) with no resolution benefit, since
            # cross-resource merging (e.g. FooDB↔HMDB by InChIKey) already
            # happens through the shared compound entity.
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.FOODB, value=f('member_compound_public_id')),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=MoleculeAnnotationsCv.CONCENTRATION_MEAN,
                    value=f('member_content'),
                    unit=f('member_unit'),
                ),
                CV(
                    term=MoleculeAnnotationsCv.CONCENTRATION_MIN,
                    value=f('member_min'),
                    unit=f('member_unit'),
                ),
                CV(
                    term=MoleculeAnnotationsCv.CONCENTRATION_MAX,
                    value=f('member_max'),
                    unit=f('member_unit'),
                ),
                CV(term=MoleculeAnnotationsCv.EXPERIMENTAL_METHOD, value=f('member_method')),
            ),
        ),
    ),
)


# =============================================================================
# Compound Schema — each unique FooDB compound projected ONCE (with synonyms)
# =============================================================================

compounds_schema = EntityBuilder(
    entity_type=EntityTypeCv.CHEMICAL,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.FOODB, value=f('public_id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.CAS, value=f('cas_number', extract='cas')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=f('moldb_inchikey')),
        CV(term=IdentifierNamespaceCv.SMILES, value=f('moldb_smiles')),
        CV(term=IdentifierNamespaceCv.CHEBI, value=f('chebi', extract='chebi')),
        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('kegg')),
        CV(term=IdentifierNamespaceCv.IUPAC_NAME, value=f('moldb_iupac')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('synonyms')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('moldb_mono_mass')),
        CV(term=MoleculeAnnotationsCv.COMPOUND_CLASS, value=f('klass')),
        CV(term=MoleculeAnnotationsCv.COMPOUND_SUBCLASS, value=f('subklass')),
    ),
)


# =============================================================================
# Resource definition
# =============================================================================

resource = Resource(
    config,
    compounds=Dataset(
        download=download_csv,
        mapper=compounds_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener, data_type='compounds', **kwargs
        ),
    ),
    foods=Dataset(
        download=download_csv,
        mapper=foods_schema,
        raw_parser=_raw,
    ),
)


__all__ = ['config', 'resource', 'download_csv']
