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

from biolink_model.datamodel.model import (
    ChemicalEntity,
    Food,
    slots,
)
from omnipath_core.naming import Namespace

from pypath.inputs_v2._measurements import measurement as _measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.foodb import MEMBER_DELIMITER, _raw
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
)

config = ResourceConfig(
    id=ResourceCv.FOODB,
    name='FooDB',
    url='https://foodb.ca/',
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    primary_category='foods',
    description="FooDB is the world's largest and most comprehensive resource on food constituents, chemistry and biology. It provides information on both macronutrients and micronutrients, including compounds that give foods their taste, colour, texture and aroma.",
)
download_csv = Download(
    url='https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz',
    filename='foodb_2020_4_7_csv.tar.gz',
    subfolder='foodb',
    large=True,
    ext='tar',
    needed=[
        'foodb_2020_04_07_csv/Food.csv',
        'foodb_2020_04_07_csv/Content.csv',
        'foodb_2020_04_07_csv/Compound.csv',
        'foodb_2020_04_07_csv/CompoundSynonym.csv',
        'foodb_2020_04_07_csv/CompoundExternalDescriptor.csv',
    ],
)
f = FieldConfig(
    delimiter=MEMBER_DELIMITER,
    preserve_indices=True,
    extract={
        'taxid': '^([1-9]\\d*)$',
        'cas': '(\\d{1,7}-\\d{2}-\\d)',
        'chebi': '^(?:CHEBI:)?(\\d+)$',
    },
)
foods_schema = EntityBuilder(
    entity_type=Food,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.FOODB, value=f('public_id')),
        CV(term=Namespace.NAME, value=f('name')),
        CV(term=Namespace.SYNONYM, value=f('name_scientific')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=f(
                'ncbi_taxonomy_id',
                extract='taxid',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:'),
            ),
        ),
        CV(term=slots.description, value=f('description')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.FOODB, value=f('member_compound_public_id'))
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term='PATO:0000033',
                    value=lambda row: [
                        _measurement(
                            v,
                            units[i] if i < len(units) else None,
                            'member_content',
                        )
                        for units in [
                            str(row.get('member_unit') or '').split(
                                MEMBER_DELIMITER
                            )
                        ]
                        for i, v in enumerate(
                            str(row.get('member_content') or '').split(
                                MEMBER_DELIMITER
                            )
                        )
                    ],
                ),
                CV(
                    term='PATO:0000033',
                    value=lambda row: [
                        _measurement(
                            v,
                            units[i] if i < len(units) else None,
                            'member_min',
                        )
                        for units in [
                            str(row.get('member_unit') or '').split(
                                MEMBER_DELIMITER
                            )
                        ]
                        for i, v in enumerate(
                            str(row.get('member_min') or '').split(
                                MEMBER_DELIMITER
                            )
                        )
                    ],
                ),
                CV(
                    term='PATO:0000033',
                    value=lambda row: [
                        _measurement(
                            v,
                            units[i] if i < len(units) else None,
                            'member_max',
                        )
                        for units in [
                            str(row.get('member_unit') or '').split(
                                MEMBER_DELIMITER
                            )
                        ]
                        for i, v in enumerate(
                            str(row.get('member_max') or '').split(
                                MEMBER_DELIMITER
                            )
                        )
                    ],
                ),
            ),
            predicate=slots.has_part,
        )
    ),
)
compounds_schema = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.FOODB, value=f('public_id')),
        CV(term=Namespace.CAS, value=f('cas_number', extract='cas')),
        CV(term=Namespace.INCHIKEY, value=f('moldb_inchikey')),
        CV(term=Namespace.SMILES, value=f('moldb_smiles')),
        CV(term=Namespace.CHEBI, value=f('chebi', extract='chebi')),
        CV(term=Namespace.KEGG, value=f('kegg')),
        CV(term=Namespace.NAME, value=f('name')),
        CV(term=Namespace.SYNONYM, value=f('moldb_iupac')),
        CV(term=Namespace.SYNONYM, value=f('synonyms')),
    ),
    annotations=AnnotationsBuilder(
        CV(term='chemrof:monoisotopic_mass', value=f('moldb_mono_mass'))
    ),
)
resource = Resource(
    config,
    compounds=Dataset(
        download=download_csv,
        mapper=compounds_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener, data_type='compounds', **kwargs
        ),
    ),
    foods=Dataset(download=download_csv, mapper=foods_schema, raw_parser=_raw),
)
__all__ = ['config', 'resource', 'download_csv']
