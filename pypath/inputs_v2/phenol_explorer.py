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

import math
import re

from biolink_model.datamodel.model import (
    ChemicalEntity,
    Food,
    QuantityValue,
    slots,
)
from omnipath_core.measurements import Measurement
from omnipath_core.naming import Namespace

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phenol_explorer import MEMBER_DELIMITER, _raw
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


def _measurement(value, unit=None, source_field=None, comparator=None):
    """Keep numeric source observations, units and comparison bounds together."""
    if value is None:
        return None
    match = re.fullmatch(
        '\\s*(<=|>=|<|>|=|~)?\\s*([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)\\s*',
        str(value),
    )
    if match is None:
        return None
    original_comparator = comparator or match.group(1)
    if original_comparator not in {None, '<', '>', '=', '<=', '>=', '~', '≈'}:
        return None
    relation = {'<': 'less_than', '>': 'greater_than', '=': 'equal_to'}.get(
        original_comparator
    )
    number = float(match.group(2))
    if not math.isfinite(number):
        return None
    quantity = QuantityValue(
        has_numeric_value=number,
        has_unit=unit or None,
        has_binary_relation=relation,
    )
    return Measurement(
        quantity=quantity,
        source_field=source_field,
        comparator=original_comparator,
    )


config = ResourceConfig(
    id=ResourceCv.PHENOL_EXPLORER,
    name='Phenol-Explorer',
    url='http://phenol-explorer.eu/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='24103452',
    primary_category='foods',
    description='Phenol-Explorer is the first comprehensive database on polyphenol content in foods. It contains data on the content and composition of polyphenols and other bioactive compounds in foods, with detailed classification of both compounds and food sources.',
)
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
f = FieldConfig(
    delimiter=MEMBER_DELIMITER,
    preserve_indices=True,
    extract={'chebi': '^(?:CHEBI:)?(\\d+)$'},
)
foods_schema = EntityBuilder(
    entity_type=Food,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.PHENOL_EXPLORER_FOOD, value=f('id')),
        CV(term=Namespace.NAME, value=f('name')),
        CV(term=Namespace.SYNONYM, value=f('scientific_name')),
    ),
    annotations=AnnotationsBuilder(),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.PHENOL_EXPLORER_COMPOUND,
                    value=f('member_compound_id'),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f('member_chebi', extract='chebi'),
                ),
                CV(term=Namespace.PUBCHEM, value=f('member_pubchem')),
                CV(term=Namespace.CAS, value=f('member_cas')),
                CV(term=Namespace.SMILES, value=f('member_smiles')),
                CV(term=Namespace.NAME, value=f('member_compound_name')),
                CV(term=Namespace.SYNONYM, value=f('member_synonyms')),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=slots.has_chemical_formula, value=f('member_formula')),
                CV(term='chemrof:mass', value=f('member_molecular_weight')),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.publications,
                    value=f(
                        'member_pubmed',
                        transform=lambda v: [
                            'PMID:' + p.strip().removeprefix('PMID:')
                            for p in str(v).split(';')
                            if p.strip().removeprefix('PMID:').isdigit()
                        ],
                    ),
                ),
                CV(
                    term='PATO:0000033',
                    value=lambda row: [
                        _measurement(
                            v,
                            units[i] if i < len(units) else None,
                            'member_mean',
                        )
                        for units in [
                            str(row.get('member_units') or '').split(
                                MEMBER_DELIMITER
                            )
                        ]
                        for i, v in enumerate(
                            str(row.get('member_mean') or '').split(
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
                            str(row.get('member_units') or '').split(
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
                            str(row.get('member_units') or '').split(
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
                CV(
                    term='PATO:0000033',
                    value=lambda row: [
                        _measurement(
                            v, units[i] if i < len(units) else None, 'member_sd'
                        )
                        for units in [
                            str(row.get('member_units') or '').split(
                                MEMBER_DELIMITER
                            )
                        ]
                        for i, v in enumerate(
                            str(row.get('member_sd') or '').split(
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
