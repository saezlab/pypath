"""
Parse Phenol-Explorer data and emit Entity records.

Phenol-Explorer is a comprehensive database on polyphenol content in foods.
This module creates Food entities with polyphenol compounds as members,
where each compound member includes full identifiers (SMILES, CHEBI, etc.)
and concentration/measurement annotations.

The food is associated with its food group and subgroup and its botanical
family; each compound with its class, subclass and aglycones, terms of
Phenol-Explorer's classification (``terms``). A measurement keeps its number
of samples (``supporting_study_size``), of data points and its method.

Data sources:
- http://phenol-explorer.eu/
"""

from __future__ import annotations

from biolink_model.datamodel.model import (
    ChemicalEntity,
    Food,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._class_terms import (
    Level,
    class_associations,
    iter_class_terms,
    map_class_term,
    member_class_associations,
    member_rows,
    plain,
)
from pypath.inputs_v2._measurements import measurement as _measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phenol_explorer import MEMBER_DELIMITER, _raw
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    AssociationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
)

config = ResourceConfig(
    id=ResourceCv.PHENOL_EXPLORER,
    name='Phenol-Explorer',
    url='http://phenol-explorer.eu/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='24103452',
    topics=(Topic.NUTRITION,),
    content=(Content.FOOD_COMPOSITION,),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL,),
    summary='Polyphenol constituents of foods with identifiers and concentration measurements.',
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
def _no_dash(value):
    """Phenol-Explorer writes "-" for a compound without an aglycone."""
    return None if str(value or '').strip() == '-' else plain(value)


FOOD_LEVELS = [
    Level('food_group', 'phenol_explorer_food_group'),
    Level('food_subgroup', 'phenol_explorer_food_subgroup'),
]
BOTANICAL_FAMILY = [Level('botanical_family', 'phenol_explorer_botanical_family')]
COMPOUND_LEVELS = [
    Level('member_compound_class', 'phenol_explorer_compound_class'),
    Level('member_compound_subclass', 'phenol_explorer_compound_subclass'),
]
# A compound's aglycones ("Quercetin", "(-)-Epicatechin // Gallic acid") group
# its glycosides; they are terms of their own, not placed under a subclass.
AGLYCONE = [Level('member_aglycones', 'phenol_explorer_aglycone', _no_dash)]
AGLYCONE_SEPARATOR = ' // '


def _member_field(field):
    """A per-member field, one value per member."""
    return lambda row: str(row.get(field) or '').split(MEMBER_DELIMITER)


def _terms(opener, max_records=None, **kwargs):
    rows = list(_raw(
        opener,
        download_compounds=download_compounds,
        download_compounds_structures=download_compounds_structures,
        download_composition=download_composition,
        **kwargs,
    ))
    yield from iter_class_terms(rows, FOOD_LEVELS)
    yield from iter_class_terms(rows, BOTANICAL_FAMILY)
    yield from iter_class_terms(member_rows(rows, [l.field for l in COMPOUND_LEVELS]), COMPOUND_LEVELS)
    yield from iter_class_terms(
        member_rows(rows, [AGLYCONE[0].field], split=AGLYCONE_SEPARATOR), AGLYCONE
    )


foods_schema = EntityBuilder(
    entity_type=Food,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.PHENOL_EXPLORER_FOOD, value=f('id')),
        CV(term=Namespace.NAME, value=f('name')),
        CV(term=Namespace.SYNONYM, value=f('scientific_name')),
    ),
    annotations=AnnotationsBuilder(),
    associations=AssociationsBuilder(*class_associations(FOOD_LEVELS + BOTANICAL_FAMILY)),
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
            entity_associations=AssociationsBuilder(
                *member_class_associations(COMPOUND_LEVELS, MEMBER_DELIMITER),
                *member_class_associations(AGLYCONE, MEMBER_DELIMITER, AGLYCONE_SEPARATOR),
            ),
            annotations=AnnotationsBuilder(
                # n: samples, N: data points the mean is computed from.
                CV(term=slots.supporting_study_size, value=_member_field('member_n')),
                CV(term='phenol_explorer:data_point_count', value=_member_field('member_N')),
                CV(
                    term=slots.supporting_study_method_description,
                    value=_member_field('member_experimental_method'),
                ),
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
    terms=Dataset(download=download_foods, mapper=map_class_term, raw_parser=_terms),
)
__all__ = ['config', 'resource']
