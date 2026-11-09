"""Parse PhytoHub data and emit Entity records.

PhytoHub is a freely available database containing detailed information about dietary phytochemicals and their human and animal metabolites.
It includes (poly)phenols, terpenoids, alkaloids and other plant secondary metabolites present in commonly consumed foods,
together with their chemical structures, physico-chemical properties, dietary sources, known and in silico predicted metabolites, and intake biomarkers.

Data sources:
- https://phytohub.eu/
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

from biolink_model.datamodel.model import (
    ChemicalEntity,
    Food,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, Topic

from pypath.inputs_v2._class_terms import (
    Level,
    class_associations,
    iter_class_terms,
    map_class_term,
    member_class_associations,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phytohub import (
    MEMBER_DELIMITER,
    _raw,
    _raw_compounds,
)
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
    id=ResourceCv.PHYTOHUB,
    name='PhytoHub',
    url='https://phytohub.eu/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.REGULAR,
    topics=(Topic.NUTRITION,),
    content=(Content.FOOD_COMPOSITION,),
    summary='Dietary phytochemicals, food sources, and their metabolites.',
    description=(
        'PhytoHub is a freely available database containing detailed information about dietary phytochemicals and their human and animal metabolites. '
    ),
)

download_compounds = Download(
    url='https://phytohub.eu/entries',
    filename='PhytoHub Compounds.html',
    subfolder='phytohub',
    large=True,
    ext='html',
)
download_foods = Download(
    url='https://phytohub.eu/food_sources',
    filename='PhytoHub Foods.html',
    subfolder='phytohub',
    large=True,
    ext='html',
)

f = FieldConfig(
    delimiter=MEMBER_DELIMITER,
    preserve_indices=True,
    extract={
        'cas': r'(\d{1,7}-\d{2}-\d)',
    },
)

FOOD_LEVELS = [
    Level('food_group', 'phytohub_food_group'),
    Level('food_subgroup', 'phytohub_food_subgroup'),
]
PHYTOCHEMICAL_LEVELS = [
    Level('Phytochemical family', 'phytohub_phytochemical_family'),
    Level('Phytochemical class', 'phytohub_phytochemical_class'),
    Level('Phytochemical subclass', 'phytohub_phytochemical_subclass'),
]
METABOLITE_LEVELS = [
    Level('Metabolite family', 'phytohub_metabolite_family'),
    Level('Metabolite class', 'phytohub_metabolite_class'),
    Level('Metabolite subclass', 'phytohub_metabolite_subclass'),
]
MEMBER_CLASS_LEVELS = [
    Level('member_phytochemical_family', 'phytohub_phytochemical_family'),
    Level('member_phytochemical_class', 'phytohub_phytochemical_class'),
    Level('member_phytochemical_subclass', 'phytohub_phytochemical_subclass'),
    Level('member_metabolite_family', 'phytohub_metabolite_family'),
    Level('member_metabolite_class', 'phytohub_metabolite_class'),
    Level('member_metabolite_subclass', 'phytohub_metabolite_subclass'),
]
CLASS_LEVELS = [*FOOD_LEVELS, *PHYTOCHEMICAL_LEVELS, *METABOLITE_LEVELS]


foods_schema = EntityBuilder(
    entity_type=Food,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.PHYTOHUB_FOOD, value=f('id')),
        CV(term=Namespace.NAME, value=f('name')),
        CV(term=Namespace.SYNONYM, value=f('synonyms')),
    ),
    annotations=AnnotationsBuilder(
        CV(term='phytohub:intake_biomarkers', value=f('intake_biomarkers')),
    ),
    associations=AssociationsBuilder(
        *class_associations(FOOD_LEVELS),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=ChemicalEntity,
            identifiers=IdentifiersBuilder(
                CV(term=Namespace.PHYTOHUB_COMPOUND, value=f('member_compound_id')),
                CV(term=Namespace.NAME, value=f('member_compound_name')),
                CV(term=Namespace.SYNONYM, value=f('member_synonyms')),
                CV(term=Namespace.SYNONYM, value=f('member_systematic_name')),
                CV(term=Namespace.SYNONYM, value=f('member_iupac_name')),
                CV(term=Namespace.CAS, value=f('member_cas', extract='cas')),
                CV(term=Namespace.INCHIKEY, value=f('member_inchi_key')),
                CV(term=Namespace.SMILES, value=f('member_smiles')),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=slots.has_chemical_formula, value=f('member_formula')),
                CV(term='chemrof:mass', value=f('member_average_mass')),
                CV(term='chemrof:monoisotopic_mass', value=f('member_monoisotopic_mass')),
                CV(term='phytohub:solubility', value=f('member_solubility')),
                CV(term='phytohub:predicted_logp', value=f('member_predicted_logp')),
                CV(term='phytohub:predicted_logs', value=f('member_predicted_logs')),
                CV(term='phytohub:pka_basic', value=f('member_pka_basic')),
                CV(term='phytohub:pka_acidic', value=f('member_pka_acidic')),
            ),
            entity_associations=AssociationsBuilder(
                *member_class_associations(MEMBER_CLASS_LEVELS, MEMBER_DELIMITER),
            ),
            predicate=slots.has_part,
        )
    ),
)

compounds_schema = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.PHYTOHUB_COMPOUND, value=f('Identifier')),
        CV(term=Namespace.NAME, value=f('Chemical Name')),
        CV(term=Namespace.SYNONYM, value=f('Synonyms')),
        CV(term=Namespace.SYNONYM, value=f('Systematic name')),
        CV(term=Namespace.SYNONYM, value=f('IUPAC Name')),
        CV(term=Namespace.CAS, value=f('CAS Number', extract='cas')),
        CV(term=Namespace.INCHIKEY, value=f('InChI Key')),
        CV(term=Namespace.SMILES, value=f('Smiles')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_chemical_formula, value=f('Formula')),
        CV(term='chemrof:mass', value=f('Average mass')),
        CV(term='chemrof:monoisotopic_mass', value=f('Monoisotopic mass')),
        CV(term='phytohub:solubility', value=f('Solubility')),
        CV(term='phytohub:predicted_logp', value=f('Predicted alogps logp')),
        CV(term='phytohub:predicted_logs', value=f('Predicted alogps logs')),
        CV(term='phytohub:pka_basic', value=f('Pka strongest basic')),
        CV(term='phytohub:pka_acidic', value=f('Pka strongest acidic')),
        CV(term='phytohub:metabolites_precursor', value=f('Metabolites / Precursor')),
        CV(term='phytohub:intake_biomarkers', value=f('Intake biomarkers')),
    ),
    associations=AssociationsBuilder(
        *class_associations([*PHYTOCHEMICAL_LEVELS, *METABOLITE_LEVELS]),
    ),
)


def _classification_rows(
    opener: object,
    download_compounds: Download,
    **kwargs: object,
) -> Iterator[dict[str, Any]]:
    for row in _raw(opener, download_compounds=download_compounds, **kwargs):
        yield {level.field: row.get(level.field) for level in FOOD_LEVELS}

    for row in _raw_compounds(download_compounds.open(), **kwargs):
        for levels in (PHYTOCHEMICAL_LEVELS, METABOLITE_LEVELS):
            yield {level.field: row.get(level.field) for level in levels}


resource = Resource(
    config,
    foods=Dataset(
        download=download_foods,
        mapper=foods_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener,
            download_compounds=download_compounds,
            **kwargs,
        ),
    ),
    compounds=Dataset(
        download=download_compounds,
        mapper=compounds_schema,
        raw_parser=lambda opener, **kwargs: _raw_compounds(opener, **kwargs),
    ),
    terms=Dataset(
        download=download_foods,
        mapper=map_class_term,
        raw_parser=lambda opener, max_records=None, **kwargs: iter_class_terms(
            _classification_rows(
                opener,
                download_compounds=download_compounds,
                **kwargs,
            ),
            CLASS_LEVELS,
            max_records=max_records,
        ),
    ),
)

__all__ = ['config', 'resource', 'download_compounds', 'download_foods']
