"""
Parse Reactome BioPAX data and emit Entity records.

This module converts Reactome BioPAX data into Entity records using the
declarative schema pattern.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    ParticipantMetadataCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.reactome import _raw


config = ResourceConfig(
    id=ResourceCv.REACTOME,
    name='Reactome',
    url='https://reactome.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='34788843',
    description=(
        'Reactome is a free, open-source, curated and peer-reviewed '
        'pathway database.'
    ),
)


# --------------------------------------------------------------------------- #
# Processing Helpers
# --------------------------------------------------------------------------- #

_MISSING_VALUE = "__MISSING__"

entity_type_map = {v.value: v for v in EntityTypeCv}
role_map = {
    'reactant': BiologicalRoleCv.REACTANT,
    'product': BiologicalRoleCv.PRODUCT,
    'template': BiologicalRoleCv.TEMPLATE,
    'controller': BiologicalRoleCv.CONTROLLER,
    'controlled': BiologicalRoleCv.CONTROLLED,
}

f = FieldConfig(
    delimiter=';',
    map={
        'entity_type': lambda value: (
            EntityTypeCv.PHYSICAL_ENTITY
            if not value or value == _MISSING_VALUE
            else entity_type_map.get(value, EntityTypeCv.PHYSICAL_ENTITY)
        ),
        'role': lambda value: (
            None
            if not value or value == _MISSING_VALUE
            else role_map.get(value)
        ),
        'missing': lambda value: '' if not value or value == _MISSING_VALUE else value,
        'split': lambda value: [] if not value or value == _MISSING_VALUE else [
            item for item in value.split(';') if item
        ],
    },
)


reactions_schema = EntityBuilder(
    entity_type=f('entity_type', map=entity_type_map),
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('display_name')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('synonyms')),
        CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('reactome_stable_id')),
        CV(term=IdentifierNamespaceCv.REACTOME_ID, value=f('reactome_id')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed')),
        CV(term=MoleculeAnnotationsCv.EC_NUMBER, value=f('ec_number')),
        CV(term=InteractionMetadataCv.CONVERSION_DIRECTION, value=f('direction')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=f('participant_entity_type', delimiter='||', map='entity_type'),
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.NAME, value=f('participant_display_name', delimiter='||', map='missing')),
                CV(term=IdentifierNamespaceCv.SYNONYM, value=f('participant_synonyms', delimiter='||', map='split')),
                CV(
                    term=IdentifierNamespaceCv.REACTOME_STABLE_ID,
                    value=f('participant_reactome_stable_id', delimiter='||', map='split'),
                ),
                CV(term=IdentifierNamespaceCv.UNIPROT, value=f('participant_uniprot', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.CHEBI, value=f('participant_chebi', delimiter='||', map='split')),
                CV(
                    term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                    value=f('participant_pubchem_compound', delimiter='||', map='split'),
                ),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('participant_kegg', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('participant_go', delimiter='||', map='split')),
            ),
            annotations=AnnotationsBuilder(
                CV(term=f('participant_role', delimiter='||', map='role')),
                CV(
                    term=ParticipantMetadataCv.STOICHIOMETRY,
                    value=f('participant_stoichiometry', delimiter='||', map='missing'),
                ),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term=IdentifierNamespaceCv.NCBI_TAX_ID,
                    value=f('participant_ncbi_tax_id', delimiter='||', map='missing'),
                ),
            ),
        ),
    ),
)


pathways_schema = EntityBuilder(
    entity_type=EntityTypeCv.PATHWAY,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('display_name')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('synonyms')),
        CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('reactome_stable_id')),
        CV(term=IdentifierNamespaceCv.REACTOME_ID, value=f('reactome_id')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('go')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('ncbi_tax_id')),
        CV(term=MoleculeAnnotationsCv.DESCRIPTION, value=f('description')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=f('component_entity_type', delimiter='||', map='entity_type'),
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.NAME, value=f('component_display_name', delimiter='||', map='missing')),
                CV(
                    term=IdentifierNamespaceCv.REACTOME_STABLE_ID,
                    value=f('component_reactome_stable_id', delimiter='||', map='split'),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.PATHWAY_COMPONENT),
                CV(
                    term=ParticipantMetadataCv.STEP_ORDER,
                    value=f('component_step_order', delimiter='||', map='missing'),
                ),
            ),
        ),
    ),
)


controls_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('display_name')),
        CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('reactome_stable_id')),
        CV(term=IdentifierNamespaceCv.REACTOME_ID, value=f('reactome_id')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('go')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.CONTROL_TYPE, value=f('control_type')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=f('participant_entity_type', delimiter='||', map='entity_type'),
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.NAME, value=f('participant_display_name', delimiter='||', map='missing')),
                CV(term=IdentifierNamespaceCv.SYNONYM, value=f('participant_synonyms', delimiter='||', map='split')),
                CV(
                    term=IdentifierNamespaceCv.REACTOME_STABLE_ID,
                    value=f('participant_reactome_stable_id', delimiter='||', map='split'),
                ),
                CV(term=IdentifierNamespaceCv.UNIPROT, value=f('participant_uniprot', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.CHEBI, value=f('participant_chebi', delimiter='||', map='split')),
                CV(
                    term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                    value=f('participant_pubchem_compound', delimiter='||', map='split'),
                ),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('participant_kegg', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('participant_go', delimiter='||', map='split')),
            ),
            annotations=AnnotationsBuilder(
                CV(term=f('participant_role', delimiter='||', map='role')),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term=IdentifierNamespaceCv.NCBI_TAX_ID,
                    value=f('participant_ncbi_tax_id', delimiter='||', map='missing'),
                ),
            ),
        ),
    ),
)


download=Download(
            url='https://reactome.org/download/current/biopax.zip',
            filename='reactome_biopax.zip',
            subfolder='reactome',
            large=True,
            ext='zip',
            default_mode='rb',
        )

resource = Resource(
    config,
    reactions=Dataset(
        download=download,
        mapper=reactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='reactions', **kwargs),
    ),
    pathways=Dataset(
        download=download,
        mapper=pathways_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='pathways', **kwargs),
    ),
    controls=Dataset(
        download=download,
        mapper=controls_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='controls', **kwargs),
    ),
)
