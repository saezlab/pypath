"""
Parse Reactome BioPAX data and emit Entity records.

This module converts Reactome BioPAX data into Entity records using the
 declarative schema pattern.
"""

from __future__ import annotations

from functools import partial
import re

from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    ControlEffectCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    PathwayAnnotationsCv,
    ReactionAnnotationsCv,
    OntologyCv,
    ParticipantMetadataCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import EntityRef, OntologyRelation
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembersFromList,
    MembershipBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.reactome import _raw


config = ResourceConfig(
    id=ResourceCv.REACTOME,
    name='Reactome',
    url='https://reactome.org/',
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='34788843',
    primary_category='pathways',
    annotation_ontologies=(OntologyCv.GENE_ONTOLOGY, OntologyCv.REACTOME_PATHWAYS),
    description=(
        'Reactome is a free, open-source, curated and peer-reviewed '
        'pathway database.'
    ),
)


_MISSING_VALUE = '__MISSING__'

entity_type_map = {v.value: v for v in EntityTypeCv}
role_map = {
    'reactant': BiologicalRoleCv.REACTANT,
    'product': BiologicalRoleCv.PRODUCT,
    'template': BiologicalRoleCv.TEMPLATE,
    'controller': BiologicalRoleCv.CONTROLLER,
    'controlled': BiologicalRoleCv.CONTROLLED,
}


def _control_effect(row):
    control_type = (row.get('control_type') or '').upper()
    control_class = (row.get('control_class') or '').upper()
    if 'CATALYSIS' in control_class or 'CATALYSIS' in control_type:
        return ControlEffectCv.CATALYSIS
    if 'ACTIV' in control_type or 'POSITIVE' in control_type:
        return ControlEffectCv.ACTIVATION
    if 'INHIB' in control_type or 'NEGATIVE' in control_type:
        return ControlEffectCv.INHIBITION
    if control_type or control_class:
        return ControlEffectCv.MODULATION
    return ControlEffectCv.UNKNOWN


_TAXON_SCOPED_ENTITY_TYPES = {
    EntityTypeCv.PROTEIN,
    EntityTypeCv.GENE,
    EntityTypeCv.RNA,
    EntityTypeCv.DNA,
}
_REACTOME_PATHWAY_ONTOLOGY_ID = 'reactome_pathways'


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
        'split_chebi': lambda value: [] if not value or value == _MISSING_VALUE else [
            match.group(1)
            for item in value.split(';')
            if (match := re.fullmatch(r'(?:CHEBI:)?(\d+)', item.strip()))
        ],
    },
)


def _participant_taxon_ids(row):
    entity_types = f(
        'participant_entity_type',
        delimiter='||',
        map='entity_type',
        preserve_indices=True,
    ).extract(row)
    taxon_ids = f(
        'participant_ncbi_tax_id',
        delimiter='||',
        map='missing',
        preserve_indices=True,
    ).extract(row)
    size = max(len(entity_types), len(taxon_ids))
    values = []
    for idx in range(size):
        entity_type = entity_types[idx] if idx < len(entity_types) else None
        taxon_id = taxon_ids[idx] if idx < len(taxon_ids) else None
        values.append(
            taxon_id if entity_type in _TAXON_SCOPED_ENTITY_TYPES and taxon_id else None
        )
    return values


def _taxon_id_for_entity_type(entity_type, taxon_id):
    mapped_entity_type = (
        entity_type
        if isinstance(entity_type, EntityTypeCv)
        else entity_type_map.get(entity_type, EntityTypeCv.PHYSICAL_ENTITY)
    )
    return taxon_id if mapped_entity_type in _TAXON_SCOPED_ENTITY_TYPES and taxon_id else None


def _entity_taxon_id(entity_type_field, taxon_field):
    def extractor(row):
        entity_type = row.get(entity_type_field)
        taxon_id = f(taxon_field, map='missing').extract(row)
        return _taxon_id_for_entity_type(entity_type, taxon_id)

    return extractor


def _controller_member_taxon_ids(row):
    entity_types = f(
        'controller_member_entity_type',
        delimiter='||',
        map='entity_type',
        preserve_indices=True,
    ).extract(row)
    taxon_ids = f(
        'controller_member_ncbi_tax_id',
        delimiter='||',
        map='missing',
        preserve_indices=True,
    ).extract(row)
    size = max(len(entity_types), len(taxon_ids))
    values = []
    for idx in range(size):
        entity_type = entity_types[idx] if idx < len(entity_types) else None
        taxon_id = taxon_ids[idx] if idx < len(taxon_ids) else None
        values.append(_taxon_id_for_entity_type(entity_type, taxon_id))
    return values


def _split_values(value: object, *, delimiter: str = ';') -> list[str]:
    if not value or value == _MISSING_VALUE:
        return []
    return [
        item
        for item in (part.strip() for part in str(value).split(delimiter))
        if item and item != _MISSING_VALUE
    ]


def _pathway_ref(stable_id: str) -> EntityRef:
    return EntityRef(
        type=EntityTypeCv.PATHWAY,
        identifier_type=IdentifierNamespaceCv.REACTOME_STABLE_ID,
        identifier=stable_id,
    )


def _pathway_association(object_identifier) -> AssociationBuilder:
    return AssociationBuilder(
        object_entity_type=EntityTypeCv.PATHWAY,
        object_identifier_type=IdentifierNamespaceCv.REACTOME_STABLE_ID,
        object_identifier=object_identifier,
    )


def _pathway_associations(object_identifier) -> AssociationsBuilder:
    return AssociationsBuilder(_pathway_association(object_identifier))


def _cv_term_association(object_identifier) -> AssociationBuilder:
    return AssociationBuilder(
        object_entity_type=EntityTypeCv.CV_TERM,
        object_identifier_type=IdentifierNamespaceCv.CV_TERM_ACCESSION,
        object_identifier=object_identifier,
    )


def _cv_term_associations(object_identifier) -> AssociationsBuilder:
    return AssociationsBuilder(_cv_term_association(object_identifier))


def _combined_associations(
    *associations: AssociationBuilder,
) -> AssociationsBuilder:
    return AssociationsBuilder(*associations)


def _pathway_ontology_relations(row: dict) -> list[OntologyRelation]:
    relations: list[OntologyRelation] = []
    seen: set[tuple[str, str]] = set()
    for parent_group in _split_values(row.get('parent_pathway_reactome_stable_id'), delimiter='||'):
        for parent_id in _split_values(parent_group):
            key = ('part_of', parent_id)
            if key in seen:
                continue
            seen.add(key)
            relations.append(
                OntologyRelation(
                    predicate='part_of',
                    object=_pathway_ref(parent_id),
                    ontology_id=_REACTOME_PATHWAY_ONTOLOGY_ID,
                )
            )
    return relations


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
        CV(term=ReactionAnnotationsCv.EC_NUMBER, value=f('ec_number')),
        CV(term=InteractionMetadataCv.CONVERSION_DIRECTION, value=f('direction')),
    ),
    associations=_pathway_associations(
        f('pathway_term_accession', map='split')
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
                CV(term=IdentifierNamespaceCv.CHEBI, value=f('participant_chebi', delimiter='||', map='split_chebi')),
                CV(
                    term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                    value=f('participant_pubchem_compound', delimiter='||', map='split'),
                ),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('participant_kegg', delimiter='||', map='split')),
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
                    value=_participant_taxon_ids,
                ),
            ),
            entity_associations=_combined_associations(
                _pathway_association(
                    f('participant_pathway_term_accession', delimiter='||', map='split')
                ),
                _cv_term_association(
                    f('participant_go', delimiter='||', map='split')
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
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.CONTROL_TYPE, value=f('control_type')),
        CV(term=InteractionMetadataCv.CONTROL_EFFECT, value=f(_control_effect)),
        CV(term=f('causal_statement')),
    ),
    associations=_combined_associations(
        _pathway_association(
            f('pathway_term_accession', map='split')
        ),
        _cv_term_association(f('go')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=f('controller_entity_type', map='entity_type'),
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.NAME, value=f('controller_display_name', map='missing')),
                    CV(term=IdentifierNamespaceCv.SYNONYM, value=f('controller_synonyms', map='split')),
                    CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('controller_reactome_stable_id', map='split')),
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('controller_uniprot', map='split')),
                    CV(term=IdentifierNamespaceCv.CHEBI, value=f('controller_chebi', map='split_chebi')),
                    CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('controller_pubchem_compound', map='split')),
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('controller_kegg', map='split')),
                ),
                annotations=AnnotationsBuilder(
                    CV(
                        term=IdentifierNamespaceCv.NCBI_TAX_ID,
                        value=_entity_taxon_id('controller_entity_type', 'controller_ncbi_tax_id'),
                    ),
                ),
                associations=_combined_associations(
                    _pathway_association(
                        f('controller_pathway_term_accession', map='split')
                    ),
                    _cv_term_association(f('controller_go', map='split')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.CONTROLLER),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=f('controlled_entity_type', map='entity_type'),
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.NAME, value=f('controlled_display_name', map='missing')),
                    CV(term=IdentifierNamespaceCv.SYNONYM, value=f('controlled_synonyms', map='split')),
                    CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('controlled_reactome_stable_id', map='split')),
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('controlled_uniprot', map='split')),
                    CV(term=IdentifierNamespaceCv.CHEBI, value=f('controlled_chebi', map='split_chebi')),
                    CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('controlled_pubchem_compound', map='split')),
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('controlled_kegg', map='split')),
                ),
                annotations=AnnotationsBuilder(
                    CV(
                        term=IdentifierNamespaceCv.NCBI_TAX_ID,
                        value=_entity_taxon_id('controlled_entity_type', 'controlled_ncbi_tax_id'),
                    ),
                ),
                associations=_combined_associations(
                    _pathway_association(
                        f('controlled_pathway_term_accession', map='split')
                    ),
                    _cv_term_association(f('controlled_go', map='split')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.CONTROLLED),
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
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('ncbi_tax_id')),
        CV(term=PathwayAnnotationsCv.DESCRIPTION, value=f('definition')),
        CV(term=PathwayAnnotationsCv.COMMENT, value=f('comments')),
    ),
    associations=_cv_term_associations(f('go')),
    ontology_relations=_pathway_ontology_relations,
)


control_groups_schema = EntityBuilder(
    entity_type=f('controller_entity_type', map='entity_type'),
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('controller_display_name', map='missing')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('controller_synonyms', map='split')),
        CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('controller_reactome_stable_id', map='split')),
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('controller_uniprot', map='split')),
        CV(term=IdentifierNamespaceCv.CHEBI, value=f('controller_chebi', map='split_chebi')),
        CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('controller_pubchem_compound', map='split')),
        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('controller_kegg', map='split')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=IdentifierNamespaceCv.NCBI_TAX_ID,
            value=_entity_taxon_id('controller_entity_type', 'controller_ncbi_tax_id'),
        ),
    ),
    associations=_combined_associations(
        _pathway_association(
            f('controller_pathway_term_accession', map='split')
        ),
        _cv_term_association(f('controller_go', map='split')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=f('controller_member_entity_type', delimiter='||', map='entity_type'),
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.NAME, value=f('controller_member_display_name', delimiter='||', map='missing')),
                CV(term=IdentifierNamespaceCv.SYNONYM, value=f('controller_member_synonyms', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('controller_member_reactome_stable_id', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.UNIPROT, value=f('controller_member_uniprot', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.CHEBI, value=f('controller_member_chebi', delimiter='||', map='split_chebi')),
                CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('controller_member_pubchem_compound', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('controller_member_kegg', delimiter='||', map='split')),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term=IdentifierNamespaceCv.NCBI_TAX_ID,
                    value=_controller_member_taxon_ids,
                ),
            ),
            entity_associations=_combined_associations(
                _pathway_association(
                    f('controller_member_pathway_term_accession', delimiter='||', map='split')
                ),
                _cv_term_association(
                    f('controller_member_go', delimiter='||', map='split')
                ),
            ),
        ),
    ),
)

download = Download(
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
        raw_parser=partial(_raw, data_type='reactions'),
    ),
    pathways=Dataset(
        download=download,
        mapper=pathways_schema,
        raw_parser=partial(_raw, data_type='pathways'),
    ),
    controls=Dataset(
        download=download,
        mapper=controls_schema,
        raw_parser=partial(_raw, data_type='controls'),
    ),
    control_groups=Dataset(
        download=download,
        mapper=control_groups_schema,
        raw_parser=partial(_raw, data_type='control_groups'),
    ),
)
