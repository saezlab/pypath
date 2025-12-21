"""
Parse Reactome BioPAX data and emit Entity records.

This module converts Reactome BioPAX data into Entity records using the
declarative schema pattern.
"""

from __future__ import annotations

from pypath.internals.silver_schema import Entity, Identifier, Annotation, Membership
from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    MembershipRoleCv,
    MoleculeAnnotationsCv,
    ParticipantMetadataCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.inputs_v2.base import Dataset, Resource, ResourceConfig
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
# Helper Functions
# --------------------------------------------------------------------------- #

def _build_member_entities(participants: list[dict]) -> list[Entity]:
    """
    Build Entity objects for reaction participants.
    Handles both simple entities and families with members.
    """
    from pypath.internals.silver_schema import Membership

    role_mapping = {
        'reactant': BiologicalRoleCv.REACTANT,
        'product': BiologicalRoleCv.PRODUCT,
        'template': BiologicalRoleCv.TEMPLATE,
        'controller': BiologicalRoleCv.CONTROLLER,
        'controlled': BiologicalRoleCv.CONTROLLED,
        'pathway_component': BiologicalRoleCv.PATHWAY_COMPONENT,
        'member': MembershipRoleCv.MEMBER_OF,  # For family members
    }

    entities = []
    for p in participants:
        entity_type_str = p.get('entity_type', 'physical_entity')
        type_mapping = {v.value: v for v in EntityTypeCv}
        entity_type = type_mapping.get(entity_type_str, EntityTypeCv.PHYSICAL_ENTITY)

        identifiers = []
        annotations = []

        if p.get('display_name'):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=p['display_name']))

        if p.get('synonyms'):
            for syn in p['synonyms'].split(';'):
                if syn:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=syn))

        if p.get('reactome_stable_id'):
            for rid in p['reactome_stable_id'].split(';'):
                if rid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))

        if p.get('uniprot'):
            for uid in p['uniprot'].split(';'):
                if uid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.UNIPROT, value=uid))

        if p.get('chebi'):
            for cid in p['chebi'].split(';'):
                if cid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.CHEBI, value=cid))

        if p.get('pubchem_compound'):
            for pcid in p['pubchem_compound'].split(';'):
                if pcid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=pcid))

        if p.get('kegg'):
            for kid in p['kegg'].split(';'):
                if kid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.KEGG_COMPOUND, value=kid))

        if p.get('go'):
            for goid in p['go'].split(';'):
                if goid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=goid))

        role_str = p.get('role', '')
        role_cv = role_mapping.get(role_str)
        if role_cv:
            annotations.append(Annotation(term=role_cv))

        if p.get('stoichiometry'):
            annotations.append(Annotation(term=ParticipantMetadataCv.STOICHIOMETRY, value=p['stoichiometry']))

        if p.get('ncbi_tax_id'):
            annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=p['ncbi_tax_id']))

        # Check if this is a family with members
        if p.get('is_family') and p.get('members'):
            from pypath.internals.silver_schema import Membership

            # Build member entities first
            family_membership = []
            for member in p['members']:
                member_identifiers = []
                member_annotations = []

                if member.get('display_name'):
                    member_identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=member['display_name']))

                if member.get('synonyms'):
                    for syn in member['synonyms'].split(';'):
                        if syn:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=syn))

                if member.get('reactome_stable_id'):
                    for rid in member['reactome_stable_id'].split(';'):
                        if rid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))

                if member.get('uniprot'):
                    for uid in member['uniprot'].split(';'):
                        if uid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.UNIPROT, value=uid))

                if member.get('chebi'):
                    for cid in member['chebi'].split(';'):
                        if cid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.CHEBI, value=cid))

                if member.get('pubchem_compound'):
                    for pcid in member['pubchem_compound'].split(';'):
                        if pcid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=pcid))

                if member.get('kegg'):
                    for kid in member['kegg'].split(';'):
                        if kid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.KEGG_COMPOUND, value=kid))

                if member.get('go'):
                    for goid in member['go'].split(';'):
                        if goid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=goid))

                # Add member role
                member_annotations.append(Annotation(term=MembershipRoleCv.MEMBER_OF))

                if member.get('ncbi_tax_id'):
                    member_annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=member['ncbi_tax_id']))

                # Get member entity type
                member_type_str = member.get('entity_type', 'physical_entity')
                member_type = type_mapping.get(member_type_str, EntityTypeCv.PHYSICAL_ENTITY)

                # Create member entity
                member_entity = Entity(
                    type=member_type,
                    identifiers=member_identifiers if member_identifiers else None,
                    annotations=member_annotations if member_annotations else None
                )

                family_membership.append(Membership(
                    member=member_entity,
                    is_parent=False,
                    annotations=member_annotations if member_annotations else None
                ))

            # Create family entity with members
            family_entity = Entity(
                type=entity_type,
                identifiers=identifiers if identifiers else None,
                annotations=annotations if annotations else None,
                membership=family_membership if family_membership else None
            )

            entities.append((family_entity, annotations if annotations else None))
        else:
            # Simple entity without members
            entities.append((
                Entity(type=entity_type, identifiers=identifiers if identifiers else None, annotations=annotations if annotations else None),
                annotations if annotations else None,
            ))

    return entities


# --------------------------------------------------------------------------- #
# Mapper Functions
# --------------------------------------------------------------------------- #

def _map_reaction_record(record: dict) -> Entity:
    from pypath.internals.silver_schema import Membership

    identifiers = []
    if record.get('display_name'):
        identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=record['display_name']))
    if record.get('synonyms'):
        for syn in record['synonyms'].split(';'):
            if syn:
                identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=syn))
    if record.get('reactome_stable_id'):
        for rid in record['reactome_stable_id'].split(';'):
            if rid:
                identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))
    if record.get('reactome_id'):
        for rid in record['reactome_id'].split(';'):
            if rid:
                identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_ID, value=rid))

    annotations = []
    if record.get('pubmed'):
        for pmid in record['pubmed'].split(';'):
            if pmid:
                annotations.append(Annotation(term=IdentifierNamespaceCv.PUBMED, value=pmid))
    if record.get('ec_number'):
        annotations.append(Annotation(term=MoleculeAnnotationsCv.EC_NUMBER, value=record['ec_number']))
    if record.get('direction'):
        annotations.append(Annotation(term=InteractionMetadataCv.CONVERSION_DIRECTION, value=record['direction']))

    membership = []
    all_participants = record.get('reactants', []) + record.get('products', []) + record.get('templates', [])

    for entity, member_annotations in _build_member_entities(all_participants):
        membership.append(Membership(
            member=entity,
            is_parent=False,
            annotations=member_annotations,
        ))

    type_mapping = {v.value: v for v in EntityTypeCv}
    entity_type = type_mapping.get(record.get('entity_type', ''), EntityTypeCv.REACTION)

    return Entity(
        type=entity_type,
        identifiers=identifiers if identifiers else None,
        annotations=annotations if annotations else None,
        membership=membership if membership else None,
    )


def _map_pathway_record(record: dict) -> Entity:
    from pypath.internals.silver_schema import Membership

    identifiers = []
    if record.get('display_name'):
        identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=record['display_name']))
    if record.get('synonyms'):
        for syn in record['synonyms'].split(';'):
            if syn:
                identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=syn))
    if record.get('reactome_stable_id'):
        for rid in record['reactome_stable_id'].split(';'):
            if rid:
                identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))
    if record.get('reactome_id'):
        for rid in record['reactome_id'].split(';'):
            if rid:
                identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_ID, value=rid))
    if record.get('go'):
        for goid in record['go'].split(';'):
            if goid:
                identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=goid))

    annotations = []
    if record.get('pubmed'):
        for pmid in record['pubmed'].split(';'):
            if pmid:
                annotations.append(Annotation(term=IdentifierNamespaceCv.PUBMED, value=pmid))
    if record.get('ncbi_tax_id'):
        annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=record['ncbi_tax_id']))
    if record.get('description'):
        annotations.append(Annotation(term=MoleculeAnnotationsCv.DESCRIPTION, value=record['description']))

    membership = []
    for comp in record.get('components', []):
        type_str = comp.get('type', '')
        if 'Pathway' in type_str:
            comp_type = EntityTypeCv.PATHWAY
        elif 'BiochemicalReaction' in type_str:
            comp_type = EntityTypeCv.REACTION
        elif 'Catalysis' in type_str:
            comp_type = EntityTypeCv.CATALYSIS
        elif 'Control' in type_str:
            comp_type = EntityTypeCv.CONTROL
        else:
            comp_type = EntityTypeCv.INTERACTION

        comp_identifiers = []
        if comp.get('display_name'):
            comp_identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=comp['display_name']))
        if comp.get('reactome_stable_id'):
            for rid in comp['reactome_stable_id'].split(';'):
                if rid:
                    comp_identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))

        comp_annotations = [Annotation(term=BiologicalRoleCv.PATHWAY_COMPONENT)]
        if comp.get('step_order') is not None:
            comp_annotations.append(Annotation(term=ParticipantMetadataCv.STEP_ORDER, value=str(comp['step_order'])))

        membership.append(Membership(
            member=Entity(type=comp_type, identifiers=comp_identifiers if comp_identifiers else None),
            is_parent=False,
            annotations=comp_annotations,
        ))

    return Entity(
        type=EntityTypeCv.PATHWAY,
        identifiers=identifiers if identifiers else None,
        annotations=annotations if annotations else None,
        membership=membership if membership else None,
    )


def _map_control_record(record: dict) -> Entity:
    from pypath.internals.silver_schema import Membership

    identifiers = []
    if record.get('display_name'):
        identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=record['display_name']))
    elif record.get('controller', {}).get('display_name'):
        identifiers.append(Identifier(
            type=IdentifierNamespaceCv.NAME,
            value=f"{record['control_class']} by {record['controller']['display_name']}",
        ))

    if record.get('reactome_stable_id'):
        for rid in record['reactome_stable_id'].split(';'):
            if rid:
                identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))
    if record.get('reactome_id'):
        for rid in record['reactome_id'].split(';'):
            if rid:
                identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_ID, value=rid))
    if record.get('go'):
        for goid in record['go'].split(';'):
            if goid:
                identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=goid))

    annotations = []
    if record.get('control_type'):
        annotations.append(Annotation(term=InteractionMetadataCv.CONTROL_TYPE, value=record['control_type']))

    membership = []

    controller = record.get('controller', {})
    if controller:
        controller_identifiers = []
        controller_annotations = [Annotation(term=BiologicalRoleCv.CONTROLLER)]

        if controller.get('display_name'):
            controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=controller['display_name']))

        if controller.get('synonyms'):
            for syn in controller['synonyms'].split(';'):
                if syn:
                    controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=syn))

        if controller.get('reactome_stable_id'):
            for rid in controller['reactome_stable_id'].split(';'):
                if rid:
                    controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))

        if controller.get('uniprot'):
            for uid in controller['uniprot'].split(';'):
                if uid:
                    controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.UNIPROT, value=uid))

        if controller.get('pubchem_compound'):
            for pcid in controller['pubchem_compound'].split(';'):
                if pcid:
                    controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=pcid))

        if controller.get('kegg'):
            for kid in controller['kegg'].split(';'):
                if kid:
                    controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.KEGG_COMPOUND, value=kid))

        if controller.get('go'):
            for goid in controller['go'].split(';'):
                if goid:
                    controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=goid))

        if controller.get('ncbi_tax_id'):
            controller_annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=controller['ncbi_tax_id']))

        type_mapping = {v.value: v for v in EntityTypeCv}
        controller_type = type_mapping.get(controller.get('entity_type', ''), EntityTypeCv.PHYSICAL_ENTITY)

        if controller_identifiers:
            membership.append(Membership(
                member=Entity(type=controller_type, identifiers=controller_identifiers, annotations=controller_annotations if controller_annotations else None),
                is_parent=False,
                annotations=[Annotation(term=BiologicalRoleCv.CONTROLLER)],
            ))

    controlled = record.get('controlled', {})
    if controlled:
        controlled_identifiers = []
        if controlled.get('display_name'):
            controlled_identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=controlled['display_name']))
        if controlled.get('reactome_stable_id'):
            for rid in controlled['reactome_stable_id'].split(';'):
                if rid:
                    controlled_identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))

        controlled_type_str = controlled.get('type', '')
        if 'BiochemicalReaction' in controlled_type_str:
            controlled_type = EntityTypeCv.REACTION
        elif 'Pathway' in controlled_type_str:
            controlled_type = EntityTypeCv.PATHWAY
        else:
            controlled_type = EntityTypeCv.INTERACTION

        if controlled_identifiers:
            membership.append(Membership(
                member=Entity(type=controlled_type, identifiers=controlled_identifiers),
                is_parent=False,
                annotations=[Annotation(term=BiologicalRoleCv.CONTROLLED)],
            ))

    entity_type = EntityTypeCv.INTERACTION

    return Entity(
        type=entity_type,
        identifiers=identifiers if identifiers else None,
        annotations=annotations if annotations else None,
        membership=membership if membership else None,
    )


resource = Resource(
    config,
    reactions=Dataset(
        download=None,
        mapper=_map_reaction_record,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='reactions', **kwargs),
    ),
    pathways=Dataset(
        download=None,
        mapper=_map_pathway_record,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='pathways', **kwargs),
    ),
    controls=Dataset(
        download=None,
        mapper=_map_control_record,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='controls', **kwargs),
    ),
)

reactions = resource.reactions
pathways = resource.pathways
controls = resource.controls
