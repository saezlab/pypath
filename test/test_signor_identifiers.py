"""SIGNOR source crosswalk and native builder contract tests (no downloads)."""
import pytest
from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core import Namespace
from omnipath_core.biolink import entity_type
from pypath.inputs_v2.signor import (
    complexes_schema, phenotypes_schema, stimuli_schema, protein_families_schema,
    interactions_schema, _infer_signor_interactor_type,
)
from pypath.internals.tabular_builder import CV, EntityBuilder, IdentifiersBuilder, AnnotationsBuilder, MembersFromList, MembershipBuilder


def interaction_row(causal='MI:2236'):
    return {
        '\ufeff#ID(s) interactor A': 'uniprotkb:P04637',
        'ID(s) interactor B': 'uniprotkb:P38398',
        'Type(s) interactor A': 'psi-mi:"MI:0326"(protein)',
        'Type(s) interactor B': 'psi-mi:"MI:0326"(protein)',
        'Taxid interactor A': 'taxid:9606(human)',
        'Taxid interactor B': 'taxid:9606(human)',
        'Causal statement': f'psi-mi:"{causal}"(source label)',
        'Interaction identifier(s)': 'SIGNOR:SIGNOR-1',
        'Publication Identifier(s)': 'pubmed:123456',
    }


@pytest.mark.parametrize('builder,cls,id,name_field', [
    (complexes_schema, model.MacromolecularComplex, 'SIGNOR-C1', 'COMPLEX NAME'),
    (protein_families_schema, model.ProteinFamily, 'SIGNOR-PF1', 'PROT. FAMILY NAME'),
    (phenotypes_schema, model.PhenotypicFeature, 'SIGNOR-PH1', 'PHENOTYPE NAME'),
    (stimuli_schema, model.ExposureEvent, 'SIGNOR-ST1', 'STIMULUS NAME'),
])
def test_entity_identifiers_and_attributes(builder, cls, id, name_field):
    entity = builder.build({'SIGNOR ID': id, name_field: 'example', 'LIST OF ENTITIES': 'P05412,P01100'})
    assert entity.type == entity_type(cls)
    assert entity.identifiers[0].type == Namespace.SIGNOR
    assert entity.identifiers[0].value == id
    assert all(i.type != 'name' for i in entity.identifiers)
    assert any(a.term == 'name' and a.value == 'example' for a in entity.annotations)


@pytest.mark.parametrize('code,direction,aspect', [
    ('MI:2235', 'increased', None), ('MI:2236', 'increased', 'activity'),
    ('MI:2237', 'increased', 'abundance'), ('MI:2238', 'increased', 'expression'),
    ('MI:2239', 'increased', 'stability'), ('MI:2240', 'decreased', None),
    ('MI:2241', 'decreased', 'activity'), ('MI:2242', 'decreased', 'abundance'),
    ('MI:2243', 'decreased', 'stability'), ('MI:2244', 'decreased', 'expression'),
])
def test_causal_qualifiers_are_source_accession_driven(code, direction, aspect):
    relation = interactions_schema.build(interaction_row(code))
    assert relation.predicate == 'affects'
    annotations = {a.term: a.value for a in relation.annotations}
    assert 'original_predicate' not in annotations
    assert annotations['object_direction_qualifier'] == direction
    assert annotations.get('object_aspect_qualifier') == aspect
    assert annotations['publications'] == 'PMID:123456'
    assert relation.subject.annotations[0].value == 'NCBITaxon:9606'
    assert relation.subject.identifiers[0].type == Namespace.UNIPROT


def test_unknown_source_semantics_fail_instead_of_falling_back():
    with pytest.raises(ValueError, match='causal statement'):
        interactions_schema.build(interaction_row('MI:9999'))
    row = interaction_row()
    row['Type(s) interactor A'] = 'psi-mi:"MI:9999"(unknown)'
    with pytest.raises(ValueError, match='interactor type'):
        interactions_schema.build(row)


@pytest.mark.parametrize('identifier,cls', [
    ('uniprotkb:URS0000123', model.RNAProduct),
    ('uniprotkb:DB00091', model.ChemicalEntity),
    ('signor:SIGNOR-FP1', model.ProteinFamily),
])
def test_documented_missing_type_repairs(identifier, cls):
    assert _infer_signor_interactor_type({'ID(s) interactor A': identifier}, 'A') is cls
    with pytest.raises(ValueError):
        _infer_signor_interactor_type({'ID(s) interactor A': 'unknown:XYZ'}, 'A')


def test_builder_accepts_native_enum_values_in_nested_members():
    builder = EntityBuilder(
        entity_type=model.MacromolecularComplex,
        identifiers=IdentifiersBuilder(CV(term=Namespace.SIGNOR, value='SIGNOR-C1')),
        membership=MembershipBuilder(MembersFromList(
            entity_type=model.Protein,
            identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value='P04637')),
            annotations=AnnotationsBuilder(CV(term=slots.object_direction_qualifier,
                                             value=model.DirectionQualifierEnum.increased)),
        )),
    )
    member = builder.build({}).membership[0]
    assert member.member.type == 'protein'
    assert member.annotations[0].value == 'increased'


def test_legacy_types_are_not_adapted():
    from pypath.internals.cv_terms import EntityTypeCv
    for invalid in [EntityTypeCv.PROTEIN, 'nonsense', 42]:
        with pytest.raises(ValueError):
            EntityBuilder(entity_type=invalid)
