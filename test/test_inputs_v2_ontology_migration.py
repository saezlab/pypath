"""Source semantics survive ontology/RNA mapping and the Silver boundary."""
import io
from types import SimpleNamespace

import pytest
from biolink_model.datamodel.model import ChemicalEntity, MolecularActivity, slots
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import Entity, Identifier, Membership
from pypath.inputs_v2 import chemont, go, hpo, mirbase, mondo, omnipath_ontology, psi_mi
from pypath.internals.tabular_builder import EntityBuilder, IdentifiersBuilder, CV, Member, MembersFromList, MembershipBuilder, FieldConfig


def extract(record):
    result = SilverExtractor('fixture', 'terms')
    result.process_record(record, raw_payload={}, row_id='1', row_number=1)
    return result


@pytest.mark.parametrize('module,prefix,namespace', [
    (chemont, 'CHEMONTID', 'chemont'), (go, 'GO', 'go'),
    (hpo, 'HP', 'hpo'), (mondo, 'MONDO', 'mondo'), (psi_mi, 'MI', 'mi'),
])
def test_ontology_hierarchy_and_unknown_relations(module, prefix, namespace):
    row = {'id': f'{prefix}:1', 'name': 'a concept', 'is_a': [f'{prefix}:2'],
           'synonyms': ['another name'], 'alt_ids': [f'{prefix}:3'],
           'relationships': [{'type': 'unknown_relation', 'target': 'OTHER:4'}]}
    result = extract(module.terms_schema(row))
    assert {e.entity_type for e in result.entities.values()} == {'ontology_class'}
    assert {r.predicate for r in result.relations} == {'subclass_of'}
    assert {e.namespace for e in result.entities.values()} == {namespace}
    source = next(e for e in result.entities.values() if e.identifier == f'{prefix}:1')
    assert ('name', 'a concept') in {(i['ns'], i['id']) for i in source.identifiers}
    assert ('synonym', 'another name') in {(i['ns'], i['id']) for i in source.identifiers}
    assert source.label == 'a concept'
    assert not any(a['term'] in {'name', 'synonym'} for a in source.annotations)
    assert not any(a['term'] == 'description'
                   for e in result.entities.values() for a in e.annotations)
    assert module.terms_schema({**row, 'is_obsolete': True}) is None
    assert module.terms_schema({}) is None


def test_hpo_uses_gene_identity_and_explicit_ontology_namespace():
    result = extract(hpo.annotations_schema({'ncbi_gene_id': '1', 'gene_symbol': 'A1BG', 'hpo_id': 'HP:1'}))
    assert {e.entity_type for e in result.entities.values()} == {'gene', 'ontology_class'}
    assert {r.predicate for r in result.relations} == {'associated_with'}
    assert next(e for e in result.entities.values() if e.entity_type == 'gene').taxon == '9606'


def test_mondo_obsolete_and_unknown_relationships_do_not_emit():
    text = '[Term]\nid: MONDO:1\nrelationship: RO:0004001 HGNC:1\nis_obsolete: true\n\n[Term]\nid: MONDO:2\nrelationship: RO:0004001 HGNC:2\n'
    rows = list(mondo.iter_gene_disease_associations(SimpleNamespace(result=io.StringIO(text))))
    assert [r['mondo_id'] for r in rows] == ['MONDO:2']
    result = extract(mondo.gene_disease_association_to_entity(rows[0]))
    assert {r.predicate for r in result.relations} == {'associated_with'}
    assert mondo.gene_disease_association_to_entity({**rows[0], 'relation': 'unknown'}) is None


def test_mirna_maturation_has_correct_direction_and_distinct_identifiers():
    result = extract(mirbase.matures_schema({'mirbase_mat': 'MIMAT1', 'name': 'mature', 'precursors': ['MI1']}))
    relation, = result.relations
    assert relation.predicate == 'derives_from'
    assert result.entities[relation.subject_entity_key].namespace == 'mirbase_mature'
    assert result.entities[relation.object_entity_key].namespace == 'mirbase_precursor'
    assert not extract(mirbase.precursors_schema({'mirbase_pre': 'MI1'})).relations


def test_legacy_vocabulary_does_not_invent_psi_mi_parent():
    assert all(row['is_a'] != 'MI:0000' for row in omnipath_ontology._extract_om_terms())
    result = extract(omnipath_ontology.terms_schema({'accession': 'OM:1', 'name': 'archived term'}))
    assert not result.relations


@pytest.mark.parametrize('predicate', [slots.has_input, slots.has_output, slots.enabled_by])
def test_nested_membership_accepts_explicit_native_predicate(predicate):
    child = Entity(type=ChemicalEntity, identifiers=[Identifier(Namespace.CHEBI, '1')])
    parent = Entity(type=MolecularActivity, identifiers=[Identifier(Namespace.RHEA, '1')],
                    membership=[Membership(child, predicate=predicate)])
    relation, = extract(parent).relations
    assert relation.predicate == predicate.name.replace(' ', '_')


def test_member_builders_preserve_predicate_and_reject_unknown():
    f = FieldConfig()
    chemical = EntityBuilder(entity_type=ChemicalEntity, identifiers=IdentifiersBuilder(CV(term=Namespace.CHEBI, value=f('id'))))
    for member in [Member(entity=chemical, predicate=slots.has_input), MembersFromList(entity_type=ChemicalEntity, identifiers=chemical.identifiers, predicate=slots.has_input)]:
        builder = EntityBuilder(entity_type=MolecularActivity, identifiers=IdentifiersBuilder(CV(term=Namespace.RHEA, value='1')), membership=MembershipBuilder(member))
        assert {r.predicate for r in extract(builder({'id': '2'})).relations} == {'has_input'}
    with pytest.raises(ValueError):
        Member(entity=chemical, predicate='invented')


def test_dynamic_member_predicates_keep_alignment_and_skip_unspecified():
    f = FieldConfig()
    builder = EntityBuilder(entity_type=MolecularActivity,
        identifiers=IdentifiersBuilder(CV(term=Namespace.RHEA, value='1')),
        membership=MembershipBuilder(MembersFromList(
            entity_type=ChemicalEntity,
            identifiers=IdentifiersBuilder(CV(term=Namespace.CHEBI, value=f('ids'))),
            predicate=f('roles', preserve_indices=True, transform=lambda value: {'input': slots.has_input, 'output': slots.has_output}.get(value)),
        )))
    result = extract(builder({'ids': ['1', '2', '3'], 'roles': ['input', None, 'output']}))
    assert {r.predicate for r in result.relations} == {'has_input', 'has_output'}
    assert len(result.entities) == 3
    output = next(r for r in result.relations if r.predicate == 'has_output')
    assert result.entities[output.object_entity_key].identifier == '3'


def test_mirbase_raw_mature_rows_link_back_to_precursor(monkeypatch):
    monkeypatch.setattr(mirbase, '_precursor_to_matures', lambda: {'MI1': ['MIMAT1']})
    monkeypatch.setattr(mirbase, 'mirbase_mirna_mature', lambda _: [(None, 'mature', 'precursor', 'MIMAT1')])
    row, = mirbase._matures_raw()
    assert row['precursors'] == ['MI1']
    assert {r.predicate for r in extract(mirbase.matures_schema(row)).relations} == {'derives_from'}
