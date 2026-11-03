"""Execution contracts for reusable row mappings."""

import pytest
from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from pypath.internals.tabular_builder import (
    AnnotationsBuilder, AssociationBuilder, Column, ColumnCache, CV, EntityBuilder, IdentifiersBuilder,
)


def test_pair_columns_are_aligned_and_evaluated_once_per_row():
    calls = []

    def pairs(row):
        calls.append(row)
        return row['pairs']

    builder = IdentifiersBuilder(CV.from_pairs(pairs))
    first = {'pairs': [(Namespace.UNIPROT, 'P1'), (Namespace.ENTREZ, None),
                       (Namespace.CHEBI, '42')]}
    second = {'pairs': [(Namespace.ENTREZ, '123')]}
    assert [(x.type, x.value) for x in builder.build(first)] == [
        (Namespace.UNIPROT, 'P1'), (Namespace.CHEBI, '42'),
    ]
    assert [(x.type, x.value) for x in builder.build(second)] == [(Namespace.ENTREZ, '123')]
    assert calls == [first, second]


def test_empty_shared_cache_reuses_an_extraction_across_builders():
    calls = []
    column = Column(lambda row: calls.append(row) or row['value'])
    identifiers = IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=column))
    annotations = AnnotationsBuilder(CV(term=slots.name, value=column))
    row = {'value': 'P1'}
    cache = ColumnCache()
    assert identifiers.build(row, cache)[0].value == 'P1'
    assert annotations.build(row, cache)[0].value == 'P1'
    assert calls == [row]


def test_static_and_dynamic_terms_still_deduplicate_identically():
    builder = AnnotationsBuilder(
        CV(term=slots.name, value=Column('name')),
        CV(term=lambda row: slots.name, value=Column('name')),
    )
    for name in ['one', 'two', 'one']:
        result = builder.build({'name': name})
        assert len(result) == 1
        assert result[0].value == name


@pytest.mark.parametrize('value', [None, '', 'name', 0, False,
                                  ['a', None, '', 'a', 'b'],
                                  [['a', 'b'], ['c']], ('a', 'b')])
def test_static_fast_path_matches_dynamic_broadcasting(value):
    static = AnnotationsBuilder(CV(term=slots.name, value=lambda row: row['value']))
    dynamic = AnnotationsBuilder(CV(term=lambda row: slots.name, value=lambda row: row['value']))
    row = {'value': value}
    assert static.build(row) == dynamic.build(row)


def test_presence_annotations_match_dynamic_terms():
    static = AnnotationsBuilder(CV(term=slots.name))
    dynamic = AnnotationsBuilder(CV(term=lambda row: slots.name))
    assert static.build({}) == dynamic.build({})




















def test_entity_type_and_identifier_share_row_local_extraction():
    calls = []
    source = Column(lambda row: calls.append(1) or row['kind'])
    builder = EntityBuilder(
        entity_type=source,
        identifiers=IdentifiersBuilder(CV(term=Namespace.SIGNOR, value=source)),
    )
    assert builder({'kind': 'biolink:Protein'}) is not None
    assert calls == [1]




def test_relation_emission_reuses_the_same_mapping_and_validation():
    from pypath.internals.tabular_builder import RelationBuilder

    endpoint = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=Column('id'))),
    )
    relation = RelationBuilder(subject=endpoint, predicate=slots.affects, object=endpoint)
    row = {'id': 'P1'}
    assert relation.build(row, emit=dict) == relation.build(row)._asdict()
    invalid = RelationBuilder(subject=endpoint, predicate=lambda row: 'invalid-predicate', object=endpoint)
    with pytest.raises(ValueError):
        invalid.build(row, emit=dict)


def test_associations_deduplicate_identifiers_under_a_slot_predicate():
    from biolink_model.datamodel.model import OntologyClass

    builder = AssociationBuilder(
        predicate=slots.associated_with,
        object_entity_type=OntologyClass,
        object_identifier_type=Namespace.GO,
        object_identifier=Column('go', delimiter=';'),
    )
    associations = builder.build({'go': 'GO:0001;GO:0002;GO:0001'}, ColumnCache())
    assert [a.object.identifier for a in associations] == ['GO:0001', 'GO:0002']
    assert {(a.predicate, a.object.type, a.object.identifier_type) for a in associations} == {
        ('associated_with', 'ontology_class', 'go'),
    }
    assert builder.build({'go': ''}, ColumnCache()) == []


def test_measurements_deduplicate_by_their_fields():
    from pypath.inputs_v2._measurements import measurement

    values = [measurement('5', 'nM', 'Ki (nM)'), measurement('5', 'nM', 'Ki (nM)'),
              measurement('7', 'nM', 'Ki (nM)'), measurement('5', 'nM', 'Kd (nM)')]
    builder = AnnotationsBuilder(CV(term=slots.has_quantitative_value, value=lambda row: values))
    assert len(builder.build({})) == 3
