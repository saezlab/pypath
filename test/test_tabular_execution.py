"""Execution contracts for reusable row mappings."""

import pytest
from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from pypath.internals.tabular_builder import (
    AnnotationsBuilder, Column, ColumnCache, CV, EntityBuilder, IdentifiersBuilder,
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


def test_flat_entity_cache_tracks_fields_and_isolates_mutation():
    calls = []

    def identifier(row):
        calls.append(row['id'])
        return row['id']

    builder = EntityBuilder(
        entity_type=Protein, cache_by=('id', 'name'), cache_size=2,
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=identifier)),
        annotations=AnnotationsBuilder(CV(term=slots.name, value=Column('name'))),
    )
    first = builder({'id': 'P1', 'name': 'first'})
    first.identifiers.clear()
    first.annotations.clear()
    repeated = builder({'id': 'P1', 'name': 'first'})
    assert repeated.identifiers[0].value == 'P1'
    assert repeated.annotations[0].value == 'first'
    assert calls == ['P1']
    assert builder({'id': 'P1', 'name': 'changed'}).annotations[0].value == 'changed'
    builder({'id': 'P2', 'name': 'second'})
    builder({'id': 'P1', 'name': 'first'})
    assert calls == ['P1', 'P1', 'P2', 'P1']  # bounded cache evicted first


def test_flat_cache_distinguishes_missing_null_and_scalar_types():
    def identifier(row):
        return 'missing' if 'id' not in row else f'{type(row["id"]).__name__}:{row["id"]}'

    builder = EntityBuilder(
        entity_type=Protein, cache_by=('id',),
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=identifier)),
    )
    assert [builder(row).identifiers[0].value for row in [{}, {'id': None}, {'id': 0}, {'id': False}]] == [
        'missing', 'NoneType:None', 'int:0', 'bool:False',
    ]


def test_column_memoization_is_automatic_bounded_and_isolated():
    calls = []
    column = Column('id', transform=lambda value: calls.append(value) or value.upper(), cache_size=2)
    first = column.extract({'id': 'a'})
    first.clear()
    assert column.extract({'id': 'a'}) == ['A']
    column.extract({'id': 'b'})
    column.extract({'id': 'a'})  # refresh LRU
    column.extract({'id': 'c'})
    column.extract({'id': 'b'})
    assert calls == ['a', 'b', 'c', 'b']


def test_column_bypasses_mutable_inputs_and_outputs():
    calls = []
    column = Column('id', transform=lambda value: calls.append(value) or value)
    row = {'id': ['a']}
    assert column.extract(row) == ['a']
    row['id'].append('b')
    assert column.extract(row) == ['a', 'b']
    assert calls == ['a', 'a', 'b']
    mutable = Column('id', transform=lambda value: [value])
    first = mutable.extract({'id': 'a'})
    first[0].append('changed')
    assert mutable.extract({'id': 'a'}) == [['a']]


def test_pair_column_reuses_cell_transformation_across_different_rows():
    calls = []
    pairs = Column('ids', transform=lambda value: calls.append(value) or ((Namespace.UNIPROT, value),))
    builder = IdentifiersBuilder(CV.from_pairs(pairs))
    assert builder.build({'ids': 'P1', 'other': 'a'}) == builder.build({'ids': 'P1', 'other': 'b'})
    assert calls == ['P1']


def test_entity_infers_fields_and_checks_dynamic_type_each_time(monkeypatch):
    from biolink_model.datamodel.model import ChemicalEntity

    calls = []
    builder = EntityBuilder(
        entity_type=lambda row: Protein if row['kind'] == 'protein' else ChemicalEntity,
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=Column('id'))),
        annotations=AnnotationsBuilder(CV(term=slots.name, value=Column('name'))),
    )
    original = builder._build
    monkeypatch.setattr(builder, '_build', lambda *args: calls.append(1) or original(*args))
    first = builder({'id': 'P1', 'name': 'one', 'kind': 'protein', 'evidence': 'a'})
    first.annotations.clear()
    repeated = builder({'id': 'P1', 'name': 'one', 'kind': 'protein', 'evidence': 'b'})
    assert repeated.annotations[0].value == 'one'
    assert len(calls) == 1
    changed_type = builder({'id': 'P1', 'name': 'one', 'kind': 'chemical'})
    assert changed_type.type != repeated.type
    assert len(calls) == 2
    assert builder({'id': 'P1', 'name': 'two', 'kind': 'protein'}).annotations[0].value == 'two'
    assert len(calls) == 3
    with pytest.raises(KeyError):
        builder({'id': 'P1', 'name': 'one'})


def test_opaque_callback_cache_uses_whole_row_including_presence_and_order():
    calls = []

    def identifier(row):
        calls.append(1)
        return repr(list(row.items()))

    builder = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=identifier)),
    )
    rows = [{}, {'id': None}, {'id': False}, {'id': 0}, {'id': -0.0}, {'id': 0.0},
            {'a': 1, 'b': 2}, {'b': 2, 'a': 1}]
    results = [builder(row) for row in rows]
    assert len({result.identifiers[0].value for result in results}) == len(rows)
    assert [builder(row) for row in rows] == results
    assert len(calls) == len(rows)


def test_disabling_caches_reexecutes_transformations():
    calls = []
    column = Column('id', transform=lambda value: calls.append(value) or value, cache_size=0)
    builder = EntityBuilder(
        entity_type=Protein, cache_size=0,
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=column)),
    )
    assert builder({'id': 'a'}) == builder({'id': 'a'})
    assert calls == ['a', 'a']


def test_mutable_opaque_row_is_not_cached():
    builder = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=lambda row: row['ids'][0])),
    )
    row = {'ids': ['P1']}
    assert builder(row).identifiers[0].value == 'P1'
    row['ids'][0] = 'P2'
    assert builder(row).identifiers[0].value == 'P2'


def test_entity_type_and_identifier_share_row_local_extraction():
    calls = []
    source = Column(lambda row: calls.append(1) or row['kind'])
    builder = EntityBuilder(
        entity_type=source,
        identifiers=IdentifiersBuilder(CV(term=Namespace.SIGNOR, value=source)),
    )
    assert builder({'kind': 'biolink:Protein'}) is not None
    assert calls == [1]


def test_nested_entities_do_not_reuse_complete_results():
    calls = []
    builder = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=Column('id'))),
        ontology_relations=lambda row: calls.append(row) or [],
    )
    assert builder({'id': 'P1'}) == builder({'id': 'P1'})
    assert len(calls) == 2


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
