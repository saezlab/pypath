"""Source contract of the NRF2ome input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor

from pypath.inputs_v2.nrf2ome import resource

# Lines of nrf2ome.csv: an undirected link, an inhibition, a directed link
# without effect and the link of a protein to itself (undirected).
LINES = (
    'O00221,O15379,1,0,0,18241676,\n'
    'O00221,Q16236,0,1,-1,18241676,\n'
    'O15525,Q16236,1,1,0,18585411,\n'
    'Q9ULX9,Q9ULX9,1,0,0,12490281,\n'
)


def rows(text=LINES):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(resource.interactions._raw_parser(opener))


def extract(row):
    record = resource.interactions.mapper(row)
    assert record is not None
    out = SilverExtractor('nrf2ome', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    return out.relations[0]


def terms(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_reader_names_the_columns():
    parsed = rows()
    assert parsed[1] == {
        'source': 'O00221',
        'target': 'Q16236',
        'is_direct': '0',
        'is_directed': '1',
        'effect': '-1',
        'references': '18241676',
        'trailing': '',
    }


def test_reader_rejects_a_line_with_another_number_of_fields_and_an_empty_file():
    with pytest.raises(ValueError, match='fields'):
        rows('O00221,O15379,1,0,0,18241676,9,\n')
    with pytest.raises(ValueError, match='no rows'):
        rows('')


def test_an_inhibition_is_a_decrease_with_its_reference():
    relation = extract(rows()[1])
    assert relation.predicate == 'affects'
    assert terms(relation, 'object_direction_qualifier') == ['decreased']
    assert terms(relation, 'publications') == ['PMID:18241676']


def test_a_directed_link_without_effect_has_no_direction_qualifier():
    relation = extract(rows()[2])
    assert relation.predicate == 'affects'
    assert terms(relation, 'object_direction_qualifier') == []
    assert terms(relation, 'object_aspect_qualifier') == []


def test_an_undirected_link_is_an_interaction():
    assert extract(rows()[0]).predicate == 'interacts_with'


def test_a_link_to_itself_is_kept():
    relation = extract(rows()[3])
    assert relation.predicate == 'interacts_with'
    assert relation.subject_entity_key == relation.object_entity_key


def test_an_unknown_effect_raises():
    with pytest.raises(ValueError, match='Unknown NRF2ome effect'):
        resource.interactions.mapper({**rows()[0], 'effect': '2'})
