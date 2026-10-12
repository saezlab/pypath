"""Source contract of the Autophagy Regulatory Network input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.arn import resource

# Lines of arn_curated.csv: a directed link without effect, a stimulation with
# four references, an inhibition, an undirected link and an undirected link
# with a sign. The last has a reference that is not a PubMed ID added to it.
LINES = (
    'O15350,O95352,0,1,0,0,1,19001857,\n'
    'O75385,O75143,1,1,1,1,1,20562859:19225151:19287211:19211835,\n'
    'O15350,Q9P2Y5,0,1,-1,0,1,19001857,\n'
    'O75385,Q6ZNE5,1,0,0,1,1,23685627,\n'
    'O95166,O14641,1,0,-1,1,0,0123:20639871,\n'
)


def rows(text=LINES):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(resource.interactions._raw_parser(opener))


def extract(row):
    record = resource.interactions.mapper(row)
    assert record is not None
    out = SilverExtractor('arn', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    return out.relations[0]


def terms(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_reader_names_the_columns():
    parsed = rows()
    assert len(parsed) == 5
    assert parsed[1]['source'] == 'O75385'
    assert parsed[1]['effect'] == '1'
    assert parsed[1]['references'].count(':') == 3


def test_reader_rejects_a_line_with_another_number_of_fields_and_an_empty_file():
    with pytest.raises(ValueError, match='fields'):
        rows('O15350,O95352,0,1,0\n')
    with pytest.raises(ValueError, match='no rows'):
        rows('')


def test_a_directed_link_without_effect_has_no_sign():
    relation = extract(rows()[0])
    assert relation.predicate == 'affects'
    assert terms(relation, 'object_direction_qualifier') == []
    assert terms(relation, 'object_aspect_qualifier') == []


def test_stimulation_and_inhibition_are_signed_effects():
    parsed = rows()
    assert terms(extract(parsed[1]), 'object_direction_qualifier') == ['increased']
    assert terms(extract(parsed[2]), 'object_direction_qualifier') == ['decreased']


def test_references_are_split_on_colons():
    assert terms(extract(rows()[1]), 'publications') == [
        'PMID:20562859',
        'PMID:19225151',
        'PMID:19287211',
        'PMID:19211835',
    ]


def test_an_undirected_link_is_an_interaction_and_its_sign_is_not_carried():
    for row in rows()[3:]:
        relation = extract(row)
        assert relation.predicate == 'interacts_with'
        assert terms(relation, 'object_direction_qualifier') == []


def test_only_a_pubmed_id_is_a_publication():
    # "0123" has a leading zero.
    assert terms(extract(rows()[4]), 'publications') == ['PMID:20639871']


def test_uniprot_accessions_and_the_taxon():
    record = resource.interactions.mapper(rows()[0])
    assert (Namespace.UNIPROT, 'O15350') in {(i.type, i.value) for i in record.subject.identifiers}


def test_an_unknown_flag_raises():
    row = {**rows()[0], 'effect': '2'}
    with pytest.raises(ValueError, match='Unknown ARN effect'):
        resource.interactions.mapper(row)
    row = {**rows()[0], 'is_directed': '2'}
    with pytest.raises(ValueError, match='is_directed'):
        resource.interactions.mapper(row)
