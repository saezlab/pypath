"""Source contract of the Laudanna Lab signaling network input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.laudanna import resource

# Lines of SignalingFlow.EA: its title, two edges, an edge with a second line
# of its own and a line without the second node.
FLOW = (
    'Signaling_flow\n'
    'F2R (pp) PAWR = state_change\n'
    'ATP8A2 (pp) CASP9 = state_change\n'
    'HSPA4 (pp) HSPA4 = complex_formation\n'
    'PKMYT1 (pp)  = state_change\n'
)
# Lines of SignalingDirection.EA.
DIRECTION = (
    'SignalingDirection\n'
    'EDNRA (pp) NBN = activation\n'
    'TGFB1 (pp) TGFB1 = inhibition\n'
    'DEK (pp) UBE3A = docking\n'
)


def rows(dataset, text):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(resource.datasets()[dataset]._raw_parser(opener))


def extract(dataset, row):
    record = resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor('laudanna', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return out


def terms(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_module_has_a_dataset_per_file():
    assert list(resource.datasets()) == ['signaling_flow', 'signaling_direction']


def test_reader_yields_the_edges_and_the_line_without_a_target():
    parsed = rows('signaling_flow', FLOW)
    assert parsed[0] == {'source': 'F2R', 'target': 'PAWR', 'value': 'state_change'}
    assert parsed[-1] == {'source': 'PKMYT1', 'target': None, 'value': 'state_change'}
    assert len(parsed) == 4


def test_reader_rejects_a_file_with_another_title():
    with pytest.raises(ValueError, match='title'):
        rows('signaling_flow', DIRECTION)


def test_reader_rejects_an_unknown_value_and_a_malformed_line():
    with pytest.raises(ValueError, match='Unknown value'):
        rows('signaling_flow', 'Signaling_flow\nA (pp) B = phosphorylation\n')
    with pytest.raises(ValueError, match='malformed'):
        rows('signaling_flow', 'Signaling_flow\nA - B = state_change\n')
    with pytest.raises(ValueError, match='no edges'):
        rows('signaling_flow', 'Signaling_flow\n')


def test_state_change_is_an_effect_without_sign_or_aspect():
    out = extract('signaling_flow', rows('signaling_flow', FLOW)[0])
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert terms(relation, 'object_direction_qualifier') == []
    assert terms(relation, 'object_aspect_qualifier') == []
    assert {a['value'] for a in relation.annotations if a['term'] == 'in_taxon'} == {
        'NCBITaxon:9606'
    }


def test_complex_formation_is_an_interaction():
    out = extract('signaling_flow', rows('signaling_flow', FLOW)[2])
    assert out.relations[0].predicate == 'interacts_with'


def test_activation_and_inhibition_are_signed_effects():
    parsed = rows('signaling_direction', DIRECTION)
    up = extract('signaling_direction', parsed[0]).relations[0]
    down = extract('signaling_direction', parsed[1]).relations[0]
    assert (up.predicate, terms(up, 'object_direction_qualifier')) == ('affects', ['increased'])
    assert (down.predicate, terms(down, 'object_direction_qualifier')) == ('affects', ['decreased'])
    # The file does not say what is changed.
    assert terms(up, 'object_aspect_qualifier') == []


def test_docking_has_no_sign():
    relation = extract('signaling_direction', rows('signaling_direction', DIRECTION)[2]).relations[0]
    assert relation.predicate == 'interacts_with'
    assert terms(relation, 'object_direction_qualifier') == []


def test_gene_symbols_identify_the_proteins():
    record = resource.signaling_direction.mapper(rows('signaling_direction', DIRECTION)[0])
    assert (Namespace.GENESYMBOL, 'EDNRA') in {
        (i.type, i.value) for i in record.subject.identifiers
    }


def test_a_line_without_a_target_gives_no_relation():
    row = rows('signaling_flow', FLOW)[-1]
    assert resource.signaling_flow.mapper(row) is None


def test_a_value_outside_the_crosswalk_raises_in_the_mapper():
    with pytest.raises(ValueError, match='Unknown Laudanna'):
        resource.signaling_direction.mapper({'source': 'A', 'target': 'B', 'value': 'activity'})
