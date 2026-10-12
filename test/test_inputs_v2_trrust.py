"""Source-specific Biolink contract for the TRRUST input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import trrust as mod

# Lines copied from trrust_rawdata.human.tsv.
HUMAN = (
    'AATF\tBAX\tRepression\t22909821\n'
    'ABL1\tBAX\tActivation\t11753601\n'
    'E2F1\tE2F1\tUnknown\t7958836\n'
    'STAT3\tVEGFA\tUnknown\t23104276\n'
    'STAT3\tVEGFA\tActivation\t11960372;15735682;17169602;17533737;19284568;19880237;21278793;21418908\n'
    'STAT3\tVEGFA\tRepression\t17533737;19880237;22313388\n'
    'TP53\tFLT1\tUnknown\t1724219e\n'
)
# Lines copied from trrust_rawdata.mouse.tsv.
MOUSE = (
    'Trp53\tCdkn1a\tRepression\t15467370\n'
    'Gata1\tSpi1\tRepression\t19491391\n'
    'Gata1\tSpi1\tRepression\t19491391\n'
)


def rows(dataset='human_interactions', text=HUMAN):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(mod.resource.datasets()[dataset]._raw_parser(opener))


def extract(row, dataset='human_interactions'):
    record = mod.resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor('trrust', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_module_has_one_dataset_per_species():
    assert list(mod.resource.datasets()) == [
        'human_interactions',
        'mouse_interactions',
    ]


def test_reader_names_the_columns_and_adds_the_species():
    parsed = rows()
    assert len(parsed) == 7
    assert parsed[0] == {
        'tf': 'AATF',
        'target': 'BAX',
        'mode': 'Repression',
        'pmids': '22909821',
        'taxon_id': '9606',
    }
    assert {r['taxon_id'] for r in rows('mouse_interactions', MOUSE)} == {'10090'}


def test_reader_rejects_an_empty_file():
    with pytest.raises(ValueError, match='no rows'):
        rows(text='')


def test_reader_rejects_a_line_with_another_number_of_fields():
    with pytest.raises(ValueError, match='fields'):
        rows(text='AATF\tBAX\tRepression\n')


def test_activation_is_increased_expression():
    _, out = extract(rows()[1])
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'object_aspect_qualifier') == ['expression']
    assert annotations(relation, 'publications') == ['PMID:11753601']
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'protein',
        Namespace.GENESYMBOL,
        'ABL1',
    )
    assert (target.entity_type, target.namespace, target.identifier) == (
        'gene',
        Namespace.GENESYMBOL,
        'BAX',
    )
    assert subject.taxon == target.taxon == '9606'


def test_repression_is_decreased_expression():
    _, out = extract(rows()[0])
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_direction_qualifier') == ['decreased']
    assert annotations(relation, 'object_aspect_qualifier') == ['expression']


@pytest.mark.parametrize('index', [2, 3, 6])
def test_unknown_mode_gets_no_direction(index):
    record, out = extract(rows()[index])
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')
    assert all(
        'qualifier' not in slot_name(a.term)
        for a in record.annotations or ()
    )


def test_conflicting_modes_of_a_pair_stay_separate_statements():
    out = SilverExtractor('trrust', 'human_interactions')
    stat3_vegfa = rows()[3:6]
    for number, row in enumerate(stat3_vegfa, 1):
        record = mod.resource.human_interactions.mapper(row)
        out.process_record(record, row, f'fixture:{number}', number)
    assert len(out.relations) == 3
    assert [
        annotations(relation, 'object_direction_qualifier')
        for relation in out.relations
    ] == [[], ['increased'], ['decreased']]
    # Each statement keeps its own articles.
    assert annotations(out.relations[2], 'publications') == [
        'PMID:17533737',
        'PMID:19880237',
        'PMID:22313388',
    ]
    assert len(out.entities) == 2


def test_self_regulation_keeps_protein_and_gene_apart():
    _, out = extract(rows()[2])
    relation = out.relations[0]
    assert relation.subject_entity_key != relation.object_entity_key
    assert {e.entity_type for e in out.entities.values()} == {'protein', 'gene'}
    assert {e.identifier for e in out.entities.values()} == {'E2F1'}


def test_mouse_file_gives_mouse_participants():
    parsed = rows('mouse_interactions', MOUSE)
    _, out = extract(parsed[0], 'mouse_interactions')
    relation = out.relations[0]
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert (subject.identifier, target.identifier) == ('Trp53', 'Cdkn1a')
    assert subject.taxon == target.taxon == '10090'
    assert annotations(relation, 'object_direction_qualifier') == ['decreased']


@pytest.mark.parametrize('numbers', ['0123456;123456789', '0123456', '123456789'])
def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(numbers):
    _, out = extract(dict(rows()[0], pmids=numbers))
    assert not annotations(out.relations[0], 'publications')


def test_malformed_pubmed_id_is_not_a_publication():
    _, out = extract(rows()[6])
    assert not annotations(out.relations[0], 'publications')


@pytest.mark.parametrize('missing', ['tf', 'target'])
def test_row_without_a_participant_emits_no_relation(missing):
    row = dict(rows()[0], **{missing: ''})
    assert mod.resource.human_interactions.mapper(row) is None


@pytest.mark.parametrize('mode', ['activation', 'Inhibition', '', None])
def test_unknown_mode_value_is_rejected(mode):
    with pytest.raises(ValueError, match='Unknown TRRUST mode of regulation'):
        mod.resource.human_interactions.mapper(dict(rows()[0], mode=mode))
