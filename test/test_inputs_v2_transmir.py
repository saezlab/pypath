"""Source-specific Biolink contract for the TransmiR input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import transmir as mod
from pypath.inputs_v2.parsers.transmir import COLUMNS, iter_literature

# Lines copied from db_literaturev3.txt of TransmiR v3.0 (tab separated, no
# header, the TSS and site columns empty).
LINES = {
    'activation': 'TFEB\thsa-mir-16\t\t\tActivation\t30591554\tliterature\thuman endothelial cells\tH.sapiens',
    'repression': 'SOX9\thsa-mir-96\t\t\tRepression\t37717559\tliterature\tHEK-293T\tH.sapiens',
    'regulation': 'TCF4\thsa-mir-484\t\t\tRegulation\t28921929\tliterature\tHCT-116\tH.sapiens',
    'activation_feedback': 'AHR\thsa-mir-106b\t\t\tActivation(feedback)\t24798859\tliterature\t\tH.sapiens',
    'repression_feedback': 'AR\thsa-mir-190a\t\t\tRepression(feedback)\t26314494\tliterature\t\tH.sapiens',
    'double_negative': 'LIN28\thsa-let-7a\t\t\tRegulation(Double-Negative Feedback Loop)\t21045151\tliterature\t\tH.sapiens',
    'loop': 'PSMD9\thsa-mir-221\t\t\tregulatory loop\t20203269\tliterature\t\tH.sapiens',
    'activation_loop': 'SRSF1\thsa-mir-221\t\t\tActivation(a negative feedback loop)\t20385090\tliterature\t\tH.sapiens',
    'plant': 'AGL15\tcsi-MIR156a\t\t\tActivation\t35661206\tliterature\tembryogenic callus\tC.sinensis',
    'herpesvirus': 'HHV8gp76\thsa-mir-146a\t\t\tActivation\t20023696\tliterature\t\tHuman herpesvirus 8 type M',
    'ebv': 'PDLIM7\tmir-146a\t\t\tActivation\t18057241\tliterature\t\tEpstein-Barr virus (EBV)',
    'space': 'Smad5 \tmmu-mir-23a\t\t\tActivation\t25476774\tliterature\t\tM.musculus',
    'alpha': 'HIF-1α\trno-mir-322\t\t\tActivation\t37566526\tliterature\tH9c2\tR.norvegicus',
    'other_species_mirna': 'SALL4\tmmu-mir-146a\t\t\tActivation\t31143524\tliterature\tHepa1-6\tH.sapiens',
}


def parse(*keys, text=None):
    text = text if text is not None else ''.join(LINES[k] + '\n' for k in keys)
    return list(iter_literature(SimpleNamespace(result=io.StringIO(text))))


def row(key):
    return parse(key)[0]


def extract(key):
    parsed = row(key)
    record = mod.interactions_schema(parsed)
    assert record is not None
    out = SilverExtractor('transmir', 'literature')
    out.process_record(record, parsed, 'fixture:1', 1)
    return record, out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_module_has_the_literature_dataset():
    assert list(mod.resource.datasets()) == ['literature']


def test_parser_names_the_columns_and_strips_cells():
    parsed = row('space')
    assert tuple(parsed) == COLUMNS
    assert parsed == {
        'tf': 'Smad5',
        'mirna': 'mmu-mir-23a',
        'tss': '',
        'site': '',
        'action': 'Activation',
        'pmid': '25476774',
        'evidence': 'literature',
        'tissue': '',
        'species': 'M.musculus',
    }
    assert row('alpha')['tf'] == 'HIF-1α'


@pytest.mark.parametrize(
    'text',
    ['', 'TFEB\thsa-mir-16\tActivation\t30591554\n'],
)
def test_parser_rejects_an_empty_file_or_a_wrong_table(text):
    with pytest.raises(ValueError):
        parse(text=text)


def test_activation_is_increased_expression_of_the_microrna():
    record, out = extract('activation')
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'object_aspect_qualifier') == ['expression']
    assert annotations(relation, 'publications') == ['PMID:30591554']
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'protein',
        Namespace.NAME,
        'TFEB',
    )
    assert (target.entity_type, target.namespace, target.identifier) == (
        'microrna',
        Namespace.NAME,
        'hsa-mir-16',
    )
    assert subject.taxon == target.taxon == '9606'


def test_repression_is_decreased_expression_of_the_microrna():
    _, out = extract('repression')
    relation = out.relations[0]
    assert annotations(relation, 'object_direction_qualifier') == ['decreased']
    assert annotations(relation, 'object_aspect_qualifier') == ['expression']


@pytest.mark.parametrize(
    'key,direction',
    [
        ('activation_feedback', 'increased'),
        ('repression_feedback', 'decreased'),
        ('activation_loop', 'increased'),
    ],
)
def test_remark_in_brackets_does_not_change_the_direction(key, direction):
    _, out = extract(key)
    relation = out.relations[0]
    assert annotations(relation, 'object_direction_qualifier') == [direction]
    assert annotations(relation, 'object_aspect_qualifier') == ['expression']


@pytest.mark.parametrize('key', ['regulation', 'double_negative', 'loop'])
def test_regulation_without_a_sign_is_a_bare_affects(key):
    record, out = extract(key)
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')
    assert all('qualifier' not in str(a.term) for a in record.annotations)


def test_every_observed_action_type_is_in_the_crosswalk():
    assert set(mod._ACTIONS) == {
        'Activation',
        'Activation(feedback)',
        'Activation(a negative feedback loop)',
        'Activation(negative regulatory loop)',
        'Repression',
        'Repression(feedback)',
        'Regulation',
        'Regulation(feedback)',
        'Regulation(Double-Negative Feedback Loop)',
        'auto-regulatory feedback circuit',
        'autoregulatory negative feedback loop(feedback)',
        'A negative regulatory loop(feedback)',
        'regulatory loop',
        'regulatory network',
    }


def test_unknown_action_type_is_rejected():
    with pytest.raises(ValueError, match='Unknown TransmiR action type'):
        mod.interactions_schema(dict(row('activation'), action='Inhibition'))


def test_species_of_the_regulation_is_the_taxon_of_both_participants():
    _, out = extract('plant')
    relation = out.relations[0]
    assert out.entities[relation.subject_entity_key].taxon == '2711'
    assert out.entities[relation.object_entity_key].taxon == '2711'


def test_species_column_is_repeated_even_if_the_microrna_prefix_differs():
    # A human regulation of a mouse microRNA name: the source says human.
    _, out = extract('other_species_mirna')
    relation = out.relations[0]
    assert out.entities[relation.object_entity_key].taxon == '9606'


@pytest.mark.parametrize('key', ['herpesvirus', 'ebv'])
def test_virus_regulation_gets_no_taxon(key):
    record, out = extract(key)
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert out.entities[relation.subject_entity_key].taxon is None
    assert out.entities[relation.object_entity_key].taxon is None


def test_unknown_species_is_rejected():
    with pytest.raises(ValueError, match='Unknown TransmiR species'):
        mod.interactions_schema(dict(row('activation'), species='P.troglodytes'))


def test_names_are_kept_as_written_without_an_accession():
    record, out = extract('alpha')
    relation = out.relations[0]
    assert out.entities[relation.subject_entity_key].identifier == 'HIF-1α'
    # No miRBase accession is made up from a name.
    assert all(i.type == Namespace.NAME for i in record.object.identifiers)


@pytest.mark.parametrize('missing', ['tf', 'mirna'])
def test_row_without_a_participant_emits_no_relation(missing):
    assert mod.interactions_schema(dict(row('activation'), **{missing: ''})) is None


@pytest.mark.parametrize('number', ['0', '0123456', '123456789'])
def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(number):
    parsed = dict(row('activation'), pmid=number)
    out = SilverExtractor('transmir', 'literature')
    out.process_record(mod.interactions_schema(parsed), parsed, 'fixture:1', 1)
    assert not annotations(out.relations[0], 'publications')


def test_publication_that_is_not_a_pubmed_id_is_left_out():
    _, out = extract('activation')
    parsed = dict(row('activation'), pmid='not available')
    out = SilverExtractor('transmir', 'literature')
    out.process_record(
        mod.interactions_schema(parsed), parsed, 'fixture:1', 1
    )
    assert not annotations(out.relations[0], 'publications')
