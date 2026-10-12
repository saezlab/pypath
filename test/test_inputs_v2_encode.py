"""Source-specific Biolink contract for the ENCODE-Nets input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import encode as mod
from pypath.inputs_v2.parsers import encode as parser

# Lines copied from enets2.Proximal_filtered.txt (space separated),
# enets3.Distal.txt (tab separated, Windows line ends) and
# enets10.TF-miRNA.txt (space separated).
PROXIMAL = (
    'CTCF proximal_filtered MYBPC3\n'
    'MEF2A proximal_filtered MEF2A\n'
    'PRDM1 proximal_filtered APOL2\n'
)
DISTAL = (
    'ATF3\tdistal\tCD7\r\n'
    'POL2\tdistal\tA2ML1\r\n'
    'POU2F2\tdistal\tPOU2F2\r\n'
)
TF_MIRNA = (
    'BDP1 (TF-miRNA) hsa-miR-19b-1*\n'
    'ESR1 (TF-miRNA) hsa-miR-23b*\n'
    'hsa-miR-137 (miRNA-TF) TCF12\n'
    'hsa-miR-124 (miRNA-TF) BCL11A\n'
)


def rows(dataset, text):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(mod.resource.datasets()[dataset]._raw_parser(opener))


def extract(dataset, row):
    record = mod.resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor('encode', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_module_has_one_dataset_per_network():
    assert list(mod.resource.datasets()) == [
        'proximal_filtered',
        'distal',
        'tf_mirna',
    ]


def test_reader_names_the_columns_of_each_file():
    assert rows('proximal_filtered', PROXIMAL)[0] == {
        'source': 'CTCF',
        'edge_type': 'proximal_filtered',
        'target': 'MYBPC3',
    }
    # The tab separated file has Windows line ends.
    assert rows('distal', DISTAL)[0] == {
        'source': 'ATF3',
        'edge_type': 'distal',
        'target': 'CD7',
    }
    assert len(rows('proximal_filtered', PROXIMAL)) == 3
    assert len(rows('distal', DISTAL)) == 3


def test_reader_of_the_microrna_file_keeps_the_factor_to_microrna_edges():
    parsed = rows('tf_mirna', TF_MIRNA)
    assert [row['source'] for row in parsed] == ['BDP1', 'ESR1']
    assert all(row['edge_type'] == '(TF-miRNA)' for row in parsed)


def test_microrna_to_factor_edges_can_be_read_apart():
    opener = SimpleNamespace(result=io.StringIO(TF_MIRNA))
    other = list(
        parser.iter_edges(
            opener,
            delimiter=' ',
            edge_type=parser.MIRNA_TF,
            other_edge_types=(parser.TF_MIRNA,),
        )
    )
    assert [(r['source'], r['target']) for r in other] == [
        ('hsa-miR-137', 'TCF12'),
        ('hsa-miR-124', 'BCL11A'),
    ]


@pytest.mark.parametrize(
    'dataset,text',
    [
        ('proximal_filtered', ''),
        ('proximal_filtered', 'CTCF proximal_filtered\n'),
        ('proximal_filtered', 'CTCF distal MYBPC3\n'),
        ('distal', 'ATF3 distal CD7\n'),
        ('tf_mirna', 'BDP1 (miRNA-TF) hsa-miR-19b-1*\n'),
        ('tf_mirna', 'BDP1 (TF-TF) ESR1\n'),
    ],
)
def test_reader_rejects_an_empty_file_a_wrong_table_or_edge_type(dataset, text):
    with pytest.raises(ValueError):
        rows(dataset, text)


def test_proximal_edge_is_binding_from_chip_seq_not_an_effect():
    record, out = extract('proximal_filtered', rows('proximal_filtered', PROXIMAL)[0])
    relation = out.relations[0]
    # Binding near a promoter shows no effect on expression, so not `affects`.
    assert relation.predicate == 'interacts_with'
    assert annotations(relation, 'encode:edge_type') == ['proximal_filtered']
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')
    assert all('qualifier' not in str(a.term) for a in record.annotations)
    assert annotations(relation, 'has_evidence_of_type') == ['ECO:0000229']
    assert annotations(relation, 'publications') == ['PMID:22955619']
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'protein',
        Namespace.GENESYMBOL,
        'CTCF',
    )
    assert (target.entity_type, target.namespace, target.identifier) == (
        'gene',
        Namespace.GENESYMBOL,
        'MYBPC3',
    )
    assert subject.taxon == target.taxon == '9606'


def test_distal_edge_does_not_claim_chip_seq_evidence():
    record, out = extract('distal', rows('distal', DISTAL)[0])
    relation = out.relations[0]
    # Binding at a distal element linked to the gene by chromatin features and
    # expression: a binding, not an effect.
    assert relation.predicate == 'interacts_with'
    assert annotations(relation, 'encode:edge_type') == ['distal']
    assert not annotations(relation, 'has_evidence_of_type')
    assert annotations(relation, 'publications') == ['PMID:22955619']
    target = out.entities[relation.object_entity_key]
    assert (target.entity_type, target.identifier) == ('gene', 'CD7')


def test_factor_names_are_kept_as_written():
    # POL2 is RNA polymerase II, listed with the factors of the file.
    _, out = extract('distal', rows('distal', DISTAL)[1])
    subject = out.entities[out.relations[0].subject_entity_key]
    assert (subject.entity_type, subject.identifier) == ('protein', 'POL2')


def test_factor_regulating_its_own_gene_is_kept():
    _, out = extract('proximal_filtered', rows('proximal_filtered', PROXIMAL)[1])
    relation = out.relations[0]
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert (subject.entity_type, target.entity_type) == ('protein', 'gene')
    assert subject.identifier == target.identifier == 'MEF2A'


def test_factor_microrna_edge_targets_a_microrna_by_name():
    record, out = extract('tf_mirna', rows('tf_mirna', TF_MIRNA)[0])
    relation = out.relations[0]
    assert relation.predicate == 'interacts_with'
    assert annotations(relation, 'encode:edge_type') == ['(TF-miRNA)']
    assert not annotations(relation, 'object_direction_qualifier')
    assert annotations(relation, 'has_evidence_of_type') == ['ECO:0000229']
    target = out.entities[relation.object_entity_key]
    assert (target.entity_type, target.namespace, target.identifier) == (
        'microrna',
        Namespace.NAME,
        'hsa-miR-19b-1*',
    )
    assert target.taxon == '9606'
    # The name is not turned into a miRBase accession.
    assert all(i.type == Namespace.NAME for i in record.object.identifiers)


@pytest.mark.parametrize('dataset', ['proximal_filtered', 'distal', 'tf_mirna'])
@pytest.mark.parametrize('missing', ['source', 'target'])
def test_row_without_a_participant_emits_no_relation(dataset, missing):
    row = {'source': 'CTCF', 'edge_type': 'x', 'target': 'MYBPC3'}
    row[missing] = ''
    assert mod.resource.datasets()[dataset].mapper(row) is None


def test_no_edge_of_any_network_is_an_effect_on_expression():
    for dataset, text in [
        ('proximal_filtered', PROXIMAL),
        ('distal', DISTAL),
        ('tf_mirna', TF_MIRNA),
    ]:
        for row in rows(dataset, text):
            record, out = extract(dataset, row)
            assert {r.predicate for r in out.relations} == {'interacts_with'}
            assert all('qualifier' not in str(a.term) for a in record.annotations)
