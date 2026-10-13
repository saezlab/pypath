"""Phobius: parser on rows of the real file, mappers and extractor."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import phobius
from pypath.inputs_v2.parsers import phobius as phobius_parser

HEADER = 'Source\tUniprot_AC\tUniprot_ID\tTM\tSP\tPREDICTION'
# Rows of phobius_human_proteome.tsv as written (the file has CRLF line ends).
# A transmembrane protein with a signal peptide.
DQB1 = 'sp\tP01920\tDQB1_HUMAN\t1\tY\tn15-25c32/33o231-251i'
# A signal peptide only.
FA9 = 'sp\tP00740\tFA9_HUMAN\t0\tY\tn13-24c28/29o'
# Neither prediction.
ARY1 = 'sp\tP18440\tARY1_HUMAN\t0\t0\to'
# Helices only, with the entry name split by a tab after "_H" (the TM column
# is padded with spaces).
CP2D7 = 'sp\tA0A087X1C5\tCP2D7_H\tUMAN       2\t0\to6-23i302-325o'
# A signal peptide only, split after "_HU".
TVA2 = 'sp\tA0A0B4J234\tTVA2_HU\tMAN        0\tY\tn6-16c25/26o'
# Three helices and a signal peptide, split after "_H".
TM271 = (
    'sp\tA0A286YF58\tTM271_H\tUMAN       3\tY\tn8-20c28/29o52-72i120-141o218-239i'
)


def table(*rows):
    text = '\r\n'.join([HEADER, *rows]) + '\r\n'
    return SimpleNamespace(result=text)


def proteins(*rows):
    return list(phobius_parser.iter_proteins(table(*rows)))


def predictions(*rows):
    return list(phobius_parser.iter_predictions(table(*rows)))


def extract(mapper, rows, dataset):
    out = SilverExtractor('phobius', dataset)
    for index, row in enumerate(rows):
        out.process_record(mapper(row), row, f'fixture:{index}', index)
    return out


def statements(out):
    return {
        (
            out.entities[r.subject_entity_key].identifier,
            r.predicate,
            out.entities[r.object_entity_key].namespace,
            out.entities[r.object_entity_key].identifier,
        )
        for r in out.relations
    }


def test_parser_reads_rows_of_the_file():
    rows = proteins(DQB1, FA9, ARY1)
    assert rows[0] == {
        'accession': 'P01920', 'entry_name': 'DQB1_HUMAN', 'helices': 1,
        'signal_peptide': True, 'topology': 'n15-25c32/33o231-251i', 'repaired': False,
    }
    assert [(r['helices'], r['signal_peptide']) for r in rows] == [
        (1, True), (0, True), (0, False),
    ]


def test_parser_repairs_an_entry_name_split_by_a_tab():
    rows = proteins(CP2D7, TVA2, TM271, ARY1)
    assert [(r['entry_name'], r['helices'], r['repaired']) for r in rows] == [
        ('CP2D7_HUMAN', 2, True),
        ('TVA2_HUMAN', 0, True),
        ('TM271_HUMAN', 3, True),
        ('ARY1_HUMAN', 0, False),
    ]
    # The prediction of a repaired row is read from its own columns.
    assert rows[1]['signal_peptide'] and rows[1]['topology'] == 'n6-16c25/26o'


def test_parser_raises_on_a_row_the_rule_does_not_repair():
    # A fourth column that is neither a number nor the rest of the name.
    with pytest.raises(ValueError, match='repair rule'):
        proteins('sp\tP18440\tARY1_HUMAN\tx\t0\to')
    # The tail does not complete "HUMAN".
    with pytest.raises(ValueError, match='repair rule'):
        proteins('sp\tA0A087X1C5\tCP2D7_H\tMOUSE       2\t0\to6-23i302-325o')
    # A row with the wrong number of fields.
    with pytest.raises(ValueError, match='fields'):
        proteins('sp\tP18440\tARY1_HUMAN\t0\t0')
    with pytest.raises(ValueError, match='lacks expected columns'):
        list(phobius_parser.iter_proteins(
            SimpleNamespace(result='Source\tUniprot_AC\tTM\r\nsp\tP18440\t0\r\n')
        ))
    with pytest.raises(ValueError, match='no rows'):
        proteins()


def test_parser_raises_on_a_protein_that_is_not_human():
    with pytest.raises(ValueError, match='_HUMAN'):
        proteins('sp\tP18440\tARY1_MOUSE\t0\t0\to')


def test_parser_raises_when_tm_is_not_the_number_of_helices():
    with pytest.raises(ValueError, match='number of helices'):
        proteins(DQB1.replace('\t1\tY\t', '\t2\tY\t'))
    with pytest.raises(ValueError, match='number of helices'):
        proteins(ARY1.replace('\t0\t0\t', '\t1\t0\t'))


def test_parser_raises_when_sp_does_not_agree_with_the_topology():
    with pytest.raises(ValueError, match='SP'):
        proteins(FA9.replace('\t0\tY\t', '\t0\t0\t'))
    with pytest.raises(ValueError, match='SP'):
        proteins(ARY1.replace('\t0\t0\t', '\t0\tY\t'))
    with pytest.raises(ValueError, match='SP value'):
        proteins(ARY1.replace('\t0\t0\t', '\t0\tN\t'))


def test_parser_raises_on_an_unknown_topology_and_a_repeated_accession():
    with pytest.raises(ValueError, match='topology'):
        proteins(ARY1.replace('\to', '\tx'))
    with pytest.raises(ValueError, match='twice'):
        proteins(ARY1, ARY1)


def test_predictions_give_one_row_per_statement_and_none_for_neither():
    rows = predictions(DQB1, FA9, ARY1, CP2D7)
    assert [(r['accession'], r['term']) for r in rows] == [
        ('P01920', 'transmembrane'),
        ('P01920', 'signal_peptide'),
        ('P00740', 'signal_peptide'),
        ('A0A087X1C5', 'transmembrane'),
    ]


def test_a_transmembrane_protein_with_a_signal_peptide_has_two_statements():
    out = extract(phobius.map_prediction, predictions(DQB1), 'predictions')
    assert statements(out) == {
        ('P01920', 'associated_with', 'phobius_prediction', 'transmembrane'),
        ('P01920', 'associated_with', 'phobius_prediction', 'signal_peptide'),
    }
    protein = next(e for e in out.entities.values() if e.entity_type == 'protein')
    assert (protein.namespace, protein.identifier) == (Namespace.UNIPROT, 'P01920')
    assert {i['id'] for i in protein.identifiers if i['ns'] == Namespace.UNIPROT_ENTRY} == {
        'DQB1_HUMAN'
    }
    # The entry names end in _HUMAN and the file is a run on the human proteome.
    assert protein.taxon == '9606'
    by_object = {
        out.entities[r.object_entity_key].identifier: {
            a['term']: (a['quantity']['has_numeric_value'] if a.get('quantity') else a['value'])
            for a in r.annotations if a['term'].startswith('phobius:')
        }
        for r in out.relations
    }
    # The count is on the transmembrane statement only; the topology on both.
    assert by_object == {
        'transmembrane': {
            # A count: a measurement, not the text '1'.
            'phobius:transmembrane_helices': 1,
            'phobius:topology': 'n15-25c32/33o231-251i',
        },
        'signal_peptide': {'phobius:topology': 'n15-25c32/33o231-251i'},
    }


def test_helices_only_and_signal_peptide_only():
    out = extract(phobius.map_prediction, predictions(CP2D7, FA9), 'predictions')
    assert statements(out) == {
        ('A0A087X1C5', 'associated_with', 'phobius_prediction', 'transmembrane'),
        ('P00740', 'associated_with', 'phobius_prediction', 'signal_peptide'),
    }
    for relation in out.relations:
        terms = {a['term'] for a in relation.annotations}
        if out.entities[relation.object_entity_key].identifier == 'signal_peptide':
            assert 'phobius:transmembrane_helices' not in terms
        else:
            assert 'phobius:transmembrane_helices' in terms
    # Every annotation has a value.
    assert all(a['value'] or a['quantity'] for r in out.relations for a in r.annotations)


def test_a_repaired_row_is_mapped_with_its_whole_entry_name():
    out = extract(phobius.map_prediction, predictions(TM271), 'predictions')
    protein = next(e for e in out.entities.values() if e.entity_type == 'protein')
    assert {i['id'] for i in protein.identifiers if i['ns'] == Namespace.UNIPROT_ENTRY} == {
        'TM271_HUMAN'
    }
    helices = [
        a['quantity'] for r in out.relations for a in r.annotations
        if a['term'] == 'phobius:transmembrane_helices'
    ]
    assert [(h['has_numeric_value'], h['has_unit'], h['source_field']) for h in helices] == [
        (3, None, 'TM')
    ]


def test_a_protein_with_neither_prediction_gives_no_statement():
    assert predictions(ARY1) == []
    out = extract(phobius.map_prediction, predictions(ARY1, FA9), 'predictions')
    assert {s[0] for s in statements(out)} == {'P00740'}
    assert 'ARY1_HUMAN' not in {
        i['id'] for e in out.entities.values() for i in e.identifiers
    }


def test_terms_dataset_has_names_and_descriptions():
    rows = list(phobius_parser.iter_terms())
    out = extract(phobius.map_term, rows, 'terms')
    by_id = {(e.namespace, e.identifier): e for e in out.entities.values()}
    assert set(by_id) == {
        ('phobius_prediction', 'transmembrane'),
        ('phobius_prediction', 'signal_peptide'),
    }
    transmembrane = by_id[('phobius_prediction', 'transmembrane')]
    assert transmembrane.entity_type == 'ontology_class'
    assert transmembrane.label == 'predicted transmembrane protein'
    assert by_id[('phobius_prediction', 'signal_peptide')].label == 'predicted signal peptide'
    assert all(
        e.annotations and e.annotations[0]['term'] == 'description'
        for e in by_id.values()
    )
    with pytest.raises(KeyError):
        phobius.map_term({'term': 'cytoplasmic'})


def test_statement_objects_are_the_entities_of_the_terms_dataset():
    term_out = extract(phobius.map_term, phobius_parser.iter_terms(), 'terms')
    out = extract(phobius.map_prediction, predictions(DQB1, FA9, CP2D7), 'predictions')
    term_keys = {
        (e.entity_type, e.namespace, e.identifier) for e in term_out.entities.values()
    }
    objects = {
        (o.entity_type, o.namespace, o.identifier)
        for o in (out.entities[r.object_entity_key] for r in out.relations)
    }
    assert objects == term_keys
    with pytest.raises(ValueError, match='Unknown Phobius term'):
        phobius.map_prediction(dict(predictions(FA9)[0], term='cytoplasmic'))


def test_resource_datasets():
    assert set(phobius.resource.datasets()) == {'terms', 'predictions'}
    assert phobius.config.pubmed == '15111065'
