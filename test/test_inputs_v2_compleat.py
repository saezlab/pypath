"""COMPLEAT complexes: reader, parsed rows and Biolink mapping."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.compleat import COMPLEAT_NAMESPACE, DATABASES, map_complex, resource
from pypath.inputs_v2.parsers.compleat import COLUMNS, iter_compleat_complexes

# Rows of the file, as it is written: CRLF line ends, members separated by
# spaces with one at the end, and a last line ``null``.
ROWS = [
    ('HC1', '5', 'Literature', '13S condensin complex',
     '13S condensin complex;condensin;Condensin I complex;condensin complex;nuclear condensin complex',
     '', 'Affinity Capture-Western', 'nuclear condensin complex', 'CYC2008', 'Yeast', '10811823',
     '9918 10592 64151 23397 10051 '),
    ('HC2543', '3', 'Literature', 'ITGAV-ITGB5-PLAUR complex', 'ITGAV-ITGB5-PLAUR complex', '',
     'MI:0019- coimmunoprecipitation', 'ITGAV-ITGB5-PLAUR complex', 'CORUM', 'Human', '10537314',
     '3685 3693 5329 '),
    ('HC2501', '2', 'Literature', 'USP1-UAF1 complex', 'USP1-UAF1 complex;USP1/UAF1', '',
     'Flag-HA-tagged USP1', 'USP1/UAF1', 'PINdb', 'Human', '18082604', '7398 57599 '),
    # This row has its PubMed ID and method in the wrong columns.
    ('HC3339', '2', 'Literature', 'DNMT3L-DNMT3A', 'DNMT3L-DNMT3A', '15105426', 'complex',
     'DNMT3L-DNMT3A', 'CORUM', 'Human', 'MI:0071- molecular sieving', '1788 29947 '),
    ('HC1332', '7', 'Literature', 'mRNA processing',
     'mRNA processing;mRNA metabolic process;RNA processing',
     '', 'Predicted by HT data (cluster was obtained from CYC2008)', 'HT_SC_Cluster37', 'CYC2008',
     'Yeast', '', '11044 23517 84186 57466 64282 3183 23064 '),
    ('HC4831', '15', 'Predicted', 'protein transport',
     'protein transport;cellular membrane fusion;Golgi vesicle transport',
     '', 'NetworkBlast', '', 'Predicted from PPI network', 'Human & Yeast', '',
     '8674 6845 9482 4905 8773 8673 9341 8417 143187 6293 6810 8775 8675 7316 10228 '),
    ('HC6106', '5', 'Predicted', 'positive regulation of neurogenesis',
     'positive regulation of neurogenesis;protein heterooligomerization', '', 'CFinder', '',
     'Predicted from PPI network', 'Human', '', '22938 3987 55742 3611 988 '),
]


def text(rows=ROWS, last='null'):
    lines = ['\t'.join(row) + '\r\n' for row in rows]
    return ''.join(lines) + (last + '\r\n' if last else '')


def opener(content=None):
    return SimpleNamespace(result={'compleat.tsv': content if content is not None else text()})


def parsed(source=None, content=None):
    return {row['compleat_id']: row for row in iter_compleat_complexes(opener(content), source=source)}


def extract(row, dataset='literature_complexes'):
    record = map_complex(row)
    assert record is not None
    out = SilverExtractor('compleat', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def complexes(out):
    return [e for e in out.entities.values() if e.entity_type == 'macromolecular_complex']


def test_module_imports_with_two_datasets():
    assert set(resource.datasets()) == {'literature_complexes', 'predicted_complexes'}
    assert resource.config.name == 'COMPLEAT'


def test_reader_gives_every_complex_and_skips_the_closing_line():
    rows = parsed()
    assert list(rows) == ['HC1', 'HC2543', 'HC2501', 'HC3339', 'HC1332', 'HC4831', 'HC6106']
    assert set(rows['HC1']) >= set(COLUMNS)
    assert rows['HC1']['members'] == ['9918', '10592', '64151', '23397', '10051']


def test_reader_selects_a_dataset_by_the_source_column():
    assert list(parsed('Predicted')) == ['HC4831', 'HC6106']
    assert list(parsed('Literature')) == ['HC1', 'HC2543', 'HC2501', 'HC3339', 'HC1332']
    with pytest.raises(ValueError, match='Unknown COMPLEAT source'):
        list(iter_compleat_complexes(opener(), source='Curated'))


def test_reader_decodes_publications_and_psimi_methods_only():
    rows = parsed()
    assert rows['HC2543']['publications'] == ['10537314']
    assert rows['HC2543']['detection_methods'] == ['MI:0019']
    # a BioGRID experimental system, free text and a prediction method are no PSI-MI term
    for accession in ('HC1', 'HC2501', 'HC4831'):
        assert rows[accession]['detection_methods'] == []
    # the malformed row: nothing is read from the shifted columns
    assert rows['HC3339']['publications'] == []
    assert rows['HC3339']['detection_methods'] == []


def test_reader_rejects_malformed_files():
    with pytest.raises(ValueError, match='no complexes'):
        list(iter_compleat_complexes(opener(text([]))))
    with pytest.raises(ValueError, match='too few fields'):
        list(iter_compleat_complexes(opener(text(last='HC9\t3'))))
    with pytest.raises(ValueError, match='size does not match'):
        bad = [('HC1', '4', *ROWS[0][2:])]
        list(iter_compleat_complexes(opener(text(bad))))
    with pytest.raises(ValueError, match='Unknown COMPLEAT source'):
        bad = [(*ROWS[0][:2], 'Curated', *ROWS[0][3:])]
        list(iter_compleat_complexes(opener(text(bad))))
    with pytest.raises(ValueError, match='malformed members'):
        bad = [(*ROWS[0][:11], '9918 TP53 ')]
        list(iter_compleat_complexes(opener(text(bad))))
    with pytest.raises(ValueError):
        list(iter_compleat_complexes(opener('<html><body>Error</body></html>')))


def test_curated_complex_has_accession_name_members_and_evidence():
    record, out = extract(parsed()['HC2543'])
    (cplx,) = complexes(out)
    assert cplx.namespace == COMPLEAT_NAMESPACE and cplx.identifier == 'HC2543'
    assert [i['ns'] for i in cplx.identifiers] == [COMPLEAT_NAMESPACE, Namespace.NAME]
    assert cplx.label == 'ITGAV-ITGB5-PLAUR complex'
    terms = {(a['term'], a['value']) for a in cplx.annotations}
    assert ('supporting_data_source', 'CORUM') in terms
    assert ('has_evidence_of_type', 'MI:0019') in terms
    assert ('publications', 'PMID:10537314') in terms
    members = [e for e in out.entities.values() if e.entity_type == 'protein']
    assert sorted(e.identifier for e in members) == ['3685', '3693', '5329']
    assert all(e.namespace == Namespace.ENTREZ and e.taxon == '9606' for e in members)
    assert len(out.relations) == 3 and {r.predicate for r in out.relations} == {'has_member'}


def test_predicted_complex_has_no_name_or_literature_evidence():
    _, out = extract(parsed()['HC4831'], 'predicted_complexes')
    (cplx,) = complexes(out)
    assert [i['ns'] for i in cplx.identifiers] == [COMPLEAT_NAMESPACE]
    assert cplx.annotations == []
    assert len(out.relations) == 15


def test_database_that_is_not_a_database_gives_no_data_source():
    for row in (parsed()['HC6106'], {**parsed()['HC1332'], 'database': 'Additional publication'}):
        _, out = extract(row)
        (cplx,) = complexes(out)
        assert not [a for a in cplx.annotations if a['term'] == 'supporting_data_source']
    # a cluster of high-throughput data, which the source lists as literature: the name and
    # the database are as the source gives them
    _, out = extract(parsed()['HC1332'])
    (cplx,) = complexes(out)
    assert cplx.label == 'HT_SC_Cluster37'
    assert ('supporting_data_source', 'CYC2008') in {(a['term'], a['value']) for a in cplx.annotations}
    assert DATABASES['Additional publication'] is None


def test_malformed_row_keeps_the_complex_without_an_invented_reference():
    _, out = extract(parsed()['HC3339'])
    (cplx,) = complexes(out)
    assert [a for a in cplx.annotations if a['term'] in ('publications', 'has_evidence_of_type')] == []


def test_equal_members_stay_two_complexes_with_their_accessions():
    first = dict(parsed()['HC2543'])
    second = {**first, 'compleat_id': 'HC99999', 'name': 'another'}
    out = SilverExtractor('compleat', 'literature_complexes')
    for i, row in enumerate((first, second)):
        out.process_record(map_complex(row), row, f'fixture:{i}', i)
    assert sorted(e.identifier for e in complexes(out)) == ['HC2543', 'HC99999']
    assert len(out.relations) == 6


def test_unknown_database_is_an_error_and_there_are_no_pairs():
    row = {**parsed()['HC1'], 'database': 'MIPS'}
    with pytest.raises(ValueError, match='Unknown COMPLEAT database'):
        map_complex(row)
    _, out = extract(parsed()['HC4831'])
    assert {r.predicate for r in out.relations} == {'has_member'}
    assert map_complex({**parsed()['HC1'], 'members': []}) is None
