"""Source contract of the CollecTRI2 input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.collectri2 import resource
from pypath.inputs_v2.parsers.collectri2 import iter_pairs

COMMENT = '#: :filename=/home/mvazque2/.rbbt/var/jobs/ExTRI/pairs/Default.tsv#:namespace=Hsa/oct2016#:type=:double'
# Header of CollecTRI2.tsv.gz.
HEADER = '\t'.join(
    [
        '#TF:TG',
        'Transcription Factor (Associated Gene Name)',
        'Target Gene (Associated Gene Name)',
        '[ExTRI2] present', '[ExTRI2] PMID', '[ExTRI2] Sign', '[ExTRI2] Transcription Factor Type',
        '[HTRI] present', '[HTRI] Technique', '[HTRI] PMID', '[HTRI] Confidence',
        '[TRRUST] present', '[TRRUST] Regulation', '[TRRUST] PMID',
        '[TFactS] present', '[TFactS] Sign', '[TFactS] Species', '[TFactS] Source', '[TFactS] PMID', '[TFactS] Confidence',
        '[GOA] present', '[GOA] Sign', '[GOA] PMID',
        '[IntAct] present', '[IntAct] PMID', '[IntAct] Method ID',
        '[SIGNOR] present', '[SIGNOR] Effect', '[SIGNOR] Sign', '[SIGNOR] PMID',
        '[CytReg] present', '[CytReg] Assay type', '[CytReg] species', '[CytReg] Activation/Repression', '[CytReg] PMID', '[CytReg] Year of publication',
        '[GEREDB] present', '[GEREDB] Effect', '[GEREDB] PMID',
        '[NTNU Curated] present', '[NTNU Curated] Sign', '[NTNU Curated] PMID',
        '[Pavlidis2021] present', '[Pavlidis2021] PMID', '[Pavlidis2021] MoR',
        '[DoRothEA_A] present', '[DoRothEA_A] PMID', '[DoRothEA_A] Directed', '[DoRothEA_A] Effect',
        'Lambert', 'Lovering', 'GO:0003700', 'GO:0140223', 'GO:0003712', 'TFClass', 'Auto-regulation',
    ]
)
COLUMNS = HEADER.split('\t')
# The non-empty cells of lines of CollecTRI2.tsv.gz, with the weight of the
# pair in CollecTRI2_regulome.tsv.gz.
PAIRS = {
    # A sign of one of three items, in two sources.
    'majority up': (
        {'[ExTRI2] present': 'ExTRI2', '[ExTRI2] PMID': '22206000|22206000|22206000', '[ExTRI2] Sign': 'UNKNOWN|UP|UNKNOWN',
         '[ExTRI2] Transcription Factor Type': 'dbTF', '[NTNU Curated] present': 'NTNU Curated',
         '[NTNU Curated] Sign': '+|', '[NTNU Curated] PMID': '22206000|22206000'},
        ('MEOX1', 'CDKN1A', '1'),
    ),
    # Two PubMed IDs say down, the third says nothing.
    'majority down': (
        {'[ExTRI2] present': 'ExTRI2', '[ExTRI2] PMID': '23451207|31602225|31538102|31602225',
         '[ExTRI2] Sign': 'UNKNOWN|DOWN|DOWN|DOWN', '[ExTRI2] Transcription Factor Type': 'dbTF',
         '[NTNU Curated] present': 'NTNU Curated', '[NTNU Curated] PMID': '23451207'},
        ('KLF4', 'CDH2', '-1'),
    ),
    # One article says up and another down; the weight is -1.
    'tie': (
        {'[ExTRI2] present': 'ExTRI2', '[ExTRI2] PMID': '28649987|28649987|19801450', '[ExTRI2] Sign': 'UP|UP|DOWN',
         '[ExTRI2] Transcription Factor Type': 'dbTF', '[GEREDB] present': 'GEREDB', '[GEREDB] Effect': 'unknown',
         '[GEREDB] PMID': '18332871'},
        ('PPARG', 'PDGFRB', '-1'),
    ),
    # No sign anywhere; the weight is -1 for the pair.
    'unsigned negative': (
        {'[ExTRI2] present': 'ExTRI2', '[ExTRI2] PMID': '21408158', '[ExTRI2] Sign': 'UNKNOWN',
         '[ExTRI2] Transcription Factor Type': 'dbTF', '[NTNU Curated] present': 'NTNU Curated',
         '[NTNU Curated] PMID': '21408158'},
        ('NR1D1', 'NUPR1', '-1'),
    ),
    # No sign anywhere; the weight is the default.
    'unsigned default': (
        {'[ExTRI2] present': 'ExTRI2', '[ExTRI2] PMID': '35388756', '[ExTRI2] Sign': 'UNKNOWN',
         '[ExTRI2] Transcription Factor Type': 'coTF candidate'},
        ('CHCHD2', 'GNPTG', '1'),
    ),
    # The transcription factor is a dimer.
    'complex': (
        {'[ExTRI2] present': 'ExTRI2', '[ExTRI2] PMID': '23974100|23974100|23974100|23974100',
         '[ExTRI2] Sign': 'UNKNOWN|UNKNOWN|UNKNOWN|UNKNOWN', '[ExTRI2] Transcription Factor Type': 'dbTF'},
        ('NFKB', 'PLK4', '1'),
    ),
}


def main_text(pairs=PAIRS):
    lines = [COMMENT, HEADER]
    for cells, (tf, target, _) in pairs.values():
        row = {**cells, '#TF:TG': f'{tf}:{target}', COLUMNS[1]: tf, COLUMNS[2]: target}
        lines.append('\t'.join(row.get(column, '') for column in COLUMNS))
    return '\n'.join(lines) + '\n'


def regulome_text(pairs=PAIRS):
    lines = [
        '#: :filename=/home/mvazque2/.rbbt/var/jobs/ExTRI2/regulome/Default.tsv#:type=:list',
        '#ID\tsource\ttarget\tweight\tTF Type\tAuto-regulation',
    ]
    for number, (_, (tf, target, weight)) in enumerate(pairs.values(), 2):
        lines.append(f'{number}\t{tf}\t{target}\t{weight}\tdbTF\t')
    return '\n'.join(lines) + '\n'


class Download:
    def __init__(self, text):
        self.text = text

    def open(self, force_refresh=False):
        return SimpleNamespace(result=io.StringIO(self.text))


def parse(main=None, regulome=None):
    return list(
        iter_pairs(
            SimpleNamespace(result=io.StringIO(main or main_text())),
            regulome=Download(regulome or regulome_text()),
        )
    )


PARSED = dict(zip(PAIRS, parse()))


def extract(name):
    row = PARSED[name]
    record = resource.interactions.mapper(row)
    assert record is not None
    out = SilverExtractor('collectri2', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out.relations[0]


def terms(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_module_wires_the_two_files():
    assert list(resource.datasets()) == ['interactions']


def test_reader_counts_the_distinct_evidence_of_each_sign():
    assert (PARSED['majority up']['up'], PARSED['majority up']['down']) == (1, 0)
    # The same article twice says down once.
    assert (PARSED['majority down']['up'], PARSED['majority down']['down']) == (0, 2)
    assert (PARSED['tie']['up'], PARSED['tie']['down']) == (1, 1)
    assert (PARSED['unsigned default']['up'], PARSED['unsigned default']['down']) == (0, 0)


def test_reader_says_how_the_weight_relates_to_the_evidence():
    assert {name: row['sign_decision'] for name, row in PARSED.items()} == {
        'majority up': 'pmid_majority',
        'majority down': 'pmid_majority',
        'tie': 'pmid_tie',
        'unsigned negative': 'no_signed_evidence',
        'unsigned default': 'no_signed_evidence',
        'complex': 'no_signed_evidence',
    }
    contradicted = {**PAIRS, 'majority up': (PAIRS['majority up'][0], ('MEOX1', 'CDKN1A', '-1'))}
    assert parse(main_text(contradicted), regulome_text(contradicted))[0]['sign_decision'] == 'contradicted'


def test_reader_lists_the_sources_and_the_pubmed_ids_of_the_pair():
    row = PARSED['majority up']
    assert row['sources'] == ['ExTRI2', 'NTNU Curated']
    assert row['pmids'] == ['22206000']
    assert row['tf_type'] == 'dbTF'


def test_reader_rejects_a_defect_of_the_files():
    with pytest.raises(ValueError, match='comment'):
        parse(main_text().split('\n', 1)[1])
    with pytest.raises(ValueError, match='lacks expected columns'):
        parse(main_text().replace('[SIGNOR] Sign', '[SIGNOR] Sig'))
    with pytest.raises(ValueError, match='Unknown CollecTRI2 sign'):
        parse(main_text().replace('UNKNOWN|UP|UNKNOWN', 'UNKNOWN|MAYBE|UNKNOWN'))
    with pytest.raises(ValueError, match='evidence items and'):
        parse(main_text().replace('UNKNOWN|UP|UNKNOWN', 'UNKNOWN|UP'))
    with pytest.raises(ValueError, match='lacks the pair'):
        parse(regulome=regulome_text({k: v for k, v in PAIRS.items() if k != 'tie'}))
    with pytest.raises(ValueError, match='twice'):
        parse(main_text() + main_text().split('\n')[2] + '\n')


def test_a_majority_of_signed_evidence_sets_the_expression_effect():
    _, up = extract('majority up')
    _, down = extract('majority down')
    assert up.predicate == 'affects'
    assert terms(up, 'object_direction_qualifier') == ['increased']
    assert terms(up, 'object_aspect_qualifier') == ['expression']
    assert terms(down, 'object_direction_qualifier') == ['decreased']


@pytest.mark.parametrize('name', ['tie', 'unsigned negative', 'unsigned default', 'complex'])
def test_a_tie_a_default_or_an_inferred_sign_gives_a_bare_effect(name):
    _, relation = extract(name)
    assert relation.predicate == 'affects'
    assert terms(relation, 'object_direction_qualifier') == []
    assert terms(relation, 'object_aspect_qualifier') == []


def test_a_weight_against_the_evidence_is_not_carried():
    row = {**PARSED['majority up'], 'sign_decision': 'contradicted'}
    out = SilverExtractor('collectri2', 'interactions')
    out.process_record(resource.interactions.mapper(row), row, 'fixture:1', 1)
    assert terms(out.relations[0], 'object_direction_qualifier') == []


def test_the_sources_tell_text_mining_from_curation():
    _, relation = extract('majority up')
    assert terms(relation, 'supporting_data_source') == ['ExTRI2', 'NTNU.Curated']
    _, relation = extract('tie')
    assert terms(relation, 'supporting_data_source') == ['ExTRI2', 'GEREDB']


def test_publications_are_the_pubmed_ids_of_the_evidence():
    _, relation = extract('majority down')
    assert terms(relation, 'publications') == ['PMID:23451207', 'PMID:31602225', 'PMID:31538102']


def test_a_dimer_is_a_complex_and_a_gene_symbol_a_protein():
    record, _ = extract('complex')
    assert {(i.type, i.value) for i in record.subject.identifiers} == {(Namespace.NAME, 'NFKB')}
    record, _ = extract('majority up')
    assert {(i.type, i.value) for i in record.subject.identifiers} == {(Namespace.GENESYMBOL, 'MEOX1')}
    assert {(i.type, i.value) for i in record.object.identifiers} == {(Namespace.GENESYMBOL, 'CDKN1A')}


def test_an_unknown_weight_decision_or_source_raises():
    with pytest.raises(ValueError, match='Unknown CollecTRI2 weight'):
        resource.interactions.mapper({**PARSED['majority up'], 'weight': '0'})
    with pytest.raises(ValueError, match='sign decision'):
        resource.interactions.mapper({**PARSED['majority up'], 'sign_decision': 'regulon'})
    with pytest.raises(ValueError, match='Unknown CollecTRI2 source'):
        resource.interactions.mapper({**PARSED['majority up'], 'sources': ['ExTRI']})
