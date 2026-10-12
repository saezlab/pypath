"""Source-specific Biolink contract for the ORegAnno input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import oreganno as mod
from pypath.inputs_v2.parsers import oreganno as parser

HEADER = (
    'ORegAnno_ID\tSpecies\tOutcome\tType\tGene_Symbol\tGene_ID\tGene_Source\t'
    'Regulatory_Element_Symbol\tRegulatory_Element_ID\tRegulatory_Element_Source\t'
    'dbSNP_ID\tPMID\tDataset\tBuild\tStrand\tChr\tStart\tEnd\n'
)
TFBS = 'TRANSCRIPTION FACTOR BINDING SITE'


def line(oreg, species, outcome, type_, gene, gene_id, gene_src, tf, tf_id, tf_src,
         pmid, dataset, build='hg19'):
    return '\t'.join(
        [oreg, species, outcome, type_, gene, gene_id, gene_src, tf, tf_id, tf_src,
         'N/A', pmid, dataset, build, '1', 'chr1', '10', '20']
    ) + '\n'


# Lines copied from ORegAnno_Combined_2016.01.19.tsv (the three genome builds of
# the first record are on three lines).
POS = 'POSITIVE OUTCOME'
LINES = [
    line('OREG0000010', 'Homo sapiens', POS, TFBS, 'MACROD1', 'ENSG00000133315', 'ENSEMBL', 'ESR1', 'ENSG00000091831', 'ENSEMBL', '15691879', 'N/A', 'hg18'),
    line('OREG0000010', 'Homo sapiens', POS, TFBS, 'MACROD1', 'ENSG00000133315', 'ENSEMBL', 'ESR1', 'ENSG00000091831', 'ENSEMBL', '15691879', 'N/A', 'hg19'),
    line('OREG0000010', 'Homo sapiens', POS, TFBS, 'MACROD1', 'ENSG00000133315', 'ENSEMBL', 'ESR1', 'ENSG00000091831', 'ENSEMBL', '15691879', 'N/A', 'hg38'),
    line('OREG0004530', 'Homo sapiens', POS, TFBS, 'CACNA2D2', '9254', 'NCBI', 'REST', '5978', 'NCBI', 'N/A', 'NRSF/REST ChIPSeq sites'),
    line('OREG0054735', 'Escherichia coli K-12', POS, TFBS, 'alkA', 'NA', 'USER DEFINED', 'ada', 'ECK120012595', 'RegulonDB', '23203884', 'RegulonDB'),
    line('OREG0002907', 'Homo sapiens', POS, TFBS, 'TFF1', '7031', 'NCBI', 'FOXA1/FOXA2', 'FOXA1/FOXA2', 'N/A', 'N/A', 'N/A'),
    line('OREG0001817', 'Caenorhabditis elegans', POS, TFBS, 'ges-1', 'WBGene00001578', 'ENSEMBL', 'elt-2', '181250', 'NCBI', ' 7649372', 'N/A'),
    # Left out: no gene, negative outcome, no factor, unspecified factor,
    # a microRNA, another record type.
    line('OREG0003048', 'Homo sapiens', POS, TFBS, 'N/A', 'N/A', 'N/A', 'STAT1', '6772', 'NCBI', 'N/A', 'STAT1 experiment-derived sites'),
    line('OREG0000011', 'Homo sapiens', 'NEGATIVE OUTCOME', TFBS, 'NR3C1', 'ENSG00000113580', 'ENSEMBL', 'N/A', 'N/A', 'N/A', '11996936', 'N/A'),
    line('OREG0000012', 'Homo sapiens', POS, TFBS, 'NR3C1', 'ENSG00000113580', 'ENSEMBL', 'N/A', 'N/A', 'N/A', '11996936', 'N/A'),
    line('OREG0000448', 'Drosophila melanogaster', POS, TFBS, 'Sgs4', 'FBgn0003374', 'ENSEMBL', 'Unspecified', 'Unspecified', 'N/A', '9819059', 'FLYREG'),
    line('OREG0069321', 'Homo sapiens', POS, TFBS, 'POT1', 'ENSG00000128513', 'ENSEMBL', 'hsa-mir-141', 'MI0000457', 'miRBase', '23874542', 'N/A'),
    line('OREG0000013', 'Homo sapiens', POS, 'REGULATORY REGION', 'CYP2A6', '1548', 'NCBI', 'N/A', 'N/A', 'N/A', '14974084', 'N/A'),
]
TEXT = HEADER + ''.join(LINES)
KEPT = ['OREG0000010', 'OREG0004530', 'OREG0054735', 'OREG0002907', 'OREG0001817']


def read(text=TEXT, function=None):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list((function or mod.resource.interactions._raw_parser)(opener))


def extract(row):
    record = mod.resource.interactions.mapper(row)
    assert record is not None
    out = SilverExtractor('oreganno', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def ids(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def test_parser_keeps_the_binding_sites_of_a_named_factor_and_gene_once():
    rows = read()
    assert [r['oreganno_id'] for r in rows] == KEPT
    assert rows[0] == {
        'oreganno_id': 'OREG0000010',
        'species': 'Homo sapiens',
        'tf_symbol': 'ESR1',
        'tf_source': 'ENSEMBL',
        'tf_id': 'ENSG00000091831',
        'gene_symbol': 'MACROD1',
        'gene_source': 'ENSEMBL',
        'gene_id': 'ENSG00000133315',
        'pmid': '15691879',
        'dataset': None,
    }


def test_every_other_line_has_one_reason_to_be_left_out():
    reasons = [parser.exclusion_reason(r) for r in read(function=parser.iter_lines)]
    assert reasons.count(None) == 7
    assert reasons[7:] == [
        'no target gene',
        'outcome is not positive',
        'no transcription factor',
        'transcription factor unspecified',
        'regulatory element is a microRNA',
        'not a transcription factor binding site',
    ]


def test_an_identifier_without_a_source_is_not_an_accession():
    row = read()[3]
    assert (row['tf_symbol'], row['tf_id'], row['tf_source']) == (
        'FOXA1/FOXA2',
        None,
        None,
    )


def test_parser_rejects_an_empty_file_and_a_missing_column():
    with pytest.raises(ValueError, match='no rows'):
        read(HEADER)
    with pytest.raises(ValueError, match='Outcome'):
        read(HEADER.replace('Outcome', 'Result') + LINES[0])


def test_parser_rejects_a_short_line_except_at_the_end():
    short = 'OREG0000011\tHomo sapiens\tN/A\n'
    with pytest.raises(ValueError, match='number of fields'):
        read(HEADER + short + LINES[0])
    # The rescued file ends in an incomplete line.
    assert [r['oreganno_id'] for r in read(TEXT + short)] == KEPT


def test_binding_site_is_interacts_with_and_has_no_effect():
    _, out = extract(read()[0])
    relation = out.relations[0]
    assert relation.predicate == 'interacts_with'
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')
    assert annotations(relation, 'publications') == ['PMID:15691879']
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert subject.entity_type == 'protein'
    assert ids(subject) == {
        (Namespace.GENESYMBOL, 'ESR1'),
        (Namespace.ENSG, 'ENSG00000091831'),
    }
    assert target.entity_type == 'gene'
    assert ids(target) == {
        (Namespace.GENESYMBOL, 'MACROD1'),
        (Namespace.ENSG, 'ENSG00000133315'),
    }
    assert subject.taxon == target.taxon == '9606'


def test_ncbi_identifier_is_an_entrez_gene_id_and_dataset_is_a_data_source():
    _, out = extract(read()[1])
    relation = out.relations[0]
    subject = out.entities[relation.subject_entity_key]
    assert (Namespace.ENTREZ, '5978') in ids(subject)
    assert annotations(relation, 'supporting_data_source') == [
        'NRSF/REST ChIPSeq sites'
    ]
    assert not annotations(relation, 'publications')


def test_species_other_than_human_get_their_own_taxon():
    _, out = extract(read()[2])
    relation = out.relations[0]
    subject = out.entities[relation.subject_entity_key]
    assert subject.taxon == '83333'
    # The user defined gene identifier of RegulonDB records is not an accession.
    assert {n for n, _ in ids(out.entities[relation.object_entity_key])} == {
        Namespace.GENESYMBOL
    }
    assert (Namespace.GENESYMBOL, 'ada') in ids(subject)
    assert ('regulondb', 'ECK120012595') in {
        (str(getattr(n, 'value', n)), i) for n, i in ids(subject)
    }


@pytest.mark.parametrize('number', ['0123456', '123456789'])
def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(number):
    _, out = extract(dict(read()[4], pmid=number))
    assert annotations(out.relations[0], 'publications') == []


def test_publication_with_a_leading_blank_is_read():
    _, out = extract(read()[4])
    assert annotations(out.relations[0], 'publications') == ['PMID:7649372']


def test_unknown_species_and_dataset_are_rejected():
    row = read()[0]
    with pytest.raises(ValueError, match='Unknown ORegAnno species'):
        mod.resource.interactions.mapper(dict(row, species='Mus spretus'))
    with pytest.raises(ValueError, match='Unknown ORegAnno dataset'):
        mod.resource.interactions.mapper(dict(row, dataset='Other'))


@pytest.mark.parametrize('missing', ['tf_symbol', 'gene_symbol'])
def test_row_without_a_participant_emits_no_relation(missing):
    role = missing.split('_')[0]
    row = dict(read()[0], **{missing: None, f'{role}_id': None})
    assert mod.resource.interactions.mapper(row) is None
