"""HuRI and the other CCSB datasets: reader, mapper and extractor on source lines."""

from collections import Counter
import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import huri
from pypath.inputs_v2.parsers.huri import iter_screen

DB = 'comment:vector name: pDEST-DB|comment:centromeric vector|comment:yeast strain: '
AD = 'comment:vector name: pDEST-AD|comment:centromeric vector|comment:yeast strain: '


def mitab(
    id_a, id_b, alt_a, alt_b, orf_a, orf_b, method, author, publication,
    confidence='-', roles=True, vector_a='-', vector_b='-', screens='1',
):
    """A line of the PSI-MI TAB files; only these columns vary in the source."""
    return '\t'.join([
        id_a, id_b, alt_a, alt_b,
        f'human orfeome collection:{orf_a}(author assigned name)',
        f'human orfeome collection:{orf_b}(author assigned name)',
        method, author, publication,
        'taxid:9606(Homo Sapiens)', 'taxid:9606(Homo Sapiens)',
        'psi-mi:MI:0915(physical association)', '-', '-', confidence, '-', '-', '-',
        'psi-mi:MI:0496(bait)' if roles else '-',
        'psi-mi:MI:0496(prey)' if roles else '-',
        'psi-mi:MI:0326(protein)', 'psi-mi:MI:0326(protein)', '-', '-', '-',
        vector_a, vector_b, f'comment:Found in screens {screens}.',
        'taxid:4932(Saccharomyces cerevisiae)', '-', '2019/10/16',
        '-', '-', '-', '-', '-',
        'gal4 dna binding domain:n-n (DB domain (n-terminal))',
        'gal4 activation domain:n-n (AD domain (n-terminal))', '-', '-',
        'psi-mi:MI1180(partial DNA sequence identification)',
        'psi-mi:MI1180(partial DNA sequence identification)',
    ]) + '\n'


# HuRI.psi line 1: a primary screen hit with an author score and no PubMed ID.
HURI_LINE = mitab(
    'uniprotkb:P0DP25', 'uniprotkb:A0A087WXN0',
    'ensembl:ENST00000291295.13|ensembl:ENSP00000291295.8|ensembl:ENSG00000160014.16',
    'ensembl:ENST00000612316.4|ensembl:ENSP00000481156.1|ensembl:ENSG00000130518.16',
    1, 56859, 'psi-mi:MI:1112(two hybrid prey pooling approach)',
    'Luck et al.(2019)', 'unassigned1304',
    confidence='author score: 0.899116271005', vector_a=DB + 'Y8930',
    vector_b='comment:vector name: pQZ213|comment:2u vector|comment:yeast strain: Y8800',
    screens='4',
)
# Yu-11.psi line 2: a validation assay, without bait and prey roles.
YU_LINE = mitab(
    'uniprotkb:P62166-1', 'uniprotkb:Q86UW9-1',
    'ensembl:ENST00000372398.5|ensembl:ENSP00000361475.3|ensembl:ENSG00000107130.9',
    'ensembl:ENST00000324432.9|ensembl:ENSP00000322885.5|ensembl:ENSG00000091073.19',
    4378, 5426, 'psi-mi:MI:1356(validated two hybrid)',
    'Yu et al.(2011)', 'pubmed:21516116', roles=False,
)
# H-I-05.psi line 1: the prey ORF has no UniProt accession.
RUAL_LINE = mitab(
    'uniprotkb:Q96HB5-1', 'ensembl:ENSP00000493357.1',
    'ensembl:ENST00000606812.5|ensembl:ENSP00000475676.1|ensembl:ENSG00000147144.12',
    'ensembl:ENST00000641098.1|ensembl:ENSP00000493357.1|ensembl:ENSG00000105443.14',
    4535, 3723, 'psi-mi:MI:0397(two hybrid array)',
    'Rual et al.(2005)', 'pubmed:16189514',
    vector_a=DB + 'MaV203', vector_b=AD + 'MaV103',
)
# H-I-05.psi: SMN1 and SMN2 encode the same protein, which pairs with itself.
SMN = (
    'ensembl:ENST00000380707.8|ensembl:ENSP00000370083.4|ensembl:ENSG00000172062.16|'
    'ensembl:ENST00000380743.8|ensembl:ENSP00000370119.4|ensembl:ENSG00000205571.13'
)
SMN_LINE = mitab(
    'uniprotkb:Q16637-1', 'uniprotkb:Q16637-1', SMN, SMN, 1063, 1063,
    'psi-mi:MI:1112(two hybrid prey pooling approach)',
    'Rual et al.(2005)', 'pubmed:16189514',
    vector_a=DB + 'MaV203', vector_b=AD + 'MaV103',
)
# Lit-BM.tsv lines 1 and 3; the second is a self interaction.
PAIRS = 'ENSG00000001167\tENSG00000066136\nENSG00000002016\tENSG00000002016\n'


def read(text, dataset='huri'):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(huri.resource.datasets()[dataset]._raw_parser(opener))


def extract(rows, dataset='huri'):
    out = SilverExtractor('huri', dataset)
    for index, row in enumerate(rows):
        record = huri.resource.datasets()[dataset].mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def identifiers(entity):
    return [(i.type, i.value) for i in entity.identifiers]


def values(relation, term, scope='relation'):
    return [
        a['value'] for a in relation.annotations
        if a['term'] == term and a['scope'] == scope
    ]


def test_module_import_defines_the_datasets_without_downloading():
    assert list(huri.resource.datasets()) == [
        'huri', 'hi_ii_14', 'yang_2016', 'yu_2011', 'venkatesan_2009',
        'hi_i_05', 'test_space_screens_2019', 'lit_bm',
    ]
    assert huri.config.id.name == 'HURI'


def test_reader_names_the_columns_of_the_headerless_files():
    rows = read(HURI_LINE + YU_LINE)
    assert len(rows) == 2
    assert list(rows[0]) == list(huri._MITAB_COLUMNS)
    assert rows[0]['id_a'] == 'uniprotkb:P0DP25'
    assert rows[0]['aliases_b'] == 'human orfeome collection:56859(author assigned name)'
    assert rows[0]['annotations_interaction'] == 'comment:Found in screens 4.'
    assert rows[1]['publication'] == 'pubmed:21516116'
    assert read(PAIRS, 'lit_bm') == [
        {'ensembl_gene_a': 'ENSG00000001167', 'ensembl_gene_b': 'ENSG00000066136'},
        {'ensembl_gene_a': 'ENSG00000002016', 'ensembl_gene_b': 'ENSG00000002016'},
    ]


def test_reader_rejects_a_file_with_another_layout():
    with pytest.raises(ValueError, match='Line 1 has 1 fields, expected 42'):
        read('<!DOCTYPE html>\n<html>\n')
    with pytest.raises(ValueError, match='Line 1 has 42 fields, expected 2'):
        read(HURI_LINE, 'lit_bm')
    assert read('') == []


def test_screen_row_keeps_source_identifiers_method_and_score():
    row = read(HURI_LINE)[0]
    relation = huri.resource.huri.mapper(row)
    # The row states "physical association" (MI:0915): not a direct interaction.
    assert slot_name(relation.predicate) == 'interacts_with'
    assert identifiers(relation.subject) == [
        (Namespace.UNIPROT, 'P0DP25'),
        (Namespace.ENST, 'ENST00000291295.13'),
        (Namespace.ENSP, 'ENSP00000291295.8'),
        (Namespace.ENSG, 'ENSG00000160014.16'),
    ]
    assert identifiers(relation.object)[0] == (Namespace.UNIPROT, 'A0A087WXN0')

    out = extract([row])
    (observed,) = out.relations
    assert observed.predicate == 'interacts_with'
    assert values(observed, 'has_evidence_of_type') == ['MI:1112']
    # The file's placeholder stands for the paper the portal names for HuRI.
    assert values(observed, 'publications') == ['PMID:32296183']
    (score,) = [
        a['quantity'] for a in observed.annotations
        if a['term'] == 'has_confidence_score'
    ]
    assert score['has_numeric_value'] == pytest.approx(0.899116271005)
    assert score['source_field'] == 'author score'
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}
    # The file states the interaction type of the row (physical association),
    # which is also the source predicate kept on the relation, and the screens
    # of the pair.
    assert values(observed, 'original_predicate') == ['MI:0915']
    assert values(observed, 'huri:screen') == ['4']
    # The yeast host is not the taxon of anything.
    assert not values(observed, 'in_taxon')


def test_prey_label_wins_over_the_bait_accession_the_files_repeat():
    (observed,) = extract(read(HURI_LINE)).relations
    assert values(observed, 'psi_mi:experimental_role', 'subject') == ['MI:0496']
    assert values(observed, 'psi_mi:experimental_role', 'object') == ['MI:0498']


@pytest.mark.parametrize('reference', ['pubmed:0123456', 'pubmed:123456789'])
def test_pubmed_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(
    reference,
):
    (row,) = read(YU_LINE, 'yu_2011')
    (observed,) = extract([dict(row, publication=reference)], 'yu_2011').relations
    assert values(observed, 'publications') == []


def test_validation_row_has_pubmed_and_no_invented_role_or_score():
    (observed,) = extract(read(YU_LINE), 'yu_2011').relations
    assert values(observed, 'has_evidence_of_type') == ['MI:1356']
    assert values(observed, 'publications') == ['PMID:21516116']
    terms = {a['term'] for a in observed.annotations}
    assert 'psi_mi:experimental_role' not in terms
    assert 'has_confidence_score' not in terms
    # Every row states the interaction type, which is not a direct interaction.
    assert observed.predicate == 'interacts_with'
    assert values(observed, 'original_predicate') == ['MI:0915']
    assert values(observed, 'huri:screen') == ['1']


def test_isoform_and_ensembl_protein_stay_as_the_source_gives_them():
    relation = huri.resource.hi_i_05.mapper(read(RUAL_LINE)[0])
    assert identifiers(relation.subject)[0] == (Namespace.UNIPROT, 'Q96HB5-1')
    assert relation.subject.molecular_form['isoform_identifier'] == {
        'ns': 'uniprot', 'id': 'Q96HB5-1',
    }
    # No UniProt accession is looked up for an ORF the source gives as ENSP.
    assert identifiers(relation.object) == [
        (Namespace.ENSP, 'ENSP00000493357.1'),
        (Namespace.ENST, 'ENST00000641098.1'),
        (Namespace.ENSG, 'ENSG00000105443.14'),
    ]


def test_orf_matching_two_genes_keeps_both_and_self_interaction_is_kept():
    relation = huri.resource.hi_i_05.mapper(read(SMN_LINE)[0])
    genes = [v for ns, v in identifiers(relation.subject) if ns == Namespace.ENSG]
    assert genes == ['ENSG00000172062.16', 'ENSG00000205571.13']
    out = extract(read(SMN_LINE), 'hi_i_05')
    assert len(out.entities) == 1
    (observed,) = out.relations
    assert observed.subject_entity_key == observed.object_entity_key


def test_same_pair_in_three_assays_is_three_observations_of_one_statement():
    lines = [
        HURI_LINE,
        HURI_LINE.replace('MI:1112(two hybrid prey pooling approach)', 'MI:0397(two hybrid array)'),
    ]
    out = extract(read(''.join(lines)))
    assert len(out.relations) == 2
    assert len({r.relation_key for r in out.relations}) == 1
    assert [values(r, 'has_evidence_of_type') for r in out.relations] == [
        ['MI:1112'], ['MI:0397'],
    ]


def test_rows_that_must_not_become_interactions():
    row = read(HURI_LINE)[0]
    mapper = huri.resource.huri.mapper
    # A negative result, in any spelling the format does not define as positive.
    assert mapper(dict(row, negative='true')) is None
    assert mapper(dict(row, negative='yes')) is None
    assert mapper(dict(row, negative='false')) is not None
    # A participant that is not a protein, or has no usable identifier.
    assert mapper(dict(row, interactor_type_b='psi-mi:MI:0328(small molecule)')) is None
    assert mapper(dict(row, id_b='-')) is None
    assert mapper(dict(row, id_b='human orfeome collection:56859')) is None


def test_unknown_values_are_left_out_not_defaulted():
    row = dict(
        read(HURI_LINE)[0],
        publication='unassigned9999',
        experimental_role_b='psi-mi:MI:0499(unspecified role)',
        detection_method='two hybrid',
        confidence='intact-miscore:0.5',
    )
    (observed,) = extract([row]).relations
    assert not values(observed, 'publications')
    assert not values(observed, 'has_evidence_of_type')
    assert not values(observed, 'psi_mi:experimental_role', 'object')
    assert 'has_confidence_score' not in {a['term'] for a in observed.annotations}


def test_gene_pairs_are_ensembl_genes_and_only_the_screens_have_a_publication():
    rows = read(PAIRS, 'lit_bm')
    literature = extract(rows, 'lit_bm')
    # Lit-BM.tsv states no interaction type for a pair; the source's definition
    # of the set (a binary assay among the evidence) is not a statement of the row.
    assert [r.predicate for r in literature.relations] == ['interacts_with'] * 2
    # Lit-BM.tsv names no method and no publication for a pair.
    assert {a['term'] for r in literature.relations for a in r.annotations} == {
        'in_taxon',
    }
    assert {
        (i['ns'], i['id']) for e in literature.entities.values() for i in e.identifiers
    } == {
        (Namespace.ENSG, 'ENSG00000001167'),
        (Namespace.ENSG, 'ENSG00000002016'),
        (Namespace.ENSG, 'ENSG00000066136'),
    }

    screens = extract(rows, 'test_space_screens_2019')
    assert [r.predicate for r in screens.relations] == ['interacts_with'] * 2
    assert values(screens.relations[0], 'publications') == ['PMID:32296183']


def test_gene_pair_with_a_malformed_identifier_is_not_emitted():
    mapper = huri.resource.lit_bm.mapper
    assert mapper({'ensembl_gene_a': 'ENSG00000001167', 'ensembl_gene_b': ''}) is None
    assert mapper({'ensembl_gene_a': 'ENSG00000001167', 'ensembl_gene_b': 'TP53'}) is None


# HuRI.psi line 58: a pair of the test space screens, listed in HuRI.psi but not
# part of the published map.
TEST_SPACE_LINE = mitab(
    'ensembl:ENSP00000476265.1', 'uniprotkb:Q8TC71-1',
    'ensembl:ENST00000609374.5|ensembl:ENSP00000476265.1|ensembl:ENSG00000167110.17',
    'ensembl:ENST00000295213.8|ensembl:ENSP00000295213.4|ensembl:ENSG00000163071.10',
    14862, 1340, 'psi-mi:MI:1356(validated two hybrid)',
    'Luck et al.(2019)', 'unassigned1304', roles=False, screens='8,3,12',
).replace('Found in screens', 'Found in test space screens')


def test_screens_of_a_pair_are_attribute_values_as_the_file_numbers_them():
    rows = read(mitab(
        'uniprotkb:Q16637-1', 'uniprotkb:Q16637-1', SMN, SMN, 1063, 1063,
        'psi-mi:MI:1112(two hybrid prey pooling approach)',
        'Rual et al.(2005)', 'pubmed:16189514', screens='1,4,7',
    ) + mitab(
        'uniprotkb:Q16637-1', 'uniprotkb:Q16637-1', SMN, SMN, 1063, 1063,
        'psi-mi:MI:1112(two hybrid prey pooling approach)',
        'Rual et al.(2005)', 'pubmed:16189514', screens='',
    ))
    first, second = extract(rows).relations
    assert values(first, 'huri:screen') == ['1', '4', '7']
    # "Found in screens ." (HI-II-14.psi has 1,744 such rows) names none.
    assert values(second, 'huri:screen') == []


def test_pairs_of_the_test_space_screens_are_not_in_the_huri_dataset():
    rows = read(HURI_LINE + TEST_SPACE_LINE + HURI_LINE, 'hi_ii_14')
    assert len(rows) == 3
    # The row is still read by the parser of the other screens ...
    assert huri.parse_screens(rows[1]['annotations_interaction']) == (['8', '3', '12'], True)
    assert huri.parse_screens(rows[0]['annotations_interaction']) == (['4'], False)
    # ... and left out of the HuRI dataset, with the exclusion counted.
    stats = Counter()
    kept = list(iter_screen(
        SimpleNamespace(result=io.StringIO(HURI_LINE + TEST_SPACE_LINE + HURI_LINE)),
        exclude_test_space=True, stats=stats,
    ))
    assert [row['id_a'] for row in kept] == ['uniprotkb:P0DP25'] * 2
    assert dict(stats) == {'rows': 3, 'test_space_screens_excluded': 1}
    opener = SimpleNamespace(result=io.StringIO(TEST_SPACE_LINE))
    assert list(huri.resource.huri._raw_parser(opener)) == []
    opener = SimpleNamespace(result=io.StringIO(TEST_SPACE_LINE))
    assert len(list(huri.resource.hi_ii_14._raw_parser(opener))) == 1


def test_unknown_interaction_annotation_is_an_error():
    for cell in ['-', 'comment:Found in screen 4.', 'comment:Found in test screens 4.']:
        with pytest.raises(ValueError, match='Unexpected interaction annotation'):
            huri.parse_screens(cell)
    line = HURI_LINE.replace('comment:Found in screens 4.', 'comment:another comment')
    with pytest.raises(ValueError, match='Unexpected interaction annotation'):
        list(huri.resource.huri._raw_parser(SimpleNamespace(result=io.StringIO(line))))


def test_only_a_row_that_states_a_direct_interaction_type_is_direct():
    row = read(HURI_LINE)[0]
    mapper = huri.resource.huri.mapper
    # "physical association" and "association" say that the proteins interact.
    for cell in ['psi-mi:MI:0915(physical association)', 'psi-mi:MI:0914(association)']:
        assert slot_name(mapper(dict(row, interaction_type=cell)).predicate) == 'interacts_with'
    # "direct interaction" (MI:0407) is a statement of the row.
    direct = mapper(dict(row, interaction_type='psi-mi:MI:0407(direct interaction)'))
    assert slot_name(direct.predicate) == 'directly_physically_interacts_with'
    (observed,) = extract([dict(row, interaction_type='psi-mi:MI:0407(direct interaction)')]).relations
    assert values(observed, 'original_predicate') == ['MI:0407']


@pytest.mark.parametrize('cell', ['-', '', 'psi-mi:MI:0403(colocalization)', 'physical association'])
def test_unknown_interaction_type_is_an_error_not_a_default(cell):
    row = dict(read(HURI_LINE)[0], interaction_type=cell)
    with pytest.raises(ValueError, match='Unknown PSI-MI interaction type'):
        huri.resource.huri.mapper(row)
