"""Source contract of the ncRDeathDB non-coding RNA-gene input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.ncrdeathdb import resource
from pypath.inputs_v2.parsers.ncrdeathdb import (
    SOURCE_HEADER,
    iter_homology_predictions,
    iter_interactions,
)


def line(*values):
    assert len(values) == len(SOURCE_HEADER)
    return '\t'.join(values)


# Rows of allNcRNACelldeathData.txt. Long cells are shortened.
LNCRNA = line(
    '465', 'HOTAIR', 'http://lncrna.com/hotair/', '', 'lncRNA', 'PTEN', 'apoptosis',
    'down', 'Homo sapiens', 'laryngeal squamous cell cancer', '23141928', '9606',
    '5728', '10q23del|BZS|CWS1', 'HGNC:9588|MIM:601728|HPRD:03431', '10', '10q23.3',
    'phosphatase and tensin homolog', 'protein-coding',
    'phosphatase and tensin homolog', '',
    'SiRNA-mediated knockdown of HOTAIR led to reduced invasion and increased apoptosis.',
)
LNCRNA_WITHOUT_GENE = line(
    '49', 'BC200', 'http://lncrna.com/BC200/', '', 'lncRNA', '', 'necrosis', '',
    'Homo sapiens', '', '25645334', '9606', '', '', '', '', '', '', '', '', '',
    'TUG1, BC200 and MIR156HG are down regulated when necrosis is induced.',
)
SNORNA = line(
    '2438', 'E2', '93A1D1FEA2961365482571F8000697C3?OpenDocument', '', 'snoRNA',
    'CEBPA', 'apoptosis', 'up', 'Homo sapiens', 'myeloid precursor cells', '11753385',
    '9606', '1050', 'C/EBP-alpha|CEBP', 'HGNC:1833|MIM:116897|Ensembl:ENSG00000245848',
    '19', '19q13.1', 'CCAAT/enhancer binding protein (C/EBP), alpha', 'protein-coding',
    'CCAAT/enhancer binding protein (C/EBP), alpha', '',
    'Ectopic expression of hnRNP E2 in myeloid precursor cells down-regulates C/EBPalpha.',
)
MIRNA_LIST = line(
    '218', 'miR-34a', 'hsa-miR-34a-5p,hsa-miR-34a-3p', 'MIMAT0000255,MIMAT0004557',
    'miRNA', 'ATG9A', 'autophagy', 'down', 'Homo sapiens', '', '22902544', '9606',
    '79065', 'APG9L1|MGD3208|mATG9', 'HGNC:22408|MIM:612204|Ensembl:ENSG00000198925',
    '2', '2q35', 'autophagy related 9A', 'protein-coding', 'autophagy related 9A', '',
    'miR-34a has recently been identified as a regulator of autophagy.',
)
# The category cell has a trailing space.
MIRNA_ONE_PRODUCT = line(
    '309', 'miR-17-5p', 'hsa-miR-17-5p', 'MIMAT0000070', 'miRNA ', 'BECN1', 'autophagy',
    'down', 'Homo sapiens', 'lung cancer cells A549', '24755562', '9606', '8678',
    'ATG6|VPS30|beclin1', 'HGNC:1034|MIM:604378|Ensembl:ENSG00000126581', '17', '17q21',
    'beclin 1, autophagy related', 'protein-coding', 'beclin 1, autophagy related', '',
    'Overexpression of miR-17-5p reduced beclin1 expression.',
)
# A list in the family notation, which is no list of miRBase names.
MIRNA_FAMILY_NOTATION = line(
    '2810', 'gga-miR-17-5p', 'miR-17/17-5p/20ab/20b-5p/93/106ab/427/518a-3p/519d',
    'MIMAT0001114', 'miRNA', 'HIF1A', 'autophagy', 'down', 'Gallus gallus', '',
    'prediction', '9031', '374177', '-', 'CGNC:8999|Ensembl:ENSGALG00000011870', '5', '-',
    'hypoxia inducible factor 1, alpha subunit', 'protein-coding',
    'hypoxia inducible factor 1, alpha subunit', '',
    'gga-miR-17-5p predicted by miR-20a (Homo sapiens), and HIF1A predicted by HIF1A '
    '(Homo sapiens) according to Pubmed ID (25485521)',
)
UNKNOWN_CATEGORY = MIRNA_ONE_PRODUCT.replace('miRNA ', 'circRNA')


def table(*lines):
    return SimpleNamespace(result=io.StringIO('\r\n'.join(lines) + '\r\n'))


def parse(text, reader=iter_interactions):
    (row,) = list(reader(table('\t'.join(SOURCE_HEADER), text)))
    return row


def extract(text, dataset='interactions', reader=iter_interactions):
    row = parse(text, reader)
    record = resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor('ncrdeathdb', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def annotations(relation, scope='relation'):
    return sorted(
        (a['term'], a['value']) for a in relation.annotations if a['scope'] == scope
    )


def test_lncrna_is_a_noncoding_rna_named_by_the_symbol():
    relation, rna, gene = extract(LNCRNA)
    assert relation.predicate == 'interacts_with'
    assert rna.entity_type == 'noncoding_rna_product'
    assert identifiers(rna) == {(Namespace.NAME, 'HOTAIR')}
    assert rna.taxon == '9606'
    assert gene.entity_type == 'gene'
    assert (gene.namespace, gene.identifier) == ('entrez', '5728')
    assert (Namespace.HGNC, '9588') in identifiers(gene)
    assert all(ns != Namespace.ENSG for ns, _ in identifiers(gene))
    assert annotations(relation) == [
        ('description', 'SiRNA-mediated knockdown of HOTAIR led to reduced invasion and increased apoptosis.'),
        ('publications', 'PMID:23141928'),
    ]


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    for number in ('0123456', '123456789'):
        text = LNCRNA.replace('23141928', number)
        assert parse(text)['pmids'] == []
        relation, _, _ = extract(text)
        assert [a for a in annotations(relation) if a[0] == 'publications'] == []


def test_the_link_in_the_mature_id_cell_of_a_lncrna_is_not_an_identifier():
    _, rna, _ = extract(LNCRNA)
    assert all('lncrna.com' not in i['id'] for i in rna.identifiers)
    assert parse(LNCRNA).get('mature_names') is None


def test_snorna_is_a_noncoding_rna_product_too():
    _, rna, gene = extract(SNORNA)
    assert rna.entity_type == 'noncoding_rna_product'
    assert (gene.namespace, gene.identifier) == ('entrez', '1050')


def test_microrna_list_is_not_expanded_and_has_no_accession():
    _, rna, _ = extract(MIRNA_LIST)
    assert rna.entity_type == 'microrna'
    assert identifiers(rna) == {(Namespace.NAME, 'miR-34a')}


def test_microrna_with_one_product_is_the_mirbase_mature_microrna():
    relation, rna, _ = extract(MIRNA_ONE_PRODUCT)
    assert (rna.namespace, rna.identifier) == ('mirbase_mature', 'MIMAT0000070')
    assert identifiers(rna) == {
        (Namespace.MIRBASE_MATURE, 'MIMAT0000070'),
        (Namespace.NAME, 'hsa-miR-17-5p'),
        (Namespace.SYNONYM, 'miR-17-5p'),
    }


def test_pathway_tissue_and_expression_direction_are_not_on_the_relation():
    relation, _, _ = extract(MIRNA_ONE_PRODUCT)
    values = {v for _, v in annotations(relation)}
    assert 'autophagy' not in values
    assert 'down' not in values
    assert 'lung cancer cells A549' not in values
    assert relation.predicate == 'interacts_with'


def test_row_without_a_gene_gives_no_relation():
    row = parse(LNCRNA_WITHOUT_GENE)
    assert resource.datasets()['interactions'].mapper(row) is None


def test_homology_prediction_is_its_own_dataset_without_publication():
    with pytest.raises(ValueError, match='no curated observations'):
        parse(MIRNA_FAMILY_NOTATION)
    row = parse(MIRNA_FAMILY_NOTATION, iter_homology_predictions)
    assert row['prediction'] is True
    assert row['pmids'] == []
    relation, rna, gene = extract(
        MIRNA_FAMILY_NOTATION, 'homology_predictions', iter_homology_predictions
    )
    # The family notation and the accession next to it are not asserted.
    assert identifiers(rna) == {(Namespace.NAME, 'gga-miR-17-5p')}
    assert rna.taxon == '9031'
    assert [t for t, _ in annotations(relation)] == ['description']
    assert (gene.namespace, gene.identifier) == ('entrez', '374177')
    assert gene.taxon == '9031'


def test_curated_reader_leaves_predictions_out():
    rows = list(
        iter_interactions(table('\t'.join(SOURCE_HEADER), LNCRNA, MIRNA_FAMILY_NOTATION))
    )
    assert [r['id'] for r in rows] == ['465']


def test_unknown_category_and_organism_are_errors():
    mapper = resource.datasets()['interactions'].mapper
    with pytest.raises(ValueError, match='Unknown RNA category'):
        mapper(parse(UNKNOWN_CATEGORY))
    with pytest.raises(ValueError, match='Unknown organism'):
        mapper(dict(parse(LNCRNA), organism='Ovis aries', tax_id='9940'))


def test_reader_strips_the_category_and_reads_missing_value_marks_as_none():
    row = parse(MIRNA_FAMILY_NOTATION, iter_homology_predictions)
    assert row['rna_category'] == 'miRNA'
    assert row['gene_synonyms'] == []
    assert parse(MIRNA_ONE_PRODUCT)['rna_category'] == 'miRNA'


def test_reader_rejects_a_table_with_changed_columns():
    with pytest.raises(ValueError, match='header differs'):
        list(iter_interactions(table('\t'.join(SOURCE_HEADER[:-1] + ('Other',)), LNCRNA)))


def test_reader_rejects_an_unexpected_pmid_and_an_empty_download():
    bad = LNCRNA.replace('23141928', 'pending')
    with pytest.raises(ValueError, match='PMID'):
        list(iter_interactions(table('\t'.join(SOURCE_HEADER), bad)))
    with pytest.raises(ValueError, match='is empty'):
        list(iter_interactions(table()))
