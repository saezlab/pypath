"""Source contract of the miRDeathDB microRNA-target input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.mirdeathdb import resource
from pypath.inputs_v2.parsers.mirdeathdb import COLUMNS, iter_interactions


def line(*values):
    assert len(values) == len(COLUMNS)
    return '\t'.join(values)


# Rows of miRDeathDB_all_data.txt. Long cells are shortened.
SINGLE_PRODUCT = line(
    'miR-1', 'hsa-miR-1', 'MIMAT0000416', 'FOXP1', 'apoptosis', 'up', 'human',
    'lung cancer cell lines', '18818206', '9606', '27086', '12CC4|QRF1|hFKH1B',
    'HGNC:3823|MIM:605515|Ensembl:ENSG00000114861|HPRD:18518|M2D:hsa-miR-1',
    '3', '3p14.1', 'forkhead box P1', 'protein-coding', 'forkhead box P1',
    'forkhead box protein P1',
)
TWO_PRODUCTS = line(
    'miR-221', '"hsa-miR-221-5p,hsa-miR-221-3p"', '"MIMAT0004568,MIMAT0000278"',
    'PTEN', 'apoptosis', 'down', 'human', 'Gastric cancer', '20618998', '9606',
    '5728', '10q23del|BZS|DEC', 'HGNC:9588|MIM:601728|Ensembl:ENSG00000171862',
    '10', '10q23.3', 'phosphatase and tensin homolog', 'protein-coding',
    'phosphatase and tensin homolog', 'MMAC1',
)
VIRAL_PRODUCT = line(
    'miR-BART5', 'ebv-miR-BART5', 'MIMAT0003413', 'BBC3', 'apoptosis', 'down',
    'human', 'EBV-nasopharyngeal carcinoma', '18838543', '9606', '27113',
    'JFY-1|JFY1|PUMA', 'HGNC:17868|MIM:605854|Ensembl:ENSG00000105327',
    '19', '19q13.3-q13.4', 'BCL2 binding component 3', 'protein-coding',
    'BCL2 binding component 3', 'bcl-2-binding component 3',
)
NO_GENE = line(
    'miR-31', '"mmu-miR-31-5p,mmu-miR-31-3p"', '"MIMAT0000538,MIMAT0004634"',
    'NA', 'autophagy', 'up', 'mouse', 'caveolin-1(-/-) null stromal cells',
    '20861672', '10090', '', '', 'M2D:hsa-miR-31', '', '', '', '', '', '',
)
SYMBOL_ONLY_GENE = line(
    'miR-24a', '"xtr-miR-24a-5p,xtr-miR-24a-3p"', '"MIMAT0003653,MIMAT0003654 "',
    'APAF1', 'apoptosis', 'down', 'Xenopus', 'neural retina', '19372388',
    '262014', '', '', 'M2D:hsa-miR-24', '', '', '', '', '', '',
)
NO_PMID = line(
    'miR-30b', '"hsa-miR-30b-5p,hsa-miR-30b-3p"', '"MIMAT0000420,MIMAT0004589 "',
    'ATG12', 'autophagy', 'down', 'human', 'Helicobacter pylori infected cells',
    '', '9606', '9140', 'APG12|APG12L', 'HGNC:588|MIM:609608', '5', '5q21-q22',
    'ATG12 autophagy related 12 homolog (S. cerevisiae)', 'protein-coding',
    'ATG12 autophagy related 12 homolog (S. cerevisiae)', 'autophagy-related protein 12',
)


def table(*lines):
    return SimpleNamespace(result=io.StringIO('\r\n'.join(lines) + '\r\n'))


def parse(text):
    (row,) = list(iter_interactions(table('\t'.join(COLUMNS), text)))
    return row


def extract(text):
    row = parse(text)
    record = resource.datasets()['interactions'].mapper(row)
    assert record is not None
    out = SilverExtractor('mirdeathdb', 'interactions')
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


def test_single_product_symbol_is_the_mirbase_mature_microrna():
    relation, mirna, gene = extract(SINGLE_PRODUCT)
    assert relation.predicate == 'interacts_with'
    assert mirna.entity_type == 'microrna'
    assert (mirna.namespace, mirna.identifier) == ('mirbase_mature', 'MIMAT0000416')
    assert identifiers(mirna) == {
        (Namespace.MIRBASE_MATURE, 'MIMAT0000416'),
        (Namespace.NAME, 'hsa-miR-1'),
        (Namespace.SYNONYM, 'miR-1'),
    }
    assert mirna.taxon == '9606'
    assert gene.entity_type == 'gene'
    assert (gene.namespace, gene.identifier) == ('entrez', '27086')
    assert identifiers(gene) == {
        (Namespace.ENTREZ, '27086'),
        (Namespace.HGNC, '3823'),
        (Namespace.ENSG, 'ENSG00000114861'),
        (Namespace.GENESYMBOL, 'FOXP1'),
        (Namespace.NAME, 'forkhead box P1'),
        (Namespace.SYNONYM, '12CC4'),
        (Namespace.SYNONYM, 'QRF1'),
        (Namespace.SYNONYM, 'hFKH1B'),
    }
    assert gene.taxon == '9606'
    assert annotations(relation) == [('publications', 'PMID:18818206')]


def test_symbol_with_two_products_is_not_expanded_to_either():
    relation, mirna, _ = extract(TWO_PRODUCTS)
    # The source names miR-221, not one of its arms: no accession, no pair per arm.
    assert identifiers(mirna) == {(Namespace.NAME, 'miR-221')}
    assert mirna.taxon == '9606'
    assert all(
        i['ns'] != 'mirbase_mature' for i in mirna.identifiers
    )


def test_pathway_tissue_and_expression_direction_are_not_on_the_relation():
    relation, _, _ = extract(TWO_PRODUCTS)
    values = {v for _, v in annotations(relation)}
    assert 'apoptosis' not in values
    assert 'down' not in values
    assert 'Gastric cancer' not in values
    assert relation.predicate == 'interacts_with'
    assert {a['term'] for a in relation.annotations} <= {
        'publications',
        'in_taxon',
    }


def test_viral_microrna_keeps_its_accession_but_not_the_host_taxon():
    _, mirna, gene = extract(VIRAL_PRODUCT)
    assert (mirna.namespace, mirna.identifier) == ('mirbase_mature', 'MIMAT0003413')
    assert mirna.taxon is None
    assert gene.taxon == '9606'


def test_gene_without_entrez_id_is_identified_by_symbol():
    _, _, gene = extract(SYMBOL_ONLY_GENE)
    assert (gene.namespace, gene.identifier) == ('genesymbol', 'APAF1')
    assert gene.taxon == '262014'


def test_row_without_a_target_gene_gives_no_relation():
    row = parse(NO_GENE)
    assert row['gene_symbol'] is None
    assert resource.datasets()['interactions'].mapper(row) is None


@pytest.mark.parametrize('number', ['0123456', '123456789'])
def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(number):
    text = SINGLE_PRODUCT.replace('18818206', number)
    assert parse(text)['pmids'] == []
    relation, _, _ = extract(text)
    assert annotations(relation) == []


def test_row_without_pmid_has_no_publication():
    relation, _, _ = extract(NO_PMID)
    assert annotations(relation) == []


def test_unknown_organism_is_an_error():
    row = dict(parse(SINGLE_PRODUCT), organism='zebrafish', tax_id='7955')
    with pytest.raises(ValueError, match='Unknown organism'):
        resource.datasets()['interactions'].mapper(row)


def test_reader_splits_the_product_lists():
    row = parse(TWO_PRODUCTS)
    assert row['mature_names'] == ['hsa-miR-221-5p', 'hsa-miR-221-3p']
    assert row['mature_accessions'] == ['MIMAT0004568', 'MIMAT0000278']
    assert row['mirbase_mature'] is None
    assert row['mature_organism_codes'] == ['hsa']
    assert row['pmids'] == ['20618998']
    assert row['expression_direction'] == 'down'


def test_reader_rejects_a_table_with_changed_columns():
    with pytest.raises(ValueError, match='header differs'):
        list(iter_interactions(table('\t'.join(COLUMNS[:-1] + ('Other',)), SINGLE_PRODUCT)))


def test_reader_rejects_a_non_numeric_pmid():
    bad = SINGLE_PRODUCT.replace('18818206', 'prediction')
    with pytest.raises(ValueError, match='PMID'):
        list(iter_interactions(table('\t'.join(COLUMNS), bad)))


def test_reader_rejects_an_empty_download():
    with pytest.raises(ValueError, match='no rows'):
        list(iter_interactions(table('\t'.join(COLUMNS))))
