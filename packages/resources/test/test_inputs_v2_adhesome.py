"""Source contract of the Adhesome input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.adhesome import resource
from pypath.inputs_v2.parsers.adhesome import (
    INTERACTION,
    WEBSITE_COLUMNS,
    iter_interactions,
    iter_website_pairs,
    join_components,
    uniprot_accessions,
)


def component(name, synonyms, gene, uniprot, function, category, classification):
    return {
        'Protein name': name,
        'Synonyms': synonyms,
        'Gene name': gene,
        'swiss prot ID': uniprot,
        'Function': function,
        'Category': category,
        'Clasification': classification,
        'Description': None,
    }


def interaction(source, code, target, pmid):
    # The spreadsheet reader returns the two numeric columns as floats.
    return {'Source': source, INTERACTION: code, 'target': target, 'PMID': pmid}


# Rows of Supplementary Table S1 of Zaidel-Bar et al. 2007.
COMPONENTS = [
    component('Abl', 'c-abl, p150', 'ABL1', 'P00519', 'Tyrosine Kinase', 'Tyrosine Kinase', 'intrinsic'),
    component('paxillin', None, 'PXN', 'P49023', 'adaptor', 'adaptor', 'intrinsic'),
    component('talin', None, 'TLN1', 'Q9Y490', 'Cytoskeletal', 'adaptor', 'intrinsic'),
    component('Calpain', 'calpain-1, Calpain-2, muCANP', 'CAPN1, CAPN2', 'P07384, P17655', 'Cysteine protease', 'Protease', 'associated'),
    component('Abi', 'ABI1, Abl interactor 1', 'ABI1, ABI2, ABI3', 'Q8IZP0, Q9NYB9, Q9P2A4', 'adaptor', 'adaptor', 'associated'),
    component('CIN85', 'CIN85, SH3KBP1', 'SH3KBP1', 'Q96B97', 'adaptor', 'adaptor', 'intrinsic'),
    component('PKC', None, 'PRKCA', 'P17252, Q05655, Q15139, Q02156, Q04759', 'serine/threonine kinase', 'serine/threonine kinase', 'intrinsic'),
    component('DAG', 'Diacylglycerol', None, None, 'lipid', 'lipid', 'associated'),
    component('[Ca+]', 'calcium ions', None, None, None, None, None),
]
# Rows of Supplementary Table S2.
INTERACTIONS = [
    interaction('ABL', 1.0, 'PAXILLIN', 9603926.0),
    interaction('CALPAIN', 2.0, 'TALIN', 15448700.0),
    interaction('ABI', 3.0, 'CIN85', 10858458.0),
    interaction('DAG', 1.0, 'PKC', 10579907.0),
    interaction('[Ca+]', 1.0, 'CALPAIN', 6100833.0),
]
ACTIVATION, INHIBITION, BINDING, LIPID, CALCIUM = join_components(
    INTERACTIONS, COMPONENTS
)


def extract(dataset, source_row):
    record = resource.datasets()[dataset].mapper(source_row)
    assert record is not None
    out = SilverExtractor('adhesome', dataset)
    out.process_record(record, source_row, 'fixture:1', 1)
    relation = next(r for r in out.relations if r.predicate != 'has_member')
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
        out,
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def annotations(relation, scope='relation'):
    return sorted(
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def members(out, group):
    return {
        (member.entity_type, member.identifier, member.taxon)
        for relation in out.relations
        if relation.predicate == 'has_member'
        and relation.subject_entity_key == group.entity_key
        for member in [out.entities[relation.object_entity_key]]
    }


def test_join_reads_the_components_of_an_interaction():
    assert ACTIVATION == {
        'source': 'ABL',
        'target': 'PAXILLIN',
        'interaction': '1',
        'pmid': '9603926',
        'source_name': 'Abl',
        'source_synonyms': 'c-abl, p150',
        'source_gene': 'ABL1',
        'source_uniprot': 'P00519',
        'source_function': 'Tyrosine Kinase',
        'source_category': 'Tyrosine Kinase',
        'source_classification': 'intrinsic',
        'target_name': 'paxillin',
        'target_synonyms': None,
        'target_gene': 'PXN',
        'target_uniprot': 'P49023',
        'target_function': 'adaptor',
        'target_category': 'adaptor',
        'target_classification': 'intrinsic',
    }


def test_activation_is_an_increasing_effect_of_the_source_on_the_target():
    relation, source, target, _ = extract('interactions', ACTIVATION)
    assert relation.predicate == 'affects'
    assert identifiers(source) == {
        (Namespace.UNIPROT, 'P00519'),
        (Namespace.GENESYMBOL, 'ABL1'),
        (Namespace.NAME, 'Abl'),
    }
    assert identifiers(target) == {
        (Namespace.UNIPROT, 'P49023'),
        (Namespace.GENESYMBOL, 'PXN'),
        (Namespace.NAME, 'paxillin'),
    }
    assert (source.namespace, source.identifier) == ('uniprot', 'P00519')
    assert {source.entity_type, target.entity_type} == {'protein'}
    assert {source.taxon, target.taxon} == {'9606'}
    assert annotations(relation) == [
        ('object_direction_qualifier', 'increased'),
        ('publications', 'PMID:9603926'),
    ]


def test_inhibition_by_a_group_of_proteins_is_a_decreasing_effect():
    relation, source, target, out = extract('interactions', INHIBITION)
    assert relation.predicate == 'affects'
    assert annotations(relation) == [
        ('object_direction_qualifier', 'decreased'),
        ('publications', 'PMID:15448700'),
    ]
    # Calpain has two accessions: the component is either protein.
    assert source.entity_type == 'protein_family'
    assert identifiers(source) == {(Namespace.NAME, 'Calpain')}
    assert members(out, source) == {
        ('protein', 'P07384', '9606'),
        ('protein', 'P17655', '9606'),
    }
    assert (target.entity_type, target.identifier) == ('protein', 'Q9Y490')


def test_binding_has_no_direction():
    relation, source, target, out = extract('interactions', BINDING)
    assert relation.predicate == 'interacts_with'
    assert annotations(relation) == [('publications', 'PMID:10858458')]
    assert len(members(out, source)) == 3
    assert (Namespace.GENESYMBOL, 'SH3KBP1') in identifiers(target)


def test_no_aspect_and_no_component_classification_is_asserted():
    relation, _, _, _ = extract('interactions', ACTIVATION)
    assert {a['term'] for a in relation.annotations} == {
        'in_taxon',
        'object_direction_qualifier',
        'publications',
    }
    values = {a['value'] for a in relation.annotations}
    assert not values & {'intrinsic', 'adaptor', 'Tyrosine Kinase'}


def test_lipid_and_calcium_are_chemicals_without_a_taxon():
    relation, lipid, kinases, out = extract('interactions', LIPID)
    assert relation.predicate == 'affects'
    assert lipid.entity_type == 'chemical_entity'
    assert identifiers(lipid) == {
        (Namespace.NAME, 'DAG'),
        (Namespace.SYNONYM, 'Diacylglycerol'),
    }
    assert lipid.taxon is None
    assert annotations(relation, 'subject') == []
    # The gene name of a group (PRKCA) names one of its five proteins only.
    assert identifiers(kinases) == {(Namespace.NAME, 'PKC')}
    assert len(members(out, kinases)) == 5
    _, calcium, _, _ = extract('interactions', CALCIUM)
    assert calcium.entity_type == 'chemical_entity'
    assert identifiers(calcium) == {
        (Namespace.NAME, '[Ca+]'),
        (Namespace.SYNONYM, 'calcium ions'),
    }


def test_unknown_interaction_code_is_rejected():
    mapper = resource.datasets()['interactions'].mapper
    with pytest.raises(ValueError, match='interaction code'):
        mapper(dict(ACTIVATION, interaction='4'))
    with pytest.raises(ValueError, match='interaction code'):
        mapper(dict(ACTIVATION, interaction=None))


def test_no_relation_for_a_component_the_table_does_not_describe():
    (row,) = join_components(
        [interaction('ABL', 1.0, 'TENSIN', 9211900.0)], COMPONENTS
    )
    assert row['target'] == 'TENSIN' and row['target_name'] is None
    assert resource.datasets()['interactions'].mapper(row) is None
    # A protein component without an accession is not identified either.
    mapper = resource.datasets()['interactions'].mapper
    assert mapper(dict(ACTIVATION, target_uniprot=None)) is None


def test_join_rejects_tables_with_changed_columns():
    with pytest.raises(ValueError, match='swiss prot ID'):
        list(join_components(INTERACTIONS, [{'Protein name': 'Abl'}]))
    with pytest.raises(ValueError, match='PMID'):
        list(join_components([{'Source': 'ABL', 'target': 'PAXILLIN'}], COMPONENTS))
    with pytest.raises(ValueError, match='no rows'):
        list(join_components([], COMPONENTS))


def website(lines):
    text = '\n'.join([','.join(WEBSITE_COLUMNS), *lines]) + '\n'
    return SimpleNamespace(result=io.StringIO(text))


# Lines of adhesome.org's interactions.csv.
WEBSITE_LINES = [
    # The 2007 network, in both orientations: read from the publication.
    '"ABL1","CBL","+","Activation","7925293","adhesome.org"',
    '"CBL","ABL1","+","Activation","7925293","adhesome.org"',
    '"ITGB1","LYN","0","Binding","14593208","wangz10_ppi"',
    '"LYN","ITGB1","0","Binding","14593208","wangz10_ppi"',
    '"PTPN1","INSR","_","Dephosphorylation","14722096","axel_ppi"',
    '"INSR","PTPN1","0","Binding","14722096","axel_ppi"',
    '"ITGB1","LYN","0","colocalization","20159608","axel_ppi"',
    '"LYN","ITGB1","0","colocalization","20159608","axel_ppi"',
    '"PTK2B","SRC","0","Unknown","12345678","axel_ppi"',
    '"SRC","PTK2B","0","Unknown","12345678","axel_ppi"',
    '"CRK","ELMO1","0","Binding","0","wangz10_ppi"',
    '"ELMO1","CRK","0","Binding","0","wangz10_ppi"',
]


def test_website_reader_merges_mirrored_rows_and_leaves_out_the_2007_network():
    assert list(iter_website_pairs(website(WEBSITE_LINES))) == [
        {
            'component_a': 'ITGB1',
            'component_b': 'LYN',
            'PMID': '14593208',
            'types': ['Binding'],
            'effects': ['0'],
            'data_sources': ['wangz10_ppi'],
        },
        {
            'component_a': 'INSR',
            'component_b': 'PTPN1',
            'PMID': '14722096',
            'types': ['Binding', 'Dephosphorylation'],
            'effects': ['0', '_'],
            'data_sources': ['axel_ppi'],
        },
        {
            'component_a': 'ITGB1',
            'component_b': 'LYN',
            'PMID': '20159608',
            'types': ['colocalization'],
            'effects': ['0'],
            'data_sources': ['axel_ppi'],
        },
        {
            'component_a': 'PTK2B',
            'component_b': 'SRC',
            'PMID': '12345678',
            'types': ['Unknown'],
            'effects': ['0'],
            'data_sources': ['axel_ppi'],
        },
        {
            'component_a': 'CRK',
            'component_b': 'ELMO1',
            'PMID': '0',
            'types': ['Binding'],
            'effects': ['0'],
            'data_sources': ['wangz10_ppi'],
        },
    ]


def test_website_reader_rejects_a_table_with_changed_columns():
    lines = SimpleNamespace(result=io.StringIO('Source,Target,PMID\nA,B,1\n'))
    with pytest.raises(ValueError, match='data_source'):
        list(iter_website_pairs(lines))
    with pytest.raises(ValueError, match='no rows'):
        list(iter_website_pairs(website([])))


BOUND, SIGNED, COLOCALIZED, UNKNOWN, NO_PMID = iter_website_pairs(
    website(WEBSITE_LINES)
)


def test_website_interaction_is_between_gene_symbols_with_its_publication():
    relation, first, second, _ = extract('website_interactions', BOUND)
    assert relation.predicate == 'interacts_with'
    assert identifiers(first) == {(Namespace.GENESYMBOL, 'ITGB1')}
    assert identifiers(second) == {(Namespace.GENESYMBOL, 'LYN')}
    assert {first.entity_type, second.entity_type} == {'protein'}
    assert {first.taxon, second.taxon} == {'9606'}
    assert annotations(relation) == [('publications', 'PMID:14593208')]


def test_website_effect_gets_no_direction():
    # The file does not say which of the two is the phosphatase or the
    # regulator, so neither the sign nor the modification is asserted.
    relation, _, _, _ = extract('website_interactions', SIGNED)
    assert relation.predicate == 'interacts_with'
    assert annotations(relation) == [('publications', 'PMID:14722096')]


def test_website_colocalization_is_not_an_interaction():
    relation, _, _, _ = extract('website_interactions', COLOCALIZED)
    assert relation.predicate == 'colocalizes_with'
    assert annotations(relation) == [('publications', 'PMID:20159608')]


def test_website_interaction_of_unknown_kind_is_not_emitted():
    mapper = resource.datasets()['website_interactions'].mapper
    assert mapper(UNKNOWN) is None
    with pytest.raises(ValueError, match='interaction type'):
        mapper(dict(BOUND, types=['Ubiquitination']))
    with pytest.raises(ValueError, match='interaction type'):
        mapper(dict(BOUND, types=[]))


@pytest.mark.parametrize(
    'pmid', ['0', '12345678', '94349924', '', '0123456', '123456789']
)
def test_number_that_is_not_a_pubmed_id_is_not_a_publication(pmid):
    relation, _, _, _ = extract('website_interactions', dict(NO_PMID, PMID=pmid))
    assert relation.predicate == 'interacts_with'
    assert annotations(relation) == []


@pytest.mark.parametrize('pmid', ['0123456', '123456789'])
def test_2007_network_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(
    pmid,
):
    relation, _, _, _ = extract('interactions', dict(ACTIVATION, pmid=pmid))
    assert ('publications', 'PMID:9603926') in annotations(
        extract('interactions', ACTIVATION)[0]
    )
    assert not [term for term, _ in annotations(relation) if term == 'publications']


def test_accession_cell_is_split_and_only_accessions_are_kept():
    assert uniprot_accessions('P00519') == ['P00519']
    assert uniprot_accessions('P07384, P17655') == ['P07384', 'P17655']
    assert uniprot_accessions('P07384, null,  ') == ['P07384']
    assert uniprot_accessions(None) == []
    assert uniprot_accessions('') == []


def test_dataset_joins_the_component_table_from_its_download():
    from pypath.inputs_v2 import adhesome

    raw_parser = adhesome.resource.interactions._raw_parser
    assert raw_parser.func is iter_interactions
    assert raw_parser.keywords == {'components_download': adhesome.components_download}
    with pytest.raises(ValueError, match='component table'):
        list(iter_interactions(None))
