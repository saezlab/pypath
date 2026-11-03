"""Source contract of the Cancer Cell Map input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.cancercellmap import resource
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.parsers.cancercellmap import (
    EDGE_COLUMNS,
    NODE_COLUMNS,
    join_nodes,
)

# Rows of cell-map-node-attributes.txt.
NODES = (
    ('101318', 'MAP3K4', 'Q5VTT9:Q5VTT8:A6H8W0:B7ZLD3:Q9Y6R4:Q92612:Q9H408:B9EG75', '4216', 'NOT_SPECIFIED', 'protein', '9606'),
    ('140665', 'AXIN1', 'Q96S28:Q86YA7:Q4TT27:O15169:Q8WVW6:Q4TT26', '8312', 'NOT_SPECIFIED', 'protein', '9606'),
    ('101498', 'MET', 'Q12875:O60366:P08581:Q9UDX7:Q9UPL8', '4233', 'NOT_SPECIFIED', 'protein', '9606'),
    ('5035', 'ITGB4', 'Q9UIQ4:O15339:P16144:O14690:O14691:O15340:O15341', '3691', 'NOT_SPECIFIED', 'protein', '9606'),
    ('101582', 'MAP2K2', 'P36507', '5605', 'NOT_SPECIFIED', 'protein', '9606'),
    ('14516', 'EGFR', 'Q9H2C9:O00688:Q92795:Q9UMD7:Q9UMD8:Q9UMG5:O00732:Q14225:Q9GZX1:Q68GS5:P06268:P00533:Q9H3C9:Q9BZS2', '1956', 'NOT_SPECIFIED', 'protein', '9606'),
    ('32635', 'CSK', 'Q2M3N2:Q6FGZ6:P41240', '1445', 'NOT_SPECIFIED', 'protein', '9606'),
    ('103193', 'TGFBR2', 'Q99474:P37173:Q6DKT6:B4DTV5:Q15580', '7048', 'NOT_SPECIFIED', 'protein', '9606'),
    # A record of "Mammalia" without a gene symbol.
    ('476621', 'NOT_SPECIFIED', 'C0H5W9', '58533', 'NOT_SPECIFIED', 'protein', '40674'),
    ('123691', 'KIT', 'P10721:Q9UM99:Q99662', '3815', 'NOT_SPECIFIED', 'protein', '9606'),
    # A GenBank protein accession in the UniProt column and nothing else.
    ('543350', 'NOT_SPECIFIED', 'BAC32181', 'NOT_SPECIFIED', 'NOT_SPECIFIED', 'protein', '10090'),
    ('52090', 'AR', 'P10275:Q9UD95:A2RUN2:B1AKD7', '367','NOT_SPECIFIED', 'protein', '9606'),
    ('542700', 'NOT_SPECIFIED', 'NOT_SPECIFIED', 'NOT_SPECIFIED', 'NOT_SPECIFIED', 'protein', '9606'),
    ('102124', 'LYN', 'A0AVQ5:P07948', '4067', 'NOT_SPECIFIED', 'protein', '9606'),
    ('14194', 'DOK1', 'Q53TY2:Q99704:O43204:Q9UHG6', '1796', 'NOT_SPECIFIED', 'protein', '9606'),
    ('103317', 'TFDP2', 'Q8WU88:Q6R754:D3DNG1:Q13331:Q14187:Q14188', '7029', 'NOT_SPECIFIED', 'protein', '9606'),
    ('14354', 'E2F4', 'Q12991:Q16254:B5BU56:A6NGR8:Q15328', '1874', 'NOT_SPECIFIED', 'protein', '9606'),
)
# Rows of cell-map-edge-attributes.txt.
INTERACTS = ('101318', 'INTERACTS_WITH', '140665', 'MAP3K4', 'AXIN1', 'CELL_MAP', 'PUBMED:15262978;')
INTERACTS_MIRROR = ('140665', 'INTERACTS_WITH', '101318', 'AXIN1', 'MAP3K4', 'CELL_MAP', 'PUBMED:15262978;')
STATE_CHANGE = ('101498', 'STATE_CHANGE', '5035', 'MET', 'ITGB4', 'CELL_MAP', 'PUBMED:11733063;')
STATE_CHANGE_SELF = ('101582', 'STATE_CHANGE', '101582', 'MAP2K2', 'MAP2K2', 'CELL_MAP', 'PUBMED:8226933;')
NO_PMID = ('14516', 'INTERACTS_WITH', '32635', 'EGFR', 'CSK', 'CELL_MAP', 'null:null;PUBMED:Literature;')
MAMMALIA = ('103193', 'INTERACTS_WITH', '476621', 'TGFBR2', 'NOT_SPECIFIED', 'CELL_MAP', 'PUBMED:11279102;')
GENBANK_ONLY = ('123691', 'INTERACTS_WITH', '543350', 'KIT', 'NOT_SPECIFIED', 'CELL_MAP', 'PUBMED:12177062;')
NO_IDENTIFIER = ('52090', 'INTERACTS_WITH', '542700', 'AR', 'NOT_SPECIFIED', 'CELL_MAP', 'PUBMED:12941817;')
SAME_COMPONENT = ('102124', 'IN_SAME_COMPONENT', '14194', 'LYN', 'DOK1', 'CELL_MAP', 'NOT_SPECIFIED')
REACTS = ('103317', 'REACTS_WITH', '14354', 'TFDP2', 'E2F4', 'CELL_MAP', 'PUBMED:9144196;')

mapper = resource.datasets()['interactions'].mapper


def table(columns, lines):
    text = '\n'.join('\t'.join(line) for line in (columns, *lines)) + '\n'
    return SimpleNamespace(result=io.StringIO(text))


def read(*edges, nodes=NODES):
    return list(
        join_nodes(
            iter_tsv(table(EDGE_COLUMNS, edges)),
            iter_tsv(table(NODE_COLUMNS, nodes)),
        )
    )


def extract(source_row):
    record = mapper(source_row)
    assert record is not None
    out = SilverExtractor('cancercellmap', 'interactions')
    out.process_record(record, source_row, 'fixture:1', 1)
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
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def test_reader_adds_the_node_attributes_of_both_participants():
    (row,) = read(INTERACTS)
    assert row['UNIPROT_ACCESSION_B'] == 'Q96S28:Q86YA7:Q4TT27:O15169:Q8WVW6:Q4TT26'
    assert row['ENTREZ_GENE_ID_A'] == '4216'
    assert row['NCBI_TAX_ID_A'] == '9606'
    assert row['NODE_TYPE_B'] == 'protein'
    assert row['MIRROR'] is False


def test_interacts_with_is_an_undirected_interaction_with_its_publication():
    relation, first, second = extract(read(INTERACTS)[0])
    assert relation.predicate == 'interacts_with'
    assert identifiers(first) == {
        (Namespace.UNIPROT, 'Q5VTT9'),
        (Namespace.UNIPROT, 'Q5VTT8'),
        (Namespace.UNIPROT, 'A6H8W0'),
        (Namespace.UNIPROT, 'B7ZLD3'),
        (Namespace.UNIPROT, 'Q9Y6R4'),
        (Namespace.UNIPROT, 'Q92612'),
        (Namespace.UNIPROT, 'Q9H408'),
        (Namespace.UNIPROT, 'B9EG75'),
        (Namespace.ENTREZ, '4216'),
        (Namespace.GENESYMBOL, 'MAP3K4'),
    }
    assert (Namespace.GENESYMBOL, 'AXIN1') in identifiers(second)
    assert (Namespace.ENTREZ, '8312') in identifiers(second)
    assert {first.entity_type, second.entity_type} == {'protein'}
    assert {first.taxon, second.taxon} == {'9606'}
    assert annotations(relation) == [('publications', 'PMID:15262978')]


def test_state_change_is_directed_without_sign_or_aspect():
    relation, controller, changed = extract(read(STATE_CHANGE)[0])
    assert relation.predicate == 'affects'
    assert (Namespace.GENESYMBOL, 'MET') in identifiers(controller)
    assert (Namespace.GENESYMBOL, 'ITGB4') in identifiers(changed)
    # The rule covers any change of state, so no qualifier is asserted.
    assert {a['term'] for a in relation.annotations} == {
        'in_taxon',
        'publications',
    }
    assert annotations(relation) == [('publications', 'PMID:11733063')]


def test_state_change_of_a_protein_by_itself_is_kept():
    relation, controller, changed = extract(read(STATE_CHANGE_SELF)[0])
    assert relation.predicate == 'affects'
    assert relation.subject_entity_key == relation.object_entity_key
    assert identifiers(controller) == {
        (Namespace.UNIPROT, 'P36507'),
        (Namespace.ENTREZ, '5605'),
        (Namespace.GENESYMBOL, 'MAP2K2'),
    }


def test_second_orientation_of_an_undirected_pair_is_not_emitted():
    first, second = read(INTERACTS, INTERACTS_MIRROR)
    assert (first['MIRROR'], second['MIRROR']) == (False, True)
    assert mapper(first) is not None
    assert mapper(second) is None


def test_same_pair_with_another_publication_is_not_a_mirror():
    other = (*INTERACTS_MIRROR[:6], 'PUBMED:9305638;')
    assert [row['MIRROR'] for row in read(INTERACTS, other)] == [False, False]


def test_reversed_state_change_is_a_statement_of_its_own():
    reverse = ('5035', 'STATE_CHANGE', '101498', 'ITGB4', 'MET', 'CELL_MAP', 'PUBMED:11733063;')
    rows = read(STATE_CHANGE, reverse)
    assert [row['MIRROR'] for row in rows] == [False, False]
    assert all(mapper(row) is not None for row in rows)


@pytest.mark.parametrize('edge', [SAME_COMPONENT, REACTS])
def test_complex_and_conversion_projections_are_not_interactions(edge):
    (row,) = read(edge)
    assert mapper(row) is None


def test_placeholder_reference_is_not_a_publication():
    relation, _, _ = extract(read(NO_PMID)[0])
    assert annotations(relation) == []


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    edge = (*INTERACTS[:6], 'PUBMED:0123456;PUBMED:123456789;PUBMED:15262978;')
    relation, _, _ = extract(read(edge)[0])
    assert annotations(relation) == [('publications', 'PMID:15262978')]


def test_record_keeps_the_taxon_the_source_gives_it():
    relation, _, second = extract(read(MAMMALIA)[0])
    assert identifiers(second) == {
        (Namespace.UNIPROT, 'C0H5W9'),
        (Namespace.ENTREZ, '58533'),
    }
    assert second.taxon == '40674'
    assert annotations(relation, 'object') == [('in_taxon', 'NCBITaxon:40674')]


@pytest.mark.parametrize('edge', [NO_IDENTIFIER, GENBANK_ONLY])
def test_no_relation_for_a_record_without_an_identifier(edge):
    (row,) = read(edge)
    assert mapper(row) is None


def test_pathway_commons_record_id_is_not_an_identifier():
    _, first, second = extract(read(INTERACTS)[0])
    values = {value for _, value in identifiers(first) | identifiers(second)}
    assert not values & {'101318', '140665'}


def test_unknown_interaction_type_gets_no_default():
    (row,) = read(INTERACTS)
    with pytest.raises(ValueError, match='CO_CONTROL'):
        mapper(dict(row, INTERACTION_TYPE='CO_CONTROL'))


def test_unknown_node_type_is_rejected():
    (row,) = read(INTERACTS)
    with pytest.raises(ValueError, match='small_molecule'):
        mapper(dict(row, NODE_TYPE_B='small_molecule'))


def test_reader_rejects_an_edge_without_its_node():
    with pytest.raises(ValueError, match='140665'):
        read(INTERACTS, nodes=NODES[:1])


def test_reader_rejects_tables_with_changed_columns():
    with pytest.raises(ValueError, match='INTERACTION_PUBMED_ID'):
        list(
            join_nodes(
                iter_tsv(table(EDGE_COLUMNS[:6], [INTERACTS[:6]])),
                iter_tsv(table(NODE_COLUMNS, NODES)),
            )
        )
    with pytest.raises(ValueError, match='NCBI_TAX_ID'):
        list(
            join_nodes(
                iter_tsv(table(EDGE_COLUMNS, [INTERACTS])),
                iter_tsv(table(NODE_COLUMNS[:6], [node[:6] for node in NODES])),
            )
        )


def test_reader_rejects_empty_downloads():
    with pytest.raises(ValueError, match='edge table has no rows'):
        read()
    with pytest.raises(ValueError, match='node table has no rows'):
        read(INTERACTS, nodes=())
