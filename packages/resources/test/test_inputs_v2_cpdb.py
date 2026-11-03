"""ConsensusPathDB human protein interactions: parser, mapper and extractor on source lines."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import cpdb
from pypath.inputs_v2.parsers import cpdb as parser

HEADER = (
    '#  ConsensusPathDB (version 35) list of human protein interactions\n'
    '#  source_databases\tinteraction_publications\t'
    'interaction_participants__uniprot_entry\tinteraction_participants__uniprot_id\t'
    'interaction_participants__genename\tinteraction_participants__hgnc_id\t'
    'interaction_participants__entrez_gene\tinteraction_participants__ensembl_gene\t'
    'interaction_confidence\n'
)
# Lines copied from ConsensusPathDB_human_PPI.gz, release 35.
PAIR = (
    'Reactome\t\tSUMF2_HUMAN,SUMF1_HUMAN\tQ8NBJ7,Q8NBK3\tSUMF2,SUMF1\t'
    'HGNC:20415,HGNC:20376\t25870,285362\tENSG00000129103,ENSG00000144455\t'
    '0.99805100000000002\n'
)
# "unassigned1304" is no PubMed ID.
PAIR_MINT = (
    'CORUM,MINT,Biogrid\tunassigned1304,22939629,11113176\tT2EA_HUMAN,T2EB_HUMAN\t'
    'P29083,P29084\tGTF2E1,GTF2E2\tHGNC:4650,HGNC:4651\t2960,2961\t'
    'ENSG00000153767,ENSG00000197265\t0.97316999999999998\n'
)
SELF = (
    'Reactome,PhosphoPOINT,HPRD\t22307056,12840032\tMK03_HUMAN\tP27361\tMAPK3\t'
    'HGNC:6877\t5595\tENSG00000102882\tNA\n'
)
COMPLEX = (
    'Reactome\t18627462,15485916,23645155\tDPH1_HUMAN,DPH2_HUMAN,DPH3_HUMAN\t'
    'Q9BZG8,Q9BQC3,Q96FX2\tDPH1,DPH2,DPH3\tHGNC:3003,HGNC:3004,HGNC:27717\t'
    '1801,1802,285381\tENSG00000108963,ENSG00000132768,ENSG00000154813\tNA\n'
)
# The second participant is two UniProt entries; it is not known which of
# them interacts.
AMBIGUOUS_PAIR = (
    'Reactome\t21325273\tVENTX_HUMAN,ARF_HUMAN.CDN2A_HUMAN\tO95231,Q8N726.P42771\t'
    'VENTX,CDKN2A\tHGNC:13639,HGNC:1787\t27287,1029\t'
    'ENSG00000151650,ENSG00000147889\tNA\n'
)
AMBIGUOUS_COMPLEX = (
    'Reactome\t17599906,17461779\tMTOR_HUMAN,SIN1_HUMAN,PRR5_HUMAN.RHG08_HUMAN,'
    'LST8_HUMAN,RICTR_HUMAN\tP42345,Q9BPZ7,P85299.P85298,Q9BVC4,Q6R327\t'
    'MTOR,MAPKAP1,ARHGAP8.PRR5,MLST8,RICTOR\t'
    'HGNC:3942,HGNC:18752,HGNC:31682.HGNC:677,HGNC:24825,HGNC:28611\t'
    '2475,79109,23779.55615,64223,253260\t'
    'ENSG00000198793,ENSG00000119487,ENSG00000186654.ENSG00000241484,'
    'ENSG00000167965,ENSG00000164327\tNA\n'
)
# One protein with several genes (the protein is unambiguous, its genes are
# not) and publication tokens that are not PubMed IDs.
PROTEIN_OF_TWO_GENES = (
    'IntAct\tcccc,29477555,888800000002\tCCL4_HUMAN,HTRA3_HUMAN\tP13236,P83110\t'
    'CCL4,HTRA3\tHGNC:10630,HGNC:30406\t388372.6351,94031\t'
    'ENSG00000275302.ENSG00000277943,ENSG00000170801\t0.115468\n'
)
TEXT = HEADER + ''.join(
    [PAIR, PAIR_MINT, SELF, COMPLEX, AMBIGUOUS_PAIR, AMBIGUOUS_COMPLEX, PROTEIN_OF_TWO_GENES]
)


def read(text=TEXT):
    return list(parser.iter_interactions(SimpleNamespace(result=io.StringIO(text))))


def map_rows(dataset, rows):
    return [cpdb.resource.datasets()[dataset].mapper(row) for row in rows]


def extract(dataset, rows):
    out = SilverExtractor('cpdb', dataset)
    for index, row in enumerate(rows):
        record = cpdb.resource.datasets()[dataset].mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def values(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def identifiers(entity):
    return [(i.type, i.value) for i in entity.identifiers]


def test_module_import_defines_the_two_datasets_of_the_one_file():
    assert list(cpdb.resource.datasets()) == ['interactions', 'complexes']
    urls = {d.download.url for d in cpdb.resource.datasets().values()}
    assert urls == {'http://cpdb.molgen.mpg.de/download/ConsensusPathDB_human_PPI.gz'}


def test_parser_keeps_the_participant_columns_aligned_and_the_dots_as_alternatives():
    rows = read()
    assert len(rows) == 7
    assert rows[0] == {
        'source_databases': ['Reactome'],
        'publications': [],
        'participants': [
            {
                'uniprot_entry': ['SUMF2_HUMAN'],
                'uniprot_id': ['Q8NBJ7'],
                'gene_name': ['SUMF2'],
                'hgnc_id': ['HGNC:20415'],
                'entrez_gene': ['25870'],
                'ensembl_gene': ['ENSG00000129103'],
            },
            {
                'uniprot_entry': ['SUMF1_HUMAN'],
                'uniprot_id': ['Q8NBK3'],
                'gene_name': ['SUMF1'],
                'hgnc_id': ['HGNC:20376'],
                'entrez_gene': ['285362'],
                'ensembl_gene': ['ENSG00000144455'],
            },
        ],
        'confidence': '0.99805100000000002',
    }
    assert rows[2]['confidence'] is None
    assert rows[4]['participants'][1]['uniprot_entry'] == ['ARF_HUMAN', 'CDN2A_HUMAN']
    assert rows[1]['publications'] == ['unassigned1304', '22939629', '11113176']


def test_parser_raises_on_a_changed_header_an_empty_table_and_unaligned_columns():
    with pytest.raises(ValueError, match='lacks expected columns'):
        read(HEADER.replace('interaction_confidence', 'confidence') + PAIR)
    with pytest.raises(ValueError, match='no rows'):
        read(HEADER)
    unaligned = PAIR.replace('Q8NBJ7,Q8NBK3', 'Q8NBJ7')
    with pytest.raises(ValueError, match='not aligned'):
        read(HEADER + unaligned)


def test_pair_is_an_interaction_not_a_direct_one_with_sources_publications_and_score():
    out = extract('interactions', read()[:2])
    first, second = out.relations
    assert first.predicate == 'interacts_with'
    assert values(first, 'publications') == []
    assert values(first, 'supporting_data_source') == ['Reactome']
    (score,) = [a for a in first.annotations if a['term'] == 'has_confidence_score']
    assert score['quantity']['has_numeric_value'] == pytest.approx(0.998051)
    assert score['quantity']['source_field'] == 'interaction_confidence'
    # Every database the line names, as the file spells it; "unassigned1304"
    # is no publication.
    assert values(second, 'supporting_data_source') == ['CORUM', 'MINT', 'Biogrid']
    assert values(second, 'publications') == ['PMID:22939629', 'PMID:11113176']
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}
    subject = cpdb.resource.interactions.mapper(read()[0]).subject
    assert identifiers(subject) == [
        (Namespace.UNIPROT, 'Q8NBJ7'),
        (Namespace.UNIPROT_ENTRY, 'SUMF2_HUMAN'),
        (Namespace.GENESYMBOL, 'SUMF2'),
        (Namespace.HGNC, '20415'),
        (Namespace.ENTREZ, '25870'),
        (Namespace.ENSG, 'ENSG00000129103'),
    ]


def test_line_with_one_protein_is_a_self_interaction_without_a_score():
    (relation,) = extract('interactions', [read()[2]]).relations
    assert relation.predicate == 'interacts_with'
    assert relation.subject_entity_key == relation.object_entity_key
    assert not values(relation, 'has_confidence_score')
    assert values(relation, 'supporting_data_source') == [
        'Reactome',
        'PhosphoPOINT',
        'HPRD',
    ]


def test_line_with_more_than_two_proteins_is_a_complex_never_pairs():
    rows = read()
    out = extract('complexes', [rows[3]])
    (complex_,) = [e for e in out.entities.values() if e.entity_type == 'macromolecular_complex']
    # The file gives the complex no accession and no name.
    assert not complex_.identifiers
    assert [r.predicate for r in out.relations] == ['has_member'] * 3
    assert not [r for r in out.relations if r.predicate == 'interacts_with']
    assert {(a['term'], a['value']) for a in complex_.annotations} == {
        ('supporting_data_source', 'Reactome'),
        ('publications', 'PMID:18627462'),
        ('publications', 'PMID:15485916'),
        ('publications', 'PMID:23645155'),
    }
    # The interactions dataset does not map it, the complexes dataset maps
    # nothing else.
    assert cpdb.resource.interactions.mapper(rows[3]) is None
    assert [cpdb.resource.complexes.mapper(r) is not None for r in rows] == [
        False, False, False, True, False, False, False
    ]


def test_participant_of_several_uniprot_entries_gives_no_record():
    pair, complex_ = read()[4], read()[5]
    assert cpdb.exclusion_reason(pair) == 'participant with several UniProt entries'
    assert cpdb.exclusion_reason(complex_) == 'participant with several UniProt entries'
    # Neither a pair of the identified protein nor a complex of the others.
    assert cpdb.resource.interactions.mapper(pair) is None
    assert cpdb.resource.complexes.mapper(complex_) is None
    assert not extract('interactions', [pair]).relations
    assert not extract('complexes', [complex_]).relations


def test_gene_identifiers_are_used_only_where_the_protein_has_one_gene():
    (relation,) = extract('interactions', [read()[6]]).relations
    first = cpdb.resource.interactions.mapper(read()[6]).subject
    # Two genes (and Entrez and Ensembl IDs) for CCL4: only the protein's own
    # identifiers and its single gene name and HGNC ID stay.
    assert identifiers(first) == [
        (Namespace.UNIPROT, 'P13236'),
        (Namespace.UNIPROT_ENTRY, 'CCL4_HUMAN'),
        (Namespace.GENESYMBOL, 'CCL4'),
        (Namespace.HGNC, '10630'),
    ]
    # "cccc" and a placeholder number are no PubMed IDs.
    assert values(relation, 'publications') == ['PMID:29477555']


def test_unknown_source_database_raises():
    row = read()[0]
    row = {**row, 'source_databases': ['Reactome', 'NewDB']}
    with pytest.raises(ValueError, match='Unknown ConsensusPathDB source database'):
        cpdb.resource.interactions.mapper(row)
