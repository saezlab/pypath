"""HINT binary and co-complex tables: shared reader, mapper and extractor on source lines."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import hint

HEADER = (
    'Uniprot_A\tUniprot_B\tGene_A\tGene_B\tpmid:method:quality:type\ttaxid\t'
    'high_quality\tin_pdb_source\tPDB_IDs\n'
)
# Lines copied from HomoSapiens_binary_hq.txt, release 2026-07: several
# publications and methods, a PDB entry as the only reference, two PubMed IDs
# run together, a chain of a UniProt entry, and a protein without a gene name.
BINARY = HEADER + (
    'A0A024R2I8\tF1D8Q5\tNR1A2\tNR2B1\t15604093:0018:HT:binary|19211732:0018:LC:binary|'
    '7746322:0018:LC:binary|9171239:0096:LC:binary|9368056:0096:LC:binary|'
    '9415406:0047:LC:binary\t9606\tTrue\tFalse\t\n'
    'A0A075B6N1\tP01848\tTRBV19\tTRAC\tPDB_5NQK:0114:LC:binary\t9606\tTrue\tTrue\t5NQK-B:A\n'
    'P42345\tQ8N122-1\tMTOR\tRPTOR\t12150925:0492:LC:binary|12150925:0493:LC:binary|'
    '12150926:0492:LC:binary|12150926:0493:LC:binary|15066126:0492:LC:binary|'
    '15066126:0493:LC:binary|16354680:0492:LC:binary|16354680:0493:LC:binary|'
    '1683716516603397:0492:LC:binary|1683716516603397:0493:LC:binary\t9606\tTrue\tFalse\t\n'
    'O15303-PRO_0000012934\tO60894\tGRM6\tRAMP1\t39083597:0947:HT:binary\t9606\tTrue\tFalse\t\n'
    'A0A087WZY1\tO60883\t\tGPR37L1\t32296183:0397:HT:binary|32296183:1112:HT:binary|'
    '32296183:1356:HT:binary\t9606\tTrue\tFalse\t\n'
)
# Lines copied from HomoSapiens_cocomp_hq.txt; the first cites a placeholder
# number of one of HINT's sources instead of a PubMed ID.
COCOMPLEX = HEADER + (
    'A0A1W2PQ72\tQ9NQB0\tMSANTD7\tTCF7L2\t888800001068:0004:HT:co-complex\t9606\tTrue\tFalse\t\n'
    'A0A087WWS6\tO94842\tIFNA13\tTOX4\t28514442:0004:HT:co-complex|'
    '33961781:0004:HT:co-complex\t9606\tTrue\tFalse\t\n'
)


def read(text, dataset='human_binary'):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(hint.resource.datasets()[dataset]._raw_parser(opener))


def extract(rows, dataset='human_binary'):
    out = SilverExtractor('hint', dataset)
    for index, row in enumerate(rows):
        record = hint.resource.datasets()[dataset].mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def identifiers(entity):
    return [(i.type, i.value) for i in entity.identifiers]


def values(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_module_import_defines_the_two_human_datasets():
    assert list(hint.resource.datasets()) == ['human_binary', 'human_cocomplex']
    binary, cocomplex = (d.download.url for d in hint.resource.datasets().values())
    assert binary.endswith('/2026-07/HomoSapiens_binary_hq.txt')
    assert cocomplex.endswith('/2026-07/HomoSapiens_cocomp_hq.txt')


def test_shared_reader_yields_the_source_columns():
    rows = read(BINARY)
    assert len(rows) == 5
    assert rows[1] == {
        'Uniprot_A': 'A0A075B6N1',
        'Uniprot_B': 'P01848',
        'Gene_A': 'TRBV19',
        'Gene_B': 'TRAC',
        'pmid:method:quality:type': 'PDB_5NQK:0114:LC:binary',
        'taxid': '9606',
        'high_quality': 'True',
        'in_pdb_source': 'True',
        'PDB_IDs': '5NQK-B:A',
    }


def test_binary_pair_is_a_direct_interaction_with_every_publication_and_method():
    out = extract(read(BINARY)[:1])
    (relation,) = out.relations
    assert relation.predicate == 'directly_physically_interacts_with'
    assert values(relation, 'publications') == [
        'PMID:15604093', 'PMID:19211732', 'PMID:7746322',
        'PMID:9171239', 'PMID:9368056', 'PMID:9415406',
    ]
    assert values(relation, 'has_evidence_of_type') == ['MI:0018', 'MI:0096', 'MI:0047']
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}
    record = hint.resource.human_binary.mapper(read(BINARY)[0])
    assert identifiers(record.subject) == [
        (Namespace.UNIPROT, 'A0A024R2I8'),
        (Namespace.GENESYMBOL, 'NR1A2'),
    ]


def test_cocomplex_pair_is_not_asserted_to_be_direct():
    out = extract(read(COCOMPLEX), 'human_cocomplex')
    assert [r.predicate for r in out.relations] == ['interacts_with'] * 2
    assert values(out.relations[1], 'publications') == ['PMID:28514442', 'PMID:33961781']
    assert values(out.relations[1], 'has_evidence_of_type') == ['MI:0004']


def test_references_that_are_not_pubmed_ids_do_not_become_publications():
    rows = read(BINARY)
    pdb, joined = extract(rows[1:3]).relations
    # A PDB entry in the reference position: the method stays, no PMID is made.
    assert values(pdb, 'publications') == []
    assert values(pdb, 'has_evidence_of_type') == ['MI:0114']
    # Two PubMed IDs the source ran together are not split by guessing.
    assert values(joined, 'publications') == [
        'PMID:12150925', 'PMID:12150926', 'PMID:15066126', 'PMID:16354680',
    ]
    (placeholder,) = extract(read(COCOMPLEX)[:1], 'human_cocomplex').relations
    assert values(placeholder, 'publications') == []
    assert values(placeholder, 'has_evidence_of_type') == ['MI:0004']


def test_isoform_and_chain_are_kept_as_the_form_of_the_participant():
    rows = read(BINARY)
    mapper = hint.resource.human_binary.mapper
    isoform = mapper(rows[2]).object
    assert identifiers(isoform)[0] == (Namespace.UNIPROT, 'Q8N122-1')
    assert isoform.molecular_form['isoform_identifier'] == {'ns': 'uniprot', 'id': 'Q8N122-1'}
    chain = mapper(rows[3]).subject
    assert identifiers(chain) == [
        (Namespace.UNIPROT, 'O15303'),
        (Namespace.GENESYMBOL, 'GRM6'),
    ]
    assert chain.molecular_form['sequence_identifiers'] == [
        {'ns': 'uniprot', 'id': 'O15303-PRO_0000012934'},
    ]
    # An entry accession names no particular sequence.
    assert mapper(rows[0]).subject.molecular_form is None


def test_missing_gene_name_is_not_filled_in():
    record = hint.resource.human_binary.mapper(read(BINARY)[4])
    assert identifiers(record.subject) == [(Namespace.UNIPROT, 'A0A087WZY1')]
    assert identifiers(record.object) == [
        (Namespace.UNIPROT, 'O60883'),
        (Namespace.GENESYMBOL, 'GPR37L1'),
    ]


def test_rows_that_must_not_become_interactions():
    row = read(BINARY)[0]
    mapper = hint.resource.human_binary.mapper
    # A gene name alone, or something that is not a UniProt accession.
    assert mapper(dict(row, Uniprot_B='')) is None
    assert mapper(dict(row, Uniprot_A='NR1A2')) is None
    assert mapper(dict(row, Uniprot_A='ENSP00000493357')) is None
    assert mapper(row) is not None


def test_unknown_evidence_layout_and_taxon_are_left_out():
    row = dict(read(BINARY)[0], taxid='human')
    row['pmid:method:quality:type'] = '15604093|two hybrid|15604093:18:HT:binary'
    record = hint.resource.human_binary.mapper(row)
    assert slot_name(record.predicate) == 'directly_physically_interacts_with'
    assert not record.annotations
    assert not record.subject.annotations


def test_quality_flags_of_the_evidence_items_are_kept_per_relation():
    rows = read(BINARY)
    mixed, pdb, long_list, ht_only, three_methods = extract(rows).relations
    # Real rows: one high-throughput and several literature curated items ...
    assert values(mixed, 'hint:quality') == ['High-throughput', 'Literature curated']
    # ... only literature curated items ...
    assert values(pdb, 'hint:quality') == ['Literature curated']
    assert values(long_list, 'hint:quality') == ['Literature curated']
    # ... and only high-throughput ones; the flag is not repeated for each item.
    assert values(ht_only, 'hint:quality') == ['High-throughput']
    assert values(three_methods, 'hint:quality') == ['High-throughput']
    cocomplex = extract(read(COCOMPLEX), 'human_cocomplex').relations
    assert [values(r, 'hint:quality') for r in cocomplex] == [['High-throughput']] * 2


def test_unknown_quality_flag_raises():
    row = dict(read(BINARY)[0])
    row['pmid:method:quality:type'] += '|9415406:0047:MQ:binary'
    with pytest.raises(ValueError, match="Unknown HINT quality flag: 'MQ'"):
        hint.resource.human_binary.mapper(row)
