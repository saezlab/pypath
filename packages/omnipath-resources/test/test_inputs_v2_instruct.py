"""INstruct: reader, mapper and extractor on lines of sapiens.sin."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import instruct
from pypath.inputs_v2.parsers.instruct import COLUMNS, parse_locations

HEADER = '\t'.join(name for name, _ in COLUMNS) + '\n'
# Lines copied from sapiens.sin: a protein kinase C with itself, a domain that
# occurs three times in one protein and four times in the other, a line that
# cites a consortium instead of a PubMed ID, and one without PDB entries.
KINASE = (
    'Q05513\tQ05513\t5590\t5590\tPRKCZ\tPRKCZ\tPF00130\tPF00069\tC1_1\tPkinase\t(131, 183)\t(252, 518)\t'
    '3pfq\t11078718;15665819;17203073;18650932;19920073;21900206;21911421;9748166\n'
)
REPEATS = (
    'Q12860\tQ12860\t1272\t1272\tCNTN1\tCNTN1\tPF00041\tPF07679\tfn3\tI-set\t'
    '(810, 893);(605, 693);(708, 796)\t(242, 327);(331, 408);(505, 600);(413, 501)\t'
    '2cqv;2jlk;2jll;2nzi;3mtr\t12700241;7628014\n'
)
CONSORTIUM = (
    'P29317\tP21709\t1969\t2041\tEPHA2\tEPHA1\tPF00536\tPF00536\tSAM_1\tSAM_1\t(902, 966)\t(911, 975)\t'
    '1b4f;1pk1;1pk3;2kiv;2qkq;3h8m;3hil;3kka\tStructural Genomics Consortium\n'
)
NO_PDB = (
    'P25054\tQ96RT1\t324\t55914\tAPC\tERBB2IP\tPF05937\tPF00595\tEB1_binding\tPDZ\t(2670, 2843)\t(1321, 1407)\t'
    '\t15123239;20936779\n'
)
TEXT = HEADER + KINASE + REPEATS + CONSORTIUM + NO_PDB


def read(text=TEXT):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(instruct.resource.interactions._raw_parser(opener))


def extract(rows):
    out = SilverExtractor('instruct', 'interactions')
    for index, row in enumerate(rows):
        record = instruct.resource.interactions.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def values(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def identifiers(entity):
    return [(i.type, i.value) for i in entity.identifiers]


def test_module_import_defines_the_dataset_without_downloading():
    assert list(instruct.resource.datasets()) == ['interactions']
    # The original address is gone; the copy on the rescued-data server is read.
    assert instruct.resource.interactions.download.url == (
        'https://rescued.omnipathdb.org/instruct-sapiens.sin'
    )


def test_locations_are_read_as_pairs_of_residue_numbers():
    assert parse_locations('(131, 183)') == [(131, 183)]
    assert parse_locations('(810, 893);(605, 693)') == [(810, 893), (605, 693)]
    for cell in ('131-183', '(0, 5)', '(9, 3)', ''):
        with pytest.raises(ValueError, match='unexpected domain location'):
            parse_locations(cell)


def test_reader_yields_lists_for_locations_entries_and_publications():
    rows = read()
    assert len(rows) == 4
    assert rows[1]['locations_a'] == [(810, 893), (605, 693), (708, 796)]
    assert rows[1]['pdb_ids'] == ['2cqv', '2jlk', '2jll', '2nzi', '3mtr']
    assert rows[2]['publications'] == ['Structural Genomics Consortium']
    assert rows[3]['pdb_ids'] == []
    assert rows[0]['pfam_b'] == 'PF00069' and rows[0]['domain_b'] == 'Pkinase'


def test_reader_rejects_a_table_that_is_not_the_expected_one():
    with pytest.raises(ValueError, match='lacks expected columns'):
        read('ProtA[uniprot]\tProtB[uniprot]\nQ05513\tQ05513\n')
    with pytest.raises(ValueError, match='no rows'):
        read(HEADER)
    with pytest.raises(ValueError, match='not a Pfam accession'):
        read(HEADER + KINASE.replace('PF00130', 'C1_1'))
    with pytest.raises(ValueError, match='not a PDB entry'):
        read(HEADER + KINASE.replace('3pfq', 'pdb'))


def test_pair_is_a_plain_interaction_with_its_pubmed_ids_and_pdb_entries():
    out = extract(read()[:1])
    (relation,) = out.relations
    # The file states no interaction type, so only the broad predicate holds.
    assert relation.predicate == 'interacts_with'
    assert values(relation, 'publications') == [
        'PMID:11078718', 'PMID:15665819', 'PMID:17203073', 'PMID:18650932',
        'PMID:19920073', 'PMID:21900206', 'PMID:21911421', 'PMID:9748166',
    ]
    assert values(relation, 'instruct:pdb_id') == ['3pfq']
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}


def test_participant_has_uniprot_entrez_and_symbol():
    record = instruct.resource.interactions.mapper(read()[2])
    assert identifiers(record.subject) == [
        (Namespace.UNIPROT, 'P29317'),
        (Namespace.ENTREZ, '1969'),
        (Namespace.GENESYMBOL, 'EPHA2'),
    ]


def test_domain_is_a_region_of_each_protein_with_its_pfam_accession():
    record = instruct.resource.interactions.mapper(read()[0])
    assert record.subject.molecular_form['regions'] == [
        {
            'type': 'interacting domain',
            'identifier': {'ns': 'pfam', 'id': 'PF00130'},
            'position': 131,
            'end_position': 183,
            'coordinate_reference': {
                'identifier': None,
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': 'C1_1',
        }
    ]
    (kinase,) = record.object.molecular_form['regions']
    assert kinase['identifier'] == {'ns': 'pfam', 'id': 'PF00069'}
    assert (kinase['position'], kinase['end_position']) == (252, 518)


def test_every_location_of_a_repeated_domain_is_kept():
    record = instruct.resource.interactions.mapper(read()[1])
    assert [
        (r['position'], r['end_position']) for r in record.subject.molecular_form['regions']
    ] == [(810, 893), (605, 693), (708, 796)]
    assert len(record.object.molecular_form['regions']) == 4


def test_text_that_is_not_a_pubmed_id_is_not_a_publication():
    (relation,) = extract(read()[2:3]).relations
    assert values(relation, 'publications') == []
    assert len(values(relation, 'instruct:pdb_id')) == 8


def test_row_without_pdb_entries_still_gives_the_interaction():
    (relation,) = extract(read()[3:4]).relations
    assert values(relation, 'instruct:pdb_id') == []
    assert values(relation, 'publications') == ['PMID:15123239', 'PMID:20936779']
    terms = {a['term'] for a in relation.annotations}
    assert 'has_evidence_of_type' not in terms
