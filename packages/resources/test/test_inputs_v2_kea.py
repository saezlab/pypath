"""KEA kinase-substrate sites: parser, mapper and extractor on source lines."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import kea
from pypath.inputs_v2.parsers.kea import iter_kea

# Lines copied from the source file, including its one malformed line, a
# kinase with a trailing space and references that are not PubMed IDs.
SOURCE = (
    'ABL1\tABL1_Y226\t10964922\tphosphosite\n'
    'ABL1\tABL1_Y226\t12522270\thprd\n'
    'ABL1\tABL1_Y245\t10964922\tkinexus\n'
    'ABL1\tACAP2_Y742\t17570479\tnetworkin\n'
    'CSF1R\tCSF1R_Y809\tN.N.\tphosphoelm\n'
    'MAPK1\tPTPRR_Y44\tCST-PHOSPHOSITE-5539963\tkinexus\n'
    'PRKACA\tSOX9_S181\t\tkinexus\n'
    'EGFR \tGAB1_Y285\t9890893\tkinexus\n'
    'PRKAR1A\tHLA-A_S336\t6333425\thprd\n'
    'NME2\tNME2_H118\t1851158\thprd\n'
    'NA\tMAP2K4  AKT1_IERLRTHSIESSGK1L\t\tkinexus\n'
)


def parse(text=SOURCE):
    return list(iter_kea(SimpleNamespace(result=io.StringIO(text))))


def extract(rows):
    out = SilverExtractor('kea', 'enzyme_substrate')
    for index, row in enumerate(rows):
        record = kea.resource.enzyme_substrate.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_parser_reads_site_reference_and_source_database():
    rows = parse()
    assert len(rows) == 11
    assert rows[0] == {
        'enzyme': 'ABL1',
        'substrate_site': 'ABL1_Y226',
        'reference': '10964922',
        'source_database': 'phosphosite',
        'substrate': 'ABL1',
        'residue': 'Y',
        'position': 226,
        'pmid': '10964922',
    }
    # A kinase with a trailing space and a hyphenated gene symbol.
    assert (rows[7]['enzyme'], rows[7]['substrate']) == ('EGFR', 'GAB1')
    assert (rows[8]['substrate'], rows[8]['residue'], rows[8]['position']) == (
        'HLA-A',
        'S',
        336,
    )
    # Other reference labels stay in the row but are not PubMed IDs.
    assert [row['pmid'] for row in rows[4:7]] == [None, None, None]
    assert rows[5]['reference'] == 'CST-PHOSPHOSITE-5539963'


def test_parser_rejects_changed_layout():
    with pytest.raises(ValueError, match='Line 1 has 3 fields, expected 4'):
        parse('ABL1\tABL1_Y226\t10964922\n')
    with pytest.raises(ValueError, match='Line 1 has 1 fields, expected 4'):
        parse('<html>\n<head><title>404 Not Found</title></head>\n')


def test_kinase_increases_substrate_phosphorylation_at_the_site():
    record = kea.resource.enzyme_substrate.mapper(parse()[0])
    assert slot_name(record.predicate) == 'affects'
    assert record.subject.identifiers == [(Namespace.GENESYMBOL, 'ABL1')]
    assert record.object.identifiers == [(Namespace.GENESYMBOL, 'ABL1')]
    assert record.subject.molecular_form is None
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'Y',
            'position': 226,
            'end_position': 226,
            'coordinate_reference': {
                'identifier': None,
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': 'ABL1_Y226',
        }
    ]

    out = extract(parse()[:1])
    (relation,) = out.relations
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'publications') == ['PMID:10964922']
    assert annotations(relation, 'supporting_data_source') == ['phosphosite']
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}


def test_sites_stay_distinct_and_evidence_lines_share_a_statement():
    out = extract(parse()[:3])
    assert len(out.relations) == 3
    # Y226 from two databases is one statement; Y245 is another.
    keys = [r.relation_key for r in out.relations]
    assert keys[0] == keys[1] != keys[2]
    sites = {
        out.entities[r.object_entity_key].molecular_form['modifications'][0]['position']
        for r in out.relations
    }
    assert sites == {226, 245}
    assert [annotations(r, 'supporting_data_source')[0] for r in out.relations] == [
        'phosphosite',
        'hprd',
        'kinexus',
    ]


def test_source_residue_letter_is_kept_as_reported():
    out = extract([parse()[9]])
    site = next(
        e.molecular_form for e in out.entities.values() if e.molecular_form
    )['modifications'][0]
    assert (site['residue'], site['position']) == ('H', 118)


def test_no_relation_without_kinase_or_site_and_no_invented_values():
    rows = parse()
    assert kea.resource.enzyme_substrate.mapper(rows[10]) is None
    assert kea.resource.enzyme_substrate.mapper(dict(rows[0], enzyme=None)) is None

    out = extract(rows)
    assert len(out.relations) == 10
    for relation in out.relations:
        # Phosphorylation is not an effect on activity.
        assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
        assert annotations(relation, 'object_direction_qualifier') == ['increased']
    # "N.N.", a CST label and an empty reference are not publications.
    assert [annotations(r, 'publications') for r in out.relations[4:7]] == [[], [], []]


def test_unknown_source_database_is_rejected():
    with pytest.raises(ValueError, match='Unknown KEA source database'):
        kea.resource.enzyme_substrate.mapper(
            dict(parse()[0], source_database='newdb')
        )


def test_only_well_formed_pubmed_ids_become_publications():
    # Rows copied from the source: "0" (kinexus), a 6-digit and a 9-digit
    # number are not PubMed IDs of 7 or 8 digits.
    rows = parse(
        'ATM\tRAD17_S635\t0\tkinexus\n'
        'LCK\tVAV1_Y174\t898671\tphosphoelm\n'
        'MAPK14\tATF2_T71\t602896\tphosphosite\n'
        'HIPK2\tRUNX1_S249\t169175071\thprd\n'
        'ABL1\tABL1_Y226\t10964922\tphosphosite\n'
        'ABL1\tACAP2_Y742\t1234567\tnetworkin\n'
        'ABL1\tACAP2_Y743\t01234567\tnetworkin\n'
    )
    assert [row['pmid'] for row in rows] == [
        None,
        None,
        None,
        None,
        '10964922',
        '1234567',
        None,
    ]
    # The reference stays in the parsed row as written.
    assert rows[3]['reference'] == '169175071'
    out = extract(rows)
    # The site relation remains; it has no publication.
    assert len(out.relations) == 7
    assert [annotations(r, 'publications') for r in out.relations] == [
        [],
        [],
        [],
        [],
        ['PMID:10964922'],
        ['PMID:1234567'],
        [],
    ]
