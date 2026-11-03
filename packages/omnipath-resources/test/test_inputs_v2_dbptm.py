"""dbPTM enzyme-substrate sites: reader, mapper and extractor on source rows."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import dbptm
from pypath.inputs_v2.parsers import dbptm as dbptm_parser

HEADER = '\t'.join(dbptm_parser.COLUMNS) + '\n'
# Rows copied from the rescued table: sites without an enzyme, two enzymes in
# one cell, a group name, autocatalysis, a former symbol, a bacterial enzyme
# with a short window, alternatives and a row without references.
SOURCE = HEADER + (
    'P30443\t\t110\tN\t104\t116\tTLRGYYNQSEDGS\tn-linked glycosylation\tSwiss-Prot\t19159218\n'
    'P63104\tPKA;PKB\t58\tS\t52\t64\tVGARRSSWRVVSS\tphosphorylation\tSwiss-Prot\t11956222;12865427;15883165;16376338\n'
    'P63104\t\t58\tS\t52\t64\tVGARRSSWRVVSS\tphosphorylation\tHPRD\t16376338\n'
    'P31749\tPDPK1\t308\tT\t302\t314\tDGATMKTFCGTPE\tphosphorylation\tSwiss-Prot\t8978681;9512493;15718470;18456494;20333297;20481595\n'
    'P31749\tTNK2\t176\tY\t170\t182\tKATGRYYAMKILK\tphosphorylation\tSwiss-Prot\t20333297\n'
    'Q13315\tQ13315\t367\tS\t361\t373\tTRSLEISQSYTTT\tphosphorylation\tSwiss-Prot\t16858402;17525332\n'
    'Q8TDC3\tLKB1\t189\tT\t183\t195\tGDSLLETSCGSPH\tphosphorylation\tSwiss-Prot\t14976552;18339622\n'
    'Q04206\tPCAF;EP300\t122\tK\t116\t128\tLGIQCVKKRDLEQ\tacetylation\tSwiss-Prot\t12419806\n'
    'P03372\tPRMT1\t260\tR\t254\t266\tGIRKDRRGGRMLK\tmethylation\tSwiss-Prot\t18657504\n'
    'Q07065\tZDHHC2\t100\tC\t94\t106\tASSSASCSRRLGR\tpalmitoylation\tSwiss-Prot\t8314870;19144824\n'
    'P36507\tYersinia yopJ\t222\tV\t216\t228\tSMANSFVGTR\tacetylation\tSwiss-Prot\t-\n'
    'Q92974\tMAPK1 or MAPK3\t679\tT\t673\t685\tGVELLLTPREPAL\tphosphorylation\tSwiss-Prot\t18211802;19690332\n'
    'P08047\tMAPK1 AND MAPK3\t453\tT\t447\t459\tGPIIIRTPTVGPN\tphosphorylation\tSwiss-Prot\t11904305;14593115;16377629\n'
    'P25963\tTyr\t42\tY\t36\t48\tSMKDEEYEQMVKE\tphosphorylation\tSwiss-Prot\t8797825\n'
    'P08195\t\t147\tK\t141\t153\tNGLVKIKVAEDEA\tubiquitylation\tSwiss-Prot\t18781797\n'
)


def parse(text=SOURCE):
    return list(dbptm_parser.iter_enzyme_sites(SimpleNamespace(result=io.StringIO(text))))


def extract(rows):
    out = SilverExtractor('dbptm', 'enzyme_substrate')
    for index, row in enumerate(rows):
        record = dbptm.resource.enzyme_substrate.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def site(out, relation):
    return out.entities[relation.object_entity_key].molecular_form['modifications'][0]


def test_reader_yields_one_row_per_enzyme_and_skips_sites_without_one():
    rows = parse()
    # 12 of the 15 rows name an enzyme; three cells name two.
    assert [(r['substrate'], r['enzyme'], r['position']) for r in rows] == [
        ('P63104', 'PKA', 58),
        ('P63104', 'PKB', 58),
        ('P31749', 'PDPK1', 308),
        ('P31749', 'TNK2', 176),
        ('Q13315', 'Q13315', 367),
        ('Q8TDC3', 'LKB1', 189),
        ('Q04206', 'PCAF', 122),
        ('Q04206', 'EP300', 122),
        ('P03372', 'PRMT1', 260),
        ('Q07065', 'ZDHHC2', 100),
        ('P36507', 'Yersinia yopJ', 222),
        ('Q92974', 'MAPK1 or MAPK3', 679),
        ('P08047', 'MAPK1', 453),
        ('P08047', 'MAPK3', 453),
        ('P25963', 'Tyr', 42),
    ]
    assert rows[0]['kinase'] == 'PKA;PKB'
    assert rows[0]['pmids'] == ['11956222', '12865427', '15883165', '16376338']
    assert (rows[0]['residue'], rows[0]['instance']) == ('S', 'VGARRSSWRVVSS')
    # "-" stands for no reference.
    assert rows[10]['pmids'] == []
    # In a window shorter than 13 the seventh letter is not the site.
    assert (rows[10]['resaa'], rows[10]['residue']) == ('V', None)


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    text = HEADER + (
        'P08195\tPDPK1\t147\tK\t141\t153\tNGLVKIKVAEDEA\tphosphorylation\tSwiss-Prot\t'
        '0123456;123456789;8978681\n'
    )
    (row,) = parse(text)
    assert row['pmids'] == ['8978681']


def test_reader_sorts_each_enzyme_label_into_one_kind():
    kinds = {
        row['enzyme']: [
            kind
            for kind in ('uniprot', 'genesymbol', 'name', 'group')
            if row[f'enzyme_{kind}']
        ]
        for row in parse()
    }
    assert kinds == {
        'PKA': ['group'],
        'PKB': ['group'],
        'PDPK1': ['genesymbol'],
        'TNK2': ['genesymbol'],
        'Q13315': ['uniprot'],
        'LKB1': ['name'],
        'PCAF': ['name'],
        'EP300': ['genesymbol'],
        'PRMT1': ['genesymbol'],
        'ZDHHC2': ['genesymbol'],
        'Yersinia yopJ': ['name'],
        'MAPK1 or MAPK3': [],
        'MAPK1': ['genesymbol'],
        'MAPK3': ['genesymbol'],
        'Tyr': [],
    }
    # The label sets do not overlap.
    sets = (
        dbptm_parser.GENE_SYMBOLS,
        dbptm_parser.GROUP_NAMES,
        dbptm_parser.PROTEIN_NAMES,
        frozenset(dbptm_parser.GENE_SYMBOL_SPELLINGS),
        dbptm_parser.BACTERIAL_NAMES,
        dbptm_parser.NOT_ONE_ENZYME,
        frozenset(dbptm_parser.JOINED_LABELS),
    )
    assert len(frozenset().union(*sets)) == sum(len(labels) for labels in sets)


def test_reader_rejects_changed_layout_and_unknown_values():
    with pytest.raises(ValueError, match='columns'):
        parse('substrate\tenzyme\tresnum\nP31749\tTNK2\t176\n')
    # A table that has the header and no rows is not an empty dataset.
    with pytest.raises(ValueError, match='no rows'):
        parse(HEADER)
    with pytest.raises(ValueError, match='HTML'):
        parse('<html>\n<head><title>404 Not Found</title></head>\n')
    row = 'P31749\t{}\t176\tY\t170\t182\tKATGRYYAMKILK\t{}\t{}\t20333297\n'
    with pytest.raises(ValueError, match='enzyme label'):
        parse(HEADER + row.format('NEWKINASE', 'phosphorylation', 'Swiss-Prot'))
    # A modification type that has no reviewed aspect is not phosphorylation.
    with pytest.raises(ValueError, match='modification type'):
        parse(HEADER + row.format('TNK2', 'ubiquitylation', 'Swiss-Prot'))
    with pytest.raises(ValueError, match='source database'):
        parse(HEADER + row.format('TNK2', 'phosphorylation', 'NewDB'))


def test_label_in_another_letter_case_of_one_gene_symbol_is_that_gene():
    text = HEADER + (
        'Q14155\tCaMK1\t694\tS\t688\t700\tKPERKPSDEEFAS\tphosphorylation\tSwiss-Prot\t18184567;18669648\n'
        'P54646\tLKB1;CaMKK2\t172\tT\t166\t178\tDGEFLRTSCGSPN\tphosphorylation\tSwiss-Prot\t15980064\n'
        'P25963\tCaMK2\t42\tY\t36\t48\tSMKDEEYEQMVKE\tphosphorylation\tSwiss-Prot\t8797825\n'
    )
    rows = parse(text)
    camk1, lkb1, camkk2, camk2 = rows
    assert (camk1['enzyme_genesymbol'], camk1['enzyme_name']) == ('CAMK1', 'CaMK1')
    assert (camkk2['enzyme_genesymbol'], camkk2['enzyme_name']) == ('CAMKK2', 'CaMKK2')
    # A former name stays a name, and CaMK2 (several genes) a group.
    assert (lkb1['enzyme_genesymbol'], lkb1['enzyme_name']) == (None, 'LKB1')
    assert (camk2['enzyme_genesymbol'], camk2['enzyme_group']) == (None, 'CaMK2')

    out = extract(rows)
    subjects = [out.entities[r.subject_entity_key] for r in out.relations]
    assert [(e.entity_type, e.namespace, e.identifier, e.taxon) for e in subjects] == [
        ('protein', Namespace.GENESYMBOL, 'CAMK1', '9606'),
        ('protein', Namespace.NAME, 'LKB1', '9606'),
        ('protein', Namespace.GENESYMBOL, 'CAMKK2', '9606'),
        ('protein_family', Namespace.NAME, 'CaMK2', None),
    ]


def test_kinase_increases_substrate_phosphorylation_at_the_site():
    row = parse()[3]
    record = dbptm.resource.enzyme_substrate.mapper(row)
    assert slot_name(record.predicate) == 'affects'
    assert record.subject.identifiers == [(Namespace.GENESYMBOL, 'TNK2')]
    assert record.object.identifiers == [(Namespace.UNIPROT, 'P31749')]
    assert record.subject.molecular_form is None
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'Y',
            'position': 176,
            'end_position': 176,
            'coordinate_reference': {
                'identifier': {'ns': 'uniprot', 'id': 'P31749'},
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': 'KATGRYYAMKILK',
        }
    ]

    out = extract([row])
    (relation,) = out.relations
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'publications') == ['PMID:20333297']
    assert annotations(relation, 'supporting_data_source') == ['Swiss-Prot']
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}


def test_other_modifications_have_their_own_aspect_and_residue_term():
    rows = parse()
    out = extract([rows[7], rows[8], rows[9]])
    assert [
        (
            annotations(r, 'object_aspect_qualifier'),
            annotations(r, 'object_direction_qualifier'),
            site(out, r)['term'],
            site(out, r)['residue'],
            site(out, r)['position'],
        )
        for r in out.relations
    ] == [
        (['acetylation'], ['increased'], 'acetylated residue', 'K', 122),
        (['methylation'], ['increased'], 'methylated residue', 'R', 260),
        (['palmitoylation'], ['increased'], 'palmitoylated residue', 'C', 100),
    ]


def test_sites_and_enzymes_stay_distinct():
    out = extract(parse()[:4])
    assert len(out.relations) == 4
    assert len({r.relation_key for r in out.relations}) == 4
    # PKA and PKB at one site of P63104; two sites of P31749.
    assert out.relations[0].object_entity_key == out.relations[1].object_entity_key
    assert out.relations[2].object_entity_key != out.relations[3].object_entity_key
    assert [site(out, r)['position'] for r in out.relations] == [58, 58, 308, 176]


def test_enzyme_identifier_follows_the_kind_of_its_label():
    rows = parse()
    out = extract([rows[0], rows[4], rows[5], rows[10]])
    subjects = [out.entities[r.subject_entity_key] for r in out.relations]
    assert [(e.entity_type, e.namespace, e.identifier, e.taxon) for e in subjects] == [
        # A kinase family name is neither a gene symbol nor one protein.
        ('protein_family', Namespace.NAME, 'PKA', None),
        # Autocatalysis: the table has the substrate's accession as the enzyme.
        ('protein', Namespace.UNIPROT, 'Q13315', '9606'),
        # A former name is kept, not replaced by the gene symbol STK11.
        ('protein', Namespace.NAME, 'LKB1', '9606'),
        # A bacterial enzyme is not human.
        ('protein', Namespace.NAME, 'Yersinia yopJ', None),
    ]
    # The short window gives a position but no residue letter.
    short = site(out, out.relations[3])
    assert (short['residue'], short['position'], short['description']) == (
        None,
        222,
        'SMANSFVGTR',
    )
    assert annotations(out.relations[3], 'publications') == []


def test_no_relation_without_one_enzyme_and_no_invented_values():
    rows = parse()
    mapper = dbptm.resource.enzyme_substrate.mapper
    # Alternatives and "Tyr" do not name one enzyme.
    assert mapper(rows[11]) is None
    assert mapper(rows[14]) is None

    out = extract(rows)
    assert len(out.relations) == 13
    for relation in out.relations:
        # A modification is not an effect on activity, and only the row's own
        # modification type is asserted.
        assert annotations(relation, 'object_direction_qualifier') == ['increased']
        (aspect,) = annotations(relation, 'object_aspect_qualifier')
        assert aspect in {'phosphorylation', 'acetylation', 'methylation', 'palmitoylation'}
    # No site of the rows without an enzyme reaches the output.
    assert 'P30443' not in {e.identifier for e in out.entities.values()}
    assert 'P08195' not in {e.identifier for e in out.entities.values()}


def test_modification_crosswalk_covers_the_types_the_parser_accepts():
    assert set(dbptm._MODIFICATIONS) == dbptm_parser.MODIFICATION_TYPES
