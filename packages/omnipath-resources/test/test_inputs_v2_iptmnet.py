"""iPTMnet enzyme-substrate sites: readers, score join, mapper and extractor on source lines."""

from collections import Counter
import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import iptmnet
from pypath.inputs_v2.parsers import iptmnet as parser

YEAST = "Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)"
# Lines copied from ptm.txt, including one without an enzyme, one whose
# enzyme is only a label, cells with a trailing space and a value that is
# no protein identifier.
PTM = (
    'PHOSPHORYLATION\thprd\tO00401\tWASL\tHomo sapiens (Human)\tY175\tP00519\tABL1\tin vitro;in vivo\t16199863\n'
    'PHOSPHORYLATION\tnpro\tO00401\tWASL\tHomo sapiens (Human)\tY175\tP00519\tABL1\t\t16199863\n'
    'PHOSPHORYLATION\tnpro\tO00401\tWASL\tHomo sapiens (Human)\tY256\tP00519\tABL1\t\t16199863\n'
    'ACETYLATION\tpro\tO15265-1\tATXN7\tHomo sapiens (Human)\tK257\tPR:Q92793\t\t\t19955365\n'
    'ACETYLATION\tpro\tO00716-1\tE2F3\tHomo sapiens (Human)\tK168\tPR:000027108\t\t\t10753885\n'
    f'PHOSPHORYLATION\tpgrd\tD6VTK4\tSTE2\t{YEAST}\tS342\tYNL154C\tYCK2\t\t19828733\n'
    'PHOSPHORYLATION\tnpro\tO43524\tFOXO3\tHomo sapiens (Human)\tS253\tP11309-2\tPIM1\t\t18593906\n'
    'O-GLYCOSYLATION\trlim+\tP34707\tskn-1\tCaenorhabditis elegans\tS470\tO18158\togt-1\t\t28272406\n'
    'UBIQUITINATION\tunip\tQ15386\tUBE3C\tHomo sapiens (Human)\tK903\tQ15386\tUBE3C\t\t32039437\n'
    'SUMOYLATION\tpro\tP42224-1\tSTAT1\tHomo sapiens (Human)\tK703\tPR:000015829\t\t\t12855578\n'
    'METHYLATION\tunip\tP13906\tbspRIM\tLysinibacillus sphaericus (Bacillus sphaericus)\tC181\tP13906\tbspRIM\t\t8065896\n'
    'PHOSPHORYLATION\trlim+\tQ14653\tIRF3\tHomo sapiens (Human)\tS339\tP0C731\tBGLF4\t\t19052084\n'
    'PHOSPHORYLATION\trlim+\tP08057\tTNNI3\tBos taurus (Bovine)\tS25 \tP00516\tPRKG1\t\t2373082\n'
    'PHOSPHORYLATION\trlim+\tP16386\tLIPE\tBos taurus (Bovine)\tS552\tP00517 \t\t\t3345839\n'
    'PHOSPHORYLATION\trlim+\tQ9NB22\t\t\tS189\tQ16513\tPKN2\t\t11112322\n'
    'ACETYLATION\tunip\tO60566\tBUB1B\tHomo sapiens (Human)\tK250\t\tPCAF\t\t19407811\n'
    'PHOSPHORYLATION\tpro\tP35222-1\tCTNNB1\tHomo sapiens (Human)\tT393\tGO:0005956\t\t\t12700239,19565571\n'
    'ACETYLATION\tiedb\tA8TX70\tCOL6A5\tHomo sapiens (Human)\tK652\t\t\tiedb:1797193\t34357683\n'
    'S-NITROSYLATION\tsno\tO75369\tFLNB\tHomo sapiens (Human)\tC2556\t\t\tin vitro\t2212679\n'
)
# Lines copied from score.txt. The table has the isoform-less accession of
# PIM1 and no line for the yeast kinase.
SCORE = (
    'O00401\tY175\tP00519\tPhosphorylation\t3\n'
    'O00401\tY256\tP00519\tPhosphorylation\t4\n'
    'O15265-1\tK257\tPR:Q92793\tAcetylation\t2\n'
    'O00716-1\tK168\tPR:000027108\tAcetylation\t2\n'
    'O43524\tS253\tP11309\tPhosphorylation\t3\n'
    'P34707\tS470\tO18158\tO-Glycosylation\t2\n'
    'Q15386\tK903\tQ15386\tUbiquitination\t2\n'
    'P42224-1\tK703\tPR:000015829\tSumoylation\t2\n'
    'P13906\tC181\tP13906\tMethylation\t2\n'
    'P08057\tS25\tP00516\tPhosphorylation\t2\n'
    'O60566\tK250\t\tAcetylation\t3\n'
    'A8TX70\tK652\t\tAcetylation\t1\n'
    'O75369\tC2556\t\tS-Nitrosylation\t0\n'
)


def opener(text):
    return SimpleNamespace(result=io.StringIO(text))


def parse(ptm=PTM, score=SCORE):
    return list(
        parser.join_scores(
            parser.iter_table(opener(ptm), parser.PTM_COLUMNS),
            parser.iter_table(opener(score), parser.SCORE_COLUMNS),
        )
    )


def extract(rows):
    out = SilverExtractor('iptmnet', 'enzyme_substrate')
    for index, row in enumerate(rows):
        record = iptmnet.resource.enzyme_substrate.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term, scope='relation'):
    return [
        a['value']
        for a in relation.annotations
        if a['term'] == term and a['scope'] == scope
    ]


def score(relation):
    return [
        a['quantity']['has_numeric_value']
        for a in relation.annotations
        if a['term'] == 'has_confidence_score'
    ]


def site(out, relation):
    return out.entities[relation.object_entity_key].molecular_form['modifications'][0]


def one(row):
    out = extract([row])
    (relation,) = out.relations
    return (
        out,
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    row = dict(parse()[0], pmid='0123456,123456789,16199863')
    _, relation, _, _ = one(row)
    assert annotations(relation, 'publications') == ['PMID:16199863']


def test_reader_names_the_columns_of_a_headerless_table():
    rows = list(parser.iter_table(opener(PTM), parser.PTM_COLUMNS))
    assert len(rows) == 19
    assert rows[0] == {
        'ptm_type': 'PHOSPHORYLATION',
        'source': 'hprd',
        'substrate_UniProtAC': 'O00401',
        'substrate_genename': 'WASL',
        'organism': 'Homo sapiens (Human)',
        'site': 'Y175',
        'enzyme_UniProtAC': 'P00519',
        'enzyme_genename': 'ABL1',
        'note': 'in vitro;in vivo',
        'pmid': '16199863',
    }
    # Cells are stripped: "S25 " and "P00517 " in the source.
    assert rows[12]['site'] == 'S25'
    assert rows[13]['enzyme_UniProtAC'] == 'P00517'


def test_reader_rejects_a_changed_layout():
    with pytest.raises(ValueError, match='fields'):
        list(parser.iter_table(opener('O00401\tY175\tP00519\n'), parser.SCORE_COLUMNS))
    with pytest.raises(ValueError, match='fields'):
        list(
            parser.iter_table(
                opener('<html>\n<head><title>403 Forbidden</title></head>\n'),
                parser.PTM_COLUMNS,
            )
        )
    with pytest.raises(ValueError, match='no rows'):
        parse(ptm='')
    with pytest.raises(ValueError, match='two scores'):
        parse(score=SCORE + 'O00401\tY175\tP00519\tPhosphorylation\t1\n')


def test_join_keeps_rows_with_an_enzyme_and_matches_scores_exactly():
    stats = Counter()
    rows = list(
        parser.join_scores(
            parser.iter_table(opener(PTM), parser.PTM_COLUMNS),
            parser.iter_table(opener(SCORE), parser.SCORE_COLUMNS),
            stats,
        )
    )
    # The two sites without an enzyme are not enzyme-substrate rows, and the
    # row whose enzyme is a Gene Ontology term is malformed.
    assert len(rows) == 16
    assert dict(stats) == {
        'rows': 19,
        'sites_without_enzyme': 2,
        'enzyme_rows': 17,
        'malformed_enzyme_accession': 1,
        'enzyme_label_only': 1,
    }
    assert all(row['enzyme_UniProtAC'] or row['enzyme_genename'] for row in rows)
    assert [row['score'] for row in rows[:5]] == ['3', '3', '4', '2', '2']
    by_enzyme = {row['enzyme_UniProtAC']: row for row in rows}
    # No score line names the yeast kinase, and the one for PIM1 has the
    # accession without the isoform: neither row gets a score.
    assert by_enzyme['YNL154C']['score'] is None
    assert by_enzyme['P11309-2']['score'] is None
    # The score of a site without an enzyme accession is not the score of a
    # relation with an enzyme label.
    assert next(r for r in rows if r['enzyme_genename'] == 'PCAF')['score'] is None
    assert next(r for r in rows if r['enzyme_genename'] == 'PRKG1')['score'] == '2'


def test_malformed_enzyme_rows_are_counted_not_silently_skipped(caplog):
    row = 'PHOSPHORYLATION\tpro\t{}\tCTNNB1\tHomo sapiens (Human)\t{}\t{}\t\t\t12700239\n'
    text = (
        row.format('P35222-1', 'T393', 'GO:0005956')
        + row.format('P35222-1', 'T393', 'Q845TI6')
        + row.format('', 'T393', 'P00519')
        + row.format('P35222-1', '', 'P00519')
        + row.format('P35222-1', 'T393', 'P00519')
    )
    stats = Counter()
    with caplog.at_level('WARNING'):
        rows = list(
            parser.join_scores(
                parser.iter_table(opener(text), parser.PTM_COLUMNS),
                parser.iter_table(opener(''), parser.SCORE_COLUMNS),
                stats,
            )
        )
    assert [r['enzyme_UniProtAC'] for r in rows] == ['P00519']
    assert dict(stats) == {
        'rows': 5,
        'enzyme_rows': 5,
        'malformed_enzyme_accession': 2,
        'malformed_substrate_accession': 1,
        'malformed_site': 1,
    }
    assert len(caplog.records) == 4


def test_kinase_increases_substrate_phosphorylation_at_the_site():
    row = parse()[0]
    record = iptmnet.resource.enzyme_substrate.mapper(row)
    assert slot_name(record.predicate) == 'affects'
    assert record.subject.identifiers == [
        (Namespace.UNIPROT, 'P00519'),
        (Namespace.GENESYMBOL, 'ABL1'),
    ]
    assert record.object.identifiers == [
        (Namespace.UNIPROT, 'O00401'),
        (Namespace.GENESYMBOL, 'WASL'),
    ]
    assert record.subject.molecular_form is None
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'Y',
            'position': 175,
            'end_position': 175,
            'coordinate_reference': {
                'identifier': {'ns': 'uniprot', 'id': 'O00401'},
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': None,
        }
    ]

    out, relation, enzyme, substrate = one(row)
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'publications') == ['PMID:16199863']
    assert annotations(relation, 'supporting_data_source') == ['hprd']
    assert score(relation) == [3.0]
    assert (enzyme.entity_type, substrate.entity_type) == ('protein', 'protein')
    # The organism column is the substrate's; the enzyme gets none from it.
    assert (enzyme.taxon, substrate.taxon) == (None, '9606')


def test_sites_stay_distinct_and_source_rows_share_a_statement():
    out = extract(parse()[:3])
    assert len(out.relations) == 3
    # Y175 from HPRD and neXtProt is one statement; Y256 is another.
    keys = [r.relation_key for r in out.relations]
    assert keys[0] == keys[1] != keys[2]
    assert [site(out, r)['position'] for r in out.relations] == [175, 175, 256]
    assert [annotations(r, 'supporting_data_source') for r in out.relations] == [
        ['hprd'],
        ['npro'],
        ['npro'],
    ]
    assert [score(r) for r in out.relations] == [[3.0], [3.0], [4.0]]


def test_ptm_types_map_to_their_own_aspect():
    # The first row of each type in the fixture.
    rows = {row['ptm_type']: row for row in reversed(parse())}
    expected = {
        'PHOSPHORYLATION': ('phosphorylation', 'phosphorylated residue', 'Y'),
        'ACETYLATION': ('acetylation', 'acetylated residue', 'K'),
        'O-GLYCOSYLATION': ('o_linked_glycosylation', 'O-glycosylated residue', 'S'),
        'UBIQUITINATION': ('ubiquitination', 'ubiquitinylated lysine', 'K'),
        'SUMOYLATION': ('sumoylation', 'sumoylated lysine', 'K'),
        'METHYLATION': ('methylation', 'methylated residue', 'C'),
    }
    for ptm_type, (aspect, term, residue) in expected.items():
        out, relation, _, _ = one(rows[ptm_type])
        assert annotations(relation, 'object_aspect_qualifier') == [aspect]
        assert annotations(relation, 'object_direction_qualifier') == ['increased']
        assert (site(out, relation)['term'], site(out, relation)['residue']) == (term, residue)


def test_enzyme_identifier_systems():
    rows = {row['enzyme_UniProtAC']: row for row in parse()}

    # The Protein Ontology names the protein of a UniProt entry by its
    # accession: PR:Q92793 is UniProt Q92793.
    _, relation, enzyme, substrate = one(rows['PR:Q92793'])
    assert (enzyme.entity_type, enzyme.namespace, enzyme.identifier) == (
        'protein',
        Namespace.UNIPROT,
        'Q92793',
    )
    assert {(i['ns'], i['id']) for i in enzyme.identifiers} == {(Namespace.UNIPROT, 'Q92793')}
    assert score(relation) == [2.0]
    # The substrate isoform is kept, and the position refers to it.
    assert substrate.identifier == 'O15265-1'
    assert substrate.molecular_form['isoform_identifier'] == {
        'ns': 'uniprot',
        'id': 'O15265-1',
    }

    # An isoform term of the Protein Ontology is the UniProt isoform.
    _, _, enzyme, _ = one(dict(rows['PR:Q92793'], enzyme_UniProtAC='PR:Q92793-1'))
    assert enzyme.identifier == 'Q92793-1'
    assert enzyme.molecular_form['isoform_identifier'] == {'ns': 'uniprot', 'id': 'Q92793-1'}

    # Another Protein Ontology term can be a class, a family or a complex.
    _, _, enzyme, _ = one(rows['PR:000027108'])
    assert (enzyme.entity_type, enzyme.namespace, enzyme.identifier, enzyme.taxon) == (
        'named_thing',
        'pr',
        'PR:000027108',
        None,
    )

    # A yeast kinase named by its gene's systematic name.
    _, relation, enzyme, substrate = one(rows['YNL154C'])
    assert enzyme.entity_type == 'protein'
    assert {(i['ns'], i['id']) for i in enzyme.identifiers} == {
        ('ordered_locus_name', 'YNL154C'),
        (Namespace.GENESYMBOL, 'YCK2'),
    }
    assert (enzyme.taxon, substrate.taxon) == ('559292', '559292')
    assert score(relation) == []

    # The enzyme isoform is kept as the source gives it.
    _, relation, enzyme, _ = one(rows['P11309-2'])
    assert enzyme.identifier == 'P11309-2'
    assert enzyme.molecular_form['isoform_identifier'] == {'ns': 'uniprot', 'id': 'P11309-2'}
    assert score(relation) == []

    # A viral kinase of a human substrate does not become human.
    _, _, enzyme, substrate = one(rows['P0C731'])
    assert (enzyme.taxon, substrate.taxon) == (None, '9606')


def test_species_is_the_source_organism_or_absent():
    rows = {row['substrate_UniProtAC']: row for row in parse()}
    _, _, _, substrate = one(rows['P34707'])
    assert substrate.taxon == '6239'
    _, _, _, substrate = one(rows['P13906'])
    assert substrate.taxon == '1421'
    # The source leaves the organism of this substrate empty.
    _, _, _, substrate = one(rows['Q9NB22'])
    assert substrate.taxon is None
    # An organism that is not in the crosswalk gets no taxon.
    _, _, _, substrate = one(dict(rows['O00401'], organism='Homo sapiens neanderthalensis'))
    assert substrate.taxon is None


def test_no_relation_without_a_named_enzyme_and_no_invented_values():
    rows = parse()
    mapper = iptmnet.resource.enzyme_substrate.mapper
    # The row whose enzyme is a Gene Ontology term is left out by the reader.
    assert 'GO:0005956' not in {r['enzyme_UniProtAC'] for r in rows}
    assert mapper(dict(rows[0], enzyme_UniProtAC='GO:0005956', enzyme_genename='')) is None
    assert mapper(dict(rows[0], enzyme_UniProtAC='', enzyme_genename='')) is None
    assert mapper(dict(rows[0], substrate_UniProtAC='')) is None

    out = extract(rows)
    assert len(out.relations) == 16
    for relation in out.relations:
        # A modification is not an effect on activity.
        assert annotations(relation, 'object_aspect_qualifier') != ['activity']
        assert annotations(relation, 'object_direction_qualifier') == ['increased']

    # A PTM type without a Biolink aspect is not put under another one.
    assert mapper(dict(rows[0], ptm_type='S-NITROSYLATION')) is None
    assert mapper(dict(rows[0], ptm_type='C-GLYCOSYLATION')) is None
    with pytest.raises(ValueError, match='Unknown iPTMnet PTM type'):
        mapper(dict(rows[0], ptm_type='NEDDYLATION'))
    with pytest.raises(ValueError, match='Unknown iPTMnet source database'):
        mapper(dict(rows[0], source='newdb'))


def test_site_that_is_not_a_letter_and_a_number_keeps_its_text():
    record = iptmnet.resource.enzyme_substrate.mapper(
        dict(parse()[0], site='S147 (or T146')
    )
    (modification,) = record.object.molecular_form['modifications']
    assert modification['description'] == 'S147 (or T146'
    assert (modification['residue'], modification['position']) == (None, None)



# Lines copied from ptm.txt: enzymes given only by a label (human, a dog
# substrate, HPRD with a note), and HPRD rows with each kind of note.
LABELS_AND_NOTES = (
    'PHOSPHORYLATION\tunip\tA6NNZ2\tTUBB8B\tHomo sapiens (Human)\tS172\t\tCDK1\t\t16371510\n'
    'PHOSPHORYLATION\tunip\tO00141\tSGK1\tHomo sapiens (Human)\tT369\t\tPKA\t\t11096081\n'
    'PHOSPHORYLATION\tunip\tP12637\tCASQ2\tCanis lupus familiaris (Dog) (Canis familiaris)\tS397\t\tCK2\t\t1985907\n'
    'PHOSPHORYLATION\thprd\tQ92934\tBAD\tHomo sapiens (Human)\tS75\t\tPAK7\tin vitro;in vivo\t10521512,10837486\n'
    'PHOSPHORYLATION\thprd\tA8K161\tCAMK2A\tHomo sapiens (Human)\tT286\tQ9UQM7\tCAMK2A\tin vitro\t9452388\n'
    'PHOSPHORYLATION\thprd\tO00763-1\tACACB\tHomo sapiens (Human)\tS222\tP54646\tPRKAA2\tin vivo\t12413941\n'
    'PHOSPHORYLATION\thprd\tA8K161\tCAMK2A\tHomo sapiens (Human)\tT305\tQ9UQM7\tCAMK2A\tin vitro;in vivo\t12408851\n'
    'ACETYLATION\tpro\tP04637-1\tTP53\tHomo sapiens (Human)\tK120\tPR:Q9H7Z6\t\t\t17189187\n'
)


def test_enzyme_given_only_by_a_label_is_named_by_the_label():
    stats = Counter()
    rows = list(
        parser.join_scores(
            parser.iter_table(opener(LABELS_AND_NOTES), parser.PTM_COLUMNS),
            parser.iter_table(opener(''), parser.SCORE_COLUMNS),
            stats,
        )
    )
    assert stats['enzyme_label_only'] == 4
    cdk1, pka, ck2, pak7 = rows[:4]

    # PAK7 is only a previous symbol of PAK5 in HGNC, and the row is HPRD's:
    # the label as written, no accession, no gene symbol, no class, no taxon.
    out, relation, enzyme, substrate = one(pak7)
    assert (enzyme.entity_type, enzyme.namespace, enzyme.identifier, enzyme.taxon) == (
        'named_thing',
        Namespace.NAME,
        'PAK7',
        None,
    )
    assert [(i['ns'], i['id']) for i in enzyme.identifiers] == [(Namespace.NAME, 'PAK7')]
    assert (substrate.identifier, substrate.taxon) == ('Q92934', '9606')
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(relation, 'publications') == ['PMID:10521512', 'PMID:10837486']
    assert annotations(relation, 'supporting_data_source') == ['hprd']
    assert score(relation) == []

    # CDK1 of a UniProt annotation of a human substrate (the label is the
    # approved symbol) and the labels PKA and CK2 (families) are others.
    out = extract([cdk1, pka, ck2, pak7])
    assert [out.entities[r.subject_entity_key].identifier for r in out.relations] == [
        'CDK1',
        'PKA',
        'CK2',
        'PAK7',
    ]
    # The substrate's organism is the dog's; the family has no taxon.
    _, _, enzyme, substrate = one(ck2)
    assert (enzyme.entity_type, enzyme.taxon, substrate.taxon) == ('protein_family', None, '9615')


# Lines copied from ptm.txt: enzymes given only by a label, of UniProt
# annotations.
LABEL_ROWS = {
    # Approved symbol, human substrate.
    'CDK1 human': 'PHOSPHORYLATION\tunip\tA6NNZ2\tTUBB8B\tHomo sapiens (Human)\tS172\t\tCDK1\t\t16371510\n',
    # The same label, mouse and yeast substrates.
    'CDK1 mouse': 'PHOSPHORYLATION\tunip\tE9Q8I9\tFry\tMus musculus (Mouse)\tT2516\t\tCDK1\t\t22753416\n',
    'CDK1 yeast': f'PHOSPHORYLATION\tunip\tP38859\tDNA2\t{YEAST}\tS17\t\tCDK1\t\t21841787\n',
    # The approved symbol in another letter case.
    'CaMK1 human': 'PHOSPHORYLATION\tunip\tO43432\tEIF4G3\tHomo sapiens (Human)\tS1156\t\tCaMK1\t\t14507913\n',
    'CaMK1 rat': 'PHOSPHORYLATION\tunip\tO08679\tMark2\tRattus norvegicus (Rat)\tS91\t\tCaMK1\t\t17442826\n',
    'MAPk14 human': 'PHOSPHORYLATION\tunip\tQ9C000\tNLRP1\tHomo sapiens (Human)\tS168\t\tMAPk14\t\t35857590\n',
    # An alias and previous symbol only (of KAT2B), and an approved symbol
    # that is also an alias of another gene (PDK1 of PDPK1).
    'PCAF human': 'ACETYLATION\tunip\tO60566\tBUB1B\tHomo sapiens (Human)\tK250\t\tPCAF\t\t19407811\n',
    'PDK1 human': 'PHOSPHORYLATION\tunip\tP08559\tPDHA1\tHomo sapiens (Human)\tS232\t\tPDK1\t\t11486000\n',
    # Families.
    'PKC human': 'PHOSPHORYLATION\tunip\tO00213\tAPBB1\tHomo sapiens (Human)\tS459\t\tPKC\t\t36250347\n',
    'PKC mouse': 'PHOSPHORYLATION\tunip\tO35717\tSocs2\tMus musculus (Mouse)\tS52\t\tPKC\t\t31578312\n',
    # Not classified: a bacterial toxin, a bacterial protein that matches the
    # symbol HPR in letter case only, a name of one gene that is no symbol,
    # and a class of kinases.
    'TcdB human': 'O-GLYCOSYLATION\tunip\tP08134\tRHOC\tHomo sapiens (Human)\tT37\t\tTcdB\t\t24905543\n',
    'HPr bacillus': 'PHOSPHORYLATION\tunip\tO31644\tmanR\tBacillus subtilis (strain 168)\tH334\t\tHPr\t\t23551403\n',
    'GSK3-beta human': 'PHOSPHORYLATION\tunip\tP04150\tNR3C1\tHomo sapiens (Human)\tS404\t\tGSK3-beta\t\t18838540\n',
    'Tyr-kinases human': 'PHOSPHORYLATION\tunip\tP03372\tESR1\tHomo sapiens (Human)\tY537\t\tTyr-kinases\t\t7539106\n',
}


def label_enzyme(name):
    (row,) = parse(LABEL_ROWS[name], '')
    _, _, enzyme, substrate = one(row)
    return enzyme, substrate


def enzyme_ids(enzyme):
    return [(i['ns'], i['id']) for i in enzyme.identifiers]


def test_label_that_is_a_gene_symbol_of_a_human_substrate_is_that_gene():
    enzyme, substrate = label_enzyme('CDK1 human')
    assert (enzyme.entity_type, enzyme.namespace, enzyme.identifier) == (
        'protein',
        Namespace.GENESYMBOL,
        'CDK1',
    )
    assert enzyme_ids(enzyme) == [(Namespace.GENESYMBOL, 'CDK1'), (Namespace.NAME, 'CDK1')]
    # UniProt's "by CDK1" on a human entry is the human kinase.
    assert (enzyme.taxon, substrate.taxon) == ('9606', '9606')
    # No accession is assigned.
    assert Namespace.UNIPROT not in {ns for ns, _ in enzyme_ids(enzyme)}

    # The row has no score: the score table does not name its enzyme.
    out = extract(parse(LABEL_ROWS['CDK1 human'], ''))
    (relation,) = out.relations
    assert score(relation) == []
    assert annotations(relation, 'supporting_data_source') == ['unip']


def test_same_label_of_a_non_human_substrate_stays_a_name():
    for name, taxon in (('CDK1 mouse', '10090'), ('CDK1 yeast', '559292')):
        enzyme, substrate = label_enzyme(name)
        assert (enzyme.entity_type, enzyme.namespace, enzyme.identifier) == (
            'named_thing',
            Namespace.NAME,
            'CDK1',
        )
        assert enzyme_ids(enzyme) == [(Namespace.NAME, 'CDK1')]
        assert (enzyme.taxon, substrate.taxon) == (None, taxon)
    # The human and the mouse row are not one entity.
    out = extract(parse(LABEL_ROWS['CDK1 human'] + LABEL_ROWS['CDK1 mouse'], ''))
    keys = {r.subject_entity_key for r in out.relations}
    assert len(keys) == 2
    # An HPRD or other source's row is not taken as a UniProt annotation.
    (row,) = parse(LABEL_ROWS['CDK1 human'], '')
    enzyme = one(dict(row, source='hprd'))[2]
    assert (enzyme.entity_type, enzyme.taxon) == ('named_thing', None)


def test_label_that_is_a_gene_symbol_in_another_letter_case_has_the_symbol_and_the_label():
    enzyme, _ = label_enzyme('CaMK1 human')
    assert (enzyme.entity_type, enzyme.taxon) == ('protein', '9606')
    assert enzyme_ids(enzyme) == [(Namespace.GENESYMBOL, 'CAMK1'), (Namespace.NAME, 'CaMK1')]
    enzyme, _ = label_enzyme('MAPk14 human')
    assert enzyme_ids(enzyme) == [(Namespace.GENESYMBOL, 'MAPK14'), (Namespace.NAME, 'MAPk14')]
    # Rat substrate: a name only.
    enzyme, substrate = label_enzyme('CaMK1 rat')
    assert (enzyme.entity_type, enzyme.taxon, substrate.taxon) == ('named_thing', None, '10116')
    assert enzyme_ids(enzyme) == [(Namespace.NAME, 'CaMK1')]


def test_label_that_is_only_an_alias_or_a_previous_symbol_stays_a_name():
    # PCAF is an alias and previous symbol of KAT2B, not an approved symbol.
    enzyme, substrate = label_enzyme('PCAF human')
    assert (enzyme.entity_type, enzyme.taxon, substrate.taxon) == ('named_thing', None, '9606')
    assert enzyme_ids(enzyme) == [(Namespace.NAME, 'PCAF')]
    # PDK1 is an alias of PDPK1 but the approved symbol of the pyruvate
    # dehydrogenase kinase 1, the kinase of PDHA1: the approved symbol decides.
    enzyme, _ = label_enzyme('PDK1 human')
    assert enzyme_ids(enzyme) == [(Namespace.GENESYMBOL, 'PDK1'), (Namespace.NAME, 'PDK1')]
    assert not (set(iptmnet.GENE_SYMBOL_SPELLINGS) & iptmnet.GENE_SYMBOLS)
    assert not {'PCAF', 'LKB1', 'IKKB', 'CDC2', 'PAK7'} & iptmnet.GENE_SYMBOLS


def test_label_that_names_a_family_is_a_protein_family_of_any_organism():
    for name in ('PKC human', 'PKC mouse'):
        enzyme, _ = label_enzyme(name)
        assert (enzyme.entity_type, enzyme.namespace, enzyme.identifier, enzyme.taxon) == (
            'protein_family',
            Namespace.NAME,
            'PKC',
            None,
        )
        assert enzyme_ids(enzyme) == [(Namespace.NAME, 'PKC')]
    # A family label is a name of a group, and is not expanded to genes.
    out = extract(parse(LABEL_ROWS['PKC human'] + LABEL_ROWS['PKC mouse'], ''))
    assert len({r.subject_entity_key for r in out.relations}) == 1
    assert not set(iptmnet.GROUP_NAMES) & (
        iptmnet.GENE_SYMBOLS | set(iptmnet.GENE_SYMBOL_SPELLINGS)
    )


def test_label_in_no_table_is_unchanged():
    for name, label in (
        ('TcdB human', 'TcdB'),
        ('GSK3-beta human', 'GSK3-beta'),
        ('Tyr-kinases human', 'Tyr-kinases'),
        # Matches HPR (haptoglobin-related protein) in letter case only.
        ('HPr bacillus', 'HPr'),
    ):
        enzyme, _ = label_enzyme(name)
        assert (enzyme.entity_type, enzyme.namespace, enzyme.identifier, enzyme.taxon) == (
            'named_thing',
            Namespace.NAME,
            label,
            None,
        )
        assert enzyme_ids(enzyme) == [(Namespace.NAME, label)]


def test_enzyme_with_an_accession_is_not_touched_by_the_label_tables():
    # The gene name of a row with an accession is the source's, as before.
    (row,) = parse(PTM.splitlines(True)[0], '')
    enzyme = one(row)[2]
    assert enzyme_ids(enzyme) == [(Namespace.UNIPROT, 'P00519'), (Namespace.GENESYMBOL, 'ABL1')]
    assert enzyme.taxon is None
    # Even a row of source unip with a human substrate and a group name there.
    enzyme = one(dict(row, enzyme_genename='PKA', source='unip'))[2]
    assert (enzyme.entity_type, enzyme_ids(enzyme)) == (
        'protein',
        [(Namespace.UNIPROT, 'P00519'), (Namespace.GENESYMBOL, 'PKA')],
    )


def test_note_is_the_evidence_type_of_the_row():
    rows = {(r['substrate_genename'], r['site']): r for r in parse(LABELS_AND_NOTES, '')}
    expected = {
        ('TUBB8B', 'S172'): [],
        ('CAMK2A', 'T286'): ['ECO:0000181'],
        ('ACACB', 'S222'): ['ECO:0000178'],
        ('CAMK2A', 'T305'): ['ECO:0000181', 'ECO:0000178'],
        # A label-only row has its note as well.
        ('BAD', 'S75'): ['ECO:0000181', 'ECO:0000178'],
    }
    for key, terms in expected.items():
        _, relation, _, _ = one(rows[key])
        assert annotations(relation, 'has_evidence_of_type') == terms
    mapper = iptmnet.resource.enzyme_substrate.mapper
    with pytest.raises(ValueError, match='Unknown iPTMnet note'):
        mapper(dict(rows['CAMK2A', 'T286'], note='in silico'))
    with pytest.raises(ValueError, match='Unknown iPTMnet note'):
        mapper(dict(rows['CAMK2A', 'T286'], note='in vitro;ex vivo'))


def test_protein_ontology_term_of_a_uniprot_entry_is_that_entry():
    rows = {r['enzyme_UniProtAC']: r for r in parse(LABELS_AND_NOTES, '')}
    _, _, enzyme, _ = one(rows['PR:Q9H7Z6'])
    assert (enzyme.entity_type, enzyme.namespace, enzyme.identifier) == (
        'protein',
        Namespace.UNIPROT,
        'Q9H7Z6',
    )
    # A PRO identifier of the form PR:000012345 is not an accession.
    _, _, enzyme, _ = one(dict(rows['PR:Q9H7Z6'], enzyme_UniProtAC='PR:000027108'))
    assert (enzyme.entity_type, enzyme.namespace, enzyme.identifier) == (
        'named_thing',
        'pr',
        'PR:000027108',
    )
    # The part after PR: must be a UniProt accession.
    assert parser.enzyme_kind('PR:Q9H7Z6') == 'pro_protein'
    assert parser.enzyme_kind('PR:Q9H7Z') is None
    assert parser.enzyme_kind('PR:q9h7z6') is None
