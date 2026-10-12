"""ProtMapper phosphosite annotations: table join, mapper and extractor on source rows."""

from collections import Counter
import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import protmapper
from pypath.inputs_v2.parsers import protmapper as parser
from pypath.inputs_v2.parsers.base import iter_csv

# Rows copied from data/export.csv of the version 4 supplement.
ANNOTATIONS = (
    'ID,CTRL_NS,CTRL_ID,CTRL_GENE_NAME,CTRL_IS_KINASE,TARGET_UP_ID,'
    'TARGET_GENE_NAME,TARGET_RES,TARGET_POS,SOURCES,BELIEF\n'
    '2235,UP,P51955,NEK2,True,P51955,NEK2,S,171,"psp,signor",1\n'
    '1704,UP,P51955,NEK2,True,P51955,NEK2,S,241,"psp,signor",1\n'
    '35689,FPLX,GSK3,,True,Q13516,OLIG2,S,10,reach,0.6071839959530149\n'
    '3705,HGNC,33510,HOTAIR,False,P67809,YBX1,S,102,reach,0.6113459482108936\n'
    '484,UP,P05696,Prkca,False,P10636-5,MAPT,S,365,psp,1\n'
    '35375,UP,Q61846,Melk,False,Q9UER7,DAXX,T,578,"reach,sparser",0.7251124837678647\n'
    '13583,UP,Q99683,MAP3K5,True,P45379,TNNT2,S,189,hprd,1\n'
    '35462,FPLX,GSK3,,True,Q13191,CBLB,S,480,reach,0.49164949642864\n'
    '71,UP,Q14289,PTK2B,True,P08631,HCK,Y,411,pid,1\n'
)
# Their rows of data/evidences.csv (the first of two for annotation 3705).
EVIDENCES = (
    'ID,SOURCE,PMID,DBID,TEXT,DESCRIPTION,ORIG_UP_ID,ORIG_RES,ORIG_POS,'
    'MAPPED_UP_ID,MAPPED_RES,MAPPED_POS\n'
    '2235,psp,,http://www.phosphosite.org/phosphosite.owl#Catalysis_12574512_729,,VALID,P51955,S,171,,,\n'
    '2235,signor,17197699,SIGNOR-151755,"Enzymatic activity, induced;",VALID,P51955,S,171,,,\n'
    '1704,psp,,http://www.phosphosite.org/phosphosite.owl#Catalysis_12574500_729,,VALID,P51955,S,241,,,\n'
    '1704,signor,17197699,SIGNOR-151767,"Enzymatic activity, inhibited;",VALID,P51955,S,241,,,\n'
    '35689,reach,28355568,,Olig2 is Phosphorylated by GSK3 at S10.,VALID,Q13516,S,10,,,\n'
    '3705,reach,34266873,,HOTAIR promotes the Ser102 phosphorylation and nuclear translocation of YBX1.,VALID,P67809,S,102,,,\n'
    '484,psp,,http://www.phosphosite.org/phosphosite.owl#Catalysis_9759_2451,,VALID,P10636-5,S,365,,,\n'
    '35375,reach,34365141,,MPK38 phosphorylated DAXX at Thr578.,VALID,Q9UER7,T,578,,,\n'
    '35375,sparser,34365141,,MPK38 phosphorylated DAXX at Thr578.,VALID,Q9UER7,T,578,,,\n'
    '13583,hprd,12819028,http://hprd.org/ptms?hprd_id=01844&isoform_id=01844_4&isoform_name=Isoform_4,,INFERRED_MOUSE_SITE,P45379,S,192,P45379,S,189\n'
    '35462,reach,,,The GSK3 kinase phosphorylates Cbl-b at the S476 and S480 sites to stabilize it.,VALID,Q13191,S,480,,,\n'
    '71,pid,,http://pathwaycommons.org/pc12/Catalysis_05a9c813710014e7db42136d01d0f884,,VALID,P09769,Y,412,,,\n'
)


def archive(annotations=ANNOTATIONS, evidences=EVIDENCES):
    return SimpleNamespace(
        result={
            parser.ANNOTATIONS_FILE: io.StringIO(annotations),
            parser.EVIDENCES_FILE: io.StringIO(evidences),
        }
    )


def parse(annotations=ANNOTATIONS, evidences=EVIDENCES):
    return list(parser.iter_evidences(archive(annotations, evidences)))


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    row = parse()[1]
    assert annotations(extract([row]).relations[0], 'publications') == ['PMID:17197699']
    for number in ('0123456', '123456789'):
        out = extract([dict(row, PMID=number)])
        assert annotations(out.relations[0], 'publications') == []


def by_source(rows, annotation, source):
    return next(r for r in rows if r['ID'] == annotation and r['SOURCE'] == source)


def extract(rows):
    out = SilverExtractor('protmapper', 'enzyme_substrate')
    for index, row in enumerate(rows):
        record = protmapper.resource.enzyme_substrate.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term, scope='relation'):
    return [
        a['value']
        for a in relation.annotations
        if a['term'] == term and a['scope'] == scope
    ]


def site(out, relation):
    return out.entities[relation.object_entity_key].molecular_form['modifications'][0]


def test_parser_yields_one_row_per_evidence_with_its_annotation():
    rows = parse()
    assert len(rows) == 12
    assert rows[1] == {
        'ID': '2235',
        'CTRL_NS': 'UP',
        'CTRL_ID': 'P51955',
        'CTRL_GENE_NAME': 'NEK2',
        'CTRL_IS_KINASE': 'True',
        'TARGET_UP_ID': 'P51955',
        'TARGET_GENE_NAME': 'NEK2',
        'TARGET_RES': 'S',
        'TARGET_POS': '171',
        'SOURCES': 'psp,signor',
        'BELIEF': '1',
        'SOURCE': 'signor',
        'PMID': '17197699',
        'DBID': 'SIGNOR-151755',
        'TEXT': 'Enzymatic activity, induced;',
        'DESCRIPTION': 'VALID',
        'ORIG_UP_ID': 'P51955',
        'ORIG_RES': 'S',
        'ORIG_POS': '171',
        'MAPPED_UP_ID': '',
        'MAPPED_RES': '',
        'MAPPED_POS': '',
    }
    # The rows keep the order of the evidence table.
    assert [row['ID'] for row in rows[:5]] == ['2235', '2235', '1704', '1704', '35689']
    assert [row['SOURCE'] for row in rows if row['ID'] == '35375'] == ['reach', 'sparser']


def test_parser_rejects_a_changed_archive():
    with pytest.raises(ValueError, match='lacks'):
        list(
            parser.iter_evidences(
                SimpleNamespace(result={parser.ANNOTATIONS_FILE: io.StringIO(ANNOTATIONS)})
            )
        )
    with pytest.raises(ValueError, match='lacks columns'):
        parse(annotations=ANNOTATIONS.replace('BELIEF', 'SCORE'))
    with pytest.raises(ValueError, match='lacks columns'):
        parse(evidences=EVIDENCES.replace('SOURCE,', 'DATABASE,'))
    with pytest.raises(ValueError, match='unknown annotation'):
        parse(evidences=EVIDENCES + '99999,psp,,,,VALID,P51955,S,171,,,\n')
    with pytest.raises(ValueError, match='HTML'):
        parse(annotations='<html>\n<head><title>429 Too Many Requests</title></head>\n')
    with pytest.raises(ValueError, match='no rows'):
        parse(evidences=EVIDENCES.splitlines()[0] + '\n')


def test_annotation_without_evidence_is_kept_without_a_source():
    rows = parse(evidences=''.join(EVIDENCES.splitlines(keepends=True)[:-1]))
    assert len(rows) == 12
    assert (rows[-1]['ID'], rows[-1]['SOURCE'], rows[-1]['PMID']) == ('71', None, None)
    (relation,) = extract(rows[-1:]).relations
    assert annotations(relation, 'supporting_data_source') == []
    assert annotations(relation, 'publications') == []
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']


def test_controller_increases_target_phosphorylation_at_the_site():
    row = by_source(parse(), '2235', 'signor')
    record = protmapper.resource.enzyme_substrate.mapper(row)
    assert slot_name(record.predicate) == 'affects'
    assert record.subject.identifiers == [
        (Namespace.UNIPROT, 'P51955'),
        (Namespace.NAME, 'NEK2'),
    ]
    assert record.object.identifiers == [
        (Namespace.UNIPROT, 'P51955'),
        (Namespace.NAME, 'NEK2'),
    ]
    assert record.subject.molecular_form is None
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'S',
            'position': 171,
            'end_position': 171,
            'coordinate_reference': {
                'identifier': {'ns': 'uniprot', 'id': 'P51955'},
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': None,
        }
    ]

    out = extract([row])
    (relation,) = out.relations
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'publications') == ['PMID:17197699']
    assert annotations(relation, 'supporting_data_source') == ['signor']
    (score,) = [
        a['quantity'] for a in relation.annotations if a['term'] == 'has_confidence_score'
    ]
    assert (score['has_numeric_value'], score['source_field']) == (1.0, 'BELIEF')
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    # The target is human; the source gives no species for the controller.
    assert out.entities[relation.object_entity_key].taxon == '9606'
    assert out.entities[relation.subject_entity_key].taxon is None


def test_sites_stay_distinct_and_evidence_rows_share_a_statement():
    rows = [row for row in parse() if row['ID'] in {'2235', '1704'}]
    out = extract(rows)
    assert len(out.relations) == 4
    # S171 from two sources is one statement; S241 is another.
    keys = [r.relation_key for r in out.relations]
    assert keys[0] == keys[1] != keys[2] == keys[3]
    assert [site(out, r)['position'] for r in out.relations] == [171, 171, 241, 241]
    # Each evidence keeps its own source and publication.
    assert [annotations(r, 'supporting_data_source') for r in out.relations] == [
        ['psp'],
        ['signor'],
        ['psp'],
        ['signor'],
    ]
    assert [annotations(r, 'publications') for r in out.relations] == [
        [],
        ['PMID:17197699'],
        [],
        ['PMID:17197699'],
    ]
    assert len(out.payloads) == 4


def test_controller_identifier_systems():
    rows = parse()

    out = extract([by_source(rows, '35689', 'reach')])
    (relation,) = out.relations
    family = out.entities[relation.subject_entity_key]
    # A FamPlex family or complex: no more specific class, no species.
    assert (family.entity_type, family.namespace, family.identifier, family.taxon) == (
        'named_thing',
        'fplx',
        'GSK3',
        None,
    )
    assert annotations(relation, 'supporting_data_source') == ['reach']
    (score,) = [
        a['quantity'] for a in relation.annotations if a['term'] == 'has_confidence_score'
    ]
    assert score['has_numeric_value'] == pytest.approx(0.6071839959530149)

    out = extract([by_source(rows, '3705', 'reach')])
    gene = out.entities[out.relations[0].subject_entity_key]
    assert (gene.entity_type, gene.namespace, gene.identifier, gene.taxon) == (
        'gene',
        Namespace.HGNC,
        '33510',
        '9606',
    )
    assert gene.label == 'HOTAIR'

    with pytest.raises(ValueError, match='Unknown ProtMapper controller namespace'):
        protmapper.resource.enzyme_substrate.mapper(dict(rows[0], CTRL_NS='MESH'))


def test_isoform_target_keeps_its_accession_and_coordinates():
    out = extract([by_source(parse(), '484', 'psp')])
    (relation,) = out.relations
    target = out.entities[relation.object_entity_key]
    assert target.identifier == 'P10636-5'
    assert target.molecular_form['isoform_identifier'] == {'ns': 'uniprot', 'id': 'P10636-5'}
    modification = site(out, relation)
    assert (modification['residue'], modification['position']) == ('S', 365)
    assert modification['coordinate_reference']['identifier'] == {
        'ns': 'uniprot',
        'id': 'P10636-5',
    }
    # A rat kinase: the controller does not get the target's species.
    controller = out.entities[relation.subject_entity_key]
    assert (controller.identifier, controller.label, controller.taxon) == (
        'P05696',
        'Prkca',
        None,
    )


def test_site_is_the_normalized_one_not_the_one_an_evidence_first_reported():
    rows = parse()
    out = extract([by_source(rows, '13583', 'hprd')])
    (mapped,) = out.relations
    # HPRD reported S192 (mouse numbering); the annotation is S189.
    assert (site(out, mapped)['residue'], site(out, mapped)['position']) == ('S', 189)
    assert annotations(mapped, 'publications') == ['PMID:12819028']


# Annotation 67 (CSK, HCK Y411) has the evidence of PhosphoSitePlus and, three
# times, a Pathway Commons (pid) catalysis whose own site is Y412 of FGR
# (P09769). The same pid rows are also the evidence of annotation 29680
# (CSK, FGR Y412). Annotation 187 has a Reactome row at Y577, not Y576, and
# is left with no other evidence; annotation 23020 has a pid row at T399 of
# isoform 2 of its target, which is no site of the canonical sequence given.
MISMATCH_ANNOTATIONS = (
    'ID,CTRL_NS,CTRL_ID,CTRL_GENE_NAME,CTRL_IS_KINASE,TARGET_UP_ID,'
    'TARGET_GENE_NAME,TARGET_RES,TARGET_POS,SOURCES,BELIEF\n'
    '67,UP,P41240,CSK,True,P08631,HCK,Y,411,"pid,psp",1\n'
    '29680,UP,P41240,CSK,True,P09769,FGR,Y,412,pid,1\n'
    '187,UP,P12931,SRC,True,Q05397,PTK2,Y,576,reactome,1\n'
    '23020,UP,P17252,PRKCA,True,P21731,TBXA2R,T,337,pid,1\n'
)
MISMATCH_EVIDENCES = (
    'ID,SOURCE,PMID,DBID,TEXT,DESCRIPTION,ORIG_UP_ID,ORIG_RES,ORIG_POS,'
    'MAPPED_UP_ID,MAPPED_RES,MAPPED_POS\n'
    '67,psp,,http://www.phosphosite.org/phosphosite.owl#Catalysis_3520_1240,,VALID,P08631,Y,411,,,\n'
    '67,pid,11489946,http://pathwaycommons.org/pc12/Catalysis_b9b6e1d9bf57ab1b68c90e4e7855eba3,,VALID,P09769,Y,412,,,\n'
    '67,pid,14967142,http://pathwaycommons.org/pc12/Catalysis_b9b6e1d9bf57ab1b68c90e4e7855eba3,,VALID,P09769,Y,412,,,\n'
    '67,pid,16537894,http://pathwaycommons.org/pc12/Catalysis_b9b6e1d9bf57ab1b68c90e4e7855eba3,,VALID,P09769,Y,412,,,\n'
    '29680,pid,11489946,http://pathwaycommons.org/pc12/Catalysis_b9b6e1d9bf57ab1b68c90e4e7855eba3,,VALID,P09769,Y,412,,,\n'
    '29680,pid,14967142,http://pathwaycommons.org/pc12/Catalysis_b9b6e1d9bf57ab1b68c90e4e7855eba3,,VALID,P09769,Y,412,,,\n'
    '29680,pid,16537894,http://pathwaycommons.org/pc12/Catalysis_b9b6e1d9bf57ab1b68c90e4e7855eba3,,VALID,P09769,Y,412,,,\n'
    '187,reactome,,http://pathwaycommons.org/pc12/Catalysis_d2b42fb082776b086feea15a8b0e744a,,VALID,Q05397,Y,577,,,\n'
    '23020,pid,9054415,http://pathwaycommons.org/pc12/Catalysis_701c923f96237790a27ca744b759d41b,,INFERRED_ALTERNATIVE_ISOFORM,P21731,T,399,P21731-2,T,399\n'
)


def test_evidence_of_another_site_is_not_attached_to_the_annotation():
    stats = Counter()
    rows = list(
        parser.join_evidences(
            iter_csv(SimpleNamespace(result=io.StringIO(MISMATCH_ANNOTATIONS))),
            iter_csv(SimpleNamespace(result=io.StringIO(MISMATCH_EVIDENCES))),
            stats,
        )
    )
    assert dict(stats) == {
        'annotations': 4,
        'evidence_rows': 9,
        'evidence_attached': 4,
        'evidence_dropped': 5,
        'evidence_dropped_protein': 3,
        'evidence_dropped_position': 1,
        'evidence_dropped_isoform': 1,
        'annotations_with_dropped_evidence': 3,
        'annotations_without_evidence': 2,
        'annotations_with_only_dropped_evidence': 2,
    }
    # The rows of the sites the evidence states are kept, and the two
    # annotations left without evidence follow, each once.
    assert [(r['ID'], r['SOURCE']) for r in rows] == [
        ('67', 'psp'),
        ('29680', 'pid'),
        ('29680', 'pid'),
        ('29680', 'pid'),
        ('187', None),
        ('23020', None),
    ]

    out = extract(rows)
    by_site = {}
    for relation in out.relations:
        modification = site(out, relation)
        key = (out.entities[relation.object_entity_key].identifier, modification['position'])
        by_site.setdefault(key, []).append(relation)
    # HCK Y411 keeps the PhosphoSitePlus row and none of the FGR Y412 evidence.
    (hck,) = by_site['P08631', 411]
    assert annotations(hck, 'supporting_data_source') == ['psp']
    assert annotations(hck, 'publications') == []
    # FGR Y412 has it.
    assert [annotations(r, 'publications') for r in by_site['P09769', 412]] == [
        ['PMID:11489946'],
        ['PMID:14967142'],
        ['PMID:16537894'],
    ]
    # An annotation without any evidence at its site is still a statement of
    # the annotation table, without publications or source.
    for key in (('Q05397', 576), ('P21731', 337)):
        (bare,) = by_site[key]
        assert annotations(bare, 'supporting_data_source') == []
        assert annotations(bare, 'publications') == []
        assert annotations(bare, 'object_aspect_qualifier') == ['phosphorylation']
    assert ('Q05397', 577) not in by_site and ('P21731', 399) not in by_site


def test_evidence_site_is_the_mapped_one_when_there_is_one():
    row = {
        'ORIG_UP_ID': 'P45379', 'ORIG_RES': 'S', 'ORIG_POS': '192',
        'MAPPED_UP_ID': 'P45379', 'MAPPED_RES': 'S', 'MAPPED_POS': '189',
    }  # fmt: skip
    assert parser.evidence_site(row) == ('P45379', 'S', '189')
    unmapped = {**row, 'MAPPED_UP_ID': '', 'MAPPED_RES': '', 'MAPPED_POS': ''}
    assert parser.evidence_site(unmapped) == ('P45379', 'S', '192')
    with pytest.raises(ValueError, match='partly empty'):
        parser.evidence_site({**row, 'MAPPED_POS': ''})
    with pytest.raises(ValueError, match='no site'):
        parser.evidence_site({**unmapped, 'ORIG_RES': ''})
    target = {'TARGET_UP_ID': 'P45379', 'TARGET_RES': 'S', 'TARGET_POS': '189'}
    assert parser.site_mismatch(target, row) is None
    assert parser.site_mismatch(target, unmapped) == 'position'
    assert parser.site_mismatch({**target, 'TARGET_RES': 'T'}, row) == 'residue'
    assert parser.site_mismatch({**target, 'TARGET_UP_ID': 'P45379-2'}, row) == 'isoform'
    assert parser.site_mismatch({**target, 'TARGET_UP_ID': 'P08631'}, row) == 'protein'


def test_kinase_flag_and_site_mapping_are_attributes_of_the_relation():
    rows = parse()
    out = extract(
        [
            by_source(rows, '13583', 'hprd'),
            by_source(rows, '3705', 'reach'),
            by_source(rows, '2235', 'signor'),
        ]
    )
    mouse, text_mined, valid = out.relations
    assert annotations(mouse, 'protmapper:ctrl_is_kinase') == ['true']
    assert annotations(mouse, 'protmapper:site_mapping') == ['INFERRED_MOUSE_SITE']
    # HOTAIR is not a kinase; the predicate and qualifiers stay those of the others.
    assert annotations(text_mined, 'protmapper:ctrl_is_kinase') == ['false']
    assert annotations(text_mined, 'protmapper:site_mapping') == ['VALID']
    assert annotations(text_mined, 'object_direction_qualifier') == ['increased']
    assert annotations(valid, 'protmapper:site_mapping') == ['VALID']
    # An annotation without evidence has the kinase flag and no site mapping.
    (bare,) = extract(parse(evidences=EVIDENCES.splitlines()[0] + '\n71,pid,,x,,VALID,P09769,Y,412,,,\n')[-1:]).relations
    assert annotations(bare, 'protmapper:ctrl_is_kinase') == ['true']
    assert annotations(bare, 'protmapper:site_mapping') == []

    mapper = protmapper.resource.enzyme_substrate.mapper
    row = rows[0]
    with pytest.raises(ValueError, match='CTRL_IS_KINASE'):
        mapper(dict(row, CTRL_IS_KINASE='Maybe'))
    with pytest.raises(ValueError, match='site mapping'):
        mapper(dict(row, DESCRIPTION='mapped by hand'))


def test_nothing_is_asserted_beyond_phosphorylation():
    rows = parse()
    out = extract(rows)
    assert len(out.relations) == 12
    for relation in out.relations:
        # Not an effect on activity, whatever the evidence text says
        # ("Enzymatic activity, inhibited;").
        assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
        assert annotations(relation, 'object_direction_qualifier') == ['increased']
        assert not annotations(relation, 'description')
        assert not annotations(relation, 'supporting_text')
    # Text mining evidence without a PubMed ID has no publication.
    (unreferenced,) = [r for r in out.relations if site(out, r)['position'] == 480]
    assert annotations(unreferenced, 'publications') == []
    assert annotations(unreferenced, 'supporting_data_source') == ['reach']

    mapper = protmapper.resource.enzyme_substrate.mapper
    assert mapper(dict(rows[0], TARGET_UP_ID='')) is None
    assert mapper(dict(rows[0], CTRL_ID='', CTRL_GENE_NAME='')) is None
    with pytest.raises(ValueError, match='Unknown ProtMapper evidence source'):
        mapper(dict(rows[0], SOURCE='trips'))


def test_malformed_site_is_rejected():
    mapper = protmapper.resource.enzyme_substrate.mapper
    row = parse()[0]
    with pytest.raises(ValueError, match='Malformed ProtMapper site'):
        mapper(dict(row, TARGET_POS=''))
    with pytest.raises(ValueError, match='Malformed ProtMapper site'):
        mapper(dict(row, TARGET_RES='Ser'))
