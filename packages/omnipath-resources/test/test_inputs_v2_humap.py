"""hu.MAP 1, 2 and 3 complexes: parser, mapper and extractor on rows of the files."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import humap
from pypath.inputs_v2.parsers.humap import iter_humap1, iter_humap_table

# Lines of the three files (hu.MAP 1 clusters.txt, humap2_complexes_20200809.txt
# and hu.MAP3.0_complexes_wConfidenceScores_total15326_wGenenames_20240922.csv).
V1 = (
    '153129 10670 64121\n'
    '10055 ENSG00000149531 79912 9518\n'
    'GFP_AEQVI 155660 7062\n'
)
V2 = (
    'HuMAP2_ID,Confidence,Uniprot_ACCs,genenames\n'
    'HuMAP2_00000,3,Q9BQS8 O95900,FYCO1 TRUB2\n'
    'HuMAP2_00006,4,P07199 Q96C28 O43296  Q5CZA5,CENPB ZNF707 ZNF264 ZNF678 ZNF805\n'
    'HuMAP2_00117,4,O75385 Q9Y3C0 Q9UM21 O75143 ,ULK1 WASHC3 MGAT4A ATG13 WASH6P\n'
    'HuMAP2_00482,4,Q9P2F5 Q9GZY0 O14640 O14641 Q9GZY0,STOX2 NXF2 DVL1 DVL2 NXF2\n'
    'HuMAP2_02799,3,Q8TAT6 Q92890,NPLOC4 UFD1\n'
    'HuMAP2_05154,1, ,LILRB3 LILRA6\n'
    'HuMAP2_06192,1,Q8TAT6 Q92890,NPLOC4 UFD1\n'
)
V3 = (
    'HuMAP3_ID,ComplexConfidence,Uniprot_ACCs,genenames\n'
    'huMAP3_00000.1,1,P20963 Q9NWV4 Q9UGQ2,CD247 CACFD1 CZIB\n'
    'huMAP3_00075.1,1,Q8WTP9 Q96JM2 Q96KQ7 Q9H2P0 Q9H582 Q9H9B1 Q9UEU5,'
    'ZNF462 EHMT2 ADNP ZNF644 EHMT1 XAGE3 GAGE2D; GAGE8\n'
    'huMAP3_02280.1,1,645345 P62891 Q5VTL8 Q9NQ29,RPL39 LUC7L PRPF38B\n'
    'huMAP3_08127.1,3,645345 O75596,CLEC3A\n'
)


def opener(text):
    return SimpleNamespace(result=io.StringIO(text))


def parse(version, text):
    if version == 1:
        return list(iter_humap1(opener(text)))
    return list(iter_humap_table(opener(text), version=version))


def extract(version, text):
    dataset = f'humap{version}_complexes'
    out = SilverExtractor('humap', dataset)
    for index, row in enumerate(parse(version, text)):
        record = getattr(humap.resource, dataset).mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def identifiers(entity):
    # An entity seen in several rows repeats its identifiers.
    return list(dict.fromkeys((i['ns'], i['id']) for i in entity.identifiers))


def terms(entity, term):
    return [a['value'] for a in entity.annotations if a['term'] == term]


def complexes(out):
    return {
        identifiers(e)[0][1]: e
        for e in out.entities.values()
        if e.entity_type == 'macromolecular_complex'
    }


def members(out, entity):
    return [
        identifiers(out.entities[r.object_entity_key])
        for r in out.relations
        if r.subject_entity_key == entity.entity_key
    ]


def test_version_1_complexes_are_identified_by_their_position_in_the_file():
    rows = parse(1, V1)
    assert [row['complex'] for row in rows] == ['0', '1', '2']
    assert rows[0]['confidence'] is None
    assert [m['entrez'] for m in rows[0]['members']] == ['153129', '10670', '64121']
    assert rows[1]['members'][1] == {
        'uniprot': None, 'entrez': None, 'ensembl': 'ENSG00000149531', 'entry': None,
        'genename': None,
    }
    assert rows[2]['members'][0]['entry'] == 'GFP_AEQVI'


def test_version_1_complex_has_entrez_members_and_no_confidence():
    out = extract(1, V1)
    entity = complexes(out)['0']
    assert identifiers(entity) == [('humap1', '0')]
    assert members(out, entity) == [
        [(Namespace.ENTREZ, '153129')],
        [(Namespace.ENTREZ, '10670')],
        [(Namespace.ENTREZ, '64121')],
    ]
    assert terms(entity, 'in_taxon') == ['NCBITaxon:9606']
    assert terms(entity, 'has_confidence_level') == []
    assert members(out, complexes(out)['1'])[1] == [(Namespace.ENSG, 'ENSG00000149531')]


def test_version_1_jellyfish_protein_has_no_human_taxon():
    out = extract(1, V1)
    (gfp,) = [
        e for e in out.entities.values() if (Namespace.UNIPROT_ENTRY, 'GFP_AEQVI') in identifiers(e)
    ]
    assert gfp.taxon is None
    human = [e for e in out.entities.values() if (Namespace.ENTREZ, '7062') in identifiers(e)]
    assert human[0].taxon == '9606'


def test_version_2_pairs_accessions_with_gene_names_by_position():
    rows = parse(2, V2)
    assert [m['uniprot'] for m in rows[0]['members']] == ['Q9BQS8', 'O95900']
    assert [m['genename'] for m in rows[0]['members']] == ['FYCO1', 'TRUB2']
    # A position without an accession keeps its gene name.
    by_name = {m['genename']: m['uniprot'] for m in rows[1]['members']}
    assert by_name == {
        'CENPB': 'P07199', 'ZNF707': 'Q96C28', 'ZNF264': 'O43296',
        'ZNF678': None, 'ZNF805': 'Q5CZA5',
    }
    assert rows[1]['confidence'] == '4'


def test_version_2_member_listed_twice_is_one_member():
    row = parse(2, V2)[3]
    assert [m['uniprot'] for m in row['members']] == ['Q9P2F5', 'Q9GZY0', 'O14640', 'O14641']


def test_version_2_complex_has_identifier_confidence_and_members():
    out = extract(2, V2)
    entity = complexes(out)['HuMAP2_00000']
    assert identifiers(entity) == [('humap2', 'HuMAP2_00000')]
    assert terms(entity, 'has_confidence_level') == ['High (3)']
    assert members(out, entity) == [
        [(Namespace.UNIPROT, 'Q9BQS8'), (Namespace.GENESYMBOL, 'FYCO1')],
        [(Namespace.UNIPROT, 'O95900'), (Namespace.GENESYMBOL, 'TRUB2')],
    ]
    # No stoichiometry is claimed by listing a protein once.
    assert not [a for r in out.relations for a in r.annotations]


def test_version_2_protein_without_accession_is_identified_by_gene_name():
    out = extract(2, V2)
    assert members(out, complexes(out)['HuMAP2_05154']) == [
        [(Namespace.GENESYMBOL, 'LILRB3')],
        [(Namespace.GENESYMBOL, 'LILRA6')],
    ]
    assert terms(complexes(out)['HuMAP2_05154'], 'has_confidence_level') == ['Extremely High (1)']


def test_complexes_with_the_same_members_stay_distinct_with_their_own_confidence():
    out = extract(2, V2)
    first, second = complexes(out)['HuMAP2_02799'], complexes(out)['HuMAP2_06192']
    assert first.entity_key != second.entity_key
    assert members(out, first) == members(out, second)
    assert terms(first, 'has_confidence_level') == ['High (3)']
    assert terms(second, 'has_confidence_level') == ['Extremely High (1)']
    # The two proteins are shared, not duplicated.
    assert len([e for e in out.entities.values() if e.entity_type == 'protein' and 'NPLOC4' in str(identifiers(e))]) == 1


def test_version_3_members_are_the_accessions_and_gene_names_are_not_read():
    rows = parse(3, V3)
    assert [m['uniprot'] for m in rows[1]['members']] == [
        'Q8WTP9', 'Q96JM2', 'Q96KQ7', 'Q9H2P0', 'Q9H582', 'Q9H9B1', 'Q9UEU5',
    ]
    assert all(m['genename'] is None for row in rows for m in row['members'])
    assert rows[0]['complex'] == 'huMAP3_00000.1'
    assert rows[0]['confidence'] == '1'


def test_version_3_confidence_scale_is_its_own():
    out = extract(3, V3)
    assert terms(complexes(out)['huMAP3_00000.1'], 'has_confidence_level') == [
        'Extremely High (1)'
    ]
    assert terms(complexes(out)['huMAP3_08127.1'], 'has_confidence_level') == ['High (3)']
    # The same number means another level in version 2.
    assert humap._CONFIDENCE[3]['4'] != humap._CONFIDENCE[2]['4']


def test_version_3_number_in_the_accession_column_is_an_entrez_gene_id():
    out = extract(3, V3)
    assert members(out, complexes(out)['huMAP3_08127.1']) == [
        [(Namespace.ENTREZ, '645345')],
        [(Namespace.UNIPROT, 'O75596')],
    ]


def test_versions_are_not_mixed():
    assert set(humap.resource.datasets()) == {
        'humap1_complexes', 'humap2_complexes', 'humap3_complexes',
    }
    urls = {
        name: dataset.download.url
        for name, dataset in humap.resource.datasets().items()
    }
    assert urls['humap1_complexes'].startswith('https://hu1.')
    assert 'humap2' in urls['humap2_complexes']
    assert 'humap3' in urls['humap3_complexes']
    assert complexes(extract(2, V2)).keys().isdisjoint(complexes(extract(3, V3)))


def test_no_threshold_and_no_pairs():
    out = extract(3, V3)
    assert len(complexes(out)) == 4
    assert {r.predicate for r in out.relations} == {'has_member'}


def test_parsers_raise_on_changed_header_empty_file_and_unknown_member():
    with pytest.raises(ValueError, match='lacks columns'):
        parse(2, 'ID,Level,Accessions,names\nx,1,P20963,A\n')
    with pytest.raises(ValueError, match='no complexes'):
        parse(3, 'HuMAP3_ID,ComplexConfidence,Uniprot_ACCs,genenames\n')
    with pytest.raises(ValueError, match='no complexes'):
        parse(1, '')
    with pytest.raises(ValueError, match='unknown member'):
        parse(1, '153129 not-an-id\n')
    with pytest.raises(ValueError, match='gene names'):
        parse(2, 'HuMAP2_ID,Confidence,Uniprot_ACCs,genenames\nHuMAP2_00000,3,Q9BQS8 O95900,FYCO1\n')


def test_unknown_confidence_level_raises_and_complex_without_members_gives_none():
    row = parse(3, V3)[0]
    with pytest.raises(ValueError, match='confidence'):
        humap.map_complex({**row, 'confidence': '9'}, version=3)
    assert humap.map_complex({**row, 'members': []}, version=3) is None
