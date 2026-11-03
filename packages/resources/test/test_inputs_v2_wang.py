"""Source contract of the Wang Lab Human Signaling Network input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.parsers.wang import iter_network
from pypath.inputs_v2.wang import resource

# Lines of WangLab_HumanSignalingNetwork_v6-format2.csv: its head, the notes,
# a few of its nodes and links.
HEAD = [
    '#By: Combination/Codes/NewSig/Combination/Codes/NewSig/Update_Edwin.pl On: Thu Jan 16 2014  at 10:16:47',
    '',
    'Huamn signaling netwiork, Version 6, Updated: January, 2014',
    '',
    'Citations:',
    '',
    'Cui et al., A map of human cancer signaling, Mol. Syst. Biol., 3:152 , 2007',
    '',
]
NOTES = [
    'Notes:,,,,,',
    'GF,Growth Factor,,,,',
    'NA,Not Available,,,,',
    'TF,Transcription Factor,,,,',
    'Cm,Membrane,,,,',
    'Cy,Intracellular,,,',
    'Nu,Neucleus,,,',
    '',
    '+,activation,,,',
    '_,inhibition,,,',
    '0,physical interaction,,,',
    '',
    'List of nodes and their annotation:,,,,',
    '',
]
NODES = [
    'Node Number in the  network,Gene_Name,Entrez Gene ID,Molecular Function,Cellular Location',
    '1453,EDNRA,1909,Multiple Domain Transmembrane Proteins,Cm',
    '1226,NBN,4683,TF,Cy',
    '2746,MDM2,4193,Adapter,Cy',
    '5080,TP53,7157,TF,Cy',
    '5080,TP53,7157,TF,Nu',
    '6295,GRB2,2885,SP,Cm',
    '6295,GRB2,2885,SP,Cy',
    '3762,SOS1,6654,Adapter,Cm',
    '3762,SOS1,6654,Adapter,Cy',
    '3762,SOS1,6654,Generics,Cy',
    '6104,SRC,6714,Kinase,Cm',
    '4918,ADCY1,107,Cyclase,Cy',
    '1352,CAMP,820,Molecule,Cy',
    '5665,PRKACA,5566,Kinase,Cy',
    '10,CIAPIN1,57019,NA,NA',
    '1317,PKC,50818,Kinase,Cm',
    '1317,PKC,50818,Kinase,Cy',
    '3487,ACTIN,280979,Cytoskeleton,Cy',
    '3533,CYP3A,1574,NA,NA',
    '4661,EGFR,1956,Receptor,Cm',
    '2654,EZR,7430,Generics,Cy',
    '5649,NR1I3,9970,TF,Nu',
    '3429,SLC10A2,6555,NA,NA',
    '3201,RAD23B,5887,NA,NA',
    '',
    '',
]
LINKS = [
    'From (Node Number),To (Node Number),Link_Type,,',
    '1453,1226,+',
    '2746,5080,-',
    '6295,3762,0',
    '6104,6104,0',
    '4918,1352,+',
    '1352,5665,+',
    '1317,4661,-',
    '3487,2654,0',
    '5649,3533,+',
    '3429,3201,0',
    '3201,3429,0',
]


def network(*parts):
    lines = [line for part in parts for line in part]
    # The file has CRLF line ends.
    return SimpleNamespace(result=io.StringIO('\r\n'.join(lines) + '\r\n'))


(
    ACTIVATION,
    INHIBITION,
    PHYSICAL,
    SELF,
    TO_CAMP,
    FROM_CAMP,
    FROM_PKC,
    FROM_ACTIN,
    TO_CYP3A,
    NEUTRAL,
    NEUTRAL_REVERSED,
) = iter_network(network(HEAD, NOTES, NODES, LINKS))


def extract(source_row):
    record = resource.datasets()['interactions'].mapper(source_row)
    assert record is not None
    out = SilverExtractor('wang', 'interactions')
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


def test_parser_joins_a_link_with_its_nodes():
    assert ACTIVATION == {
        'source_node': '1453',
        'source_symbol': 'EDNRA',
        'source_entrez': '1909',
        'source_functions': ['Multiple Domain Transmembrane Proteins'],
        'source_locations': ['Cm'],
        'target_node': '1226',
        'target_symbol': 'NBN',
        'target_entrez': '4683',
        'target_functions': ['TF'],
        'target_locations': ['Cy'],
        'link_type': '+',
        'mirror': False,
    }


def test_parser_collects_the_lines_of_a_node():
    assert PHYSICAL['source_functions'] == ['SP']
    assert PHYSICAL['source_locations'] == ['Cm', 'Cy']
    assert PHYSICAL['target_functions'] == ['Adapter', 'Generics']
    assert PHYSICAL['target_locations'] == ['Cm', 'Cy']
    assert INHIBITION['target_locations'] == ['Cy', 'Nu']
    # The inhibition is read as the links write it.
    assert INHIBITION['link_type'] == '-'


def test_parser_rejects_a_file_without_its_tables():
    with pytest.raises(ValueError, match='no links'):
        list(iter_network(network(HEAD, NOTES, NODES)))
    with pytest.raises(ValueError, match='no links'):
        list(iter_network(network(HEAD, NOTES, NODES, LINKS[:1])))
    with pytest.raises(ValueError, match='no node table'):
        list(iter_network(network(HEAD, NOTES, LINKS)))
    with pytest.raises(ValueError, match='no links'):
        list(iter_network(SimpleNamespace(result=io.StringIO(''))))


def test_parser_rejects_notes_that_define_the_link_types_differently():
    notes = [line.replace('+,activation', '+,inhibition') for line in NOTES]
    with pytest.raises(ValueError, match="link type '\\+'"):
        list(iter_network(network(HEAD, notes, NODES, LINKS)))
    notes = [line for line in NOTES if not line.startswith('0,')]
    with pytest.raises(ValueError, match="link type '0'"):
        list(iter_network(network(HEAD, notes, NODES, LINKS)))


def test_parser_rejects_inconsistent_nodes_and_links():
    with pytest.raises(ValueError, match='unknown node 99'):
        list(iter_network(network(HEAD, NOTES, NODES, [*LINKS, '1453,99,+'])))
    nodes = [*NODES[:2], '1453,EDN1,1906,GF,Ex', *NODES[2:]]
    with pytest.raises(ValueError, match='two genes'):
        list(iter_network(network(HEAD, NOTES, nodes, LINKS)))
    with pytest.raises(ValueError, match='Malformed Wang network link'):
        list(iter_network(network(HEAD, NOTES, NODES, [*LINKS, '1453,1226'])))


def test_activation_is_an_increasing_effect_of_the_source_on_the_target():
    relation, source, target = extract(ACTIVATION)
    assert relation.predicate == 'affects'
    assert identifiers(source) == {
        (Namespace.ENTREZ, '1909'),
        (Namespace.GENESYMBOL, 'EDNRA'),
    }
    assert identifiers(target) == {
        (Namespace.ENTREZ, '4683'),
        (Namespace.GENESYMBOL, 'NBN'),
    }
    assert (source.namespace, source.identifier) == ('entrez', '1909')
    assert {source.entity_type, target.entity_type} == {'protein'}
    assert {source.taxon, target.taxon} == {'9606'}
    assert annotations(relation) == [('object_direction_qualifier', 'increased')]


@pytest.mark.parametrize('link_type', ['-', '_'])
def test_inhibition_is_a_decreasing_effect(link_type):
    relation, source, target = extract(dict(INHIBITION, link_type=link_type))
    assert relation.predicate == 'affects'
    assert (source.identifier, target.identifier) == ('4193', '7157')
    assert annotations(relation) == [('object_direction_qualifier', 'decreased')]


def test_physical_interaction_has_no_direction():
    relation, source, target = extract(PHYSICAL)
    assert relation.predicate == 'interacts_with'
    assert (source.identifier, target.identifier) == ('2885', '6654')
    assert annotations(relation) == []


def test_link_asserts_no_aspect_function_location_or_publication():
    relation, _, _ = extract(INHIBITION)
    assert {a['term'] for a in relation.annotations} == {
        'in_taxon',
        'object_direction_qualifier',
    }
    assert annotations(relation, 'object') == [('in_taxon', 'NCBITaxon:9606')]


def test_link_of_a_node_to_itself_is_kept():
    relation, source, target = extract(SELF)
    assert relation.predicate == 'interacts_with'
    assert source.entity_key == target.entity_key


def test_unknown_link_type_is_rejected():
    mapper = resource.datasets()['interactions'].mapper
    for link_type in ('?', 'Pos', '', None):
        with pytest.raises(ValueError, match='link type'):
            mapper(dict(ACTIVATION, link_type=link_type))


def test_no_relation_for_a_small_molecule_under_a_gene_id():
    # ``CAMP`` with the function ``Molecule`` is cyclic AMP; Entrez Gene 820
    # is the cathelicidin gene.
    mapper = resource.datasets()['interactions'].mapper
    assert TO_CAMP['target_functions'] == ['Molecule']
    assert mapper(TO_CAMP) is None
    assert mapper(FROM_CAMP) is None


@pytest.mark.parametrize('role', ['source', 'target'])
def test_no_relation_without_a_gene_id(role):
    mapper = resource.datasets()['interactions'].mapper
    for missing in ('0', '', 'NA', None):
        assert mapper(dict(ACTIVATION, **{f'{role}_entrez': missing})) is None


def test_no_relation_for_the_names_of_protein_families():
    # ``PKC`` is a family of kinases and ``ACTIN`` of cytoskeletal proteins,
    # ``CYP3A`` a cytochrome P450 subfamily; the file files them under 50818
    # (the dystonia 10 locus), 280979 (the beta-actin gene of the cow) and
    # 1574 (a withdrawn record).
    mapper = resource.datasets()['interactions'].mapper
    assert FROM_PKC['source_functions'] == ['Kinase']
    for row in (FROM_PKC, FROM_ACTIN, TO_CYP3A):
        assert mapper(row) is None
    # Their neighbours are genes and are mapped in their other links.
    assert mapper(dict(FROM_PKC, source_entrez='1956')) is not None


def test_physical_interaction_listed_both_ways_is_one_relation():
    assert (NEUTRAL['source_node'], NEUTRAL['target_node']) == ('3429', '3201')
    assert (NEUTRAL_REVERSED['source_node'], NEUTRAL_REVERSED['target_node']) == (
        '3201',
        '3429',
    )
    assert NEUTRAL['mirror'] is False and NEUTRAL_REVERSED['mirror'] is True
    # The two listings carry the same nodes and the same type.
    assert NEUTRAL['link_type'] == NEUTRAL_REVERSED['link_type'] == '0'
    mapper = resource.datasets()['interactions'].mapper
    assert mapper(NEUTRAL) is not None
    assert mapper(NEUTRAL_REVERSED) is None
    # A node linked to itself is not a mirror; a signed link is never one.
    assert SELF['mirror'] is False and PHYSICAL['mirror'] is False
    assert FROM_PKC['mirror'] is False


def test_signed_links_in_both_directions_are_both_kept():
    links = [*LINKS[:1], '1453,1226,+', '1226,1453,+', '1453,1226,0', '1226,1453,0']
    rows = list(iter_network(network(HEAD, NOTES, NODES, links)))
    assert [r['mirror'] for r in rows] == [False, False, False, True]
    mapper = resource.datasets()['interactions'].mapper
    assert [mapper(r) is None for r in rows] == [False, False, False, True]
