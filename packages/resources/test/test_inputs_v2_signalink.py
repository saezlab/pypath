"""SignaLink 3: parsers, mappers and extractor on source rows."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import signalink
from pypath.inputs_v2.parsers import signalink as parser

NODE_HEADER = 'id\tname\talt_accession\ttax_id\tpathways\taliases\ttopology\n'
# Rows copied from the nodes file.
NODES = [
    (125, 'Uniprot:F8VW24', '9606', 'Toll-like receptor', ''),
    (1154, 'Uniprot:H0YK56', '9606', '', ''),
    (6121, 'Uniprot:O60674', '9606', 'JAK/STAT', 'Mediator'),
    (6577, 'Uniprot:P13528', '6239', 'WNT/Wingless', 'Transcription factor'),
    (6536, 'Uniprot:P17486', '6239', 'WNT/Wingless', 'Transcription factor'),
    (69, 'Uniprot:Q06124', '9606', 'JAK/STAT|Receptor tyrosine kinase|T-cell receptor', 'Co-factor|Scaffold'),
    (895, 'Uniprot:Q15759', '9606', 'Receptor tyrosine kinase|Toll-like receptor', 'Endocytosis related|Mediator'),
    (6937, 'Uniprot:Q7SX76', '7955', '', ''),
    (2091, 'Uniprot:Q8IW41', '9606', 'Receptor tyrosine kinase', 'Mediator'),
    (6938, 'Uniprot:Q8JG38', '7955', 'RTK', ''),
    (6641, 'Uniprot:Q9BZR9', '9606', 'JAK/STAT', 'Co-factor'),
    (1339, 'Uniprot:Q9Y6X2', '9606', 'JAK/STAT|TGF', 'Co-factor'),
]
NODES_TEXT = NODE_HEADER + ''.join(
    f'{i}\t{name}\t\ttaxid:{tax}\t{pathways}\t\t{functions}\n'
    for i, name, tax, pathways, functions in NODES
)
EDGE_HEADER = (
    'interactor_a_node_name\tinteractor_b_node_name\tinteraction_detection_method\t'
    'first_author\tpublication_ids\tinteraction_types\tsource_db\t'
    'interaction_identifiers\tconfidence_scores\tlayer\n'
)
# Rows copied from the layer 0 file.
EDGES = [
    ('Uniprot:H0YK56', 'Uniprot:F8VW24', 'MI:0007(anti tag coimmunoprecipitation)',
     'pubmed:12860405|pubmed:23180781',
     'is_direct:false|is_directed:true|MI:0915(physical association)', 'InnateDB'),
    ('Uniprot:Q15759', 'Uniprot:Q8IW41', 'MI:0423(in-gel kinase assay)',
     'pubmed:23180781|pubmed:23331499|pubmed:9628874',
     'is_direct:false|is_direct:true|is_directed:true|MI:0217(phosphorylation reaction)|MI:0624(stimulation)',
     'InnateDB|SLKv2.1'),
    ('Uniprot:O60674', 'Uniprot:Q06124', '',
     'pubmed:21071413|pubmed:23331499|pubmed:8995399',
     'is_direct:true|is_directed:true|MI:0624(stimulation)', 'SLKv2.0|SLKv2.1'),
    ('Uniprot:P17486', 'Uniprot:P13528', '', 'pubmed:10049362|pubmed:23331499',
     'is_direct:false|is_directed:true|MI:0623(inhibition)|MI:0624(stimulation)',
     'SLKv2.0|SLKv2.1'),
    ('Uniprot:Q9BZR9', 'Uniprot:Q9Y6X2', '', 'pubmed:20516148',
     'is_direct:true|is_directed:true|MI:0220(ubiquitination reaction)|MI:0623(inhibition)',
     'SLKv3.0'),
    ('Uniprot:Q8JG38', 'Uniprot:Q7SX76', '', 'pubmed:20802085',
     'is_direct:true|is_directed:true', 'SignaFish'),
]


def edges_text(edges=EDGES):
    return EDGE_HEADER + ''.join(
        f'{a}\t{b}\t{method}\t\t{pubs}\t{types}\t{sources}\t\t\t0\n'
        for a, b, method, pubs, types, sources in edges
    )


class Source:
    def __init__(self, text):
        self.text = text

    def open(self, **_kwargs):
        return SimpleNamespace(result=self.text)


def interactions(edges=EDGES, nodes=NODES_TEXT):
    return list(
        parser.iter_interactions(
            Source(edges_text(edges)).open(), nodes_download=Source(nodes)
        )
    )


def nodes():
    return list(parser.iter_nodes(Source(NODES_TEXT).open()))


def extract(dataset, rows):
    out = SilverExtractor('signalink', dataset)
    mapper = signalink.resource.datasets()[dataset].mapper
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def relation_attributes(relation):
    result = {}
    for a in relation.annotations:
        if a['scope'] == 'relation':
            result.setdefault(a['term'], []).append(a['value'])
    return result


def test_nodes_parser_keeps_the_annotated_nodes():
    rows = nodes()
    assert [r['accession'] for r in rows] == [
        'F8VW24', 'O60674', 'P13528', 'P17486', 'Q06124', 'Q15759', 'Q8IW41',
        'Q8JG38', 'Q9BZR9', 'Q9Y6X2',
    ]
    q06124 = rows[4]
    assert q06124 == {
        'name': 'Uniprot:Q06124',
        'database': 'Uniprot',
        'accession': 'Q06124',
        'taxon': '9606',
        'pathways': ['JAK/STAT', 'Receptor tyrosine kinase', 'T-cell receptor'],
        'functions': ['Co-factor', 'Scaffold'],
    }


def test_terms_are_the_pathways_and_functions_as_written():
    rows = list(parser.iter_terms(Source(NODES_TEXT).open()))
    kinds = {}
    for row in rows:
        kinds.setdefault(row['kind'], []).append(row['term'])
    # The source writes "WNT/Wingless" and the pathway "RTK" next to "Receptor
    # tyrosine kinase"; both are kept as they are.
    assert 'WNT/Wingless' in kinds['pathway'] and 'RTK' in kinds['pathway']
    assert kinds['function'] == [
        'Mediator', 'Transcription factor', 'Co-factor', 'Scaffold', 'Endocytosis related'
    ]
    out = extract('terms', rows)
    assert {e.entity_type for e in out.entities.values()} == {'ontology_class'}
    assert {e.namespace for e in out.entities.values()} == {
        'signalink_pathway',
        'signalink_function',
    }
    assert not out.relations


def test_annotation_associates_a_protein_with_its_pathways_and_functions():
    rows = nodes()
    out = extract('annotations', rows)
    assert {r.predicate for r in out.relations} == {'associated_with'}
    q06124 = next(k for k, e in out.entities.items() if e.identifier == 'Q06124')
    objects = {
        (out.entities[r.object_entity_key].namespace, out.entities[r.object_entity_key].identifier)
        for r in out.relations
        if r.subject_entity_key == q06124
    }
    assert objects == {
        ('signalink_pathway', 'JAK/STAT'),
        ('signalink_pathway', 'Receptor tyrosine kinase'),
        ('signalink_pathway', 'T-cell receptor'),
        ('signalink_function', 'Co-factor'),
        ('signalink_function', 'Scaffold'),
    }
    proteins = [e for e in out.entities.values() if e.entity_type == 'protein']
    assert {e.taxon for e in proteins} == {'9606', '6239', '7955'}
    # A node with neither pathway nor function gives no statement.
    assert not any(e.identifier in ('H0YK56', 'Q7SX76') for e in out.entities.values())


def test_objects_are_the_same_entities_as_the_terms():
    terms = extract('terms', list(parser.iter_terms(Source(NODES_TEXT).open())))
    annotations = extract('annotations', nodes())
    identity = lambda out, keys: {  # noqa: E731
        (out.entities[k].entity_type, out.entities[k].namespace, out.entities[k].identifier)
        for k in keys
    }
    objects = identity(annotations, {r.object_entity_key for r in annotations.relations})
    assert objects == identity(terms, terms.entities)


def test_signed_pair_is_an_effect_of_the_first_node_on_the_second():
    rows = interactions()
    stim = next(r for r in rows if r['a']['accession'] == 'O60674')
    record = signalink.map_interaction(stim)
    assert slot_name(record.predicate) == 'affects'
    assert record.subject.identifiers == [(Namespace.UNIPROT, 'O60674')]
    assert record.object.identifiers == [(Namespace.UNIPROT, 'Q06124')]
    out = extract('interactions', [stim])
    attrs = relation_attributes(out.relations[0])
    assert attrs['object_direction_qualifier'] == ['increased']
    assert 'object_aspect_qualifier' not in attrs
    assert attrs['publications'] == ['PMID:21071413', 'PMID:23331499', 'PMID:8995399']
    assert attrs['supporting_data_source'] == ['SLKv2.0', 'SLKv2.1']
    assert attrs['original_predicate'] == ['MI:0624']
    assert attrs['signalink:is_direct'] == ['true']


def test_pair_with_both_signs_gives_one_relation_for_each():
    rows = [r for r in interactions() if r['a']['accession'] == 'P17486']
    assert sorted(r['direction'] for r in rows) == ['decreased', 'increased']
    out = extract('interactions', rows)
    assert sorted(
        relation_attributes(r)['object_direction_qualifier'][0] for r in out.relations
    ) == ['decreased', 'increased']
    # Both are effects of the same pair; neither is a direct statement.
    assert {relation_attributes(r)['signalink:is_direct'][0] for r in out.relations} == {'false'}


def test_enzymatic_term_gives_the_aspect_and_its_own_direction_only_where_it_agrees():
    rows = interactions()
    phos = next(r for r in rows if r['a']['accession'] == 'Q15759')
    assert (phos['direction'], phos['aspect']) == ('increased', 'phosphorylation')
    # "is_direct:false|is_direct:true": one of the merged statements is direct.
    assert phos['is_direct'] is True
    ubiq = [r for r in rows if r['a']['accession'] == 'Q9BZR9']
    assert sorted((r['direction'], r['aspect']) for r in ubiq) == [
        ('decreased', None),
        ('increased', 'ubiquitination'),
    ]


def test_pair_without_an_effect_is_an_interaction_and_without_a_type_an_association():
    rows = interactions()
    phys = next(r for r in rows if r['a']['accession'] == 'H0YK56')
    assert not phys['affects']
    record = signalink.map_interaction(phys)
    # Not direct and a physical association: it says only that they interact.
    assert slot_name(record.predicate) == 'interacts_with'
    assert 'object_direction_qualifier' not in relation_attributes(
        extract('interactions', [phys]).relations[0]
    )
    untyped = next(r for r in rows if r['a']['accession'] == 'Q8JG38')
    # No type, but the row says the pair is direct (the SignaFish rows): a direct
    # interaction, like a typed pair that says so.
    assert untyped['interaction_types'] == [] and untyped['is_direct']
    assert (
        slot_name(signalink.map_interaction(untyped).predicate)
        == 'directly_physically_interacts_with'
    )
    # No type and no statement of directness: only that they are joined.
    undirect = {**untyped, 'is_direct': False}
    assert slot_name(signalink.map_interaction(undirect).predicate) == 'associated_with'
    direct_binding = {**phys, 'is_direct': True}
    assert (
        slot_name(signalink.map_interaction(direct_binding).predicate)
        == 'directly_physically_interacts_with'
    )


def test_rows_it_cannot_read_are_errors():
    bad_term = [(*EDGES[0][:4], 'is_direct:false|is_directed:true|MI:9999(x)', 'InnateDB')]
    with pytest.raises(ValueError, match='interaction type'):
        interactions(bad_term)
    undirected = [(*EDGES[0][:4], 'is_direct:false|MI:0915(physical association)', 'InnateDB')]
    with pytest.raises(ValueError, match='not directed'):
        interactions(undirected)
    unknown_node = [('Uniprot:P00000', *EDGES[0][1:])]
    with pytest.raises(ValueError, match='edge not understood'):
        interactions(unknown_node)
    bad_flag = [(*EDGES[0][:4], 'is_direct:maybe|is_directed:true', 'InnateDB')]
    with pytest.raises(ValueError, match='flag'):
        interactions(bad_flag)
    with pytest.raises(ValueError, match='no rows'):
        interactions([])
    with pytest.raises(ValueError):
        list(parser.iter_nodes(Source(NODE_HEADER).open()))
    # A source database outside the crosswalk is not dropped.
    row = interactions()[0]
    with pytest.raises(ValueError, match='Unknown SignaLink source'):
        signalink.map_interaction({**row, 'sources': ['Elsewhere']})


def test_no_pathway_membership_is_derived_for_the_interaction_partners():
    out = extract('interactions', interactions())
    # The interactions dataset says nothing of pathways or functions.
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert not [e for e in out.entities.values() if e.namespace.startswith('signalink_')]


def test_effect_terms_follow_the_psi_mi_vocabulary():
    # PSI-MI CV: MI:2243 down-regulates quantity by repression (of expression),
    # MI:2244 by destabilization (of stability); signor.py reads them the same way.
    assert parser.EFFECT_TERMS['MI:2243'] == ('decreased', 'expression')
    assert parser.EFFECT_TERMS['MI:2244'] == ('decreased', 'stability')
    assert parser.EFFECT_TERMS['MI:2239'] == ('increased', 'stability')
    assert parser.EFFECT_TERMS['MI:2242'] == ('decreased', 'abundance')


def test_cleavage_increases_the_cleavage_of_the_second_node():
    # Row of the layer 0 file: "protein cleavage" with the generic interaction type.
    edges = [
        (
            'Uniprot:O60674',
            'Uniprot:Q06124',
            '',
            'pubmed:8995399',
            'is_direct:false|is_directed:true|MI:0190(interaction type)|MI:0570(protein cleavage)',
            'SLKv2.1',
        )
    ]
    (row,) = interactions(edges)
    assert (row['direction'], row['aspect']) == ('increased', 'cleavage')
    out = extract('interactions', [row])
    attrs = relation_attributes(out.relations[0])
    assert out.relations[0].predicate == 'affects'
    assert attrs['object_direction_qualifier'] == ['increased']
    assert attrs['object_aspect_qualifier'] == ['cleavage']


def test_obsolete_other_modification_is_no_effect():
    # Row of the layer 0 file: InnateDB's only type for the pair is the obsolete
    # residue modification MI:0179 ("other modification"), which names no effect.
    nodes_text = NODE_HEADER + (
        '1\tUniprot:U3KQ93\t\ttaxid:9606\t\t\t\n'
        '2\tUniprot:P63165\t\ttaxid:9606\t\t\t\n'
    )
    edges = [
        (
            'Uniprot:U3KQ93',
            'Uniprot:P63165',
            'MI:0415(enzymatic study)',
            'pubmed:10187858|pubmed:23180781',
            'is_direct:false|is_directed:true|MI:0179(other modification)',
            'InnateDB',
        )
    ]
    (row,) = interactions(edges, nodes_text)
    assert not row['affects'] and row['direction'] is None
    record = signalink.map_interaction(row)
    assert slot_name(record.predicate) == 'interacts_with'
    assert relation_attributes(extract('interactions', [row]).relations[0])[
        'original_predicate'
    ] == ['MI:0179']


def test_the_vocabulary_dataset_is_called_terms():
    assert list(signalink.resource.datasets()) == ['terms', 'annotations', 'interactions']
