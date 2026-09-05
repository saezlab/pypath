"""Regression cases from the resource-by-resource web review."""
import pytest
from test.test_inputs_v2_chemical_migration import project
from test.test_inputs_v2_pathway_migration import extract


def test_foodb_list_columns_preserve_empty_zero_and_context():
    out = project('foodb', 'foods_schema', {
        'public_id': 'FOOD1',
        'member_compound_public_id': ['FDB1', 'FDB2', 'FDB3'],
        'member_compound_name': ['Alpha', 'Beta', 'Gamma'],
        'member_content': [None, '0', '3.5'],
        'member_unit': ['mg/g', 'ug/g', 'mg/100g'],
        'member_citation': ['', 'PMID:123', ''],
        'member_food_part': ['', 'leaf', ''],
        'member_method': ['', 'HPLC', ''],
    })
    assert len(out.relations) == 3
    quantities = [[a['quantity'] for a in r.annotations if a.get('quantity')] for r in out.relations]
    assert quantities[0] == []
    assert quantities[1][0]['has_numeric_value'] == 0
    assert quantities[1][0]['has_unit'] == 'ug/g'
    assert quantities[2][0]['has_numeric_value'] == 3.5
    assert {a['value'] for a in out.relations[1].annotations if a.get('value')} >= {'PMID:123', 'Food part: leaf', 'HPLC'}
    assert {'Alpha', 'Beta', 'Gamma'} <= {i['id'] for e in out.entities.values() for i in e.identifiers}


@pytest.mark.parametrize('module,id_field', [('metatlas', 'human_gem_reaction_id'), ('recon3d', 'bigg_reaction_id')])
def test_transport_keeps_compartment_on_each_participant(module, id_field):
    out = extract(module, 'reactions_schema', {id_field: 'R1', 'reactants': 'h2o:c:2', 'products': 'h2o:e:1'})
    assert {(r.predicate, a['value']) for r in out.relations for a in r.annotations if a['term'] == 'biopax:cellularLocation'} == {('has_input', 'c'), ('has_output', 'e')}


def test_reactome_physical_state_survives_parser_and_mapping():
    from rdflib import Graph, URIRef, Literal, RDF
    from pypath.inputs_v2.parsers import reactome as parser
    g = Graph()
    p, ref, loc, feature, vocab, site = [URIRef('https://example.org/' + x) for x in ('protein', 'ref', 'loc', 'feature', 'vocab', 'site')]
    for triple in [(p, RDF.type, parser.BP.Protein), (p, parser.BP.entityReference, ref),
                   (p, parser.BP.cellularLocation, loc), (loc, parser.BP['term'], Literal('cytosol')),
                   (p, parser.BP.feature, feature), (feature, parser.BP.modificationType, vocab),
                   (vocab, parser.BP['term'], Literal('phosphoserine')), (feature, parser.BP.featureLocation, site),
                   (site, parser.BP.sequencePosition, Literal(259))]:
        g.add(triple)
    participant = parser._extract_participant_data(g, p, 'reactant', {}, {}, {})
    assert participant['compartment'] == 'cytosol'
    assert participant['modification'] == 'phosphoserine; position 259'
    out = extract('reactome', 'reactions_schema', {
        'entity_type': 'reaction', 'reactome_stable_id': 'R-HSA-123',
        'participant_entity_type': 'protein', 'participant_uniprot': 'P12345',
        'participant_role': 'reactant', 'participant_source_physical_entity': participant['source_physical_entity'],
        'participant_compartment': participant['compartment'], 'participant_modification': participant['modification'],
    })
    annotations = {a['term']: a['value'] for r in out.relations for a in r.annotations}
    assert annotations['biopax:cellularLocation'] == 'cytosol'
    assert annotations['biopax:feature'] == 'phosphoserine; position 259'
    assert annotations['source_record_urls'] == str(p)


def test_cellchat_role_does_not_rename_reference_protein():
    from pypath.inputs_v2.cellchat import map_cellchat_cofactor_group, map_cellchat_cofactor_interaction
    entity = map_cellchat_cofactor_group({'name': 'ACTIVIN antagonist', 'genes': ['FST'], 'taxon_id': '9606'})
    assert 'ACTIVIN antagonist' not in {i.value for i in entity.identifiers}
    relation = map_cellchat_cofactor_interaction({'cofactor_gene': 'FST', 'cofactor_group': 'ACTIVIN antagonist', 'cofactor_role': 'antagonist', 'effect': 'inhibition', 'target_name': 'INHBA', 'target_genes': ['INHBA'], 'taxon_id': '9606'})
    from omnipath_core.biolink import annotation_term, annotation_value
    assert {annotation_term(a.term): annotation_value(a.term, a.value) for a in relation.annotations}['original_subject'] == 'ACTIVIN antagonist'
    assert {annotation_term(a.term): annotation_value(a.term, a.value) for a in relation.annotations}['object_direction_qualifier'] == 'decreased'


def test_bindingdb_construct_is_observation_context():
    from pypath.inputs_v2.bindingdb import interactions_schema
    record = interactions_schema({'BindingDB MonomerID': '1', 'Target Name': 'Dimer of polyprotein [501-599]', 'UniProt (SwissProt) Primary ID of Target Chain 1': 'P03367', 'Ki (nM)': '0.99'})
    assert 'Dimer of polyprotein [501-599]' not in {i.value for i in record.object.identifiers}
    assert {a.term: a.value for a in record.annotations}['original_object'] == 'Dimer of polyprotein [501-599]'


def test_unspecified_diagram_edges_do_not_claim_interaction():
    from pypath.inputs_v2.wikipathways import wikipathways_predicate
    from omnipath_core.biolink import predicate
    assert predicate(wikipathways_predicate({'interaction_types': ['Conversion']})) == 'related_to'
    assert predicate(wikipathways_predicate({'interaction_types': ['Binding']})) == 'physically_interacts_with'
