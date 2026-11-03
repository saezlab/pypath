"""Source fixture tests; no network calls or resource rebuilds."""

import pytest
from biolink_model.datamodel.model import Protein
from omnipath_core import Namespace
from pypath.inputs_v2 import intact, signor
from pypath.internals.tabular_builder import CV, EntityBuilder, IdentifiersBuilder


def mitab_row():
    return {
        '\ufeff#ID(s) interactor A': 'uniprotkb:P04637-2',
        '#ID(s) interactor A': 'uniprotkb:P04637-2',
        'ID(s) interactor B': 'uniprotkb:P38398-PRO_0000000123',
        'Alt. ID(s) interactor A': 'refseq:NM_000546.6|uniprotkb:P04637-3',
        'Type(s) interactor A': 'psi-mi:"MI:0326"(protein)',
        'Type(s) interactor B': 'psi-mi:"MI:0326"(protein)',
        'Taxid interactor A': 'taxid:9606(human)',
        'Taxid interactor B': 'taxid:9606(human)',
        'Causal statement': 'psi-mi:"MI:2236"(up-regulates activity)',
        'Interaction identifier(s)': 'intact:EBI-123',
        'Feature(s) interactor A': 'phosphorylated residue:15-15|mutation:175-175(p.R175H)',
        'Feature(s) interactor B': 'MOD:00696:27-27|phosphorylated residue:?-?',
    }


@pytest.mark.parametrize('builder', [signor.interactions_schema, intact.interactions_schema])
def test_pilot_sources_preserve_paired_isoform_chain_ptm_variant(builder):
    relation = builder(mitab_row())
    assert relation.subject.type == relation.object.type == 'protein'
    subject, object = relation.subject.molecular_form, relation.object.molecular_form
    assert subject['isoform_identifier'] == {'ns': 'uniprot', 'id': 'P04637-2'}
    assert subject['sequence_identifiers'] == [{'ns': 'uniprot', 'id': 'P04637-2'}]
    assert object['sequence_identifiers'][0]['id'] == 'P38398-PRO_0000000123'
    # General identifiers normalize the chain; occurrence retains its identity.
    assert relation.object.identifiers[0].value == 'P38398'
    ptm = subject['modifications'][0]
    assert ptm['term'] == 'phosphorylated residue'
    assert ptm['position'] == ptm['end_position'] == 15
    assert ptm['coordinate_reference']['identifier']['id'] == 'P04637-2'
    variant = subject['variants'][0]
    assert (variant['reference'], variant['position'], variant['alternate']) == ('R', 175, 'H')
    assert object['modifications'][1]['position'] is None
    assert object['modifications'][1]['description'] == 'phosphorylated residue:?-?'
    assert subject['protein_entity_key'] is None


def test_plain_protein_and_gene_crossreferences_do_not_assert_canonical_sequence():
    row = mitab_row()
    row['#ID(s) interactor A'] = row['\ufeff#ID(s) interactor A'] = 'uniprotkb:P04637'
    row['Feature(s) interactor A'] = '-'
    assert intact.interactions_schema(row).subject.molecular_form is None
    row['Type(s) interactor A'] = 'psi-mi:"MI:0250"(gene)'
    row['#ID(s) interactor A'] = row['\ufeff#ID(s) interactor A'] = 'entrezgene/locuslink:7157'
    row['Alt. ID(s) interactor A'] = 'uniprotkb:P04637-2'
    subject = intact.interactions_schema(row).subject
    assert subject.type == 'gene'
    assert subject.molecular_form is None


def test_unknown_sequence_and_fuzzy_feature_range_are_not_guessed():
    row = mitab_row()
    row['#ID(s) interactor A'] = row['\ufeff#ID(s) interactor A'] = 'uniprotkb:P04637'
    row['Feature(s) interactor A'] = 'phosphorylated residue:5..7-9..10'
    feature = intact.interactions_schema(row).subject.molecular_form['modifications'][0]
    assert feature['position'] is feature['end_position'] is None
    assert feature['coordinate_reference']['identifier'] is None
    assert feature['description'].endswith('5..7-9..10')


def test_builder_keeps_each_rows_form_and_rejects_resolved_keys():
    builder = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=lambda row: row['id'])),
        molecular_form=lambda row: {'isoform_identifier': {'ns': 'uniprot', 'id': row['isoform']}},
    )
    first = builder({'id': 'P04637', 'isoform': 'P04637-2'})
    second = builder({'id': 'P04637', 'isoform': 'P04637-3'})
    assert first.molecular_form != second.molecular_form
    first.molecular_form['isoform_identifier']['id'] = 'mutated'
    repeated = builder({'id': 'P04637', 'isoform': 'P04637-2'})
    assert repeated.molecular_form['isoform_identifier']['id'] == 'P04637-2'
    builder.molecular_form = {'protein_entity_key': 'forged'}
    with pytest.raises(ValueError, match='resolution'):
        builder({'id': 'P04637'})


def test_signor_nested_member_isoforms_are_not_merged_into_canonical_members():
    complex = signor.complexes_schema({
        'SIGNOR ID': 'SIGNOR-C1',
        'LIST OF ENTITIES': 'P04637-2,P04637-3,P38398',
    })
    members = [membership.member for membership in complex.membership]
    assert members[0].molecular_form['isoform_identifier']['id'] == 'P04637-2'
    assert members[1].molecular_form['isoform_identifier']['id'] == 'P04637-3'
    assert members[2].molecular_form is None


def test_bindingdb_retains_each_reported_chain_isoform_in_its_complex():
    from pypath.inputs_v2.bindingdb import _target
    from omnipath_core.biolink import entity_type
    row = {
        'Target Name': 'test complex',
        'Number of Protein Chains in Target (>1 implies a multichain complex)': '2',
        'UniProt (SwissProt) Primary ID of Target Chain 1': 'P04637-2',
        'UniProt (TrEMBL) Primary ID of Target Chain 2': 'P38398-PRO_0000000123',
    }
    complex = _target(row)
    assert entity_type(complex.type) == 'macromolecular_complex'
    members = [member.member for member in complex.membership]
    assert members[0].identifiers[0].value == 'P04637-2'
    assert members[0].molecular_form['isoform_identifier']['id'] == 'P04637-2'
    assert members[1].molecular_form['sequence_identifiers'][0]['id'] == 'P38398-PRO_0000000123'
    assert complex.molecular_form is None


def test_reactome_retains_exact_features_before_tabular_flattening():
    from rdflib import Graph, Literal, URIRef
    from rdflib.namespace import RDF
    from omnipath_core import Entity, Identifier
    from pypath.inputs_v2 import reactome
    from pypath.inputs_v2.parsers.reactome import BP, _extract_participant_data, _flatten_participants

    graph = Graph()
    molecule, reference = URIRef('urn:protein'), URIRef('urn:reference')
    graph.add((molecule, RDF.type, BP.Protein))
    graph.add((molecule, BP.entityReference, reference))
    for i, status in enumerate(['EQUAL', 'GREATER-THAN']):
        feature, vocabulary, site = (URIRef(f'urn:{name}{i}') for name in ('feature', 'mod', 'site'))
        graph.add((molecule, BP.feature, feature))
        graph.add((feature, BP.modificationType, vocabulary))
        graph.add((vocabulary, BP['term'], Literal('phosphorylated residue | source')))
        graph.add((feature, BP.featureLocation, site))
        graph.add((site, BP.sequencePosition, Literal(15 + i)))
        graph.add((site, BP.positionStatus, Literal(status)))
        if i == 0:
            xref = URIRef('urn:modxref')
            graph.add((vocabulary, BP.xref, xref))
            graph.add((xref, BP.db, Literal('PSI-MOD')))
            graph.add((xref, BP.id, Literal('MOD:00696')))
    index = {str(reference): {'entity': Entity('protein', [Identifier('uniprot', 'P04637-2')])}}
    participant = _extract_participant_data(graph, molecule, 'reactant', index, {}, {})
    row = {'entity_type': 'reaction', 'reactome_stable_id': 'R-HSA-1',
           **_flatten_participants([participant])}
    assert row['participant_molecular_form'].count('||') == 0
    member = reactome.reactions_schema(row).membership[0].member
    form = member.molecular_form
    assert form['isoform_identifier']['id'] == 'P04637-2'
    assert form['modifications'][0]['position'] == 15
    assert form['modifications'][0]['coordinate_reference']['identifier']['id'] == 'P04637-2'
    assert form['modifications'][1]['position'] is None
    assert 'GREATER-THAN' in form['modifications'][1]['description']
    assert form['modifications'][0]['term'] == 'MOD:00696'
    assert ' | source' in form['modifications'][0]['description']


def test_reactome_gene_type_stays_gene_and_unannotated_positions_stay_unknown():
    from rdflib import Graph, Literal, URIRef
    from pypath.inputs_v2.parsers.reactome import BP, PHYSICAL_ENTITY_TYPE_MAP, _participant_molecular_form
    assert PHYSICAL_ENTITY_TYPE_MAP['gene'] == 'gene'
    graph = Graph()
    molecule, feature, vocabulary, site = [URIRef('urn:' + name) for name in ('p', 'f', 'v', 's')]
    graph.add((molecule, BP.feature, feature))
    graph.add((feature, BP.modificationType, vocabulary))
    graph.add((vocabulary, BP['term'], Literal('phosphorylated residue')))
    graph.add((feature, BP.featureLocation, site))
    graph.add((site, BP.sequencePosition, Literal(22)))
    form = _participant_molecular_form(graph, molecule, {'entity_type': 'protein', 'uniprot': 'P04637'})
    assert form['modifications'][0]['position'] is None
    assert form['modifications'][0]['coordinate_reference']['identifier'] is None
    assert _participant_molecular_form(graph, molecule, {'entity_type': 'gene', 'uniprot': 'P04637-2'}) is None


def test_unknown_mitab_feature_label_does_not_acquire_guessed_ptm_semantics():
    row = mitab_row()
    row['Feature(s) interactor A'] = 'test phosphorylation experiment:15-15'
    form = intact.interactions_schema(row).subject.molecular_form
    assert form['modifications'] is None
    assert form['isoform_identifier']['id'] == 'P04637-2'
