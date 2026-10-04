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


def test_builder_cache_preserves_distinct_forms_even_with_explicit_cache_by():
    calls = []
    builder = EntityBuilder(
        entity_type=Protein, cache_by=('id',),
        identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT,
                                         value=lambda row: calls.append(row['id']) or row['id'])),
        molecular_form=lambda row: {'isoform_identifier': {'ns': 'uniprot', 'id': row['isoform']}},
    )
    first = builder({'id': 'P04637', 'isoform': 'P04637-2'})
    second = builder({'id': 'P04637', 'isoform': 'P04637-3'})
    assert first.molecular_form != second.molecular_form
    first.molecular_form['isoform_identifier']['id'] = 'mutated'
    repeated = builder({'id': 'P04637', 'isoform': 'P04637-2'})
    assert repeated.molecular_form['isoform_identifier']['id'] == 'P04637-2'
    assert calls == ['P04637', 'P04637']
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
