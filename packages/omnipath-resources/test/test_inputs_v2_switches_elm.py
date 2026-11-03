"""Source contract of the switches.ELM input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor

from pypath.inputs_v2.parsers.switches_elm import iter_switches
from pypath.inputs_v2.switches_elm import resource

HEADER = (
    'Switch ID\tStatus\tInteraction ID\tIntramolecular\tID A\tBindingsite A ID\t'
    'Bindingsite A Start\tBindingsite A End\tID B\tBindingsite B ID\t'
    'Bindingsite B Start\tBindingsite B End\tAffected interactor\tSwitch type\t'
    'Switch subtype\tSwitch mechanism\tSwitch direction\tSwitch outcome direction\t'
    'Switch outcome\tModification\tModification sites\tModifying enzymes\tEffector\t'
    'Cell cycle phase\tLocalisation\tPathway\tPMID'
)
# Lines of switches.ELM-v1.txt (the real file ends its lines with a bare
# carriage return).
ROWS = {
    # A phosphotyrosine induces the binding of an SH2 domain.
    'induction': 'SWTI000001\tActive\tINTI000001\t\tUNIPROT:O43561\tELM:LIG_SH2_STAT5\t161\t164\tUNIPROT:P19174\tPFAM:PF00017\t550\t639\tID A\tBinary\tPhysicochemical compatibility\tPTM-dependent\tReversible\tPositive\tInduction\tMOD:00696\tY161\t\t\t\tGO:0009898\tKEGG:hsa04660;KEGG:hsa04064;KEGG:hsa04650;REACTOME:202403;REACTOME:2172127\tPMID: 20610546;PMID: 15065860',
    # The phosphorylation abolishes the binding.
    'abrogation': 'SWTI000201\tActive\tINTI000201\t\tUNIPROT:P05067\tELM:LIG_PTB_Apo_2\t756\t763\tUNIPROT:O00213\tPFAM:PF00640\t543\t664\tID A\tBinary\tPhysicochemical compatibility\tPTM-dependent\tReversible\tNegative\tAbrogation\tMOD:00696\tT743\t\t\t\t\t\tPMID: 18833287',
    # The modified protein is the one with the domain.
    'object': 'SWTI000009\tActive\tINTI000009\t\tUNIPROT:P11117\tELM:TRG_ENDOCYTIC_2\t413\t416\tUNIPROT:Q96CW1\tPFAM:PF00928\t159\t435\tID B\tUncategorised\tUncategorised\tPTM-dependent\tReversible\tPositive\tInduction\tMOD:00696\tT156\tUNIPROT:Q2M2I8\t\t\t\t\tPMID: 11877457',
    # Two motifs and two sites; the references have no space.
    'two sites': 'SWTI000317\tActive\tINTI000367\t\tUNIPROT:Q9P0K1\tELM:LIG_14-3-3_3\t831;854\t836;859\tUNIPROT:P31946\tPFAM:PF00244\t5\t238\tID A\tSpecificity\tMotif hiding\tPTM-dependent\tReversible\tPositive\tInduction\tMOD:00696\tS834;S857\t\t\t\t\t\tPMID: 16868027',
    # An isoform loses the motif; the domain occurs four times.
    'splicing': 'SWTI000499\tActive\tINTI000601\t\tUNIPROT:P97288-3\tELM:LIG_PDZ_Class_2\t368\t371\tUNIPROT:Q63ZW7\tPFAM:PF00595\t134;248;365;555\t221;328;453;641\tID A\tBinary\tPre-translational\tAlternative splicing\tIrreversible\tNegative\tAbrogation\t\t\t\t\t\t\t\tPMID:15466885',
    # An effector binds the participant B.
    'effector': 'SWTI000212\tActive\tINTI000212\t\tUNIPROT:P19814\tELM:TRG_ENDOCYTIC_2\t350\t353\tUNIPROT:P84092\tPFAM:PF00928\t1\t137\tID B\tBinary\tAllostery\tEffector binding\tReversible\tPositive\tInduction\t\t\t\tCHEBI:18348\t\t\t\tPMID: 20603002',
    # A proteolytic cleavage names no residue modification.
    'cleavage': 'SWTI000527\tActive\tINTI000630\t\tUNIPROT:P55211\tELM:LIG_BIR_III_2\t315\t319\tUNIPROT:P98170\tPFAM:PF00653\t268\t331\tID A\tBinary\tPhysicochemical compatibility\tPTM-dependent\tIrreversible\tPositive\tInduction\tProteolytic cleavage\tD315\t\t\t\t\t\tPMID: 12620238',
    # A small molecule, a complex formed by another interaction, and a switch
    # of the table that is inactive.
    'small molecule': 'SWTI000325\tActive\tINTI000383\t\tCHEBI:16618\t\t\t\tUNIPROT:P27986\tPFAM:PF00017\t624\t698\t\tSpecificity\tDomain hiding\t\tReversible\t\t\t\t\t\t\t\t\t\tPMID: 8521499',
    'complex': 'SWTI000298\tActive\tINTI000327\t\tUNIPROT:O43854\tELM:LIG_RGD\t96\t98\tINTI000326\t\t\t\t\tPre-assembly\tComposite binding site formation\t\tReversible\tPositive\tInduction\t\t\t\t\t\tGO:0031012\t\tPMID: 9420328',
    'inactive': 'SWTI000659\tInactive\tINTI000758\t\tUNIPROT:O43613\tELM:MOD_TYR_ITSM\t79\t86\tUNIPROT:Q06124\tPFAM:PF00017\t112\t197\tID A\tAvidity-sensing\t\tPTM-dependent\tReversible\tPositive\tInduction\tMOD:00696\tY83\t\t\t\tGO:0009898\t\tPMID: 19661287',
}
TABLE = '\r'.join([HEADER, *ROWS.values()]) + '\r'


def parse(table=TABLE):
    return {
        name: row
        for name, row in zip(
            ROWS, iter_switches(SimpleNamespace(result=io.StringIO(table)))
        )
    }


PARSED = parse()


def extract(name):
    row = PARSED[name]
    record = resource.interactions.mapper(row)
    assert record is not None, name
    out = SilverExtractor('switches_elm', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out.relations[0]


def terms(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_reader_decodes_a_line_of_a_file_with_carriage_returns():
    row = PARSED['induction']
    assert row['switch_id'] == 'SWTI000001'
    assert row['a'] == {
        'kind': 'uniprot',
        'id': 'O43561',
        'binding_sites': [{'id': 'ELM:LIG_SH2_STAT5', 'start': 161, 'end': 164}],
    }
    assert row['b']['binding_sites'] == [{'id': 'PFAM:PF00017', 'start': 550, 'end': 639}]
    assert row['affected'] == 'A'
    assert row['modification_sites'] == [{'residue': 'Y', 'position': 161}]
    assert row['pmids'] == ['20610546', '15065860']
    assert len(PARSED) == len(ROWS)


def test_reader_pairs_repeated_domains_and_motifs_and_keeps_other_participants():
    assert [s['start'] for s in PARSED['splicing']['b']['binding_sites']] == [134, 248, 365, 555]
    assert [s['id'] for s in PARSED['splicing']['b']['binding_sites']] == ['PFAM:PF00595'] * 4
    assert [s['start'] for s in PARSED['two sites']['a']['binding_sites']] == [831, 854]
    assert PARSED['small molecule']['a']['kind'] == 'chebi'
    assert PARSED['complex']['b'] == {'kind': 'interaction', 'id': 'INTI000326', 'binding_sites': []}


def test_reader_raises_on_a_defect_of_the_table():
    with pytest.raises(ValueError, match='lacks expected columns'):
        parse(HEADER.replace('PMID', 'Papers') + '\r')
    with pytest.raises(ValueError, match='no rows'):
        parse(HEADER + '\r')
    with pytest.raises(ValueError, match='fields'):
        parse(HEADER + '\r' + ROWS['induction'].replace('\t\t', '\t', 1) + '\r')
    unreadable = ROWS['induction'].replace('\tY161\t', '\tTyr161\t')
    with pytest.raises(ValueError, match='modification site'):
        parse(HEADER + '\r' + unreadable + '\r')
    unpaired = ROWS['induction'].replace('\t161\t164\t', '\t161;170\t164\t')
    with pytest.raises(ValueError, match='start and'):
        parse(HEADER + '\r' + unpaired + '\r')


def test_the_interaction_is_between_two_proteins_with_their_binding_sites():
    record, relation = extract('abrogation')
    assert relation.predicate == 'interacts_with'
    regions = record.subject.molecular_form['regions']
    assert regions == [
        {
            'type': 'binding site',
            'identifier': {'ns': 'elm', 'id': 'LIG_PTB_Apo_2'},
            'position': 756,
            'end_position': 763,
            'coordinate_reference': {
                'identifier': None,
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': None,
        }
    ]
    assert record.object.molecular_form['regions'][0]['identifier'] == {'ns': 'pfam', 'id': 'PF00640'}
    assert terms(relation, 'publications') == ['PMID:18833287']
    # The table says nothing of the effect of the motif interaction itself.
    assert terms(relation, 'object_direction_qualifier') == []


def test_a_modification_that_induces_the_interaction_is_on_the_bound_form():
    record, _ = extract('induction')
    (modification,) = record.subject.molecular_form['modifications']
    assert modification['term'] == 'MOD:00696'
    assert (modification['residue'], modification['position']) == ('Y', 161)
    assert record.object.molecular_form['modifications'] is None


def test_the_modification_is_on_the_participant_that_the_table_names():
    record, _ = extract('object')
    assert record.subject.molecular_form['modifications'] is None
    (modification,) = record.object.molecular_form['modifications']
    assert (modification['residue'], modification['position']) == ('T', 156)


def test_all_the_modified_residues_of_a_switch_are_on_the_form():
    record, _ = extract('two sites')
    assert [m['position'] for m in record.subject.molecular_form['modifications']] == [834, 857]


def test_a_modification_that_abrogates_the_interaction_is_not_on_the_bound_form():
    record, _ = extract('abrogation')
    assert record.subject.molecular_form.get('modifications') is None


def test_an_effector_a_splice_form_and_a_cleavage_give_no_modification():
    for name in ('effector', 'splicing', 'cleavage'):
        record, _ = extract(name)
        assert record.subject.molecular_form.get('modifications') is None
        assert record.object.molecular_form.get('modifications') is None


def test_an_isoform_names_its_sequence_as_the_reference_of_the_positions():
    record, _ = extract('splicing')
    region = record.subject.molecular_form['regions'][0]
    assert region['coordinate_reference']['identifier'] == {'ns': 'uniprot', 'id': 'P97288-3'}
    assert record.object.molecular_form['regions'][1]['position'] == 248


def test_only_a_pair_of_proteins_of_an_active_switch_gives_a_relation():
    for name in ('small molecule', 'complex', 'inactive'):
        assert resource.interactions.mapper(PARSED[name]) is None


def test_an_unknown_status_outcome_or_mechanism_raises():
    with pytest.raises(ValueError, match='status'):
        resource.interactions.mapper({**PARSED['induction'], 'status': 'Pending'})
    with pytest.raises(ValueError, match='outcome'):
        resource.interactions.mapper({**PARSED['induction'], 'outcome': 'Reversal'})
    with pytest.raises(ValueError, match='mechanism'):
        resource.interactions.mapper({**PARSED['induction'], 'switch_mechanism': 'Cleavage'})
