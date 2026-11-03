"""Source contract of the scConnect ligand-receptor input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import PARTICIPANT_ROLE

from pypath.inputs_v2 import scconnect as mod
from pypath.inputs_v2.parsers import scconnect as parser

FIELDS = (
    'interaction_id', 'target', 'target_gene_symbol', 'target_uniprot',
    'target_ensembl_gene_id', 'target_species', 'ligand', 'ligand_gene_symbol',
    'ligand_species', 'type', 'action', 'endogenous', 'pubmed_id',
)


def row(*values):
    return dict(zip(FIELDS, values, strict=True))


# Rows of GTP_tables_clean/2020-5/interactions.csv, with the columns read
# and the ligand's type, action and endogenous flag.
PEPTIDE = row(
    '2019', 'apelin receptor', 'APLNR', 'P35414', 'ENSG00000134817', 'Human',
    'apelin-36', 'APLN', 'Human', 'Agonist', 'Full agonist', 't',
    '10777510|11336787|12603839|12939143',
)
NEUROTRANSMITTER = row(
    '31', '5-HT1A receptor', 'HTR1A', 'P08908', 'ENSG00000178394', 'Human',
    '5-hydroxytryptamine', '', '', 'Agonist', 'Full agonist', 't',
    '1386736|9760039|10431754|9550290|9205951|15628665',
)
RAT_NEUROTRANSMITTER = row(
    '6060', 'D2 receptor', 'Drd2', 'P61169', 'ENSRNOG00000008428', 'Rat',
    'dopamine', '', '', 'Agonist', 'Full agonist', 't', '1975644|10617689',
)
# The ligand's genes are not in the order of its species.
UNORDERED_SPECIES = row(
    '2018', 'apelin receptor', 'APLNR', 'P35414', 'ENSG00000134817', 'Human',
    'apelin-13', 'APLN|Apln|Apln', 'Human|Mouse|Rat', 'Agonist',
    'Full agonist', 't', '10777510|12603839|12939143',
)
SAME_SYMBOL = row(
    '1388', 'ACKR3', 'Ackr3', 'O89039', 'ENSRNOG00000019622', 'Rat', 'α-CGRP',
    'Calca|Calca', 'Mouse|Rat', 'Agonist', 'Agonist', 'f', '8554605',
)
OTHER_SPECIES = row(
    '8123', 'GAL1 receptor', 'GALR1', 'P47211', 'ENSG00000166573', 'Human',
    'galanin', 'GAL', 'Pig', 'Agonist', 'Full agonist', 'f',
    '9808667|9880084|7524088',
)
SEVERAL_TARGET_GENES = row(
    '10550', 'Interleukin-2 receptor', 'IL2RG|IL2RA|IL2RB',
    'P31785|P01589|P14784',
    'ENSG00000100385|ENSG00000134460|ENSG00000147168', 'Human', 'IL-2', 'IL2',
    'Human', 'Agonist', 'Agonist', 't', '',
)
HETEROMERIC_RECEPTOR = row(
    '912', '5-HT3AB', 'HTR3B|HTR3A', 'P46098|O95264',
    'ENSG00000166736|ENSG00000149305', 'Human', '5-hydroxytryptamine', '', '',
    'Agonist', 'Agonist', 't', '9950429|10521471|16081679|10825381',
)
# Two rat paralogues, not the subunits of a complex.
PARALOGUES = row(
    '2632', 'AT1 receptor', 'Agtr1a|Agtr1b', 'P29089|P25095',
    'ENSRNOG00000010640|ENSRNOG00000018346', 'Rat', 'angiotensin A', 'AGT',
    'Human', 'Agonist', 'Partial agonist', 't', '17138938',
)
SEVERAL_LIGAND_GENES = row(
    '7682', 'FSH receptor', 'FSHR', 'P23945', 'ENSG00000170820', 'Human',
    'FSH', 'CGA|FSHB', 'Human', 'Agonist', 'Agonist', 't', '',
)
SYNTHETIC = row(
    '0', '12S-LOX', 'ALOX12', 'P18054', 'ENSG00000108839', 'Human', 'ML355',
    '', '', 'Inhibitor', 'Inhibition', 'f', '24393039',
)
NO_TARGET_GENE = row(
    '3649', 'Calmodulin-dependent kinase II', '', '', '', 'Human',
    'calmodulin', 'CALM2|CALM3|CALM1', 'Human', 'Activator', 'Activation',
    't', '',
)
NO_TARGET = row(
    '19872', '', '', '', '', 'Human', 'BAFF', 'TNFSF13B', 'Human', 'None',
    'Binding', 'f', '',
)


# Rows of the file whose action is not an agonist's.
ANTAGONIST = row(
    '4304', 'CCR1', 'CCR1', 'P32246', 'ENSG00000163823', 'Human', 'CCL18',
    'CCL18', 'Human', 'Antagonist', 'Antagonist', 'f', '23951310',
)
INVERSE_AGONIST = row(
    '12250', 'MC1 receptor', 'Mc1r', 'Q01727', 'ENSMUSG00000074037', 'Mouse',
    'agouti', 'a', 'Mouse', 'Antagonist', 'Inverse agonist', 't',
    '7935841|9454589',
)
INHIBITION = row(
    '16742', 'Proton-coupled Amino acid Transporter 1', 'SLC36A1', 'Q7Z2H8',
    'ENSG00000123643', 'Human', '5-hydroxytryptamine', '', '', 'Inhibitor',
    'Inhibition', 'f', '16126914',
)
PARTIAL_AGONIST = row(
    '1756', 'α2B-adrenoceptor', 'ADRA2B', 'P18089', 'ENSG00000222040', 'Human',
    '(-)-noradrenaline', '', '', 'Agonist', 'Partial agonist', 't',
    '9824686|9605427',
)
POSITIVE = row(
    '15136', 'P2X7', 'P2RX7', 'Q99572', 'ENSG00000089041', 'Human', 'LL-37',
    'CAMP', 'Human', 'Allosteric modulator', 'Positive', 'f', '18765670',
)
BINDING = row(
    '6532', 'epidermal growth factor receptor', 'EGFR', 'P00533',
    'ENSG00000146648', 'Human', 'amphiregulin', 'AREG', 'Human', 'None',
    'Binding', 't', '',
)
NO_ACTION = row(
    '16', '4-1BB', 'TNFRSF9', 'Q07011', 'ENSG00000049249', 'Human',
    '4-1BB ligand', 'TNFSF9', 'Human', 'None', 'None', 't', '',
)
UNKNOWN_ACTION = row(
    '4465', 'CCRL2', 'CCRL2', 'O00421', 'ENSG00000121797', 'Human', 'CCL19',
    'CCL19', 'Human', 'Agonist', 'Unknown', 't', '20002784',
)


def extract(source_row):
    record = mod.resource.datasets()['interactions'].mapper(source_row)
    assert record is not None
    out = SilverExtractor('scconnect', 'interactions')
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


def test_peptide_pair_is_ligand_to_receptor_with_source_and_publications():
    relation, ligand, receptor = extract(PEPTIDE)
    # A full agonist: an effect of the ligand on the target.
    assert relation.predicate == 'affects'
    assert identifiers(ligand) == {(Namespace.GENESYMBOL, 'APLN')}
    assert identifiers(receptor) == {
        (Namespace.GENESYMBOL, 'APLNR'),
        (Namespace.UNIPROT, 'P35414'),
        (Namespace.ENSG, 'ENSG00000134817'),
        (Namespace.NAME, 'apelin receptor'),
    }
    assert {ligand.entity_type, receptor.entity_type} == {'protein'}
    assert {ligand.taxon, receptor.taxon} == {'9606'}
    assert annotations(relation) == [
        ('causal_mechanism_qualifier', 'agonism'),
        ('object_direction_qualifier', 'increased'),
        ('publications', 'PMID:10777510'),
        ('publications', 'PMID:11336787'),
        ('publications', 'PMID:12603839'),
        ('publications', 'PMID:12939143'),
        ('supporting_data_source', 'Guide to Pharmacology'),
    ]


def test_pair_asserts_no_ligand_type_or_endogenous_flag():
    relation, ligand, _ = extract(PEPTIDE)
    # The action is carried as qualifiers; the ligand type ('Agonist') and the
    # endogenous flag ('t') are not.
    assert {
        a['term'] for a in relation.annotations if a['scope'] == 'relation'
    } == {
        'publications',
        'supporting_data_source',
        'causal_mechanism_qualifier',
        'object_direction_qualifier',
    }
    assert all(
        a['value'] not in {'Agonist', 'Full agonist', 't'}
        for a in relation.annotations
    )
    # The peptide's name is not a name of its precursor's gene product.
    assert 'apelin-36' not in {value for _, value in identifiers(ligand)}


@pytest.mark.parametrize(
    'source_row,predicate,mechanism,direction',
    [
        (ANTAGONIST, 'affects', 'antagonism', 'decreased'),
        (INVERSE_AGONIST, 'affects', 'inverse_agonism', 'decreased'),
        (INHIBITION, 'affects', 'inhibition', 'decreased'),
        (PARTIAL_AGONIST, 'affects', 'partial_agonism', 'increased'),
        (PEPTIDE, 'affects', 'agonism', 'increased'),
        (SAME_SYMBOL, 'affects', 'agonism', 'increased'),
        (POSITIVE, 'affects', 'positive_modulation', None),
        # As in guidetopharma.py, binding is an interaction, not an effect.
        (BINDING, 'interacts_with', 'binding', None),
        (NO_ACTION, 'interacts_with', None, None),
        (UNKNOWN_ACTION, 'interacts_with', None, None),
    ],
)
def test_action_is_the_mechanism_and_direction_of_guidetopharma(
    source_row, predicate, mechanism, direction
):
    relation, _, _ = extract(source_row)
    assert relation.predicate == predicate
    found = {
        a['term']: a['value']
        for a in relation.annotations
        if a['term'] in {'causal_mechanism_qualifier', 'object_direction_qualifier'}
    }
    assert found.get('causal_mechanism_qualifier') == mechanism
    assert found.get('object_direction_qualifier') == direction


def test_antagonists_differ_from_agonists_of_the_same_pair():
    # Before the action was carried, these were identical to an agonist's row.
    agonist = dict(ANTAGONIST, action='Agonist')
    assert extract(agonist)[0].relation_key != extract(ANTAGONIST)[0].relation_key


def test_unknown_action_is_rejected_for_every_row():
    for source_row in (PEPTIDE, SYNTHETIC):
        with pytest.raises(ValueError, match="Unknown scConnect action 'Blocker'"):
            mod.map_scconnect_interaction(dict(source_row, action='Blocker'))
    with pytest.raises(ValueError, match="Unknown scConnect action ''"):
        mod.map_scconnect_interaction(dict(PEPTIDE, action=''))


def test_every_action_of_the_file_is_in_the_table():
    assert set(mod.ACTIONS) == {
        'Activation', 'Agonist', 'Full agonist', 'Partial agonist',
        'Biased agonist', 'Antagonist', 'Inhibition', 'Feedback inhibition',
        'Irreversible inhibition', 'Inverse agonist',
        'Voltage-dependent inhibition', 'Positive', 'Negative',
        'Potentiation', 'Pore blocker', 'Binding', 'Mixed', 'Biphasic',
        'Competitive', 'Non-competitive', 'Neutral', 'Slows inactivation',
        'Irreversible agonist', 'Unknown', 'None',
    }  # fmt: skip


@pytest.mark.parametrize('source_row', [PEPTIDE, NEUROTRANSMITTER, INHIBITION])
def test_ligand_column_gives_the_ligand_role_and_the_target_has_none(source_row):
    relation, ligand, receptor = extract(source_row)
    assert [
        v for t, v in annotations(relation, 'subject') if t == PARTICIPANT_ROLE
    ] == ['ligand']
    # The table's column is `target`, not receptor: targets include enzymes
    # and transporters (INHIBITION is serotonin on the amino acid transporter
    # SLC36A1), so the target is given no role.
    assert all(t != PARTICIPANT_ROLE for t, _ in annotations(relation, 'object'))
    assert all(t != PARTICIPANT_ROLE for t, _ in annotations(relation))


def test_neurotransmitter_is_a_chemical_not_its_enzymes():
    relation, ligand, receptor = extract(NEUROTRANSMITTER)
    assert ligand.entity_type == 'chemical_entity'
    assert identifiers(ligand) == {(Namespace.NAME, '5-hydroxytryptamine')}
    assert ligand.taxon is None
    # The source's ligand column makes it a ligand; it has no species.
    assert annotations(relation, 'subject') == [(PARTICIPANT_ROLE, 'ligand')]
    assert (Namespace.GENESYMBOL, 'HTR1A') in identifiers(receptor)
    assert len(annotations(relation)) == 9
    # scConnect scores serotonin by TPH2, DDC and SLC18A2; none is a ligand.
    out = SilverExtractor('scconnect', 'interactions')
    out.process_record(
        mod.map_scconnect_interaction(NEUROTRANSMITTER), NEUROTRANSMITTER, 'f:1', 1
    )
    assert {e.entity_type for e in out.entities.values()} == {
        'chemical_entity',
        'protein',
    }
    assert len(out.entities) == 2


def test_species_is_the_targets_and_is_not_translated():
    _, ligand, receptor = extract(RAT_NEUROTRANSMITTER)
    assert receptor.taxon == '10116'
    assert identifiers(receptor) == {
        (Namespace.GENESYMBOL, 'Drd2'),
        (Namespace.UNIPROT, 'P61169'),
        (Namespace.ENSG, 'ENSRNOG00000008428'),
        (Namespace.NAME, 'D2 receptor'),
    }
    assert ligand.taxon is None


def test_same_symbol_in_every_species_identifies_the_ligand_gene():
    _, ligand, receptor = extract(SAME_SYMBOL)
    assert identifiers(ligand) == {(Namespace.GENESYMBOL, 'Calca')}
    assert {ligand.taxon, receptor.taxon} == {'10116'}


@pytest.mark.parametrize(
    'source_row,reason',
    [
        (UNORDERED_SPECIES, 'ligand gene of the target species is not identified'),
        (OTHER_SPECIES, 'ligand is from another species than the target'),
        (SEVERAL_TARGET_GENES, 'target is named by several genes'),
        (HETEROMERIC_RECEPTOR, 'target is named by several genes'),
        (PARALOGUES, 'target is named by several genes'),
        (SEVERAL_LIGAND_GENES, 'ligand is named by several genes'),
        (SYNTHETIC, 'ligand is not a peptide with a gene or a scConnect molecule'),
        (NO_TARGET_GENE, 'target has no gene'),
        (NO_TARGET, 'target has no gene'),
        (dict(PEPTIDE, target_species='Monkey', ligand_species='Monkey'),
         'target species has no taxon'),
    ],
)
def test_rows_that_yield_no_relation(source_row, reason):
    assert mod.exclusion(source_row) == reason
    assert mod.map_scconnect_interaction(source_row) is None


def test_genes_joined_in_a_cell_are_never_a_complex():
    # Neither file read by scConnect says whether such genes are subunits.
    for source_row in (SEVERAL_TARGET_GENES, SEVERAL_LIGAND_GENES, PARALOGUES):
        assert mod.map_scconnect_interaction(source_row) is None
    for source_row in (PEPTIDE, NEUROTRANSMITTER, SAME_SYMBOL):
        out = SilverExtractor('scconnect', 'interactions')
        record = mod.map_scconnect_interaction(source_row)
        out.process_record(record, source_row, 'fixture:1', 1)
        assert [r.predicate for r in out.relations] == ['affects']
        assert 'macromolecular_complex' not in {
            e.entity_type for e in out.entities.values()
        }


def test_pair_without_or_with_malformed_pmids():
    relation, _, _ = extract(dict(PEPTIDE, pubmed_id=''))
    assert annotations(relation) == [
        ('causal_mechanism_qualifier', 'agonism'),
        ('object_direction_qualifier', 'increased'),
        ('supporting_data_source', 'Guide to Pharmacology'),
    ]
    # The cell holds PubMed IDs only, so a short number is one (the shared rule
    # is 1 to 8 digits); a leading zero, nine digits and text are not.
    relation, _, _ = extract(dict(PEPTIDE, pubmed_id='10777510|10777510|953|n/a'))
    assert [v for t, v in annotations(relation) if t == 'publications'] == [
        'PMID:10777510',
        'PMID:953',
    ]
    relation, _, _ = extract(
        dict(PEPTIDE, pubmed_id='0123456|123456789|10777510')
    )
    assert [v for t, v in annotations(relation) if t == 'publications'] == [
        'PMID:10777510'
    ]


def table(lines):
    return SimpleNamespace(result=io.StringIO('\n'.join(lines) + '\n'))


def test_reader_reads_quoted_cells_and_keeps_every_column():
    header = ','.join(FIELDS)
    line = (
        '63,5-HT1A receptor,HTR1A,P08908,ENSG00000178394,Human,'
        '"compound 34 [PMID: 20866075], racemate",,,Agonist,Agonist,f,20866075'
    )
    (record,) = parser.iter_interactions(table([header, line]))
    assert record['ligand'] == 'compound 34 [PMID: 20866075], racemate'
    assert record['target_uniprot'] == 'P08908'
    assert record['action'] == 'Agonist'


def test_reader_rejects_a_table_with_changed_columns():
    with pytest.raises(ValueError, match='ligand_gene_symbol'):
        list(parser.iter_interactions(table(['ligand,target', 'IL-2,IL2RA'])))


def test_reader_rejects_an_empty_download():
    with pytest.raises(ValueError, match='no rows'):
        list(parser.iter_interactions(table([','.join(FIELDS)])))
