"""Source contract of the CellCall ligand-receptor input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import PARTICIPANT_ROLE

from pypath.inputs_v2.cellcall import resource
from pypath.inputs_v2.parsers.cellcall import COLUMNS, iter_pairs


HEADER = '\t'.join(COLUMNS)
# Lines of new_ligand_receptor_TFs.txt (human core set).
CCL26_CCR1_LINES = [
    '10344\t1230\t2309\thsa04062_5\thsa04062\tCCL26\tCCR1\tFOXO3',
    '10344\t1230\t3551\thsa04062_5\thsa04062\tCCL26\tCCR1\tIKBKB',
]
TGFB1_COMPLEX_LINES = [
    '7040\t7046,7048\t4087\t'
    'hsa04926_1,hsa05200_1,hsa05212_1,hsa05225_1,hsa05226_1\t'
    'hsa04926,hsa05200,hsa05212,hsa05225,hsa05226\tTGFB1\tTGFBR1,TGFBR2\tSMAD2',
    '7040\t7046,7048\t4089\thsa05200_2,hsa05212_2,hsa05225_2,hsa05226_2\t'
    'hsa05200,hsa05212,hsa05225,hsa05226\tTGFB1\tTGFBR1,TGFBR2\tSMAD4',
]
# Lines 373, 382 and 2802 of new_ligand_receptor_TFs_homology.txt (mouse core
# set) are the same row; the pathway identifiers keep the human prefix.
GH_PRLR_LINE = '14599\t19116\t17869\thsa04630_4\thsa04630\tGh\tPrlr\tMyc'
IFNG_IFNGR1_LINE = (
    '15978\t15979\t12912\thsa04380_6\thsa04380\tIfng\tIfngr1\tCreb1'
)

# The same pairs as the reader yields them.
PROTEIN_PAIR = {
    'Ligand_ID': '10344',
    'Ligand_Symbol': 'CCL26',
    'Receptor_ID': '1230',
    'Receptor_Symbol': 'CCR1',
    'taxon_id': '9606',
    'transcription_factors': [
        {
            'TF_ID': '2309',
            'TF_Symbol': 'FOXO3',
            'Pathway': 'hsa04062_5',
            'pathway_ID': 'hsa04062',
        },
        {
            'TF_ID': '3551',
            'TF_Symbol': 'IKBKB',
            'Pathway': 'hsa04062_5',
            'pathway_ID': 'hsa04062',
        },
    ],
}
COMPLEX_PAIR = {
    'Ligand_ID': '7040',
    'Ligand_Symbol': 'TGFB1',
    'Receptor_ID': '7046,7048',
    'Receptor_Symbol': 'TGFBR1,TGFBR2',
    'taxon_id': '9606',
    'transcription_factors': [
        {
            'TF_ID': '4087',
            'TF_Symbol': 'SMAD2',
            'Pathway': 'hsa04926_1,hsa05200_1,hsa05212_1,hsa05225_1,hsa05226_1',
            'pathway_ID': 'hsa04926,hsa05200,hsa05212,hsa05225,hsa05226',
        },
        {
            'TF_ID': '4089',
            'TF_Symbol': 'SMAD4',
            'Pathway': 'hsa05200_2,hsa05212_2,hsa05225_2,hsa05226_2',
            'pathway_ID': 'hsa05200,hsa05212,hsa05225,hsa05226',
        },
    ],
}
MOUSE_PAIR = {
    'Ligand_ID': '14599',
    'Ligand_Symbol': 'Gh',
    'Receptor_ID': '19116',
    'Receptor_Symbol': 'Prlr',
    'taxon_id': '10090',
    'transcription_factors': [
        {
            'TF_ID': '17869',
            'TF_Symbol': 'Myc',
            'Pathway': 'hsa04630_4',
            'pathway_ID': 'hsa04630',
        },
    ],
}
# The tight junction proteins pair with themselves in the source.
SELF_PAIR = {
    'Ligand_ID': '100506658',
    'Ligand_Symbol': 'OCLN',
    'Receptor_ID': '100506658',
    'Receptor_Symbol': 'OCLN',
    'taxon_id': '9606',
    'transcription_factors': [
        {
            'TF_ID': '2626',
            'TF_Symbol': 'GATA4',
            'Pathway': 'hsa04530_3,hsa04530_5',
            'pathway_ID': 'hsa04530',
        },
    ],
}


def table(lines):
    # The source files have Windows line ends.
    return SimpleNamespace(result=io.StringIO('\r\n'.join(lines) + '\r\n', newline=''))


def extract(row, dataset='human_core_interactions'):
    record = resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor('cellcall', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def pair(out):
    (relation,) = [r for r in out.relations if r.predicate != 'has_member']
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def test_datasets_are_one_per_set_and_species():
    assert set(resource.datasets()) == {
        'human_core_interactions',
        'human_extended_interactions',
        'mouse_core_interactions',
        'mouse_extended_interactions',
    }


def test_pair_is_ligand_to_receptor_with_symbols_and_gene_ids():
    record, out = extract(PROTEIN_PAIR)
    relation, ligand, receptor = pair(out)
    assert [r.predicate for r in out.relations] == ['interacts_with']
    assert identifiers(ligand) == {
        (Namespace.GENESYMBOL, 'CCL26'),
        (Namespace.ENTREZ, '10344'),
    }
    assert identifiers(receptor) == {
        (Namespace.GENESYMBOL, 'CCR1'),
        (Namespace.ENTREZ, '1230'),
    }
    assert record.subject.identifiers[0].type == Namespace.GENESYMBOL
    assert {ligand.entity_type, receptor.entity_type} == {'protein'}
    assert {ligand.taxon, receptor.taxon} == {'9606'}


def test_multi_subunit_receptor_is_a_complex_of_the_listed_subunits():
    _, out = extract(COMPLEX_PAIR)
    relation, ligand, receptor = pair(out)
    assert relation.predicate == 'interacts_with'
    assert ligand.entity_type == 'protein'
    assert receptor.entity_type == 'macromolecular_complex'
    # The source has no accession for the complex, only the subunit list.
    assert identifiers(receptor) == {(Namespace.NAME, 'TGFBR1,TGFBR2')}
    assert receptor.taxon == '9606'
    members = [r for r in out.relations if r.predicate == 'has_member']
    assert {m.subject_entity_key for m in members} == {receptor.entity_key}
    # Symbols and gene IDs are paired by their position in the two cells.
    assert [
        identifiers(out.entities[m.object_entity_key]) for m in members
    ] == [
        {(Namespace.GENESYMBOL, 'TGFBR1'), (Namespace.ENTREZ, '7046')},
        {(Namespace.GENESYMBOL, 'TGFBR2'), (Namespace.ENTREZ, '7048')},
    ]
    assert {
        (out.entities[m.object_entity_key].entity_type,
         out.entities[m.object_entity_key].taxon)
        for m in members
    } == {('protein', '9606')}


def roles(relation):
    return sorted(
        (a['scope'], a['value'])
        for a in relation.annotations
        if a['term'] == PARTICIPANT_ROLE
    )


@pytest.mark.parametrize('dataset,source_row', [
    ('human_core_interactions', PROTEIN_PAIR),
    ('mouse_core_interactions', MOUSE_PAIR),
])
def test_ligand_and_receptor_columns_give_the_participant_roles(dataset, source_row):
    _, out = extract(source_row, dataset)
    relation, ligand, receptor = pair(out)
    assert roles(relation) == [('object', 'receptor'), ('subject', 'ligand')]
    assert (Namespace.GENESYMBOL, source_row['Ligand_Symbol']) in identifiers(ligand)
    assert (Namespace.GENESYMBOL, source_row['Receptor_Symbol']) in identifiers(receptor)


def test_complex_receptor_has_the_role_and_its_subunits_have_none():
    _, out = extract(COMPLEX_PAIR)
    relation, ligand, receptor = pair(out)
    assert receptor.entity_type == 'macromolecular_complex'
    # The source lists the complex on the receptor side, so the complex is the
    # receptor; the subunits are members of it and carry no role.
    assert roles(relation) == [('object', 'receptor'), ('subject', 'ligand')]
    members = [r for r in out.relations if r.predicate == 'has_member']
    assert len(members) == 2
    assert all(
        a['term'] != PARTICIPANT_ROLE for m in members for a in m.annotations
    )
    assert all(
        a['term'] != PARTICIPANT_ROLE
        for e in out.entities.values()
        for a in e.annotations
    )


def test_missing_subunit_symbol_does_not_shift_the_gene_ids():
    row = dict(COMPLEX_PAIR, Receptor_Symbol=',TGFBR2')
    record = resource.datasets()['human_core_interactions'].mapper(row)
    assert [
        [(i.type, i.value) for i in m.member.identifiers]
        for m in record.object.membership
    ] == [
        [(Namespace.ENTREZ, '7046')],
        [(Namespace.GENESYMBOL, 'TGFBR2'), (Namespace.ENTREZ, '7048')],
    ]


def test_complex_is_not_expanded_into_subunit_pairs():
    _, out = extract(COMPLEX_PAIR)
    interactions = [r for r in out.relations if r.predicate == 'interacts_with']
    assert len(interactions) == 1
    assert len(out.entities) == 4


def test_mouse_pair_carries_mouse_taxon_and_gene_ids():
    _, out = extract(MOUSE_PAIR, 'mouse_core_interactions')
    _, ligand, receptor = pair(out)
    assert identifiers(ligand) == {
        (Namespace.GENESYMBOL, 'Gh'),
        (Namespace.ENTREZ, '14599'),
    }
    assert identifiers(receptor) == {
        (Namespace.GENESYMBOL, 'Prlr'),
        (Namespace.ENTREZ, '19116'),
    }
    assert {ligand.taxon, receptor.taxon} == {'10090'}


def test_pair_asserts_no_direction_transcription_factor_or_pathway():
    _, out = extract(PROTEIN_PAIR)
    relation, _, _ = pair(out)
    # Only the participants' species and roles. The downstream transcription
    # factors and KEGG pathways of the pair are neither annotations nor
    # participants.
    assert sorted(
        (a['scope'], a['term'], a['value']) for a in relation.annotations
    ) == [
        ('object', PARTICIPANT_ROLE, 'receptor'),
        ('object', 'in_taxon', 'NCBITaxon:9606'),
        ('subject', PARTICIPANT_ROLE, 'ligand'),
        ('subject', 'in_taxon', 'NCBITaxon:9606'),
    ]
    values = {i['id'] for e in out.entities.values() for i in e.identifiers}
    assert not values & {'FOXO3', '2309', 'IKBKB', '3551'}
    assert len(out.entities) == 2


def test_transcription_factors_stay_in_the_evidence_payload():
    _, out = extract(PROTEIN_PAIR)
    (payload,) = out.payloads
    assert '"TF_Symbol":"FOXO3"' in payload['payload_json']
    assert '"pathway_ID":"hsa04062"' in payload['payload_json']


def test_self_pair_of_the_source_is_kept():
    _, out = extract(SELF_PAIR)
    relation, ligand, receptor = pair(out)
    assert relation.subject_entity_key == relation.object_entity_key
    assert identifiers(ligand) == {
        (Namespace.GENESYMBOL, 'OCLN'),
        (Namespace.ENTREZ, '100506658'),
    }


def test_pair_listed_in_both_orientations_keeps_each_row_orientation():
    # Claudin and nectin pairs appear as A-B and as B-A in the source; each
    # row keeps its own ligand as subject and receptor as object.
    forward = dict(PROTEIN_PAIR)
    reverse = dict(
        PROTEIN_PAIR,
        Ligand_ID=PROTEIN_PAIR['Receptor_ID'],
        Ligand_Symbol=PROTEIN_PAIR['Receptor_Symbol'],
        Receptor_ID=PROTEIN_PAIR['Ligand_ID'],
        Receptor_Symbol=PROTEIN_PAIR['Ligand_Symbol'],
    )
    _, out_forward = extract(forward)
    _, out_reverse = extract(reverse)
    (rel_f, subj_f, obj_f) = pair(out_forward)
    (rel_r, subj_r, obj_r) = pair(out_reverse)
    assert (Namespace.GENESYMBOL, 'CCL26') in identifiers(subj_f)
    assert (Namespace.GENESYMBOL, 'CCR1') in identifiers(obj_f)
    assert (Namespace.GENESYMBOL, 'CCR1') in identifiers(subj_r)
    assert (Namespace.GENESYMBOL, 'CCL26') in identifiers(obj_r)


@pytest.mark.parametrize('role', ['Ligand', 'Receptor'])
def test_no_relation_without_a_participant(role):
    mapper = resource.datasets()['human_core_interactions'].mapper
    row = dict(PROTEIN_PAIR, **{f'{role}_ID': '', f'{role}_Symbol': ''})
    assert mapper(row) is None


def test_reader_yields_one_row_per_pair_with_its_transcription_factors():
    lines = [HEADER, *CCL26_CCR1_LINES, *TGFB1_COMPLEX_LINES]
    assert list(iter_pairs(table(lines), taxon_id='9606')) == [
        PROTEIN_PAIR,
        COMPLEX_PAIR,
    ]


def test_reader_drops_exact_repeats_and_joins_rows_that_are_apart():
    lines = [HEADER, GH_PRLR_LINE, IFNG_IFNGR1_LINE, GH_PRLR_LINE]
    rows = list(iter_pairs(table(lines), taxon_id='10090'))
    assert rows[0] == MOUSE_PAIR
    assert [(r['Ligand_Symbol'], r['Receptor_Symbol']) for r in rows] == [
        ('Gh', 'Prlr'),
        ('Ifng', 'Ifngr1'),
    ]


def test_reader_rejects_subunit_ids_that_do_not_match_the_symbols():
    line = '7040\t7046,7048\t4087\thsa04926_1\thsa04926\tTGFB1\tTGFBR1\tSMAD2'
    with pytest.raises(ValueError, match='unaligned subunits'):
        list(iter_pairs(table([HEADER, line]), taxon_id='9606'))


def test_reader_rejects_a_table_with_changed_columns():
    lines = ['Ligand_ID\tReceptor_ID\tTF_ID', '10344\t1230\t2309']
    with pytest.raises(ValueError, match='Ligand_Symbol'):
        list(iter_pairs(table(lines), taxon_id='9606'))


def test_reader_rejects_an_empty_download():
    with pytest.raises(ValueError, match='no rows'):
        list(iter_pairs(table([HEADER]), taxon_id='9606'))
