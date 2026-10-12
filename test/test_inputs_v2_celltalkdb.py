"""Source contract of the CellTalkDB ligand-receptor input."""

import io
from types import SimpleNamespace

import pandas as pd
import pytest
import rdata
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import PARTICIPANT_ROLE

from pypath.inputs_v2.celltalkdb import resource
from pypath.inputs_v2.parsers.celltalkdb import COLUMNS, iter_lr_pairs


# Rows of human_lr_pair.rds and mouse_lr_pair.rds as the parser yields them.
HUMAN_ROW = {
    'lr_pair': 'CX3CL1_CX3CR1',
    'ligand_gene_symbol': 'CX3CL1',
    'receptor_gene_symbol': 'CX3CR1',
    'ligand_gene_id': '6376',
    'receptor_gene_id': '1524',
    'ligand_ensembl_protein_id': 'ENSP00000006053',
    'receptor_ensembl_protein_id': 'ENSP00000351059',
    'ligand_ensembl_gene_id': 'ENSG00000006210',
    'receptor_ensembl_gene_id': 'ENSG00000168329',
    'evidence': '10699158,24218476',
    'taxon_id': '9606',
}
# The source has the literal text NA for this ligand's Ensembl gene.
HUMAN_ROW_NO_ENSG = {
    'lr_pair': 'THPO_MPL',
    'ligand_gene_symbol': 'THPO',
    'receptor_gene_symbol': 'MPL',
    'ligand_gene_id': '7066',
    'receptor_gene_id': '4352',
    'ligand_ensembl_protein_id': 'ENSP00000204615',
    'receptor_ensembl_protein_id': 'ENSP00000361548',
    'ligand_ensembl_gene_id': None,
    'receptor_ensembl_gene_id': 'ENSG00000117400',
    'evidence': '10495959',
    'taxon_id': '9606',
}
MOUSE_ROW = {
    'lr_pair': 'Wnt3_Fzd6',
    'ligand_gene_symbol': 'Wnt3',
    'receptor_gene_symbol': 'Fzd6',
    'ligand_gene_id': '22415',
    'receptor_gene_id': '14368',
    'ligand_ensembl_protein_id': 'ENSMUSP00000000127',
    'receptor_ensembl_protein_id': 'ENSMUSP00000022906',
    'ligand_ensembl_gene_id': 'ENSMUSG00000000125',
    'receptor_ensembl_gene_id': 'ENSMUSG00000022297',
    'evidence': '19901330',
    'taxon_id': '10090',
}


def rds(tmp_path, rows, columns):
    """An opener over an RDS file holding one data frame."""
    path = tmp_path / 'lr_pair.rds'
    rdata.write_rds(path, pd.DataFrame(rows, columns=columns, dtype=object))
    return SimpleNamespace(result=io.BytesIO(path.read_bytes()))


def extract(dataset, row):
    record = resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor('celltalkdb', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def annotations(relation, scope):
    return sorted(
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def test_datasets_are_one_per_species():
    assert set(resource.datasets()) == {
        'human_interactions',
        'mouse_interactions',
    }


def test_human_pair_is_ligand_to_receptor_with_source_identifiers():
    _, out = extract('human_interactions', HUMAN_ROW)
    (relation,) = out.relations
    ligand = out.entities[relation.subject_entity_key]
    receptor = out.entities[relation.object_entity_key]
    assert relation.predicate == 'interacts_with'
    assert relation.upstream_id == 'CX3CL1_CX3CR1'
    assert identifiers(ligand) == {
        (Namespace.GENESYMBOL, 'CX3CL1'),
        (Namespace.ENTREZ, '6376'),
        (Namespace.ENSG, 'ENSG00000006210'),
        (Namespace.ENSP, 'ENSP00000006053'),
    }
    assert identifiers(receptor) == {
        (Namespace.GENESYMBOL, 'CX3CR1'),
        (Namespace.ENTREZ, '1524'),
        (Namespace.ENSG, 'ENSG00000168329'),
        (Namespace.ENSP, 'ENSP00000351059'),
    }
    assert {ligand.entity_type, receptor.entity_type} == {'protein'}
    assert {ligand.taxon, receptor.taxon} == {'9606'}
    assert annotations(relation, 'relation') == [
        ('publications', 'PMID:10699158'),
        ('publications', 'PMID:24218476'),
    ]


def test_mouse_pair_carries_mouse_taxon_and_identifiers():
    _, out = extract('mouse_interactions', MOUSE_ROW)
    (relation,) = out.relations
    ligand = out.entities[relation.subject_entity_key]
    receptor = out.entities[relation.object_entity_key]
    assert (Namespace.GENESYMBOL, 'Wnt3') in identifiers(ligand)
    assert (Namespace.ENSP, 'ENSMUSP00000022906') in identifiers(receptor)
    assert {ligand.taxon, receptor.taxon} == {'10090'}
    assert ('publications', 'PMID:19901330') in annotations(
        relation, 'relation'
    )


def test_pair_asserts_no_direction_role_or_isoform():
    record, out = extract('human_interactions', HUMAN_ROW)
    (relation,) = out.relations
    # A ligand-receptor pair is an interaction: no effect, sign or aspect.
    assert {a['term'] for a in relation.annotations} == {
        'in_taxon',
        'publications',
        PARTICIPANT_ROLE,
    }
    assert {
        a['term'] for a in relation.annotations if a['scope'] == 'relation'
    } == {'publications'}
    # The roles are annotations of the participants, not of the relation.
    assert annotations(relation, 'subject') == [
        (PARTICIPANT_ROLE, 'ligand'),
        ('in_taxon', 'NCBITaxon:9606'),
    ]
    assert annotations(relation, 'object') == [
        (PARTICIPANT_ROLE, 'receptor'),
        ('in_taxon', 'NCBITaxon:9606'),
    ]
    # The Ensembl protein is a cross-reference of the gene symbol, not a claim
    # that this one sequence is the participant.
    assert record.subject.identifiers[0].type == Namespace.GENESYMBOL
    assert all(e.molecular_form is None for e in out.entities.values())


@pytest.mark.parametrize('dataset,row', [
    ('human_interactions', HUMAN_ROW),
    ('mouse_interactions', MOUSE_ROW),
])
def test_ligand_and_receptor_columns_give_the_participant_roles(dataset, row):
    record, out = extract(dataset, row)
    (relation,) = out.relations
    ligand = out.entities[relation.subject_entity_key]
    receptor = out.entities[relation.object_entity_key]
    assert (Namespace.GENESYMBOL, row['ligand_gene_symbol']) in identifiers(ligand)
    assert (Namespace.GENESYMBOL, row['receptor_gene_symbol']) in identifiers(receptor)
    assert [v for t, v in annotations(relation, 'subject') if t == PARTICIPANT_ROLE] == ['ligand']
    assert [v for t, v in annotations(relation, 'object') if t == PARTICIPANT_ROLE] == ['receptor']
    # The role is not on the relation and does not split a protein in two.
    assert PARTICIPANT_ROLE not in {a['term'] for a in relation.annotations if a['scope'] == 'relation'}
    assert len(out.entities) == 2


def test_missing_cross_reference_is_not_filled_in():
    _, out = extract('human_interactions', HUMAN_ROW_NO_ENSG)
    (relation,) = out.relations
    ligand = out.entities[relation.subject_entity_key]
    assert identifiers(ligand) == {
        (Namespace.GENESYMBOL, 'THPO'),
        (Namespace.ENTREZ, '7066'),
        (Namespace.ENSP, 'ENSP00000204615'),
    }


@pytest.mark.parametrize(
    'evidence,publications',
    [
        # NXPH3_NRXN2: the second token has nine digits.
        ('26198319,321961154', ['PMID:26198319']),
        # HLA-E_KLRC2: 91 is a fragment, not the 1975 paper with that PMID.
        ('9486650,91,26697006', ['PMID:26697006', 'PMID:9486650']),
        # IL12A_IL12RB2: the only token is a fragment.
        ('943050', []),
        # A leading zero or nine digits is not the form of a PubMed ID.
        ('0123456,123456789,10699158', ['PMID:10699158']),
    ],
)
def test_evidence_token_that_is_not_an_intact_pmid_is_left_out(
    evidence, publications
):
    row = dict(HUMAN_ROW, evidence=evidence)
    _, out = extract('human_interactions', row)
    assert [
        value for _, value in annotations(out.relations[0], 'relation')
    ] == publications


@pytest.mark.parametrize(
    'blank',
    [
        ('ligand_gene_symbol', 'ligand_gene_id', 'ligand_ensembl_gene_id',
         'ligand_ensembl_protein_id'),
        ('receptor_gene_symbol', 'receptor_gene_id', 'receptor_ensembl_gene_id',
         'receptor_ensembl_protein_id'),
    ],
)
def test_no_relation_without_a_participant(blank):
    row = dict(HUMAN_ROW, **{column: None for column in blank})
    assert resource.datasets()['human_interactions'].mapper(row) is None


def test_reader_keeps_every_column_and_adds_the_species(tmp_path):
    rows = [
        [HUMAN_ROW[column] for column in COLUMNS],
        [
            'THPO_MPL', 'THPO', 'MPL', '7066', '4352', 'ENSP00000204615',
            'ENSP00000361548', 'NA', 'ENSG00000117400', '10495959',
        ],
    ]
    opener = rds(tmp_path, rows, list(COLUMNS))
    assert list(iter_lr_pairs(opener, taxon_id='9606')) == [
        HUMAN_ROW,
        HUMAN_ROW_NO_ENSG,
    ]


def test_reader_rejects_a_table_with_changed_columns(tmp_path):
    opener = rds(
        tmp_path,
        [['CX3CL1_CX3CR1', 'CX3CL1', 'CX3CR1']],
        ['lr_pair', 'ligand', 'receptor'],
    )
    with pytest.raises(ValueError, match='ligand_gene_symbol'):
        list(iter_lr_pairs(opener, taxon_id='9606'))


def test_reader_rejects_an_empty_download():
    with pytest.raises(ValueError, match='no rows'):
        list(iter_lr_pairs(SimpleNamespace(result=None), taxon_id='9606'))
