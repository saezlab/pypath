"""Source contract of the iTALK ligand-receptor input."""

import io
from types import SimpleNamespace

import pandas as pd
import pytest
import rdata
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import PARTICIPANT_ROLE

from pypath.inputs_v2.italk import resource
from pypath.inputs_v2.parsers.italk import COLUMNS, iter_pairs


def row(*values):
    return dict(zip(COLUMNS, values, strict=True))


# Rows of the ``database`` data frame in LR_database.rda.
GROWTH_FACTOR_ROW = row(
    'EGF_EGFR', 'EGF', 'epidermal growth factor', 'EGFR',
    'epidermal growth factor receptor', 'growth factor',
)
# The pair name and the participant names use aliases here.
CHECKPOINT_ROW = row(
    'PD-L1_PD-1', 'CD274', 'PD-L1', 'PDCD1', 'PD-1', 'checkpoint'
)
NO_NAME_ROW = row('KRTAP4-12_DLK2', 'KRTAP4-12', None, 'DLK2', None, 'other')
# A viral receptor: the source has no approved symbol for it.
VIRAL_RECEPTOR_ROW = row('CCL2_US28', 'CCL2', 'CCL2', None, 'US28', 'cytokine')


def extract(source_row):
    record = resource.datasets()['interactions'].mapper(source_row)
    assert record is not None
    out = SilverExtractor('italk', 'interactions')
    out.process_record(record, source_row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def test_pair_is_ligand_to_receptor_by_approved_symbol():
    relation, ligand, receptor = extract(GROWTH_FACTOR_ROW)
    assert relation.predicate == 'interacts_with'
    assert relation.upstream_id == 'EGF_EGFR'
    assert identifiers(ligand) == {
        (Namespace.GENESYMBOL, 'EGF'),
        (Namespace.NAME, 'epidermal growth factor'),
    }
    assert identifiers(receptor) == {
        (Namespace.GENESYMBOL, 'EGFR'),
        (Namespace.NAME, 'epidermal growth factor receptor'),
    }
    assert {ligand.entity_type, receptor.entity_type} == {'protein'}
    assert {ligand.taxon, receptor.taxon} == {'9606'}


def test_alias_is_a_name_and_the_approved_symbol_identifies_the_protein():
    relation, ligand, receptor = extract(CHECKPOINT_ROW)
    assert relation.upstream_id == 'PD-L1_PD-1'
    assert (ligand.namespace, ligand.identifier) == ('genesymbol', 'CD274')
    assert (receptor.namespace, receptor.identifier) == ('genesymbol', 'PDCD1')
    assert (Namespace.NAME, 'PD-L1') in identifiers(ligand)
    assert (Namespace.NAME, 'PD-1') in identifiers(receptor)


def test_pair_asserts_no_direction_role_or_classification():
    relation, _, _ = extract(CHECKPOINT_ROW)
    # Only the participants' species and roles; the classification has no
    # faithful slot and a ligand-receptor pair states no effect on the
    # receptor. Nothing is annotated on the relation itself.
    assert sorted(
        (a['scope'], a['term'], a['value']) for a in relation.annotations
    ) == [
        ('object', PARTICIPANT_ROLE, 'receptor'),
        ('object', 'in_taxon', 'NCBITaxon:9606'),
        ('subject', PARTICIPANT_ROLE, 'ligand'),
        ('subject', 'in_taxon', 'NCBITaxon:9606'),
    ]


def test_ligand_and_receptor_columns_give_the_participant_roles():
    relation, ligand, receptor = extract(GROWTH_FACTOR_ROW)
    assert (ligand.identifier, receptor.identifier) == ('EGF', 'EGFR')
    roles = [
        (a['scope'], a['value'])
        for a in relation.annotations
        if a['term'] == PARTICIPANT_ROLE
    ]
    assert sorted(roles) == [('object', 'receptor'), ('subject', 'ligand')]


def test_missing_name_adds_no_identifier():
    _, ligand, receptor = extract(NO_NAME_ROW)
    assert identifiers(ligand) == {(Namespace.GENESYMBOL, 'KRTAP4-12')}
    assert identifiers(receptor) == {(Namespace.GENESYMBOL, 'DLK2')}


def test_no_relation_for_a_receptor_without_an_approved_symbol():
    mapper = resource.datasets()['interactions'].mapper
    # The alias alone (``US28``) must not become a human protein.
    assert mapper(VIRAL_RECEPTOR_ROW) is None


def test_no_relation_without_a_ligand_symbol():
    mapper = resource.datasets()['interactions'].mapper
    assert mapper(dict(GROWTH_FACTOR_ROW, **{'Ligand.ApprovedSymbol': None})) is None


def rda(tmp_path, rows, columns, name='database'):
    """An opener over an .rda file holding one named data frame."""
    path = tmp_path / 'LR_database.rda'
    frame = pd.DataFrame(rows, columns=columns, dtype=object)
    rdata.write_rda(path, {name: frame})
    return SimpleNamespace(result=io.BytesIO(path.read_bytes()))


def test_reader_yields_the_rows_of_the_database_object(tmp_path):
    rows = [GROWTH_FACTOR_ROW, VIRAL_RECEPTOR_ROW]
    opener = rda(
        tmp_path, [list(r.values()) for r in rows], list(COLUMNS)
    )
    assert list(iter_pairs(opener)) == rows


def test_reader_rejects_a_table_with_changed_columns(tmp_path):
    opener = rda(tmp_path, [['EGF', 'EGFR']], ['ligand', 'receptor'])
    with pytest.raises(ValueError, match='Pair.Name'):
        list(iter_pairs(opener))


def test_reader_rejects_an_empty_download():
    with pytest.raises(ValueError, match='no rows'):
        list(iter_pairs(SimpleNamespace(result=None)))
