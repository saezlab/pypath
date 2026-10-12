"""Source contract of the EMBRACE ligand-receptor input."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import PARTICIPANT_ROLE

from pypath.inputs_v2.embrace import resource
from pypath.inputs_v2.parsers.embrace import COLUMNS, SHEET, iter_pairs


def row(*values):
    return dict(zip(COLUMNS, values, strict=True))


# Rows of the ``expressed_pairs`` sheet of Table S11.
PAIR_ROW = row('A2m', 'Lrp1', 1, 0, 0, 0, 1, 0, 1, 0)
# Microglial Apoe and neural Ldlr: the pair the paper confirmed by experiment.
APOE_ROW = row('Apoe', 'Ldlr', 0, 0, 1, 0, 1, 0, 0, 0)
# Neither partner is expressed in any of the four cell types.
NOT_EXPRESSED_ROW = row('Guca2a', 'Gucy2c', 0, 0, 0, 0, 0, 0, 0, 0)
# A mouse gene without a human-style symbol.
RIKEN_ROW = row('1700013F07Rik', 'Plscr4', 0, 0, 0, 1, 0, 1, 0, 0)


def extract(source_row):
    record = resource.datasets()['interactions'].mapper(source_row)
    assert record is not None
    out = SilverExtractor('embrace', 'interactions')
    out.process_record(record, source_row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def test_pair_is_ligand_to_receptor_by_mouse_symbol():
    relation, ligand, receptor = extract(PAIR_ROW)
    assert relation.predicate == 'interacts_with'
    assert identifiers(ligand) == {(Namespace.GENESYMBOL, 'A2m')}
    assert identifiers(receptor) == {(Namespace.GENESYMBOL, 'Lrp1')}
    assert {ligand.entity_type, receptor.entity_type} == {'protein'}


def test_participants_keep_the_source_species():
    # The symbols are neither translated to human orthologs nor upper-cased.
    _, ligand, receptor = extract(APOE_ROW)
    assert {ligand.taxon, receptor.taxon} == {'10090'}
    assert (ligand.identifier, receptor.identifier) == ('Apoe', 'Ldlr')
    _, ligand, _ = extract(RIKEN_ROW)
    assert identifiers(ligand) == {(Namespace.GENESYMBOL, '1700013F07Rik')}


def test_pair_asserts_no_direction_role_or_cell_type():
    relation, _, _ = extract(APOE_ROW)
    # Only the participants' species and roles: the expression flags have no
    # faithful slot and a ligand-receptor pair states no effect on the
    # receptor. Nothing is annotated on the relation itself.
    assert sorted(
        (a['scope'], a['term'], a['value']) for a in relation.annotations
    ) == [
        ('object', PARTICIPANT_ROLE, 'receptor'),
        ('object', 'in_taxon', 'NCBITaxon:10090'),
        ('subject', PARTICIPANT_ROLE, 'ligand'),
        ('subject', 'in_taxon', 'NCBITaxon:10090'),
    ]


def test_ligand_and_receptor_columns_give_the_participant_roles():
    relation, ligand, receptor = extract(PAIR_ROW)
    assert (ligand.identifier, receptor.identifier) == ('A2m', 'Lrp1')
    roles = [
        (a['scope'], a['value'])
        for a in relation.annotations
        if a['term'] == PARTICIPANT_ROLE
    ]
    assert sorted(roles) == [('object', 'receptor'), ('subject', 'ligand')]


def test_expression_flags_do_not_change_the_relation():
    expressed, _, _ = extract(PAIR_ROW)
    silent, _, _ = extract(dict(PAIR_ROW, N_ligand=0, N_receptor=0, M_receptor=0))
    assert expressed.relation_key == silent.relation_key
    # A pair whose partners are expressed nowhere is still a pair of the table.
    relation, ligand, receptor = extract(NOT_EXPRESSED_ROW)
    assert (ligand.identifier, receptor.identifier) == ('Guca2a', 'Gucy2c')


@pytest.mark.parametrize('role', ['ligand', 'receptor'])
def test_no_relation_without_a_participant_symbol(role):
    mapper = resource.datasets()['interactions'].mapper
    assert mapper(dict(PAIR_ROW, **{f'{role}_symbol': None})) is None


def workbook(rows, sheet=SHEET):
    """An opener over an .xlsx file with one sheet."""
    book = openpyxl.Workbook()
    table = book.active
    table.title = sheet
    for values in rows:
        table.append(values)
    content = io.BytesIO()
    book.save(content)
    return SimpleNamespace(result=io.BytesIO(content.getvalue()))


def test_reader_yields_the_rows_of_the_sheet():
    rows = [PAIR_ROW, NOT_EXPRESSED_ROW, RIKEN_ROW]
    opener = workbook([list(COLUMNS), *(list(r.values()) for r in rows)])
    assert list(iter_pairs(opener)) == rows


def test_reader_rejects_a_table_with_changed_columns():
    opener = workbook([['ligand', 'receptor'], ['A2m', 'Lrp1']])
    with pytest.raises(ValueError, match='ligand_symbol'):
        list(iter_pairs(opener))


def test_reader_rejects_a_workbook_without_the_sheet():
    opener = workbook([list(COLUMNS)], sheet='Sheet1')
    with pytest.raises(KeyError):
        list(iter_pairs(opener))


def test_reader_rejects_an_empty_download():
    with pytest.raises(ValueError, match='no rows'):
        list(iter_pairs(SimpleNamespace(result=None)))
    with pytest.raises(ValueError, match='no rows'):
        list(iter_pairs(workbook([list(COLUMNS)])))
