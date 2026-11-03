"""Source contract of the PhosphoPOINT kinase interactome input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.parsers.phosphopoint import COLUMNS, iter_pairs
from pypath.inputs_v2.phosphopoint import resource

HEADER = 'Kinase;Entrez GeneID;Gene;Entrez GeneID;Category'
# Lines of phosphopoint.csv, one of each category.
INTERACTING_PROTEIN = 'AKT1;207;PHLPP;23239;Category 1'
INTERACTING_PHOSPHOPROTEIN = 'AKT1;207;ESR2;2100;Category 2'
SUBSTRATE = 'AKT2;208;MAP3K5;4217;Category 3'
INTERACTING_SUBSTRATE = 'AKT1;207;ESR1;2099;Category 4'
# The kinase with itself.
SELF = 'AKT1;207;AKT1;207;Category 2'
# A spreadsheet turned the symbol SEPT3 into a date.
DATE_SYMBOL = 'PRKG1;5592;03-Sep;55964;Category 4'

datasets = resource.datasets()


def table(*lines, header=HEADER):
    text = '\n'.join((header, *lines)) + '\n'
    return SimpleNamespace(result=io.BytesIO(text.encode('ascii')))


def read(line):
    (row,) = iter_pairs(table(line))
    return row


def extract(dataset, source_row):
    record = datasets[dataset].mapper(source_row)
    assert record is not None
    out = SilverExtractor('phosphopoint', dataset)
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


def test_reader_keeps_the_gene_ids_of_both_proteins():
    assert read(INTERACTING_PHOSPHOPROTEIN) == dict(
        zip(COLUMNS, ('AKT1', '207', 'ESR2', '2100', 'Category 2'), strict=True)
    )


def test_substrate_is_phosphorylated_by_the_kinase():
    relation, kinase, substrate = extract('enzyme_substrate', read(SUBSTRATE))
    assert relation.predicate == 'affects'
    assert identifiers(kinase) == {
        (Namespace.ENTREZ, '208'),
        (Namespace.GENESYMBOL, 'AKT2'),
    }
    assert identifiers(substrate) == {
        (Namespace.ENTREZ, '4217'),
        (Namespace.GENESYMBOL, 'MAP3K5'),
    }
    assert (kinase.namespace, kinase.identifier) == ('entrez', '208')
    assert {kinase.entity_type, substrate.entity_type} == {'protein'}
    assert {kinase.taxon, substrate.taxon} == {'9606'}
    # Phosphorylation goes up; nothing is said about activity, and no site.
    assert annotations(relation) == [
        ('object_aspect_qualifier', 'phosphorylation'),
        ('object_direction_qualifier', 'increased'),
        ('phosphopoint:category', 'Category 3'),
    ]
    assert substrate.molecular_form is None


@pytest.mark.parametrize(
    'line', [INTERACTING_PROTEIN, INTERACTING_PHOSPHOPROTEIN, SELF]
)
def test_interacting_protein_is_an_interaction_and_not_a_substrate(line):
    row = read(line)
    relation, _, _ = extract('interactions', row)
    assert relation.predicate == 'interacts_with'
    # The category the table states is kept: Category 1 and 2 differ.
    assert annotations(relation) == [('phosphopoint:category', row['Category'])]
    # An interacting phospho-protein is not thereby a substrate of the kinase.
    assert datasets['enzyme_substrate'].mapper(row) is None


def test_categories_stay_distinct_on_the_relations():
    interacting = extract('interactions', read(INTERACTING_PROTEIN))[0]
    phospho = extract('interactions', read(INTERACTING_PHOSPHOPROTEIN))[0]
    both = extract('interactions', read(INTERACTING_SUBSTRATE))[0]
    assert [annotations(r) for r in (interacting, phospho, both)] == [
        [('phosphopoint:category', 'Category 1')],
        [('phosphopoint:category', 'Category 2')],
        [('phosphopoint:category', 'Category 4')],
    ]
    # A pair in both datasets has its category on both relations.
    assert ('phosphopoint:category', 'Category 4') in annotations(
        extract('enzyme_substrate', read(INTERACTING_SUBSTRATE))[0]
    )


def test_substrate_without_interaction_is_not_an_interaction():
    assert datasets['interactions'].mapper(read(SUBSTRATE)) is None


def test_interacting_substrate_is_in_both_datasets():
    row = read(INTERACTING_SUBSTRATE)
    relation, kinase, substrate = extract('enzyme_substrate', row)
    assert relation.predicate == 'affects'
    assert (Namespace.ENTREZ, '2099') in identifiers(substrate)
    relation, _, _ = extract('interactions', row)
    assert relation.predicate == 'interacts_with'


def test_kinase_with_itself_is_kept():
    relation, _, _ = extract('interactions', read(SELF))
    assert relation.subject_entity_key == relation.object_entity_key


def test_symbol_turned_into_a_date_is_not_an_identifier():
    _, _, substrate = extract('enzyme_substrate', read(DATE_SYMBOL))
    assert identifiers(substrate) == {(Namespace.ENTREZ, '55964')}


@pytest.mark.parametrize('dataset', ['enzyme_substrate', 'interactions'])
def test_unknown_category_gets_no_default(dataset):
    row = dict(read(SUBSTRATE), Category='Category 5')
    with pytest.raises(ValueError, match='Category 5'):
        datasets[dataset].mapper(row)


@pytest.mark.parametrize('role', ['Kinase', 'Gene'])
def test_no_relation_without_a_participant(role):
    row = dict(read(INTERACTING_SUBSTRATE), **{role: '', f'{role} Entrez GeneID': ''})
    assert datasets['enzyme_substrate'].mapper(row) is None
    assert datasets['interactions'].mapper(row) is None


def test_reader_rejects_a_table_with_changed_columns():
    with pytest.raises(ValueError, match='unexpected header'):
        list(iter_pairs(table(SUBSTRATE, header=HEADER.replace('Gene;', 'Substrate;'))))
    with pytest.raises(ValueError, match='Line 1 has 3 fields, expected 5'):
        list(iter_pairs(table(SUBSTRATE, header='Kinase;Gene;Category')))


def test_reader_rejects_a_row_with_missing_fields():
    with pytest.raises(ValueError, match='Line 2 has 4 fields, expected 5'):
        list(iter_pairs(table('AKT2;208;MAP3K5;Category 3')))


def test_reader_rejects_an_empty_download():
    with pytest.raises(ValueError, match='no rows'):
        list(iter_pairs(table()))
