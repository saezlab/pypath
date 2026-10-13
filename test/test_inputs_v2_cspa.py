"""Source contract of the CSPA input."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.cspa import (
    CELL_TYPE_NAMESPACE,
    LIST_ID,
    LIST_NAMESPACE,
    resource,
)
from pypath.inputs_v2.parsers.cspa import (
    PROTEIN_COLUMNS,
    iter_detections,
    iter_human_proteins,
    iter_mouse_proteins,
    iter_terms,
)

HUMAN_HEADER = (*PROTEIN_COLUMNS, 'UniProt Cell surface')
MOUSE_HEADER = (*PROTEIN_COLUMNS, 'UniProt cell surface')
# Rows of Table A (human) and Table B (mouse) of S2_File.xlsx.
NEP = ['Human', 'P08473', '1 - high confidence', 'Neprilysin', 'NEP_HUMAN', 'CD10', 4311, 'MME', 1, 29, 1, 1, 0, 'yes']
AOC3 = ['Human', 'Q16853', '3 - unspecific', 'Membrane primary amine oxidase', 'AOC3_HUMAN', 'no', 8639, 'AOC3', 1, 3, 0, 0, 0, '#N/A']
MOUSE_NEP = ['Mouse', 'Q61391', '1 - high confidence', 'Neprilysin', 'NEP_MOUSE', 'CD10', 17380, 'Mme', 1, 13, 1, 1, 0, 'yes']
ENV1 = ['Mouse', 'P10404', '1 - high confidence', 'MLV-related proviral Env polyprotein [Cleaved into: Surface protein', 'ENV1_MOUSE', 'no', None, '0', 1, 19, 1, 2, 0, 'yes']
YA003 = ['Mouse', 'Q9D7S0', '1 - high confidence', 'Uncharacterized protein UNQ511/PRO1026 homolog', 'YA003_MOUSE', 'no', 70163, 0, 0.9999, 1, 0, 0, 1, 'yes']
HLA_B38 = ['Human', 'Q95365', '1 - high confidence', 'HLA class I histocompatibility antigen, B-38 alpha chain', '1B38_HUMAN', 'no', 0, 'HLA-B', 1, 2, 1, 1, 0, 'yes']
TVA3 = ['Human', 'P01737', '3 - unspecific', 'T-cell receptor alpha chain V region PY14', 'TVA3_HUMAN', 'no', 0, '0', 1, 1, 0, 0, 0, '#N/A']
# Rows of Table_A of S1_File.xlsx: organism, experiment tag, accession, CD,
# Entrez ID, symbol, entry name, cell types, probability, unique peptides.
ANO9_TCELLS = ['Human', 'CD4+, CD25- Tcells', 'A1A5B4', 'no', 338440, 'ANO9', 'ANO9_HUMAN', 8, 1, 2]
ANO9_A431 = ['Human', 'A-431', 'A1A5B4', 'no', 338440, 'ANO9', 'ANO9_HUMAN', 8, 0.9935, 1]
FAM171 = ['Mouse', 'J1_EB', 'A2A699', 'no', 217219, 'Fam171a2', 'F1712_MOUSE', 7, 0.9968, 1]
DETECTION_HEADER = (
    'organism', 'Experiment tag', 'ID_link', 'CD', 'ENTREZ ac', 'ENTREZ gene symbol',
    'UP_entry_name', 'count detection in different cell types', 'protein probability',
    'num unique peps', None, None, None,
)


def book(sheets):
    workbook = openpyxl.Workbook()
    workbook.active.title = 'Cover'
    for name, header, rows in sheets:
        sheet = workbook.create_sheet(name)
        sheet.append(list(header))
        for row in rows:
            sheet.append(row)
        # The mouse sheet is blank below its proteins.
        sheet.append([None] * len(header))
    buffer = io.BytesIO()
    workbook.save(buffer)
    return SimpleNamespace(result=io.BytesIO(buffer.getvalue()))


def proteins_file():
    return book(
        [('Table A', HUMAN_HEADER, [NEP, AOC3]), ('Table B', MOUSE_HEADER, [MOUSE_NEP, ENV1, YA003])]
    )


def detections_file(rows=(ANO9_TCELLS, ANO9_A431, FAM171)):
    return book([('Table_A', DETECTION_HEADER, list(rows))])


def process(dataset, record):
    out = SilverExtractor('cspa', dataset)
    out.process_record(resource.datasets()[dataset].mapper(record), record, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        out.entities[relation.subject_entity_key],
        relation,
        out.entities[relation.object_entity_key],
        out,
    )


def attributes(relation):
    return {a['term']: a['value'] for a in relation.annotations if a['scope'] == 'relation'}


def numbers(relation):
    """The Measurements of a statement: term -> (number, column of the file)."""
    return {
        a['term']: (a['quantity']['has_numeric_value'], a['quantity']['source_field'])
        for a in relation.annotations
        if a['scope'] == 'relation' and a['quantity']
    }


def test_human_protein_is_associated_with_the_atlas_with_its_category():
    nep, aoc3 = list(iter_human_proteins(proteins_file()))
    protein, relation, listed, _ = process('human_proteins', nep)
    assert relation.predicate == 'associated_with'
    assert (protein.entity_type, protein.namespace, protein.identifier, protein.taxon) == (
        'protein', 'uniprot', 'P08473', '9606'
    )
    assert {(i['ns'], i['id']) for i in protein.identifiers} == {
        (Namespace.UNIPROT, 'P08473'),
        (Namespace.UNIPROT_ENTRY, 'NEP_HUMAN'),
        (Namespace.GENESYMBOL, 'MME'),
        (Namespace.ENTREZ, '4311'),
        (Namespace.SYNONYM, 'CD10'),
    }
    assert (listed.entity_type, listed.namespace, listed.identifier) == (
        'ontology_class', LIST_NAMESPACE, LIST_ID
    )
    # Flags are 'true'/'false'; numbers are Measurements (their text value is empty).
    assert attributes(relation) == {
        'has_confidence_level': '1 - high confidence',
        'has_confidence_score': '',
        'cspa:n_cell_types': '',
        'cspa:phobius_transmembrane_predicted': 'true',
        'cspa:predicted_transmembrane_domains': '',
        'cspa:gpi': 'false',
        'cspa:uniprot_cell_surface': 'true',
    }
    assert numbers(relation) == {
        'has_confidence_score': (1.0, 'protein probability'),
        'cspa:n_cell_types': (29.0, 'count detection in different cell types'),
        'cspa:predicted_transmembrane_domains': (1.0, '# predicted TM domains'),
    }
    (score,) = [a for a in relation.annotations if a['term'] == 'has_confidence_score']
    assert score['quantity']['has_numeric_value'] == 1.0
    assert score['quantity']['source_field'] == 'protein probability'


def test_an_unspecific_protein_is_in_the_list_with_that_category_and_no_cd_name():
    _, aoc3 = list(iter_human_proteins(proteins_file()))
    protein, relation, _, _ = process('human_proteins', aoc3)
    # The list does not call the protein a cell surface protein, so the statement
    # is an association with the category on it.
    assert relation.predicate == 'associated_with'
    assert attributes(relation)['has_confidence_level'] == '3 - unspecific'
    # UniProt's "#N/A" (no annotation) is no statement.
    assert 'cspa:uniprot_cell_surface' not in attributes(relation)
    assert Namespace.SYNONYM not in {i['ns'] for i in protein.identifiers}


def test_mouse_proteins_are_mouse_and_a_zero_symbol_is_no_symbol():
    nep, env1, ya003 = list(iter_mouse_proteins(proteins_file()))
    protein, _, _, _ = process('mouse_proteins', nep)
    assert (protein.namespace, protein.identifier, protein.taxon) == ('uniprot', 'Q61391', '10090')
    for record in (env1, ya003):
        assert record['gene'] is None
    protein, _, _, _ = process('mouse_proteins', env1)
    assert {i['ns'] for i in protein.identifiers} == {Namespace.UNIPROT, Namespace.UNIPROT_ENTRY}
    assert env1['entrez'] is None


def test_zero_is_no_gene_id_and_no_symbol():
    # The sheets write 0 for a protein that has no Entrez gene or symbol.
    book_ = book([('Table A', HUMAN_HEADER, [HLA_B38, TVA3])])
    hla, tva = list(iter_human_proteins(book_))
    assert (hla['entrez'], hla['gene']) == (None, 'HLA-B')
    assert (tva['entrez'], tva['gene']) == (None, None)
    protein, _, _, _ = process('human_proteins', tva)
    assert {(i['ns'], i['id']) for i in protein.identifiers} == {
        (Namespace.UNIPROT, 'P01737'),
        (Namespace.UNIPROT_ENTRY, 'TVA3_HUMAN'),
    }


def test_parser_skips_the_blank_rows_below_the_table():
    assert len(list(iter_mouse_proteins(proteins_file()))) == 3


def test_detection_is_expression_in_a_cell_type_with_probability_and_peptides():
    rows = list(iter_detections(detections_file()))
    assert [r['cell_type'] for r in rows] == ['CD4+, CD25- Tcells', 'A-431', 'J1_EB']
    protein, relation, cell_type, _ = process('detections', rows[1])
    assert relation.predicate == 'expressed_in'
    assert (protein.namespace, protein.identifier, protein.taxon) == ('uniprot', 'A1A5B4', '9606')
    assert (cell_type.entity_type, cell_type.namespace, cell_type.identifier) == (
        'ontology_class', CELL_TYPE_NAMESPACE, 'A-431'
    )
    assert numbers(relation) == {
        'has_confidence_score': (0.9935, 'protein probability'),
        'cspa:unique_peptides': (1.0, 'num unique peps'),
    }
    protein, _, _, _ = process('detections', rows[2])
    assert protein.taxon == '10090'


def test_the_vocabulary_dataset_is_called_terms_and_the_list_is_cspa_list():
    assert 'terms' in resource.datasets()
    assert LIST_NAMESPACE == 'cspa_list'
    _, relation, listed, _ = process('human_proteins', list(iter_human_proteins(proteins_file()))[0])
    assert listed.namespace == 'cspa_list'


def test_terms_are_the_atlas_and_the_cell_types_with_names():
    terms = list(iter_terms(detections_file()))
    assert terms[0] == {'kind': 'list'}
    assert [t['cell_type'] for t in terms[1:]] == ['CD4+, CD25- Tcells', 'A-431', 'J1_EB']
    out = SilverExtractor('cspa', 'terms')
    out.process_record(resource.terms.mapper(terms[2]), terms[2], 'fixture:1', 1)
    (entity,) = out.entities.values()
    assert (entity.entity_type, entity.namespace, entity.identifier) == (
        'ontology_class', CELL_TYPE_NAMESPACE, 'A-431'
    )
    assert {(i['ns'], i['id']) for i in entity.identifiers} == {
        (CELL_TYPE_NAMESPACE, 'A-431'),
        (Namespace.NAME, 'A-431'),
    }
    out = SilverExtractor('cspa', 'terms')
    out.process_record(resource.terms.mapper(terms[0]), terms[0], 'fixture:1', 1)
    (entity,) = out.entities.values()
    assert (entity.namespace, entity.identifier) == (LIST_NAMESPACE, LIST_ID)
    assert (Namespace.NAME, 'Cell Surface Protein Atlas') in {
        (i['ns'], i['id']) for i in entity.identifiers
    }


def test_objects_have_the_type_namespace_and_identifier_of_their_terms():
    terms = list(iter_terms(detections_file()))
    for detection in iter_detections(detections_file()):
        _, _, cell_type, _ = process('detections', detection)
        term = next(t for t in terms if t.get('cell_type') == detection['cell_type'])
        out = SilverExtractor('cspa', 'terms')
        out.process_record(resource.terms.mapper(term), term, 'fixture:1', 1)
        (defined,) = out.entities.values()
        assert (defined.entity_type, defined.namespace, defined.identifier) == (
            cell_type.entity_type, cell_type.namespace, cell_type.identifier
        )


def test_cell_types_are_not_given_ontology_accessions():
    terms = list(iter_terms(detections_file()))
    for term in terms[1:]:
        entity = resource.terms.mapper(term)
        assert {i.type for i in entity.identifiers} == {CELL_TYPE_NAMESPACE, Namespace.NAME}


@pytest.mark.parametrize(
    'change',
    [
        {2: '4 - maybe'},  # unknown category
        {12: 2},  # a flag that is not 0 or 1
        {13: 'maybe'},  # unknown UniProt cell surface value
        {0: 'Mouse'},  # a mouse protein on the human sheet
    ],
)
def test_parser_raises_on_a_value_it_does_not_know(change):
    row = list(NEP)
    for index, value in change.items():
        row[index] = value
    with pytest.raises(ValueError):
        list(iter_human_proteins(book([('Table A', HUMAN_HEADER, [row])])))


def test_parser_raises_on_missing_columns_repeats_and_empty_tables():
    with pytest.raises(ValueError, match='lacks expected columns'):
        list(iter_human_proteins(book([('Table A', HUMAN_HEADER[:-1], [NEP[:-1]])])))
    with pytest.raises(ValueError, match='occurs twice'):
        list(iter_human_proteins(book([('Table A', HUMAN_HEADER, [NEP, NEP])])))
    with pytest.raises(ValueError, match='occurs twice'):
        list(iter_detections(detections_file([ANO9_A431, ANO9_A431])))
    with pytest.raises(ValueError):
        list(iter_detections(detections_file([])))
    other = list(FAM171)
    other[0] = 'Rat'
    with pytest.raises(ValueError, match='organism'):
        list(iter_detections(detections_file([other])))


def test_a_cell_type_name_is_in_one_organism_only():
    other = list(FAM171)
    other[1] = 'A-431'
    with pytest.raises(ValueError, match='two organisms'):
        list(iter_terms(detections_file([ANO9_A431, other])))


def test_unknown_kind_of_term_raises():
    with pytest.raises(ValueError, match='Unknown kind'):
        resource.terms.mapper({'kind': 'tissue'})
