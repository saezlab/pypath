"""Source contract of the Surfaceome input."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.parsers.surfaceome import COLUMNS, SHEET, iter_proteins, iter_terms
from pypath.inputs_v2.surfaceome import LIST_ID, LIST_NAMESPACE, resource

# Rows of the sheet "in silico surfaceome only" of table S3 of Bausch-Fluck et
# al. 2018, in the order of COLUMNS: entry name, accession, gene, GeneID, label,
# label source, training set, score, FPR class, TM domains, signal peptide,
# topology, topology source, Almen main class, Almen sub-class.
S12A8 = ['S12A8_HUMAN', 'A0AV02', 'SLC12A8', 84561, 'surface', 'machine learning', None, 0.8363, 1, 13, 0, 'NC:1-37;TM:38-60;CY:61-71;TM:72-93;NC:94-98;TM:99-116;CY:117-122;TM:123-142;NC:143-153;TM:154-174;CY:175-183;TM:184-204;NC:205-232;TM:233-254;CY:255-266;TM:267-289;NC:290-307;TM:308-330;CY:331-359;TM:360-377;NC:378-382;TM:383-403;CY:404-592;TM:593-616;NC:617-621;TM:622-646;CY:647-714', 'phobius', 'Transporters', 'SLC;APC;SLC12']
ESYT3 = ['ESYT3_HUMAN', 'A0FGR9', 'ESYT3', 83850, 'surface', 'pos. trainingset', 'surface', 0.0465, None, 2, 0, 'CY:1-27;TM:28-46;NC:47-51;TM:52-72;CY:73-886', 'uniprot', 'Unclassified', 'Unclassified']
JUNO = ['JUNO_HUMAN', 'A6ND01', 'IZUMO1R', 390243, 'surface', 'GPI (UniProt)', None, 0.6267, '2 (nonTM)', 0, 1, 'SP:1-19;NC:20-250', 'uniprot', 'Unclassified', 'Unclassified']
TM11D = ['TM11D_HUMAN', 'O60235', 'TMPRSS11D', 9407, 'surface', 'pos. trainingset', 'surface', 0.546, 3, 1, 0, 'CY:1-17;TM:18-42;NC:43-418', 'uniprot', 'Enzymes', 3.4]
QSOX1 = ['QSOX1_HUMAN', 'O00391', 'QSOX1', 5768, 'surface', 'machine learning', None, 0.6248, 2, 1, 1, 'SP:1-29;NC:30-705;TM:706-731;CY:732-747', 'uniprot - CSPA based corrected Nterm - CSPA confirmed', 'Enzymes', '1.8.3']
T131L = ['T131L_HUMAN', 'A2VDJ0', 'KIAA0922', 23240, 'surface', 'machine learning', None, 0.7784, 1, 1, 1, 'SP:1-40;NC:41-869;TM:870-891;CY:892-1609', 'uniprot - CSPA confirmed', None, None]
UPK3L = ['UPK3L_HUMAN', 'B0FP48', 'UPK3BL', '100134938;102724187', 'surface', 'machine learning', None, 0.8842, 1, 1, 1, 'SP:1-33;NC:34-202;TM:203-228;CY:229-263', 'uniprot', 'Miscellaneous', 'StructuralAndAdhesion;UPK3']
GSG1M = ['GSG1M_HUMAN', 'A8MUP6', None, None, 'surface', 'machine learning', None, 0.6747, 2, 3, 1, 'SP:1-28;NC:29-120;TM:121-144;CY:145-158;TM:159-184;NC:185-203;TM:204-225;CY:226-235', 'phobius', 'Unclassified', 'Unclassified']
ROWS = [S12A8, ESYT3, JUNO, TM11D, QSOX1, T131L, UPK3L, GSG1M]


def workbook(rows, header=COLUMNS, sheet=SHEET):
    book = openpyxl.Workbook()
    book.active.title = 'SurfaceomeMasterTable'
    table = book.create_sheet(sheet)
    table.append(['table S3: List of predicted surfaceome proteins'])
    table.append(list(header))
    for row in rows:
        table.append(row)
    buffer = io.BytesIO()
    book.save(buffer)
    return SimpleNamespace(result=io.BytesIO(buffer.getvalue()))


def records(rows):
    return list(iter_proteins(workbook(rows)))


def statements(row):
    """The statements of a protein: (object entity, relation, subject)."""
    result = []
    for record in records([row]):
        out = SilverExtractor('surfaceome', 'proteins')
        out.process_record(resource.proteins.mapper(record), record, 'fixture:1', 1)
        (relation,) = out.relations
        result.append(
            (
                out.entities[relation.object_entity_key],
                relation,
                out.entities[relation.subject_entity_key],
            )
        )
    return result


def attributes(relation):
    return {a['term']: a['value'] for a in relation.annotations if a['scope'] == 'relation'}


def numbers(relation):
    """The Measurements of a statement: term -> (number, column of the file)."""
    return {
        a['term']: (a['quantity']['has_numeric_value'], a['quantity']['source_field'])
        for a in relation.annotations
        if a['scope'] == 'relation' and a['quantity']
    }


def test_protein_is_associated_with_the_surfaceome_with_its_prediction():
    (listed, relation, protein), (almen, almen_relation, _) = statements(S12A8)
    assert relation.predicate == 'associated_with'
    assert (listed.entity_type, listed.namespace, listed.identifier) == (
        'ontology_class', LIST_NAMESPACE, LIST_ID
    )
    assert (protein.entity_type, protein.namespace, protein.identifier, protein.taxon) == (
        'protein', 'uniprot', 'A0AV02', '9606'
    )
    assert {(i['ns'], i['id']) for i in protein.identifiers} == {
        (Namespace.UNIPROT, 'A0AV02'),
        (Namespace.UNIPROT_ENTRY, 'S12A8_HUMAN'),
        (Namespace.GENESYMBOL, 'SLC12A8'),
        (Namespace.ENTREZ, '84561'),
    }
    terms = attributes(relation)
    assert terms.pop('surfaceome:topology').startswith('NC:1-37;TM:38-60;')
    # The flag is 'false'; the score and the number of domains are Measurements
    # (their text value is empty).
    assert terms == {
        'surfaceome:score': '',
        'surfaceome:label_source': 'machine learning',
        'surfaceome:transmembrane_domains': '',
        'surfaceome:signal_peptide': 'false',
        'surfaceome:topology_source': 'phobius',
        'has_confidence_level': 'class 1: false positive rate 1%',
    }
    assert numbers(relation) == {
        'surfaceome:score': (0.8363, 'MachineLearning score'),
        'surfaceome:transmembrane_domains': (13.0, 'TM domains'),
    }
    # The class of Almen et al. is a second statement without attributes.
    assert (almen.namespace, almen.identifier) == ('almen2009_class', 'SLC;APC;SLC12')
    assert almen_relation.predicate == 'associated_with'
    assert attributes(almen_relation) == {}


def test_the_vocabulary_dataset_is_called_terms_and_the_list_is_surfaceome_list():
    assert 'terms' in resource.datasets()
    assert LIST_NAMESPACE == 'surfaceome_list'
    (listed, _, _), _ = statements(S12A8)
    assert listed.namespace == 'surfaceome_list'


def test_training_set_and_gpi_rows_keep_how_they_got_their_label():
    (_, training, _), _ = statements(ESYT3)
    assert attributes(training)['surfaceome:label_source'] == 'pos. trainingset'
    assert attributes(training)['surfaceome:training_set'] == 'surface'
    # No false positive rate class: no confidence level.
    assert 'has_confidence_level' not in attributes(training)
    (_, gpi, _), _ = statements(JUNO)
    assert attributes(gpi)['surfaceome:label_source'] == 'GPI (UniProt)'
    assert attributes(gpi)['surfaceome:signal_peptide'] == 'true'
    assert attributes(gpi)['has_confidence_level'] == 'class 2 (nonTM): false positive rate 5%'


def test_unclassified_proteins_are_in_the_main_class_unclassified():
    _, (almen, _, _) = statements(ESYT3)
    assert (almen.namespace, almen.identifier) == ('almen2009_class', 'Unclassified')


def test_enzymes_are_in_the_main_class_and_their_ec_classes():
    # Excel has stored an EC class as a number (3.4); another is text.
    got = {(t.namespace, t.identifier) for t, _, _ in statements(TM11D)[1:]}
    assert got == {('almen2009_class', 'Enzymes'), ('ec', 'EC:3.4')}
    got = {(t.namespace, t.identifier) for t, _, _ in statements(QSOX1)[1:]}
    assert got == {('almen2009_class', 'Enzymes'), ('ec', 'EC:1.8.3')}


def test_a_protein_without_a_class_has_only_the_list_statement():
    assert len(statements(T131L)) == 1


def test_several_entrez_ids_are_all_kept():
    ((_, _, protein), _) = statements(UPK3L)
    assert {i['id'] for i in protein.identifiers if i['ns'] == Namespace.ENTREZ} == {
        '100134938', '102724187'
    }


def test_a_protein_without_gene_symbol_or_entrez_has_the_accession_and_entry_name():
    ((_, _, protein), _) = statements(GSG1M)
    assert {(i['ns'], i['id']) for i in protein.identifiers} == {
        (Namespace.UNIPROT, 'A8MUP6'),
        (Namespace.UNIPROT_ENTRY, 'GSG1M_HUMAN'),
    }


def test_terms_have_the_list_and_the_classes_with_names_and_parents():
    terms = list(iter_terms(workbook(ROWS)))
    assert terms[0] == {'kind': 'surfaceome'}
    classes = {t['id']: t for t in terms if t['kind'] == 'class'}
    assert classes['StructuralAndAdhesion;UPK3']['parent'] == 'StructuralAndAdhesion'
    assert classes['StructuralAndAdhesion']['parent'] == 'Miscellaneous'
    assert classes['SLC;APC;SLC12']['name'] == 'SLC > APC > SLC12'
    listed = resource.terms.mapper(terms[0])
    out = SilverExtractor('surfaceome', 'terms')
    out.process_record(listed, terms[0], 'fixture:1', 1)
    (entity,) = out.entities.values()
    assert (entity.entity_type, entity.namespace, entity.identifier) == (
        'ontology_class', LIST_NAMESPACE, LIST_ID
    )
    assert (Namespace.NAME, 'In silico human surfaceome') in {
        (i['ns'], i['id']) for i in entity.identifiers
    }


def test_object_of_each_statement_is_a_term_of_the_terms_dataset():
    terms = list(iter_terms(workbook(ROWS)))
    defined = {LIST_ID} | {t['id'] for t in terms if t['kind'] == 'class'}
    for row in ROWS:
        for obj, _, _ in statements(row):
            if obj.namespace in (LIST_NAMESPACE, 'almen2009_class'):
                assert obj.identifier in defined
                assert obj.entity_type == 'ontology_class'


@pytest.mark.parametrize(
    'change',
    [
        {4: 'nonsurface'},  # a protein the sheet does not list as surface
        {5: 'guess'},  # unknown label source
        {10: 2},  # unknown signal peptide flag
        {13: None},  # a sub-class without a main class
        {8: '4'},  # unknown false positive rate class
    ],
)
def test_parser_and_mapper_raise_on_a_value_they_do_not_know(change):
    row = list(S12A8)
    for index, value in change.items():
        row[index] = value
    with pytest.raises(ValueError):
        for record in records([row]):
            resource.proteins.mapper(record)


def test_parser_raises_on_a_repeated_accession_and_on_missing_columns():
    with pytest.raises(ValueError, match='occurs twice'):
        records([S12A8, S12A8])
    with pytest.raises(ValueError, match='lacks expected columns'):
        list(iter_proteins(workbook([S12A8], header=COLUMNS[:-1])))
    with pytest.raises(ValueError):
        records([])


def test_classes_are_mapped_by_the_shared_module_and_equal_those_of_almen2009():
    from pypath.inputs_v2 import _almen_classes, almen2009, surfaceome

    # No import from one resource module to the other.
    assert surfaceome.map_class is _almen_classes.map_class
    assert almen2009.map_class is _almen_classes.map_class
    assert surfaceome.term_entity is _almen_classes.term_entity
    # The same class row gives the same entity in both resources.
    for row in (
        {'id': 'Enzymes', 'name': 'Enzymes', 'label': 'Enzymes', 'parent': None},
        {
            'id': 'SLC;APC;SLC12', 'name': 'SLC > APC > SLC12', 'label': 'SLC12',
            'parent': 'SLC;APC',
        },
    ):
        assert surfaceome.resource.terms.mapper({**row, 'kind': 'class'}) == (
            almen2009.resource.terms.mapper(row)
        )

