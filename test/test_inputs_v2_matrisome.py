"""Matrisome master lists: parser, mappers and extractor on source rows."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import entity_type, slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import matrisome
from pypath.inputs_v2.parsers import matrisome as parser

HUMAN_TITLE = ['Homo sapiens matrisome masterlist']
HUMAN_HEADER = [
    'Matrisome Division', 'Matrisome Category', 'Gene Symbol', 'Gene Name',
    'Synonyms', 'HGNC_IDs', 'HGNC_IDs Links', 'UniProt_IDs', 'Refseq_IDs', 'Notes',
]
# Rows copied from the human list; the HGNC ids are number cells.
HUMAN_ROWS = [
    ['Core matrisome', 'ECM Glycoproteins', 'ABI3BP',
     'ABI family, member 3 (NESH) binding protein',
     'FLJ41743|FLJ41754|NESHBP|TARSH', 17265, 17265,
     'B4DSV9:D3YTG3:Q7Z7G0', 'NP_056244.2:XP_005247340.1', 'None'],
    ['Core matrisome', 'Collagens', 'COL4A1', 'collagen, type IV, alpha 1',
     'arresten|CAV', 2202, 2202, 'P02462', 'NP_001836.2', 'Basement Membrane '],
    ['Matrisome-associated', 'ECM-affiliated Proteins', 'SEMA3A',
     'sema domain, immunoglobulin domain (Ig), short basic domain, secreted, (semaphorin) 3A',
     'Hsema-I|SEMAD|SEMAL', 10723, 10723, 'Q14563', 'NP_006073.2', 'Semaphorin'],
    ['Matrisome-associated', 'Secreted Factors', 'XCL2', 'chemokine (C motif) ligand 2',
     'SCM-1b|SCM1B|SCYC2', 10646, 10646, 'Q9UBD3', 'NP_003166.1', 'None'],
]
MOUSE_HEADER = [
    'Matrisome Division', 'Matrisome Category', 'Gene Symbol', 'Gene Name',
    'Synonyms', 'MGI_IDs', 'MGI_IDs Links', 'UniProt_IDs', 'Refseq_IDs', 'Notes',
]
# Rows copied from the mouse list; the first has no synonyms (the sheet writes the
# text "None") and the last has an MGI id with a trailing space in the source.
MOUSE_ROWS = [
    ['Core matrisome', 'ECM Glycoproteins', '5430419D17Rik',
     'RIKEN cDNA 5430419D17 gene', 'None', 'MGI:1918645', 'MGI:1918645', 'Q8BZE1',
     'NP_780375.2', 'None'],
    ['Matrisome-associated', 'ECM Regulators', 'Ctsl', 'cathepsin L',
     'AA409190|Cat L|CTSL', 'MGI:88562 ', 'MGI:88562', 'P06797', 'NP_034114.1', 'Cathepsin'],
]


def workbook(header, rows, title=HUMAN_TITLE):
    book = openpyxl.Workbook()
    sheet = book.active
    for row in (title, header, *rows):
        sheet.append(row)
    data = io.BytesIO()
    book.save(data)
    return SimpleNamespace(result=io.BytesIO(data.getvalue()))


def human(rows=HUMAN_ROWS, header=HUMAN_HEADER):
    return workbook(header, rows)


def mouse(rows=MOUSE_ROWS, header=MOUSE_HEADER):
    return workbook(header, rows, ['Mus musculus matrisome masterlist'])


class FakeDownload:
    def __init__(self, opener):
        self._opener = opener

    def open(self, **_kwargs):
        return self._opener


def extract(dataset, rows):
    out = SilverExtractor('matrisome', dataset)
    mapper = matrisome.resource.datasets()[dataset].mapper
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def terms_rows():
    return list(parser.iter_terms(human(), mouse_download=FakeDownload(mouse())))


def test_human_parser_reads_the_rows_below_the_title():
    rows = list(parser.iter_genes(human()))
    assert len(rows) == 4
    assert rows[0] == {
        'division': 'Core matrisome',
        'category': 'ECM Glycoproteins',
        'symbol': 'ABI3BP',
        'name': 'ABI family, member 3 (NESH) binding protein',
        'synonyms': ['FLJ41743', 'FLJ41754', 'NESHBP', 'TARSH'],
        'gene_id': '17265',
        'taxon': '9606',
        'uniprot': ['B4DSV9', 'D3YTG3', 'Q7Z7G0'],
        'refseq': ['NP_056244.2', 'XP_005247340.1'],
        'note': None,
    }
    # A trailing space of the comment is not part of it.
    assert rows[1]['note'] == 'Basement Membrane'


def test_mouse_parser_reads_mgi_ids_and_no_synonyms():
    rows = list(parser.iter_genes(mouse(), species='mouse'))
    assert rows[0]['gene_id'] == 'MGI:1918645'
    assert rows[0]['synonyms'] == []
    assert rows[0]['taxon'] == '10090'
    assert rows[1]['gene_id'] == 'MGI:88562'


@pytest.mark.parametrize(
    'broken',
    [
        lambda rows: [],
        lambda rows: [*rows, rows[0]],
        lambda rows: [[*rows[0][:5], 'x', 'x', *rows[0][7:]]],
        lambda rows: [[*rows[0][:2], None, *rows[0][3:]]],
        lambda rows: [[rows[0][0], None, *rows[0][2:]]],
    ],
)
def test_parser_raises_on_rows_it_cannot_read(broken):
    with pytest.raises(ValueError):
        list(parser.iter_genes(human(broken(HUMAN_ROWS))))


def test_parser_raises_when_a_column_is_missing():
    header = [*HUMAN_HEADER[:5], 'HGNC', *HUMAN_HEADER[6:]]
    with pytest.raises(ValueError, match='lacks expected columns'):
        list(parser.iter_genes(human(header=header)))


def test_terms_are_divisions_and_categories_with_their_parent():
    rows = terms_rows()
    assert [(r['kind'], r['term'], r['parent']) for r in rows] == [
        ('division', 'Core matrisome', None),
        ('division', 'Matrisome-associated', None),
        ('category', 'ECM Glycoproteins', 'Core matrisome'),
        ('category', 'Collagens', 'Core matrisome'),
        ('category', 'ECM-affiliated Proteins', 'Matrisome-associated'),
        ('category', 'Secreted Factors', 'Matrisome-associated'),
        ('category', 'ECM Regulators', 'Matrisome-associated'),
    ]
    category = matrisome.map_term(rows[3])
    assert entity_type(category.type) == 'ontology_class'
    assert category.identifiers == [
        ('matrisome_category', 'Collagens'),
        (Namespace.NAME, 'Collagens'),
    ]
    assert [
        (slot_name(r.predicate), r.object.identifier_type, r.object.identifier)
        for r in category.ontology_relations
    ] == [('subclass_of', 'matrisome_division', 'Core matrisome')]
    division = matrisome.map_term(rows[0])
    assert division.ontology_relations is None

    out = extract('terms', rows)
    assert len(out.entities) == 7 + 2
    assert {r.predicate for r in out.relations} == {'subclass_of'}
    assert len(out.relations) == 5


def test_a_category_in_two_divisions_is_an_error():
    bad = [
        HUMAN_ROWS[0],
        ['Matrisome-associated', 'ECM Glycoproteins', 'XCL2', 'x', '-', 1, 1, 'Q9UBD3', 'x', 'None'],
    ]
    with pytest.raises(ValueError, match='two divisions'):
        list(parser.iter_terms(human(bad), mouse_download=FakeDownload(mouse())))


def test_human_gene_is_associated_with_its_category():
    rows = list(parser.iter_genes(human()))
    record = matrisome.map_gene(rows[1])
    assert slot_name(record.predicate) == 'associated_with'
    assert record.subject.identifiers == [
        (Namespace.HGNC, '2202'),
        (Namespace.GENESYMBOL, 'COL4A1'),
        (Namespace.NAME, 'collagen, type IV, alpha 1'),
        (Namespace.SYNONYM, 'arresten'),
        (Namespace.SYNONYM, 'CAV'),
    ]
    assert record.object.identifier == 'Collagens'

    out = extract('human', rows)
    assert len(out.relations) == 4
    assert {r.predicate for r in out.relations} == {'associated_with'}
    genes = [e for e in out.entities.values() if e.entity_type == 'gene']
    assert {e.taxon for e in genes} == {'9606'}
    notes = sorted(
        (a['value'] for r in out.relations for a in r.annotations if a['term'] == 'matrisome:note')
    )
    assert notes == ['Basement Membrane', 'Semaphorin']


def test_mouse_gene_is_identified_by_mgi_id_and_symbol():
    rows = list(parser.iter_genes(mouse(), species='mouse'))
    record = matrisome.map_gene(rows[0])
    assert record.subject.identifiers[0].type is Namespace.MGI
    assert record.subject.identifiers == [
        (Namespace.MGI, 'MGI:1918645'),
        (Namespace.GENESYMBOL, '5430419D17Rik'),
        (Namespace.NAME, 'RIKEN cDNA 5430419D17 gene'),
    ]
    out = extract('mouse', rows)
    genes = [e for e in out.entities.values() if e.entity_type == 'gene']
    assert {e.taxon for e in genes} == {'10090'}


def test_objects_are_the_same_entities_as_the_terms():
    terms = extract('terms', terms_rows())
    genes = extract('human', list(parser.iter_genes(human())))
    mouse_genes = extract('mouse', list(parser.iter_genes(mouse(), species='mouse')))
    identity = lambda out, keys: {  # noqa: E731
        (out.entities[k].entity_type, out.entities[k].namespace, out.entities[k].identifier)
        for k in keys
    }
    term_identity = identity(terms, terms.entities)
    for out in (genes, mouse_genes):
        objects = identity(out, {r.object_entity_key for r in out.relations})
        assert objects
        assert objects <= term_identity
        assert {i[:2] for i in objects} == {('ontology_class', 'matrisome_category')}


def test_nothing_is_asserted_for_the_division_or_from_the_comment():
    genes = extract('human', list(parser.iter_genes(human())))
    objects = {genes.entities[r.object_entity_key].namespace for r in genes.relations}
    # One statement per gene, to the category. The comment is an attribute and
    # no term.
    assert objects == {'matrisome_category'}
    assert len(genes.relations) == 4
    assert not any(e.identifier == 'Semaphorin' for e in genes.entities.values())
    # The UniProt and RefSeq accessions of the list are not identifiers of the gene.
    assert not any(
        e.namespace in (Namespace.UNIPROT, Namespace.REFSEQ_PROTEIN)
        for e in genes.entities.values()
    )


def test_the_text_none_is_not_a_synonym():
    # Row of the human list (AMBN): the sheet writes "None" in the synonyms and in
    # the notes of a gene that has neither.
    ambn = ['Core matrisome', 'ECM Glycoproteins', 'AMBN',
            'ameloblastin (enamel matrix protein)', 'None', 452, 452, 'Q9NP70',
            'NP_057603.1', 'None']
    (row,) = parser.iter_genes(human([ambn]))
    assert row['synonyms'] == [] and row['note'] is None
    record = matrisome.map_gene(row)
    assert record.subject.identifiers == [
        (Namespace.HGNC, '452'),
        (Namespace.GENESYMBOL, 'AMBN'),
        (Namespace.NAME, 'ameloblastin (enamel matrix protein)'),
    ]
    out = extract('human', [row])
    assert not any(
        i['id'] == 'None' for e in out.entities.values() for i in e.identifiers
    )


def test_the_vocabulary_dataset_is_called_terms():
    assert list(matrisome.resource.datasets()) == ['terms', 'human', 'mouse']
