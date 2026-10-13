"""Phosphatome: reader, mappers and extractor on rows of Table S1 of Chen et al. 2017."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import KnowledgeOrigin

from pypath.inputs_v2 import phosphatome
from pypath.inputs_v2.parsers import phosphatome as parser

S1 = 'S1-HumanPhosphatases'
S1_GROUPS = [None, None, 'Classification', None, None, None] + [None] * 15 + ['Substrate type', None, None, None]
S1_HEADER = [
    'Entrez Gene Symbol', 'Entrez Gene ID', 'Fold', 'Family', 'Subfamily',
    'Included in Protein Phosphatome?', 'Gene Synonym', 'Description',
    'Other Designations', 'OMIM ID', 'Chromosomal locus',
    'Andersen et al. (2001, 2004)', 'Alonso et al. (2004)', 'Cohen (2004)',
    'Tonks (2006)', 'Lu & Wang (2008)', 'Li et al. (2013) ',
    'Liberti et al. (2013)', 'Wang et al. (2013)', 'Seifried et al. (2013)',
    'Hatzihristidis et al. (2014)', 'Protein substrates?',
    'Non-protein substrates?', 'Catalytic activity?', 'Note',
]


def s1_row(symbol, entrez, fold, family, subfamily, protein, substrates, other, activity):
    row = dict.fromkeys(S1_HEADER)
    row.update(
        {
            'Entrez Gene Symbol': symbol, 'Entrez Gene ID': entrez, 'Fold': fold,
            'Family': family, 'Subfamily': subfamily,
            'Included in Protein Phosphatome?': protein,
            'Protein substrates?': substrates, 'Non-protein substrates?': other,
            'Catalytic activity?': activity,
        }
    )
    return [row[c] for c in S1_HEADER]


# Rows of S1 as the sheet has them (the Entrez Gene ID is a number).
S1_ROWS = [
    s1_row('ALPI', 248, 'AP', 'AP', 'AP', 'Yes', 'Yes', 'Yes', 'Yes'),
    s1_row('ALPPL2', 251, 'AP', 'AP', 'AP', 'Yes', 'No', 'No', 'Unknown'),
    s1_row('CDKN3', 1033, 'CC1', 'DSP', 'CDKN3', 'Yes', 'Yes', 'No', 'Yes'),
    s1_row('DUSP3', 1845, 'CC1', 'DSP', 'DSP3', 'Yes', 'Yes', 'No', 'Yes'),
    s1_row('STYX', 6815, 'CC1', 'DSP', 'STYX', 'Yes', 'No', 'No', 'Verified pseudophosphatase'),
    s1_row('PTPRQ', 374462, 'CC1', 'PTP', 'PTPRB', 'Yes', 'No', 'Yes', 'Yes'),
    s1_row('INPP4A', 3631, 'INPP4', 'INPP4', None, 'No', 'No', 'Yes', 'Yes'),
    s1_row('TMEM55A', 55529, 'TMEM55', 'TMEM55', None, 'No', 'No', 'No', 'N/A'),
]
S2_HEADER = [
    'Fold', 'No. of human phosphatases in the fold',
    'No. of human phosphatases have no known substrate',
    'No. of human phosphatases with protein substrates',
    '% of human phosphatases with protein substrates',
    '% of human phosphatases with substrate that have a protein substrate',
    'SCOP ID and name',
]
S2_ROWS = [
    ['CC1', 106, 31, 56, 0.5283018867924528, 0.7466666666666667, 'c.45: (Phosphotyrosine protein) phosphatases II'],
    ['AP', 4, 1, 2, 0.5, 0.6666666666666666, 'c.76: Alkaline phosphatase-like'],
    ['TMEM55', 2, 2, 0, 0, 'N/A', 'N/A'],
    ['INPP4', 2, 0, 0, 0, 0, 'N/A'],
]


def workbook(s1_rows=S1_ROWS, s1_header=S1_HEADER, s2_rows=S2_ROWS):
    book = openpyxl.Workbook()
    for name, rows in ((S1, [S1_GROUPS, s1_header, *s1_rows]), ('S2-Fold', [S2_HEADER, *s2_rows])):
        sheet = book.create_sheet(name)
        for row in rows:
            sheet.append(row)
    del book['Sheet']
    data = io.BytesIO()
    book.save(data)
    # An archive member is the opener's only value.
    return SimpleNamespace(result={'aag1796_Tables S1 to S23.xlsx': io.BytesIO(data.getvalue())})


def extract(dataset, rows):
    out = SilverExtractor('phosphatome', dataset)
    mapper = phosphatome.resource.datasets()[dataset].mapper
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def tree():
    return {t['name']: t for t in parser.iter_phosphatase_tree(workbook())}


def test_reader_builds_the_fold_family_subfamily_tree():
    folds = tree()
    assert list(folds) == ['AP', 'CC1', 'INPP4', 'TMEM55']
    # S2 gives the SCOP name of a fold, unless it says N/A.
    assert folds['CC1']['description'] == 'c.45: (Phosphotyrosine protein) phosphatases II'
    assert folds['INPP4']['description'] is None
    ap_family = folds['AP']['children'][0]
    assert (ap_family['level'], ap_family['name']) == ('family', 'AP')
    assert ap_family['children'][0]['members'] == [
        {'entrez': '248', 'symbol': 'ALPI'},
        {'entrez': '251', 'symbol': 'ALPPL2'},
    ]
    # A gene without a subfamily is a member of its family.
    assert folds['INPP4']['children'][0]['children'] == []
    assert folds['INPP4']['children'][0]['members'] == [{'entrez': '3631', 'symbol': 'INPP4A'}]
    dsp = next(c for c in folds['CC1']['children'] if c['name'] == 'DSP')
    assert [c['name'] for c in dsp['children']] == ['CDKN3', 'DSP3', 'STYX']


def test_reader_rejects_bad_input():
    with pytest.raises(ValueError, match='Table lacks expected columns'):
        list(parser.iter_phosphatase_tree(workbook(s1_header=['Symbol', *S1_HEADER[1:]])))
    unknown_fold = [s1_row('X', 1, 'NEW', 'F', None, 'Yes', 'No', 'No', 'Yes')]
    with pytest.raises(ValueError, match='not in sheet S2'):
        list(parser.iter_phosphatase_tree(workbook(s1_rows=unknown_fold)))
    maybe = [s1_row('X', 1, 'AP', 'AP', None, 'Maybe', 'No', 'No', 'Yes')]
    with pytest.raises(ValueError, match='Unknown Phosphatome value'):
        list(parser.iter_phosphatase_classes(workbook(s1_rows=maybe)))
    with pytest.raises(ValueError, match='Malformed Phosphatome gene row'):
        list(parser.iter_phosphatase_classes(workbook(s1_rows=[s1_row('X', 'n/a', 'AP', 'AP', None, 'Yes', 'No', 'No', 'Yes')])))


def test_families_are_nested_and_a_fold_carries_its_scop_name():
    folds = list(parser.iter_phosphatase_tree(workbook()))
    out = extract('families', folds)
    entities = {(e.namespace, e.identifier): e for e in out.entities.values()}
    # The fold AP, its family AP and its subfamily AP are three entities.
    for key in (
        ('phosphatome_fold', 'AP'),
        ('phosphatome_family', 'AP'),
        ('phosphatome_subfamily', 'AP'),
    ):
        assert entities[key].entity_type == 'gene_family'
    names = {key: e for key, e in out.entities.items()}
    edges = {
        (
            (names[r.subject_entity_key].namespace, names[r.subject_entity_key].identifier),
            (names[r.object_entity_key].namespace, names[r.object_entity_key].identifier),
        )
        for r in out.relations
    }
    assert (('phosphatome_fold', 'AP'), ('phosphatome_family', 'AP')) in edges
    assert (('phosphatome_family', 'AP'), ('phosphatome_subfamily', 'AP')) in edges
    assert (('phosphatome_subfamily', 'AP'), (Namespace.ENTREZ, '248')) in edges
    assert (('phosphatome_family', 'INPP4'), (Namespace.ENTREZ, '3631')) in edges
    # No gene is a member of a fold.
    assert not any(a == ('phosphatome_fold', 'AP') and b[0] == Namespace.ENTREZ for a, b in edges)
    record = phosphatome.map_family(folds[1])
    assert [(a.term, a.value) for a in record.annotations if 'description' in str(a.term)] == [
        (record.annotations[1].term, 'SCOP c.45: (Phosphotyrosine protein) phosphatases II')
    ]


def test_the_datasets_are_families_terms_and_class_annotations():
    assert list(phosphatome.resource.datasets()) == ['families', 'terms', 'class_annotations']


def test_terms_are_the_classes_that_state_something():
    rows = list(parser.iter_class_values(workbook()))
    terms = [phosphatome.map_term(r) for r in rows]
    stated = {(t.identifiers[0].type, t.identifiers[0].value, t.identifiers[1].value) for t in terms if t}
    assert stated == {
        ('phosphatome_set', 'Protein phosphatome', 'Protein phosphatome'),
        ('phosphatome_substrate_type', 'Protein substrates', 'Protein substrates'),
        ('phosphatome_substrate_type', 'Non-protein substrates', 'Non-protein substrates'),
        ('phosphatome_catalytic_activity', 'Yes', 'Catalytically active'),
        ('phosphatome_catalytic_activity', 'Verified pseudophosphatase', 'Verified pseudophosphatase'),
    }
    # A "No", an "Unknown" and an "N/A" are not classes.
    assert [(r['column'], r['value']) for r, t in zip(rows, terms) if t is None] == [
        ('Protein substrates?', 'No'),
        ('Non-protein substrates?', 'No'),
        ('Catalytic activity?', 'Unknown'),
        ('Included in Protein Phosphatome?', 'No'),
        ('Catalytic activity?', 'N/A'),
    ]
    with pytest.raises(ValueError, match='Unknown Phosphatome value'):
        phosphatome.map_term({'column': 'Catalytic activity?', 'value': 'Maybe'})


def test_genes_are_associated_with_the_classes_and_the_objects_are_the_terms():
    classes = list(parser.iter_phosphatase_classes(workbook()))
    out = extract('class_annotations', classes)
    terms = extract('terms', parser.iter_class_values(workbook()))
    term_keys = {(e.entity_type, e.namespace, e.identifier) for e in terms.entities.values()}
    statements = {}
    for r in out.relations:
        subject = out.entities[r.subject_entity_key]
        obj = out.entities[r.object_entity_key]
        assert r.predicate == 'associated_with'
        # The object of a statement is the term of the terms dataset.
        assert (obj.entity_type, obj.namespace, obj.identifier) in term_keys
        statements.setdefault(subject.identifier, set()).add(obj.identifier)
    assert statements['248'] == {
        'Protein phosphatome', 'Protein substrates', 'Non-protein substrates', 'Yes'
    }
    assert statements['6815'] == {'Protein phosphatome', 'Verified pseudophosphatase'}
    assert statements['3631'] == {'Non-protein substrates', 'Yes'}
    # "Unknown" catalytic activity is no statement; ALPPL2 is only a protein phosphatase.
    assert statements['251'] == {'Protein phosphatome'}
    # TMEM55A is in no class at all: no record.
    assert '55529' not in statements
    assert phosphatome.map_classes(classes[-1]) is None


def test_knowledge_origin_is_computationally_predicted():
    assert tuple(phosphatome.config.knowledge_origin) == (KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,)
