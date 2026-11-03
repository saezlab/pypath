"""OncoTree terms: codes, names, hierarchy, tissue, main type, cross-references, history."""

import copy
import io
import json
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import intogen, oncotree
from pypath.inputs_v2.parsers.oncotree import iter_terms

# A cut of https://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2025_10_03: the
# root, the lung branch down to LUAD, SPCC (a code with a previous one), and a
# node with two references to one database.
TREE = {
    'TISSUE': {
        'code': 'TISSUE', 'color': None, 'name': 'Tissue', 'mainType': None,
        'externalReferences': {'UMLS': ['C0040300'], 'NCI': ['C12801']},
        'tissue': None, 'parent': None, 'history': [], 'level': 0,
        'revocations': [], 'precursors': [],
        'children': {
            'LUNG': {
                'code': 'LUNG', 'color': 'Blue', 'name': 'Lung', 'mainType': 'Non-Small Cell Lung Cancer',
                'externalReferences': {'UMLS': ['C0024109'], 'NCI': ['C12468']},
                'tissue': 'Lung', 'parent': 'TISSUE', 'history': [], 'level': 1,
                'revocations': [], 'precursors': [],
                'children': {
                    'NSCLC': {
                        'code': 'NSCLC', 'color': 'Gainsboro', 'name': 'Non-Small Cell Lung Cancer',
                        'mainType': 'Non-Small Cell Lung Cancer',
                        'externalReferences': {'UMLS': ['C0007131'], 'NCI': ['C2926']},
                        'tissue': 'Lung', 'parent': 'LUNG', 'history': [], 'level': 2,
                        'revocations': [], 'precursors': [],
                        'children': {
                            'LUAD': {
                                'code': 'LUAD', 'color': 'Gainsboro', 'name': 'Lung Adenocarcinoma',
                                'mainType': 'Non-Small Cell Lung Cancer',
                                'externalReferences': {'UMLS': ['C0152013'], 'NCI': ['C3512']},
                                'tissue': 'Lung', 'parent': 'NSCLC', 'history': [], 'level': 3,
                                'revocations': [], 'precursors': [], 'children': {},
                            },
                            'SPCC': {
                                'code': 'SPCC', 'color': 'Gainsboro',
                                'name': 'Spindle Cell Carcinoma of the Lung',
                                'mainType': 'Non-Small Cell Lung Cancer',
                                'externalReferences': {'UMLS': ['C1708784'], 'NCI': ['C45541']},
                                'tissue': 'Lung', 'parent': 'NSCLC', 'history': ['SpCC'], 'level': 3,
                                'revocations': [], 'precursors': [], 'children': {},
                            },
                            'SRCCR': {
                                'code': 'SRCCR', 'color': 'Gainsboro', 'name': 'Signet Ring Cell Carcinoma',
                                'mainType': 'Non-Small Cell Lung Cancer',
                                'externalReferences': {'UMLS': ['C0279654', 'C1707436'], 'NCI': ['C9168']},
                                'tissue': 'Lung', 'parent': 'NSCLC', 'history': [], 'level': 3,
                                'revocations': [], 'precursors': [], 'children': {},
                            },
                            'NOREF': {
                                'code': 'NOREF', 'color': 'Gainsboro', 'name': 'A node without references',
                                'mainType': 'Non-Small Cell Lung Cancer', 'externalReferences': {},
                                'tissue': 'Lung', 'parent': 'NSCLC', 'history': ['an old name'], 'level': 3,
                                'revocations': [], 'precursors': [], 'children': {},
                            },
                        },
                    },
                },
            },
        },
    },
}


def opener(tree):
    return SimpleNamespace(result=io.StringIO(json.dumps(tree)))


def terms(tree=TREE):
    return {row['code']: row for row in iter_terms(opener(tree))}


def extract(row):
    record = oncotree.map_term(row)
    extractor = SilverExtractor('oncotree', 'terms')
    extractor.process_record(record, raw_payload={}, row_id='1', row_number=1)
    return record, extractor


def term_entity(extractor, identifier):
    (entity,) = (e for e in extractor.entities.values() if e.identifier == identifier and e.label)
    return entity


def test_parser_yields_each_node_after_its_parent_with_the_root():
    rows = list(iter_terms(opener(TREE)))
    assert [row['code'] for row in rows] == ['TISSUE', 'LUNG', 'NSCLC', 'LUAD', 'SPCC', 'SRCCR', 'NOREF']
    assert rows[0]['parent'] is None
    assert {row['code']: row['parent'] for row in rows[1:]} == {
        'LUNG': 'TISSUE', 'NSCLC': 'LUNG', 'LUAD': 'NSCLC', 'SPCC': 'NSCLC', 'SRCCR': 'NSCLC', 'NOREF': 'NSCLC',
    }
    assert rows[0]['main_type'] is None and rows[0]['tissue'] is None


def test_term_has_code_name_main_type_tissue_references_and_parent():
    row = terms()['LUAD']
    record, extracted = extract(row)
    entity = term_entity(extracted, 'LUAD')
    # the code as written, in the namespace IntOGen cites it in
    assert (entity.entity_type, entity.namespace) == ('ontology_class', 'oncotree')
    assert entity.label == 'Lung Adenocarcinoma'
    assert {(i['ns'], i['id']) for i in entity.identifiers} == {
        ('oncotree', 'LUAD'), ('name', 'Lung Adenocarcinoma'),
    }
    annotations = {(format_term(a.term), a.value) for a in record.annotations}
    assert annotations == {
        ('xref', 'UMLS:C0152013'),
        ('xref', 'NCIT:C3512'),
        ('oncotree:main_type', 'Non-Small Cell Lung Cancer'),
        ('oncotree:tissue', 'Lung'),
    }
    (relation,) = extracted.relations
    assert relation.predicate == 'subclass_of'
    assert extracted.entities[relation.object_entity_key].identifier == 'NSCLC'
    assert {e.namespace for e in extracted.entities.values()} == {'oncotree'}


def test_root_has_no_parent_and_no_main_type_or_tissue():
    record, extracted = extract(terms()['TISSUE'])
    assert term_entity(extracted, 'TISSUE').label == 'Tissue'
    assert not extracted.relations
    assert {(format_term(a.term), a.value) for a in record.annotations} == {
        ('xref', 'UMLS:C0040300'), ('xref', 'NCIT:C12801'),
    }


def test_previous_code_is_an_identifier_of_the_term():
    _, extracted = extract(terms()['SPCC'])
    entity = term_entity(extracted, 'SPCC')
    assert {(i['ns'], i['id']) for i in entity.identifiers} == {
        ('oncotree', 'SPCC'), ('oncotree', 'SpCC'), ('name', 'Spindle Cell Carcinoma of the Lung'),
    }


def test_history_that_is_not_a_plain_code_is_not_an_identifier():
    _, extracted = extract(terms()['NOREF'])
    entity = term_entity(extracted, 'NOREF')
    assert {(i['ns'], i['id']) for i in entity.identifiers} == {
        ('oncotree', 'NOREF'), ('name', 'A node without references'),
    }
    assert not entity.annotations or not any(a['term'] == 'xref' for a in entity.annotations)


def test_several_references_to_one_database_are_all_kept():
    record = oncotree.map_term(terms()['SRCCR'])
    assert {a.value for a in record.annotations if format_term(a.term) == 'xref'} == {
        'UMLS:C0279654', 'UMLS:C1707436', 'NCIT:C9168',
    }


def test_term_is_the_entity_intogen_cites_by_code():
    row = terms()['LUAD']
    cited = SilverExtractor('intogen', 'drivers')
    cited.process_record(
        intogen.cancer_type_object({'cancer_type': 'LUAD'}), raw_payload={}, row_id='1', row_number=1,
    )
    (ours,) = (e for e in extract(row)[1].entities.values() if e.identifier == 'LUAD')
    (theirs,) = cited.entities.values()
    assert (ours.entity_type, ours.namespace, ours.identifier) == (
        theirs.entity_type, theirs.namespace, theirs.identifier,
    )
    assert ours.namespace == Namespace.ONCOTREE.value


def test_parent_that_is_not_the_node_it_is_under_raises():
    tree = copy.deepcopy(TREE)
    tree['TISSUE']['children']['LUNG']['children']['NSCLC']['parent'] = 'BREAST'
    with pytest.raises(ValueError, match='NSCLC'):
        terms(tree)


def test_node_whose_parent_is_absent_raises():
    tree = copy.deepcopy(TREE)
    tree['LUNG'] = tree['TISSUE']['children'].pop('LUNG')  # LUNG at the top, still naming TISSUE
    with pytest.raises(ValueError, match='LUNG'):
        terms(tree)


def test_repeated_code_raises():
    tree = copy.deepcopy(TREE)
    nscl = tree['TISSUE']['children']['LUNG']['children']['NSCLC']['children']
    nscl['LUNG'] = {**nscl['LUAD'], 'code': 'LUNG', 'history': []}
    with pytest.raises(ValueError, match='repeated'):
        terms(tree)


def test_node_without_a_name_raises():
    tree = copy.deepcopy(TREE)
    tree['TISSUE']['children']['LUNG']['name'] = ''
    with pytest.raises(ValueError, match='no name'):
        terms(tree)


def test_previous_code_that_is_a_live_code_raises():
    tree = copy.deepcopy(TREE)
    tree['TISSUE']['children']['LUNG']['children']['NSCLC']['children']['SPCC']['history'] = ['LUAD']
    with pytest.raises(ValueError, match='LUAD'):
        terms(tree)


def test_unknown_reference_database_raises():
    tree = copy.deepcopy(TREE)
    tree['TISSUE']['externalReferences'] = {'MESH': ['D000001']}
    with pytest.raises(ValueError, match='MESH'):
        terms(tree)


def test_empty_download_and_empty_record_give_nothing():
    with pytest.raises(ValueError, match='empty'):
        list(iter_terms(SimpleNamespace(result=io.StringIO(''))))
    assert oncotree.map_term({}) is None


def test_resource_declares_the_classification():
    assert oncotree.config.id == ResourceCv.ONCOTREE
    assert oncotree.config.license == LicenseCV.CC_BY_4_0
    assert oncotree.config.pubmed == '33625877'
    assert oncotree.resource.terms.mapper is oncotree.map_term
    assert 'version=oncotree_2025_10_03' in oncotree.resource.terms.download.url
    assert Namespace.ONCOTREE.value == 'oncotree'
