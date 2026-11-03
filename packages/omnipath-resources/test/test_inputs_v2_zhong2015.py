"""Source contract of the Zhong2015 input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.parsers.zhong2015 import iter_cam_types, iter_genes
from pypath.inputs_v2.zhong2015 import CAM_TYPES, TYPE_NAMESPACE, resource

# Rows of Supplementary Table S1: the file has four more empty columns and
# trailing spaces in the names.
HEADER = 'symbol\tentrezID\tCAM type\tInter Pro\tGO\tNCBI\tSum\tname\t\t\t\t'
CDH1 = 'CDH1\t999\ti\t1\t1\t1\t3\tcadherin 1\t\t\t\t'
CDH13 = 'CDH13\t1012\ti\t1\t1\t1\t3\tcadherin 13 \t\t\t\t'
ITGA6 = 'ITGA6\t3655\tm\t1\t0\t1\t2\tintegrin alpha 6 \t\t\t\t'
ESAM = 'ESAM\t90952\ttj\t0\t0\t1\t1\tendothelial cell adhesion molecule \t\t\t\t'
CD24 = 'CD24\t100000000\tc\t0\t1\t1\t2\tCD24 molecule \t\t\t\t'


def opener(*lines, header=HEADER):
    return SimpleNamespace(result=io.StringIO('\n'.join([header, *lines, ''])))


def genes(*lines, header=HEADER):
    return list(iter_genes(opener(*lines, header=header)))


def extract(record, dataset='genes'):
    out = SilverExtractor('zhong2015', dataset)
    out.process_record(resource.datasets()[dataset].mapper(record), record, 'fixture:1', 1)
    return out


def test_parser_reads_genes_with_their_type_and_searches():
    assert genes(CDH13, ESAM) == [
        {
            'symbol': 'CDH13', 'entrez': '1012', 'cam_type': 'i', 'name': 'cadherin 13',
            'interpro': True, 'go': True, 'ncbi': True, 'sources': 3,
        },
        {
            'symbol': 'ESAM', 'entrez': '90952', 'cam_type': 'tj',
            'name': 'endothelial cell adhesion molecule',
            'interpro': False, 'go': False, 'ncbi': True, 'sources': 1,
        },
    ]


def test_the_placeholder_entrez_id_is_dropped():
    # 100000000 stands for three different genes.
    (record,) = genes(CD24)
    assert record['entrez'] is None
    out = extract(record)
    (relation,) = out.relations
    gene = out.entities[relation.subject_entity_key]
    assert {(i['ns'], i['id']) for i in gene.identifiers} == {
        (Namespace.GENESYMBOL, 'CD24'),
        (Namespace.NAME, 'CD24 molecule'),
    }


def test_gene_is_associated_with_its_type():
    (record,) = genes(ITGA6)
    out = extract(record)
    (relation,) = out.relations
    gene = out.entities[relation.subject_entity_key]
    term = out.entities[relation.object_entity_key]
    assert relation.predicate == 'associated_with'
    assert (gene.entity_type, gene.namespace, gene.identifier, gene.taxon) == (
        'gene', 'entrez', '3655', '9606'
    )
    assert {(i['ns'], i['id']) for i in gene.identifiers} == {
        (Namespace.ENTREZ, '3655'),
        (Namespace.GENESYMBOL, 'ITGA6'),
        (Namespace.NAME, 'integrin alpha 6'),
    }
    assert (term.entity_type, term.namespace, term.identifier) == (
        'ontology_class', TYPE_NAMESPACE, 'm'
    )
    attributes = {a['term']: a['value'] for a in relation.annotations if a['scope'] == 'relation'}
    # Flags are 'true'/'false'; the number of searches is a Measurement (its text
    # value is empty).
    assert attributes == {
        'zhong2015:found_by_interpro': 'true',
        'zhong2015:found_by_go': 'false',
        'zhong2015:found_by_ncbi': 'true',
        'zhong2015:n_searches': '',
    }
    (searches,) = [a for a in relation.annotations if a['term'] == 'zhong2015:n_searches']
    assert searches['quantity']['has_numeric_value'] == 2.0
    assert searches['quantity']['source_field'] == 'Sum'


def test_the_vocabulary_dataset_is_called_terms():
    assert 'terms' in resource.datasets()
    assert 'cam_types' not in resource.datasets()


def test_types_are_terms_with_names_and_the_statement_objects_are_those_terms():
    types = list(iter_cam_types(opener(CDH1, CDH13, ITGA6, ESAM)))
    assert [t['cam_type'] for t in types] == ['i', 'm', 'tj']
    for record, name in zip(types, ('iCAM', 'matrix adhesion', 'tight junction')):
        out = extract(record, 'terms')
        (term,) = out.entities.values()
        assert (term.entity_type, term.namespace) == ('ontology_class', TYPE_NAMESPACE)
        assert term.identifier == record['cam_type']
        assert (Namespace.NAME, name) in {(i['ns'], i['id']) for i in term.identifiers}
        assert not out.relations
    (record,) = genes(ITGA6)
    out = extract(record)
    (relation,) = out.relations
    obj = out.entities[relation.object_entity_key]
    assert (obj.entity_type, obj.namespace, obj.identifier) == (
        'ontology_class', TYPE_NAMESPACE, 'm'
    )


def test_all_eight_codes_have_a_name():
    assert set(CAM_TYPES) == {'i', 'm', 'ag', 'aj', 'c', 'fa', 'tj', 'my'}


def test_no_class_is_inferred_from_the_gene_name():
    # The table's one type is the only class; the cadherin is an iCAM here
    # because the table says so, not because its name says cadherin.
    (record,) = genes(CDH1)
    out = extract(record)
    assert len(out.relations) == 1


def test_unknown_type_flag_sum_and_header_raise():
    (unknown,) = genes(CDH1.replace('\ti\t', '\tzz\t'))
    with pytest.raises(ValueError, match='Unknown Zhong2015 CAM type'):
        resource.genes.mapper(unknown)
    with pytest.raises(ValueError, match='Inter Pro'):
        genes(CDH1.replace('\t1\t1\t1\t3', '\t2\t1\t1\t3'))
    with pytest.raises(ValueError, match='Sum'):
        genes(CDH1.replace('\t3\tcadherin', '\t2\tcadherin'))
    with pytest.raises(ValueError, match='Entrez'):
        genes(CDH1.replace('999', 'n/a'))
    with pytest.raises(ValueError, match='lacks expected columns'):
        genes(CDH1, header=HEADER.replace('CAM type', 'type'))
    with pytest.raises(ValueError, match='occurs twice'):
        genes(CDH1, CDH1)
    with pytest.raises(ValueError):
        genes()
