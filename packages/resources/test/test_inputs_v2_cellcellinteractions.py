"""Cell-Cell Interaction Database: parser, mappers and extractor on source rows."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import entity_type, slot_name
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import KnowledgeOrigin

from pypath.inputs_v2 import cellcellinteractions as cci
from pypath.inputs_v2.parsers import cellcellinteractions as parser

# Lines copied from protein_types.txt.
TABLE = (
    'hgnc_symbol\tcurated_protein_type\n'
    'SBSPON\tECM/Receptor/Ligand\n'
    'A2M\tLigand\n'
    'ITGB1\tReceptor\n'
    'COL1A1\tECM\n'
    'EFEMP2\tECM/Ligand\n'
)


def opener(text=TABLE):
    return SimpleNamespace(result=text)


def extract(dataset, rows):
    out = SilverExtractor('cellcellinteractions', dataset)
    mapper = cci.resource.datasets()[dataset].mapper
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def test_the_vocabulary_dataset_is_called_terms():
    assert 'terms' in cci.resource.datasets()
    assert 'roles' not in cci.resource.datasets()


def test_parser_splits_the_roles_of_a_cell():
    rows = list(parser.iter_genes(opener()))
    assert rows[0] == {
        'symbol': 'SBSPON',
        'roles': ['ECM', 'Receptor', 'Ligand'],
        'protein_type': 'ECM/Receptor/Ligand',
    }
    assert [r['roles'] for r in rows[1:]] == [
        ['Ligand'],
        ['Receptor'],
        ['ECM'],
        ['ECM', 'Ligand'],
    ]
    assert [r['role'] for r in parser.iter_roles(opener())] == [
        'ECM',
        'Receptor',
        'Ligand',
    ]


@pytest.mark.parametrize(
    'text',
    [
        '',
        'hgnc_symbol\tcurated_protein_type\n',
        'symbol\ttype\nA2M\tLigand\n',
        'hgnc_symbol\tcurated_protein_type\nA2M\tLigand\nA2M\tReceptor\n',
        'hgnc_symbol\tcurated_protein_type\nA2M\tLigand/Ligand\n',
        'hgnc_symbol\tcurated_protein_type\nA2M\t\n',
    ],
)
def test_parser_rejects_a_table_that_is_not_the_protein_types(text):
    with pytest.raises(ValueError):
        list(parser.iter_genes(opener(text)))


def test_role_terms_have_names_and_definitions():
    rows = list(parser.iter_roles(opener()))
    out = extract('terms', rows)
    terms = {e.identifier: e for e in out.entities.values()}
    assert set(terms) == {'ECM', 'Receptor', 'Ligand'}
    assert {e.entity_type for e in terms.values()} == {'ontology_class'}
    assert {e.namespace for e in terms.values()} == {'cellcellinteractions_role'}
    record = cci.map_role({'role': 'ECM'})
    assert record.identifiers[1] == (Namespace.NAME, 'Extracellular matrix')
    assert not out.relations


def test_gene_is_associated_with_the_term_of_each_role():
    rows = list(parser.iter_genes(opener()))
    record = cci.map_gene(rows[0])
    assert entity_type(record.type) == 'gene'
    assert record.identifiers == [(Namespace.GENESYMBOL, 'SBSPON')]

    out = extract('genes', rows)
    assert len(out.relations) == 3 + 1 + 1 + 1 + 2
    assert {r.predicate for r in out.relations} == {'associated_with'}
    genes = [e for e in out.entities.values() if e.entity_type == 'gene']
    assert {e.taxon for e in genes} == {'9606'}


def test_objects_are_the_same_entities_as_the_terms():
    terms = extract('terms', list(parser.iter_roles(opener())))
    genes = extract('genes', list(parser.iter_genes(opener())))
    term_identity = {
        (e.entity_type, e.namespace, e.identifier) for e in terms.entities.values()
    }
    object_identity = {
        (e.entity_type, e.namespace, e.identifier)
        for e in (genes.entities[r.object_entity_key] for r in genes.relations)
    }
    assert term_identity == object_identity == {
        ('ontology_class', 'cellcellinteractions_role', role)
        for role in ('ECM', 'Receptor', 'Ligand')
    }


def test_nothing_is_derived_from_the_combined_labels():
    genes = extract('genes', list(parser.iter_genes(opener())))
    # "ECM/Ligand" is two roles, not a third term, and no gene gets a role
    # it does not list.
    objects = {genes.entities[r.object_entity_key].identifier for r in genes.relations}
    assert objects == {'ECM', 'Receptor', 'Ligand'}
    a2m = [
        r
        for r in genes.relations
        if genes.entities[r.subject_entity_key].identifier == 'A2M'
    ]
    assert [genes.entities[r.object_entity_key].identifier for r in a2m] == ['Ligand']


def test_unknown_role_raises():
    with pytest.raises(ValueError, match='Unknown'):
        cci.map_gene({'symbol': 'X', 'roles': ['Adhesion'], 'protein_type': 'Adhesion'})
    with pytest.raises(ValueError, match='Unknown'):
        cci.map_role({'role': 'Adhesion'})


def test_knowledge_origin_is_computationally_predicted():
    assert tuple(cci.config.knowledge_origin) == (KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,)
