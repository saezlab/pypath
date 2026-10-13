"""Source contract of the PROGENy model input."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor

from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import progeny
from pypath.inputs_v2.parsers import progeny as parser
from pypath.inputs_v2.progeny import DESCRIPTIONS, PATHWAY, resource

# Rows of model_human_full in model_human_full.rda, as iter_rdata yields them.
RFC2_EGFR = {'gene': 'RFC2', 'pathway': 'EGFR', 'weight': 1.4706466171905381, 'p.value': 0.0016552751715941223}
CBX6_EGFR = {'gene': 'CBX6', 'pathway': 'EGFR', 'weight': -0.6755073377638218, 'p.value': 0.017641119340513132}
TNFAIP3_NFKB = {'gene': 'TNFAIP3', 'pathway': 'NFkB', 'weight': 28.373752313465076, 'p.value': 1.867140427866091e-71}
RFC2_WNT = {'gene': 'RFC2', 'pathway': 'WNT', 'weight': 0.465622632141, 'p.value': 0.704579815833}
RFC2_MAPK = {'gene': 'RFC2', 'pathway': 'MAPK', 'weight': 2.277341180050157, 'p.value': 3.324844332731664e-11}
RFC2_TRAIL = {'gene': 'RFC2', 'pathway': 'Trail', 'weight': -0.8006773132589663, 'p.value': 0.6240861307828072}
SOURCE_ROWS = [RFC2_EGFR, CBX6_EGFR, TNFAIP3_NFKB, RFC2_WNT, RFC2_MAPK, RFC2_TRAIL]


@pytest.fixture
def model(monkeypatch):
    def read(rows):
        monkeypatch.setattr(parser, '_read', lambda opener: iter(rows))

    return read


def parsed(rows, model):
    model(rows)
    return list(parser.iter_model(SimpleNamespace(result=None)))


def extract(row, dataset='responses'):
    record = resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor('progeny', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return out


def relation_of(row):
    out = extract(row)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def quantities(relation):
    return {
        a['term']: a['quantity']['has_numeric_value']
        for a in relation.annotations
        if a['quantity']
    }


def test_parser_renames_the_p_value_column(model):
    rows = parsed(SOURCE_ROWS, model)
    assert rows[0] == {'gene': 'RFC2', 'pathway': 'EGFR', 'weight': 1.4706466171905381,
                       'p_value': 0.0016552751715941223}
    assert len(rows) == len(SOURCE_ROWS)


def test_parser_rejects_missing_values_and_empty_tables(model):
    with pytest.raises(ValueError, match='lacks'):
        parsed([dict(RFC2_EGFR, **{'p.value': None})], model)
    with pytest.raises(ValueError, match='lacks'):
        parsed([{'gene': 'RFC2', 'pathway': 'EGFR', 'weight': 1.0}], model)
    with pytest.raises(ValueError, match='no rows'):
        parsed([], model)


def test_pathways_are_listed_once(model):
    model(SOURCE_ROWS)
    rows = list(parser.iter_pathways(SimpleNamespace(result=None)))
    assert [r['pathway'] for r in rows] == ['EGFR', 'NFkB', 'WNT', 'MAPK', 'Trail']


def test_gene_is_associated_with_the_pathway_with_weight_and_p_value(model):
    row = parsed([RFC2_EGFR], model)[0]
    relation, gene, pathway = relation_of(row)
    assert relation.predicate == 'associated_with'
    assert (gene.entity_type, gene.namespace, gene.identifier, gene.taxon) == (
        'gene', 'genesymbol', 'RFC2', '9606',
    )
    assert (pathway.entity_type, pathway.namespace, pathway.identifier) == (
        'ontology_class', PATHWAY, 'EGFR',
    )
    assert quantities(relation) == {
        'progeny:weight': 1.4706466171905381,
        'progeny:p_value': 0.0016552751715941223,
    }


def test_negative_weight_and_tiny_p_value_are_kept(model):
    rows = parsed([CBX6_EGFR, TNFAIP3_NFKB], model)
    assert quantities(relation_of(rows[0])[0])['progeny:weight'] == -0.6755073377638218
    assert quantities(relation_of(rows[1])[0]) == {
        'progeny:weight': 28.373752313465076,
        'progeny:p_value': 1.867140427866091e-71,
    }


def test_nothing_is_filtered_by_significance_or_rank(model):
    # A p-value of 0.70 and a weight near zero are statements like any other.
    (row,) = parsed([RFC2_WNT], model)
    assert quantities(relation_of(row)[0])['progeny:p_value'] == 0.704579815833
    mapped = [resource.responses.mapper(r) for r in parsed(SOURCE_ROWS, model)]
    assert all(record is not None for record in mapped)


def test_a_pathway_of_the_model_is_one_entity_with_a_description():
    out = extract({'pathway': 'EGFR'}, 'terms')
    (term,) = (e for e in out.entities.values() if e.label)
    assert (term.entity_type, term.namespace, term.identifier, term.label) == (
        'ontology_class', PATHWAY, 'EGFR', 'EGFR',
    )
    assert out.relations == []
    assert any(
        a['term'] == 'description' and a['value'].startswith('regulates growth, survival')
        for a in term.annotations
    )


def test_statement_object_is_the_entity_of_the_terms_dataset(model):
    _, _, pathway = relation_of(parsed([RFC2_MAPK], model)[0])
    out = extract({'pathway': 'MAPK'}, 'terms')
    (term,) = (e for e in out.entities.values() if e.label)
    assert (pathway.entity_type, pathway.namespace, pathway.identifier) == (
        term.entity_type, term.namespace, term.identifier,
    )


def test_every_pathway_of_the_model_has_a_description():
    assert len(DESCRIPTIONS) == 14
    assert set(DESCRIPTIONS) == {
        'Androgen', 'EGFR', 'Estrogen', 'Hypoxia', 'JAK-STAT', 'MAPK', 'NFkB',
        'p53', 'PI3K', 'TGFb', 'TNFa', 'Trail', 'VEGF', 'WNT',
    }


def test_unknown_pathway_raises(model):
    with pytest.raises(ValueError, match='Unknown PROGENy pathway'):
        progeny.map_response(dict(parsed([RFC2_EGFR], model)[0], pathway='Notch'))
    with pytest.raises(ValueError, match='Unknown PROGENy pathway'):
        progeny.map_pathway({'pathway': 'Notch'})


def test_nothing_else_is_asserted(model):
    relation, _, _ = relation_of(parsed([RFC2_TRAIL], model)[0])
    assert {a['term'] for a in relation.annotations} == {
        'progeny:weight', 'progeny:p_value', 'in_taxon',
    }
    # No direction or effect is stated between the gene and the pathway.
    assert relation.predicate not in ('affects', 'regulates')


def test_resource_declares_its_datasets_and_license():
    assert progeny.config.id == ResourceCv.PROGENY
    assert progeny.config.license == LicenseCV.APACHE_2_0
    assert set(resource.datasets()) == {'terms', 'responses'}
