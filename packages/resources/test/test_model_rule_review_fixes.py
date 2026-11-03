"""Regression coverage for source model identity and Boolean gene rules."""

import pytest

from omnipath_build.silver import SilverExtractor
from pypath.internals.silver_schema import format_term

from pypath.inputs_v2 import metatlas, recon3d
from pypath.inputs_v2.parsers import metatlas as human_parser
from pypath.inputs_v2.parsers import recon3d as recon_parser


@pytest.mark.parametrize('rule,expected', [
    ('1 and 2', [['1', '2']]),
    ('1 or 2', [['1'], ['2']]),
    ('1 and (2 or 3)', [['1', '2'], ['1', '3']]),
    ('(1 or 2) and (3 or 4)', [['1', '3'], ['1', '4'], ['2', '3'], ['2', '4']]),
    ('2_AT1 and 1_AT2', [['1', '2']]),
    ('', []),
])
def test_boolean_rules(rule, expected):
    assert recon_parser._parse_gene_rule(rule) == expected


@pytest.mark.parametrize('module', [metatlas, recon3d])
def test_named_logical_groups_and_stable_composition(module):
    group = module.enzyme_complexes_schema({'complex_subunits': '2||1'})
    reordered = module.enzyme_complexes_schema({'complex_subunits': '1||2||1'})
    assert group.identifiers == reordered.identifiers
    assert format_term(group.type) == 'named_thing'
    assert all(format_term(m.member.type) == 'gene' for m in group.membership)
    assert all(format_term(m.predicate) == 'has_member' for m in group.membership)
    assert module.enzyme_complexes_schema({'complex_subunits': ''}) is None


@pytest.mark.parametrize('module,id_field', [(metatlas, 'human_gem_reaction_id'), (recon3d, 'bigg_reaction_id')])
def test_reaction_views_and_gene_clause_extraction(module, id_field):
    row = {id_field: 'reaction1', 'gene_rule_clauses': [['1', '2'], ['3']], 'reactants': 'm:c:1', 'products': 'm:e:1'}
    outputs = []
    for schema in [module.reactions_schema, module.transport_reactions_schema]:
        entity = schema(row)
        associations = [m for m in entity.membership if format_term(m.predicate) == 'associated_with']
        assert len(associations) == 2
        assert format_term(associations[0].member.type) == 'named_thing'
        assert len(associations[0].member.membership) == 2
        assert format_term(associations[1].member.type) == 'gene'
        output = SilverExtractor(module.config.id.value, 'fixture')
        output.process_record(entity, row, 'fixture:1', 1)
        outputs.append(output)
    assert outputs[0].entities.keys() == outputs[1].entities.keys()
    assert len(outputs[0].relations) == 6  # 2 chemical, 2 alternatives, 2 required genes


def test_parsers_retain_rule_bounds_and_compartments(monkeypatch):
    monkeypatch.setattr(human_parser, '_metabolite_xrefs', lambda: {})
    data = {'compartments': {'c': 'cytosol', 'e': 'external'}, 'reactions': [{
        'id': 'r', 'gene_reaction_rule': '1 and (2 or 3)',
        'lower_bound': -50, 'upper_bound': 100, 'metabolites': {},
    }]}
    for parser in [human_parser._reactions, recon_parser._parse_reactions]:
        row = next(parser(data))
        assert row['gene_reaction_rule'] == '1 and (2 or 3)'
        assert row['gene_rule_clauses'] == [['1', '2'], ['1', '3']]
        assert row['lower_bound'] == -50 and row['upper_bound'] == 100
        assert row['compartments'] == data['compartments']


def test_logical_group_identity_is_model_scoped():
    row = {'complex_subunits': '1||2'}
    assert metatlas.enzyme_complexes_schema(row).identifiers != recon3d.enzyme_complexes_schema(row).identifiers


@pytest.mark.parametrize('rule', ['1 and', '(1 or 2', '1 2', '1 or )'])
def test_malformed_gene_rules_do_not_silently_flatten(rule):
    with pytest.raises(ValueError):
        recon_parser._parse_gene_rule(rule)


def test_human_gem_uses_explicit_compartment_context(monkeypatch):
    monkeypatch.setattr(human_parser, '_metabolite_xrefs', lambda: {})
    data = {
        'metabolites': [{'id': 'MAM1cyt', 'compartment': 'cyt'}],
        'reactions': [{'id': 'r', 'metabolites': {'MAM1cyt': -1}, 'lower_bound': 0, 'upper_bound': 2}],
    }
    assert next(human_parser._reactions(data))['reactants'] == 'MAM1:cyt:1'
