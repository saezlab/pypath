"""Metabolic Atlas: Human-GEM and the animal GEMs in one resource."""

from omnipath_build.silver import SilverExtractor
from omnipath_resources.inputs_v2 import metatlas
from omnipath_resources.inputs_v2.parsers import metatlas as parsers

MOUSE = next(model for model in metatlas.MODELS if model.name == 'Mouse-GEM')
RAT = next(model for model in metatlas.MODELS if model.name == 'Rat-GEM')


def test_every_model_has_its_datasets_and_human_gem_keeps_its_names():
    names = {name for name in dir(metatlas.resource) if not name.startswith('_')}
    assert {'metabolites', 'reactions', 'enzyme_complexes'} <= names
    for prefix in ('mouse', 'rat', 'zebrafish', 'fruitfly', 'worm'):
        assert {f'{prefix}_metabolites', f'{prefix}_reactions', f'{prefix}_enzyme_complexes'} <= names
    assert metatlas.HUMAN.gene_rule_namespace == 'human_gem.gpr_and'
    assert metatlas.download.url.endswith('/Human-GEM/main/model/Human-GEM.yml')


def _genes(group):
    return {
        (str(m.member.identifiers[0].type), m.member.identifiers[0].value, m.member.annotations[0].value)
        for m in group.membership
    }


def test_a_model_s_genes_carry_its_identifiers_and_taxon():
    group = metatlas._enzyme_complexes_schema({'complex_subunits': 'Acaa1a||Myo5b'}, MOUSE)
    assert str(group.identifiers[0].type) == 'mouse_gem.gpr_and'
    assert _genes(group) == {
        ('genesymbol', 'Acaa1a', 'NCBITaxon:10090'),
        ('genesymbol', 'Myo5b', 'NCBITaxon:10090'),
    }
    human = metatlas.enzyme_complexes_schema({'complex_subunits': 'ENSG1||ENSG2'})
    assert str(human.identifiers[0].type) == 'human_gem.gpr_and'
    assert {gene[0] for gene in _genes(human)} == {'ensg'}


def test_the_same_gene_rule_of_two_models_is_two_groups():
    row = {'complex_subunits': 'Acaa1a||Myo5b'}
    mouse = metatlas._enzyme_complexes_schema(row, MOUSE)
    rat = metatlas._enzyme_complexes_schema(row, RAT)
    assert (str(mouse.identifiers[0].type), mouse.identifiers[0].value) != (
        str(rat.identifiers[0].type),
        rat.identifiers[0].value,
    )


def test_a_model_reads_its_own_cross_references(monkeypatch):
    seen = []

    def metabolite_xrefs(model):
        seen.append(('metabolites', model))
        return {'MAM1': {'chebi': '15422'}}

    def reaction_xrefs(model):
        seen.append(('reactions', model))
        return {'MAR1': {'rhea': '13065'}}

    monkeypatch.setattr(parsers, '_other_metabolite_xrefs', metabolite_xrefs)
    monkeypatch.setattr(parsers, '_other_reaction_xrefs', reaction_xrefs)
    data = {
        'metabolites': [{'id': 'MAM1c', 'compartment': 'c'}],
        'reactions': [{
            'id': 'MAR1',
            'metabolites': {'MAM1c': -1},
            'lower_bound': 0,
            'upper_bound': 1000,
            'gene_reaction_rule': 'Acaa1a or Acaa1b',
            'eccodes': '2.3.1.16',
        }],
    }
    (reaction,) = parsers._reactions(data, 'Mouse-GEM')
    assert reaction['rhea'] == '13065'
    assert reaction['reactant_chebi'] == '15422'
    assert reaction['gene_rule_clauses'] == [['Acaa1a'], ['Acaa1b']]
    (metabolite,) = parsers._metabolites(data, 'Mouse-GEM')
    assert metabolite['chebi'] == '15422'
    assert set(seen) == {('metabolites', 'Mouse-GEM'), ('reactions', 'Mouse-GEM')}


def test_mouse_reaction_gene_rules_are_mouse_genes(monkeypatch):
    monkeypatch.setattr(parsers, '_other_metabolite_xrefs', lambda model: {})
    row = {
        'human_gem_reaction_id': 'MAR1',
        'name': 'reaction',
        'reactants': '',
        'products': '',
        'gene_rule_clauses': [['Acaa1a']],
    }
    schema = metatlas._reactions_schema(MOUSE)
    ex = SilverExtractor(metatlas.config.id.value, 'mouse_reactions')
    ex.process_record(metatlas._map_with_structures(schema, row, model=MOUSE), row, 'mouse_reactions:1', 1)
    genes = {
        (str(i['ns']), i['id'])
        for entity in ex.entities.values()
        for i in entity.identifiers
        if str(i['ns']) == 'genesymbol'
    }
    assert genes == {('genesymbol', 'Acaa1a')}
