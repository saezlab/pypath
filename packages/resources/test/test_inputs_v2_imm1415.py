"""iMM1415: the mouse BiGG model's reactions, gene rules and genes, read as Recon3D's are."""

from omnipath_build.silver import SilverExtractor
from omnipath_resources.inputs_v2 import imm1415, recon3d


def _extract(mapper, row, dataset):
    ex = SilverExtractor(imm1415.config.id.value, dataset)
    ex.process_record(mapper(row), row, f'{dataset}:1', 1)
    return ex


def test_reactions_genes_and_groups_are_mouse():
    assert {'metabolites', 'reactions', 'enzyme_complexes', 'genes'} <= {
        name for name in dir(imm1415.resource) if not name.startswith('_')
    }
    gene = imm1415.genes_schema({'entrez_id': '12846', 'name': 'Comt'})
    assert {(a.term, a.value) for a in gene.annotations} >= {('in_taxon', 'NCBITaxon:10090')}
    group = imm1415.enzyme_complexes_schema({'complex_subunits': '11946||11947'})
    assert str(group.identifiers[0].type) == 'imm1415.gpr_and'
    assert {m.member.annotations[0].value for m in group.membership} == {'NCBITaxon:10090'}


def test_a_gene_rule_of_recon3d_and_imm1415_is_two_groups():
    row = {'complex_subunits': '11946||11947'}
    mouse, human = imm1415.enzyme_complexes_schema(row), recon3d.enzyme_complexes_schema(row)
    assert str(mouse.identifiers[0].type) != str(human.identifiers[0].type)


def test_reaction_links_mouse_genes():
    row = {
        'bigg_reaction_id': 'HEX1',
        'name': 'hexokinase',
        'reactants': '',
        'products': '',
        'gene_rule_clauses': [['15275']],
    }
    ex = _extract(imm1415.reactions_schema, row, 'reactions')
    genes = [entity for entity in ex.entities.values() if str(entity.entity_type) == 'gene']
    assert [(i['id'], {a['value'] for a in e.annotations or []}) for e in genes for i in e.identifiers] == [
        ('15275', {'NCBITaxon:10090'})
    ]
