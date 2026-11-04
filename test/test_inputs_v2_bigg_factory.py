"""The shared BiGG factory, proven on synthetic models rather than any real one."""

from pypath.inputs_v2 import _bigg
from pypath.inputs_v2.parsers import _bigg as bigg_parser
from pypath.internals.silver_schema import format_term

_SYNTHETIC_MODEL = {
    'compartments': {'c': 'cytosol', 'e': 'extracellular space'},
    'metabolites': [
        {'id': 'atp_c', 'name': 'ATP', 'formula': 'C10H16N5O13P3', 'charge': -4, 'annotation': {}},
        {'id': 'h2o_c', 'name': 'Water', 'formula': 'H2O', 'charge': 0, 'annotation': {}},
    ],
    'reactions': [
        {
            'id': 'RX1',
            'name': 'synthetic reaction',
            'subsystem': 'Test metabolism',
            'lower_bound': -1000.0,
            'upper_bound': 1000.0,
            'gene_reaction_rule': '10 and 20',
            'metabolites': {'atp_c': -1, 'h2o_c': 1},
            'annotation': {'sbo': ['SBO:0000176']},
        },
    ],
    'genes': [
        {'id': '10', 'name': 'geneA'},
        {'id': '20', 'name': 'geneB'},
    ],
}


def test_orphan_filter_controls_whether_reaction_metabolites_are_listed_standalone():
    with_filter = list(bigg_parser._parse_metabolites(_SYNTHETIC_MODEL, orphan_filter=True))
    without_filter = list(bigg_parser._parse_metabolites(_SYNTHETIC_MODEL, orphan_filter=False))
    # Both synthetic metabolites participate in the one reaction, so the
    # orphan-only pass yields nothing; the unfiltered pass yields both.
    assert with_filter == []
    assert {m['bigg_metabolite_id'] for m in without_filter} == {'atp', 'h2o'}


def test_non_default_taxon_reaches_genes_and_gpr_catalysts():
    genes = _bigg.genes_schema(taxon='NCBITaxon:7227')
    gene_entity = genes({'entrez_id': '10', 'name': 'geneA'})
    assert [a.value for a in gene_entity.annotations if a.term == 'in_taxon'] == ['NCBITaxon:7227']

    reactions = _bigg.reactions_schema(taxon='NCBITaxon:7227', gpr_prefix='synthetic.gpr_and')
    row = {
        'bigg_reaction_id': 'RX1', 'reactants': 'atp:c:1', 'products': 'h2o:e:1',
        'gene_rule_clauses': [['10', '20']],
    }
    entity = reactions(row)
    gpr_group = next(
        m.member for m in entity.membership
        if format_term(m.predicate) == 'associated_with'
    )
    catalyst_taxa = {
        a.value
        for member in gpr_group.membership
        for a in member.member.annotations
        if format_term(a.term) == 'in_taxon'
    }
    assert catalyst_taxa == {'NCBITaxon:7227'}


def test_two_models_keep_gpr_identifiers_apart_for_the_same_gene_rule_shape():
    row = {'complex_subunits': '10||20'}
    a = _bigg.enzyme_complexes_schema(taxon='NCBITaxon:9606', gpr_prefix='model_a.gpr_and')(row)
    b = _bigg.enzyme_complexes_schema(taxon='NCBITaxon:9606', gpr_prefix='model_b.gpr_and')(row)
    assert a.identifiers != b.identifiers
