"""A KEGG reference map belongs to no organism."""

from omnipath_resources.inputs_v2.parsers import kegg

REACTION = """\
ENTRY       R00001                      Reaction
NAME        a reaction
DEFINITION  A <=> B
EQUATION    C00001 <=> C00002
PATHWAY     rn00010  Glycolysis / Gluconeogenesis
            rn01100  Metabolic pathways
///
"""


def pathway_rows(organism):
    records = kegg._build_records(
        {'get_reaction': REACTION.splitlines(keepends=True)}, organism
    )
    return {row['pathway_id']: row for row in records['pathways']}


def test_reference_maps_found_through_reactions_carry_no_taxon():
    rows = pathway_rows('mmu')
    assert set(rows) == {'rn00010', 'rn01100'}
    assert {row['taxon_id'] for row in rows.values()} == {''}


def test_the_organism_chosen_for_a_build_does_not_reach_reference_maps():
    assert pathway_rows('hsa') == pathway_rows('mmu')
