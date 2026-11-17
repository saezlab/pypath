"""What the resolver does with the alanines and the EGFR orthologs.

These pin the outcome on a small identity library, so a change to the voting or the
matching that merges or splits them shows up here.
"""

from contextlib import closing

import pytest

from library_fixture import (
    BETA_ALANINE,
    D_ALANINE,
    DL_ALANINE,
    L_ALANINE,
    build_fixture_library,
)
from omnipath_resolver.canonical.match import LibraryMatcher
from test_canonical import obs


@pytest.fixture
def matcher(tmp_path):
    with closing(LibraryMatcher(build_fixture_library(tmp_path, "benchmarks"))) as matcher:
        yield matcher


def match(matcher, *observations):
    keyed = {f"o{i}": item for i, item in enumerate(observations)}
    result = matcher.match(keyed)
    return [result[key] for key in keyed]


def chemical(namespace, identifier, *extra):
    return obs("chemical_entity", namespace, identifier, list(extra))


def test_stereoisomers_of_alanine_stay_distinct(matcher):
    l, d, dl = match(
        matcher,
        chemical("chebi", "CHEBI:16977"),
        chemical("chebi", "CHEBI:15570"),
        chemical("chebi", "CHEBI:16449"),
    )
    keys = [(m.canonical_namespace, m.canonical_identifier) for m in (l, d, dl)]
    assert keys == [("inchikey", L_ALANINE), ("inchikey", D_ALANINE), ("inchikey", DL_ALANINE)]
    assert len(set(keys)) == 3
    assert [m.label.lower() for m in (l, d, dl)] == ["l-alanine", "d-alanine", "alanine"]


def test_the_alanines_share_a_connectivity_block_and_beta_alanine_does_not(matcher):
    found = match(
        matcher,
        chemical("inchikey", L_ALANINE),
        chemical("inchikey", D_ALANINE),
        chemical("inchikey", DL_ALANINE),
        chemical("inchikey", BETA_ALANINE),
    )
    blocks = [m.canonical_identifier.split("-")[0] for m in found]
    assert blocks[0] == blocks[1] == blocks[2]
    assert blocks[3] != blocks[0]
    assert len({m.canonical_identifier for m in found}) == 4


def test_the_residue_is_not_merged_into_an_alanine(matcher):
    [residue] = match(matcher, chemical("chebi", "CHEBI:90000"))
    assert residue.canonical_namespace != "inchikey"
    assert residue.canonical_identifier not in {L_ALANINE, D_ALANINE, DL_ALANINE, BETA_ALANINE}


def test_a_bare_name_does_not_pick_one_alanine(matcher):
    [named] = match(matcher, chemical("kegg", "C99999", ("name", "alanine")))
    assert named.canonical_identifier not in {L_ALANINE, D_ALANINE, BETA_ALANINE}


@pytest.mark.parametrize(
    ("namespace", "identifier", "taxon"),
    [("uniprot", "P00533", "9606"), ("entrez", "1956", "9606"), ("genesymbol", "EGFR", "9606")],
)
def test_human_egfr_is_one_gene_by_every_identifier(matcher, namespace, identifier, taxon):
    [gene] = match(matcher, obs("protein", namespace, identifier, taxon=taxon))
    assert (gene.canonical_namespace, gene.canonical_identifier) == ("entrez", "1956")
    assert gene.label == "EGFR"
    assert gene.taxon == "9606"


def test_each_ortholog_references_its_own_gene(matcher):
    human, mouse, rat = match(
        matcher,
        obs("protein", "genesymbol", "EGFR", taxon="9606"),
        obs("protein", "uniprot", "Q01279", taxon="10090"),
        obs("protein", "entrez", "24329", taxon="10116"),
    )
    assert [m.canonical_identifier for m in (human, mouse, rat)] == ["1956", "13649", "24329"]
    assert {m.taxon for m in (human, mouse, rat)} == {"9606", "10090", "10116"}


def test_a_symbol_is_not_carried_across_taxa(matcher):
    [mouse] = match(matcher, obs("protein", "genesymbol", "Egfr", taxon="10090"))
    assert mouse.canonical_identifier == "13649"
    [rat] = match(matcher, obs("protein", "genesymbol", "Egfr", taxon="10116"))
    assert rat.canonical_identifier == "24329"
