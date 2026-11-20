"""A lipid name Goslin cannot parse is marked on the occurrence that fell back to its own identity."""

from contextlib import closing

import pytest

from library_fixture import build_fixture_library
from omnipath_build.writer import _resolution_annotations
from omnipath_resolver import observations
from omnipath_resolver.canonical.match import LibraryMatcher
from test_canonical import obs


@pytest.fixture
def matcher(tmp_path, monkeypatch):
    monkeypatch.setenv("OMNIPATH_GOSLIN_CACHE", str(tmp_path / "goslin"))
    monkeypatch.setattr(observations, "_LIPID_CACHE", None)
    observations.lipid_name_result.cache_clear()
    with closing(LibraryMatcher(build_fixture_library(tmp_path))) as matcher:
        yield matcher
    observations.lipid_name_result.cache_clear()


def status_of(matcher, *observations_):
    keyed = {f"o{i}": item for i, item in enumerate(observations_)}
    result = matcher.match(keyed)
    return [getattr(result[key], "lipid_name_status", None) for key in keyed]


def lipid(namespace, identifier, *names):
    return obs("chemical_entity", namespace, identifier, [("name", name) for name in names])


def test_an_unparsable_name_from_a_lipid_resource_is_marked(matcher):
    [status] = status_of(matcher, lipid("lipidmaps", "LM0000000", "not a lipid name at all"))
    assert status == "unparsed"


def test_a_swisslipids_identifier_makes_the_name_a_candidate(matcher):
    candidate = obs(
        "chemical_entity",
        "hmdb",
        "HMDB9999999",
        [("swisslipids", "SLM:000000000"), ("name", "garbled lipid")],
    )
    [status] = status_of(matcher, candidate)
    assert status == "unparsed"


def test_a_parsable_name_is_not_marked(matcher):
    [status] = status_of(matcher, lipid("lipidmaps", "LM0000001", "PC 16:0/18:1"))
    assert status is None


def test_a_name_from_another_resource_is_not_a_candidate(matcher):
    [status] = status_of(matcher, obs("chemical_entity", "kegg", "C99999", [("name", "Mystery")]))
    assert status is None


def test_a_lipid_without_a_name_is_not_marked(matcher):
    [status] = status_of(matcher, obs("chemical_entity", "lipidmaps", "LM0000002"))
    assert status is None


def test_the_marker_reaches_the_evidence_annotations(matcher):
    target = matcher.match({"a": lipid("lipidmaps", "LM0000000", "not a lipid name")})["a"]
    annotations = _resolution_annotations([target])
    assert {(a["term"], a["value"]) for a in annotations} == {
        ("omnipath:lipid_name_status", "unparsed")
    }


def test_a_lipid_node_without_a_structure_is_marked(matcher):
    from omnipath_resolver.canonical.match import Node
    from omnipath_resolver.canonical.policy import get_policy

    observed = {"swisslipids": ["SLM:000000999"], "name": ["garbled lipid"]}
    node = Node("n1", "9606", "swisslipids", "SLM:000000999", None, ())
    policy = get_policy("chemical_entity")
    marked = matcher._build(None, policy, observed, node, "node")
    assert marked.lipid_name_status == "unparsed"
    parsed = matcher._build(None, policy, observed | {"name": ["PC 16:0/18:1"]}, node, "node")
    assert parsed.lipid_name_status is None
    structure = Node("n2", "9606", "inchikey", "AAAAAAAAAAAAAA-BBBBBBBBBB-N", None, ())
    assert matcher._build(None, policy, observed, structure, "node").lipid_name_status is None
