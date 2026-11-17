"""``/sets/metsigdb``: the contract of cycles 012 and 013 over the subset of a served release."""

import pytest
from fastapi.testclient import TestClient

from metsigdb_served import SETS, key, serve
from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.server import create_app

DEFAULT = [
    "entity",
    "metabolite_label",
    "inchikey",
    "hmdb",
    "pubchem",
    "chebi",
    "kegg",
    "resource",
    "set",
    "set_label",
    "set_type",
    "set_sub_type",
    "set_size",
]
ALL = [
    "entity",
    "metabolite_label",
    "metabolite_entity_type",
    "metabolite_structure_key",
    "inchikey",
    "smiles",
    "hmdb",
    "pubchem",
    "chebi",
    "kegg",
    "has_structure",
    "exclusion_classes",
    "resource",
    "set",
    "set_entity",
    "set_label",
    "set_type",
    "set_sub_type",
    "organism",
    "set_size",
    "set_context",
    "provenance_source",
    "provenance_record",
    "build_id",
]
TOTAL_ROWS = sum(len(s[6]) for s in SETS)


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    root = serve(tmp_path_factory.mktemp("data"))
    engine = ReleaseQueries(root)
    with TestClient(create_app(queries=engine)) as client:
        yield client
    engine.close()


def get(client, path="/sets/metsigdb", **params):
    return client.get(path, params=params)


def rows(client, **params):
    response = get(client, **params)
    assert response.status_code == 200, response.text
    return response.json()["rows"]


def members(client, **params):
    return {(r["resource"], r["set"], r["entity"]) for r in rows(client, limit=1000, **params)}


# ----------------------------------------------------------------- shape and paging


def test_a_page_carries_the_default_fields_and_no_total_unless_asked(client):
    body = get(client).json()
    assert body["count"] == len(body["rows"]) == TOTAL_ROWS
    assert body["has_more"] is False
    assert body["total"] is None
    assert list(body["rows"][0]) == DEFAULT


def test_total_is_the_size_of_the_whole_result(client):
    body = get(client, limit=2, total="true").json()
    assert (body["count"], body["has_more"], body["total"]) == (2, True, TOTAL_ROWS)


def test_paging_neither_repeats_nor_skips_a_row(client):
    seen = []
    for offset in range(0, TOTAL_ROWS + 3, 3):
        seen += [
            (r["resource"], r["set"], r["entity"]) for r in rows(client, limit=3, offset=offset)
        ]
    assert seen == sorted(seen)
    assert len(seen) == len(set(seen)) == TOTAL_ROWS


def test_rows_come_in_the_order_of_the_key(client):
    found = [(r["resource"], r["set"], r["entity"]) for r in rows(client)]
    assert found == sorted(found)


@pytest.mark.parametrize("params", [{"limit": 0}, {"limit": 100_001}, {"offset": -1}])
def test_a_page_outside_its_bounds_is_refused(client, params):
    assert get(client, **params).status_code in (400, 422)


def test_a_parameter_the_contract_does_not_publish_is_refused_by_name(client):
    response = get(client, hmdb="HMDB0000538")
    assert response.status_code == 400
    assert "hmdb" in response.text


def test_the_old_names_of_set_and_entity_are_refused(client):
    for name in ("set_source_id", "metabolite_entity_id"):
        assert get(client, **{name: "x"}).status_code == 400


# ----------------------------------------------------------------- vocabularies


def test_vocabularies_match_without_regard_to_case(client):
    expected = members(client, resource="WikiPathways")
    assert expected and members(client, resource="wikipathways") == expected
    assert members(client, set_type="PATHWAY") == members(client, set_type="pathway")
    assert members(client, set_sub_type="Overview_Map") == members(
        client, set_sub_type="overview_map"
    )


@pytest.mark.parametrize(
    "params", [{"resource": "smpdb"}, {"set_type": "tissue"}, {"set_sub_type": "nonsense"}]
)
def test_a_value_outside_a_vocabulary_is_refused_whatever_its_case(client, params):
    response = get(client, **params)
    assert response.status_code == 400
    assert next(iter(params.values())) in response.text


def test_several_values_are_given_repeated_or_separated_by_commas(client):
    both = members(client, resource="Reactome") | members(client, resource="KEGG")
    assert members(client, resource="Reactome,KEGG") == both
    assert {r["resource"] for r in rows(client, resource=["Reactome", "KEGG"])} == {
        "Reactome",
        "KEGG",
    }


def test_organism_matches_explicit_values_only(client):
    found = rows(client, organism=9606, fields="organism")
    assert found and {r["organism"] for r in found} == {9606}
    assert all(r["set"] != "rn00010" for r in found)


# ----------------------------------------------------------------- sets and entities


def test_set_selects_by_the_sets_own_id(client):
    found = rows(client, set="R-HSA-70171")
    assert {r["resource"] for r in found} == {"Reactome"} and len(found) == 3


def test_a_macdb_set_never_reads_as_an_entity_identifier(client):
    assert {r["set"] for r in rows(client, set="MACDB:15378")} == {"MACDB:15378"}
    # The bare number is a ChEBI or PubChem identifier, and no metabolite carries it.
    assert rows(client, entity="15378") == []


@pytest.mark.parametrize(
    "value",
    [
        "ZKHQWZAMYRWXGA-KQYNXXCUSA-J",  # InChIKey
        "HMDB0000538",
        "hmdb0000538",
        "CHEBI:30616",
        "chebi:30616",
        "C00002",
        "5957",  # PubChem, a bare integer
        "30616",  # ChEBI, a bare integer
        key("atp"),  # the release's entity key
    ],
)
def test_entity_recognises_the_identifier_a_metabolite_is_given_by(client, value):
    found = rows(client, entity=value, fields="entity")
    assert found and {r["entity"] for r in found} == {key("atp")}


def test_a_bare_integer_is_looked_for_as_chebi_and_as_pubchem(client):
    # 962 is water's PubChem id; 17234 is glucose's ChEBI id.
    assert {r["entity"] for r in rows(client, entity="962,17234", fields="entity")} == {
        key("water"),
        key("glucose"),
    }


def test_an_identifier_no_rule_claims_is_refused_rather_than_answered_empty(client):
    response = get(client, entity="not an identifier")
    assert response.status_code == 400
    assert "not an identifier" in response.text


def test_a_metabolite_is_found_by_the_structure_key_it_shares_with_its_variants(client):
    found = rows(client, metabolite_structure_key="ZKHQWZAMYRWXGA", fields="entity")
    assert found and {r["entity"] for r in found} == {key("atp")}


# ----------------------------------------------------------------- fields


def test_fields_add_to_the_default(client):
    row = rows(client, fields="smiles,organism", limit=1)[0]
    assert list(row) == [f for f in ALL if f in DEFAULT + ["smiles", "organism"]]


def test_fields_may_be_repeated(client):
    row = rows(client, fields=["smiles", "organism"], limit=1)[0]
    assert "smiles" in row and "organism" in row


def test_fields_all_returns_every_published_field(client):
    row = rows(client, fields="all", limit=1)[0]
    assert list(row) == ALL


def test_an_unknown_field_is_refused_rather_than_dropped(client):
    response = get(client, fields="smiles,hmdbb")
    assert response.status_code == 400
    assert "hmdbb" in response.text


def test_json_fields_arrive_as_json_and_missing_values_as_null(client):
    row = next(
        r for r in rows(client, resource="ClassyFire", fields="all") if r["entity"] == key("atp")
    )
    assert row["set_context"] == {"assignment": "direct", "depth": 0}
    assert row["provenance_record"] == {"source": "x", "dataset": "d", "row_id": "1"}
    assert row["exclusion_classes"] == ["cofactors", "cofactors_extended"]
    assert row["has_structure"] is True
    mystery = next(
        r
        for r in rows(client, resource="ClassyFire", fields="all")
        if r["entity"] == key("mystery")
    )
    assert mystery["inchikey"] is None and mystery["has_structure"] is False


# ----------------------------------------------------------------- exclusion


def excluded(client, value):
    return members(client, exclude_metabolites=value)


def test_nothing_is_excluded_by_default_or_by_none(client):
    everything = members(client)
    assert excluded(client, "none") == everything and len(everything) == TOTAL_ROWS


def test_a_class_leaves_out_the_metabolites_in_it(client):
    left = excluded(client, "cofactors")
    assert not {m for m in left if m[2] == key("atp")}
    assert len(left) < TOTAL_ROWS


def test_all_leaves_out_every_class_including_the_extended_cofactors(client):
    left = excluded(client, "all")
    gone = {key("atp"), key("water"), key("sodium")}
    assert not {m for m in left if m[2] in gone}
    assert {m[2] for m in left} == {key("glucose"), key("mystery")}


def test_classes_may_be_combined_and_are_matched_without_regard_to_case(client):
    assert excluded(client, "Ions,small_molecules") == excluded(client, "ions") & excluded(
        client, "small_molecules"
    )


def test_the_extended_cofactor_class_is_its_own_class(client):
    # ATP is in both; leaving out only the extended class still removes it.
    assert not {m for m in excluded(client, "cofactors_extended") if m[2] == key("atp")}


def test_total_and_has_more_describe_the_filtered_result(client):
    body = get(client, exclude_metabolites="all", limit=2, total="true").json()
    full = len(excluded(client, "all"))
    assert body["total"] == full and body["has_more"] is (full > 2)


@pytest.mark.parametrize("value", ["metabolites", "none,ions", "all,ions", "all,none"])
def test_an_exclusion_that_names_nothing_or_contradicts_itself_is_refused(client, value):
    assert get(client, exclude_metabolites=value).status_code == 400


# ----------------------------------------------------------------- grouping


def test_group_returns_sets_grouped_by_resource_with_their_members(client):
    body = get(client, group="true").json()
    assert [g["resource"] for g in body["groups"]] == sorted({s[0] for s in SETS})
    reactome = next(g for g in body["groups"] if g["resource"] == "Reactome")
    glycolysis = next(s for s in reactome["sets"] if s["set"] == "R-HSA-70171")
    assert glycolysis["set_size"] == glycolysis["returned"] == 3 == len(glycolysis["members"])
    assert glycolysis["set_label"] == "Glycolysis"
    # The set side is written once, not on each member.
    assert "set_label" not in glycolysis["members"][0] and "resource" not in glycolysis


def test_under_grouping_the_counts_are_in_sets(client):
    body = get(client, group="true", limit=2, total="true").json()
    assert body["count"] == 2 and body["has_more"] is True and body["total"] == len(SETS)


def test_a_set_is_never_cut_by_the_row_ceiling(client, monkeypatch):
    from omnipath_api.routers import metsigdb

    monkeypatch.setattr(metsigdb, "ROW_CEILING", 4)
    body = get(client, group="true", limit=100).json()
    sizes = [len(s["members"]) for g in body["groups"] for s in g["sets"]]
    assert body["has_more"] is True
    assert sum(sizes) <= 4 or len(sizes) == 1
    assert all(len(s["members"]) == s["returned"] for g in body["groups"] for s in g["sets"])


def test_a_set_larger_than_the_ceiling_is_returned_whole_and_alone(client, monkeypatch):
    from omnipath_api.routers import metsigdb

    monkeypatch.setattr(metsigdb, "ROW_CEILING", 1)
    body = get(client, group="true", limit=100).json()
    assert body["count"] == 1 and body["has_more"] is True
    assert body["groups"][0]["sets"][0]["returned"] >= 1


def test_a_set_that_loses_every_member_is_still_returned_with_none(client):
    body = get(client, group="true", set="R-HSA-1", exclude_metabolites="cofactors").json()
    (only,) = body["groups"][0]["sets"]
    assert only["returned"] == 0 and only["members"] == [] and only["set_size"] == 1


def test_grouped_member_fields_follow_the_requested_fields(client):
    body = get(client, group="true", set="R-HSA-70171", fields="smiles,organism").json()
    (one,) = body["groups"][0]["sets"]
    assert (
        "organism" in one and "smiles" in one["members"][0] and "organism" not in one["members"][0]
    )


# ----------------------------------------------------------------- the set listing


def test_the_set_listing_has_no_members(client):
    body = get(client, "/sets/metsigdb/sets", total="true").json()
    assert body["total"] == len(SETS) == body["count"]
    glycolysis = next(s for s in body["sets"] if s["set"] == "R-HSA-70171")
    assert glycolysis == {
        "resource": "Reactome",
        "set": "R-HSA-70171",
        "set_entity": key("R-HSA-70171"),
        "set_label": "Glycolysis",
        "set_type": "pathway",
        "set_sub_type": None,
        "organism": 9606,
        "set_size": 3,
    }


def test_the_set_listing_takes_the_set_filters_and_pages(client):
    pathways = get(client, "/sets/metsigdb/sets", set_type="pathway").json()["sets"]
    assert {s["set_type"] for s in pathways} == {"pathway"}
    page = get(client, "/sets/metsigdb/sets", limit=2, offset=1).json()
    assert page["count"] == 2 and page["has_more"] is True


def test_the_set_listing_refuses_a_member_level_filter(client):
    assert get(client, "/sets/metsigdb/sets", entity="HMDB0000538").status_code == 400


# ----------------------------------------------------------------- public prefix


def test_the_route_answers_under_the_public_prefix(client):
    assert get(client, "/api/sets/metsigdb", limit=1).status_code == 200


# ----------------------------------------------------------------- the published description


def test_the_description_holds_the_coverage_measured_when_the_subset_was_built(client):
    schema = client.get("/openapi.json").json()
    description = schema["paths"]["/sets/metsigdb"]["get"]["description"]
    assert "| Reactome | 1 | 58% |" in description  # the share the fixture's build measured
    # What the cycle 012 text said about a build that no longer exists is gone.
    assert "cannot return a Reactome row" not in description
    assert "Deprecated" not in description
