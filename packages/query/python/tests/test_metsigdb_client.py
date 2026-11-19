"""MetSigDB through the clients: the hosted API answers, the local installation does not."""

import pytest
from fastapi.testclient import TestClient

from omnipath_api.queries import ReleaseQueries
from omnipath_api.release_database import build_served_database, publish_newest
from omnipath_api.server import create_app
from omnipath_client import ClientError, connect
from test_client import service  # noqa: F401  (the fixture)
from omnipath_core.fixtures import (
    METSIGDB_SETS,
    metsigdb_key,
    write_dataset,
    write_metsigdb_subsets,
)


@pytest.fixture
def api(tmp_path):
    root = write_dataset(tmp_path / "data")
    publish_newest(root, "1")
    write_metsigdb_subsets(root)
    build_served_database(root, "1", keep=None, memory_limit="1GB", threads=1)
    engine = ReleaseQueries(root)
    with TestClient(create_app(queries=engine)) as http:
        with connect(mode="api", api_url="http://testserver", http_client=http) as client:
            yield client
    engine.close()


def test_a_page_of_memberships_answers_as_the_route_does(api):
    page = api.metsigdb(resource="reactome", limit=2, total=True)
    assert (page["count"], page["has_more"], page["total"]) == (2, True, 4)
    assert {row["resource"] for row in page["rows"]} == {"Reactome"}


def test_filters_may_be_lists_and_booleans_reach_the_api_as_booleans(api):
    page = api.metsigdb(
        entity=["HMDB0000538", "C00031"],
        exclude_metabolites=["cofactors"],
        fields=["smiles"],
        total=True,
    )
    assert {row["entity"] for row in page["rows"]} == {metsigdb_key("glucose")}
    assert "smiles" in page["rows"][0]


def test_grouping_returns_the_sets_with_their_members(api):
    page = api.metsigdb(group=True, set="R-HSA-70171")
    (one,) = page["groups"][0]["sets"]
    assert one["returned"] == 3 and len(one["members"]) == 3


def test_the_sets_are_listed_without_their_members(api):
    listing = api.metsigdb_sets(total=True)
    assert listing["total"] == len(METSIGDB_SETS)
    assert "members" not in listing["sets"][0]


def test_a_value_the_api_refuses_is_a_value_error_with_its_reason(api):
    with pytest.raises(ValueError, match="smpdb"):
        api.metsigdb(resource="smpdb")


def test_a_filter_nobody_publishes_is_refused_before_a_request(api):
    with pytest.raises(TypeError, match="hmdb"):
        api.metsigdb(hmdb="HMDB0000538")


def test_a_release_without_the_subset_is_a_client_error(tmp_path):
    root = write_dataset(tmp_path / "data")
    publish_newest(root, "1")
    build_served_database(root, "1", keep=None, memory_limit="1GB", threads=1)
    engine = ReleaseQueries(root)
    try:
        with TestClient(create_app(queries=engine)) as http:
            with connect(mode="api", api_url="http://testserver", http_client=http) as client:
                with pytest.raises(ClientError, match="MetSigDB"):
                    client.metsigdb()
    finally:
        engine.close()


def test_a_local_installation_says_where_the_memberships_are(service):  # noqa: F811
    service.init(["alpha"])
    with connect(service.data) as local:
        with pytest.raises(ClientError, match="API"):
            local.metsigdb()


def test_the_exclusion_table_answers_as_the_route_does(api):
    page = api.exclusion_metabolites(exclusion_class="cofactors", total=True)
    assert page["total"] == page["count"] == 1
    assert {row["entity"] for row in page["rows"]} == {metsigdb_key("atp")}
    assert {row["exclusion_class"] for row in page["rows"]} == {"cofactors"}


def test_the_exclusion_report_answers_with_its_detail_as_data(api):
    page = api.exclusion_metabolite_report(kind="class_size", total=True)
    assert page["total"] == len(page["rows"])
    assert isinstance(page["rows"][0]["detail"], dict)


def test_an_exclusion_filter_nobody_publishes_is_refused_before_a_request(api):
    with pytest.raises(TypeError, match="hmdb"):
        api.exclusion_metabolites(hmdb="HMDB0000538")
    with pytest.raises(ValueError, match="xenon"):
        api.exclusion_metabolites(exclusion_class="xenon")


def test_a_local_installation_says_where_the_exclusion_table_is(service):  # noqa: F811
    service.init(["alpha"])
    with connect(service.data) as local:
        with pytest.raises(ClientError, match="API"):
            local.exclusion_metabolites()
