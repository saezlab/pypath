"""What ``/sets/metsigdb`` answers for a release that has no subset."""

import pytest
from fastapi.testclient import TestClient

from metsigdb_served import serve
from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.server import create_app


@pytest.fixture
def bare(tmp_path):
    root = serve(tmp_path, subsets=False)
    engine = ReleaseQueries(root)
    with TestClient(create_app(queries=engine)) as client:
        yield client
    engine.close()


@pytest.mark.parametrize(
    "path, params",
    [
        ("/sets/metsigdb", {}),
        ("/sets/metsigdb", {"group": "true"}),
        ("/sets/metsigdb", {"exclude_metabolites": "cofactors"}),
        ("/sets/metsigdb/sets", {}),
    ],
)
def test_a_release_without_the_subset_answers_503_and_says_why(bare, path, params):
    response = bare.get(path, params=params)
    assert response.status_code == 503
    assert "MetSigDB" in response.json()["detail"]


def test_a_request_the_route_refuses_is_refused_before_the_subset_is_looked_for(bare):
    assert bare.get("/sets/metsigdb", params={"resource": "smpdb"}).status_code == 400
