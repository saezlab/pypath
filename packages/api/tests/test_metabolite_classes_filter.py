"""The relation filter metabolite_classes through the API: the request's filter reaches
the relation queries, and a release without the metabolite_class subset refuses it.

The fixture release (omnipath_core.fixtures.write_dataset) holds proteins only, so the
table written here names TP53: what the filter tests is the row, not the entity type."""

import duckdb
import pytest
from fastapi.testclient import TestClient

from metsigdb_served import serve
from omnipath_api.queries.release import SERVING_DIR, ReleaseQueries
from omnipath_api.queries.services.export import relation_params
from omnipath_api.server import create_app
from omnipath_core.fixtures import fixture_key

TP53 = fixture_key("P04637")


def _classify(root):
    with duckdb.connect(str(root / SERVING_DIR / "1.duckdb")) as db:
        db.execute("""CREATE TABLE metabolite_class (
            entity_key VARCHAR NOT NULL, label VARCHAR, inchikey VARCHAR, hmdb VARCHAR,
            metabolite_class VARCHAR NOT NULL, match_kind VARCHAR NOT NULL,
            via_entity_key VARCHAR, supported_by VARCHAR[] NOT NULL, build_id VARCHAR NOT NULL)""")
        db.execute(
            "INSERT INTO metabolite_class VALUES (?, 'TP53', NULL, NULL, 'hmdb', 'entity', NULL, ['test'], 'b1')",
            [TP53],
        )


@pytest.fixture(params=["bare", "classed"])
def api(request, tmp_path):
    root = serve(tmp_path, subsets=False)
    if request.param == "classed":
        _classify(root)
    engine = ReleaseQueries(root)
    with TestClient(create_app(queries=engine)) as client:
        client.kind = request.param
        yield client
    engine.close()


def _search(api, **filters):
    return api.post("/relations/search", json={"filters": filters, "limit": 100})


def _ends(response):
    return [(r["subjectEntity"]["entityPk"], r["objectEntity"]["entityPk"]) for r in response.json()["rows"]]


def test_the_filter_reaches_the_relation_queries():
    assert relation_params({"metaboliteClasses": ["hmdb"]})["metabolite_classes"] == ["hmdb"]
    assert relation_params({"metabolite_classes": ["hmdb", "chebi_role"]})["metabolite_classes"] == [
        "hmdb", "chebi_role"]
    assert relation_params({})["metabolite_classes"] is None


def test_the_filter_keeps_the_relations_with_an_end_in_the_bucket(api):
    every = _search(api)
    assert every.status_code == 200 and every.json()["total"] > 0
    response = _search(api, metaboliteClasses=["hmdb"])
    if api.kind == "bare":
        assert response.status_code == 400
        assert "no metabolite_class subset" in response.json()["detail"]
        return
    assert response.status_code == 200
    kept = _ends(response)
    assert kept and all(TP53 in ends for ends in kept)
    assert sorted(kept) == sorted(e for e in _ends(every) if TP53 in e)
    assert response.json()["total"] == len(kept)
    unknown = _search(api, metaboliteClasses=["hmbd"])
    assert unknown.status_code == 200 and unknown.json()["total"] == 0 and unknown.json()["rows"] == []
    facets = api.post("/relations/scoped-facets", json={"metaboliteClasses": ["hmbd"]})
    assert facets.status_code == 200 and facets.json() == []
