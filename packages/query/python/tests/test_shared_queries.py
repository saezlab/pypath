"""Every query answers alike wherever it runs: the Python client, the R client and the
API (its generic route and the client's API mode), on one installation of the
fixtures. Each query runs with each of its examples (``-- example:`` in its file), and
with questions about the fixtures' entities, so that answers have rows."""

import json
import subprocess

import pytest
from fastapi.testclient import TestClient

from omnipath_api.queries import ReleaseQueries
from omnipath_api.release_database import serve_newest
from omnipath_api.server import create_app
from omnipath_client import connect, steps
from omnipath_core.fixtures import fixture_key as key

from test_client import service  # noqa: F401  (the fixture)
from test_steps import R_CLIENT, R_RUNNER, _r_available, _r_env

TP53, CAFFEINE = key("P04637"), key("CHEBI:27732")
# Questions about the fixtures, on top of the files' examples (whose identifiers are a
# full release's).
FIXTURE_CASES = [
    ("entity_details", {"entity_key": TP53}),
    ("entity_core", {"entity_key": CAFFEINE}),
    ("entity_relationships", {"entity_keys": [TP53], "limit": 1}),
    ("entity_evidence", {"entity_key": TP53}),
    ("entities_by_keys", {"tokens": [TP53, CAFFEINE]}),
    ("entity_search", {"query": "TP53"}),
    ("entity_search", {"query": "tp", "limit": 1}),
    ("resolve_keys", {"tokens": ["TP53", "gene:entrez:7157"]}),
    ("relation_page", {"entities": [TP53], "limit": 1, "details": True}),
    ("relation_page", {"sources": ["beta"], "records": True}),
    ("relations", {"entities": ["TP53"]}),
    ("relation_facets", {"entities": [TP53]}),
    ("relation_counts", {}),
    ("entity_facets", {"query": "TP53"}),
    ("entity_groups", {"query": "TP53", "strategy": "auto"}),
    ("selection_scope", {"entity_pks": [TP53]}),
    ("stats_entity_types", {}),
]


def cases():
    found = [
        (name, example) for name in steps.names() for example in steps.load(name).examples
    ]
    for name, params in FIXTURE_CASES:
        assert set(params) <= set(steps.load(name).params), name
        found.append((name, params))
    return found


def local_answer(db, name, params):
    try:
        text = db.query(name, output="json", **params)
    except ValueError as exc:
        return {"error": str(exc)}
    return None if text is None else json.loads(text)


def api_answer(http, name, params):
    response = http.post(f"/query/{name}", params={"format": "json"}, json=params)
    if response.status_code == 400:
        return {"error": response.json()["detail"]}
    assert response.status_code == 200, (name, params, response.text[:300])
    return json.loads(response.text)


@pytest.fixture
def installed(service):  # noqa: F811
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    return service.data


def test_every_query_has_examples():
    assert not [name for name in steps.names() if not steps.load(name).examples]


def test_queries_answer_alike_locally_and_in_the_api(installed):
    engine = ReleaseQueries(installed)
    try:
        with TestClient(create_app(queries=engine)) as http, connect(installed) as local, connect(
            mode="api", api_url="http://testserver", http_client=http
        ) as remote:
            answered = 0
            for name, params in cases():
                outputs = steps.load(name).outputs
                if "json" in outputs:
                    want = local_answer(local, name, params)
                    assert api_answer(http, name, params) == want, (name, params)
                    if "error" in (want or {}):
                        continue
                    assert json.loads(remote.query(name, output="json", **params) or "null") == want
                    answered += bool(want)
                # A table answer comes back from the API as the same table.
                if "table" in outputs:
                    table = local.query(name, **params).to_arrow_table()
                    got = remote.query(name, **params).to_arrow_table()
                    assert got.schema.names == table.schema.names, name
                    assert got.to_pylist() == table.to_pylist(), (name, params)
                    answered += bool(table.num_rows)
            assert answered > len(FIXTURE_CASES)
    finally:
        engine.close()


@pytest.mark.skipif(not _r_available(), reason="R or its packages are not installed")
def test_r_answers_every_query_like_python(installed, tmp_path):
    found = [(name, params) for name, params in cases() if "json" in steps.load(name).outputs]
    (tmp_path / "cases.json").write_text(json.dumps([[name, params] for name, params in found]))
    process = subprocess.run(
        ["Rscript", str(R_RUNNER), str(R_CLIENT), str(installed / "omnipath.duckdb"), "package",
         str(tmp_path / "cases.json")],
        capture_output=True, text=True, env=_r_env(),
    )
    assert process.returncode == 0, process.stderr[-3000:]
    answers = json.loads(process.stdout[process.stdout.index("[") :])
    with connect(installed) as local:
        for (name, params), text in zip(found, answers):
            assert json.loads(text) == local_answer(local, name, params), (name, params)
