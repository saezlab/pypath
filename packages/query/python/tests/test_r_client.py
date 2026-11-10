"""The R client (packages/query/r) and this client share data directories and answers.

R installs the fixtures into one directory and this client into another; the R
client then answers the same questions on both, and this client on both. Answers must
agree, and neither client may rewrite the other's database. Skipped without R and the
R packages DBI, duckdb, httr2, jsonlite, pkgload and testthat.
"""

import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from omnipath_client import connect
from omnipath_core.fixtures import fixture_key as key

from test_client import service  # noqa: F401  (the fixture)

R_CLIENT = Path(__file__).resolve().parents[2] / "r"
QUESTIONS = {
    "search": [
        {"query": "TP53", "filters": {}},
        {"query": "tp", "filters": {}},
        {"query": "caffeine", "filters": {}},
        {"query": "uniprot:P04637", "filters": {}},
        {"query": "TP53", "filters": {"sources": "beta"}},
        {"query": "brien", "filters": {}},
    ],
    "relations": [
        {"entities": "TP53", "filters": {}},
        {"entities": None, "filters": {"sources": "alpha", "predicates": "affects"}},
        {"entities": ["caffeine", "O'Brien"], "filters": {"gene_mode": "both"}},
        {"entities": "gene:entrez:7157", "filters": {"taxonomy_ids": "9606"}},
    ],
    "entity": [key("P04637"), key("CHEBI:27732"), "0" * 64],
}


def r_available() -> bool:
    if not shutil.which("Rscript"):
        return False
    check = subprocess.run(
        ["Rscript", "-e", "for (p in c('DBI','duckdb','httr2','jsonlite','pkgload','testthat')) "
         "if (!requireNamespace(p, quietly = TRUE)) quit(status = 1)"],
        capture_output=True,
    )
    return check.returncode == 0


pytestmark = pytest.mark.skipif(not r_available(), reason="R or its packages are not installed")


def clean(value):
    """JSON values without nulls, so that R's and Python's renderings compare."""
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items() if v is not None}
    if isinstance(value, list):
        return [clean(v) for v in value]
    return value


def python_answers(directory):
    def as_list(value):
        return [value] if isinstance(value, str) else value

    with connect(directory) as db:
        search = [
            db.search(q["query"], **q["filters"]).to_arrow_table().to_pylist() for q in QUESTIONS["search"]
        ]
        relations = []
        for q in QUESTIONS["relations"]:
            rows = db.relations(q["entities"], **q["filters"]).to_arrow_table().to_pylist()
            relations.append(rows)
        entity = [db.entity(k) for k in QUESTIONS["entity"]]
        layer = db.sql("SELECT query_layer FROM database_info").fetchone()[0]
    return dict(search=search, relations=relations, entity=entity, query_layer=layer)


def from_r(answers):
    return dict(
        search=[json.loads(rows) for rows in answers["search"]],
        relations=[json.loads(rows) for rows in answers["relations"]],
        entity=[json.loads(value) for value in answers["entity"]],
        query_layer=answers["query_layer"],
    )


def same(python, r):
    assert python["query_layer"] == r["query_layer"]
    assert [clean(rows) for rows in python["search"]] == [clean(rows) for rows in r["search"]]
    def canonical(rows):
        return sorted(json.dumps(clean(row), sort_keys=True) for row in rows)

    for p_rows, r_rows in zip(python["relations"], r["relations"]):
        assert canonical(json.loads(json.dumps(p_rows, default=str))) == canonical(r_rows)
    assert clean(python["entity"]) == clean(r["entity"])


def test_r_and_python_share_directories_and_answers(service, tmp_path):  # noqa: F811
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    r_dir = tmp_path / "r"
    questions = tmp_path / "questions.json"
    questions.write_text(json.dumps(QUESTIONS))
    database = service.data / "omnipath.duckdb"
    before = database.stat().st_ino
    run = subprocess.run(
        ["Rscript", str(R_CLIENT / "tests" / "cross_language.R"), str(R_CLIENT), str(r_dir),
         str(service.data), str(questions)],
        capture_output=True, text=True, env={**os.environ, "R_LIBS_USER": os.environ.get(
            "R_LIBS_USER", str(Path.home() / "Library" / "R" / "library-4.6"))},
    )
    assert run.returncode == 0, run.stderr[-3000:]
    answers = json.loads(run.stdout[run.stdout.index('{"r"'):])
    # R read the Python client's database as it was.
    assert database.stat().st_ino == before
    python = python_answers(service.data)
    same(python, from_r(answers["python"]))
    same(python, from_r(answers["r"]))
    # The Python client reads the R client's database as it is, with the same answers.
    r_database = r_dir / "omnipath.duckdb"
    r_before = r_database.stat().st_ino
    same(python_answers(r_dir), from_r(answers["r"]))
    assert r_database.stat().st_ino == r_before


def r_http_available() -> bool:
    """R's curl package loads (httr2 needs it for real requests; a CRAN binary built
    for another R fails to load)."""
    check = subprocess.run(
        ["Rscript", "-e", "loadNamespace('curl')"], capture_output=True,
        env={**os.environ, "R_LIBS_USER": os.environ.get(
            "R_LIBS_USER", str(Path.home() / "Library" / "R" / "library-4.6"))},
    )
    return check.returncode == 0


@pytest.mark.skipif(not r_http_available(), reason="R's curl package does not load")
def test_r_api_mode_answers_like_local_mode(service, tmp_path):  # noqa: F811
    """R in API mode, against the API serving the Python client's directory over HTTP."""
    import socket
    import threading
    import time

    import uvicorn
    from omnipath_api.queries import ReleaseQueries
    from omnipath_api.release_database import serve_newest
    from omnipath_api.server import create_app

    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    engine = ReleaseQueries(service.data)
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    server = uvicorn.Server(uvicorn.Config(create_app(queries=engine), host="127.0.0.1", port=port, log_level="warning"))
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()
    try:
        for _ in range(100):
            if server.started:
                break
            time.sleep(0.05)
        questions = tmp_path / "questions.json"
        questions.write_text(json.dumps({**QUESTIONS, "api_url": f"http://127.0.0.1:{port}"}))
        run = subprocess.run(
            ["Rscript", str(R_CLIENT / "tests" / "cross_language.R"), str(R_CLIENT), str(tmp_path / "r"),
             str(service.data), str(questions)],
            capture_output=True, text=True, env={**os.environ, "R_LIBS_USER": os.environ.get(
                "R_LIBS_USER", str(Path.home() / "Library" / "R" / "library-4.6"))},
        )
        assert run.returncode == 0, run.stderr[-3000:]
        answers = json.loads(run.stdout[run.stdout.index('{"r"'):])["api"]
        assert all(answers["entity_identical"]) and answers["sql_refused"]
        python = python_answers(service.data)
        api = from_r({**answers, "entity": ["null"] * len(QUESTIONS["entity"]), "query_layer": python["query_layer"]})
        assert [clean(rows) for rows in python["search"]] == [clean(rows) for rows in api["search"]]
        def canonical(rows):
            return sorted(json.dumps(clean(row), sort_keys=True) for row in rows)
        for p_rows, a_rows in zip(python["relations"], api["relations"]):
            assert canonical(json.loads(json.dumps(p_rows, default=str))) == canonical(a_rows)
    finally:
        server.should_exit = True
        thread.join(10)
        engine.close()
