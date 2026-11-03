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


@pytest.mark.parametrize("materialize", [True, False])
def test_r_and_python_share_directories_and_answers(service, tmp_path, materialize):  # noqa: F811
    service.init(["alpha", "beta", "gamma"], include_evidence=True, materialize=materialize)
    r_dir = tmp_path / "r"
    questions = tmp_path / "questions.json"
    questions.write_text(json.dumps({**QUESTIONS, "materialize": materialize}))
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
