"""The step-file runner: every construct, the same answers from Python and R."""

import json
import os
from pathlib import Path
import shutil
import subprocess

import duckdb
import pytest

from omnipath_client import steps

QUERIES = Path(__file__).parent / "step_queries"
R_RUNNER = Path(__file__).parents[2] / "r" / "tests" / "run_steps.R"
R_CLIENT = Path(__file__).parents[2] / "r"


def database(path):
    db = duckdb.connect(str(path))
    db.execute(
        "CREATE TABLE item AS SELECT * FROM (VALUES (1::BIGINT, 'a', 'Alpha'), (2, 'b', 'Beta'), "
        "(3, 'c', 'Gamma'), (4, 'd', 'Delta'), (5, 'e', 'Epsilon'), (6, 'f', 'Zeta')) AS t(_eid, key, label)"
    )
    db.execute("CREATE INDEX item_eid ON item (_eid)")
    db.close()
    return duckdb.connect(str(path), read_only=True)


def run(*args, **kwargs):
    """The runner without the shared macros: the toy database has no release tables."""
    return steps.run(*args, macros=False, **kwargs)


def query(name):
    return steps.parse((QUERIES / f"{name}.sql").read_text(), name)


def fragments(name):
    return (QUERIES / "fragments" / f"{name}.sql").read_text()


@pytest.fixture
def db(tmp_path):
    connection = database(tmp_path / "toy.duckdb")
    yield connection
    connection.close()


def test_lists_index_reads_slices_and_outputs(db):
    lookup = query("lookup")
    assert json.loads(run(db, lookup, {"keys": ["b", "a", "zz"]}, output="json")) == {
        "items": [{"key": "a", "label": "Alpha"}, {"key": "b", "label": "Beta"}],
        "total": 2,
    }
    page = run(db, lookup, {"keys": ["a", "b", "c"], "offset": 1, "limit": 2})
    assert page.to_pylist() == [{"key": "b", "label": "Beta"}, {"key": "c", "label": "Gamma"}]
    # An empty list matches nothing; the skipped step keeps its columns.
    assert json.loads(run(db, lookup, {"keys": []}, output="json")) == {"items": [], "total": 0}
    # One value stands for a list of one.
    assert run(db, lookup, {"keys": "c"}).to_pylist() == [{"key": "c", "label": "Gamma"}]
    # No temporary table is left behind.
    assert db.execute("SELECT count(*) FROM duckdb_tables() WHERE temporary").fetchone() == (0,)


def test_the_index_read_uses_the_index(tmp_path):
    path = tmp_path / "big.duckdb"
    with duckdb.connect(str(path)) as db:
        db.execute("CREATE TABLE item AS SELECT range AS _eid, 'k' || range AS key, 'L' AS label FROM range(200000)")
        db.execute("CREATE INDEX item_eid ON item (_eid)")
    with duckdb.connect(str(path), read_only=True) as db:
        db.execute("SET index_scan_max_count = 100000; SET index_scan_percentage = 0.5")
        lookup = query("lookup")
        run = steps.Run(db, lookup, steps.bind(lookup, {}), None, None)
        run.values["ids"] = [17, 123456, 199999]
        # The scan type is chosen when the statement runs: EXPLAIN alone shows a sequential scan.
        plan = db.execute("EXPLAIN ANALYZE " + run.render(lookup.steps[1].sql)).fetchall()[0][1]
        assert "Type: Index Scan" in plan


def test_repeat_grows_until_done(db):
    assert json.loads(run(db, query("grow"), {"target": 5}, output="json")) == {
        "n": 8, "found": 6, "enough": True
    }
    # At most max rounds.
    assert json.loads(run(db, query("grow"), {"target": 100}, output="json")) == {
        "n": 8, "found": 6, "enough": False
    }


def test_rows_tables_and_invalid_input(db):
    pairs = query("pairs")
    assert json.loads(run(db, pairs, {"labels": ["Alpha", "Gamma"]}, output="json")) == {
        "matched": ["a", "c"], "mode": "any"
    }
    assert json.loads(run(db, pairs, {"labels": ["Alpha"], "mode": "strict"}, output="json")) == {
        "matched": [], "mode": "strict"
    }
    with pytest.raises(steps.QueryError, match="mode must be any or strict"):
        run(db, pairs, {"labels": ["Alpha"], "mode": "other"})
    with pytest.raises(TypeError, match="unknown parameter"):
        run(db, pairs, {"label": "x"})
    with pytest.raises(TypeError, match="mode must be VARCHAR"):
        run(db, pairs, {"mode": 3})


def test_fragments_are_included(db):
    assert json.loads(run(db, query("shouted"), {"limit": 2}, output="json", fragments=fragments)) == [
        "ALPHA!", "BETA!"
    ]


def test_cached_steps_are_computed_once(db):
    text = """-- query: cached
-- param: n INTEGER = 2
-- output: table rows
-- step: ids AS list cache
SELECT _eid FROM item ORDER BY _eid LIMIT {{n}}
-- step: rows AS table cache
SELECT key FROM item WHERE _eid IN ({{ids*}}) ORDER BY key
"""
    store, computed = {}, []

    def cache(key, compute):
        if key not in store:
            computed.append(key[2])
            store[key] = compute()
        return store[key]

    cached = steps.parse(text)
    for _ in range(2):
        assert run(db, cached, {"n": 2}, cache=cache, scope="v1").to_pylist() == [{"key": "a"}, {"key": "b"}]
    assert computed == ["ids", "rows"]
    run(db, cached, {"n": 3}, cache=cache, scope="v1")
    assert computed == ["ids", "rows", "ids", "rows"]


def test_parse_errors():
    with pytest.raises(SyntaxError, match="no '-- output:'"):
        steps.parse("-- step: x AS value\nSELECT 1\n")
    with pytest.raises(SyntaxError, match="without '-- end'"):
        steps.parse("-- output: json x\n-- repeat: n = 1 next 2 until x max 2\n-- step: x AS value\nSELECT 1\n")
    with pytest.raises(SyntaxError, match="unknown type"):
        steps.parse("-- param: x TEXT\n-- output: json y\n")


def _r_available():
    if not shutil.which("Rscript"):
        return False
    probe = subprocess.run(
        ["Rscript", "-e", "for (p in c('DBI','duckdb','jsonlite','pkgload')) if (!requireNamespace(p, quietly=TRUE)) quit(status=1)"],
        capture_output=True, env=_r_env(),
    )
    return probe.returncode == 0


def _r_env():
    return {**os.environ, "R_LIBS_USER": os.environ.get("R_LIBS_USER", str(Path.home() / "Library" / "R" / "library-4.6"))}


@pytest.mark.skipif(not _r_available(), reason="R or its packages are not installed")
def test_r_runner_answers_like_python(tmp_path):
    database(tmp_path / "toy.duckdb").close()
    cases = [
        (name, example)
        for name in ("lookup", "grow", "pairs", "shouted")
        for example in query(name).examples
    ]
    (tmp_path / "cases.json").write_text(json.dumps([[name, example] for name, example in cases]))
    process = subprocess.run(
        ["Rscript", str(R_RUNNER), str(R_CLIENT), str(tmp_path / "toy.duckdb"), str(QUERIES), str(tmp_path / "cases.json")],
        capture_output=True, text=True, env=_r_env(),
    )
    assert process.returncode == 0, process.stderr[-3000:]
    answers = [json.loads(text) for text in json.loads(process.stdout[process.stdout.index("[") :])]
    db = duckdb.connect(str(tmp_path / "toy.duckdb"), read_only=True)
    for (name, example), answer in zip(cases, answers):
        try:
            expected = json.loads(run(db, query(name), example, output="json", fragments=fragments))
        except steps.QueryError as exc:
            expected = {"error": str(exc)}
        assert answer == expected, (name, example)
