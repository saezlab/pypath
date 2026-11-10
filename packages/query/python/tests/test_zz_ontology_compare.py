"""TEMPORARY: the old Python ontology views against the step files (removed after)."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[3] / "api" / "tests"))

from served import served_queries  # noqa: E402


def norm(x):
    return json.loads(json.dumps(x, default=str))


def test_old_and_new_ontology(tmp_path):
    from omnipath_core.fixtures import write_dataset

    root = write_dataset(tmp_path / "data")
    api = served_queries(data_root=root)
    with api.query_scope():
        terms = [r[0] for r in api._db.execute(f"SELECT DISTINCT identifier FROM {api._served().alias}.main.entity WHERE entity_type = 'ontology_class'").fetchall()]
        keys = [r[0] for r in api._db.execute(f"SELECT DISTINCT entity_key FROM {api._served().alias}.main.entity").fetchall()]
    payloads = [{}, {"query": "a"}, {"limit": 1}, {"limit": 1, "offset": 1}, {"scoped": True}, {"entityPks": keys[:2]},
                {"termIds": terms[:1]}, {"selectionScope": {"entityPks": ["TP53"]}}, {"ontologyIds": ["go"]},
                {"prefixes": ["GO"]}, {"query": "zzz"}]
    for payload in payloads:
        for counts in (False, True):
            with api.query_scope():
                old = api._old_search_ontology(dict(payload), None, counts)
                new = api.search_ontology(dict(payload), None, counts)
            assert norm(new) == norm(old), (payload, counts)
    for ids in [terms, terms[:1] + ["nope"], keys[:3], []]:
        with api.query_scope():
            assert norm(api.get_terms_info(ids)) == norm(api._old_get_terms_info(ids)), ids
    print(len(terms), "terms")


def test_old_and_new_tree(tmp_path):
    from omnipath_core.fixtures import write_dataset

    root = write_dataset(tmp_path / "data")
    api = served_queries(data_root=root)
    with api.query_scope():
        rows = api._db.execute(f"SELECT DISTINCT identifier, label, entity_key FROM {api._served().alias}.main.entity").fetchall()
    for ident, label, key in rows + [("nope", None, None)]:
        for term in [ident, label, key, f"go|{ident}"]:
            if not term:
                continue
            for call in ("get_ontology_tree", "get_ontology_children"):
                args = ([term],) if call == "get_ontology_tree" else (term,)
                with api.query_scope():
                    try:
                        old = getattr(api, "_old_" + call)(*args)
                    except ValueError as exc:
                        old = str(exc)
                    try:
                        new = getattr(api, call)(*args)
                    except ValueError as exc:
                        new = str(exc)
                assert norm(new) == norm(old), (call, term)
