"""TEMPORARY: the old Python molecular context against the step file (removed after)."""
import json

from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.release_database import serve_newest
from test_client import service  # noqa: F401
from zz_jsondiff import diff


def norm(x):
    return json.loads(json.dumps(x, default=str))


def test_old_and_new_molecular_context(service):
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    with api.query_scope():
        rows = api._db.execute(f"SELECT DISTINCT entity_key, entity_type FROM {api._served().alias}.main.entity").fetchall()
        refs = [r[0] for r in api._db.execute(f"SELECT DISTINCT reference_entity_key FROM {api._served().alias}.main.entity WHERE reference_entity_key IS NOT NULL").fetchall()]
    cases = [(key, {}) for key, _ in rows] + [(key, {"view": "product"}) for key, t in rows if t in ("protein", "transcript")]
    cases += [("gene:" + r, {}) for r in refs] + [("missing", {}), (rows[0][0], {"limit": 1}), (rows[0][0], {"offset": 1})]
    for entity, kwargs in cases:
        with api.query_scope():
            try:
                old = norm(api._old_get_molecular_context(entity, None, **kwargs))
            except ValueError as exc:
                old = str(exc)
            try:
                new = norm(api.get_molecular_context(entity, None, **kwargs))
            except ValueError as exc:
                new = str(exc)
        assert new == old, (entity, kwargs, diff(old, new) if isinstance(old, dict) and isinstance(new, dict) else (old, new))
    print(len(cases), "cases")
