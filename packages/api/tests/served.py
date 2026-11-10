"""The API's engine over fixture resources: every published release gets its database.

Tests write resources under a data root; ``served_queries`` builds the database of each
published release that has none (publishing the newest resources as a release first
when nothing is published) and returns the engine. ``republish`` publishes the current
newest resources as a new release with its database, as a server does after builds.
"""

from __future__ import annotations

from pathlib import Path

from omnipath_api.queries.release import SERVING_DIR, ReleaseQueries
from omnipath_api.queries.store.inventory import ReleaseStore
from omnipath_api.release_database import build_served_database, publish_newest


def serve_releases(root) -> None:
    root = Path(root)
    store = ReleaseStore(root)
    if not (root / "releases").is_dir() or not store.list():
        try:
            publish_newest(root)
        except ValueError as exc:
            if "at least one resource" in str(exc):
                return  # nothing built yet: nothing is served
            if "sample builds" not in str(exc):
                raise
            publish_newest(root, partial_resources=True)
    for manifest in store.list():
        if not (root / SERVING_DIR / f"{manifest['version']}.duckdb").is_file():
            build_served_database(root, manifest["version"], keep=None, memory_limit="1GB", threads=1)


def served_queries(root="data", *args, **kwargs) -> ReleaseQueries:
    data_root = kwargs.get("data_root", root)
    serve_releases(data_root)
    return ReleaseQueries(root, *args, **kwargs) if "data_root" not in kwargs else ReleaseQueries(*args, **kwargs)


def republish(root, engine: ReleaseQueries | None = None) -> str:
    """A new release of the newest resources, with its database; the engine sees it."""
    try:
        version = publish_newest(root)
    except ValueError as exc:
        if "sample builds" not in str(exc):
            raise
        version = publish_newest(root, partial_resources=True)
    build_served_database(root, version, keep=None, memory_limit="1GB", threads=1)
    if engine is not None:
        engine.reload_resources()
    return version
