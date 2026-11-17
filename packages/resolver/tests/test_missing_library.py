"""A resolver with no identity library says so: entities then keep their own identifiers."""

import logging
from pathlib import Path
import sys
import types

import pytest

try:
    import omnipath_resolver  # noqa: F401
except ImportError:  # the native extension is not built here: the Python half is enough
    stub = types.ModuleType("omnipath_resolver._omnipath_resolver")
    stub.resolve_precomputed_batch = stub.resolve_molecular_batch = lambda *a, **k: None
    sys.modules["omnipath_resolver._omnipath_resolver"] = stub
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omnipath_resolver.resolver import EntityResolver, locate_library_dir


@pytest.fixture(autouse=True)
def no_library_in_the_environment(monkeypatch):
    monkeypatch.delenv("OMNIPATH_LIBRARY_DIR", raising=False)


def test_no_library_anywhere_is_none():
    assert locate_library_dir() is None


def test_a_resolver_without_a_library_warns_how_to_give_one(caplog):
    with caplog.at_level(logging.WARNING):
        resolver = EntityResolver()
    assert resolver.library_dir is None
    messages = " ".join(record.getMessage() for record in caplog.records)
    assert "No reference library" in messages
    assert "OMNIPATH_LIBRARY_DIR" in messages


def test_a_resolver_can_be_told_to_require_one():
    with pytest.raises(FileNotFoundError, match="OMNIPATH_LIBRARY_DIR"):
        EntityResolver(require_library=True)


def test_a_resolver_given_a_data_root_without_a_library_warns_where_it_looked(caplog, tmp_path):
    with caplog.at_level(logging.WARNING):
        EntityResolver(data_root=tmp_path)
    assert str(tmp_path) in " ".join(record.getMessage() for record in caplog.records)
