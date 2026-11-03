import sqlite3
from types import SimpleNamespace
from unittest.mock import patch

from omnipath_build.source_metadata import snapshot_source_metadata


def test_only_complete_successful_download_records_are_used(tmp_path):
    with sqlite3.connect(tmp_path / "cache.sqlite") as conn:
        conn.executescript("""
            CREATE TABLE main(id INTEGER, status INTEGER);
            CREATE TABLE attr_varchar(id INTEGER, name TEXT, value TEXT);
            CREATE TABLE attr_datetime(id INTEGER, name TEXT, value TEXT);
            INSERT INTO main VALUES (1,3),(2,4);
            INSERT INTO attr_varchar VALUES (1,'url','https://example.org/data'),
                (2,'url','https://example.org/data'),
                (2,'download_finished','2026-09-05 12:00:00');
            INSERT INTO attr_datetime VALUES (1,'download_finished','2026-05-01 12:00:00');
        """)
    ds = SimpleNamespace(
        qualified_module="pypath.inputs_v2.cellchat",
        raw_dataset=SimpleNamespace(download=SimpleNamespace(url="https://example.org/data")),
    )
    with patch(
        "omnipath_build.source_metadata.importlib.import_module", return_value=SimpleNamespace()
    ):
        metadata = snapshot_source_metadata("cellchat", [ds], tmp_path)
        assert metadata["downloaded_at"] == "2026-05-01 12:00:00"
        ds.raw_dataset.download = None
        assert snapshot_source_metadata("cellchat", [ds], tmp_path)["downloaded_at"] is None


def test_composite_and_dynamic_download_descriptors(tmp_path):
    from dataclasses import dataclass
    from omnipath_build.source_metadata import _download_urls

    @dataclass
    class Composite:
        links: object
        actions: object

    combined = Composite(
        SimpleNamespace(url="https://example.org/links"),
        SimpleNamespace(url="https://example.org/actions"),
    )
    assert _download_urls(combined) == (
        {"https://example.org/links", "https://example.org/actions"},
        True,
    )
    ds = SimpleNamespace(
        qualified_module="pypath.inputs_v2.stitch", raw_dataset=SimpleNamespace(download=combined)
    )
    with patch(
        "omnipath_build.source_metadata.importlib.import_module", return_value=SimpleNamespace()
    ):
        assert snapshot_source_metadata("stitch", [ds], tmp_path)["downloaded_at"] is None
    called = []
    combined.actions = SimpleNamespace(url=lambda: called.append(True))
    assert _download_urls(combined) == ({"https://example.org/links"}, False)
    assert not called
    assert _download_urls(object()) == (set(), False)
