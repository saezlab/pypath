"""PTFI and RaMP fetch through the shared download manager, like every other input."""

import io
import json

import pytest
import requests

from omnipath_resources.inputs_v2 import ptfi, rampdb
from omnipath_resources.share import downloads


class Manager:
    """A download manager that serves canned bodies and records what was asked."""

    def __init__(self, bodies):
        self.bodies = bodies
        self.calls = []

    def download(self, url, *args, **kwargs):
        self.calls.append((url, kwargs))
        body = self.bodies.get(url)
        return None if body is None else io.BytesIO(body.encode())


@pytest.fixture
def manager(monkeypatch):
    def install(bodies):
        found = Manager(bodies)
        monkeypatch.setattr(downloads, "get_download_manager", lambda: found)
        return found

    def no_direct_http(*args, **kwargs):
        raise AssertionError("an input must fetch through the download manager")

    monkeypatch.setattr(requests, "get", no_direct_http)
    monkeypatch.setattr(requests, "Session", no_direct_http)
    return install


LISTING = """<html><body>
<a title="ramp_db_2.5.0.sqlite.gz">ramp_db_2.5.0.sqlite.gz</a>
<a title="ramp_db_2.6.1.sqlite.gz">ramp_db_2.6.1.sqlite.gz</a>
<a title="notes.txt">notes.txt</a>
</body></html>"""


def test_the_latest_ramp_release_is_read_through_the_download_manager(manager):
    found = manager({rampdb.BASE_URL % "main": LISTING})
    version, url = rampdb.get_ramp_latest_ver()
    assert version == "2.6.1"
    assert url == rampdb.BASE_URL % "main" + "ramp_db_2.6.1.sqlite.gz"
    (call,) = found.calls
    # A directory listing must not be served from an earlier cache entry.
    assert call[1].get("force_download") is True and call[1].get("dest") is False


def test_a_ramp_listing_that_cannot_be_downloaded_is_an_error(manager):
    manager({})
    with pytest.raises(RuntimeError, match="no data"):
        rampdb.get_ramp_latest_ver()


def test_ptfi_specimens_are_fetched_through_the_download_manager(manager, tmp_path, monkeypatch):
    monkeypatch.setattr(ptfi, "DELAY_SECONDS", 0)
    ids = [101, 102]
    bodies = {
        ptfi.BASE_API_URL.format(specimen_id=f"specimen_FOODON_{i}"): json.dumps({"id": i}) for i in ids
    }
    found = manager(bodies)
    out = ptfi.download_ptfi_data(output_dir=tmp_path, specimen_ids=ids, allow_download=True)
    assert out == tmp_path
    assert json.loads((tmp_path / "specimen_FOODON_101.json").read_text()) == {"id": 101}
    assert {call[0] for call in found.calls} == set(bodies)
    assert all("User-Agent" in call[1]["headers"] for call in found.calls)


def test_a_ptfi_specimen_that_fails_does_not_stop_the_others(manager, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(ptfi, "DELAY_SECONDS", 0)
    manager(
        {
            ptfi.BASE_API_URL.format(specimen_id="specimen_FOODON_2"): json.dumps({"id": 2}),
            ptfi.BASE_API_URL.format(specimen_id="specimen_FOODON_3"): "not json",
        }
    )
    ptfi.download_ptfi_data(output_dir=tmp_path, specimen_ids=[1, 2, 3], allow_download=True)
    assert [p.name for p in sorted(tmp_path.glob("*.json"))] == ["specimen_FOODON_2.json"]
    assert "1 downloaded, 0 skipped, 2 failed" in capsys.readouterr().out


def test_an_existing_specimen_is_not_fetched_again(manager, tmp_path, monkeypatch):
    monkeypatch.setattr(ptfi, "DELAY_SECONDS", 0)
    (tmp_path / "specimen_FOODON_1.json").write_text("{}")
    found = manager({})
    ptfi.download_ptfi_data(output_dir=tmp_path, specimen_ids=[1], allow_download=True)
    assert found.calls == []
