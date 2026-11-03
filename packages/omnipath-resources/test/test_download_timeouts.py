"""Shared download timeouts also apply to legacy dm.download callers."""
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest
from pypath.share import downloads


def test_default_timeouts_and_explicit_override(monkeypatch):
    calls = []
    monkeypatch.setattr(downloads.DownloadManager, 'download', lambda self, url, *args, **kw: calls.append(kw))
    manager = object.__new__(downloads._PypathDownloadManager)
    manager.download('https://example.invalid/data')
    assert calls[-1]['connecttimeout'] == 30
    assert calls[-1]['timeout'] == 120
    assert calls[-1]['retries'] == 1
    manager.download('https://example.invalid/data', timeout=600)
    assert calls[-1]['timeout'] == 600
    monkeypatch.setenv('PYPATH_READ_TIMEOUT', '0.1')
    manager.download('https://example.invalid/data')
    assert calls[-1]['timeout'] == 0.1


def test_unresponsive_http_source_times_out(tmp_path, monkeypatch):
    class StalledHandler(BaseHTTPRequestHandler):
        def do_GET(self):
            time.sleep(1)
        def log_message(self, *args):
            pass
    server = ThreadingHTTPServer(('127.0.0.1', 0), StalledHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    monkeypatch.setenv('PYPATH_READ_TIMEOUT', '0.1')
    manager = downloads._PypathDownloadManager(path=str(tmp_path / 'cache'), config={'backend': 'requests'})
    start = time.monotonic()
    try:
        with pytest.raises(Exception):
            manager.download(f'http://127.0.0.1:{server.server_port}/stalled', dest=str(tmp_path / 'result'))
        assert time.monotonic() - start < 1
    finally:
        server.shutdown()
        server.server_close()
