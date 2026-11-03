"""Local paths and read-only HTTP locations for pinned release artifacts.

HTTP bytes are streamed or read by range; no persistent Parquet cache is created.
"""

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from functools import cache
import io
import os
import ssl
from time import perf_counter
from pathlib import Path
from urllib.parse import unquote, urlsplit
from urllib.request import Request, urlopen

import certifi

Location = str | Path
TIMEOUT = 60


@dataclass
class HTTPTransfer:
    """Actual Python HTTP response bytes; excludes DuckDB's own HTTP client."""

    requests: int = 0
    bytes_read: int = 0
    seconds: float = 0.0


_transfer: ContextVar[HTTPTransfer | None] = ContextVar("http_transfer", default=None)


@contextmanager
def measure_http_transfer():
    """Measure one validation phase without recording URLs or source contents."""
    stats = HTTPTransfer()
    token = _transfer.set(stats)
    started = perf_counter()
    try:
        yield stats
    finally:
        stats.seconds = perf_counter() - started
        _transfer.reset(token)


class _CountedResponse:
    def __init__(self, response, stats):
        self.response, self.stats = response, stats

    def read(self, *args):
        data = self.response.read(*args)
        self.stats.bytes_read += len(data)
        return data

    def __getattr__(self, name):
        return getattr(self.response, name)

    def __enter__(self):
        self.response.__enter__()
        return self

    def __exit__(self, *args):
        return self.response.__exit__(*args)


def is_remote(value: Location) -> bool:
    return isinstance(value, str) and value.startswith(("https://", "http://"))


def validate_url(value: str) -> str:
    parsed = urlsplit(value)
    if parsed.scheme != "https" and not (
        parsed.scheme == "http" and parsed.hostname in {"localhost", "127.0.0.1", "::1"}
    ):
        raise ValueError("Remote artifacts require HTTPS (HTTP is allowed on loopback for tests)")
    if not parsed.hostname or parsed.username or parsed.password or parsed.query or parsed.fragment:
        raise ValueError("Artifact URLs require a host and no credentials, query or fragment")
    path = unquote(parsed.path)
    if "\x00" in path or "\\" in path or any(part in {".", ".."} for part in path.split("/")):
        raise ValueError("Artifact URL must be a literal path without traversal")
    return value


def join_location(root: Location, *parts: str) -> Location:
    if is_remote(root):
        return validate_url(root.rstrip("/") + "/" + "/".join(parts))
    return Path(root).joinpath(*parts)


@cache
def ssl_context() -> ssl.SSLContext:
    """Certificate authorities for HTTPS: SSL_CERT_FILE where it is set, else certifi's.

    The standalone Python builds that uv installs do not find the system's
    certificates on every host (NixOS): reading a release there failed with
    CERTIFICATE_VERIFY_FAILED.
    """
    return ssl.create_default_context(cafile=None if os.environ.get("SSL_CERT_FILE") else certifi.where())


def open_http(url: str, *, method: str = "GET", headers: dict | None = None):
    request = Request(
        validate_url(url), method=method, headers={"Accept-Encoding": "identity", **(headers or {})}
    )
    response = urlopen(request, timeout=TIMEOUT, context=ssl_context())
    try:
        validate_url(response.geturl())
        if response.headers.get("Content-Encoding", "identity") != "identity":
            raise OSError("Artifact server must return uncompressed HTTP bytes")
    except Exception:
        response.close()
        raise
    stats = _transfer.get()
    if stats is not None:
        stats.requests += 1
        return _CountedResponse(response, stats)
    return response


def read_bytes(location: Location, *, limit: int = 8 * 1024 * 1024) -> bytes:
    if not is_remote(location):
        return Path(location).read_bytes()
    with open_http(location) as response:
        content = response.read(limit + 1)
    if len(content) > limit:
        raise ValueError("Release metadata exceeds the 8 MiB limit")
    return content


class HTTPRangeFile(io.RawIOBase):
    """Seekable Parquet metadata reader, requiring exact HTTP byte ranges."""

    def __init__(self, url: str):
        super().__init__()
        self.url = validate_url(url)
        with open_http(url, method="HEAD") as response:
            value = response.headers.get("Content-Length")
            if value is None or not value.isdigit():
                raise OSError("Artifact server must provide a valid Content-Length")
            self.size = int(value)
        self.position = 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=io.SEEK_SET):
        position = (
            offset
            if whence == io.SEEK_SET
            else (self.position + offset if whence == io.SEEK_CUR else self.size + offset)
        )
        if whence not in {io.SEEK_SET, io.SEEK_CUR, io.SEEK_END} or position < 0:
            raise ValueError("Invalid HTTP artifact seek")
        self.position = position
        return position

    def read(self, size=-1):
        if self.closed:
            raise ValueError("Read from closed HTTP artifact")
        size = self.size - self.position if size < 0 else min(size, self.size - self.position)
        if size <= 0:
            return b""
        start, end = self.position, self.position + size - 1
        with open_http(self.url, headers={"Range": f"bytes={start}-{end}"}) as response:
            expected = f"bytes {start}-{end}/{self.size}"
            if response.status != 206 or response.headers.get("Content-Range") != expected:
                raise OSError("Artifact server must support exact HTTP byte ranges")
            data = response.read(size + 1)
        if len(data) != size:
            raise OSError("Incomplete HTTP artifact byte range")
        self.position += size
        return data
