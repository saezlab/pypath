"""HTTP line streaming with on-the-fly gzip and first-member zip inflate.

No download manager, no session logger: open the URL, decode as bytes arrive,
stop when the caller stops iterating (the connection is closed).
"""

from __future__ import annotations

import gzip
import io
import json
import struct
import zlib
from collections.abc import Iterator
from typing import BinaryIO

import urllib3

_USER_AGENT = "omnipath-build/0.1 (+https://omnipathdb.org; hub-export)"
_HTTP = urllib3.PoolManager(
    timeout=urllib3.Timeout(connect=60.0, read=None),
    retries=urllib3.Retry(total=5, backoff_factor=1.0, redirect=5),
    headers={"User-Agent": _USER_AGENT},
)


class _Peekable(io.RawIOBase):
    """File-like that can look at the first bytes without consuming them."""

    def __init__(self, raw: BinaryIO) -> None:
        super().__init__()
        self._raw = raw
        self._buf = bytearray()

    def peek(self, size: int) -> bytes:
        self._fill(size)
        return bytes(self._buf[:size])

    def _fill(self, size: int) -> None:
        while len(self._buf) < size:
            chunk = self._raw.read(max(size - len(self._buf), 65536))
            if not chunk:
                break
            self._buf.extend(chunk)

    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int:  # type: ignore[no-untyped-def]
        wanted = len(buffer)
        self._fill(wanted)
        out = bytes(self._buf[:wanted])
        del self._buf[:wanted]
        buffer[: len(out)] = out
        return len(out)

    def close(self) -> None:
        closer = getattr(self._raw, "close", None)
        if closer is not None:
            closer()
        super().close()


class _InflatingReader(io.RawIOBase):
    """Read a zip local-file payload, inflating deflate as chunks arrive."""

    def __init__(self, raw: BinaryIO, method: int, compressed_size: int | None) -> None:
        super().__init__()
        self._raw = raw
        self._method = method
        self._remaining = compressed_size
        self._decoder = zlib.decompressobj(-zlib.MAX_WBITS) if method == 8 else None
        self._buf = bytearray()
        self._eof = False

    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int:  # type: ignore[no-untyped-def]
        wanted = len(buffer)
        while len(self._buf) < wanted and not self._eof:
            take = 65536 if self._remaining is None else min(65536, self._remaining)
            chunk = b"" if take == 0 else self._raw.read(take)
            if self._remaining is not None:
                self._remaining -= len(chunk)
            if not chunk:
                if self._decoder is not None:
                    self._buf.extend(self._decoder.flush())
                self._eof = True
                break
            if self._decoder is not None:
                self._buf.extend(self._decoder.decompress(chunk))
            else:
                self._buf.extend(chunk)
        out = bytes(self._buf[:wanted])
        del self._buf[:wanted]
        buffer[: len(out)] = out
        return len(out)


def _open_zip_first_member(raw: BinaryIO) -> BinaryIO:
    header = raw.read(30)
    if len(header) < 30 or header[:4] != b"PK\x03\x04":
        raise ValueError("Not a zip local-file header")
    _sig, _ver, flags, method, _mt, _md, _crc, compressed_size, _us, name_len, extra_len = (
        struct.unpack("<IHHHHHIIIHH", header)
    )
    raw.read(name_len)
    raw.read(extra_len)
    if method not in {0, 8}:
        raise ValueError(f"Unsupported zip compression method {method}")
    size = None if flags & 0x8 else compressed_size
    return io.BufferedReader(_InflatingReader(raw, method, size))


def open_decoded(raw: BinaryIO) -> BinaryIO:
    """Wrap a binary stream: gzip, first zip member, or pass-through."""
    peeked = raw if isinstance(raw, _Peekable) else _Peekable(raw)
    magic = peeked.peek(4)
    if magic[:2] == b"\x1f\x8b":
        return gzip.GzipFile(fileobj=peeked)  # type: ignore[return-value]
    if magic[:4] == b"PK\x03\x04":
        return _open_zip_first_member(peeked)
    return io.BufferedReader(peeked)


def iter_decoded_lines(raw: BinaryIO, encoding: str = "utf-8") -> Iterator[str]:
    decoded = open_decoded(raw)
    text = io.TextIOWrapper(decoded, encoding=encoding, errors="replace", newline="")
    try:
        yield from text
    finally:
        text.close()


def iter_url_lines(url: str, encoding: str = "utf-8") -> Iterator[str]:
    """Stream text lines from ``url``, inflating gzip/zip as bytes arrive."""
    filename = url.split("?")[0].rstrip("/").split("/")[-1]
    if filename:
        from ..discovery import setup_pypath_cache

        cache = setup_pypath_cache()
        candidates = [cache / filename, cache / "pubchem" / filename]
        for candidate in candidates:
            if candidate.is_file() and candidate.stat().st_size > 0:
                with candidate.open("rb") as f:
                    yield from iter_decoded_lines(f, encoding=encoding)
                return

    response = _HTTP.request(
        "GET",
        url,
        preload_content=False,
        redirect=True,
        decode_content=False,
    )
    try:
        yield from iter_decoded_lines(response, encoding=encoding)
    finally:
        response.release_conn()
        response.close()


def fetch_json(url: str) -> dict:
    """Small JSON GET for paginated APIs (ChEMBL). Not used for bulk files."""
    response = _HTTP.request("GET", url, headers={**_HTTP.headers, "Accept": "application/json"})
    try:
        return json.loads(response.data.decode("utf-8"))
    finally:
        response.release_conn()
