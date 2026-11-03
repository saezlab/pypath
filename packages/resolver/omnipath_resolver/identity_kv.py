"""LMDB point-lookup stores of the identity layer (spec 3a/3b), read by ``IdentityRuntime``.

Hub index ``<hub dir>/kv/`` (written by ``omnipath-build build-hub-kv``)::

    manifest.json   format, hub, source index sha, shard counts, per-env statistics
    id.<s>/         db ``id``:  part + ns + b"\\0" + identifier -> "local_id<F>tag<R>local_id<F>tag..."
                    (symbol namespaces, format v3: part + ns + b"\\0" + identifier + b"\\0" + taxon)
    rec.<s>/        db ``rec``: part + local_id -> "taxon<F>anchor<F>count<F>reviewed[<R>type<F>value...]"

``part`` is the hub index partition (two hex digits of md5 of the identifier, or of
``<hub>:<local_id>``) and the shard is its first digit (shard 0 when a hub has one shard). The
prefix makes every partition's keys a contiguous ascending block, so the writer appends partition
files as they are, with no global sort. <F> is ``FIELD`` (0x1f) and <R> is ``ROW`` (0x1e); values
are UTF-8 and parsed with ``str.split``. Empty taxon/anchor mean none; reviewed is 1 or 0.

A gene symbol is the same word in hundreds of species, and a symbol key is always taxon-scoped.
From format v3 the ``id`` keys of the symbol namespaces end with the taxon of the record
(empty for a record without one), so a scoped lookup reads the rows of its taxon and those of
the records without a taxon, and never the head of another species' record. Format v2 keys have
no taxon suffix; the runtime reads both.

Identity directory ``<identity dir>/kv/`` (written by ``omnipath-build build-identity-kv``), one
environment with the dbs ``exc`` (record_id -> [entity_id, decision, quarantined]),
``exc_members`` (entity_id -> [record_id]), ``extra`` (entity_id -> [kind, taxon, quarantined,
preferred_record]), ``gp_protein`` (protein entity -> [[entrez_id, taxon]]), ``gp_gene`` (entrez
id -> [[protein entity, taxon]]), ``lipid`` (goslin name -> inchikey) and ``cand`` (record_id
-> [entity_id]). These small stores keep msgpack values with a one-byte prefix: ``\\x00`` raw,
``\\x01`` zstd (above ``COMPRESS_OVER`` bytes).
"""

from __future__ import annotations

import hashlib
import json
import threading
from pathlib import Path

import lmdb
import msgpack
import zstandard

KV_DIR = "kv"
MANIFEST = "manifest.json"
HUB_KV_FORMAT = "omnipath-hub-kv-v3"
HUB_KV_FORMATS = ("omnipath-hub-kv-v2", HUB_KV_FORMAT)  # the formats HubKv can read
FIELD, ROW = "\x1f", "\x1e"
DECISIONS_KV_FORMAT = "omnipath-identity-kv-v1"
MAX_KEY = 511  # LMDB's default key size limit
COMPRESS_OVER = 256
RAW, ZSTD = b"\x00", b"\x01"
SHARD_COUNTS = (1, 16)
# Namespaces whose lookup rows are keyed by the taxon of their record (format v3).
SYMBOL_NAMESPACES = frozenset({"genesymbol", "genesymbol-syn"})
DECISION_DBS = ("exc", "exc_members", "extra", "gp_protein", "gp_gene", "lipid", "cand")


class KvMissing(RuntimeError):
    """A key/value store the runtime needs is absent or does not match its source."""


# ---------------------------------------------------------------------------- codec

_local = threading.local()


def _compressor():
    if not hasattr(_local, "c"):
        _local.c = zstandard.ZstdCompressor(level=3)
        _local.d = zstandard.ZstdDecompressor()
    return _local.c


def pack(obj) -> bytes:
    body = msgpack.packb(obj, use_bin_type=True)
    if len(body) > COMPRESS_OVER:
        return ZSTD + _compressor().compress(body)
    return RAW + body


def body_of(raw: bytes) -> bytes:
    prefix = raw[:1]
    if prefix == RAW:
        return raw[1:]
    if prefix == ZSTD:
        _compressor()
        return _local.d.decompress(raw[1:])
    raise ValueError("Unsupported kv value prefix")


def unpack(raw: bytes):
    return msgpack.unpackb(body_of(raw), raw=False)


# ----------------------------------------------------------------------- key helpers


def partition(value: str) -> str:
    return hashlib.md5(value.encode(), usedforsecurity=False).hexdigest()[:2]


def shard_of(value: str, shards: int) -> int:
    """Shard of ``value`` (an identifier, or a record id for ``rec``) among ``shards``."""
    if shards == 1:
        return 0
    if shards != 16:
        raise ValueError(f"Unsupported shard count {shards}")
    return int(partition(value)[0], 16)


def id_key(ns: str, identifier: str, taxon: str | None = None) -> bytes | None:
    """Key of the ``id`` db, or None when it cannot be stored (longer than LMDB's key limit).

    ``taxon`` is for the symbol namespaces of format v3 ('' is the records without a taxon).
    """
    key = partition(identifier).encode() + ns.encode() + b"\0" + identifier.encode()
    if taxon is not None:
        key += b"\0" + taxon.encode()
    return key if len(key) <= MAX_KEY else None


def rec_key(hub: str, local_id: str) -> bytes | None:
    if not local_id:
        return None
    key = partition(f"{hub}:{local_id}").encode() + local_id.encode()
    return key if len(key) <= MAX_KEY else None


def decode_ids(raw: bytes) -> list[tuple[str, str]]:
    return [tuple(row.split(FIELD, 1)) for row in raw.decode().split(ROW)]


def decode_head(raw: bytes) -> tuple:
    taxon, anchor, count, reviewed = raw.split(ROW.encode(), 1)[0].decode().split(FIELD)
    return (taxon or None, anchor or None, int(count or 0), reviewed == "1")


def decode_rows(raw: bytes) -> list[tuple[str, str]]:
    return [tuple(row.split(FIELD, 1)) for row in raw.decode().split(ROW)[1:]]


def sha256_file(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ----------------------------------------------------------------------------- reading


class Env:
    """One read-only environment: a single long-lived read transaction, named dbs opened once."""

    def __init__(self, path: Path, names):
        path = Path(path)
        if not (path / "data.mdb").is_file():
            raise KvMissing(f"Missing LMDB store {path}")
        self.path = path
        self.env = lmdb.open(
            str(path),
            readonly=True,
            lock=False,
            readahead=False,
            max_dbs=max(len(names), 1),
            subdir=True,
        )
        self.txn = self.env.begin()
        try:
            self.dbs = {name: self.env.open_db(name.encode(), txn=self.txn, create=False) for name in names}
        except lmdb.Error as exc:
            self.close()
            raise KvMissing(f"LMDB store {path} lacks a database: {exc}") from exc

    def get_many(self, name: str, keys) -> dict:
        """key -> raw value for the keys present (read in key order for page locality)."""
        db, txn, found = self.dbs[name], self.txn, {}
        for key in sorted(set(keys)):
            raw = txn.get(key, db=db)
            if raw is not None:
                found[key] = raw
        return found

    def items(self, name: str):
        yield from self.txn.cursor(db=self.dbs[name])

    def close(self):
        txn, env = getattr(self, "txn", None), getattr(self, "env", None)
        self.txn = self.env = None
        if txn is not None:
            txn.abort()
        if env is not None:
            env.close()


def _read_manifest(directory: Path, formats, what: str, build_hint: str) -> dict:
    path = directory / MANIFEST
    if not path.is_file():
        raise KvMissing(f"{what}: no key/value store at {directory} ({build_hint})")
    manifest = json.loads(path.read_text())
    formats = (formats,) if isinstance(formats, str) else tuple(formats)
    if manifest.get("format") not in formats:
        raise KvMissing(
            f"{what}: {path} has format {manifest.get('format')!r}, expected {' or '.join(formats)!r}"
        )
    return manifest


class HubKv:
    """``id`` and ``rec`` lookups of one hub index (environments are opened on first use)."""

    def __init__(self, hub: str, index_dir: Path):
        self.hub = hub
        self.index_dir = Path(index_dir)
        directory = self.index_dir / KV_DIR
        hint = f"run: omnipath-build build-hub-kv --hub {hub} --hub-index-root <root>"
        self.manifest = _read_manifest(directory, HUB_KV_FORMATS, f"hub index {hub}", hint)
        # Format v3: the rows of the symbol namespaces are keyed by taxon (``taxon_ids``).
        self.by_taxon = self.manifest["format"] == HUB_KV_FORMAT
        source = self.index_dir / "manifest.json"
        if source.is_file() and self.manifest.get("source_sha256") != sha256_file(source):
            raise KvMissing(f"hub index {hub}: kv at {directory} was built from a different index")
        self.id_shards = int(self.manifest["id_shards"])
        self.rec_shards = int(self.manifest["rec_shards"])
        for kind, count in (("id", self.id_shards), ("rec", self.rec_shards)):
            for shard in range(count):
                if not (directory / f"{kind}.{shard}" / "data.mdb").is_file():
                    raise KvMissing(f"hub index {hub}: missing {kind}.{shard} under {directory}")
        self.directory = directory
        self._envs: dict[tuple[str, int], Env] = {}
        self._lock = threading.RLock()

    def _env(self, kind: str, shard: int) -> Env:
        key = (kind, shard)
        env = self._envs.get(key)
        if env is None:
            env = self._envs[key] = Env(self.directory / f"{kind}.{shard}", (kind,))
        return env

    def opened(self) -> set[tuple[str, int]]:
        return set(self._envs)

    def ids(self, pairs) -> dict[tuple[str, str], list[tuple[str, str]]]:
        """(ns, identifier) -> [(local_id, tag)] for the pairs present.

        A symbol namespace is answered for every taxon by format v2 only; format v3 reads
        it with ``taxon_ids``.
        """
        by_shard: dict[int, dict[bytes, tuple[str, str]]] = {}
        for ns, identifier in pairs:
            if self.by_taxon and ns in SYMBOL_NAMESPACES:
                raise ValueError(f"Namespace {ns} is keyed by taxon in kv format v3")
            key = id_key(ns, identifier)
            if key is not None:
                by_shard.setdefault(shard_of(identifier, self.id_shards), {})[key] = (ns, identifier)
        found = {}
        with self._lock:
            for shard, keys in by_shard.items():
                for key, raw in self._env("id", shard).get_many("id", keys).items():
                    found[keys[key]] = decode_ids(raw)
        return found

    def taxon_ids(self, wanted) -> dict[tuple[str, str, str], list[tuple[str, str]]]:
        """(ns, identifier, taxon) -> [(local_id, tag)] of the symbol rows (format v3 only).

        The rows of the records of that taxon and those of the records without a taxon, which
        belong to every scope; the others are not read.
        """
        if not self.by_taxon:
            raise ValueError(
                f"hub index {self.hub}: kv format {self.manifest['format']} has no taxon keys"
            )
        # A key without a taxon is shared by every wanted taxon of its symbol.
        by_shard: dict[int, dict[bytes, list[tuple[str, str, str]]]] = {}
        for ns, identifier, taxon in set(wanted):
            if ns not in SYMBOL_NAMESPACES:
                raise ValueError(f"Namespace {ns} is not keyed by taxon")
            for suffix in {taxon, ""}:
                key = id_key(ns, identifier, suffix)
                if key is not None:
                    shard = shard_of(identifier, self.id_shards)
                    keys = by_shard.setdefault(shard, {})
                    keys.setdefault(key, []).append((ns, identifier, taxon))
        found: dict[tuple[str, str, str], list[tuple[str, str]]] = {}
        with self._lock:
            for shard, keys in by_shard.items():
                for key, raw in self._env("id", shard).get_many("id", keys).items():
                    for item in keys[key]:
                        found.setdefault(item, []).extend(decode_ids(raw))
        return found

    def _rec(self, local_ids, decode) -> dict:
        by_shard: dict[int, dict[bytes, str]] = {}
        for local_id in local_ids:
            key = rec_key(self.hub, local_id)
            if key is not None:
                shard = shard_of(f"{self.hub}:{local_id}", self.rec_shards)
                by_shard.setdefault(shard, {})[key] = local_id
        found = {}
        with self._lock:
            for shard, keys in by_shard.items():
                for key, raw in self._env("rec", shard).get_many("rec", keys).items():
                    found[keys[key]] = decode(raw)
        return found

    def heads(self, local_ids) -> dict[str, tuple]:
        """local_id -> (taxon, anchor, anchor_count, reviewed)."""
        return self._rec(local_ids, decode_head)

    def rows(self, local_ids) -> dict[str, list[tuple[str, str]]]:
        """local_id -> [(source_type, value)] (every by_record row of the record)."""
        return self._rec(local_ids, decode_rows)

    def close(self):
        with self._lock:
            for env in self._envs.values():
                env.close()
            self._envs.clear()


class DecisionsKv:
    """The identity decisions' small stores (one environment, six dbs)."""

    def __init__(self, identity_dir: Path, fingerprint: str | None = None):
        self.directory = Path(identity_dir) / KV_DIR
        hint = f"run: omnipath-build build-identity-kv --identity-dir {identity_dir}"
        self.manifest = _read_manifest(self.directory, DECISIONS_KV_FORMAT, "identity decisions", hint)
        if fingerprint and self.manifest.get("fingerprint") != fingerprint:
            raise KvMissing(f"identity decisions: kv at {self.directory} belongs to another fingerprint")
        self.env = Env(self.directory, DECISION_DBS)
        self._lock = threading.RLock()

    def get(self, name: str, keys) -> dict[str, object]:
        """key -> decoded value for the string keys present."""
        encoded = {}
        for key in keys:
            raw = key.encode()
            if len(raw) <= MAX_KEY:
                encoded[raw] = key
        with self._lock:
            return {encoded[k]: unpack(v) for k, v in self.env.get_many(name, encoded).items()}

    def items(self, name: str):
        with self._lock:
            return [(k.decode(), unpack(v)) for k, v in self.env.items(name)]

    def close(self):
        with self._lock:
            self.env.close()
