"""Small dictionary Parquet files shipped next to the hub maps."""

from __future__ import annotations

import json
import tempfile
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import yaml

import omnipath_core

_VOCAB_DATA = Path(omnipath_core.__file__).resolve().parent / "vocab"


def _data_file(name: str) -> Path:
    if (_VOCAB_DATA / name).is_file():
        return _VOCAB_DATA / name
    raise FileNotFoundError(f"Bundled core vocabulary is missing: {_VOCAB_DATA / name}")


def _load_yaml(name: str) -> dict[str, Any]:
    path = _data_file(name)
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def load_organisms() -> dict[str, Any]:
    return _load_yaml("organisms.yaml")


def _publish_table(table: pa.Table, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{path.name}.", dir=path.parent) as temporary:
        staged = Path(temporary) / path.name
        pq.write_table(table, staged, compression="zstd")
        if pq.read_metadata(staged).num_rows != table.num_rows:
            raise ValueError(f"Incomplete dictionary export: {path}")
        staged.replace(path)


def write_id_type(path: Path) -> int:
    rows = []
    for name, info in _load_yaml("id_types.yaml").items():
        if not isinstance(info, dict):
            continue
        rows.append(
            {
                "name": str(name),
                "label": str(info.get("label") or name),
                "entity_type": str(info.get("entity_type") or ""),
                "curie_prefix": str(info.get("curie_prefix") or ""),
            }
        )
    table = pa.Table.from_pylist(rows)
    _publish_table(table, path)
    return table.num_rows


def write_backend(path: Path, backends: Iterable[str] | None = None) -> int:
    names = sorted(set(backends or []))
    if not names:
        seen: set[str] = set()
        for info in _load_yaml("id_types.yaml").values():
            if isinstance(info, dict):
                seen.update(str(key) for key in (info.get("backends") or {}))
        names = sorted(seen)
    table = pa.Table.from_pylist([{"name": name} for name in names])
    _publish_table(table, path)
    return table.num_rows


def write_organism(path: Path) -> int:
    rows = []
    for tax_id, info in load_organisms().items():
        if not isinstance(info, dict):
            continue
        rows.append(
            {
                "ncbi_tax_id": str(tax_id),
                "common_name": str(info.get("common_name") or ""),
                "latin_name": str(info.get("latin_name") or ""),
            }
        )
    table = pa.Table.from_pylist(rows)
    _publish_table(table, path)
    return table.num_rows


def write_manifest(path: Path, files: dict[str, int], export_id: str | None = None) -> None:
    now = datetime.now(timezone.utc).isoformat()
    payload = {
        "export_id": export_id or now.replace(":", "").replace("-", ""),
        "created_at": now,
        "files": files,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".manifest-", dir=path.parent) as temporary:
        staged = Path(temporary) / path.name
        staged.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        staged.replace(path)
