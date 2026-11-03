"""Shared row contract for identifier hub Parquet files."""

from __future__ import annotations

from typing import Any

import pyarrow as pa

HUB_FIELDS = (
    "source_type",
    "source_id",
    "hub_id",
    "taxonomy_id",
    "backend",
)

HUB_SCHEMA = pa.schema([(name, pa.string()) for name in HUB_FIELDS])

CHEMICAL_TAXON = "0"


def hub_row(
    source_type: str,
    source_id: str,
    hub_id: str,
    taxonomy_id: str,
    backend: str,
) -> dict[str, str]:
    return {
        "source_type": source_type,
        "source_id": source_id,
        "hub_id": hub_id,
        "taxonomy_id": taxonomy_id,
        "backend": backend,
    }


def as_values(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [text for item in value if (text := _clean_id(item))]
    text = _clean_id(value)
    return [text] if text else []


def _clean_id(value: Any) -> str:
    text = str(value).strip()
    if not text or text.lower() in {"none", "nan", "null", "inchi=none", "inchikey=none"}:
        return ""
    return text
