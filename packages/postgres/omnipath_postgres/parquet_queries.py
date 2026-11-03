"""Read and validate published Parquet tables without expanding Python rows."""

from __future__ import annotations
from pathlib import Path
from typing import Any
from .locations import is_remote, join_location

# Match Python str.strip(), including Unicode spaces and the four control
# separators Python treats as whitespace. Identifier text itself is never trimmed.
_KEY_WHITESPACE = " \t\n\r\v\f\x1c\x1d\x1e\x1f\x85\xa0                　"


def _literal(value: str) -> str:
    if not isinstance(value, str) or "\x00" in value:
        raise ValueError("SQL literal must be a string without NUL characters")
    return "'" + value.replace("'", "''") + "'"


def _scan(directory: str | Path, filename: str) -> str:
    location = join_location(directory, filename)
    path = location if is_remote(location) else str(location.absolute())
    return f"read_parquet({_literal(path)}, hive_partitioning=false)"


def quantity_sql() -> str:
    """An annotation table row's ``quantity_*`` columns as one struct, NULL when all
    are null (the published tables store a quantity as columns)."""
    from omnipath_core.measurements import QUANTITY_STRUCT

    names = [field.name for field in QUANTITY_STRUCT]
    empty = " AND ".join(f"quantity_{n} IS NULL" for n in names)
    fields = ", ".join(f"{n} := quantity_{n}" for n in names)
    return f"CASE WHEN {empty} THEN NULL ELSE struct_pack({fields}) END"


def _bad_key(expression: str) -> str:
    return f"({expression} IS NULL OR trim({expression},{_literal(_KEY_WHITESPACE)})='')"


def _any(con: Any, table: str, condition: str) -> bool:
    return bool(con.execute(f"SELECT count(*) FROM {table} WHERE {condition}").fetchone()[0])


def validate_inputs(con: Any) -> None:
    """Reject staged inputs the projection cannot represent.

    Keys must be nonempty; evidence annotation lists (the only nested lists left)
    must not contain null entries, which would crash the reference projector; and
    each present numeric quantity must be finite before SQL JSON serialization.
    """
    if _any(con, "ap_entity_raw", _bad_key("entity_key")):
        raise ValueError("entity_key must be a nonempty string")
    keys = " OR ".join(
        _bad_key(name) for name in ("relation_key", "subject_entity_key", "object_entity_key")
    )
    if _any(con, "ap_statement_raw", keys):
        raise ValueError(
            "relation_key, subject_entity_key and object_entity_key must be nonempty strings"
        )
    null_entries = "list_contains(list_transform(item.annotations, a -> a IS NULL), true)"
    if _any(con, "ap_evidence_raw", null_entries) or _any(
        con, "ap_entity_evidence_raw", null_entries
    ):
        raise ValueError("Null nested struct in evidence annotations")
    if _any(
        con,
        "ap_annotation_raw",
        "quantity.has_numeric_value IS NOT NULL AND NOT isfinite(quantity.has_numeric_value)",
    ):
        raise ValueError("Invalid JSON value in annotations: quantity must be finite")
