"""Present compiled annotations; source payload interpretation belongs to parsers."""

from __future__ import annotations

import json
from typing import Any
from omnipath_core.naming import normalize_namespace

_EMPTY_VALUES = {"", "-", "none", "null", "na", "n/a", "."}


def unwrap_attribute_value(value: Any) -> str:
    """Return scalar string from ingest wrappers like {'canonical': false, 'value': '...' }."""
    current: Any = value
    for _ in range(3):
        if current is None:
            return ""
        if isinstance(current, dict):
            if "value" in current:
                current = current.get("value")
                continue
            return json.dumps(current, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        if isinstance(current, (list, tuple)):
            return json.dumps(list(current), ensure_ascii=False, separators=(",", ":"))
        text = str(current).strip()
        if len(text) >= 2 and text[0] == "{" and '"value"' in text:
            try:
                parsed = json.loads(text)
            except json.JSONDecodeError:
                return str(current)
            if isinstance(parsed, dict) and "value" in parsed:
                current = parsed.get("value")
                continue
        return str(current)
    return "" if current is None else str(current)


def normalize_identifier_type(ns: str, value: str = "") -> str:
    """Read the namespace provided by ingestion without guessing from the value."""
    return normalize_namespace(ns)


def _annotation(term: str, value: str | None = None) -> dict[str, Any]:
    return {"term": term, "value": value}


def _is_empty(value: Any) -> bool:
    return value is None or value == ""


def split_compiled_annotations(annotations: Any) -> dict[str, list[dict[str, Any]]]:
    """Partition structured annotations list into subject / object / relation scopes."""
    subject: list[dict[str, Any]] = []
    obj: list[dict[str, Any]] = []
    relation: list[dict[str, Any]] = []
    if not isinstance(annotations, list):
        return {"subject": subject, "object": obj, "relation": relation}
    for item in annotations:
        if not isinstance(item, dict):
            continue
        term = str(item.get("term") or "").strip()
        if not term:
            continue
        value = item.get("value")
        if _is_empty(value):
            value = None
        else:
            value = str(value)
        rec = _annotation(term, value)
        if item.get("quantity") is not None:
            rec["quantity"] = item["quantity"]
        scope = str(item.get("scope") or "relation").lower()
        if scope == "subject":
            subject.append(rec)
        elif scope == "object":
            obj.append(rec)
        else:
            relation.append(rec)
    return {"subject": subject, "object": obj, "relation": relation}
