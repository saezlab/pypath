"""Deterministic SHA-256 keys and stable hashing functions."""

from __future__ import annotations

import hashlib
import json
from typing import Any


def canonical_json(value: Any) -> str:
    """Deterministic JSON serialization."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def stable_hash(*parts: Any) -> str:
    """Deterministic SHA-256 hash for record identity."""
    digest = hashlib.sha256()
    for part in parts:
        if isinstance(part, (dict, list)):
            text = canonical_json(part)
        else:
            text = "" if part is None else str(part)
        digest.update(text.encode("utf-8"))
        digest.update(b"\x00")
    return digest.hexdigest()


def entity_key(entity_type: str, namespace: str, identifier: str) -> str:
    """Deterministic SHA-256 identity key for an entity observation."""
    return stable_hash(
        str(entity_type).strip().lower(), str(namespace).strip().lower(), str(identifier).strip()
    )


def relation_key(
    subject_key: str,
    predicate: str,
    object_key: str,
    annotations=None,
    *,
    statement_kind: str = "relation",
) -> str:
    """Identity includes all Biolink qualifiers and symmetric endpoint ordering."""
    from .biolink import statement_identity

    return relation_key_from_identity(
        statement_identity(subject_key, predicate, object_key, annotations),
        statement_kind=statement_kind,
    )


def relation_key_from_identity(identity: tuple, *, statement_kind: str = "relation") -> str:
    """Hash an already validated statement without repeating normalization."""
    if statement_kind not in {"relation", "ontology"}:
        raise ValueError(f"Unknown statement kind: {statement_kind}")
    if statement_kind == "ontology":
        return stable_hash("ontology", identity)
    subject, predicate, object, qualified = identity
    if qualified:
        return stable_hash(subject, predicate, object, qualified)
    return stable_hash(subject, predicate, object)
