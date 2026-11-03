"""Stable input contract accepted by EntityResolver and LibraryMatcher.

Extraction may attach identifiers, annotations and structural derivations. Runtime
matching consumes this record without importing source parsers or build tooling.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any


@dataclass(slots=True)
class RawEntityObservation:
    """Extracted raw entity observation before canonicalization."""

    entity_key: str
    entity_type: str
    namespace: str
    identifier: str
    taxon: str | None = None
    label: str | None = None
    identifiers: list[dict[str, Any]] = field(default_factory=list)
    annotations: list[dict[str, Any]] = field(default_factory=list)
    identity_scope: str | None = None
    structure_derivations: list[dict[str, Any]] = field(default_factory=list)
    molecular_form: dict[str, Any] | None = None
    evidence: list[dict[str, Any]] = field(default_factory=list)
