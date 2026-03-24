"""Canonical ontology schema helpers for silver-layer ontology artifacts."""
from __future__ import annotations

from typing import NamedTuple

__all__ = [
    'OntologyRelationship',
    'OntologyTerm',
    'OntologyTypedef',
    'OntologyDocument',
]


class OntologyRelationship(NamedTuple):
    type: str
    target: str
    target_name: str | None = None


class OntologyTerm(NamedTuple):
    id: str
    name: str
    definition: str | None = None
    synonyms: list[str] | None = None
    comments: list[str] | None = None
    xrefs: list[str] | None = None
    is_a: list[str] | None = None
    relationships: list[OntologyRelationship] | None = None
    is_obsolete: bool | None = None


class OntologyTypedef(NamedTuple):
    id: str
    name: str


class OntologyDocument(NamedTuple):
    ontology: str
    default_namespace: str | None = None
    remark: str | None = None
    typedefs: list[OntologyTypedef] | None = None
