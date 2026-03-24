"""Declarative helpers to turn rows into ontology terms."""
from __future__ import annotations

from typing import Any, Sequence

from pypath.internals.ontology_schema import OntologyRelationship, OntologyTerm
from pypath.internals.tabular_builder import Column, ColumnCache

__all__ = [
    'RelationshipBuilder',
    'OntologyBuilder',
]


class RelationshipBuilder:
    def __init__(
        self,
        *,
        type: str,
        target: Column | str,
        target_name: Column | str | None = None,
    ) -> None:
        self.type = type
        self.target = target if isinstance(target, Column) else Column(target)
        self.target_name = (
            target_name if isinstance(target_name, Column) or target_name is None else Column(target_name)
        )

    def build(self, row: Any, cache: ColumnCache) -> list[OntologyRelationship]:
        targets = cache.values(self.target, row)
        target_names = cache.values(self.target_name, row) if self.target_name is not None else []
        rels: list[OntologyRelationship] = []
        for idx, target in enumerate(targets):
            if not target:
                continue
            target_name = target_names[idx] if idx < len(target_names) and target_names[idx] else None
            rels.append(OntologyRelationship(type=self.type, target=str(target), target_name=target_name))
        return rels


class OntologyBuilder:
    def __init__(
        self,
        *,
        id: Column | str,
        name: Column | str,
        definition: Column | str | None = None,
        synonyms: Column | str | None = None,
        comments: Column | str | None = None,
        xrefs: Column | str | None = None,
        is_a: Column | str | None = None,
        relationships: Sequence[RelationshipBuilder] | None = None,
        is_obsolete: Column | str | None = None,
    ) -> None:
        self.id = id if isinstance(id, Column) else Column(id)
        self.name = name if isinstance(name, Column) else Column(name)
        self.definition = definition if isinstance(definition, Column) or definition is None else Column(definition)
        self.synonyms = synonyms if isinstance(synonyms, Column) or synonyms is None else Column(synonyms)
        self.comments = comments if isinstance(comments, Column) or comments is None else Column(comments)
        self.xrefs = xrefs if isinstance(xrefs, Column) or xrefs is None else Column(xrefs)
        self.is_a = is_a if isinstance(is_a, Column) or is_a is None else Column(is_a)
        self.relationships = list(relationships or [])
        self.is_obsolete = (
            is_obsolete if isinstance(is_obsolete, Column) or is_obsolete is None else Column(is_obsolete)
        )

    def __call__(self, row: Any) -> OntologyTerm | None:
        return self.build(row)

    def build(self, row: Any) -> OntologyTerm | None:
        cache = ColumnCache()
        term_ids = cache.values(self.id, row)
        names = cache.values(self.name, row)
        term_id = term_ids[0] if term_ids else None
        name = names[0] if names else None
        if not term_id or not name:
            return None

        def _list(column: Column | None) -> list[str] | None:
            if column is None:
                return None
            values = [str(v) for v in cache.values(column, row) if v not in (None, '')]
            return values or None

        def _scalar(column: Column | None) -> str | None:
            values = _list(column)
            return values[0] if values else None

        relationships: list[OntologyRelationship] = []
        for relationship in self.relationships:
            relationships.extend(relationship.build(row, cache))

        is_obsolete = None
        if self.is_obsolete is not None:
            values = cache.values(self.is_obsolete, row)
            if values:
                value = values[0]
                if isinstance(value, bool):
                    is_obsolete = value
                else:
                    is_obsolete = str(value).strip().lower() in {'1', 'true', 'yes'}

        return OntologyTerm(
            id=str(term_id),
            name=str(name),
            definition=_scalar(self.definition),
            synonyms=_list(self.synonyms),
            comments=_list(self.comments),
            xrefs=_list(self.xrefs),
            is_a=_list(self.is_a),
            relationships=relationships or None,
            is_obsolete=is_obsolete,
        )
