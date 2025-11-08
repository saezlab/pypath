#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Declarative helpers to turn tabular rows into silver-layer entities."""

from __future__ import annotations

from dataclasses import dataclass
import logging
import re
from typing import Any, Callable, Mapping, Sequence

from pypath.internals.silver_schema import (
    Annotation as SilverAnnotation,
    Entity as SilverEntity,
    Identifier as SilverIdentifier,
    Membership as SilverMembership,
)
from pypath.internals.cv_terms import (
    AnnotationTypeCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
)

logger = logging.getLogger(__name__)


@dataclass
class ParsedValue:
    """Result of parsing a single column token."""

    value: str | None
    prefix: str | None = None
    units: str | None = None
    raw: str | None = None


class ColumnCache(dict):
    """Cache parsed column values per row to avoid duplicate work."""

    def values(self, column: 'Column', row: Any) -> list[ParsedValue]:
        if column not in self:
            self[column] = column.extract(row)
        return self[column]


class Column:
    """Column definition describing how to extract values from a row."""

    def __init__(
        self,
        selector: int | str | Callable[[Any], Any],
        *,
        delimiter: str | None = None,
        processing: dict[str, str | re.Pattern] | None = None,
        cv: Any | dict[str, Any] | Callable[[ParsedValue], Any] | None = None,
        term_cv: AnnotationTypeCv | None = None,
        unit_cv: IdentifierNamespaceCv | None = None,
    ) -> None:
        self.selector = selector
        self.delimiter = delimiter
        self.processing = processing or {}
        self.cv = cv
        self.term_cv = term_cv
        self.unit_cv = unit_cv
        self._compiled_regex: dict[str, re.Pattern] = {}

    def extract(self, row: Any) -> list[ParsedValue]:
        raw_value = self._lookup(row)
        if raw_value is None:
            return []

        tokens = self._split(raw_value)
        parsed: list[ParsedValue] = []
        for token in tokens:
            text = self._normalize_token(token)
            if text in (None, '', '-'):
                continue

            prefix = self._apply_regex('extract_prefix', text)
            if prefix is not None:
                prefix = prefix.strip().lower()

            value = self._apply_regex('extract_value', text)
            value = value.strip() if isinstance(value, str) else value
            if value in (None, ''):
                value = text

            units = self._apply_regex('extract_unit', text)
            units = units.strip() if isinstance(units, str) else units

            parsed.append(
                ParsedValue(
                    value=value if isinstance(value, str) else str(value),
                    prefix=prefix,
                    units=units if isinstance(units, str) else None,
                    raw=text,
                )
            )

        return parsed

    def resolve_cv(self, value: ParsedValue) -> Any | None:
        if self.cv is None:
            return None

        if callable(self.cv):
            try:
                return self.cv(value)
            except Exception as exc:  # pragma: no cover - defensive
                logger.debug('Column cv callable failed: %s', exc)
                return None

        if isinstance(self.cv, dict):
            lowered = (value.prefix or '').lower()
            if lowered in self.cv:
                return self.cv[lowered]
            return self.cv.get('default') or self.cv.get('DEFAULT')

        return self.cv

    def _apply_regex(self, key: str, text: str) -> str | None:
        pattern = self.processing.get(key)
        if not pattern:
            return None

        if callable(pattern):
            return pattern(text)

        if isinstance(pattern, re.Pattern):
            regex = pattern
        else:
            regex = self._compiled_regex.get(key)
            if not regex:
                regex = re.compile(pattern)
                self._compiled_regex[key] = regex

        match = regex.search(text)
        if not match:
            return None

        if match.lastindex:
            return match.group(1)
        return match.group(0)

    def _lookup(self, row: Any) -> Any:
        if callable(self.selector):
            return self.selector(row)

        if isinstance(row, Mapping):
            if self.selector in row:
                return row[self.selector]

        if isinstance(self.selector, int):
            if isinstance(row, Sequence) and not isinstance(row, (str, bytes)):
                idx = self.selector
                if -len(row) <= idx < len(row):
                    return row[idx]

        return None

    def _split(self, value: Any) -> list[Any]:
        if value is None:
            return []

        if isinstance(value, (list, tuple)):
            return list(value)

        if isinstance(value, str):
            if self.delimiter:
                return value.split(self.delimiter)
            return [value]

        text = str(value)
        if self.delimiter:
            return text.split(self.delimiter)
        return [text]

    def _normalize_token(self, token: Any) -> str | None:
        if token is None:
            return None
        if isinstance(token, str):
            return token.strip().strip('"')
        return str(token).strip()


class Identifiers:
    """Collection of `Column` objects describing identifiers."""

    def __init__(self, *columns: Column) -> None:
        self.columns = columns

    def build(self, row: Any, cache: ColumnCache | None = None) -> list[SilverIdentifier]:
        cache = cache or ColumnCache()
        identifiers: list[SilverIdentifier] = []
        for column in self.columns:
            for value in cache.values(column, row):
                id_type = column.resolve_cv(value)
                literal = value.value or value.raw
                if not id_type or not literal:
                    continue
                identifiers.append(
                    SilverIdentifier(type=id_type, value=literal)
                )
        return identifiers

    def build_for_index(
        self,
        row: Any,
        index: int,
        cache: ColumnCache,
    ) -> list[SilverIdentifier]:
        identifiers: list[SilverIdentifier] = []
        for column in self.columns:
            values = cache.values(column, row)
            if index >= len(values):
                continue
            value = values[index]
            id_type = column.resolve_cv(value)
            literal = value.value or value.raw
            if not id_type or not literal:
                continue
            identifiers.append(
                SilverIdentifier(type=id_type, value=literal)
            )
        return identifiers


class Annotations:
    """Collection of annotation columns for entities or memberships."""

    def __init__(self, *columns: Column) -> None:
        self.columns = columns

    def build(self, row: Any, cache: ColumnCache | None = None) -> list[SilverAnnotation]:
        cache = cache or ColumnCache()
        annotations: list[SilverAnnotation] = []
        for column in self.columns:
            const_term = None
            if column.term_cv is not None:
                const_term = column.term_cv
            elif column.cv is not None and not isinstance(column.cv, (dict, Callable)):
                const_term = column.cv
            for value in cache.values(column, row):
                resolved_term = const_term or column.resolve_cv(value)
                if not resolved_term:
                    continue
                annotations.append(
                    SilverAnnotation(
                        term=resolved_term,
                        value=value.value or value.raw,
                        units=value.units,
                    )
                )
        return annotations

    def build_for_index(
        self,
        row: Any,
        index: int,
        cache: ColumnCache,
    ) -> list[SilverAnnotation]:
        annotations: list[SilverAnnotation] = []
        for column in self.columns:
            const_term = None
            if column.term_cv is not None:
                const_term = column.term_cv
            elif column.cv is not None and not isinstance(column.cv, (dict, Callable)):
                const_term = column.cv
            values = cache.values(column, row)
            if index >= len(values):
                continue
            value = values[index]
            resolved_term = const_term or column.resolve_cv(value)
            if not resolved_term:
                continue
            annotations.append(
                SilverAnnotation(
                    term=resolved_term,
                    value=value.value or value.raw,
                    units=value.units,
                )
            )
        return annotations


class Entities:
    """Definition of member entities derived from aligned column splits."""

    def __init__(
        self,
        *,
        entity_type: EntityTypeCv,
        entity_source: str,
        identifiers: Identifiers,
        annotations: Annotations | None = None,
        entity_annotations: Annotations | None = None,
    ) -> None:
        self.entity_type = entity_type
        self.entity_source = entity_source
        self.identifiers = identifiers
        self.membership_annotations = annotations
        self.entity_annotations = entity_annotations

    def build(self, row: Any, cache: ColumnCache) -> list[SilverMembership]:
        lengths = [len(cache.values(column, row)) for column in self.identifiers.columns]

        if self.membership_annotations:
            lengths.extend(len(cache.values(column, row)) for column in self.membership_annotations.columns)

        if self.entity_annotations:
            lengths.extend(len(cache.values(column, row)) for column in self.entity_annotations.columns)

        member_count = max(lengths) if lengths else 0
        memberships: list[SilverMembership] = []

        for index in range(member_count):
            member_identifiers = self.identifiers.build_for_index(row, index, cache)
            if not member_identifiers:
                continue

            entity_annotations = (
                self.entity_annotations.build_for_index(row, index, cache)
                if self.entity_annotations
                else None
            )
            membership_annotations = (
                self.membership_annotations.build_for_index(row, index, cache)
                if self.membership_annotations
                else None
            )

            member_entity = SilverEntity(
                source=self.entity_source,
                type=self.entity_type,
                identifiers=member_identifiers,
                annotations=entity_annotations if entity_annotations else None,
                members=None,
                is_member_of=None,
            )

            memberships.append(
                SilverMembership(
                    member=member_entity,
                    annotations=membership_annotations if membership_annotations else None,
                )
            )

        return memberships


class Members:
    """Container for one or more `Entities` definitions."""

    def __init__(self, *entities: Entities) -> None:
        self.entities = entities

    def build(self, row: Any, cache: ColumnCache) -> list[SilverMembership]:
        memberships: list[SilverMembership] = []
        for entity in self.entities:
            memberships.extend(entity.build(row, cache))
        return memberships


class Entity:
    """Declarative spec that produces `SilverEntity` records from rows."""

    def __init__(
        self,
        *,
        source: str,
        entity_type: EntityTypeCv,
        identifiers: Identifiers,
        annotations: Annotations | None = None,
        members: Members | None = None,
    ) -> None:
        self.source = source
        self.entity_type = entity_type
        self.identifiers = identifiers
        self.annotations = annotations
        self.members = members

    def __call__(self, row: Any) -> SilverEntity | None:
        return self.build(row)

    def build(self, row: Any) -> SilverEntity | None:
        cache = ColumnCache()
        identifiers = self.identifiers.build(row, cache)
        if not identifiers:
            return None

        annotations = self.annotations.build(row, cache) if self.annotations else None
        members = self.members.build(row, cache) if self.members else None

        return SilverEntity(
            source=self.source,
            type=self.entity_type,
            identifiers=identifiers,
            annotations=annotations if annotations else None,
            members=members if members else None,
            is_member_of=None,
        )


__all__ = [
    'Annotations',
    'Column',
    'Entity',
    'Entities',
    'Identifiers',
    'Members',
]
