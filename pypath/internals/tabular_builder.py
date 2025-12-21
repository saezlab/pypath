#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Declarative helpers to turn tabular rows into silver-layer entities.

New DSL based on CV / FieldConfig sources (producing Columns), replacing the
older Column+cv/term_cv setup. This module intentionally breaks the previous
configuration API.
"""

from __future__ import annotations

from dataclasses import dataclass, field
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
    EntityTypeCv,
    CvEnum,
)

logger = logging.getLogger(__name__)


class Column:
    """Column definition describing how to extract and normalize values from a row.

    Parameters
    ----------
    selector:
        - ``int``      → index into a sequence row
        - ``str``      → key into a mapping row
        - ``callable`` → function ``row -> value``
    delimiter:
        Optional delimiter used to split string-valued cells into multiple
        tokens. If not provided, the entire cell is treated as a single value.
    extract:
        Optional sequence of regex patterns or callables applied in order to
        the raw string value *before* mapping. Elements may be:
        - ``str`` or ``re.Pattern`` → the first group (if any) or full match.
        - ``callable`` ``str -> str | None``.
        If any extract step returns ``None`` or fails to match, that value is
        discarded.
    map:
        - ``dict`` mapping raw/processed values to desired outputs.
          Matching is done first on the original processed value, then on
          ``lower()`` if the dict contains lowercase keys.
        - ``callable`` ``value -> mapped_value``.
        - ``None`` to return the extracted value unchanged.
    default:
        Optional fallback value if mapping yields no result.
    """

    def __init__(
        self,
        selector: int | str | Callable[[Any], Any],
        *,
        delimiter: str | None = None,
        extract: Sequence[str | re.Pattern | Callable[[str], Any]] | None = None,
        map: Mapping[Any, Any] | Callable[[Any], Any] | None = None,
        default: Any | None = None,
    ) -> None:
        self.selector = selector
        self.delimiter = delimiter
        self.extract_steps: list[str | re.Pattern | Callable[[str], Any]] = list(extract or [])
        self.mapping = map
        self.default = default

    def extract(self, row: Any, cache: ColumnCache | None = None) -> list[Any]:
        """Extract a list of processed values from the given row."""
        raw_value = self._lookup(row)
        if raw_value is None:
            return []

        tokens = self._split(raw_value)
        out: list[Any] = []
        for token in tokens:
            text = self._normalize_token(token)
            if text in (None, "", "-"):
                continue
            processed: str | Any | None = text
            for step in self.extract_steps:
                processed = self._apply_step(step, processed)
                if processed is None:
                    break

            if processed is None:
                continue

            mapped = self._apply_mapping(processed)
            if mapped is None:
                if self.default is not None:
                    mapped = self.default
                else:
                    continue

            out.append(mapped)
        return out

    # -- internal helpers -------------------------------------------------

    def _lookup(self, row: Any) -> Any:
        if callable(self.selector):
            try:
                return self.selector(row)
            except Exception as exc:  # pragma: no cover - defensive
                logger.debug("Column selector callable failed: %s", exc)
                return None

        if isinstance(row, Mapping) and isinstance(self.selector, str):
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

    def _apply_step(
        self,
        step: str | re.Pattern | Callable[[str], Any],
        value: Any,
    ) -> Any | None:
        if value is None:
            return None

        text = str(value)

        if callable(step):
            try:
                return step(text)
            except Exception as exc:  # pragma: no cover - defensive
                logger.debug("Column extract callable failed: %s", exc)
                return None

        if isinstance(step, re.Pattern):
            regex = step
        else:
            regex = re.compile(step)

        match = regex.search(text)
        if not match:
            return None

        if match.lastindex:
            return match.group(1)
        return match.group(0)

    def _apply_mapping(self, value: Any) -> Any | None:
        if self.mapping is None:
            return value

        if callable(self.mapping):
            try:
                return self.mapping(value)
            except Exception as exc:  # pragma: no cover - defensive
                logger.debug("Column mapping callable failed: %s", exc)
                return None

        try:
            if value in self.mapping:
                return self.mapping[value]
        except TypeError:
            value_str = str(value)
            if value_str in self.mapping:
                return self.mapping[value_str]
            lower = value_str.lower()
            if lower in self.mapping:
                return self.mapping[lower]
            return None

        value_str = str(value)
        if value_str in self.mapping:
            return self.mapping[value_str]
        lower = value_str.lower()
        if lower in self.mapping:
            return self.mapping[lower]

        return None


class ColumnCache(dict):
    """Cache extracted values per source to avoid duplicate work."""

    def values(self, source: Any, row: Any) -> list[Any]:
        if source not in self:
            try:
                extracted = source.extract(row, self)
            except TypeError:
                # Fallback for sources that do not accept cache (for robustness)
                extracted = source.extract(row)
            self[source] = extracted
        return self[source]


class _ConstantSource:
    """Internal source representing a constant value."""

    __slots__ = ("value",)

    def __init__(self, value: Any) -> None:
        self.value = value

    def extract(self, row: Any, cache: ColumnCache | None = None) -> list[Any]:  # noqa: ARG002
        return [self.value]


class _CallableSource:
    """Internal source wrapping a callable row -> value(s)."""

    __slots__ = ("func",)

    def __init__(self, func: Callable[[Any], Any]) -> None:
        self.func = func

    def extract(self, row: Any, cache: ColumnCache | None = None) -> list[Any]:  # noqa: ARG002
        try:
            result = self.func(row)
        except Exception as exc:  # pragma: no cover - defensive
            logger.debug("Callable source failed: %s", exc)
            return []

        if result is None:
            return []

        if isinstance(result, (list, tuple)):
            return [v for v in result if v is not None]

        return [result]


@dataclass
class FieldConfig:
    extract: dict[str, Any] = field(default_factory=dict)
    map: dict[str, Mapping[Any, Any] | Callable[[Any], Any]] = field(default_factory=dict)
    delimiter: str | None = None

    def __call__(
        self,
        selector: int | str | Callable[[Any], Any],
        *,
        extract: str | Sequence[Any] | Callable[[str], Any] | re.Pattern | None = None,
        map: str | Mapping[Any, Any] | Callable[[Any], Any] | None = None,
        delimiter: str | None = None,
        default: Any | None = None,
    ) -> Column:
        extract_steps = self._resolve_extract(extract)
        mapping = self._resolve_mapping(map)
        return Column(
            selector,
            delimiter=delimiter if delimiter is not None else self.delimiter,
            extract=extract_steps,
            map=mapping,
            default=default,
        )

    def _resolve_extract(
        self,
        extract: str | Sequence[Any] | Callable[[str], Any] | re.Pattern | None,
    ) -> list[Any] | None:
        if extract is None:
            return None
        if isinstance(extract, str):
            extract = [extract]
        if isinstance(extract, (list, tuple)):
            steps: list[Any] = []
            for item in extract:
                if isinstance(item, str) and item in self.extract:
                    value = self.extract[item]
                    if isinstance(value, (list, tuple)):
                        steps.extend(value)
                    else:
                        steps.append(value)
                else:
                    steps.append(item)
            return steps
        return [extract]

    def _resolve_mapping(
        self,
        mapping: str | Mapping[Any, Any] | Callable[[Any], Any] | None,
    ) -> Mapping[Any, Any] | Callable[[Any], Any] | None:
        if mapping is None:
            return None
        if isinstance(mapping, str):
            return self.map[mapping]
        return mapping


class CV:
    """Declarative specification of a single CV-based field.

    This is used by both :class:`IdentifiersBuilder` and :class:`AnnotationsBuilder`.

    Parameters
    ----------
    term:
        Source for the CV term / identifier type. Required.
        Can be:
        - Constant (e.g. ``MoleculeAnnotationsCv.ENDOGENOUS``)
        - :class:`Column` (typically via :class:`FieldConfig`)
        - ``callable`` ``row -> value | list[value]``
    value:
        Optional source for the value part.
    unit:
        Optional source for the unit part (only used for annotations).
    """

    def __init__(
        self,
        *,
        term: Any,
        value: Any | None = None,
        unit: Any | None = None,
    ) -> None:
        self.term_source = self._normalize_source(term)
        self.value_source = self._normalize_source(value) if value is not None else None
        self.unit_source = self._normalize_source(unit) if unit is not None else None

    # -- source normalization ---------------------------------------------

    def _normalize_source(self, spec: Any) -> Any:
        # Already a source-like object
        if isinstance(spec, (Column, _ConstantSource, _CallableSource)):
            return spec

        # Callable row -> value(s)
        if callable(spec):
            return _CallableSource(spec)

        # Fallback: constant
        return _ConstantSource(spec)


class _BaseCvBuilder:
    """Shared implementation for identifier and annotation builders."""

    def __init__(self, silver_cls: type[SilverIdentifier] | type[SilverAnnotation], *cvs: CV) -> None:
        self.silver_cls = silver_cls
        self.cvs = cvs

    # -- public API -------------------------------------------------------

    def build(self, row: Any, cache: ColumnCache | None = None) -> list[SilverIdentifier] | list[SilverAnnotation]:
        cache = cache or ColumnCache()
        results: list[SilverIdentifier] | list[SilverAnnotation] = []
        seen: set[tuple[Any, Any, Any]] = set()

        for cv in self.cvs:
            for term, value, unit in self._expand_cv(cv, row, cache):
                if term is None:
                    continue

                key = (term, value, unit if self.silver_cls is SilverAnnotation else None)
                if key in seen:
                    continue
                seen.add(key)

                if self.silver_cls is SilverIdentifier:
                    # Identifiers ignore unit and require a non-empty value
                    if value is None or value == "":
                        continue
                    results.append(SilverIdentifier(type=term, value=value))
                else:
                    results.append(SilverAnnotation(term=term, value=value, units=unit))

        return results

    def build_for_index(
        self,
        row: Any,
        index: int,
        cache: ColumnCache,
    ) -> list[SilverIdentifier] | list[SilverAnnotation]:
        results: list[SilverIdentifier] | list[SilverAnnotation] = []
        seen: set[tuple[Any, Any, Any]] = set()

        for cv in self.cvs:
            term_vals = self._safe_extract(cv.term_source, row, cache)
            if not term_vals:
                continue

            value_vals: list[Any] | None = None
            unit_vals: list[Any] | None = None

            if cv.value_source is not None:
                vv = self._safe_extract(cv.value_source, row, cache)
                value_vals = vv if vv else None

            if cv.unit_source is not None:
                uv = self._safe_extract(cv.unit_source, row, cache)
                unit_vals = uv if uv else None

            term = self._pick_index(term_vals, index)
            if term is None:
                continue
            self._validate_term(term)

            value = self._pick_index(value_vals, index) if value_vals is not None else None
            if cv.value_source is not None and not self._has_value(value):
                continue

            unit = self._pick_index(unit_vals, index) if unit_vals is not None else None
            if not self._has_value(value):
                unit = None

            key = (term, value, unit if self.silver_cls is SilverAnnotation else None)
            if key in seen:
                continue
            seen.add(key)

            if self.silver_cls is SilverIdentifier:
                if value is None or value == "":
                    continue
                results.append(SilverIdentifier(type=term, value=value))
            else:
                results.append(SilverAnnotation(term=term, value=value, units=unit))

        return results

    def max_length(self, row: Any, cache: ColumnCache) -> int:
        """Maximum sequence length implied by the CVs for this row."""
        max_len = 0
        for cv in self.cvs:
            term_vals = self._safe_extract(cv.term_source, row, cache)
            if not term_vals:
                continue
            lengths = [len(term_vals)]

            if cv.value_source is not None:
                vv = self._safe_extract(cv.value_source, row, cache)
                if vv:
                    lengths.append(len(vv))

            if cv.unit_source is not None:
                uv = self._safe_extract(cv.unit_source, row, cache)
                if uv:
                    lengths.append(len(uv))

            if lengths:
                max_len = max(max_len, max(lengths))

        return max_len

    # -- helpers ----------------------------------------------------------

    def _expand_cv(
        self,
        cv: CV,
        row: Any,
        cache: ColumnCache,
    ) -> list[tuple[Any, Any, Any]]:
        term_vals = self._safe_extract(cv.term_source, row, cache)
        if not term_vals:
            return []

        value_vals: list[Any] = [None]
        unit_vals: list[Any] = [None]

        if cv.value_source is not None:
            vv = self._safe_extract(cv.value_source, row, cache)
            if vv:
                value_vals = vv

        if cv.unit_source is not None:
            uv = self._safe_extract(cv.unit_source, row, cache)
            if uv:
                unit_vals = uv

        max_len = max(len(term_vals), len(value_vals), len(unit_vals))
        out: list[tuple[Any, Any, Any]] = []

        for i in range(max_len):
            term = self._pick_index(term_vals, i)
            if term is None:
                continue
            self._validate_term(term)

            value = self._pick_index(value_vals, i)
            if cv.value_source is not None and not self._has_value(value):
                continue

            unit = None
            has_value = self._has_value(value)
            if has_value:
                unit = self._pick_index(unit_vals, i)
            else:
                # Ensure no dangling units without value
                unit = None

            out.append((term, value if has_value else None, unit))

        return out

    @staticmethod
    def _safe_extract(source: Any, row: Any, cache: ColumnCache) -> list[Any]:
        try:
            return cache.values(source, row)
        except Exception as exc:  # pragma: no cover - defensive
            logger.debug("Source extraction failed: %s", exc)
            return []

    @staticmethod
    def _pick_index(values: list[Any] | None, index: int) -> Any | None:
        if not values:
            return None
        if len(values) == 1:
            return values[0]
        if 0 <= index < len(values):
            return values[index]
        return None

    @staticmethod
    def _has_value(value: Any) -> bool:
        if value is None:
            return False
        if isinstance(value, str):
            return value.strip() not in ("", "-")
        return True

    @staticmethod
    def _validate_term(term: Any) -> None:  # noqa: D401 - validation hook
        """Validation hook for CV terms (currently permissive)."""
        return


class IdentifiersBuilder(_BaseCvBuilder):
    """Collection of CV specifications describing identifiers.

    Example
    -------
    >>> f = FieldConfig()
    >>> identifiers = IdentifiersBuilder(
    ...     CV(
    ...         term  = f("ID Namespace", map=id_namespace_cv_mapping),
    ...         value = f("ID"),
    ...     ),
    ... )
    """

    def __init__(self, *cvs: CV) -> None:
        super().__init__(SilverIdentifier, *cvs)


class AnnotationsBuilder(_BaseCvBuilder):
    """Collection of CV specifications describing annotations.

    Example
    -------
    >>> f = FieldConfig()
    >>> annotations = AnnotationsBuilder(
    ...     CV(
    ...         term = f("Action", extract=[r"^([A-Za-z]+)"], map=action_cv_mapping),
    ...     ),
    ...     CV(
    ...         term = f("Type", map=type_cv_mapping),
    ...     ),
    ...     CV(
    ...         term = MoleculeAnnotationsCv.ENDOGENOUS,
    ...     ),
    ...     CV(
    ...         term  = MoleculeAnnotationsCv.AFFINITY_HIGH,
    ...         value = f("Affinity"),
    ...         unit  = f("Affinity Units", map=affinity_units_cv_mapping),
    ...     ),
    ... )
    """

    def __init__(self, *cvs: CV) -> None:
        super().__init__(SilverAnnotation, *cvs)


class MembersFromList:
    """Definition of member entities derived from index-aligned sources.

    This is similar in spirit to the old implementation, but driven by CV
    specifications instead of raw Columns. For each index up to the maximum
    length implied by identifiers / annotations CVs, a member is created in
    which the i-th identifier / annotation value is used (with broadcasting
    for constant sources).
    """

    def __init__(
        self,
        *,
        entity_type: EntityTypeCv,
        identifiers: IdentifiersBuilder,
        annotations: AnnotationsBuilder | None = None,
        entity_annotations: AnnotationsBuilder | None = None,
    ) -> None:
        self.entity_type = entity_type
        self.identifiers = identifiers
        self.membership_annotations = annotations
        self.entity_annotations = entity_annotations

    def build(self, row: Any, cache: ColumnCache) -> list[SilverMembership]:
        lengths: list[int] = []

        lengths.append(self.identifiers.max_length(row, cache))

        if self.membership_annotations:
            lengths.append(self.membership_annotations.max_length(row, cache))

        if self.entity_annotations:
            lengths.append(self.entity_annotations.max_length(row, cache))

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
                type=self.entity_type,
                identifiers=member_identifiers,
                annotations=entity_annotations if entity_annotations else None,
                membership=None,
            )

            memberships.append(
                SilverMembership(
                    member=member_entity,
                    annotations=membership_annotations if membership_annotations else None,
                )
            )

        return memberships


class Member:
    """Definition of a single member: entity + membership annotations."""

    def __init__(
        self,
        *,
        entity: "EntityBuilder",
        annotations: AnnotationsBuilder | None = None,
    ) -> None:
        self.entity = entity
        self.annotations = annotations

    def build(self, row: Any, cache: ColumnCache) -> SilverMembership | None:
        member_entity = self.entity.build(row)
        if not member_entity:
            return None

        membership_annot = (
            self.annotations.build(row, cache)
            if self.annotations
            else None
        )

        return SilverMembership(
            member=member_entity,
            annotations=membership_annot if membership_annot else None,
        )


class MembershipBuilder:
    """Container for one or more `Member` or `MembersFromList` definitions.

    Accepts:
    - Member: Single member with explicit entity + membership annotations
    - MembersFromList: Index-aligned members from CV-based sources
    """

    def __init__(self, *members: Member | MembersFromList) -> None:
        self.members = members

    def build(self, row: Any, cache: ColumnCache) -> list[SilverMembership]:
        memberships: list[SilverMembership] = []
        for member_def in self.members:
            if isinstance(member_def, MembersFromList):
                memberships.extend(member_def.build(row, cache))
            elif isinstance(member_def, Member):
                membership = member_def.build(row, cache)
                if membership:
                    memberships.append(membership)
        return memberships


class EntityBuilder:
    """Declarative spec that produces `Entity` records from rows."""

    def __init__(
        self,
        *,
        entity_type: EntityTypeCv | Column | Callable[[Any], EntityTypeCv],
        identifiers: IdentifiersBuilder,
        annotations: AnnotationsBuilder | None = None,
        membership: MembershipBuilder | None = None,
    ) -> None:
        self.entity_type = entity_type
        self.identifiers = identifiers
        self.annotations = annotations
        self.membership = membership

    def __call__(self, row: Any) -> SilverEntity | None:
        return self.build(row)

    def build(self, row: Any) -> SilverEntity | None:
        cache = ColumnCache()
        identifiers = self.identifiers.build(row, cache)
        if not identifiers:
            return None

        annotations = self.annotations.build(row, cache) if self.annotations else None
        membership = self.membership.build(row, cache) if self.membership else None

        # Resolve entity_type dynamically if it's a Column / callable
        resolved_type: Any = self.entity_type
        if isinstance(self.entity_type, Column):
            values = cache.values(self.entity_type, row)
            if values:
                resolved_type = values[0]
            else:
                # If extraction failed, we cannot build this entity
                logger.debug("Entity type extraction failed for row")
                return None
        elif callable(self.entity_type):
            try:
                resolved_type = self.entity_type(row)
            except Exception as exc:  # pragma: no cover - defensive
                logger.debug("Entity type callable failed: %s", exc)
                return None

        return SilverEntity(
            type=resolved_type,
            identifiers=identifiers,
            annotations=annotations if annotations else None,
            membership=membership if membership else None,
        )


__all__ = [
    "AnnotationsBuilder",
    "Column",
    "ColumnCache",
    "CV",
    "EntityBuilder",
    "IdentifiersBuilder",
    "FieldConfig",
    "Member",
    "MembershipBuilder",
    "MembersFromList",
]
