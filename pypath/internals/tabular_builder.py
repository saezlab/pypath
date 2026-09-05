#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Declarative helpers to turn tabular rows into silver-layer entities.

New DSL based on CV / FieldConfig sources (producing Columns), replacing the
older Column+cv/term_cv setup. This module intentionally breaks the previous
configuration API.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from collections import OrderedDict
from enum import Enum
import logging
import re
from typing import Any, Callable, Mapping, Sequence

from pypath.internals.silver_schema import (
    Annotation as SilverAnnotation,
    Association as SilverAssociation,
    Entity as SilverEntity,
    EntityRef as SilverEntityRef,
    Identifier as SilverIdentifier,
    Membership as SilverMembership,
    OntologyRelation as SilverOntologyRelation,
    Relation as SilverRelation,
    format_term,
)
from omnipath_core.naming import normalize_namespace
from omnipath_core.biolink import predicate as biolink_predicate, annotation_value, annotation_term
from omnipath_core.biolink import entity_type as biolink_entity_type
from biolink_model.datamodel import model

logger = logging.getLogger(__name__)
_UNCOMPUTED = object()


def _immutable_key(value: Any) -> Any:
    """Typed keys for immutable inputs; unsupported values bypass memoization."""
    if value is None or type(value) in (str, bytes, int, bool):
        return (type(value), value)
    if type(value) is float:
        if value != value:
            raise TypeError('NaN payloads cannot be memoized')
        return (float, value.hex())
    if isinstance(value, Enum):
        return (type(value), value.name)
    if isinstance(value, type):
        return (type, value)
    if type(value) is tuple:
        return (tuple, tuple(_immutable_key(item) for item in value))
    raise TypeError('Mutable or custom value cannot be memoized')


def _clean_annotation_term(term: Any) -> str:
    return annotation_term(term)


def _canonical_entity_type(value: Any) -> str:
    return biolink_entity_type(value)


def _canonical_predicate(value: Any) -> str:
    return biolink_predicate(value)


def _validate_static_entity_type(value: Any) -> None:
    if isinstance(value, Column) or (callable(value) and not isinstance(value, type)):
        return
    if _canonical_entity_type(value) is None:
        raise ValueError(
            f'EntityBuilder entity_type must not be None; got {value!r}.'
        )


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
    cache_size:
        Maximum repeated immutable cell values retained (default 4096).
        Callbacks must be pure and definitions stable. Set to zero to disable
        cell memoization; also disable enclosing entity caches for stateful code.
    """

    def __init__(
        self,
        selector: int | str | Callable[[Any], Any],
        *,
        delimiter: str | None = None,
        extract: Sequence[str | re.Pattern | Callable[[str], Any]] | None = None,
        transform: Callable[[Any], Any] | None = None,
        map: Mapping[Any, Any] | Callable[[Any], Any] | None = None,
        default: Any | None = None,
        preserve_indices: bool = False,
        cache_size: int = 4096,
    ) -> None:
        self.selector = selector
        self.delimiter = delimiter
        self.extract_steps: list[str | re.Pattern | Callable[[str], Any]] = list(extract or [])
        self.transform = transform
        self.mapping = map
        self.default = default
        self.preserve_indices = preserve_indices
        if cache_size < 0:
            raise ValueError('cache_size must not be negative')
        self.cache_size = cache_size
        self._value_cache: OrderedDict = OrderedDict()

    def extract(self, row: Any, cache: ColumnCache | None = None) -> list[Any]:
        """Extract a list of processed values from the given row.
        
        When preserve_indices is True, empty/placeholder values are preserved
        as None in the output list to maintain index alignment across multiple
        delimited fields (needed for MembersFromList).
        """
        raw_value = self._lookup(row)
        if not self.cache_size:
            return self._extract_value(raw_value)
        try:
            key = _immutable_key(raw_value)
        except TypeError:
            return self._extract_value(raw_value)
        if key in self._value_cache:
            self._value_cache.move_to_end(key)
            return list(self._value_cache[key])
        result = self._extract_value(raw_value)
        try:
            # Only retain immutable results; the outer list is copied on hits.
            _immutable_key(tuple(result))
        except TypeError:
            return result
        self._value_cache[key] = tuple(result)
        if len(self._value_cache) > self.cache_size:
            self._value_cache.popitem(last=False)
        return result


    def _extract_value(self, raw_value: Any) -> list[Any]:
        if raw_value is None:
            return []

        tokens = self._split(raw_value)
        out: list[Any] = []
        for token in tokens:
            text = self._normalize_token(token)
            if text in (None, "", "-"):
                if self.preserve_indices:
                    out.append(None)
                continue
            processed: str | Any | None = text
            for step in self.extract_steps:
                processed = self._apply_step(step, processed)
                if processed is None:
                    break

            if processed is None:
                if self.preserve_indices:
                    out.append(None)
                continue

            processed = self._apply_transform(processed)
            if processed is None:
                if self.preserve_indices:
                    out.append(None)
                continue

            mapped = self._apply_mapping(processed)
            if mapped is None:
                if self.default is not None:
                    mapped = self.default
                elif self.preserve_indices:
                    out.append(None)
                    continue
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

    def _apply_transform(self, value: Any) -> Any | None:
        if self.transform is None:
            return value

        try:
            return self.transform(value)
        except Exception as exc:  # pragma: no cover - defensive
            logger.debug("Column transform callable failed: %s", exc)
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

    __slots__ = ("value", "_dedupe_key")

    def __init__(self, value: Any) -> None:
        self.value = value
        # Constants are definition-time values. Snapshot their comparison key
        # once instead of formatting generated LinkML objects for every row.
        self._dedupe_key = _UNCOMPUTED


    def dedupe_key(self) -> Any:
        if self._dedupe_key is _UNCOMPUTED:
            self._dedupe_key = _BaseCvBuilder._make_hashable(self.value)
        return self._dedupe_key

    def extract(self, row: Any, cache: ColumnCache | None = None) -> list[Any]:  # noqa: ARG002
        return [self.value]


class _CallableSource:
    """Internal source wrapping a callable row -> value(s)."""

    __slots__ = ("func",)

    def __init__(self, func: Callable[[Any], Any]) -> None:
        self.func = func

    def extract(self, row: Any, cache: ColumnCache | None = None) -> list[Any]:  # noqa: ARG002
        result = self.func(row)

        if result is None:
            return []

        if isinstance(result, (list, tuple)):
            return list(result)

        return [result]


class _PairColumn:
    """Project aligned term/value pairs, sharing their row-local extraction."""

    def __init__(self, source: Any, index: int) -> None:
        self.source = source
        self.index = index

    def extract(self, row: Any, cache: ColumnCache) -> list[Any]:
        return [pair[self.index] for pair in cache.values(self.source, row)]


@dataclass
class FieldConfig:
    extract: dict[str, Any] = field(default_factory=dict)
    map: dict[str, Mapping[Any, Any] | Callable[[Any], Any]] = field(default_factory=dict)
    transform: dict[str, Callable[[Any], Any]] = field(default_factory=dict)
    delimiter: str | None = None
    preserve_indices: bool = False

    def __call__(
        self,
        selector: int | str | Callable[[Any], Any],
        *,
        extract: str | Sequence[Any] | Callable[[str], Any] | re.Pattern | None = None,
        transform: str | Callable[[Any], Any] | None = None,
        map: str | Mapping[Any, Any] | Callable[[Any], Any] | None = None,
        delimiter: str | None = None,
        default: Any | None = None,
        preserve_indices: bool | None = None,
        cache_size: int = 4096,
    ) -> Column:
        extract_steps = self._resolve_extract(extract)
        transform_func = self._resolve_transform(transform)
        mapping = self._resolve_mapping(map)
        return Column(
            selector,
            delimiter=delimiter if delimiter is not None else self.delimiter,
            extract=extract_steps,
            transform=transform_func,
            map=mapping,
            default=default,
            preserve_indices=preserve_indices if preserve_indices is not None else self.preserve_indices,
            cache_size=cache_size,
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

    def _resolve_transform(
        self,
        transform: str | Callable[[Any], Any] | None,
    ) -> Callable[[Any], Any] | None:
        if transform is None:
            return None
        if isinstance(transform, str):
            return self.transform[transform]
        return transform


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
        if isinstance(spec, (Column, _ConstantSource, _CallableSource, _PairColumn)):
            return spec

        # Callable row -> value(s)
        if callable(spec) and not isinstance(spec, type):
            return _CallableSource(spec)

        # Fallback: constant
        return _ConstantSource(spec)


    @classmethod
    def from_pairs(cls, extract: Column | Callable[[Any], Sequence[tuple[Any, Any]]]) -> CV:
        """Extract aligned (term, value) pairs once per row.

        Accept a row callback, or a Column whose cell transformation returns a
        sequence of pairs. Columns also reuse immutable cell transformations
        across rows. The term and value projections share one extraction.
        """
        source = _PairsSource(extract) if isinstance(extract, Column) else _CallableSource(extract)
        return cls(term=_PairColumn(source, 0), value=_PairColumn(source, 1))


class _PairsSource:
    """Flatten pair sequences produced by a Column's cell transformation."""

    def __init__(self, column: Column) -> None:
        self.column = column

    def extract(self, row: Any, cache: ColumnCache) -> list[Any]:
        return [pair for pairs in cache.values(self.column, row) for pair in pairs]


def _source_columns(source: Any) -> set[str] | None:
    """None means an opaque row callback: conservatively depend on the whole row."""
    if source is None or isinstance(source, _ConstantSource):
        return set()
    if type(source) is Column and isinstance(source.selector, str):
        return {source.selector}
    if isinstance(source, _PairColumn):
        return _source_columns(source.source)
    if isinstance(source, _PairsSource):
        return _source_columns(source.column)
    return None


def _entity_columns(*builders: Any) -> tuple[str, ...] | None:
    columns = set()
    for builder in builders:
        if builder is None:
            continue
        if type(builder) not in (IdentifiersBuilder, AnnotationsBuilder):
            return None
        for cv in builder.cvs:
            for source in (cv.term_source, cv.value_source, cv.unit_source):
                dependencies = _source_columns(source)
                if dependencies is None:
                    return None
                columns.update(dependencies)
    return tuple(sorted(columns))


def _normalize_source(spec: Any) -> Any:
    if isinstance(spec, (Column, _ConstantSource, _CallableSource)):
        return spec
    if callable(spec) and not isinstance(spec, type):
        return _CallableSource(spec)
    return _ConstantSource(spec)


class _BaseCvBuilder:
    """Shared implementation for identifier and annotation builders."""

    def __init__(self, silver_cls: type[SilverIdentifier] | type[SilverAnnotation], *cvs: CV) -> None:
        self.silver_cls = silver_cls
        self.cvs = cvs
        self._constant_terms: dict[CV, str] = {}

    # -- public API -------------------------------------------------------

    def build(self, row: Any, cache: ColumnCache | None = None) -> list[SilverIdentifier] | list[SilverAnnotation]:
        if cache is None:
            cache = ColumnCache()
        results: list[SilverIdentifier] | list[SilverAnnotation] = []
        seen: set[tuple[Any, Any, Any]] = set()

        for cv in self.cvs:
            for term, value, unit in self._expand_cv(cv, row, cache):
                if term is None:
                    continue

                key = self._cv_dedupe_key(cv, term, value, unit if self.silver_cls is SilverAnnotation else None)
                if key in seen:
                    continue
                seen.add(key)

                if self.silver_cls is SilverIdentifier:
                    # Identifiers ignore unit and require a non-empty value
                    if value is None or value == "":
                        continue
                    clean_term = self._clean_term(cv, term)
                    results.append(SilverIdentifier(type=clean_term, value=str(value)))
                else:
                    clean_term = self._clean_term(cv, term)
                    results.append(SilverAnnotation(term=clean_term, value=annotation_value(term, value), units=str(unit) if unit is not None else None))

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

            if cv.value_source is None:
                if self.silver_cls is not SilverAnnotation:
                    continue
                key = self._cv_dedupe_key(cv, term, None, None)
                if key in seen:
                    continue
                seen.add(key)
                results.append(SilverAnnotation(term=annotation_term(term), value=None, units=None))
                continue

            value = self._pick_index(value_vals, index) if value_vals is not None else None
            value_items = self._explode_values(value)
            if not value_items:
                continue

            unit = self._pick_index(unit_vals, index) if unit_vals is not None else None
            unit_items = self._explode_values(unit) if cv.unit_source is not None else [None]
            if len(unit_items) != len(value_items):
                unit_items = [unit_items[0] if unit_items else None] * len(value_items)

            for value_item, unit_item in zip(value_items, unit_items):
                key = self._cv_dedupe_key(cv, term, value_item, unit_item if self.silver_cls is SilverAnnotation else None)
                if key in seen:
                    continue
                seen.add(key)

                if self.silver_cls is SilverIdentifier:
                    if value_item is None or value_item == "":
                        continue
                    clean_term = self._clean_term(cv, term)
                    results.append(SilverIdentifier(type=clean_term, value=str(value_item)))
                else:
                    clean_term = self._clean_term(cv, term)
                    results.append(SilverAnnotation(term=clean_term, value=annotation_value(term, value_item), units=str(unit_item) if unit_item is not None else None))

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

    def _clean_term(self, cv: CV, term: Any) -> str:
        if isinstance(cv.term_source, _ConstantSource):
            if cv in self._constant_terms:
                return self._constant_terms[cv]
        result = (
            normalize_namespace(term) or str(term)
            if self.silver_cls is SilverIdentifier
            else _clean_annotation_term(term)
        )
        if isinstance(cv.term_source, _ConstantSource):
            self._constant_terms[cv] = result
        return result

    def _expand_cv(
        self,
        cv: CV,
        row: Any,
        cache: ColumnCache,
    ) -> list[tuple[Any, Any, Any]]:
        # A static term with no units is the common scalar/list mapping case.
        # It needs neither three-column broadcasting nor temporary unit lists.
        if isinstance(cv.term_source, _ConstantSource) and cv.unit_source is None:
            term = cv.term_source.value
            if term is None:
                return []
            if cv.value_source is None:
                self._validate_term(term)
                return [(term, None, None)]
            result = []
            for value in cache.values(cv.value_source, row):
                self._validate_term(term)
                for item in self._explode_values(value):
                    if self._has_value(item):
                        result.append((term, item, None))
            return result

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
            unit = self._pick_index(unit_vals, i) if unit_vals else None

            if cv.value_source is None:
                out.append((term, None, None))
                continue

            value_items = self._explode_values(value)
            if not value_items:
                continue

            unit_items = self._explode_values(unit) if cv.unit_source is not None else [None]
            if len(unit_items) != len(value_items):
                unit_items = [unit_items[0] if unit_items else None] * len(value_items)

            for idx in range(len(value_items)):
                value_item = value_items[idx]
                unit_item = unit_items[idx] if idx < len(unit_items) else None
                if not self._has_value(value_item):
                    continue
                out.append((term, value_item, unit_item))

        return out

    @staticmethod
    def _safe_extract(source: Any, row: Any, cache: ColumnCache) -> list[Any]:
        return cache.values(source, row)

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
        if isinstance(value, (list, tuple, set)):
            return any(_BaseCvBuilder._has_value(item) for item in value)
        if isinstance(value, str):
            return value.strip() not in ("", "-")
        return True

    @staticmethod
    def _explode_values(value: Any) -> list[Any]:
        if value is None:
            return []
        if isinstance(value, (list, tuple, set)):
            return [item for item in value if _BaseCvBuilder._has_value(item)]
        if _BaseCvBuilder._has_value(value):
            return [value]
        return []

    @staticmethod
    def _cv_dedupe_key(cv: CV, term: Any, value: Any, unit: Any) -> tuple:
        term_key = (
            cv.term_source.dedupe_key()
            if isinstance(cv.term_source, _ConstantSource)
            else _BaseCvBuilder._make_hashable(term)
        )
        value_key = (
            cv.value_source.dedupe_key()
            if isinstance(cv.value_source, _ConstantSource) and value is cv.value_source.value
            else _BaseCvBuilder._make_hashable(value)
        )
        unit_key = (
            cv.unit_source.dedupe_key()
            if isinstance(cv.unit_source, _ConstantSource) and unit is cv.unit_source.value
            else _BaseCvBuilder._make_hashable(unit)
        )
        return term_key, value_key, unit_key


    @staticmethod
    def _dedupe_key(term: Any, value: Any, unit: Any) -> tuple[Any, Any, Any]:
        return (
            _BaseCvBuilder._make_hashable(term),
            _BaseCvBuilder._make_hashable(value),
            _BaseCvBuilder._make_hashable(unit),
        )

    @staticmethod
    def _make_hashable(value: Any) -> Any:
        if value is None or type(value) in (str, int, float, bool):
            return value
        try:
            hash(value)
            return value
        except TypeError:
            pass

        if hasattr(value, '_asdict'):
            return tuple(
                (key, _BaseCvBuilder._make_hashable(item))
                for key, item in value._asdict().items()
            )
        if isinstance(value, dict):
            return tuple(
                sorted(
                    (_BaseCvBuilder._make_hashable(key), _BaseCvBuilder._make_hashable(item))
                    for key, item in value.items()
                )
            )
        if isinstance(value, (list, tuple, set)):
            return tuple(_BaseCvBuilder._make_hashable(item) for item in value)

        return repr(value)

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
        entity_type: type[model.NamedThing] | Column | Callable[[Any], Any],
        predicate: Any = model.slots.has_member,
        identifiers: IdentifiersBuilder,
        annotations: AnnotationsBuilder | None = None,
        associations: "AssociationsBuilder" | None = None,
        entity_annotations: AnnotationsBuilder | None = None,
        entity_associations: "AssociationsBuilder" | None = None,
    ) -> None:
        _validate_static_entity_type(entity_type)
        self.predicate = (predicate if isinstance(predicate, Column) or callable(predicate)
                          else _canonical_predicate(predicate))
        self.entity_type = entity_type
        self.identifiers = identifiers
        self.membership_annotations = annotations
        self.membership_associations = associations
        self.entity_annotations = entity_annotations
        self.entity_associations = entity_associations

    def build(self, row: Any, cache: ColumnCache) -> list[SilverMembership]:
        lengths: list[int] = []

        lengths.append(self.identifiers.max_length(row, cache))

        if self.membership_annotations:
            lengths.append(self.membership_annotations.max_length(row, cache))

        if self.membership_associations:
            lengths.append(self.membership_associations.max_length(row, cache))

        if self.entity_annotations:
            lengths.append(self.entity_annotations.max_length(row, cache))

        if self.entity_associations:
            lengths.append(self.entity_associations.max_length(row, cache))

        member_count = max(lengths) if lengths else 0
        memberships: list[SilverMembership] = []

        for index in range(member_count):
            predicate = self._predicate_for_index(row, cache, index)
            if predicate is None:
                continue
            member_identifiers = self.identifiers.build_for_index(row, index, cache)
            if not member_identifiers:
                continue

            member_entity_type = self._resolve_entity_type(row, cache, index)
            if member_entity_type is None:
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
            entity_associations = (
                self.entity_associations.build_for_index(row, index, cache)
                if self.entity_associations
                else None
            )
            membership_associations = (
                self.membership_associations.build_for_index(row, index, cache)
                if self.membership_associations
                else None
            )

            member_entity = SilverEntity(
                type=member_entity_type,
                identifiers=member_identifiers,
                annotations=entity_annotations if entity_annotations else None,
                associations=entity_associations if entity_associations else None,
                membership=None,
            )

            memberships.append(
                SilverMembership(
                    member=member_entity,
                    predicate=predicate,
                    annotations=membership_annotations if membership_annotations else None,
                    associations=membership_associations if membership_associations else None,
                )
            )

        return memberships

    def _predicate_for_index(self, row: Any, cache: ColumnCache, index: int) -> str | None:
        values = (cache.values(self.predicate, row) if isinstance(self.predicate, Column)
                  else self.predicate(row) if callable(self.predicate) else self.predicate)
        if isinstance(values, (list, tuple)):
            values = _BaseCvBuilder._pick_index(values, index)
        return None if values is None else _canonical_predicate(values)

    def _resolve_entity_type(
        self,
        row: Any,
        cache: ColumnCache,
        index: int,
    ) -> str:
        if isinstance(self.entity_type, Column):
            values = cache.values(self.entity_type, row)
            value = _BaseCvBuilder._pick_index(values, index)
        elif isinstance(self.entity_type, type):
            value = self.entity_type
        elif callable(self.entity_type):
            try:
                value = self.entity_type(row)
            except Exception as exc:  # pragma: no cover - defensive
                logger.debug("MembersFromList entity_type callable failed: %s", exc)
                return None
            if isinstance(value, (list, tuple)):
                value = _BaseCvBuilder._pick_index(list(value), index)
        else:
            value = self.entity_type

        entity_type = _canonical_entity_type(value)
        if entity_type is None:
            logger.debug("Invalid member entity_type: %r", value)
        return entity_type


class Member:
    """Definition of a single member: entity + membership annotations."""

    def __init__(
        self,
        *,
        entity: "EntityBuilder",
        predicate: Any = model.slots.has_member,
        annotations: AnnotationsBuilder | None = None,
        associations: "AssociationsBuilder" | None = None,
    ) -> None:
        self.predicate = (predicate if isinstance(predicate, Column) or callable(predicate)
                          else _canonical_predicate(predicate))
        self.entity = entity
        self.annotations = annotations
        self.associations = associations

    def build(self, row: Any, cache: ColumnCache) -> SilverMembership | None:
        member_entity = self.entity.build(row)
        if not member_entity:
            return None

        membership_annot = (
            self.annotations.build(row, cache)
            if self.annotations
            else None
        )
        membership_associations = (
            self.associations.build(row, cache)
            if self.associations
            else None
        )

        return SilverMembership(
            member=member_entity,
            predicate=_canonical_predicate(
                cache.values(self.predicate, row)[0] if isinstance(self.predicate, Column)
                else self.predicate(row) if callable(self.predicate) else self.predicate
            ),
            annotations=membership_annot if membership_annot else None,
            associations=membership_associations if membership_associations else None,
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


class OntologyRelationBuilder:
    """Definition of ontology hierarchy edges from the built entity."""

    def __init__(
        self,
        *,
        predicate: Any,
        object_entity_type: Any,
        object_identifier_type: Any,
        object_identifier: Any,
        ontology_id: Any | None = None,
    ) -> None:
        self.predicate_source = _normalize_source(predicate)
        self.object_entity_type_source = _normalize_source(object_entity_type)
        self.object_identifier_type_source = _normalize_source(
            object_identifier_type,
        )
        self.object_identifier_source = _normalize_source(object_identifier)
        self.ontology_id_source = (
            _normalize_source(ontology_id)
            if ontology_id is not None
            else None
        )

    def build(
        self,
        row: Any,
        cache: ColumnCache,
    ) -> list[SilverOntologyRelation]:
        predicates = cache.values(self.predicate_source, row)
        object_entity_types = cache.values(self.object_entity_type_source, row)
        object_identifier_types = cache.values(
            self.object_identifier_type_source,
            row,
        )
        object_identifiers = cache.values(self.object_identifier_source, row)
        ontology_ids = (
            cache.values(self.ontology_id_source, row)
            if self.ontology_id_source is not None
            else []
        )
        relations: list[SilverOntologyRelation] = []
        seen: set[tuple[str, str, str, str, str | None]] = set()
        for index, object_identifier in enumerate(object_identifiers):
            predicate = _pick_index_or_first(predicates, index)
            object_entity_type = _pick_index_or_first(object_entity_types, index)
            object_identifier_type = _pick_index_or_first(
                object_identifier_types,
                index,
            )
            ontology_id = _pick_index_or_first(ontology_ids, index)
            if not predicate or not object_identifier:
                continue
            if object_entity_type is None or object_identifier_type is None:
                continue
            key = (
                str(predicate),
                str(object_entity_type),
                str(object_identifier_type),
                str(object_identifier),
                str(ontology_id) if ontology_id is not None else None,
            )
            if key in seen:
                continue
            seen.add(key)
            relations.append(
                SilverOntologyRelation(
                    predicate=_canonical_predicate(predicate),
                    object=SilverEntityRef(
                        type=_canonical_entity_type(object_entity_type),
                        identifier_type=normalize_namespace(object_identifier_type) or str(object_identifier_type),
                        identifier=str(object_identifier),
                    ),
                    ontology_id=str(ontology_id) if ontology_id else None,
                )
            )
        return relations


class OntologyRelationsBuilder:
    """Container for one or more ontology relation builders."""

    def __init__(self, *relations: OntologyRelationBuilder) -> None:
        self.relations = relations

    def build(
        self,
        row: Any,
        cache: ColumnCache,
    ) -> list[SilverOntologyRelation]:
        values: list[SilverOntologyRelation] = []
        for relation in self.relations:
            values.extend(relation.build(row, cache))
        return values


def _pick_index_or_first(values: list[Any], index: int) -> Any | None:
    if not values:
        return None
    if index < len(values):
        return values[index]
    return values[0]


class AssociationBuilder:
    """Definition of graph associations from the built entity to an object."""

    def __init__(
        self,
        *,
        object_entity_type: Any,
        object_identifier_type: Any,
        object_identifier: Any,
        predicate: Any | None = None,
    ) -> None:
        self.predicate_source = (
            _normalize_source(predicate)
            if predicate is not None
            else None
        )
        self.object_entity_type_source = _normalize_source(object_entity_type)
        self.object_identifier_type_source = _normalize_source(
            object_identifier_type,
        )
        self.object_identifier_source = _normalize_source(object_identifier)

    def build(
        self,
        row: Any,
        cache: ColumnCache,
    ) -> list[SilverAssociation]:
        max_len = self.max_length(row, cache)
        associations: list[SilverAssociation] = []
        seen: set[tuple[str | None, str, str, str]] = set()
        for index in range(max_len):
            for association in self.build_for_index(row, index, cache):
                key = (
                    str(association.predicate) if association.predicate else None,
                    str(association.object.type),
                    str(association.object.identifier_type),
                    str(association.object.identifier),
                )
                if key in seen:
                    continue
                seen.add(key)
                associations.append(association)
        return associations

    def build_for_index(
        self,
        row: Any,
        index: int,
        cache: ColumnCache,
    ) -> list[SilverAssociation]:
        predicates = (
            cache.values(self.predicate_source, row)
            if self.predicate_source is not None
            else []
        )
        object_entity_types = cache.values(self.object_entity_type_source, row)
        object_identifier_types = cache.values(
            self.object_identifier_type_source,
            row,
        )
        object_identifiers = cache.values(self.object_identifier_source, row)
        predicate = _pick_index_or_first(predicates, index)
        object_entity_type = _pick_index_or_first(object_entity_types, index)
        object_identifier_type = _pick_index_or_first(
            object_identifier_types,
            index,
        )
        object_identifier = _BaseCvBuilder._pick_index(object_identifiers, index)
        if object_entity_type is None or object_identifier_type is None:
            return []

        associations: list[SilverAssociation] = []
        seen: set[tuple[str | None, str, str, str]] = set()
        for object_identifier_item in _BaseCvBuilder._explode_values(
            object_identifier
        ):
            key = (
                str(predicate) if predicate else None,
                str(object_entity_type),
                str(object_identifier_type),
                str(object_identifier_item),
            )
            if key in seen:
                continue
            seen.add(key)
            associations.append(
                SilverAssociation(
                    predicate=_canonical_predicate(predicate) if predicate else None,
                    object=SilverEntityRef(
                        type=_canonical_entity_type(object_entity_type),
                        identifier_type=normalize_namespace(object_identifier_type) or str(object_identifier_type),
                        identifier=str(object_identifier_item),
                    ),
                )
            )
        return associations

    def max_length(self, row: Any, cache: ColumnCache) -> int:
        lengths = [
            len(cache.values(self.object_entity_type_source, row)),
            len(cache.values(self.object_identifier_type_source, row)),
            len(cache.values(self.object_identifier_source, row)),
        ]
        if self.predicate_source is not None:
            lengths.append(len(cache.values(self.predicate_source, row)))
        return max(lengths)


class AssociationsBuilder:
    """Container for one or more association builders."""

    def __init__(self, *associations: AssociationBuilder) -> None:
        self.associations = associations

    def build(
        self,
        row: Any,
        cache: ColumnCache,
    ) -> list[SilverAssociation]:
        values: list[SilverAssociation] = []
        seen: set[tuple[str | None, str, str, str]] = set()
        for association_builder in self.associations:
            for association in association_builder.build(row, cache):
                key = (
                    str(association.predicate) if association.predicate else None,
                    str(association.object.type),
                    str(association.object.identifier_type),
                    str(association.object.identifier),
                )
                if key in seen:
                    continue
                seen.add(key)
                values.append(association)
        return values

    def build_for_index(
        self,
        row: Any,
        index: int,
        cache: ColumnCache,
    ) -> list[SilverAssociation]:
        values: list[SilverAssociation] = []
        for association_builder in self.associations:
            values.extend(association_builder.build_for_index(row, index, cache))
        return values

    def max_length(self, row: Any, cache: ColumnCache) -> int:
        lengths = [
            association_builder.max_length(row, cache)
            for association_builder in self.associations
        ]
        return max(lengths) if lengths else 0


class EntityBuilder:
    """Declarative spec that produces `Entity` records from rows.

    Flat entities automatically reuse pure mappings, inferring dependencies
    from identifier/annotation columns. Opaque callbacks depend on the whole
    dictionary row. Dynamic types are evaluated on every row. ``cache_size=0``
    disables entity reuse; fields have their own independent cache setting.
    """

    def __init__(
        self,
        *,
        entity_type: type[model.NamedThing] | Column | Callable[[Any], type[model.NamedThing]],
        identifiers: IdentifiersBuilder | None = None,
        annotations: AnnotationsBuilder | None = None,
        associations: AssociationsBuilder | None = None,
        membership: MembershipBuilder | None = None,
        ontology_relations: OntologyRelationsBuilder
        | Callable[[Any], list[SilverOntologyRelation]]
        | None = None,
        cache_by: Sequence[str] | None = None,
        cache_size: int = 4096,
    ) -> None:
        _validate_static_entity_type(entity_type)
        self.entity_type = entity_type
        self.identifiers = identifiers
        self.annotations = annotations
        self.associations = associations
        self.membership = membership
        self.ontology_relations = ontology_relations
        if isinstance(cache_by, str):
            raise TypeError('cache_by must be a sequence of column names, not a string')
        if cache_by is not None and (associations or membership or ontology_relations):
            raise ValueError('cache_by is supported for flat entities only')
        if cache_size < 0:
            raise ValueError('cache_size must not be negative')
        self.cache_by = tuple(cache_by) if cache_by is not None else None
        self._inferred_columns = _entity_columns(identifiers, annotations)
        standard_fields = (
            (identifiers is None or type(identifiers) is IdentifiersBuilder)
            and (annotations is None or type(annotations) is AnnotationsBuilder)
        )
        self.cache_size = (
            cache_size if standard_fields and not (associations or membership or ontology_relations) else 0
        )
        self._entity_cache: OrderedDict[tuple, SilverEntity | None] = OrderedDict()

    def __call__(self, row: Any) -> SilverEntity | None:
        return self.build(row)

    def build(self, row: Any) -> SilverEntity | None:
        if not self.cache_size or type(row) is not dict:
            return self._build(row)
        # Evaluate dynamic type on every call. It may depend on arbitrary fields
        # and must still validate a row even when identifiers/annotations repeat.
        cache = ColumnCache()
        resolved_type = self._resolve_type(row, cache)
        if resolved_type is None:
            return None
        columns = self.cache_by if self.cache_by is not None else self._inferred_columns
        try:
            if columns is None:
                # Preserve iteration order for opaque callbacks as well as values.
                values = tuple((_immutable_key(k), _immutable_key(v)) for k, v in row.items())
            else:
                values = tuple((name in row, _immutable_key(row.get(name))) for name in columns)
            key = (resolved_type, values)
        except TypeError:
            return self._build(row, resolved_type, cache)
        if key in self._entity_cache:
            self._entity_cache.move_to_end(key)
            entity = self._entity_cache[key]
        else:
            entity = self._build(row, resolved_type, cache)
            self._entity_cache[key] = entity
            if len(self._entity_cache) > self.cache_size:
                self._entity_cache.popitem(last=False)
        # Flat Identifier/Annotation records contain immutable scalar values.
        # Copy their containers so consumers cannot mutate later results.
        if entity is None:
            return None
        return entity._replace(
            identifiers=list(entity.identifiers),
            annotations=list(entity.annotations) if entity.annotations is not None else None,
        )


    def _resolve_type(self, row: Any, cache: ColumnCache) -> str | None:
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
        elif isinstance(self.entity_type, type):
            resolved_type = self.entity_type
        elif callable(self.entity_type):
            resolved_type = self.entity_type(row)

        resolved_type = _canonical_entity_type(resolved_type)
        if resolved_type is None:
            logger.debug("Invalid entity_type for row")
            return None
        return resolved_type


    def _build(
        self, row: Any, resolved_type: Any = _UNCOMPUTED, cache: ColumnCache | None = None,
    ) -> SilverEntity | None:
        if cache is None:
            cache = ColumnCache()
        if resolved_type is _UNCOMPUTED:
            resolved_type = self._resolve_type(row, cache)
        if resolved_type is None:
            return None

        identifiers = self.identifiers.build(row, cache) if self.identifiers else []
        if not identifiers:
            return None

        annotations = self.annotations.build(row, cache) if self.annotations else None
        associations = self.associations.build(row, cache) if self.associations else None
        membership = self.membership.build(row, cache) if self.membership else None
        ontology_relations = self._build_ontology_relations(row, cache)

        return SilverEntity(
            type=resolved_type,
            identifiers=identifiers,
            annotations=annotations if annotations else None,
            associations=associations if associations else None,
            membership=membership if membership else None,
            ontology_relations=(
                ontology_relations if ontology_relations else None
            ),
        )

    def _build_ontology_relations(
        self,
        row: Any,
        cache: ColumnCache,
    ) -> list[SilverOntologyRelation]:
        if self.ontology_relations is None:
            return []
        if hasattr(self.ontology_relations, 'build'):
            return self.ontology_relations.build(row, cache)
        try:
            return list(self.ontology_relations(row) or [])
        except Exception as exc:  # pragma: no cover - defensive
            logger.debug("Ontology relations callable failed: %s", exc)
            return []



class RelationBuilder:
    """Declarative spec that produces `Relation` records from rows."""

    def __init__(
        self,
        *,
        subject: EntityBuilder | Callable[[Any], SilverEntity | SilverEntityRef | None],
        predicate: str | Column | Callable[[Any], str | None],
        object: EntityBuilder | Callable[[Any], SilverEntity | SilverEntityRef | None],
        category: str | Column | Callable[[Any], str | None] | None = None,
        interaction_class: str | Column | Callable[[Any], str | None] | None = None,
        identifiers: IdentifiersBuilder | None = None,
        annotations: AnnotationsBuilder | None = None,
    ) -> None:
        self.subject = subject
        self.predicate = predicate
        self.object = object
        self.category = category
        self.interaction_class = interaction_class
        self.identifiers = identifiers
        self.annotations = annotations

    def __call__(self, row: Any) -> SilverRelation | None:
        return self.build(row)

    def build(self, row: Any, *, emit: Callable[..., Any] | None = None) -> Any:
        """Evaluate the mapping, optionally sending its fields to an execution sink.

        Ordinary callers receive a SilverRelation. The optional factory allows
        executor experiments to reuse exactly the same mapping and validation.
        """
        cache = ColumnCache()

        # Build subject entity
        if isinstance(self.subject, EntityBuilder):
            subject_entity = self.subject.build(row)
        elif callable(self.subject):
            try:
                subject_entity = self.subject(row)
            except Exception as exc:  # pragma: no cover - defensive
                logger.debug("Relation subject callable failed: %s", exc)
                return None
        else:
            subject_entity = self.subject

        if subject_entity is None:
            return None

        # Build object entity
        if isinstance(self.object, EntityBuilder):
            object_entity = self.object.build(row)
        elif callable(self.object):
            try:
                object_entity = self.object(row)
            except Exception as exc:  # pragma: no cover - defensive
                logger.debug("Relation object callable failed: %s", exc)
                return None
        else:
            object_entity = self.object

        if object_entity is None:
            return None

        # Invalid or missing predicates are source modeling errors, never fallback edges.
        if isinstance(self.predicate, Column):
            pred_vals = cache.values(self.predicate, row)
            predicate_val = _canonical_predicate(pred_vals[0] if pred_vals else None)
        elif callable(self.predicate):
            predicate_val = _canonical_predicate(self.predicate(row))
        else:
            predicate_val = _canonical_predicate(self.predicate)

        # Resolve optional category
        category_val = None
        if isinstance(self.category, Column):
            cat_vals = cache.values(self.category, row)
            category_val = str(cat_vals[0]) if cat_vals else None
        elif callable(self.category):
            try:
                category_val = self.category(row)
            except Exception:  # pragma: no cover
                pass
        elif self.category is not None:
            category_val = str(self.category)

        # Resolve optional interaction_class
        class_val = None
        if isinstance(self.interaction_class, Column):
            cls_vals = cache.values(self.interaction_class, row)
            class_val = str(cls_vals[0]) if cls_vals else None
        elif callable(self.interaction_class):
            try:
                class_val = self.interaction_class(row)
            except Exception:  # pragma: no cover
                pass
        elif self.interaction_class is not None:
            class_val = str(self.interaction_class)

        identifiers = self.identifiers.build(row, cache) if self.identifiers else None
        annotations = self.annotations.build(row, cache) if self.annotations else None

        factory = emit if emit is not None else SilverRelation
        return factory(
            subject=subject_entity,
            predicate=predicate_val,
            object=object_entity,
            category=category_val,
            interaction_class=class_val,
            identifiers=identifiers if identifiers else None,
            annotations=annotations if annotations else None,
        )


__all__ = [
    "AssociationBuilder",
    "AssociationsBuilder",
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
    "RelationBuilder",
]
