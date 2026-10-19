"""Classes a source states for its molecules.

A source with a column per level of a classification (category, class,
subclass) is modelled as terms of its own vocabulary: each level a term in the
namespace of its level, or in one namespace when the identifiers are unique
across levels (LIPID MAPS codes), ``subclass_of`` the level above it in the
same row. A molecule is ``associated_with`` every level the source states for
it; the classification makes no claim of function.

ClassyFire names are the labels of ChemOnt, the ontology ClassyFire assigns.
Such a name is the ChemOnt term whose label it is, when exactly one current
term has that label; otherwise there is no statement.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
import re
from typing import Any

from biolink_model.datamodel.model import OntologyClass, slots
from omnipath_core.naming import Namespace

from pypath.inputs_v2.base import Download
from pypath.inputs_v2.parsers.hmdb import _normalize_label
from pypath.inputs_v2.parsers.obo import iter_obo
from pypath.internals.silver_schema import Entity, EntityRef, Identifier, OntologyRelation
from pypath.internals.tabular_builder import AssociationBuilder


def plain(value: Any) -> tuple[str, str] | None:
    """A class given by its name: the name is the identifier."""
    text = ' '.join(str(value or '').split())
    return (text, text) if text else None


_CODED = re.compile(r'^(?P<name>.*?)\s*\[(?P<code>[A-Za-z0-9]+)\]$')


def coded(value: Any) -> tuple[str, str] | None:
    """A class written ``Name [CODE]`` (LIPID MAPS): the code is the identifier."""
    text = ' '.join(str(value or '').split())
    match = _CODED.match(text)
    if match:
        return match['code'], match['name']
    return None


@dataclass(frozen=True)
class Level:
    """A level of a classification: the field that holds it and its namespace."""

    field: str
    namespace: str
    parse: Callable[[Any], tuple[str, str] | None] = plain


def _values(value: Any) -> list[Any]:
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple)) else [value]


def class_associations(levels: Iterable[Level]) -> list[AssociationBuilder]:
    """``associated_with`` the term of each level a row states."""
    return [
        AssociationBuilder(
            predicate=slots.associated_with,
            object_entity_type=OntologyClass,
            object_identifier_type=level.namespace,
            object_identifier=lambda row, level=level: [
                parsed[0]
                for value in _values(row.get(level.field))
                if (parsed := level.parse(value))
            ],
        )
        for level in levels
    ]


def iter_class_terms(
    rows: Iterable[dict[str, Any]],
    levels: list[Level],
    max_records: int | None = None,
) -> Iterator[dict[str, Any]]:
    """Each term of the classification once, with the level above it in the row.

    A level a row leaves empty is skipped: the next level is below the last one
    the row states.
    """
    seen = set()
    for row in rows:
        parent = None
        for level in levels:
            parsed = level.parse(row.get(level.field))
            if not parsed:
                continue
            identifier, name = parsed
            key = (level.namespace, identifier, parent)
            if key not in seen:
                seen.add(key)
                yield {
                    'namespace': level.namespace,
                    'id': identifier,
                    'name': name,
                    'parent_namespace': parent[0] if parent else None,
                    'parent': parent[1] if parent else None,
                }
                if max_records is not None and len(seen) >= max_records:
                    return
            parent = (level.namespace, identifier)


def map_class_term(row: dict[str, Any]) -> Entity:
    """A term of a source's classification, below its parent."""
    identifiers = [Identifier(type=row['namespace'], value=row['id'])]
    identifiers.append(Identifier(type=Namespace.NAME, value=row['name']))
    relations = None
    if row.get('parent'):
        relations = [
            OntologyRelation(
                predicate=slots.subclass_of,
                object=EntityRef(
                    type=OntologyClass,
                    identifier_type=row['parent_namespace'],
                    identifier=row['parent'],
                ),
                ontology_id=row['namespace'],
            )
        ]
    return Entity(type=OntologyClass, identifiers=identifiers, ontology_relations=relations)


chemont_download = Download(
    url='http://classyfire.wishartlab.com/system/downloads/1_0/chemont/ChemOnt_2_1.obo.zip',
    filename='ChemOnt_2_1.obo.zip',
    subfolder='chemont',
    large=True,
    ext='zip',
    needed=['ChemOnt_2_1.obo'],
)


def chemont_labels(opener) -> dict[str, str]:
    """The current ChemOnt term of each label that names exactly one."""
    terms: dict[str, set[str]] = {}
    for row in iter_obo(opener):
        if row.get('is_obsolete') or not row.get('name') or not row.get('id'):
            continue
        terms.setdefault(_normalize_label(str(row['name'])), set()).add(str(row['id']))
    return {label: next(iter(ids)) for label, ids in terms.items() if len(ids) == 1}


def with_chemont(raw_parser, fields: Iterable[str], target: str = 'chemont_ids'):
    """Wrap a row parser: ``target`` holds the ChemOnt terms of the ClassyFire
    names in ``fields``."""
    fields = tuple(fields)

    def parser(opener, force_refresh: bool = False, **kwargs):
        labels = chemont_labels(chemont_download.open(force_refresh=force_refresh))
        for row in raw_parser(opener, force_refresh=force_refresh, **kwargs):
            ids = []
            for field in fields:
                for value in _values(row.get(field)):
                    term = labels.get(_normalize_label(str(value or '')))
                    if term and term not in ids:
                        ids.append(term)
            row[target] = ids
            yield row

    return parser


def chemont_association(target: str = 'chemont_ids') -> AssociationBuilder:
    """``associated_with`` the ChemOnt terms :func:`with_chemont` found."""
    return AssociationBuilder(
        predicate=slots.associated_with,
        object_entity_type=OntologyClass,
        object_identifier_type=Namespace.CHEMONT,
        object_identifier=lambda row: list(row.get(target) or ()),
    )


def _member_values(row: dict[str, Any], field: str, delimiter: str) -> list[str]:
    value = row.get(field)
    if value is None or value == '':
        return []
    return list(value) if isinstance(value, (list, tuple)) else str(value).split(delimiter)


def member_class_associations(
    levels: Iterable[Level],
    delimiter: str = '||',
    split: str | None = None,
) -> list[AssociationBuilder]:
    """``associated_with`` for the members of a row whose fields list one value
    per member, joined by ``delimiter``; ``split`` separates several classes of
    one member."""

    def ids(row, level):
        out = []
        for value in _member_values(row, level.field, delimiter):
            items = value.split(split) if split else [value]
            out.append([parsed[0] for item in items if (parsed := level.parse(item))])
        return out

    return [
        AssociationBuilder(
            predicate=slots.associated_with,
            object_entity_type=OntologyClass,
            object_identifier_type=level.namespace,
            object_identifier=lambda row, level=level: ids(row, level),
        )
        for level in levels
    ]


def member_rows(
    rows: Iterable[dict[str, Any]],
    fields: Iterable[str],
    delimiter: str = '||',
    split: str | None = None,
) -> Iterator[dict[str, Any]]:
    """One row per member and class, for :func:`iter_class_terms`."""
    fields = tuple(fields)
    for row in rows:
        columns = [_member_values(row, field, delimiter) for field in fields]
        for index in range(max(map(len, columns), default=0)):
            values = [column[index] if index < len(column) else None for column in columns]
            choices = [
                (value.split(split) if split and value else [value]) for value in values
            ]
            # Several classes of one member are siblings of the same parent.
            for combination in _product(choices):
                yield dict(zip(fields, combination))


def _product(choices: list[list[Any]]) -> Iterator[tuple[Any, ...]]:
    if not choices:
        yield ()
        return
    for head in choices[0]:
        for tail in _product(choices[1:]):
            yield (head, *tail)
