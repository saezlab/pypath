"""Biolink boundary: schema-derived names, validation, and qualified statements.

Storage names are generated from schema names (lowercase, spaces to underscores).
CURIEs, schema names and generated Python primitives resolve through the same index.
No source vocabulary or source identifier repair belongs in this module.
"""

from __future__ import annotations

from functools import lru_cache
from enum import Enum
import re
from importlib.resources import files
from typing import Any

from linkml_runtime.linkml_model.meta import PermissibleValue
from linkml_runtime.utils.enumerations import EnumDefinitionImpl
from linkml_runtime.utils.schemaview import SchemaView


@lru_cache(maxsize=1)
def schema() -> SchemaView:
    return SchemaView(str(files("biolink_model").joinpath("schema/biolink_model.yaml")))


def _storage_name(name: str) -> str:
    return name.lower().replace(" ", "_")


@lru_cache(maxsize=2)
def _index(kind: str) -> dict[str, str]:
    result = {}
    elements = schema().all_classes() if kind == "class" else schema().all_slots()
    for name, element in elements.items():
        curie = schema().get_uri(element, expand=False)
        for key in (name, _storage_name(name), curie):
            result[str(key).lower()] = name
    return result


def _name(value: Any, kind: str) -> str:
    if type(value) is str:
        return _string_name(value, kind)
    if isinstance(value, type) and hasattr(value, "class_class_curie"):
        value = value.class_class_curie
    elif not isinstance(value, (str, EnumDefinitionImpl, PermissibleValue)):
        value = getattr(value, "curie", value)
    if not isinstance(value, str):
        raise ValueError(f"Expected a Biolink {kind}, got {value!r}")
    return _string_name(value, kind)


@lru_cache(maxsize=8192)
def _string_name(value: str, kind: str) -> str:
    """Cache only immutable names, never mutable generated schema objects."""
    try:
        return _index(kind)[value.strip().lower()]
    except KeyError:
        raise ValueError(
            f"Unknown Biolink {kind}: {value!r}; adapt it at the source boundary"
        ) from None


def entity_type(value: Any) -> str:
    return _storage_name(_name(value, "class"))


def slot_name(value: Any) -> str:
    return _storage_name(_name(value, "slot"))


def class_curie(value: Any) -> str:
    return str(schema().get_uri(schema().get_class(_name(value, "class")), expand=False))


def slot_curie(value: Any) -> str:
    return str(schema().get_uri(schema().get_slot(_name(value, "slot")), expand=False))


def predicate(value: Any) -> str:
    name = _name(value, "slot")
    return _predicate_name(name)


@lru_cache(maxsize=2048)
def _predicate_name(name: str) -> str:
    if "related to" not in schema().slot_ancestors(name, reflexive=True):
        raise ValueError(f"Biolink slot {name!r} is not a predicate")
    if schema().get_slot(name).deprecated:
        raise ValueError(f"Deprecated Biolink predicate {name!r}; adapt the source modeling")
    return _storage_name(name)


def enum_value(value: Any) -> str:
    if isinstance(value, PermissibleValue):
        return str(value.text)
    if isinstance(value, EnumDefinitionImpl):
        return str(value.code.text)
    return str(value)


def annotation_term(value: Any) -> str:
    """Biolink slot or an external ontology attribute CURIE, preserving its case."""
    if type(value) is str:
        return _annotation_text(value)
    return _annotation_object(value)


@lru_cache(maxsize=8192)
def _annotation_text(value: str) -> str:
    return _annotation_object(value)


def _annotation_object(value: Any) -> str:
    if isinstance(value, Enum):
        raise ValueError(
            "Annotation terms must use Biolink slots or ontology CURIEs, not legacy enums"
        )
    try:
        return slot_name(value)
    except ValueError:
        if isinstance(value, str) and re.fullmatch(r"[A-Za-z][A-Za-z0-9._-]*:[^\s:]+", value):
            if value.startswith("biolink:"):
                raise
            return value
        raise ValueError(
            f"Unknown annotation term {value!r}; use a Biolink slot or ontology CURIE"
        ) from None


def annotation_value(term: Any, value: Any) -> str | None:
    if value is None:
        return None
    from .measurements import quantity_json

    quantity = quantity_json(value)
    if quantity is not None:
        try:
            range_name, allowed = _annotation_range(_name(term, "slot"))
        except ValueError:
            pass  # published external attribute type
        else:
            if allowed is not None or range_name not in {
                "quantity value",
                "double",
                "float",
                "integer",
                "decimal",
            }:
                raise ValueError(
                    f"Quantities require a quantitative slot or external attribute type, got {term!r}"
                )
        return quantity
    text = enum_value(value)
    try:
        name = _name(term, "slot")
    except ValueError:
        return text  # external ontology attributes remain evidence, not Biolink slots
    range_name, allowed = _annotation_range(name)
    if allowed is not None and text not in allowed:
        raise ValueError(f"Invalid {range_name} value {text!r} for {name!r}")
    if name == "qualified predicate":
        return predicate(value)
    return text


@lru_cache(maxsize=2048)
def _annotation_range(name: str):
    slot = schema().induced_slot(name)
    enum = schema().get_enum(slot.range)
    return slot.range, frozenset(enum.permissible_values) if enum else None


def is_symmetric(value: Any) -> bool:
    return _symmetric_name(_name(predicate(value), "slot"))


@lru_cache(maxsize=2048)
def _symmetric_name(name: str) -> bool:
    return schema().induced_slot(name).symmetric is True


def is_descendant(value: Any, ancestor: Any) -> bool:
    return _descendant_names(_name(value, "slot"), _name(ancestor, "slot"))


@lru_cache(maxsize=8192)
def _descendant_names(name: str, ancestor: str) -> bool:
    return ancestor in schema().slot_ancestors(name, reflexive=True)


def qualifiers(annotations: list[dict] | None) -> tuple[tuple[str, str], ...]:
    result = {}
    for annotation in annotations or ():
        try:
            name = _name(annotation.get("term"), "slot")
        except ValueError:
            continue
        if not _descendant_names(name, "qualifier"):
            continue
        term = _storage_name(name)
        if annotation.get("scope", "relation") != "relation":
            raise ValueError("Biolink edge qualifiers must have relation scope")
        value = annotation_value(term, annotation.get("value"))
        if value is None:
            raise ValueError(f"Missing value for qualifier {term}")
        if term in result and result[term] != value:
            raise ValueError(f"Conflicting values for qualifier {term}")
        result[term] = value
    return tuple(sorted(result.items()))


def statement_identity(subject: str, relation: Any, object: str, annotations=None):
    relation = predicate(relation)
    qualified = qualifiers(annotations)
    # Qualifiers can distinguish subject and object even on a symmetric base edge.
    if is_symmetric(relation) and not qualified:
        subject, object = sorted((subject, object))
    return subject, relation, object, qualified


def direction_sign(annotations: list[dict] | None) -> int:
    return direction_sign_from_qualifiers(qualifiers(annotations))


def direction_sign_from_qualifiers(qualified: tuple[tuple[str, str], ...]) -> int:
    """Project direction from the validated qualifier tuple of a statement."""
    direction = dict(qualified).get("object_direction_qualifier")
    return _direction_sign(direction)


@lru_cache(maxsize=128)
def _direction_sign(direction: str | None) -> int:
    if direction is None:
        return 0
    values = schema().get_enum("DirectionQualifierEnum").permissible_values
    while values[direction].is_a:
        direction = values[direction].is_a
    # Integer signs are a serving projection, not a second biological vocabulary.
    return {"increased": 1, "decreased": -1}[direction]


def presentation_category(value: Any) -> str:
    """Documented legacy UI grouping; never used to infer biological semantics."""
    return _presentation_category_name(_name(value, "slot"))


@lru_cache(maxsize=2048)
def _presentation_category_name(name: str) -> str:
    for root, category in (
        ("has part", "membership"),
        ("part of", "membership"),
        ("has member", "membership"),
        ("member of", "membership"),
        ("subclass of", "ontology"),
        ("superclass of", "ontology"),
        ("participates in", "pathway"),
        ("affects", "interaction"),
        ("interacts with", "interaction"),
    ):
        if _descendant_names(name, root):
            return category
    return "association"


def hierarchy_direction(value: Any) -> bool | None:
    """Direction of an explicitly identified ontology axiom (child -> parent)."""
    if value == "is_a":
        return False
    try:
        if is_descendant(value, "subclass_of") or is_descendant(value, "part_of"):
            return False
        if is_descendant(value, "superclass_of") or is_descendant(value, "has_part"):
            return True
    except ValueError:
        pass
    return None
