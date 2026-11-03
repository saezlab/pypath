"""Contract tests for native Biolink primitives and qualified statement identity."""

import pytest
from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core import format_term, relation_key
from omnipath_core.biolink import (
    annotation_term,
    annotation_value,
    class_curie,
    direction_sign,
    entity_type,
    is_symmetric,
    predicate,
    presentation_category,
    qualifiers,
)


@pytest.mark.parametrize("value", [model.Protein, "protein", "biolink:Protein"])
def test_class_representations_are_identical(value):
    assert entity_type(value) == "protein"
    assert class_curie(value) == "biolink:Protein"


@pytest.mark.parametrize("value", [None, 42, "nonsense", "complex", "MI:0326"])
def test_unmodeled_classes_are_rejected(value):
    with pytest.raises(ValueError):
        entity_type(value)


@pytest.mark.parametrize(
    "value", [model.DirectionQualifierEnum.increased, model.DirectionQualifierEnum("increased")]
)
def test_native_enumeration_values(value):
    assert format_term(value) == "increased"
    assert annotation_value(slots.object_direction_qualifier, value) == "increased"


def test_annotation_and_predicate_validation():
    assert annotation_term("MI:2236") == "MI:2236"
    assert annotation_term(slots.object_direction_qualifier) == "object_direction_qualifier"
    for term in ["direction_by_guessing", "biolink:invented"]:
        with pytest.raises(ValueError):
            annotation_term(term)
    with pytest.raises(ValueError):
        annotation_value(slots.object_direction_qualifier, "positive")
    for value in ["activates", "positively_regulates", slots.description, None]:
        with pytest.raises(ValueError):
            predicate(value)


@pytest.mark.parametrize(
    "value", [slots.interacts_with, slots.associated_with, slots.directly_physically_interacts_with]
)
def test_schema_symmetry(value):
    assert is_symmetric(value)
    assert relation_key("a", value, "b") == relation_key("b", value, "a")
    assert not is_symmetric(slots.affects)


def test_qualified_edges_are_distinct_and_qualifier_order_is_irrelevant():
    direction = {"term": "object_direction_qualifier", "value": "increased"}
    aspect = {"term": "object_aspect_qualifier", "value": "activity"}
    negative = dict(direction, value="decreased")
    assert relation_key("a", slots.affects, "b", [direction, aspect]) == relation_key(
        "a", "affects", "b", [aspect, direction]
    )
    assert relation_key("a", "affects", "b", [direction]) != relation_key(
        "a", "affects", "b", [negative]
    )
    assert relation_key("a", "affects", "b", [direction]) != relation_key(
        "a", "affects", "b", [direction, aspect]
    )
    assert relation_key("a", slots.affects, "b") != relation_key("b", slots.affects, "a")
    with pytest.raises(ValueError, match="Conflicting"):
        qualifiers([direction, negative])
    with pytest.raises(ValueError, match="scope"):
        qualifiers([dict(direction, scope="subject")])


def test_legacy_display_projections_are_explicit():
    assert presentation_category(slots.has_part) == "membership"
    assert presentation_category(slots.subclass_of) == "ontology"
    assert presentation_category(slots.associated_with) == "association"
    assert presentation_category(slots.affects) == "interaction"
    assert direction_sign([{"term": "object_direction_qualifier", "value": "upregulated"}]) == 1
    assert direction_sign([{"term": "object_direction_qualifier", "value": "downregulated"}]) == -1
    assert direction_sign([{"term": "MI:2236", "value": None}]) == 0
