"""Preserve inclusive/approximate source operators without contradictory claims."""

import itertools

from biolink_model.datamodel.model import QuantityValue
import pytest

from omnipath_core.measurements import Measurement, quantity_dict, quantity_json


@pytest.mark.parametrize(
    "comparator,relation",
    [
        ("<=", None),
        ("<=", "less_than"),
        (">=", None),
        (">=", "greater_than"),
        ("~", None),
        ("≈", None),
        ("=", "equal_to"),
    ],
)
def test_comparison_round_trip(comparator, relation):
    measurement = Measurement(
        QuantityValue(has_numeric_value=10, has_binary_relation=relation),
        source_field="IC50",
        comparator=comparator,
    )
    encoded = quantity_json(measurement)
    decoded = quantity_dict(encoded)
    assert decoded["comparator"] == comparator
    assert decoded["has_binary_relation"] == relation
    assert decoded["source_field"] == "IC50"


@pytest.mark.parametrize(
    "comparator,relation",
    [
        ("<=", "greater_than"),
        (">=", "less_than"),
        ("<=", "equal_to"),
        (">=", "equal_to"),
        *itertools.product(("~", "≈"), ("equal_to", "less_than", "greater_than")),
    ],
)
def test_contradictory_comparisons_are_rejected(comparator, relation):
    measurement = Measurement(
        QuantityValue(has_numeric_value=10, has_binary_relation=relation),
        comparator=comparator,
    )
    with pytest.raises(ValueError, match="contradicts"):
        quantity_dict(measurement)
