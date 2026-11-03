"""Lossless serving projection of native Biolink quantitative attributes."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from typing import Any

from biolink_model.datamodel.model import QuantityValue
import pyarrow as pa


@dataclass(frozen=True)
class Measurement:
    """A Biolink quantity with source-field and comparator provenance.

    The source field distinguishes e.g. reported minimum/median/maximum without
    inventing a statistical interpretation. The original comparator can retain
    <=/>=, which Biolink 4.4.4's binary relation enum cannot express.
    """

    quantity: QuantityValue
    source_field: str | None = None
    comparator: str | None = None


QUANTITY_STRUCT = pa.struct(
    [
        ("has_numeric_value", pa.float64()),
        ("has_unit", pa.string()),
        ("has_unit_prefix", pa.string()),
        ("has_binary_relation", pa.string()),
        ("source_field", pa.string()),
        ("comparator", pa.string()),
    ]
)

# Non-strict comparisons retain their inclusive source operator while allowing
# the nearest Biolink relation. Approximation must not become exact equality.
COMPARATOR_RELATIONS = {
    "=": {None, "equal_to"},
    "<": {None, "less_than"},
    ">": {None, "greater_than"},
    "<=": {None, "less_than"},
    ">=": {None, "greater_than"},
    "~": {None},
    "≈": {None},
}


def quantity_dict(value: Any) -> dict[str, Any] | None:
    """Project supported values; scalar annotations remain scalar."""
    if isinstance(value, str):
        if not value.startswith('{"has_numeric_value":'):
            return None
        try:
            data = json.loads(value)
        except ValueError:
            return None
        if set(data) != set(QUANTITY_STRUCT.names):
            return None
        return _validate(data)
    measurement = value if isinstance(value, Measurement) else None
    quantity = measurement.quantity if measurement else value
    if not isinstance(quantity, QuantityValue):
        return None
    data = {name: getattr(quantity, name, None) for name in QUANTITY_STRUCT.names}
    for name in ("has_unit", "has_unit_prefix", "has_binary_relation"):
        if data[name] is not None:
            data[name] = str(data[name])
    if measurement:
        data.update(source_field=measurement.source_field, comparator=measurement.comparator)
    return _validate(data)


def _validate(data: dict[str, Any]) -> dict[str, Any]:
    number = data["has_numeric_value"]
    if number is None or isinstance(number, bool) or not math.isfinite(float(number)):
        raise ValueError("Quantitative annotations require a finite numeric value")
    data["has_numeric_value"] = float(number)
    if data["has_binary_relation"] not in (None, "less_than", "equal_to", "greater_than"):
        raise ValueError(
            "Use a native Biolink binary relation; retain source operators in comparator"
        )
    comparator = data["comparator"]
    if comparator is not None:
        if comparator not in COMPARATOR_RELATIONS:
            raise ValueError("Unsupported measurement comparison operator")
        if data["has_binary_relation"] not in COMPARATOR_RELATIONS[comparator]:
            raise ValueError("Source comparator contradicts the native binary relation")
    return data


def quantity_json(value: Any) -> str | None:
    data = quantity_dict(value)
    return (
        None
        if data is None
        else json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    )
