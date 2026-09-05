"""Source-neutral parsing of finite numeric observations and comparison bounds."""

from __future__ import annotations

import math
import re

from biolink_model.datamodel.model import QuantityValue
from omnipath_core.measurements import Measurement

_NUMBER = re.compile(
    r'\s*(<=|>=|<|>|=|~|≈)?\s*'
    r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*'
)
_COMPARATORS = {None, '<', '>', '=', '<=', '>=', '~', '≈'}


def measurement(
    value: object,
    unit: str | None = None,
    source_field: str | None = None,
    comparator: str | None = None,
) -> Measurement | None:
    """Parse a scalar without inferring its endpoint, units or biological meaning.

    Invalid, non-finite and contradictory values omit only the measurement.
    Preserve inclusive/approximate bounds without replacing them by strict ones.
    """
    if value is None:
        return None
    match = _NUMBER.fullmatch(str(value))
    if match is None:
        return None
    explicit = (
        (str(comparator).strip() or None) if comparator is not None else None
    )
    embedded = match.group(1)
    if explicit not in _COMPARATORS:
        return None
    if explicit and embedded and explicit != embedded:
        return None
    comparison = explicit or embedded
    number = float(match.group(2))
    if not math.isfinite(number):
        return None
    return Measurement(
        quantity=QuantityValue(
            has_numeric_value=number,
            has_unit=unit or None,
            has_binary_relation={
                '<': 'less_than',
                '>': 'greater_than',
                '=': 'equal_to',
            }.get(comparison),
        ),
        source_field=source_field,
        comparator=comparison,
    )
