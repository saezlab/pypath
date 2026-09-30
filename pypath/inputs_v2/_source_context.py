"""Source-event context carried by ordinary membership attributes."""

from collections.abc import Mapping
from typing import Any

from omnipath_core.source_attributes import CONVERSION_DIRECTION

from pypath.internals.tabular_builder import CV


def _direction_assertions(row: Mapping[str, Any]) -> list[list[Any]] | None:
    # One outer item broadcasts all assertions onto every participant. The inner
    # list preserves contradictory fields for downstream evidence diagnostics.
    values = [
        row[field]
        for field in ('direction', 'conversion_direction')
        if row.get(field) not in (None, '')
    ]
    return [values] if values else None


def conversion_direction_cv() -> CV:
    """Retain original direction values without interpreting or qualifying them."""
    return CV(term=CONVERSION_DIRECTION, value=_direction_assertions)
