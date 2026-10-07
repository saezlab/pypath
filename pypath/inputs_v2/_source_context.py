"""Source-event context carried by ordinary membership attributes."""

from collections.abc import Mapping
from typing import Any

from omnipath_core.source_attributes import CONVERSION_DIRECTION

from pypath.internals.tabular_builder import CV

# Compartment codes of Recon3D and Human-GEM, which share them, as GO cellular
# component names. [i] is the mitochondrial intermembrane space (Recon3D/BiGG:
# "inner mitochondrial compartment", Human-GEM: "Inner mitochondria"); [m] is
# the mitochondrion.
_COMPARTMENTS = {
    'c': 'cytosol',
    'e': 'extracellular space',
    'g': 'Golgi apparatus',
    'i': 'mitochondrial intermembrane space',
    'l': 'lysosome',
    'm': 'mitochondrion',
    'n': 'nucleus',
    'r': 'endoplasmic reticulum',
    'x': 'peroxisome',
}
# BioPAX conversionDirection values, by their spellings in the sources.
_CONVERSION_DIRECTIONS = {
    'left-to-right': 'LEFT-TO-RIGHT',
    'left_to_right': 'LEFT-TO-RIGHT',
    'right-to-left': 'RIGHT-TO-LEFT',
    'right_to_left': 'RIGHT-TO-LEFT',
    'reversible': 'REVERSIBLE',
}


def compartment_name(code: str | None) -> str | None:
    """GO cellular component name of a model compartment code; unknown codes are kept."""
    return _COMPARTMENTS.get(code, code) or None


def conversion_direction(value: Any) -> Any:
    """BioPAX spelling of a conversion direction; unknown values are kept."""
    return _CONVERSION_DIRECTIONS.get(str(value).strip().lower(), value)


def flux_direction(lower_bound: float, upper_bound: float) -> str:
    """BioPAX conversion direction from a model reaction's flux bounds."""
    if lower_bound < 0 < upper_bound:
        return 'REVERSIBLE'
    return 'RIGHT-TO-LEFT' if lower_bound < 0 else 'LEFT-TO-RIGHT'


def stoichiometric_coefficient(value: Any) -> str | None:
    """Coefficient as text, integral values without decimals; symbols are kept."""
    text = str(value if value is not None else '').strip()
    try:
        number = float(text)
    except ValueError:
        return text or None
    return str(int(number)) if number.is_integer() else text


def _direction_assertions(row: Mapping[str, Any]) -> list[list[Any]] | None:
    # One outer item broadcasts all assertions onto every participant. The inner
    # list preserves contradictory fields for downstream evidence diagnostics.
    values = [
        conversion_direction(row[field])
        for field in ('direction', 'conversion_direction')
        if row.get(field) not in (None, '')
    ]
    return [values] if values else None


def conversion_direction_cv() -> CV:
    """Retain reported direction values in BioPAX spelling, without qualifying them."""
    return CV(term=CONVERSION_DIRECTION, value=_direction_assertions)
