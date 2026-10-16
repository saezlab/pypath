"""Builders shared by the modules that emit reactions (see ``omnipath_core.reactions``).

A reaction record states its SBO type, the compartment of each participant
(a GO cellular component, as Biolink's ``anatomical_context_qualifier``) and the
set of compartments it takes place in, which keeps its compartment variants
apart after resolution. Recon3D and Human-GEM share the participant encoding
``base_id:compartment_code:coefficient``.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from typing import Any

from omnipath_core import reactions as _r
from omnipath_core.source_attributes import SOURCE_REACTION, SOURCE_REACTION_TYPE

from pypath.inputs_v2._source_context import compartment_name
from pypath.internals.tabular_builder import CV

Row = Mapping[str, Any]


def model_compartment_go(code: str | None) -> str | None:
    """GO cellular component of a Recon3D / Human-GEM compartment code."""
    return _r.go_compartment(compartment_name(code))


def stoichiometry_compartment_go(value: str | None) -> str | None:
    """GO compartment of one ``base_id:compartment:coefficient`` entry."""
    parts = (value or '').split(':')
    return model_compartment_go(parts[1]) if len(parts) > 1 else None


def model_reaction_compartments(row: Row) -> list[str]:
    """The GO compartments of a model reaction's participants, once each."""
    found = []
    for field in ('reactants', 'products'):
        for entry in str(row.get(field) or '').split('||'):
            go = stoichiometry_compartment_go(entry)
            if go and go not in found:
                found.append(go)
    return found


def reaction_type_cv(value: Callable[[Row], Any] | Any) -> CV:
    """The reaction's SBO terms, in the Biolink ``type`` slot."""
    return CV(term=_r.REACTION_TYPE, value=value)


def source_reaction_type_cv(value: Callable[[Row], Any] | Any) -> CV:
    """The source's own name for the reaction's type."""
    return CV(term=SOURCE_REACTION_TYPE, value=value)


def reaction_compartments_cv(value: Callable[[Row], Iterable[str]]) -> CV:
    """The compartments the reaction takes place in, on the reaction entity."""
    return CV(term=_r.COMPARTMENT, value=lambda row: list(value(row)) or None)


def source_reaction_cv(value: Callable[[Row], Any] | Any) -> CV:
    """The source's id of the reaction, where catalysis is stated apart from it."""
    return CV(term=SOURCE_REACTION, value=value)


def model_reaction_type(
    row: Row,
    *,
    transport: bool,
    pseudoreaction: bool,
) -> str:
    """SBO type of a genome-scale model reaction the model types by subsystem.

    A pseudoreaction (one side empty) is an exchange when it opens the
    extracellular compartment, a sink when it can run both ways and a demand
    otherwise, as in the COBRA convention Recon3D's own SBO terms follow.
    """
    if transport:
        return _r.TRANSPORT
    if not pseudoreaction:
        return _r.BIOCHEMICAL_REACTION
    compartments = model_reaction_compartments(row)
    if compartments == [_r.GO_COMPARTMENTS['extracellular space']]:
        return _r.EXCHANGE
    lower, upper = (_number(row.get(f'{bound}_bound')) for bound in ('lower', 'upper'))
    if lower is not None and upper is not None and lower < 0 < upper:
        return _r.SINK
    return _r.DEMAND


def _number(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
