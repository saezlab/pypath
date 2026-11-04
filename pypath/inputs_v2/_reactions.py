"""Builders shared by the modules that emit reactions (see ``omnipath_core.reactions``).

A reaction record states its SBO type, the compartment of each participant
(a GO cellular component, as Biolink's ``anatomical_context_qualifier``) and the
set of compartments it takes place in, which keeps its compartment variants
apart after resolution. Recon3D and Human-GEM share the participant encoding
``base_id:compartment_code:coefficient``.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
import hashlib
import json
from typing import Any

from biolink_model.datamodel.model import Gene, NamedThing, slots
from omnipath_core import reactions as _r
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import (
    SOURCE_REACTION,
    SOURCE_REACTION_TYPE,
)

from pypath.inputs_v2._source_context import compartment_name
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)
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


def _gene_rule_group(
    genes: list[str],
    *,
    identifier_type: str,
    taxon: str,
    gene_namespace: Namespace,
) -> Entity:
    """Source-scoped logical AND clause, without a physical assembly claim.

    Shared by every genome-scale model that emits a gene-protein-reaction
    (GPR) rule. ``identifier_type`` must be unique per model (e.g.
    ``'recon3d.gpr_and'``, ``'imm1415.gpr_and'``) so that two models never
    collide on an identical, unqualified gene-rule digest.

    Args:
        genes: The subunit identifiers of one AND clause.
        identifier_type: The model-scoped identifier type for the group
            entity — not a registered namespace, a model-private digest key.
        taxon: The ``NCBITaxon:...`` value every subunit gene is annotated
            with.
        gene_namespace: The identifier namespace the model's genes use
            (Entrez for Recon3D and iMM1415, Ensembl for Human-GEM).

    Returns:
        A ``NamedThing`` entity whose membership is one ``Gene`` per
        subunit, each in ``gene_namespace`` and annotated with ``taxon``.
    """
    genes = sorted(set(genes))
    digest = hashlib.sha256(json.dumps(genes, separators=(',', ':')).encode()).hexdigest()
    return Entity(
        type=NamedThing,
        identifiers=[Identifier(type=identifier_type, value=digest)],
        membership=[
            Membership(
                member=Entity(
                    type=Gene,
                    identifiers=[Identifier(type=gene_namespace, value=gene)],
                    annotations=[Annotation(term=slots.in_taxon, value=taxon)],
                ),
                predicate=slots.has_member,
            )
            for gene in genes
        ],
    )
