"""``/sets/exclusion-metabolites``: the table of metabolites that mass spectrometry does not
measure, as the exclusion subset of the served release holds it.

The same table drives ``exclude_metabolites`` on ``/sets/metsigdb``; this route returns it
itself. Read-only; the queries are the step files ``exclusion_metabolites`` and
``exclusion_metabolites_report``. A release built without the subset answers 503.
"""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Query, Request, Response

from omnipath_api import metsigdb as contract
from omnipath_api.routers.common import _call
from omnipath_api.routers.metsigdb import _listed, _refuse, _reject_unknown

router = APIRouter(prefix="/sets/exclusion-metabolites", tags=["sets"])

# The row contract of the subset this route reads (``subset_layout.ROW_CONTRACTS``).
ROW_CONTRACT = 1

# The vocabularies of the table's columns (the CHECK constraints of its definition).
ASSIGNMENT_SOURCES = (
    "metaproviz_seed",
    "project_seed",
    "override_add",
    "chebi_hierarchy",
    "chebi_role",
)
ORIGINS = ("curated", "derived")
MATCH_KINDS = ("entity", "structure_key")
REPORT_KINDS = (
    "seed_coverage",
    "class_size",
    "seed_not_derived",
    "derived_not_seeded",
    "shared_key",
    "no_key",
    "veto_applied",
    "uncommitted_inputs",
    "cofactor_reactions",
    "cofactor_cross_check",
)

# Everything a row may carry, in the order of the macro `op_exclusion_fields`.
ALL_FIELDS = (
    "entity",
    "metabolite_label",
    "metabolite_entity_type",
    "inchikey",
    "hmdb",
    "pubchem",
    "chebi",
    "kegg",
    "exclusion_class",
    "assignment_source",
    "origin",
    "match_kind",
    "seed_namespace",
    "seed_identifier",
    "evidence",
    "via_entity",
    "list_version",
    "build_id",
)
DEFAULT_FIELDS = (
    "entity",
    "metabolite_label",
    "inchikey",
    "hmdb",
    "pubchem",
    "chebi",
    "kegg",
    "exclusion_class",
    "assignment_source",
    "origin",
    "match_kind",
)

# A parameter the route does not publish would be ignored and answer another question with
# a 200, so naming one fails loudly (as on ``/sets/metsigdb``).
ROW_PARAMS = (
    "exclusion_class",
    "assignment_source",
    "origin",
    "match_kind",
    "entity",
    "fields",
    "limit",
    "offset",
    "total",
    "release",
)
REPORT_PARAMS = ("kind", "exclusion_class", "limit", "offset", "total", "release")

ROWS_DESCRIPTION = """The exclusion metabolite table: metabolites that mass spectrometry does
not measure, each in a class, with the assignment that put it there.

Classes are ``ions``, ``atoms``, ``small_molecules``, ``xenobiotics``, ``cofactors`` and
``cofactors_extended`` (the cofactors ChEBI tags with the cofactor role, beyond the reviewed
``cofactors``). A metabolite in one class through two sources has two rows. ``origin`` is
``curated`` (a seed list) or ``derived`` (a ChEBI closure). ``match_kind`` is ``entity`` for
a listed metabolite and ``structure_key`` for a charge, tautomer or stereo variant that
shares its structure key; ``via_entity`` names the listed one. Variants stay separate
entities.

Every response is paged and the order is stable. ``fields=`` adds any of the other fields to
the default, ``fields=all`` returns every one. ``entity`` accepts an entity key, an InChIKey,
or an HMDB, ChEBI, KEGG or PubChem identifier; a bare integer matches ChEBI or PubChem.

The same table filters ``/sets/metsigdb`` through ``exclude_metabolites``."""

REPORT_DESCRIPTION = """What the exclusion build measured: the size of each class and where its
members came from, how many seeds resolved, and where the curated and derived lists disagree.
One row per finding; ``detail`` is JSON."""


def _answer(request: Request, name: str, params: dict[str, Any]) -> Response:
    text = _call(
        request,
        "run_subset_query",
        "exclusion_metabolites",
        name,
        params,
        row_contract=ROW_CONTRACT,
    )
    return Response(content=text if text is not None else "null", media_type="application/json")


def _fields(requested: list[str] | None) -> tuple[str, ...]:
    names = _listed(requested)
    if not names:
        return DEFAULT_FIELDS
    if contract.ALL_KEYWORD in names:
        return ALL_FIELDS
    unknown = sorted({name for name in names if name not in ALL_FIELDS})
    if unknown:
        raise ValueError(
            f"Unsupported field: {', '.join(unknown)}. "
            f'Supported: {", ".join(ALL_FIELDS)}, or "{contract.ALL_KEYWORD}".'
        )
    wanted = set(DEFAULT_FIELDS) | set(names)
    return tuple(field for field in ALL_FIELDS if field in wanted)


@router.get("", description=ROWS_DESCRIPTION, summary="Exclusion metabolite table")
@router.get("/", include_in_schema=False)
def rows(
    request: Request,
    exclusion_class: list[str] | None = Query(default=None),
    assignment_source: list[str] | None = Query(default=None),
    origin: list[str] | None = Query(default=None),
    match_kind: list[str] | None = Query(default=None),
    entity: list[str] | None = Query(default=None),
    fields: list[str] | None = Query(default=None),
    limit: int = Query(default=contract.DEFAULT_LIMIT, ge=1, le=contract.MAX_LIMIT),
    offset: int = Query(default=0, ge=0),
    total: bool = Query(default=False),
) -> Response:
    _reject_unknown(request, ROW_PARAMS)
    try:
        params: dict[str, Any] = {
            "classes": contract.canonical_values(
                _listed(exclusion_class), contract.EXCLUSION_CLASSES, "exclusion_class"
            )
            or None,
            "sources": contract.canonical_values(
                _listed(assignment_source), ASSIGNMENT_SOURCES, "assignment_source"
            )
            or None,
            "origins": contract.canonical_values(_listed(origin), ORIGINS, "origin") or None,
            "match_kinds": contract.canonical_values(
                _listed(match_kind), MATCH_KINDS, "match_kind"
            )
            or None,
            "fields": list(_fields(fields)),
        }
        identifiers = contract.entity_parameters(_listed(entity))
    except ValueError as exc:
        raise _refuse(exc) from exc
    params |= {
        "entity_given": bool(_listed(entity)),
        **{name: values or None for name, values in identifiers.items()},
        "limit": limit,
        "offset": offset,
        "total": total,
    }
    return _answer(request, "exclusion_metabolites", params)


@router.get("/report", description=REPORT_DESCRIPTION, summary="Exclusion build report")
def report(
    request: Request,
    kind: list[str] | None = Query(default=None),
    exclusion_class: list[str] | None = Query(default=None),
    limit: int = Query(default=contract.DEFAULT_LIMIT, ge=1, le=contract.MAX_LIMIT),
    offset: int = Query(default=0, ge=0),
    total: bool = Query(default=False),
) -> Response:
    _reject_unknown(request, REPORT_PARAMS)
    try:
        params = {
            "kinds": contract.canonical_values(_listed(kind), REPORT_KINDS, "kind") or None,
            "classes": contract.canonical_values(
                _listed(exclusion_class), contract.EXCLUSION_CLASSES, "exclusion_class"
            )
            or None,
        }
    except ValueError as exc:
        raise _refuse(exc) from exc
    params |= {"limit": limit, "offset": offset, "total": total}
    return _answer(request, "exclusion_metabolites_report", params)
