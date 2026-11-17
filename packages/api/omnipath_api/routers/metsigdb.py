"""``/sets/metsigdb``: the metabolite memberships of resource-native sets (MetSigDB).

Read-only over the MetSigDB subset of the served release; the queries are the step files
``metsigdb_rows``, ``metsigdb_groups`` and ``metsigdb_sets``. This module validates the
request, bounds the page and dispatches. A release built without the subset answers 503.
"""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, HTTPException, Query, Request, Response

from omnipath_api import metsigdb as contract
from omnipath_api.routers.common import _call

router = APIRouter(prefix="/sets/metsigdb", tags=["sets"])

# What a grouped response may carry in total. `limit` counts sets, and set sizes span four
# orders of magnitude, so a page of N sets has no bounded row count without this. Sets are
# added whole until the next would cross it; a set larger than the ceiling is returned
# complete and alone, so a published set always stays reachable.
ROW_CEILING = 50_000

# Everything the route accepts. A query parameter the contract does not publish would be
# ignored and answer the wrong question with a 200: `?hmdb=HMDB00077` reads as "the
# memberships of this metabolite" and would return the unfiltered first page. Naming one has
# to fail loudly. `release` selects the release, as everywhere in the API.
MEMBERSHIP_PARAMS = (
    "resource",
    "set_type",
    "set_sub_type",
    "organism",
    "set",
    "entity",
    "metabolite_structure_key",
    "exclude_metabolites",
    "fields",
    "group",
    "limit",
    "offset",
    "total",
    "release",
)
SET_PARAMS = (
    "resource",
    "set_type",
    "set_sub_type",
    "organism",
    "set",
    "limit",
    "offset",
    "total",
    "release",
)

MEMBERSHIPS_DESCRIPTION = """Membership rows: one metabolite in one set of one resource.

Every response is paged. No filter combination returns the whole subset, and the row order
is stable, so paging through a result neither repeats nor skips a membership.

``has_more`` says whether another page follows. ``total`` is the size of the whole result
and is computed only when the request asks for it.

**Fields.** Thirteen by default. ``fields=`` adds any of the others by name, and
``fields=all`` returns every one. An unknown name is refused rather than ignored.

**Grouping.** ``group=true`` groups by resource and then by set. Under grouping ``limit``,
``offset``, ``has_more`` and ``total`` count **sets**, not rows, so no set straddles a page.
A set is never truncated: sets are added whole until the next would cross a 50,000-row
ceiling, and a set larger than the ceiling is returned complete and alone.

**Exclusion.** ``exclude_metabolites`` leaves out metabolites that mass spectrometry does
not measure: ``none`` (the default), ``all``, or a list of ``ions``, ``atoms``,
``small_molecules``, ``xenobiotics``, ``cofactors`` and ``cofactors_extended``. The last is
the cofactors that ChEBI tags with the cofactor role, beyond the reviewed ``cofactors``;
``all`` includes it. The memberships stay in the subset; the request leaves them out. It
applies before paging and counting, so ``total`` and ``has_more`` describe the filtered
result. Under grouping a set that loses every member is still returned, with an empty
``members`` list, ``returned`` of zero and its published ``set_size``.

**Identifiers.** ``entity`` accepts an entity key, an InChIKey, or an HMDB, ChEBI, KEGG or
PubChem identifier, recognised from the value itself. A bare integer matches ChEBI **or**
PubChem and returns the union, because both are bare integers. ``metabolite_structure_key``
selects the metabolites that share a connectivity-level structure (the charge, stereo and
tautomer variants of one molecule). SMILES is published and reachable through ``fields``,
but is not searchable.

Nothing is translated on the request path. Values are matched against the identifier columns
the subset publishes, and the subset publishes what each resource and the release's other
resources state about the metabolite, so coverage depends on the resource.
"""

SETS_DESCRIPTION = """The sets of the MetSigDB subset without their members: resource, id,
label, type, sub-type, organism and size, one row per set, in the order of the key. Takes the
set-level filters of ``/sets/metsigdb``."""


def coverage_table(coverage: dict[str, dict]) -> str:
    """The identifier coverage of each resource, measured when the subset was built."""
    columns = ("hmdb", "chebi", "kegg", "pubchem", "inchikey")
    lines = [
        "| resource | metabolites | " + " | ".join(columns) + " |",
        "|---|---|" + "|".join("---" for _ in columns) + "|",
    ]
    for resource in contract.RESOURCES:
        found = coverage.get(resource)
        if not found:
            continue
        shares = found.get("identifiers", {})
        lines.append(
            f"| {resource} | {found.get('metabolites', 0):,} | "
            + " | ".join(f"{round(100 * shares.get(name, 0))}%" for name in columns)
            + " |"
        )
    return "\n".join(lines)


def describe_openapi(app) -> None:
    """Publish the coverage of the release being served in the route's description.

    The table is measured when the subset is built and read from the served release each
    time the document is asked for, so it cannot go stale: a request for an identifier a
    resource's metabolites lack answers empty, and the table says which those are.
    """
    original = app.openapi

    def openapi() -> dict[str, Any]:
        schema = original()
        try:
            status = app.state.queries.subset_status()
            coverage = status["metsigdb"]["stats"]["coverage"]
            text = (
                "**Coverage** of the release being served: the share of each resource's "
                "metabolites that carry each identifier, measured when the subset was built. "
                "A query for an identifier a resource's metabolites lack returns nothing for "
                "that resource.\n\n" + coverage_table(coverage)
            )
        except Exception:
            text = "**Coverage** is not available: the served release has no MetSigDB subset."
        operation = schema.get("paths", {}).get("/sets/metsigdb", {}).get("get")
        if operation is not None:
            operation["description"] = MEMBERSHIPS_DESCRIPTION + "\n" + text
        return schema

    app.openapi = openapi


def _listed(values: list[str] | None) -> list[str]:
    """A repeated or comma-separated parameter as one list."""
    return [part.strip() for value in values or () for part in value.split(",") if part.strip()]


def _reject_unknown(request: Request, accepted: tuple[str, ...]) -> None:
    unknown = sorted(set(request.query_params) - set(accepted))
    if unknown:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Unsupported query parameter: {', '.join(unknown)}. "
                f"Supported: {', '.join(p for p in accepted if p != 'release')}."
            ),
        )


def _refuse(error: ValueError) -> HTTPException:
    return HTTPException(status_code=400, detail=str(error))


def _set_filters(
    resource: list[str] | None,
    set_type: list[str] | None,
    set_sub_type: list[str] | None,
    organism: int | None,
    set_: list[str] | None,
) -> dict[str, Any]:
    try:
        return {
            "resource": contract.canonical_values(_listed(resource), contract.RESOURCES, "resource")
            or None,
            "set_type": contract.canonical_values(_listed(set_type), contract.SET_TYPES, "set_type")
            or None,
            "set_sub_type": contract.canonical_values(
                _listed(set_sub_type), contract.SET_SUB_TYPES, "set_sub_type"
            )
            or None,
            "organism": organism,
            "sets": _listed(set_) or None,
        }
    except ValueError as exc:
        raise _refuse(exc) from exc


def _answer(request: Request, name: str, params: dict[str, Any], *, also: tuple[str, ...] = ()):
    text = _call(
        request,
        "run_subset_query",
        "metsigdb",
        name,
        params,
        row_contract=contract.ROW_CONTRACT,
        also=also,
    )
    return Response(content=text if text is not None else "null", media_type="application/json")


@router.get("", description=MEMBERSHIPS_DESCRIPTION, summary="MetSigDB membership rows")
@router.get("/", include_in_schema=False)
def memberships(
    request: Request,
    resource: list[str] | None = Query(default=None),
    set_type: list[str] | None = Query(default=None),
    set_sub_type: list[str] | None = Query(default=None),
    organism: int | None = Query(default=None),
    set: list[str] | None = Query(default=None),
    entity: list[str] | None = Query(default=None),
    metabolite_structure_key: list[str] | None = Query(default=None),
    exclude_metabolites: list[str] | None = Query(default=None),
    fields: list[str] | None = Query(default=None),
    group: bool = Query(default=False),
    limit: int = Query(default=contract.DEFAULT_LIMIT, ge=1, le=contract.MAX_LIMIT),
    offset: int = Query(default=0, ge=0),
    total: bool = Query(default=False),
) -> Response:
    _reject_unknown(request, MEMBERSHIP_PARAMS)
    params: dict[str, Any] = _set_filters(resource, set_type, set_sub_type, organism, set)
    try:
        shown = contract.resolve_fields(_listed(fields))
        excluded = contract.exclusion_classes(_listed(exclude_metabolites))
        identifiers = contract.entity_parameters(_listed(entity))
    except ValueError as exc:
        raise _refuse(exc) from exc
    params |= {
        "entity_given": bool(_listed(entity)),
        **{name: values or None for name, values in identifiers.items()},
        "structure_keys": _listed(metabolite_structure_key) or None,
        "exclude_classes": excluded or None,
        "limit": limit,
        "offset": offset,
        "total": total,
    }
    also = ("exclusion_metabolites",) if excluded else ()
    if group:
        params |= {
            # The resource is the group's own key, so a set does not repeat it.
            "set_fields": [f for f in shown if f in contract.SET_FIELDS and f != "resource"],
            "member_fields": [f for f in shown if f in contract.MEMBER_FIELDS],
            "row_ceiling": ROW_CEILING,
        }
        return _answer(request, "metsigdb_groups", params, also=also)
    params["fields"] = list(shown)
    return _answer(request, "metsigdb_rows", params, also=also)


@router.get("/sets", description=SETS_DESCRIPTION, summary="MetSigDB sets")
def sets(
    request: Request,
    resource: list[str] | None = Query(default=None),
    set_type: list[str] | None = Query(default=None),
    set_sub_type: list[str] | None = Query(default=None),
    organism: int | None = Query(default=None),
    set: list[str] | None = Query(default=None),
    limit: int = Query(default=contract.DEFAULT_LIMIT, ge=1, le=contract.MAX_LIMIT),
    offset: int = Query(default=0, ge=0),
    total: bool = Query(default=False),
) -> Response:
    _reject_unknown(request, SET_PARAMS)
    params = _set_filters(resource, set_type, set_sub_type, organism, set)
    params |= {"limit": limit, "offset": offset, "total": total}
    return _answer(request, "metsigdb_sets", params)
