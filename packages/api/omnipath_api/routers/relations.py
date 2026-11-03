from __future__ import annotations
from typing import Any
from fastapi import HTTPException, Request, Query
from fastapi.responses import Response

from omnipath_api.models import (
    FacetOption,
    RelationEvidenceResponse,
    RelationFilterOptions,
    RelationPayloadsResponse,
    RelationSearchRequest,
    RelationSearchResponse,
    ScopedRelationFacetsRequest,
)
from fastapi import APIRouter
from omnipath_api.routers.common import _call, _dump, _timed, _parse_json_param

router = APIRouter()


@router.get("/relations/filter-options", response_model=RelationFilterOptions, tags=["relations"])
def relation_filter_options(request: Request) -> dict[str, Any]:
    """Predicate/source/type enumerations for the relation explorer."""
    return _call(request, "get_relation_filter_options")


@router.post("/relations/search", response_model=RelationSearchResponse, tags=["relations"])
def search_relations(
    request: Request,
    response: Response,
    payload: RelationSearchRequest,
) -> dict[str, Any]:
    """Filtered relation search with offset pagination."""
    body = _dump(payload)
    result = _call(
        request,
        "search_relations_api",
        filters=body.get("filters") or {},
        resources=body.get("resources"),
        limit=body.get("limit", 20),
        offset=body.get("offset", 0),
    )
    return _timed(response, result)


@router.post(
    "/relations/scoped-facets",
    response_model=list[FacetOption],
    tags=["relations"],
    response_model_exclude_none=False,
)
def scoped_relation_facets(
    request: Request, payload: ScopedRelationFacetsRequest
) -> list[dict[str, Any]]:
    """The relation filter options the selection and the other groups' filters leave."""
    body = _dump(payload)
    return _call(request, "get_scoped_relation_facets", body, resources=body.get("resources"))


@router.get(
    "/relations/{relation_pk}/evidence",
    response_model=RelationEvidenceResponse,
    tags=["relations"],
)
def relation_evidence(
    request: Request, relation_pk: str, filters: str | None = Query(None)
) -> dict[str, Any]:
    """Nested evidence and annotations for one relation."""
    return _call(
        request, "get_relation_evidence", relation_pk, filters=_parse_json_param(filters, {})
    )


@router.get(
    "/relations/{relation_pk}/payloads",
    response_model=RelationPayloadsResponse,
    tags=["relations"],
)
def relation_payloads(request: Request, relation_pk: str) -> dict[str, Any]:
    """Raw source payloads attached to a relation."""
    return _call(request, "get_relation_payloads", relation_pk)


@router.get("/relations/{relation_pk}", tags=["relations"])
def get_relation(request: Request, relation_pk: str) -> dict[str, Any]:
    """Fetch one relation plus hydrated endpoints."""
    result = _call(request, "get_relation_by_pk", relation_pk)
    if not result:
        raise HTTPException(status_code=404, detail="Relation not found")
    return result
