from __future__ import annotations
from typing import Any, Literal
from fastapi import HTTPException, Query, Request
from fastapi.responses import Response

from omnipath_api.models import (
    DetailField,
    EntityGroupsRequest,
    GroupRelationshipsRequest,
    EntityDetailsResponse,
    EntityExamplesResponse,
    EntityFilterOptions,
    EntityPksRequest,
    EntityPublicIdsRequest,
    EntityResolveRequest,
    EntitySearchRequest,
    EntitySearchResponse,
    EntitiesResponse,
    FacetCount,
    ScopedEntityFacetsRequest,
)
from fastapi import APIRouter
from omnipath_api.routers.common import _call, _cursor_dict, _dump, _parse_json_param, _timed

router = APIRouter()


@router.get("/entities/search", response_model=EntitySearchResponse, tags=["entities"])
def search_entities_get(
    request: Request,
    response: Response,
    q: str = Query("", description="Identifier or label substring."),
    limit: int = Query(20, ge=1, le=1000),
    cursor: str | None = Query(None, description="JSON object `{relationCount, entityPk}`."),
    filters: str | None = Query(None, description="JSON object of search filters."),
) -> dict[str, Any]:
    """Autocomplete-style entity search used by GET callers."""
    parsed_filters = _parse_json_param(filters, {})
    parsed_cursor = _parse_json_param(cursor, None)
    result = _call(
        request,
        "search_entities_api",
        query=q,
        filters=parsed_filters,
        limit=limit,
        cursor=parsed_cursor,
    )
    result["rows"] = [
        {
            "entity_key": entity["entityPk"],
            "label": entity.get("label") or entity["canonicalIdentifier"],
            "entity_type": entity.get("entityType"),
            "namespace": entity.get("canonicalIdentifierType"),
            "identifier": entity.get("canonicalIdentifier"),
            "taxon": entity.get("taxonomyId"),
        }
        for entity in result.get("entities") or []
    ]
    result["count"] = len(result["rows"])
    return _timed(response, result)


@router.post("/entities/search", response_model=EntitySearchResponse, tags=["entities"])
def search_entities_post(
    request: Request,
    response: Response,
    payload: EntitySearchRequest,
) -> dict[str, Any]:
    """Entity search with JSON filters, cursor pagination, and optional resource scope."""
    body = _dump(payload)
    query = body.get("query") or body.get("q") or ""
    filters = body.get("filters") or {}
    if body.get("entity_types"):
        filters = {**filters, "entity_types": body.get("entity_types")}
    if body.get("shape") == "web":
        result = _call(
            request,
            "search_entities",
            query=query,
            resources=body.get("resources"),
            entity_types=body.get("entity_types"),
            taxon=body.get("taxon"),
            limit=body.get("limit", 20),
        )
    else:
        result = _call(
            request,
            "search_entities_api",
            query=query,
            filters=filters,
            resources=body.get("resources"),
            limit=body.get("limit", 20),
            cursor=_cursor_dict(body.get("cursor")),
        )
    return _timed(response, result)


@router.post("/entities/connectivity-groups", tags=["entities"], include_in_schema=False)
@router.post("/entities/groups", tags=["entities"])
def connectivity_groups(request: Request, payload: EntityGroupsRequest):
    """Group matching entities using an explicit grouping strategy.

    Counts and member pagination are scoped to the current query and filters.
    Auto groups chemicals by the first 14 InChIKey characters and biological
    entities by a shared gene reference, in one ordered page. Missing or
    conflicting references remain separate singleton results. Explicit gene
    reference and chemical connectivity strategies remain available.
    """
    return _call(request, "search_entity_groups", **_dump(payload))


@router.get("/entities/filter-options", response_model=EntityFilterOptions, tags=["entities"])
def entity_filter_options(request: Request) -> dict[str, Any]:
    """Distinct entity types, sources, and taxa for explorer filters."""
    return _call(request, "get_entity_filter_options")


@router.post("/entities/by-pks", response_model=EntitiesResponse, tags=["entities"])
def entities_by_pks(request: Request, payload: EntityPksRequest) -> dict[str, Any]:
    """Hydrate entities by primary key."""
    body = _dump(payload)
    pks = body.get("pks") or body.get("entityPks") or []
    return _call(request, "get_entities_by_pks", [str(v) for v in pks])


@router.post("/entities/by-public-ids", response_model=EntitiesResponse, tags=["entities"])
def entities_by_public_ids(request: Request, payload: EntityPublicIdsRequest) -> dict[str, Any]:
    """Hydrate entities by public identifier."""
    body = _dump(payload)
    public_ids = body.get("public_ids") or body.get("publicIds") or []
    return _call(request, "get_entities_by_public_ids", [str(v) for v in public_ids])


@router.post("/entities/resolve", response_model=EntitiesResponse, tags=["entities"])
def resolve_entities(request: Request, payload: EntityResolveRequest) -> dict[str, Any]:
    """Resolve raw identifiers to canonical entities."""
    return _call(request, "resolve_identifiers", [str(v) for v in payload.identifiers])


@router.post("/entities/scoped-facets", response_model=list[FacetCount], tags=["entities"])
def scoped_entity_facets(
    request: Request, payload: ScopedEntityFacetsRequest
) -> list[dict[str, Any]]:
    """Facet counts for the current entity selection."""
    return _call(request, "get_scoped_entity_facets", _dump(payload))


@router.get("/entities/examples", response_model=EntityExamplesResponse, tags=["entities"])
def entity_examples(request: Request) -> dict[str, Any]:
    """Curated landing suggestions from the selected release; not a search page."""
    return _call(request, "get_entity_examples")


@router.post("/entities/group-relationships", tags=["entities"])
def group_relationships(request: Request, payload: GroupRelationshipsRequest):
    keys = payload.member_keys
    if payload.group_key:
        groups = _call(
            request,
            "search_entity_groups",
            strategy=payload.strategy,
            group_key=payload.group_key,
            query=payload.query,
            filters=_dump(payload)["filters"],
            resources=payload.resources,
            include_member_keys=True,
        )["groups"]
        keys = groups[0]["entity"]["groupMemberKeys"] if groups else []
    elif not keys:
        raise HTTPException(status_code=422, detail="Provide a group or member keys")
    result = _call(
        request,
        "get_entity_relationships",
        keys,
        resources=payload.resources,
        limit=payload.limit,
        offset=payload.offset,
    )
    member_set = set(keys)
    for row in result["relationships"]:
        row["groupOutgoing"] = row["relation"]["subjectEntityPk"] in member_set
    return result


@router.get("/entities/{entity_id}/molecular-context", tags=["entities"])
def molecular_context(
    request: Request,
    entity_id: str,
    view: Literal["reference", "product"] = "reference",
    isoform_identifier: str | None = None,
    limit: int = Query(20, ge=1, le=100),
    offset: int = Query(0, ge=0),
):
    result = _call(
        request,
        "get_molecular_context",
        entity_id,
        view=view,
        isoform_identifier=isoform_identifier,
        limit=limit,
        offset=offset,
    )
    if result is None:
        raise HTTPException(status_code=404, detail="Entity not found")
    return result


@router.get("/entities/{entity_id}/relationships", tags=["entities"])
def entity_relationships(
    request: Request,
    entity_id: str,
    limit: int = Query(50, ge=1, le=100),
    offset: int = Query(0, ge=0),
):
    if not _call(request, "get_entity_core", entity_id):
        raise HTTPException(status_code=404, detail="Entity not found")
    return _call(request, "get_entity_relationships", entity_id, limit=limit, offset=offset)


@router.get("/entities/{entity_id}/evidence", tags=["entities"])
def entity_evidence(
    request: Request,
    entity_id: str,
    limit: int = Query(20, ge=1, le=100),
    offset: int = Query(0, ge=0),
):
    """A page of the entity's molecular evidence: where and how sources state it."""
    result = _call(request, "get_entity_evidence", entity_id, limit=limit, offset=offset)
    if result is None:
        raise HTTPException(status_code=404, detail="Entity not found")
    return result


@router.get("/entities/{entity_id}", response_model=EntityDetailsResponse, tags=["entities"])
def get_entity(
    request: Request,
    entity_id: str,
    includeRelationships: bool = True,
    detail_limit: int = Query(20, ge=1, le=500),
    detail_offset: int = Query(0, ge=0),
    detail_field: DetailField | None = None,
) -> dict[str, Any]:
    """Fetch one entity by its exact entity key. No identifier or text fallback.

    ``detail_field`` pages one collection (identifiers or attributes) on its own.
    """
    result = _call(
        request, "get_entity_details" if includeRelationships else "get_entity_core", entity_id
    )
    if not result:
        raise HTTPException(status_code=404, detail="Entity not found")
    from omnipath_api.queries.entity_details import page_details

    result["entity"] = page_details(result["entity"], detail_limit, detail_offset, detail_field)
    return result
