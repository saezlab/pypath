from __future__ import annotations
from typing import Any
from fastapi import Request
from fastapi.responses import Response

from omnipath_api.models import (
    ScopedRelationFacetsRequest,
)
from fastapi import APIRouter
from omnipath_api.routers.common import _call, _dump, _timed

router = APIRouter()


@router.post(
    "/facets",
    response_model=None,
    tags=["relations"],
    summary="Relation facet counts",
)
def relation_facets(
    request: Request,
    response: Response,
    payload: ScopedRelationFacetsRequest,
) -> Any:
    """Facet payload for filtered relation searches."""
    body = _dump(payload)
    if body.get("filters"):
        result = _call(
            request,
            "get_facets",
            filters=body.get("filters") or {},
            resources=body.get("resources"),
        )
        svelte = {
            "categoryCounts": result.get("categories") or {},
            "predicateCounts": result.get("predicates") or {},
            "signCounts": result.get("signs") or {},
            "elapsed_ms": result.get("elapsed_ms"),
        }
        return _timed(response, {**result, **svelte})
    return _call(request, "get_scoped_relation_facets", body, resources=body.get("resources"))
