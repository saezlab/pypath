from __future__ import annotations
from typing import Any
from fastapi import Query, Request

from omnipath_api.models import (
    OntologySearchRequest,
    OntologyTerm,
    OntologyChildrenResponse,
    OntologyTreeRequest,
    TermsRequest,
    TermsResponse,
)
from fastapi import APIRouter
from omnipath_api.routers.common import _call, _dump

router = APIRouter()


@router.get("/ontology/search", tags=["ontology"])
def ontology_search_get(
    request: Request, q: str = "", limit: int = 24, offset: int = 0
) -> list[Any]:
    return _call(request, "search_ontology", {"query": q, "limit": limit, "offset": offset})


@router.post("/ontology/search", tags=["ontology"], response_model=list[OntologyTerm])
@router.post("/ontology/scoped-search", tags=["ontology"])
def ontology_search_post(request: Request, payload: OntologySearchRequest) -> list[Any]:
    body = _dump(payload)
    body["scoped"] = request.url.path.endswith("/scoped-search")
    return _call(request, "search_ontology", body)


@router.get("/ontology/prefixes", tags=["ontology"])
def ontology_prefixes() -> dict[str, list[str]]:
    return {"prefixes": []}


@router.post("/ontology/prefix-counts", tags=["ontology"])
def ontology_prefix_counts(_payload: OntologySearchRequest) -> list[Any]:
    return []


@router.post("/ontology/ontology-id-counts", tags=["ontology"])
def ontology_id_counts(request: Request, payload: OntologySearchRequest) -> list[Any]:
    return _call(request, "search_ontology", _dump(payload), counts=True)


@router.post("/ontology/source-counts", tags=["ontology"])
def ontology_source_counts(_payload: OntologySearchRequest) -> list[Any]:
    return []


@router.get("/ontology/children", tags=["ontology"], response_model=OntologyChildrenResponse)
def ontology_children(
    request: Request,
    termId: str = Query("", description="Term identifier."),
    ontologyId: str | None = Query(None),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
) -> dict[str, Any]:
    """Direct child terms of an ontology term."""
    return _call(
        request,
        "get_ontology_children",
        term_id=termId,
        ontology_id=ontologyId,
        limit=limit,
        offset=offset,
    )


@router.post("/ontology/tree", tags=["ontology"])
def ontology_tree(request: Request, payload: OntologyTreeRequest) -> dict[str, Any]:
    """Ancestor tree for one or more ontology terms."""
    body = _dump(payload)
    term_ids = body.get("term_ids") or body.get("termIds") or []
    ontology_id = body.get("ontologyId") or body.get("ontology_id")
    return _call(request, "get_ontology_tree", term_ids=term_ids, ontology_id=ontology_id)


@router.post("/ontology/entity-ids", tags=["ontology"])
def ontology_entity_ids(_payload: OntologySearchRequest) -> dict[str, list[str]]:
    return {"entityIds": []}


@router.post("/terms", response_model=TermsResponse, tags=["ontology"])
def terms_info(request: Request, payload: TermsRequest) -> dict[str, Any]:
    """Batch lookup of ontology term metadata."""
    body = _dump(payload)
    term_ids = body.get("term_ids") or body.get("termIds") or []
    return _call(request, "get_terms_info", term_ids=term_ids)
