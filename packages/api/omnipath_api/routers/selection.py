from __future__ import annotations
from typing import Any
from fastapi import Request

from omnipath_api.models import (
    SelectionScopeRequest,
    SelectionScopeResponse,
)
from fastapi import APIRouter
from omnipath_api.routers.common import _call, _dump

router = APIRouter()


@router.post("/selection/scope", response_model=SelectionScopeResponse, tags=["selection"])
def selection_scope(request: Request, payload: SelectionScopeRequest) -> dict[str, Any]:
    """Expand the current explorer selection to associated entities."""
    return _call(request, "resolve_selection_scope", _dump(payload))
