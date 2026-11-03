from __future__ import annotations
from typing import Any
from fastapi import Request

from omnipath_api.models import (
    BuildManifest,
    StatsEntityType,
    StatsInteractionType,
    StatsSource,
)
from fastapi import APIRouter
from omnipath_api.routers.common import _call

router = APIRouter()


@router.get("/stats/sources", response_model=list[StatsSource], tags=["stats"])
def stats_sources(request: Request) -> list[dict[str, Any]]:
    return _call(request, "get_stats_sources")


@router.get("/stats/entity-types", response_model=list[StatsEntityType], tags=["stats"])
def stats_entity_types(request: Request) -> list[dict[str, Any]]:
    return _call(request, "get_stats_entity_types")


@router.get("/stats/interaction-types", response_model=list[StatsInteractionType], tags=["stats"])
def stats_interaction_types(request: Request) -> list[dict[str, Any]]:
    return _call(request, "get_stats_interaction_types")


@router.get("/stats/build-manifest", response_model=BuildManifest, tags=["stats"])
def stats_build_manifest(request: Request) -> dict[str, Any]:
    return _call(request, "get_stats_build_manifest")
