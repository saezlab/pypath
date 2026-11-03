"""Health check endpoints."""

from __future__ import annotations

import time
from typing import Any
from fastapi import APIRouter, Request

from omnipath_api.models import HealthResponse, ServingStatus

router = APIRouter(tags=["health"])


@router.get("/", response_model=HealthResponse)
def root(request: Request) -> dict[str, Any]:
    return {"status": "ok", "engine": "duckdb-parquet", "timestamp": time.time()}


@router.get("/health", response_model=HealthResponse)
def health(request: Request) -> dict[str, Any]:
    return {"status": "ok", "engine": "duckdb-parquet", "timestamp": time.time()}


@router.get("/status", response_model=ServingStatus)
def serving_status(request: Request) -> dict[str, Any]:
    """How busy the service is: query slots in use and waiting, the median response time
    of the last minute, and the container's CPU load and memory use (fractions)."""
    from omnipath_api.status import status

    pool = getattr(request.app.state.queries, "_query_pool", None)
    return status(pool, request.app.state.request_times)
