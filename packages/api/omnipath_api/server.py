"""FastAPI serving layer for Parquet-only OmniPath.

Routes live at the application root (`/health`, `/entities/search`, …). The
public URL prefix is `/api` via FastAPI `root_path`, matching the Svelte app's
`/api/docs` link and Vite/SvelteKit proxy. `/api/*` and `/app-api/*` are also
accepted as aliases so existing callers keep working.
"""

from __future__ import annotations

import threading
import time
import logging
from contextlib import asynccontextmanager
from urllib.parse import quote
from pathlib import Path
from typing import Any

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from starlette.middleware.gzip import GZipMiddleware
from fastapi.responses import JSONResponse
from starlette.exceptions import HTTPException as StarletteHTTPException

from omnipath_api.admin import AdminService, register_admin_routes
from omnipath_api.settings import Settings
from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.routers.health import router as health_router
from omnipath_api.routers.common import EngineJSONResponse
from omnipath_api.models import (
    ErrorResponse,
)

OPENAPI_TAGS = [
    {"name": "health", "description": "Liveness and engine identity."},
    {"name": "resources", "description": "Indexed Parquet resources and catalogs."},
    {"name": "entities", "description": "Entity search, lookup, and facets."},
    {"name": "relations", "description": "Relation search, evidence, and facets."},
    {"name": "ontology", "description": "Ontology navigation (partial on parquet serving)."},
    {"name": "selection", "description": "Selection-scope expansion for the explorer."},
    {"name": "stats", "description": "Aggregate counts for the resources dashboard."},
    {"name": "export", "description": "Binary slice export (Parquet, Arrow, CSV, JSON)."},
    {"name": "admin", "description": "Inspect resolver/resources and trigger build jobs."},
]


class StripPublicPrefixMiddleware:
    """Map `/api/*` and `/app-api/*` onto FastAPI routes at the application root."""

    def __init__(self, app: Any) -> None:
        self.app = app

    async def __call__(self, scope: dict[str, Any], receive: Any, send: Any) -> None:
        if scope["type"] in {"http", "websocket"}:
            path = scope.get("path") or ""
            for prefix in ("/app-api", "/api"):
                if path == prefix:
                    path = "/"
                    break
                if path.startswith(prefix + "/"):
                    path = path[len(prefix) :] or "/"
                    break
            if path != scope.get("path"):
                scope = dict(scope)
                scope["path"] = path
                raw = scope.get("raw_path")
                prefix_bytes = prefix.encode("ascii")
                scope["raw_path"] = (
                    (raw[len(prefix_bytes) :] or b"/")
                    if raw and raw.startswith(prefix_bytes)
                    else quote(path, safe="/").encode("ascii")
                )
        await self.app(scope, receive, send)


def create_app(
    queries: ParquetQueries | None = None,
    data_root: str | Path | None = None,
    admin_ops: Any | None = None,
    settings: Settings | None = None,
) -> FastAPI:
    settings = settings or Settings()
    serving_queries = queries or ParquetQueries(data_root=data_root or settings.data_root)

    @asynccontextmanager
    async def lifespan(_app):
        if settings.warm_cache:
            from omnipath_api.warm import start

            start(serving_queries)
        yield
        if hasattr(serving_queries, "close"):
            serving_queries.close()

    app = FastAPI(
        title="OmniPath Parquet API",
        description=(
            "REST API for querying OmniPath knowledge graphs directly from nested "
            "Parquet files via DuckDB. The Svelte explorer uses `/api/*` (and "
            "`/app-api/*` as an application alias). Interactive docs: `/api/docs`."
        ),
        version="0.1.0",
        root_path=settings.root_path or "/api",
        lifespan=lifespan,
        docs_url="/docs",
        redoc_url="/redoc",
        openapi_url="/openapi.json",
        default_response_class=EngineJSONResponse,
        openapi_tags=OPENAPI_TAGS,
        responses={500: {"model": ErrorResponse}},
    )
    serving_queries.read_only = settings.read_only
    app.state.settings = settings
    app.state.queries = serving_queries
    app.state.queries_lock = threading.Lock()
    from omnipath_api.status import RequestTimes

    app.state.request_times = RequestTimes()
    app.state.admin = AdminService(data_root=serving_queries.data_root, ops=admin_ops)
    register_admin_routes(app)

    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_methods=["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS"],
        allow_headers=["*"],
    )
    app.add_middleware(GZipMiddleware, minimum_size=1000, compresslevel=5)

    @app.middleware("http")
    async def time_requests(request: Request, call_next):
        started = time.perf_counter()
        response = await call_next(request)
        app.state.request_times.record(request.url.path, time.perf_counter() - started)
        return response

    app.add_middleware(StripPublicPrefixMiddleware)

    @app.exception_handler(Exception)
    async def unhandled_error(_request: Request, exc: Exception) -> JSONResponse:
        if isinstance(exc, StarletteHTTPException):
            detail = exc.detail
            payload = detail if isinstance(detail, dict) else {"error": detail, "detail": detail}
            return EngineJSONResponse(payload, status_code=exc.status_code)
        logging.getLogger(__name__).error("Unhandled request failure", exc_info=exc)
        return EngineJSONResponse({"error": "Internal server error"}, status_code=500)

    app.include_router(health_router)
    from omnipath_api.routers.resources import router as resources_router

    app.include_router(resources_router)
    from omnipath_api.routers.entities import router as entities_router

    app.include_router(entities_router)
    from omnipath_api.routers.relations import router as relations_router

    app.include_router(relations_router)
    from omnipath_api.routers.facets import router as facets_router

    app.include_router(facets_router)
    from omnipath_api.routers.ontology import router as ontology_router

    app.include_router(ontology_router)
    from omnipath_api.routers.selection import router as selection_router

    app.include_router(selection_router)
    from omnipath_api.routers.stats import router as stats_router

    app.include_router(stats_router)
    from omnipath_api.routers.export import router as export_router

    app.include_router(export_router)
    from omnipath_api.routers.index import router as index_router

    app.include_router(index_router)

    return app


def serve(
    host: str = "0.0.0.0",
    port: int = 8085,
    data_root: str | Path | None = None,
) -> None:
    import uvicorn

    resolved_root = Path(data_root).resolve() if data_root else Settings().data_root
    app = create_app(data_root=resolved_root)

    print("\n=======================================================", flush=True)
    print(f"OmniPath Parquet serving: http://{host}:{port}", flush=True)
    print("Engine: DuckDB over the resource Parquet tables", flush=True)
    print(f"Storage: {resolved_root}", flush=True)
    print(f"OpenAPI: http://{host}:{port}/api/docs", flush=True)
    print("Svelte API: /api/* and /app-api/*", flush=True)
    print("=======================================================\n", flush=True)

    uvicorn.run(app, host=host, port=port, log_level="info")


run_server = serve


if __name__ == "__main__":
    import sys

    port = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 8085
    serve(port=port)
