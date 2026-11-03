"""Authenticated administrative HTTP routes."""

from __future__ import annotations

from contextlib import nullcontext
import hmac
import json
import re
from collections.abc import Iterator
from typing import Any
from fastapi import HTTPException, Request
from fastapi.responses import StreamingResponse
from omnipath_api.admin import AdminService
from omnipath_api.jobs.build_ops import BuildUnavailable
from omnipath_api.models import AdminJobRequest, AdminJobResponse, AdminJobsResponse
from omnipath_api.queries.store.inventory import ReleaseStore
from omnipath_api.queries.store.query_pool import QueryCapacityError


def require_admin(request: Request) -> None:
    config = request.app.state.settings
    secret = config.admin_secret
    if not secret:
        raise HTTPException(
            status_code=503, detail="Administration is disabled: configure an admin secret"
        )
    provided = request.headers.get("x-admin-secret") or ""
    if not provided or not hmac.compare_digest(provided.encode(), secret.encode()):
        raise HTTPException(status_code=401, detail="Admin password required")
    if request.method not in {"GET", "HEAD", "OPTIONS"} and config.read_only:
        raise HTTPException(
            status_code=403, detail="Administration mutations are disabled in read-only mode"
        )


def register_admin_routes(app: Any) -> None:
    def _admin(request: Request) -> AdminService:
        require_admin(request)
        return request.app.state.admin

    @app.get("/admin/status", tags=["admin"])
    def admin_status(request: Request) -> dict[str, Any]:
        """Resolver files, built resources, and the current admin job."""
        return _admin(request).inspect()

    @app.get("/admin/sources", tags=["admin"])
    def admin_sources(request: Request) -> dict[str, Any]:
        """Discoverable inputs_v2 sources (slow first call)."""
        admin = _admin(request)
        return {
            "sources": admin.available_sources(),
            "build_available": admin.build_capabilities()["build_available"],
        }

    @app.post("/admin/releases", tags=["admin"], status_code=201)
    def publish_release(request: Request, payload: dict[str, Any]) -> dict[str, Any]:
        admin = _admin(request)
        try:
            with request.app.state.queries_lock:
                return ReleaseStore(admin.data_root).publish(payload)
        except FileExistsError as exc:
            raise HTTPException(status_code=409, detail="OmniPath release already exists") from exc
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.get("/admin/jobs", tags=["admin"], response_model=AdminJobsResponse)
    def admin_jobs(request: Request) -> dict[str, Any]:
        snapshot = _admin(request).inspect()
        return {"job": snapshot["job"], "jobs": snapshot["jobs"]}

    @app.post("/admin/jobs", tags=["admin"], status_code=202, response_model=AdminJobResponse)
    def start_admin_job(request: Request, payload: AdminJobRequest) -> dict[str, Any]:
        admin = _admin(request)
        body = payload.model_dump()
        action = body.pop("action")
        params = {key: value for key, value in body.items() if value is not None}
        if action == "build_resource" and not params.get("source") and not params.get("sources"):
            raise HTTPException(status_code=400, detail="source or sources is required")
        if action == "build_resources" and not params.get("sources") and not params.get("source"):
            raise HTTPException(status_code=400, detail="sources is required")
        if action in {"build_resource", "build_resources"}:
            for source in params.get("sources") or [params["source"]]:
                version = (params.get("versions") or {}).get(source, params.get("version"))
                if not isinstance(version, str) or not re.fullmatch(
                    r"[0-9]+(?:\.[0-9]+)*", version
                ):
                    raise HTTPException(
                        status_code=400, detail=f"Explicit numeric version required for {source}"
                    )
                if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", source):
                    raise HTTPException(status_code=400, detail="Invalid resource name")
                if (admin.data_root / "resources" / source / version).exists():
                    raise HTTPException(
                        status_code=409, detail=f"Resource {source}/{version} already exists"
                    )
        try:
            return admin.start(action, params)
        except BuildUnavailable as exc:
            raise HTTPException(status_code=503, detail=str(exc)) from exc
        except RuntimeError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.get("/admin/logs/{kind}/{name}", tags=["admin"])
    def admin_read_log(
        request: Request, kind: str, name: str, version: str | None = None
    ) -> dict[str, Any]:
        try:
            return _admin(request).read_log(kind, name, version)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.delete("/admin/resources/{resource}/versions/{version}", tags=["admin"])
    @app.delete("/admin/resources/{resource}/{version}", tags=["admin"])
    def delete_admin_resource_version(
        request: Request, resource: str, version: str
    ) -> dict[str, Any]:
        try:
            admin = _admin(request)
            queries = getattr(request.app.state, "queries", None)
            scope = queries.query_scope() if hasattr(queries, "query_scope") else nullcontext()
            with request.app.state.queries_lock, scope:
                res = admin.delete_resource_version(resource, version)
                if queries is not None and hasattr(queries, "reload_resources"):
                    queries.reload_resources()
                elif queries is not None and hasattr(queries, "discover_resources"):
                    queries.discover_resources()
            return res
        except QueryCapacityError as exc:
            raise HTTPException(
                status_code=503, detail=str(exc), headers={"Retry-After": "1"}
            ) from exc
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

    @app.get("/admin/jobs/{job_id}", tags=["admin"], response_model=AdminJobResponse)
    def get_admin_job(request: Request, job_id: str) -> dict[str, Any]:
        job = _admin(request).get_job(job_id)
        if job is None:
            raise HTTPException(status_code=404, detail="Unknown job")
        return job

    @app.post("/admin/jobs/{job_id}/cancel", tags=["admin"], response_model=AdminJobResponse)
    def cancel_admin_job(request: Request, job_id: str) -> dict[str, Any]:
        try:
            return _admin(request).cancel(job_id)
        except KeyError:
            raise HTTPException(status_code=404, detail="Unknown job") from None

    @app.get("/admin/jobs/{job_id}/events", tags=["admin"], response_model=None)
    def admin_job_events(request: Request, job_id: str) -> StreamingResponse:
        admin = _admin(request)
        if admin.get_job(job_id) is None:
            raise HTTPException(status_code=404, detail="Unknown job")

        def generate() -> Iterator[str]:
            for snapshot in admin.events(job_id):
                yield f"data: {json.dumps(snapshot, default=str)}\n\n"

        return StreamingResponse(
            generate(),
            media_type="text/event-stream",
            headers={
                "Cache-Control": "no-cache",
                "Connection": "keep-alive",
                "X-Accel-Buffering": "no",
            },
        )
