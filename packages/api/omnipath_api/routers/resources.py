from __future__ import annotations
from pathlib import Path
from typing import Any
from fastapi import HTTPException, Query, Request
from fastapi.responses import FileResponse

from omnipath_api.models import (
    ResourceDownloadRequest,
    ResourceFilesResponse,
    ResourcesResponse,
    ResourceTopicsResponse,
)
from fastapi import APIRouter
from omnipath_api.routers.common import _call, _zip_files

router = APIRouter()


@router.get("/releases", tags=["resources"])
def list_releases(request: Request) -> dict[str, Any]:
    store = request.app.state.queries.releases
    return {"default": store.default(), "releases": store.list()}


@router.get("/resources", response_model=ResourcesResponse, tags=["resources"])
@router.get(
    "/resources/catalog",
    response_model=ResourcesResponse,
    tags=["resources"],
    include_in_schema=False,
)
def list_resources(
    request: Request,
    shape: str = Query("web", description="Pass `svelte` for the explorer catalog shape."),
) -> dict[str, Any]:
    """List indexed resources. `shape=svelte` returns the explorer catalog."""
    svelte_shape = request.url.path.rstrip("/").endswith("/catalog") or shape == "svelte"
    if svelte_shape:
        resources = _call(request, "list_resource_catalog")
    else:
        resources = _call(request, "list_resources")
    return {"resources": resources}


@router.get("/resources/topics", response_model=ResourceTopicsResponse, tags=["resources"])
def list_resource_topics(request: Request) -> dict[str, Any]:
    """The biological topics of the release, each with the resources that have it."""
    return {"topics": _call(request, "list_resource_topics")}


@router.post("/resources/download", tags=["resources"], response_model=None)
def download_resources(request: Request, payload: ResourceDownloadRequest) -> FileResponse:
    """Zip parquet files for one or more resources. Evidence payloads are omitted unless requested."""
    resource_ids = [str(item).strip() for item in payload.resource_ids if str(item).strip()]
    if not resource_ids:
        raise HTTPException(status_code=400, detail="At least one resource_id is required")
    entries: list[tuple[str, Path]] = []
    for resource_id in resource_ids:
        files = _call(
            request,
            "iter_resource_file_paths",
            resource_id,
            include_evidence=payload.include_evidence,
        )
        if not files:
            raise HTTPException(status_code=404, detail=f"Resource '{resource_id}' not found")
        prefix = resource_id if len(resource_ids) > 1 else ""
        for name, path in files:
            entries.append((f"{prefix}/{name}" if prefix else name, path))
    zip_name = resource_ids[0] if len(resource_ids) == 1 else "omnipath-resources"
    return _zip_files(entries, f"{zip_name}.zip")


@router.get(
    "/resources/{resource_id}/files",
    response_model=ResourceFilesResponse,
    tags=["resources"],
)
def list_resource_files(
    request: Request,
    resource_id: str,
    include_evidence: bool = Query(
        False, description="Include evidence_payloads.parquet in the listing."
    ),
) -> dict[str, Any]:
    """List downloadable parquet artifacts (and column schema) for a resource."""
    result = _call(request, "list_resource_files", resource_id, include_evidence=include_evidence)
    if not result:
        raise HTTPException(status_code=404, detail=f"Resource '{resource_id}' not found")
    return result


@router.get("/resources/{resource_id}/files/{filename}", tags=["resources"], response_model=None)
def download_resource_file(request: Request, resource_id: str, filename: str) -> FileResponse:
    """Download one parquet file from a resource release."""
    path = _call(request, "get_resource_file_path", resource_id, filename)
    if not path:
        raise HTTPException(
            status_code=404, detail=f"File '{filename}' not found for resource '{resource_id}'"
        )
    stem = filename.removesuffix(".parquet")
    return FileResponse(
        path,
        media_type="application/vnd.apache.parquet",
        filename=f"{resource_id}-{stem}.parquet",
    )


@router.get("/resources/{resource_id}/download", tags=["resources"], response_model=None)
def download_resource(
    request: Request,
    resource_id: str,
    include_evidence: bool = Query(
        False, description="Include evidence_payloads.parquet in the zip."
    ),
) -> FileResponse:
    """Download resource parquet files as a zip. Evidence payloads are omitted unless requested."""
    files = _call(
        request, "iter_resource_file_paths", resource_id, include_evidence=include_evidence
    )
    if not files:
        raise HTTPException(status_code=404, detail=f"Resource '{resource_id}' not found")
    return _zip_files(files, f"{resource_id}.zip")
