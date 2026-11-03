from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import Response

from omnipath_api.models import ExportRequest
from omnipath_api.routers.common import _call, _dump

router = APIRouter()


def _attachment(data: bytes, content_type: str, name: str) -> Response:
    return Response(
        content=data,
        media_type=content_type,
        headers={
            "Content-Disposition": f'attachment; filename="{name}.parquet"',
            "Access-Control-Allow-Origin": "*",
        },
    )


@router.post("/export", tags=["export"], response_model=None)
def export_slice(request: Request, payload: ExportRequest) -> Response:
    """A filtered relation slice as Parquet: at most 100,000 relations, with their
    annotations and evidence."""
    body = _dump(payload)
    data, content_type = _call(
        request,
        "export_slice",
        filters=body.get("filters") or {},
        resources=body.get("resources"),
        limit=body["limit"],
    )
    return _attachment(data, content_type, "omnipath_relations")
