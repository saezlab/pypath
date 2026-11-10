"""The shared step-file queries by name: the same SQL the clients run locally."""

from __future__ import annotations

import io
from typing import Any, Literal

from fastapi import APIRouter, Body, HTTPException, Request
from fastapi.responses import Response

from omnipath_client import steps
from omnipath_api.routers.common import _call

router = APIRouter()
PARQUET = "application/vnd.apache.parquet"


@router.get("/query", tags=["query"])
def list_queries() -> list[dict[str, Any]]:
    """The queries, their parameters and their outputs."""
    out = []
    for name in steps.names():
        query = steps.load(name)
        out.append(
            {
                "name": name,
                "parameters": [
                    {"name": p.name, "type": p.type, "default": p.default}
                    for p in query.params.values()
                ],
                "outputs": sorted(query.outputs),
            }
        )
    return out


@router.post("/query/{name}", tags=["query"], response_model=None)
def run_named_query(
    request: Request,
    name: str,
    params: dict[str, Any] = Body(default_factory=dict),
    format: Literal["json", "parquet"] = "json",
) -> Response:
    """Run a query by name with its parameters (the JSON body): its JSON answer, or its
    table as Parquet. ``resources`` selects resources of the release (the query policy
    still applies)."""
    if name not in steps.names():
        raise HTTPException(status_code=404, detail=f"Unknown query: {name}")
    params = dict(params or {})
    resources = params.pop("resources", None)
    output = "json" if format == "json" else "table"
    if output not in steps.load(name).outputs:
        raise HTTPException(status_code=422, detail=f"{name} has no {format} answer")
    try:
        result = _call(request, "run_query", name, params, resources, output=output)
    except TypeError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    if output == "json":
        return Response(content=result if result is not None else "null", media_type="application/json")
    import pyarrow.parquet as pq

    sink = io.BytesIO()
    pq.write_table(result, sink, compression="zstd")
    return Response(content=sink.getvalue(), media_type=PARQUET)
