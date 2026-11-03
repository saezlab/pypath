from __future__ import annotations
from contextlib import nullcontext
import json
import tempfile
import zipfile
from pathlib import Path
from typing import Any
from fastapi import HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse, Response
from starlette.background import BackgroundTask
from pydantic import ValidationError
from omnipath_api.queries.store.query_pool import QueryCapacityError


class EngineJSONResponse(JSONResponse):
    def render(self, content: Any) -> bytes:
        return json.dumps(content, default=str, ensure_ascii=False).encode("utf-8")


def _parse_json_param(raw: str | None, default: Any) -> Any:
    if raw is None or raw == "":
        return default
    try:
        value = json.loads(raw)
    except (json.JSONDecodeError, TypeError) as exc:
        raise HTTPException(status_code=422, detail="Invalid JSON query parameter") from exc
    if not (value is None and default is None) and not isinstance(
        value, type(default) if default is not None else dict
    ):
        raise HTTPException(status_code=422, detail=f"Expected JSON {type(default).__name__}")
    return value


def _cursor_dict(cursor: Any) -> dict[str, Any] | None:
    if cursor is None:
        return None
    if hasattr(cursor, "model_dump"):
        return cursor.model_dump()
    if isinstance(cursor, dict):
        return cursor
    return None


def _dump(model: Any) -> dict[str, Any]:
    if hasattr(model, "model_dump"):
        return model.model_dump()
    return dict(model) if model else {}


def _call(request: Request, method: str, *args: Any, **kwargs: Any) -> Any:
    queries = request.app.state.queries
    try:
        with queries.query_scope() if hasattr(queries, "query_scope") else nullcontext():
            release = getattr(request.state, "omnipath_release", None)
            if release is None:
                release = (
                    request.query_params.get("release")
                    or request.headers.get("x-omnipath-release")
                    or request.cookies.get("omnipath_release")
                    or queries.releases.default()
                )
                request.state.omnipath_release = release
            snapshot = getattr(request.state, "inventory_snapshot", None)
            if snapshot is None and hasattr(queries, "inventory_snapshot"):
                snapshot = queries.inventory_snapshot()
                request.state.inventory_snapshot = snapshot
            scope = (
                queries.release_scope(release, snapshot=snapshot)
                if snapshot is not None
                else queries.release_scope(release)
            )
            with scope:
                return getattr(queries, method)(*args, **kwargs)
    except QueryCapacityError as exc:
        raise HTTPException(status_code=503, detail=str(exc), headers={"Retry-After": "1"}) from exc
    except ValidationError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


def _timed(response: Response, result: dict[str, Any]) -> dict[str, Any]:
    elapsed = result.get("elapsed_ms")
    if elapsed is not None:
        response.headers["X-Query-Duration-Ms"] = str(elapsed)
    return result


def _unlink(path: Path) -> None:
    path.unlink(missing_ok=True)


def _zip_files(entries: list[tuple[str, Path]], filename: str) -> FileResponse:
    tmp = tempfile.NamedTemporaryFile(prefix="omnipath_", suffix=".zip", delete=False)
    tmp_path = Path(tmp.name)
    tmp.close()
    with zipfile.ZipFile(tmp_path, "w", compression=zipfile.ZIP_STORED) as archive:
        for archive_name, source in entries:
            archive.write(source, arcname=archive_name)
    return FileResponse(
        tmp_path,
        media_type="application/zip",
        filename=filename,
        background=BackgroundTask(_unlink, tmp_path),
    )
