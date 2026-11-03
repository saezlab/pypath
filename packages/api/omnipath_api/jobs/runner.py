"""Job runner and lifecycle management."""

from __future__ import annotations

import time
import uuid
from pathlib import Path
from typing import Any


def make_stage(stage_id: str) -> dict[str, Any]:
    return {
        "id": stage_id,
        "status": "pending",
        "message": "",
        "current": None,
        "total": None,
        "counts": {},
        "started_at": None,
        "finished_at": None,
        "elapsed_seconds": None,
    }


def stage_elapsed(stage: dict[str, Any], now: float) -> float | None:
    started = stage.get("started_at")
    if started is None:
        return None
    finished = stage.get("finished_at")
    end = finished if finished is not None else now
    return round(end - started, 2)


def set_stage_status(stage: dict[str, Any], status: str, now: float | None = None) -> None:
    stamp = time.time() if now is None else now
    previous = stage.get("status")
    stage["status"] = status
    if status == "running":
        if previous != "running" or stage.get("started_at") is None:
            stage["started_at"] = stamp
        stage["finished_at"] = None
    elif status in {"done", "failed", "cancelled"}:
        if stage.get("started_at") is None:
            stage["started_at"] = stamp
        if previous == "running" or stage.get("finished_at") is None:
            stage["finished_at"] = stamp
    stage["elapsed_seconds"] = stage_elapsed(stage, stamp)


class Job:
    """Represents an asynchronous background job with progression stages and logs."""

    def __init__(self, action: str, params: dict[str, Any], stages: list[str]) -> None:
        self.id = uuid.uuid4().hex[:12]
        self.action = action
        self.params = params
        self.status = "queued"
        self.stages: list[dict[str, Any]] = [make_stage(stage_id) for stage_id in stages]
        self.logs: list[str] = []
        self.result: Any = None
        self.error: str | None = None
        self.created_at = time.time()
        self.started_at: float | None = None
        self.finished_at: float | None = None
        self.seq = 0
        self.cancel_requested = False
        self.log_path: Path | None = None

    def snapshot(self) -> dict[str, Any]:
        now = time.time()
        elapsed = 0.0
        if self.started_at is not None:
            end = self.finished_at if self.finished_at is not None else now
            elapsed = round(end - self.started_at, 2)
        stages = []
        for stage in self.stages:
            item = dict(stage)
            item["elapsed_seconds"] = stage_elapsed(item, now)
            stages.append(item)
        return {
            "id": self.id,
            "action": self.action,
            "params": self.params,
            "status": self.status,
            "stages": stages,
            "logs": list(self.logs[-80:]),
            "result": self.result,
            "error": self.error,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "elapsed_seconds": elapsed,
            "seq": self.seq,
            "cancel_requested": self.cancel_requested,
            "log_path": str(self.log_path) if self.log_path else None,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Job:
        job = cls(str(data.get("action") or ""), dict(data.get("params") or {}), [])
        job.id = str(data.get("id") or job.id)
        job.status = str(data.get("status") or "done")
        job.stages = [dict(s) for s in data.get("stages", [])]
        job.logs = [str(line) for line in data.get("logs", [])]
        job.result = data.get("result")
        job.error = data.get("error")
        job.created_at = float(data.get("created_at") or job.created_at)
        job.started_at = float(data["started_at"]) if data.get("started_at") is not None else None
        job.finished_at = (
            float(data["finished_at"]) if data.get("finished_at") is not None else None
        )
        job.seq = int(data.get("seq") or 0)
        job.cancel_requested = bool(data.get("cancel_requested"))
        log_path = data.get("log_path")
        job.log_path = Path(log_path) if log_path else None
        return job
