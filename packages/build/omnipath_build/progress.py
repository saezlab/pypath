"""Progress and cancellation hooks for long-running build jobs."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

ProgressCallback = Callable[[dict[str, Any]], None]
CancelCheck = Callable[[], bool]


class BuildCancelled(Exception):
    """Raised when an in-flight build or hub export is cancelled."""


def emit(on_progress: ProgressCallback | None, **event: Any) -> None:
    if on_progress is None:
        return
    on_progress(event)


def check_cancel(should_cancel: CancelCheck | None) -> None:
    if should_cancel is not None and should_cancel():
        raise BuildCancelled()
