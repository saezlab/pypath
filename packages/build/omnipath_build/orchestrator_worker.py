"""Isolated resource worker. Launched by orchestrator, not a public CLI."""

from __future__ import annotations

import json
import os
from pathlib import Path
import resource
import signal
import sys
import time
import traceback


def main():
    spec = json.loads(Path(sys.argv[1]).read_text())
    os.environ["OMNIPATH_BUILD_DUCKDB_THREADS"] = str(spec["cpus"])
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[name] = "1"
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, spec["allowed_cpus"])
    os.environ["TMPDIR"] = spec["tmpdir"]
    Path(spec["tmpdir"]).mkdir(parents=True, exist_ok=True)
    cancelled = False

    def cancel(*_):
        nonlocal cancelled
        cancelled = True

    signal.signal(signal.SIGTERM, cancel)
    signal.signal(signal.SIGINT, cancel)
    from .pipeline import build_resource
    from .progress import BuildCancelled

    started = time.monotonic()
    result = {}
    code = 1
    with Path(spec["events"]).open("a", buffering=1) as events:

        def progress(event):
            events.write(json.dumps(event, default=str) + "\n")

        try:
            result = build_resource(
                **spec["build"],
                progress=False,
                on_progress=progress,
                should_cancel=lambda: cancelled,
            )
            result["status"] = "success"
            code = 0
        except BuildCancelled:
            result = {"status": "cancelled", "error": "Build cancelled"}
            code = 130
        except Exception as exc:
            traceback.print_exc()
            result = {
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
                "memory_failure": type(exc).__name__ in {"MemoryError", "OutOfMemoryException"},
            }
        finally:
            usage = resource.getrusage(resource.RUSAGE_SELF)
            result.update(
                elapsed_seconds=time.monotonic() - started,
                peak_rss_bytes=int(usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)),
                cpu_seconds=usage.ru_utime + usage.ru_stime,
            )
            target = Path(spec["result"])
            temp = target.with_suffix(".tmp")
            temp.write_text(json.dumps(result, default=str, indent=2) + "\n")
            temp.replace(target)
    return code


if __name__ == "__main__":
    sys.exit(main())
