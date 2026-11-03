"""Run durable build jobs independently from serving: python -m omnipath_api.jobs.worker."""

from __future__ import annotations

import argparse
import fcntl
import logging
import signal
import threading
from contextlib import contextmanager
from pathlib import Path

from omnipath_api.admin import AdminService
from omnipath_api.settings import Settings


@contextmanager
def worker_lock(data_root):
    path = Path(data_root) / "jobs" / "worker.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Another build worker owns this data root") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def run_worker(data_root, *, ops=None, stop=None, once=False, poll_interval=1.0):
    """One supervised worker per data root; multiple serving processes may enqueue."""
    stop = stop or threading.Event()
    with worker_lock(data_root):
        service = AdminService(data_root, ops=ops)
        service.job_store.recover_interrupted()
        import time

        capabilities = service.build_capabilities(use_worker=False)
        capabilities["recorded_at"] = time.time()
        capabilities["sources"] = (
            service.ops.list_sources()
            if capabilities["pipeline_available"] and hasattr(service.ops, "list_sources")
            else []
        )
        service.job_store.publish_capabilities(capabilities)
        while not stop.is_set():
            job = service.job_store.claim()
            if job is not None:
                service.run_job(job, should_stop=stop.is_set)
            if once:
                return
            if job is None:
                stop.wait(poll_interval)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Settings().data_root)
    parser.add_argument(
        "--once", action="store_true", help="Process at most one queued job and exit"
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    stop = threading.Event()
    for signum in (signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, lambda *_: stop.set())
    run_worker(args.data_root, stop=stop, once=args.once)


if __name__ == "__main__":
    main()
