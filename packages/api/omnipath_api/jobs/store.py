"""Durable single-flight queue shared by API processes and the build worker."""

from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from pathlib import Path

from .runner import Job


class JobStore:
    def __init__(self, data_root):
        self.path = Path(data_root) / "jobs" / "queue.sqlite3"

    @contextmanager
    def connection(self, *, write=False):
        if write:
            self.path.parent.mkdir(parents=True, exist_ok=True)
        con = sqlite3.connect(
            self.path if write else self.path.as_uri() + "?mode=ro", uri=not write, timeout=30
        )
        try:
            if write:
                con.execute("PRAGMA synchronous=FULL")
                con.execute(
                    "CREATE TABLE IF NOT EXISTS jobs (id TEXT PRIMARY KEY, status TEXT NOT NULL, cancel INTEGER NOT NULL DEFAULT 0, created REAL NOT NULL, snapshot TEXT NOT NULL)"
                )
                con.execute("BEGIN IMMEDIATE")
            yield con
            if write:
                con.commit()
        except BaseException:
            con.rollback()
            raise
        finally:
            con.close()

    @staticmethod
    def _job(row):
        job = Job.from_dict(json.loads(row[0]))
        job.status = row[1]
        job.cancel_requested = bool(row[2])
        return job

    @staticmethod
    def _has_jobs(con):
        # A writer creates the file before its table; readers may observe the gap.
        return bool(
            con.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='jobs'").fetchone()
        )

    def list(self):
        if not self.path.is_file():
            return []
        with self.connection() as con:
            if not self._has_jobs(con):
                return []
            return [
                self._job(row)
                for row in con.execute(
                    "SELECT snapshot, status, cancel FROM jobs ORDER BY created DESC LIMIT 50"
                )
            ]

    def get(self, job_id):
        if not self.path.is_file():
            return None
        with self.connection() as con:
            if not self._has_jobs(con):
                return None
            row = con.execute(
                "SELECT snapshot, status, cancel FROM jobs WHERE id=?", [job_id]
            ).fetchone()
            return self._job(row) if row else None

    def enqueue(self, job):
        with self.connection(write=True) as con:
            if con.execute(
                "SELECT 1 FROM jobs WHERE status IN ('queued','running') LIMIT 1"
            ).fetchone():
                raise RuntimeError("A job is already queued or running")
            con.execute(
                "INSERT INTO jobs (id,status,cancel,created,snapshot) VALUES (?,?,?,?,?)",
                [job.id, job.status, 0, job.created_at, json.dumps(job.snapshot(), default=str)],
            )

    def save(self, job):
        with self.connection(write=True) as con:
            con.execute(
                "UPDATE jobs SET status=?,snapshot=? WHERE id=?",
                [job.status, json.dumps(job.snapshot(), default=str), job.id],
            )

    def cancel(self, job_id):
        with self.connection(write=True) as con:
            row = con.execute(
                "SELECT snapshot,status,cancel FROM jobs WHERE id=?", [job_id]
            ).fetchone()
            if row is None:
                raise KeyError(job_id)
            job = self._job(row)
            if job.status in {"queued", "running"}:
                job.cancel_requested = True
                job.seq += 1
                if job.status == "queued":
                    import time

                    job.status = "cancelled"
                    job.finished_at = time.time()
                con.execute(
                    "UPDATE jobs SET cancel=1,status=?,snapshot=? WHERE id=?",
                    [job.status, json.dumps(job.snapshot(), default=str), job_id],
                )
            return job

    def claim(self):
        with self.connection(write=True) as con:
            row = con.execute(
                "SELECT snapshot,status,cancel FROM jobs WHERE status='queued' ORDER BY created LIMIT 1"
            ).fetchone()
            if row is None:
                return None
            job = self._job(row)
            job.status = "running"
            con.execute("UPDATE jobs SET status='running' WHERE id=?", [job.id])
            return job

    def capabilities(self):
        if not self.path.is_file():
            return None
        with self.connection() as con:
            if not con.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='worker_metadata'"
            ).fetchone():
                return None
            row = con.execute(
                "SELECT value FROM worker_metadata WHERE key='capabilities'"
            ).fetchone()
            return json.loads(row[0]) if row else None

    def publish_capabilities(self, capabilities):
        with self.connection(write=True) as con:
            con.execute(
                "CREATE TABLE IF NOT EXISTS worker_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
            )
            con.execute(
                "INSERT OR REPLACE INTO worker_metadata VALUES ('capabilities',?)",
                [json.dumps(capabilities)],
            )

    def recover_interrupted(self):
        """Called only while holding the OS-owned worker lock, never by the API."""
        import time

        with self.connection(write=True) as con:
            rows = con.execute(
                "SELECT snapshot,status,cancel FROM jobs WHERE status='running'"
            ).fetchall()
            for row in rows:
                job = self._job(row)
                job.status = "failed"
                job.error = "Build worker stopped before completion; inspect outputs before submitting a retry"
                job.finished_at = time.time()
                job.seq += 1
                con.execute(
                    "UPDATE jobs SET status=?,snapshot=? WHERE id=?",
                    [job.status, json.dumps(job.snapshot(), default=str), job.id],
                )
