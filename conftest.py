"""Explicit opt-in integration tests own a disposable local PostgreSQL cluster."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import pytest


def _worker_port() -> int:
    """Per-xdist-worker port: gw0 -> 55479, gw1 -> 55480, ... (controller/serial -> 55479)."""
    worker = os.environ.get("PYTEST_XDIST_WORKER", "gw0")
    return 55479 + int(worker.removeprefix("gw") or 0)


@pytest.fixture(scope="session")
def postgres_dsn():
    # Session scope is per xdist worker process: every worker initialises its own cluster in a
    # private mkdtemp directory (the server is reachable only through a unix socket inside that
    # directory), so parallel integration runs never share data or sockets.
    port = _worker_port()
    if os.environ.get("OMNIPATH_TEST_POSTGRES") != "1":
        pytest.skip("Set OMNIPATH_TEST_POSTGRES=1 to start a disposable PostgreSQL cluster")
    candidates = [
        os.environ.get("OMNIPATH_POSTGRES_BIN", ""),
        "/opt/homebrew/opt/postgresql@16/bin",
        str(Path(shutil.which("initdb") or "/missing/initdb").parent),
    ]
    binary = next(
        (
            Path(candidate)
            for candidate in candidates
            if candidate
            and all((Path(candidate) / name).is_file() for name in ("initdb", "pg_ctl"))
        ),
        None,
    )
    if binary is None:
        pytest.fail("Set OMNIPATH_POSTGRES_BIN to a directory containing initdb and pg_ctl")
    root = Path(tempfile.mkdtemp(prefix="op-pg-", dir="/tmp"))
    root.chmod(0o700)
    cluster = root / "cluster"
    logfile = root / "postgres.log"
    environment = {**os.environ, "LC_ALL": "C"}

    def run(*args):
        result = subprocess.run(args, text=True, capture_output=True, env=environment, timeout=30)
        if result.returncode:
            log = logfile.read_text() if logfile.exists() else ""
            raise RuntimeError(f"Temporary PostgreSQL failed: {result.stderr}\n{log}")
        return result

    try:
        run(
            str(binary / "initdb"),
            "-D",
            str(cluster),
            "-U",
            "omnipath_migration",
            "--auth=trust",
            "--no-locale",
            "-E",
            "UTF8",
        )
        run(
            str(binary / "pg_ctl"),
            "-D",
            str(cluster),
            "-l",
            str(logfile),
            "-w",
            "start",
            "-o",
            f"-F -c listen_addresses='' -c unix_socket_directories='{root}' "
            f"-p {port} -c max_connections=20 -c shared_buffers=32MB",
        )
        yield f"host={root} port={port} user=omnipath_migration dbname=postgres"
    finally:
        if (cluster / "postmaster.pid").exists():
            run(str(binary / "pg_ctl"), "-D", str(cluster), "-m", "fast", "-w", "stop")
        shutil.rmtree(root)


def pytest_collection_modifyitems(config, items):
    # Database tests are integration tests whether or not they carry the marker,
    # so `-m integration` selects them and `-m "not integration"` deselects them.
    marker = pytest.mark.integration
    for item in items:
        if "postgres_dsn" in getattr(item, "fixturenames", ()):
            item.add_marker(marker)
