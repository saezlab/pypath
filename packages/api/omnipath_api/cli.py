"""Command-line interface of the OmniPath API.

``omnipath-api [--host --port --data-root]`` serves the API. ``omnipath-api
build-release-db <version>`` writes the database a published release is served from
(``serving/<version>.duckdb``) and keeps the newest few; ``omnipath-api
prune-release-dbs`` only prunes.
"""

from __future__ import annotations

import argparse
import json
import sys

from omnipath_api.settings import settings


def _release_db(args: list[str]) -> int:
    from omnipath_api.release_database import DEFAULT_KEEP, build_served_database

    parser = argparse.ArgumentParser(
        prog="omnipath-api build-release-db",
        description="Write the DuckDB database a published release is served from, then "
        "remove all but the newest --keep databases.",
    )
    parser.add_argument("version", help="a published release (releases/<version>.json)")
    parser.add_argument("--data-root", default=str(settings.data_root))
    parser.add_argument(
        "--keep",
        type=int,
        default=DEFAULT_KEEP,
        help=f"databases to keep (default {DEFAULT_KEEP}); 0 keeps all",
    )
    parser.add_argument("--memory-limit", default="4GB")
    parser.add_argument("--threads", type=int, default=None)
    parsed = parser.parse_args(args)
    record = build_served_database(
        parsed.data_root,
        parsed.version,
        keep=parsed.keep or None,
        memory_limit=parsed.memory_limit,
        threads=parsed.threads,
        progress=lambda message: print(message, file=sys.stderr, flush=True),
    )
    print(json.dumps(record, indent=2))
    return 0


def _prune(args: list[str]) -> int:
    from omnipath_api.release_database import DEFAULT_KEEP, prune_served_databases

    parser = argparse.ArgumentParser(prog="omnipath-api prune-release-dbs")
    parser.add_argument("--data-root", default=str(settings.data_root))
    parser.add_argument("--keep", type=int, default=DEFAULT_KEEP)
    parsed = parser.parse_args(args)
    print(json.dumps({"removed": prune_served_databases(parsed.data_root, parsed.keep)}))
    return 0


def main(args: list[str] | None = None) -> int:
    """CLI entrypoint."""
    args = sys.argv[1:] if args is None else list(args)
    if args and args[0] == "build-release-db":
        return _release_db(args[1:])
    if args and args[0] == "prune-release-dbs":
        return _prune(args[1:])
    from omnipath_api.server import run_server

    parser = argparse.ArgumentParser(
        prog="omnipath-api",
        description="Serve OmniPath's published releases from their DuckDB databases. "
        "Subcommands: build-release-db, prune-release-dbs.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=settings.port,
        help=f"Port to bind server (default: {settings.port})",
    )
    parser.add_argument(
        "--host", default=settings.host, help=f"Host address to bind (default: {settings.host})"
    )
    parser.add_argument(
        "--data-root",
        default=str(settings.data_root),
        help="Directory containing resources/, releases/ and serving/",
    )

    parsed = parser.parse_args(args)
    run_server(
        host=parsed.host,
        port=parsed.port,
        data_root=parsed.data_root,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
