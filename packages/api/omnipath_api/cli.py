"""Command-line interface for the OmniPath full-parquet presentation server.

Provides commands to start the embedded DuckDB query server and REST API.
"""

from __future__ import annotations

import argparse
import sys

from omnipath_api.settings import settings
from omnipath_api.server import run_server


def main(args: list[str] | None = None) -> int:
    """CLI entrypoint."""
    parser = argparse.ArgumentParser(
        prog="omnipath-api",
        description="Serve OmniPath knowledge graphs directly from the resource Parquet tables via DuckDB.",
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
        "--data-root", default=str(settings.data_root), help="Directory containing resources/ data"
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
