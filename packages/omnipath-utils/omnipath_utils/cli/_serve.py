"""Serve command for omnipath-utils web service."""

from __future__ import annotations

import logging
import argparse


def serve_cmd(args: list[str]) -> None:
    """Start the web service."""
    parser = argparse.ArgumentParser(
        description='Start omnipath-utils web service',
    )
    parser.add_argument(
        '--db-url',
        default=None,
        help='PostgreSQL connection URL',
    )
    parser.add_argument(
        '--host',
        default='0.0.0.0',
        help='Host to bind to',
    )
    parser.add_argument(
        '--port',
        type=int,
        default=8082,
        help='Port to listen on',
    )
    parser.add_argument('-v', '--verbose', action='store_true')

    opts = parser.parse_args(args)

    if opts.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    import uvicorn

    from omnipath_utils.server._app import create_app

    app = create_app(db_url=opts.db_url)
    uvicorn.run(app, host=opts.host, port=opts.port)
