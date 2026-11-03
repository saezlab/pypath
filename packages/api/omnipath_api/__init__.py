"""OmniPath REST API over the resource Parquet files of a serving root.

Serves analytical queries, autocomplete, faceted filtering and binary exports over
HTTP; the queries are ``omnipath_api.queries``.
"""

from .server import create_app, run_server, serve

__all__ = [
    "create_app",
    "run_server",
    "serve",
]
