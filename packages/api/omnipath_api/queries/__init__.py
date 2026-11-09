"""The API's queries over the release databases of a serving root.

Analytical queries, autocomplete, faceted filtering and binary exports over one DuckDB
database per published release (``serving/<version>.duckdb``), for the API's routes.
The database is the format the clients (``packages/query``) write locally too.
"""

from .release import ReleaseQueries

__all__ = ["ReleaseQueries"]
