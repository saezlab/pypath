"""The API's queries over a serving root of resource Parquet files.

Analytical queries, autocomplete, faceted filtering and binary exports directly from
the resource Parquet tables via DuckDB, for the API's routes. The clients do not use
this package: their search, relation export and entity details are the SQL macros of
``packages/query/sql``, which give the same answers (the client tests compare them).
"""

from .parquet import ParquetQueries

__all__ = ["ParquetQueries"]
