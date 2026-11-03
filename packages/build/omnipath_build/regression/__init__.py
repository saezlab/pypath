"""Resolution regression harness.

Extract pre-resolution observations (queries and votes) straight from source
extraction, resolve them with any runtime that offers ``resolve(queries, votes)``
and diff two result sets. Built resource Parquet outputs are never read: they are
themselves products of resolution.
"""
