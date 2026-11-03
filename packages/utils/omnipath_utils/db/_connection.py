"""PostgreSQL connection management."""

from __future__ import annotations

import os
import logging

_log = logging.getLogger(__name__)

DEFAULT_DB_URL = (
    'postgresql+psycopg://postgres:dev@localhost:5433/omnipath_utils'
)
SCHEMA = 'omnipath_utils'


def get_db_url() -> str:
    """Get DB URL from env or default."""
    return os.environ.get('OMNIPATH_UTILS_DB_URL', DEFAULT_DB_URL)


def get_engine(db_url: str | None = None, echo: bool = False):
    """Create SQLAlchemy engine."""
    from sqlalchemy import create_engine

    url = db_url or get_db_url()
    _log.info('Connecting to %s', url.split('@')[-1])  # hide password
    return create_engine(url, echo=echo)


def get_connection(db_url: str | None = None):
    """Get a psycopg3 connection directly (for COPY operations).

    Accepts either a psycopg or SQLAlchemy-style URL.
    The ``postgresql+psycopg://`` prefix is normalised to
    ``postgresql://`` so that plain psycopg can parse it.
    """
    import psycopg

    url = db_url or get_db_url()
    url = url.replace('postgresql+psycopg://', 'postgresql://')
    return psycopg.connect(url)


def ensure_schema(engine):
    """Create the omnipath_utils schema if it doesn't exist."""
    from sqlalchemy import text

    with engine.connect() as conn:
        conn.execute(text(f'CREATE SCHEMA IF NOT EXISTS {SCHEMA}'))
        conn.commit()
    _log.info('Schema %s ensured', SCHEMA)
