"""Products built in PostgreSQL over the committed normalized relational tables.

The exclusion metabolites and MetSigDB are not among them: they are built from the published
release as Parquet and loaded.
"""

from dataclasses import asdict, is_dataclass
import time
from psycopg2 import sql
from . import cosmos, network_views


def _search_path(conn, schema):
    with conn.cursor() as cur:
        cur.execute(sql.SQL("SET search_path = {}, public").format(sql.Identifier(schema)))


def _reported(value):
    return asdict(value) if is_dataclass(value) else value


class CommitDeferredConnection:
    """Keep imported main product writes atomic until the loader checkpoint.

    Real cursors and every operation except commit are delegated unchanged.
    Rollback must stay real so failures/cancellation undo the unfinished product.
    """

    def __init__(self, connection):
        self.connection = connection

    def commit(self):
        return None

    def __getattr__(self, name):
        return getattr(self.connection, name)


def run_product(conn, product, *, schema="public", progress=False):
    """Publish main's product schema from the already committed shared tables."""
    conn = CommitDeferredConnection(conn)
    _search_path(conn, schema)
    started = time.monotonic()
    if product == "network_views":
        stats = network_views.apply_all(
            conn,
            network_views.NETWORKS,
            registry_schema=schema,
            log=(lambda line: print(line, flush=True)) if progress else (lambda line: None),
        )
    elif product == "cosmos":
        stats = cosmos.build_cosmos_projection(
            conn,
            schema=schema,
            progress=progress,
            utils_db_url=None,
            use_published_identifiers=True,
        )
    else:
        raise ValueError(f"Unknown main product: {product}")
    return {"phase_seconds": time.monotonic() - started, "result": _reported(stats)}
