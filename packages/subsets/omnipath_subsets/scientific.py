"""Scientific product builders over committed normalized relational tables."""

from dataclasses import asdict, is_dataclass
import time
from psycopg2 import sql
from . import cosmos, exclusion_metabolites, metsigdb, network_views


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
    if product == "exclusion_metabolites":
        stats = exclusion_metabolites.build_exclusion_metabolites(conn, schema=schema)
    elif product == "network_views":
        stats = network_views.apply_all(
            conn,
            network_views.NETWORKS,
            registry_schema=schema,
            log=(lambda line: print(line, flush=True)) if progress else (lambda line: None),
        )
    elif product == "metsigdb":
        # main's original public build had an unqualified build_id lookup after
        # DDL reset the search path. Reapply the explicit target for isolation.
        metsigdb.ensure_membership_table(conn, schema=schema)
        _search_path(conn, schema)
        stamp = metsigdb.build_id(conn)
        loaded = []
        with conn.cursor() as cur:
            cur.execute("SELECT name FROM data_source")
            available = {row[0] for row in cur.fetchall()}
        for rule in metsigdb.RESOURCES:
            if rule.source_name in available and (
                not rule.hierarchy_source_name or rule.hierarchy_source_name in available
            ):
                _search_path(conn, schema)
                loaded.append(metsigdb.load_resource(conn, rule, stamp=stamp))
        with conn.cursor() as cur:
            cur.execute("SELECT count(*) FROM metsigdb_membership")
            rows = cur.fetchone()[0]
            # Main's publication requires statistics. VACUUM is deferred
            # until after the durable product checkpoint; ANALYZE stays atomic.
            cur.execute("ANALYZE metsigdb_membership")
        stats = {"build_id": stamp, "rows": rows, "resources": [_reported(row) for row in loaded]}
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
