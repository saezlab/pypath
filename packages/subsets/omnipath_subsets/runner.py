"""Selected product rebuilds over the aligned current-main PostgreSQL layout."""

from contextlib import closing
from copy import deepcopy
from time import perf_counter

import psycopg2
from psycopg2 import sql
from psycopg2.extras import Json

from dataclasses import dataclass

from .scientific import run_product
from .db import validate_schema

PRODUCTS = ("exclusion_metabolites", "metsigdb", "network_views", "cosmos")


@dataclass(frozen=True)
class BuildResult:
    release: str
    manifest_sha256: str
    schema: str
    products: dict
    phase_seconds: dict[str, float]


def _identity(connection, schema, expected=None):
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL("LOCK TABLE {}.parquet_release IN SHARE MODE").format(sql.Identifier(schema))
        )
        cursor.execute(
            sql.SQL(
                "SELECT version,manifest_sha256 FROM {}.parquet_release WHERE singleton"
            ).format(sql.Identifier(schema))
        )
        rows = cursor.fetchall()
    if len(rows) != 1:
        raise ValueError("Main subset build requires exactly one loaded release")
    identity = tuple(rows[0])
    if expected is not None and identity != expected:
        raise ValueError("Loaded release changed during main subset build")
    return identity


def _ready(connection, schema):
    with connection.cursor() as cursor:
        cursor.execute(sql.SQL("SELECT phase FROM {}.parquet_phase").format(sql.Identifier(schema)))
        phases = {row[0] for row in cursor.fetchall()}
    if not {"base", "derived"}.issubset(phases):
        raise ValueError("Main subset build requires committed base and shared derivations")


def _metadata(connection, schema):
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL("""CREATE TABLE IF NOT EXISTS {}.subset_build_metadata (
            product text PRIMARY KEY, release_id text NOT NULL,
            manifest_sha256 text NOT NULL, stats jsonb NOT NULL,
            built_at timestamptz NOT NULL DEFAULT now())""").format(sql.Identifier(schema))
        )


def _publish(connection, schema, product, identity, outcome, seconds):
    """Prepare all product markers in the same transaction as its tables."""
    _identity(connection, schema, identity)
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL("""INSERT INTO {}.subset_build_metadata
            (product,release_id,manifest_sha256,stats) VALUES(%s,%s,%s,%s)
            ON CONFLICT(product) DO UPDATE SET release_id=EXCLUDED.release_id,
            manifest_sha256=EXCLUDED.manifest_sha256,stats=EXCLUDED.stats,built_at=now()""").format(
                sql.Identifier(schema)
            ),
            (product, *identity, Json(outcome["result"])),
        )
        cursor.execute(
            sql.SQL("""INSERT INTO {}.parquet_phase(phase,seconds,result)
            VALUES(%s,%s,%s) ON CONFLICT(phase) DO UPDATE SET seconds=EXCLUDED.seconds,
            result=EXCLUDED.result,completed_at=now()""").format(sql.Identifier(schema)),
            (product, seconds, Json(outcome)),
        )
        cursor.execute(
            sql.SQL("""UPDATE {}.parquet_release SET
            phase_seconds=phase_seconds || jsonb_build_object(%s::text,%s::double precision),
            status=CASE WHEN (SELECT count(*) FROM {}.parquet_phase WHERE phase=ANY(%s))=%s
                        THEN 'complete' ELSE %s END WHERE singleton""").format(
                sql.Identifier(schema), sql.Identifier(schema)
            ),
            (product, seconds, list(PRODUCTS), len(PRODUCTS), product),
        )
    _identity(connection, schema, identity)


def _fresh_worker(connection, products):
    # MetSigDB's temp_buffers setting must precede any product temporary table.
    # This dedicated worker has never run shared derivations or an earlier run.
    if "metsigdb" in products:
        from .metsigdb.build import configure_temp_buffers

        with connection.cursor() as cursor:
            configure_temp_buffers(cursor)


def run_products(
    database_url,
    schema,
    *,
    owner,
    identity,
    products=PRODUCTS,
    checkpoint_products=True,
    resume=False,
    observer=None,
    progress=False,
    on_product_committed=None,
):
    """Run products while the caller holds the schema session advisory lock.

    Resume skips committed products; explicit rebuild executes every selected
    product. Suppress internal commits until product tables and both metadata
    markers agree. Observers run after durable commits while the lock is held.
    """

    results, timings = {}, {}
    _identity(owner, schema, identity)
    _ready(owner, schema)
    if resume:
        with owner.cursor() as cursor:
            cursor.execute(
                sql.SQL("SELECT phase FROM {}.parquet_phase").format(sql.Identifier(schema))
            )
            completed = {row[0] for row in cursor.fetchall()}
        products = tuple(product for product in products if product not in completed)
    owner.commit()  # Caller retains session lock, releases publication table lock.
    if checkpoint_products:
        for product in products:
            with closing(psycopg2.connect(database_url)) as worker:
                try:
                    _fresh_worker(worker, (product,))
                    _identity(worker, schema, identity)
                    _ready(worker, schema)
                    _metadata(worker, schema)
                    started = perf_counter()
                    if observer is not None:
                        observer("phase_start", phase=product)
                    outcome = run_product(worker, product, schema=schema, progress=progress)
                    seconds = round(perf_counter() - started, 6)
                    _publish(worker, schema, product, identity, outcome, seconds)
                    worker.commit()
                except BaseException:
                    worker.rollback()
                    raise
            results[product], timings[product] = outcome["result"], seconds
            if observer is not None:
                observer(
                    "phase_committed", phase=product, seconds=seconds, result=deepcopy(outcome)
                )
            if on_product_committed is not None:
                on_product_committed(
                    product, BuildResult(*identity, schema, deepcopy(results), timings.copy())
                )
    else:
        with closing(psycopg2.connect(database_url)) as worker:
            try:
                _fresh_worker(worker, products)
                _identity(worker, schema, identity)
                _ready(worker, schema)
                _metadata(worker, schema)
                for product in products:
                    started = perf_counter()
                    if observer is not None:
                        observer("phase_start", phase=product)
                    outcome = run_product(worker, product, schema=schema, progress=progress)
                    seconds = round(perf_counter() - started, 6)
                    _publish(worker, schema, product, identity, outcome, seconds)
                    results[product], timings[product] = outcome["result"], seconds
                worker.commit()
            except BaseException:
                worker.rollback()
                raise
    return BuildResult(*identity, schema, results, timings)


def acquire_schema_lock(connection, schema):
    with connection.cursor() as cursor:
        cursor.execute(
            "SELECT pg_try_advisory_lock(hashtextextended(%s,0))",
            ("omnipath:main-parquet:" + schema,),
        )
        if not cursor.fetchone()[0]:
            raise ValueError("Another migration owns this schema")
    connection.commit()


def build_subsets(
    database_url, schema, *, products=PRODUCTS, checkpoint_products=False, on_product_committed=None
):
    """Explicitly rebuild selected products from the committed current release."""
    validate_schema(schema)
    if not isinstance(products, (tuple, list)) or any(not isinstance(p, str) for p in products):
        raise ValueError("Products must be unique known products")
    if not products or len(set(products)) != len(products) or set(products) - set(PRODUCTS):
        raise ValueError("Products must be unique known products")
    if type(checkpoint_products) is not bool:
        raise ValueError("checkpoint_products must be a boolean")
    if on_product_committed is not None:
        if not callable(on_product_committed):
            raise ValueError("on_product_committed must be callable")
        if not checkpoint_products:
            raise ValueError("on_product_committed requires checkpoint_products=True")
    with closing(psycopg2.connect(database_url)) as owner:
        acquire_schema_lock(owner, schema)
        identity = _identity(owner, schema)
        return run_products(
            database_url,
            schema,
            owner=owner,
            identity=identity,
            products=products,
            checkpoint_products=checkpoint_products,
            on_product_committed=on_product_committed,
        )
