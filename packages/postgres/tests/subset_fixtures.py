"""A release schema holding what loading the subsets needs of a loaded release."""

from contextlib import closing
import hashlib
import json
import uuid

import duckdb
import psycopg2
from psycopg2 import sql
import pytest

from omnipath_core.fixtures import (
    METSIGDB_METABOLITES,
    METSIGDB_SETS,
    metsigdb_key,
    write_dataset,
    write_metsigdb_subsets,
)
from omnipath_postgres.projection import AP_UUID_MACRO, entity_id_sql

VERSION = "1"
MEMBERSHIP_ROWS = sum(len(s[6]) for s in METSIGDB_SETS)


def canonical_digest(manifest):
    text = json.dumps(manifest, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(text.encode()).hexdigest()


def entity_uuid(key):
    con = duckdb.connect()
    con.execute(AP_UUID_MACRO)
    return str(con.execute(f"SELECT {entity_id_sql('?')}", [key]).fetchone()[0])


@pytest.fixture
def root(tmp_path):
    root = write_dataset(tmp_path / "data")
    from omnipath_api.release_database import publish_newest

    publish_newest(root, VERSION)
    write_metsigdb_subsets(root, VERSION)
    return root


@pytest.fixture
def release(postgres_dsn, root):
    """A schema holding what loading the subsets needs of a loaded release."""
    schema = "subsets_" + uuid.uuid4().hex
    manifest = json.loads((root / "releases" / f"{VERSION}.json").read_text())
    names = {s[1] for s in METSIGDB_SETS} | set(METSIGDB_METABOLITES)
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        try:
            with connection.cursor() as cur:
                cur.execute(sql.SQL("CREATE SCHEMA {}").format(sql.Identifier(schema)))
                cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
                cur.execute(
                    """
                    CREATE TABLE parquet_release(singleton boolean PRIMARY KEY DEFAULT true,
                      version text, manifest_sha256 text, manifest jsonb, status text,
                      counts jsonb DEFAULT '{}', phase_seconds jsonb DEFAULT '{}',
                      compatibility jsonb DEFAULT '{}');
                    CREATE TABLE parquet_phase(phase text PRIMARY KEY,
                      completed_at timestamptz DEFAULT now(), seconds double precision,
                      result jsonb);
                    CREATE TABLE parquet_resource(resource text PRIMARY KEY, version text,
                      manifest_sha256 text, manifest jsonb);
                    CREATE TABLE entity(entity_id uuid PRIMARY KEY, label text,
                      canonical_identifier text);
                    CREATE TABLE entity_identifier_lookup(entity_id uuid, identifier_id uuid);
                    CREATE TABLE identifier_evidence(identifier_id uuid PRIMARY KEY,
                      identifier_type_id bigint, value text);
                    CREATE TABLE vocab_identifier_type(identifier_type_id bigint PRIMARY KEY,
                      name text);
                    """
                )
                cur.execute(
                    "INSERT INTO parquet_release(version, manifest_sha256, manifest, status) "
                    "VALUES (%s, %s, %s, 'derived')",
                    [VERSION, canonical_digest(manifest), json.dumps(manifest)],
                )
                cur.executemany(
                    "INSERT INTO parquet_phase(phase, seconds, result) VALUES (%s, 0, '{}')",
                    [("base",), ("derived",)],
                )
                cur.executemany(
                    "INSERT INTO entity VALUES (%s, %s, %s)",
                    [(entity_uuid(metsigdb_key(name)), name, name) for name in names],
                )
            connection.commit()
            yield connection, schema, postgres_dsn
        finally:
            connection.rollback()
            with connection.cursor() as cur:
                cur.execute(sql.SQL("DROP SCHEMA {} CASCADE").format(sql.Identifier(schema)))
            connection.commit()
