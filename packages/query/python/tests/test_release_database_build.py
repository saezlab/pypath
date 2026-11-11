"""The database of a published release, as the server writes it for the API."""

import json

import duckdb
import pyarrow.parquet as pq
import pytest

from omnipath_client.database import DATABASE_FORMAT, QUERY_LAYER, build_release_database
from omnipath_core.fixtures import fixture_key, write_dataset


def release(root, version="1", resources=None):
    resources = resources or {"chebi": "1", "signor": "1", "uniprot": "1"}
    (root / "releases").mkdir(exist_ok=True)
    (root / "releases" / f"{version}.json").write_text(
        json.dumps(dict(schema_version=1, version=version, resources=resources))
    )


def test_release_database_holds_merged_numbered_tables(tmp_path):
    root = write_dataset(tmp_path / "data")
    release(root)
    target = tmp_path / "serving" / "1.duckdb"
    target.parent.mkdir()
    result = build_release_database(root, "1", target, memory_limit="1GB", threads=1)
    assert result["path"] == str(target) and result["resources"] == 3
    assert not [p for p in target.parent.iterdir() if p != target]  # no staging left
    files = root / "resources"

    def published(table, column):
        return pq.ParquetDataset([files / r / "1" / f"{table}.parquet" for r in ("chebi", "signor", "uniprot")]) \
            .read(columns=[column]).column(column).to_pylist()

    with duckdb.connect(str(target), read_only=True) as db:
        info = db.execute(
            "SELECT release, database_format, query_layer, built_by FROM database_info"
        ).fetchone()
        assert info == ("1", DATABASE_FORMAT, QUERY_LAYER, "build_release_database")
        # One row per entity and per relation key, numbered in key order.
        keys = sorted(set(published("entity", "entity_key")))
        assert db.execute("SELECT entity_key FROM entity ORDER BY entity_n").fetchall() == [(k,) for k in keys]
        assert db.execute("SELECT min(entity_n), max(entity_n) FROM entity").fetchone() == (1, len(keys))
        relation_keys = sorted(set(published("relation", "relation_key")))
        assert db.execute("SELECT relation_key FROM relation ORDER BY relation_n").fetchall() == [
            (k,) for k in relation_keys
        ]
        assert db.execute(
            "SELECT resources FROM relation WHERE relation_key = ?", [fixture_key("relation")]
        ).fetchone() == (["signor"],)
        # Every evidence record, with the resource stating it.
        for table in ("entity_evidence", "relation_evidence"):
            assert db.execute(f"SELECT count(*) FROM {table}").fetchone() == (
                len(published(table, "ordinal")),
            ), table
        # No resource tables, schemas or scratch tables are left.
        tables = {r[0] for r in db.execute("SELECT table_name FROM duckdb_tables()").fetchall()}
        assert not tables & {"entity_numbers", "relation_merged", "entity_copies", "evidence_payloads"}
        assert db.execute("SELECT count(*) FROM duckdb_schemas() WHERE schema_name = 'signor'").fetchone() == (0,)
        assert db.execute("SELECT bool_or(has_evidence) FROM resources").fetchone() == (False,)
        assert {r[0] for r in db.execute("SELECT index_name FROM duckdb_indexes()").fetchall()} == {
            "entity_n", "entity_key", "entity_evidence_n", "entity_annotation_n", "relation_n", "relation_evidence_n",
        }
        # The scratch database's tables took no space in the database.
        assert db.execute("SELECT free_blocks FROM pragma_database_size()").fetchone() == (0,)


def test_release_database_refuses_incomplete_resources(tmp_path):
    root = write_dataset(tmp_path / "data")
    (root / "resources" / "uniprot" / "1" / "relation_endpoint.parquet").unlink()
    release(root)
    with pytest.raises(FileNotFoundError, match="uniprot/1 lacks relation_endpoint.parquet"):
        build_release_database(root, "1", tmp_path / "1.duckdb")
    assert not (tmp_path / "1.duckdb").exists()
