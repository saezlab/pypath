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


def test_release_database_holds_every_resource_with_global_ids(tmp_path):
    root = write_dataset(tmp_path / "data")
    release(root)
    target = tmp_path / "serving" / "1.duckdb"
    target.parent.mkdir()
    result = build_release_database(root, "1", target, memory_limit="1GB", threads=1)
    assert result["path"] == str(target) and result["resources"] == 3
    assert not [p for p in target.parent.iterdir() if p != target]  # no staging left
    with duckdb.connect(str(target), read_only=True) as db:
        info = db.execute(
            "SELECT release, database_format, query_layer, built_by FROM database_info"
        ).fetchone()
        assert info == ("1", DATABASE_FORMAT, QUERY_LAYER, "build_release_database")
        assert db.execute("SELECT * FROM op_resource_index ORDER BY idx").fetchall() == [
            ("chebi", 0), ("signor", 1), ("uniprot", 2)
        ]
        # Each table holds every resource's rows, with the global id of its rows.
        for table, column, global_id in [
            ("entity", "entity_id", "_eid"),
            ("entity_evidence", "entity_id", "_eid"),
            ("relation", "relation_id", "_rid"),
            ("relation_evidence", "relation_id", "_rid"),
        ]:
            rows = sum(
                pq.read_metadata(root / "resources" / r / "1" / f"{table}.parquet").num_rows
                for r in ("chebi", "signor", "uniprot")
            )
            assert db.execute(f"SELECT count(*) FROM {table}").fetchone() == (rows,), table
            assert db.execute(
                f"""SELECT count(*) FROM {table} t JOIN op_resource_index i ON i._resource = t._resource
                WHERE t.{global_id} <> (i.idx << 32) + t.{column}"""
            ).fetchone() == (0,), table
        # Raw payloads stay in their files on the server.
        tables = {r[0] for r in db.execute("SELECT table_name FROM duckdb_tables()").fetchall()}
        assert "evidence_payloads" not in tables
        assert db.execute("SELECT bool_or(has_evidence) FROM resources").fetchone() == (False,)
        # The derived tables and indexes.
        assert db.execute(
            "SELECT _rid FROM relation_key_lookup WHERE key_prefix = ?", [int(fixture_key("relation")[:15], 16)]
        ).fetchall() == [((1 << 32) + 0,)]
        assert {r[0] for r in db.execute("SELECT index_name FROM duckdb_indexes()").fetchall()} == {
            "relation_rid", "relation_evidence_rid"
        }
        assert db.execute("SELECT count(*) FROM signor.relation").fetchone() == (1,)


def test_release_database_refuses_incomplete_resources(tmp_path):
    root = write_dataset(tmp_path / "data")
    (root / "resources" / "uniprot" / "1" / "relation_endpoint.parquet").unlink()
    release(root)
    with pytest.raises(FileNotFoundError, match="uniprot/1 lacks relation_endpoint.parquet"):
        build_release_database(root, "1", tmp_path / "1.duckdb")
    assert not (tmp_path / "1.duckdb").exists()
