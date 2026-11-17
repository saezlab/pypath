"""The subset tables of a release are copied into the database it is served from."""

import json

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from omnipath_api.release_database import build_served_database, publish_newest
from omnipath_core import subset_layout
from omnipath_core.fixtures import write_dataset


@pytest.fixture
def root(tmp_path):
    root = write_dataset(tmp_path / "data")
    publish_newest(root, "1")
    return root


def write_product(root, product, tables):
    def writer(directory, context):
        for name, rows in tables.items():
            pq.write_table(pa.table(rows), directory / f"{name}.parquet")
        return {"tables": len(tables)}

    return subset_layout.build_product(root, "1", product, writer)


def serve(root):
    return build_served_database(root, "1", keep=None, memory_limit="1GB", threads=1)


def tables_of(root):
    with duckdb.connect(str(root / "serving" / "1.duckdb"), read_only=True) as db:
        return {row[0] for row in db.execute("SHOW TABLES").fetchall()}


def test_a_release_without_subsets_gets_no_subset_table(root):
    serve(root)
    found = tables_of(root)
    assert "subset_info" not in found
    assert not {"metsigdb_membership", "exclusion_metabolite"} & found


def test_each_present_subset_is_copied_with_its_rows(root):
    write_product(
        root, "exclusion_metabolites", {"exclusion_metabolite": {"entity_key": ["a", "b"]}}
    )
    write_product(
        root,
        "metsigdb",
        {"metsigdb_membership": {"set_source_id": ["s"]}},
    )
    serve(root)
    with duckdb.connect(str(root / "serving" / "1.duckdb"), read_only=True) as db:
        assert db.execute("SELECT entity_key FROM exclusion_metabolite ORDER BY 1").fetchall() == [
            ("a",),
            ("b",),
        ]
        assert db.execute("SELECT count(*) FROM metsigdb_membership").fetchone() == (1,)
        info = db.execute(
            "SELECT product, row_contract, partial FROM subset_info ORDER BY product"
        ).fetchall()
    assert info == [("exclusion_metabolites", 1, False), ("metsigdb", 2, False)]


def test_the_subset_is_recorded_with_the_digest_of_its_manifest(root):
    write_product(root, "exclusion_metabolites", {"exclusion_metabolite": {"entity_key": ["a"]}})
    record = serve(root)
    manifest = root / "subsets" / "1" / "exclusion_metabolites" / "subset_manifest.json"
    import hashlib

    digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
    with duckdb.connect(str(root / "serving" / "1.duckdb"), read_only=True) as db:
        assert db.execute("SELECT manifest_sha256 FROM subset_info").fetchone() == (digest,)
    assert record["subsets"] == {"exclusion_metabolites": digest}
    assert json.loads((root / "serving" / "1.json").read_text())["subsets"] == record["subsets"]


def test_a_subset_that_no_longer_matches_its_release_stops_the_build(root):
    write_product(root, "exclusion_metabolites", {"exclusion_metabolite": {"entity_key": ["a"]}})
    release = root / "releases" / "1.json"
    manifest = json.loads(release.read_text())
    manifest["created_at"] = "another day"
    release.write_text(json.dumps(manifest))
    with pytest.raises(subset_layout.SubsetError, match="release manifest"):
        serve(root)
    assert not (root / "serving" / "1.duckdb").exists()
