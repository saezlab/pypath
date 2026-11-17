"""Where a subset product lives, how it is written and when it is refused."""

import hashlib
import json

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from omnipath_core import subset_layout as layout


def publish_release(root, version="2026.11", *, partial=False, resources=None):
    manifest = {
        "schema_version": 1,
        "version": version,
        "resources": resources or {"chebi": "1.0.0"},
    }
    if partial:
        manifest["partial_resources"] = True
    path = root / "releases" / f"{version}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table_writer(name="members.parquet", rows=3):
    def write(directory, context):
        pq.write_table(pa.table({"entity_key": [f"k{i}" for i in range(rows)]}), directory / name)

    return write


def build(root, product="exclusion_metabolites", writer=None, release="2026.11"):
    return layout.build_product(root, release, product, writer or table_writer())


def test_a_product_is_written_under_its_release_and_name(tmp_path):
    publish_release(tmp_path)
    build(tmp_path)
    directory = tmp_path / "subsets" / "2026.11" / "exclusion_metabolites"
    assert sorted(path.name for path in directory.iterdir()) == [
        "members.parquet",
        "subset_manifest.json",
    ]


def test_the_manifest_pins_the_release_and_describes_each_table(tmp_path):
    digest = publish_release(tmp_path)
    build(tmp_path, writer=table_writer(rows=5))
    directory = tmp_path / "subsets" / "2026.11" / "exclusion_metabolites"
    manifest = json.loads((directory / "subset_manifest.json").read_text())
    assert manifest["schema_version"] == 1
    assert manifest["product"] == "exclusion_metabolites"
    assert manifest["release"] == "2026.11"
    assert manifest["release_manifest_sha256"] == digest
    assert manifest["depends_on"] == {}
    assert manifest["partial"] is False
    assert manifest["built_at"]
    assert manifest["code_version"]
    assert manifest["tables"] == {
        "members.parquet": {
            "rows": 5,
            "sha256": hashlib.sha256((directory / "members.parquet").read_bytes()).hexdigest(),
        }
    }


def test_the_manifest_names_the_row_contract_and_keeps_what_the_writer_reports(tmp_path):
    publish_release(tmp_path)

    def writer(directory, context):
        table_writer()(directory, context)
        return {"rows": 3}

    manifest = build(tmp_path, writer=writer)
    assert manifest["row_contract"] == layout.ROW_CONTRACTS["exclusion_metabolites"]
    assert manifest["stats"] == {"rows": 3}
    stored = json.loads(
        (tmp_path / "subsets/2026.11/exclusion_metabolites/subset_manifest.json").read_text()
    )
    assert stored["row_contract"] == manifest["row_contract"]


def test_the_writer_works_in_a_temporary_directory_that_is_gone_afterwards(tmp_path):
    publish_release(tmp_path)
    seen = {}

    def writer(directory, context):
        seen["directory"] = directory
        seen["final_exists"] = (tmp_path / "subsets/2026.11/exclusion_metabolites").exists()
        table_writer()(directory, context)

    build(tmp_path, writer=writer)
    assert seen["final_exists"] is False
    assert not seen["directory"].exists()
    leftovers = [p.name for p in (tmp_path / "subsets" / "2026.11").iterdir()]
    assert leftovers == ["exclusion_metabolites"]


def test_a_failing_writer_leaves_nothing_behind(tmp_path):
    publish_release(tmp_path)

    def writer(directory, context):
        table_writer()(directory, context)
        raise RuntimeError("the query failed")

    with pytest.raises(RuntimeError, match="the query failed"):
        build(tmp_path, writer=writer)
    assert list((tmp_path / "subsets").rglob("*")) in ([], [tmp_path / "subsets" / "2026.11"])


def test_a_product_without_a_table_is_refused_and_not_published(tmp_path):
    publish_release(tmp_path)
    with pytest.raises(layout.SubsetError, match="no table"):
        build(tmp_path, writer=lambda directory, context: None)
    assert not (tmp_path / "subsets/2026.11/exclusion_metabolites").exists()


def test_a_table_that_is_not_readable_parquet_is_refused(tmp_path):
    publish_release(tmp_path)

    def writer(directory, context):
        (directory / "members.parquet").write_text("not parquet")

    with pytest.raises(layout.SubsetError, match="members.parquet"):
        build(tmp_path, writer=writer)
    assert not (tmp_path / "subsets/2026.11/exclusion_metabolites").exists()


def test_a_second_build_for_the_same_release_and_product_is_refused(tmp_path):
    publish_release(tmp_path)
    build(tmp_path)
    directory = tmp_path / "subsets/2026.11/exclusion_metabolites"
    before = (directory / "members.parquet").read_bytes()
    with pytest.raises(layout.SubsetExists):
        build(tmp_path, writer=table_writer(rows=9))
    assert (directory / "members.parquet").read_bytes() == before


def test_an_unpublished_release_is_refused(tmp_path):
    with pytest.raises(layout.SubsetError, match="2026.11"):
        build(tmp_path)


def test_an_unknown_product_is_refused(tmp_path):
    publish_release(tmp_path)
    with pytest.raises(layout.SubsetError, match="network_views"):
        build(tmp_path, product="network_views")


def test_a_release_pinning_samples_gives_a_partial_subset(tmp_path):
    publish_release(tmp_path, partial=True)
    manifest = build(tmp_path)
    assert manifest["partial"] is True


def test_metsigdb_needs_the_exclusion_metabolites_first(tmp_path):
    publish_release(tmp_path)
    with pytest.raises(layout.SubsetError, match="exclusion_metabolites"):
        build(tmp_path, product="metsigdb")


def test_metsigdb_records_the_digest_of_the_exclusion_manifest_it_read(tmp_path):
    publish_release(tmp_path)
    build(tmp_path)
    seen = {}

    def writer(directory, context):
        seen.update(context.dependencies)
        table_writer("metsigdb_membership.parquet")(directory, context)

    manifest = build(tmp_path, product="metsigdb", writer=writer)
    exclusion = tmp_path / "subsets/2026.11/exclusion_metabolites/subset_manifest.json"
    expected = hashlib.sha256(exclusion.read_bytes()).hexdigest()
    assert manifest["depends_on"] == {"exclusion_metabolites": expected}
    assert seen["exclusion_metabolites"] == exclusion.parent


def test_a_dependency_built_for_another_release_digest_is_refused(tmp_path):
    publish_release(tmp_path)
    build(tmp_path)
    publish_release(tmp_path, resources={"chebi": "2.0.0"})
    with pytest.raises(layout.SubsetError, match="release manifest"):
        build(tmp_path, product="metsigdb", writer=table_writer("metsigdb_membership.parquet"))


def test_a_subset_is_read_only_while_its_release_manifest_is_unchanged(tmp_path):
    publish_release(tmp_path)
    build(tmp_path)
    found = layout.read_subset(tmp_path, "2026.11", "exclusion_metabolites")
    assert found.tables == {"members.parquet": found.directory / "members.parquet"}
    assert found.manifest["product"] == "exclusion_metabolites"


def test_a_subset_of_a_changed_release_manifest_is_refused_when_read(tmp_path):
    publish_release(tmp_path)
    build(tmp_path)
    publish_release(tmp_path, resources={"chebi": "2.0.0"})
    with pytest.raises(layout.SubsetError, match="release manifest"):
        layout.read_subset(tmp_path, "2026.11", "exclusion_metabolites")


def test_a_subset_with_an_altered_table_is_refused_when_read(tmp_path):
    publish_release(tmp_path)
    build(tmp_path)
    pq.write_table(
        pa.table({"entity_key": ["other"]}),
        tmp_path / "subsets/2026.11/exclusion_metabolites/members.parquet",
    )
    with pytest.raises(layout.SubsetError, match="members.parquet"):
        layout.read_subset(tmp_path, "2026.11", "exclusion_metabolites")


def test_a_release_without_a_subset_reads_as_none(tmp_path):
    publish_release(tmp_path)
    assert layout.read_subset(tmp_path, "2026.11", "metsigdb") is None


def test_present_subsets_follow_the_dependency_order(tmp_path):
    publish_release(tmp_path)
    build(tmp_path)
    build(tmp_path, product="metsigdb", writer=table_writer("metsigdb_membership.parquet"))
    found = layout.present_subsets(tmp_path, "2026.11")
    assert [subset.manifest["product"] for subset in found] == ["exclusion_metabolites", "metsigdb"]
