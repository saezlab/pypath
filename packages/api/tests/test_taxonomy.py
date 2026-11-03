"""Taxonomy labels are offline, release-pinned and never alter identity/filter IDs."""

import io
import tarfile

import pytest
from fastapi.testclient import TestClient

from omnipath_core.taxonomy import taxon_names
from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_api.queries.store.inventory import ReleaseStore
from omnipath_core.fixtures import rewrite_resource


def dump(path, scientific="Arabidopsis thaliana"):
    path.parent.mkdir(parents=True, exist_ok=True)
    files = {
        "names.dmp": f"3702 | {scientific} | | scientific name |\n9606 | Homo sapiens | | scientific name |\n9606 | human | | genbank common name |\n",
        "merged.dmp": "999 | 998 |\n998 | 3702 |\n",
    }
    with tarfile.open(path, "w:gz") as archive:
        for name, text in files.items():
            data = text.encode()
            member = tarfile.TarInfo(name)
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))


def test_names_common_scientific_merged_and_unknown(tmp_path):
    archive = tmp_path / "dump.tar.gz"
    dump(archive)
    rows = {r["taxon_id"]: r for r in taxon_names(archive, {"3702", "9606", "999", "123"})}
    assert rows["9606"]["common_name"] == "human"
    assert rows["3702"]["scientific_name"] == "Arabidopsis thaliana"
    assert rows["999"]["current_taxon_id"] == "3702"
    assert rows["999"]["scientific_name"] == "Arabidopsis thaliana"
    assert rows["123"]["scientific_name"] is None


def test_release_names_entities_facets_and_cached_relations(tmp_path):
    folder = tmp_path / "resources/test/1"
    folder.mkdir(parents=True)
    entities = [
        dict(
            entity_key="entity",
            entity_type="protein",
            namespace="test",
            identifier="p",
            label="Protein",
            taxon="3702",
            identifiers=[],
            annotations=[],
        )
    ]
    rewrite_resource(folder / "entity.parquet", entities=entities)
    relations = [
        dict(
            relation_key="relation",
            subject_entity_key="entity",
            object_entity_key="entity",
            subject_type="protein",
            object_type="protein",
            predicate="interacts_with",
            taxon="3702",
            category="interaction",
            sources=["test"],
            evidence_count=1,
        )
    ]
    rewrite_resource(folder / "relation.parquet", relations=relations)
    archive = tmp_path / "references/taxonomy/taxdump.tar.gz"
    store = ReleaseStore(tmp_path)
    store.publish(dict(schema_version=1, version="0", resources={"test": "1"}))
    manifests = []
    for version, name in [("1", "Arabidopsis thaliana"), ("2", "Updated name")]:
        dump(archive, name)
        manifests.append(
            store.publish(dict(schema_version=1, version=version, resources={"test": "1"}))
        )
    for manifest in manifests:
        taxonomy = (
            tmp_path
            / "references/taxonomy"
            / manifest["references"]["taxonomy"]["version"]
            / "taxonomy.parquet"
        )
        assert taxonomy.stat().st_mode & 0o777 == 0o644
    assert (
        manifests[0]["references"]["taxonomy"]["version"]
        != manifests[1]["references"]["taxonomy"]["version"]
    )
    engine = ParquetQueries(tmp_path)
    # Unfiltered relation facets share cached raw counts across these releases.
    for version, name in [
        ("1", "Arabidopsis thaliana"),
        ("2", "Updated name"),
        ("1", "Arabidopsis thaliana"),
        ("latest", "Updated name"),
        ("0", None),
    ]:
        with engine.release_scope(version):
            assert engine.get_entity_details("entity")["entity"]["taxonomyName"] == name
            for facets in (
                engine.get_scoped_entity_facets({}),
                engine.get_scoped_relation_facets({}),
            ):
                taxon = next(f for f in facets if f["facetName"] == "taxonomy_id")
                assert taxon["facetValue"] == "3702"
                assert taxon["facetLabel"] == name
    with TestClient(create_app(queries=engine)) as client:
        assert (
            client.get("/entities/entity?release=1").json()["entity"]["taxonomyName"]
            == "Arabidopsis thaliana"
        )
        response = client.post("/entities/scoped-facets?release=1", json={})
        assert any(row.get("facetLabel") == "Arabidopsis thaliana" for row in response.json())
    reference = manifests[0]["references"]["taxonomy"]["version"]
    (tmp_path / "references/taxonomy" / reference / "taxonomy.parquet").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="checksum"):
        store.validate_files(manifests[0])


def test_latest_lookup_covers_new_taxa_without_changing_pinned_names(tmp_path):
    from omnipath_core.taxonomy_lookup import prepare, read_lookup

    archive = tmp_path / "references/taxonomy/taxdump.tar.gz"
    dump(archive)
    path = prepare(tmp_path)
    names = read_lookup(str(path))
    assert names.get("9606") == "human"
    assert names.get("999") == "Arabidopsis thaliana"
    assert names.get("123") is None
    assert prepare(tmp_path) == path
    engine = ParquetQueries(tmp_path)
    with engine.release_scope("latest"):
        assert engine._taxonomy_names().get("3702") == "Arabidopsis thaliana"
        first_version = engine._taxonomy_cache_version()
    # A pinned scope without a taxonomy snapshot does not inherit moving labels.
    token = engine._release_scope.set({"version": "old", "resources": {}})
    try:
        assert engine._taxonomy_names().get("3702") is None
    finally:
        engine._release_scope.reset(token)
    dump(archive, "Updated name")
    prepare(tmp_path)
    with engine.release_scope("latest"):
        assert engine._taxonomy_names().get("3702") == "Updated name"
        assert engine._taxonomy_cache_version() != first_version
