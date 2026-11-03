"""The writer reads entity records from an identity snapshot (spec section 1 integration)."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[2] / "resolver" / "tests"))
from identity_snapshot import ETHANOL, WATER, build_snapshot  # noqa: E402

from omnipath_build.silver import SilverExtractor  # noqa: E402
from omnipath_build.writer import ParquetWriter, _reference_records  # noqa: E402
from tables_fixture import read_tables  # noqa: E402
from omnipath_resolver.identity_runtime import IdentityRuntime  # noqa: E402
from omnipath_resolver.resolver import EntityResolver  # noqa: E402


@pytest.fixture
def snapshot(tmp_path, monkeypatch):
    monkeypatch.setenv("OMNIPATH_IDENTITY_CACHE", str(tmp_path / "cache"))
    return build_snapshot(tmp_path)


def test_reference_records_are_fetched_in_batches(snapshot, tmp_path):
    runtime = IdentityRuntime(snapshot, cache_dir=tmp_path / "c")
    try:
        ids = [f"inchikey:{WATER}", f"inchikey:{ETHANOL}", "entrez:7157"]
        got = dict(_reference_records(runtime, ids, batch_size=2))
        assert set(got) == set(ids)
        assert ("chebi", "CHEBI:15377") in {tuple(p) for p in got[ids[0]]["identifiers"]}
    finally:
        runtime.close()


def test_writer_attaches_resolver_aliases_from_an_identity_snapshot(snapshot, tmp_path):
    item = {
        "type": "small_molecule",
        "identifiers": [{"type": "chebi", "value": "CHEBI:15377"}],
    }
    extractor = SilverExtractor("fixture", "identity")
    extractor.process_record(item, item, "standalone", 0)
    resolver = EntityResolver(library_dir=snapshot, defer_aliases=True)
    writer = ParquetWriter(tmp_path / "output", library_dir=snapshot)
    try:
        assert isinstance(resolver.matcher.runtime, IdentityRuntime)
        writer.append_observations(extractor, resolver)
        tables = read_tables(writer.close()["files"])
    finally:
        resolver.close()
    (entity,) = tables["entity"]
    assert entity["reference_entity_key"] == f"inchikey:{WATER}"
    identifiers = tables.children("entity_identifier", entity)
    aliases = {(a["ns"], a["id"]) for a in identifiers if a["source"] == "resolver"}
    assert {("inchikey", WATER), ("chebi", "CHEBI:15377"), ("hmdb", "HMDB0002111")} <= aliases
    assert ("name", "water") in aliases
