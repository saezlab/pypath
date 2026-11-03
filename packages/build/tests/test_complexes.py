"""Composition identity and referential integrity across shard boundaries."""

from omnipath_build.complexes import composition_keys
from test_partition_contract import entity, run_rows
from omnipath_resolver import EntityResolver
from omnipath_build.silver import SilverExtractor
from omnipath_build.writer import ParquetWriter
from tables_fixture import read_tables


def complex_entity(name, members=()):
    return {
        "type": "macromolecular_complex",
        "identifiers": [{"type": "signor", "value": name}],
        "membership": [{"member": m} for m in members],
    }


def test_composition_order_duplicates_nested_conflicts_cycles():
    keys = composition_keys(
        set("abcdefg"),
        {
            "a": [{"P1", "P2"}],
            "b": [{"P2", "P1"}, {"P1", "P2"}],
            "c": [{"a", "P3"}],
            "d": [{"P1", "P2", "P3"}],
            "e": [{"P1"}, {"P2"}],
            "f": [{"g"}],
            "g": [{"f"}],
        },
    )
    assert keys["a"] == keys["b"]
    assert keys["c"] == keys["d"]
    assert not set("efg") & keys.keys()


def records():
    return [
        (
            "interaction",
            {"subject": complex_entity("A"), "predicate": "interacts_with", "object": entity("P3")},
        ),
        ("a", complex_entity("A", [entity("P1"), entity("P2"), entity("P1")])),
        ("b", complex_entity("B", [entity("P2"), entity("P1")])),
        ("different", complex_entity("C", [entity("P1"), entity("P3")])),
    ]


def check_output(tables):
    entities, relations, payloads = (
        tables["entity"],
        tables["relation"],
        tables["evidence_payloads"],
    )
    complexes = [e for e in entities if e["entity_type"] == "macromolecular_complex"]
    assert len(complexes) == 2
    assert {e["namespace"] for e in complexes} == {"complex"}

    def ids(row):
        return tables.children("entity_identifier", row)

    merged = next(e for e in complexes if {"A", "B"} <= {i["id"] for i in ids(e)})
    assert all(i["ns"] == "complex" for i in ids(merged) if i["is_canonical"])
    keys = {e["entity_key"] for e in entities}
    relkeys = {r["relation_key"] for r in relations}
    assert all(
        r["subject_entity_key"] in keys and r["object_entity_key"] in keys for r in relations
    )
    assert all(
        p["relation_key"] in relkeys if p["relation_key"] else p["entity_key"] in keys
        for p in payloads
    )
    interaction = next(r for r in relations if r["predicate"] == "interacts_with")
    assert merged["entity_key"] in (
        interaction["subject_entity_key"],
        interaction["object_entity_key"],
    )


def test_composition_finalizes_across_batches_and_shards(tmp_path):
    rows = records()
    tables = run_rows(tmp_path / "all", rows, len(rows))
    check_output(tables)
    assert run_rows(tmp_path / "split", rows[::-1], 1).comparable() == tables.comparable()
    resolver = EntityResolver(tmp_path / "absent")
    final = ParquetWriter(tmp_path / "final")
    try:
        for n, (locator, record) in enumerate(rows):
            shard = ParquetWriter(tmp_path / f"shard-{n}")
            extractor = SilverExtractor("fixture", "interactions")
            extractor.process_record(record, record, locator, 0)
            shard.append_observations(extractor, resolver)
            final.import_observation_shard(shard.seal_observation_shard())
        check_output(read_tables(final.close()["files"]))
    finally:
        resolver.close()


def test_members_use_resolved_accessions(tmp_path):
    from library_fixture import build_fixture_library

    library = build_fixture_library(tmp_path / "library")
    rows = [
        ("primary", complex_entity("A", [entity("P04637")])),
        ("secondary", complex_entity("B", [entity("Q15086")])),
    ]
    tables = run_rows(tmp_path / "out", rows, 1, library)
    assert run_rows(tmp_path / "all", rows[::-1], len(rows), library).comparable() == (
        tables.comparable()
    )
    entities, relations, payloads = (
        tables["entity"],
        tables["relation"],
        tables["evidence_payloads"],
    )
    complexes = [e for e in entities if e["entity_type"] == "macromolecular_complex"]
    assert len(complexes) == 1
    assert {"A", "B"} <= {i["id"] for i in tables.children("entity_identifier", complexes[0])}
    by_key = {e["entity_key"]: e for e in entities}
    assert len(relations) == 1 and relations[0]["evidence_count"] == 2
    assert relations[0]["object_reference_entity_key"] == "entrez:7157"
    for evidence in tables.children("relation_evidence", relations[0]):
        product = by_key[evidence["object_molecular_form"]["protein_entity_key"]]
        assert (product["namespace"], product["identifier"]) == ("uniprot", "P04637")
    assert len(payloads) == 2


def test_same_gene_products_and_gene_only_members_have_distinct_compositions(tmp_path):
    from library_fixture import build_fixture_library

    library = build_fixture_library(tmp_path / "library")
    gene = {"type": "gene", "identifiers": [{"type": "entrez", "value": "7157"}]}
    rows = [
        ("primary", complex_entity("A", [entity("P04637")])),
        ("other_product", complex_entity("B", [entity("A0A0U1RQF1")])),
        ("gene_only", complex_entity("C", [gene])),
    ]
    tables = run_rows(tmp_path / "out", rows, len(rows), library)
    assert run_rows(tmp_path / "split", rows[::-1], 1, library).comparable() == tables.comparable()
    entities, relations, payloads = (
        tables["entity"],
        tables["relation"],
        tables["evidence_payloads"],
    )
    complexes = [e for e in entities if e["entity_type"] == "macromolecular_complex"]
    assert len(complexes) == 3
    by_key = {e["entity_key"]: e for e in entities}
    products = set()
    for relation in relations:
        assert relation["predicate"] == "has_member"
        assert relation["subject_entity_key"] in by_key and relation["object_entity_key"] in by_key
        assert relation["object_reference_entity_key"] == "entrez:7157"
        form = tables.children("relation_evidence", relation)[0]["object_molecular_form"]
        if relation["object_type"] == "gene":
            assert form is None
        else:
            products.add(by_key[form["protein_entity_key"]]["identifier"])
        assert relation["evidence_count"] == 1
    assert products == {"P04637", "A0A0U1RQF1"}
    relation_keys = {r["relation_key"] for r in relations}
    assert all(
        p["relation_key"] in relation_keys if p["relation_key"] else p["entity_key"] in by_key
        for p in payloads
    )
