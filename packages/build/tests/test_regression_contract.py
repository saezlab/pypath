import json
from pathlib import Path


from library_fixture import build_fixture_library
from omnipath_resolver import EntityResolver
from omnipath_build.silver import RawEntityObservation, SilverExtractor
from omnipath_build.writer import ParquetWriter
from test_partition_contract import records, entity
from test_resolver import obs
from tables_fixture import read_tables

# Nested lists of the frozen contract live in their own tables now.
NESTED = {"identifiers", "annotations", "evidence"}


def test_resolution_contract_across_batches(tmp_path):
    from dataclasses import asdict

    library = build_fixture_library(tmp_path / "library")
    fixture = json.loads((Path(__file__).parent / "fixtures/resolution_contract.json").read_text())
    rows = [RawEntityObservation(**r) for r in fixture["observations"]]
    resolver = EntityResolver(library, defer_aliases=True)
    try:
        baseline = {
            k: asdict(v)
            for k, v in resolver.resolve_entities({o.entity_key: o for o in rows}).items()
        }
        # Chemical identity under the identity rules (resolution.md): a frozen oracle.
        for observation in rows:
            if observation.entity_type == "chemical_entity":
                golden = fixture["matches"][observation.entity_key]
                assert {k: baseline[observation.entity_key][k] for k in golden} == golden
                assert baseline[observation.entity_key]["protein_identifier"] is None
        for items in (rows, rows[::-1], rows[::2]):
            for offset in range(0, len(items), 17):
                batch = {o.entity_key: o for o in items[offset : offset + 17]}
                actual = resolver.resolve_entities(batch)
                assert {k: asdict(v) for k, v in actual.items()} == {k: baseline[k] for k in batch}
        assert resolver.resolution_stats()["input_entities"] == len(rows)
    finally:
        resolver.close()


def test_flat_pipeline_matches_nested_contract_across_partitions(tmp_path):
    library = build_fixture_library(tmp_path / "library")
    rows = records() + [
        (
            "matched",
            {"subject": entity("P04637"), "predicate": "affects", "object": entity("P02340")},
        )
    ]
    rows += [
        (
            "ontology",
            {
                "type": "ontology_class",
                "identifiers": [{"type": "go", "value": "GO:1"}],
                "annotations": [{"term": "subclass_of", "value": "GO:2"}],
            },
        )
    ]
    golden = json.loads((Path(__file__).parent / "fixtures/resource_contract.json").read_text())
    baseline = None
    for name, source, size in [("one", rows, 1), ("all", rows[::-1], len(rows))]:
        writer = ParquetWriter(tmp_path / name, library_dir=library)
        resolver = EntityResolver(library, defer_aliases=True)
        try:
            for start in range(0, len(source), size):
                extractor = SilverExtractor("fixture", "interactions")
                for locator, record in source[start : start + size]:
                    extractor.process_record(
                        record, raw_payload=record, row_id=locator, row_number=0
                    )
                writer.append_observations(extractor, resolver)
            resolver.close()
            tables = read_tables(writer.close()["files"])
            if baseline is None:
                baseline = tables.comparable()
            else:
                assert tables.comparable() == baseline
            entities, relations, payloads = (
                tables["entity"],
                tables["relation"],
                tables["evidence_payloads"],
            )
            # Existing native entities and statements retain their exact frozen fields.
            native = {
                e["entity_key"]
                for e in golden[0]
                if e["identifier"] in {"P1", "P2", "P3", "P4", "GO:1"}
            }
            legacy_entities = {e["entity_key"]: e for e in golden[0] if e["entity_key"] in native}
            for row in entities:
                if row["entity_key"] in native:
                    expected = legacy_entities[row["entity_key"]]
                    expected = {k: v for k, v in expected.items() if k not in NESTED}
                    assert {k: row[k] for k in expected} == expected
            legacy_relations = {
                r["relation_key"]: r
                for r in golden[1]
                if r["subject_entity_key"] in native and r["object_entity_key"] in native
            }
            for row in relations:
                if row["relation_key"] in legacy_relations:
                    expected = legacy_relations[row["relation_key"]]
                    scalars = {k: v for k, v in expected.items() if k not in NESTED}
                    assert {k: row[k] for k in scalars} == scalars
                    evidence = tables.children("relation_evidence", row)
                    assert [
                        {k: ev[k] for k in old}
                        for ev, old in zip(evidence, expected["evidence"], strict=True)
                    ] == expected["evidence"]
            assert len(entities) == 9 and len(relations) == 6 and len(payloads) == 10
            assert {
                r["relation_key"] for r in relations if r["relation_key"] in legacy_relations
            } == set(legacy_relations)
            by_key = {e["entity_key"]: e for e in entities}
            matched = next(
                r
                for r in relations
                if tables.children("relation_evidence", r)[0]["row_id"] == "matched"
            )
            assert (
                matched["subject_reference_entity_key"],
                matched["object_reference_entity_key"],
            ) == ("entrez:7157", "entrez:22059")
            assert (
                by_key[matched["subject_entity_key"]]["namespace"],
                by_key[matched["subject_entity_key"]]["identifier"],
            ) == ("entrez", "7157")
            evidence = tables.children("relation_evidence", matched)[0]
            assert (
                by_key[evidence["subject_molecular_form"]["protein_entity_key"]]["identifier"]
                == "P04637"
            )
            assert (
                by_key[evidence["object_molecular_form"]["protein_entity_key"]]["identifier"]
                == "P02340"
            )
            assert matched["subject_type"] == matched["object_type"] == "protein"
            assert matched["evidence_count"] == 1
            assert all(
                r["subject_entity_key"] in by_key and r["object_entity_key"] in by_key
                for r in relations
            )
            relation_keys = {r["relation_key"] for r in relations}
            assert all(
                p["relation_key"] in relation_keys
                if p["relation_key"]
                else p["entity_key"] in by_key
                for p in payloads
            )
            raw_by_locator = dict(rows)
            assert all(
                json.loads(p["payload_json"]) == raw_by_locator[p["row_id"]] for p in payloads
            )
        finally:
            resolver.close()


def test_molecular_resolution_oracle_preserves_conflicts_and_source_types(tmp_path):
    library = build_fixture_library(tmp_path / "library")
    observations = {
        "gene": obs("gene", "gene", "entrez", "7157"),
        "gene_as_protein": obs("gene_as_protein", "protein", "entrez", "7157"),
        "primary": obs("primary", "protein", "uniprot", "P04637"),
        "other_product": obs("other_product", "protein", "uniprot", "A0A0U1RQF1"),
        "secondary": obs("secondary", "protein", "uniprot", "Q15086"),
        "multi": obs("multi", "protein", "uniprot", "P0DP23"),
        "conflict": obs("conflict", "protein", "uniprot", "P04637", [("entrez", "805")]),
    }
    expected = {
        "gene": ("entrez", "7157", None, "resolved", {"entrez:7157"}),
        "gene_as_protein": ("entrez", "7157", None, "resolved", {"entrez:7157"}),
        "primary": ("entrez", "7157", "P04637", "resolved", {"entrez:7157"}),
        "other_product": ("entrez", "7157", "A0A0U1RQF1", "resolved", {"entrez:7157"}),
        "secondary": ("entrez", "7157", "P04637", "resolved", {"entrez:7157"}),
        "multi": (
            "uniprot",
            "P0DP23",
            "P0DP23",
            "ambiguous",
            {"entrez:801", "entrez:805", "entrez:808"},
        ),
        "conflict": ("uniprot", "P04637", "P04637", "conflict", {"entrez:7157", "entrez:805"}),
    }
    resolver = EntityResolver(library, defer_aliases=True)
    try:
        for keys in (list(observations), list(reversed(observations))):
            for key in keys:
                targets = resolver.resolve_entity_targets({key: observations[key]})[key]
                assert len(targets) == 1
                result = targets[0]
                assert (
                    result.canonical_namespace,
                    result.canonical_identifier,
                    result.protein_identifier,
                    result.gene_mapping_status,
                    set(result.gene_candidates),
                ) == expected[key]
                assert result.entity_type == observations[key].entity_type
                assert result.matched
    finally:
        resolver.close()
