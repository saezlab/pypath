"""Serving contract and evidence preservation through the production writer."""

import json
import pyarrow.parquet as pq

from omnipath_resolver import EntityResolver
from omnipath_build.silver import SilverExtractor
from omnipath_build.writer import ParquetWriter
from omnipath_core.schema import PAYLOAD_SCHEMA, PUBLISHED_TABLES, SERVING_TABLES
from tables_fixture import read_tables
from test_partition_contract import entity, run_rows


def shape(schema):
    return [(f.name, str(f.type).replace("element:", "item:")) for f in schema]


def test_evidence_occurrences_survive_batches(tmp_path):
    forward = {
        "subject": entity("P1"),
        "predicate": "affects",
        "object": entity("P2"),
        "annotations": [
            {"term": "object_direction_qualifier", "value": "increased"},
            {"term": "description", "value": "first"},
        ],
    }
    later = {
        **forward,
        "annotations": [forward["annotations"][0], {"term": "description", "value": "later"}],
    }
    opposite = {
        **forward,
        "annotations": [{"term": "object_direction_qualifier", "value": "decreased"}],
    }
    source = [("1", forward), ("2", later), ("1", forward), ("3", opposite)]
    tables = run_rows(tmp_path / "all", source, len(source))
    for name, rows in [("split", source), ("reverse", source[::-1])]:
        assert run_rows(tmp_path / name, rows, 1).comparable() == tables.comparable()
    relations, payloads = tables["relation"], tables["evidence_payloads"]
    assert len(relations) == 2 and len(payloads) == 4
    merged = next(r for r in relations if r["sign"] == 1)
    evidence = tables.children("relation_evidence", merged)
    assert merged["evidence_count"] == len(evidence) == 3
    assert [e["row_id"] for e in evidence].count("1") == 2
    annotations = tables.children("relation_annotation", merged)
    assert {a["value"] for a in annotations if a["term"] == "description"} == {"first", "later"}
    assert any(any(a["value"] == "later" for a in e["annotations"]) for e in evidence)


def test_empty_writer_keeps_all_schemas(tmp_path):
    writer = ParquetWriter(tmp_path)
    result = writer.close()
    assert set(result["rows"].values()) == {0}
    for name, schema in {**PUBLISHED_TABLES, **SERVING_TABLES}.items():
        assert shape(pq.read_schema(result["files"][name])) == shape(schema)
    assert not (tmp_path / ".bulk").exists()


def test_aliases_and_label_merge_across_batches(tmp_path):
    source = [("first", entity("P1")), ("later", entity("P1", label="CDKN2A", aliases=("CDKN2A",)))]
    tables = run_rows(tmp_path, source, 1)
    (row,) = tables["entity"]
    assert row["label"] == "CDKN2A"
    identifiers = tables.children("entity_identifier", row)
    assert ("genesymbol", "CDKN2A") in {(i["ns"], i["id"]) for i in identifiers}
    assert row["identifier_count"] == len(identifiers)


def test_nonempty_output_schemas_and_payload_links(tmp_path):
    source = [("1", {"subject": entity("P1"), "predicate": "affects", "object": entity("P1")})]
    tables = run_rows(tmp_path, source, 1)
    entities, relations, payloads = (
        tables["entity"],
        tables["relation"],
        tables["evidence_payloads"],
    )
    assert (len(entities), len(relations), len(payloads)) == (1, 1, 1)
    assert payloads[0]["relation_key"] == relations[0]["relation_key"]
    assert relations[0]["subject_entity_key"] == entities[0]["entity_key"]
    # A self-loop counts as one relation of its entity.
    assert entities[0]["relation_count"] == 1
    for name, schema in {**PUBLISHED_TABLES, **SERVING_TABLES}.items():
        assert shape(pq.read_schema(tmp_path / f"{name}.parquet")) == shape(schema)


def test_payload_arrow_batches_preserve_shared_source_rows(tmp_path):
    writer = ParquetWriter(tmp_path)
    resolver = EntityResolver()
    extractor = SilverExtractor("fixture", "test")
    raw = {"source_text": "é" * 200}
    for i in range(1030):
        extractor.process_record(
            {"subject": entity("P1"), "predicate": "affects", "object": entity("P2")},
            raw,
            str(i),
            i,
        )
    try:
        writer.append_observations(extractor, resolver)
        writer.close()
    finally:
        resolver.close()
    file = pq.ParquetFile(writer.payload_path)
    assert file.metadata.num_row_groups == 1
    rows = file.read().to_pylist()
    assert len(rows) == 1030
    assert {r["row_id"] for r in rows} == {str(i) for i in range(1030)}
    assert all(json.loads(r["payload_json"]) == raw for r in rows)
    assert writer.metrics["payload_rows"] == 1030


def test_large_shared_payloads_are_bounded_before_arrow_conversion(tmp_path, monkeypatch):
    writer = ParquetWriter(tmp_path)
    resolver = EntityResolver()
    extractor = SilverExtractor("fixture", "test")
    import random

    raw = {"text": random.Random(42).randbytes(200000).hex() + "é"}
    for i in range(40):
        extractor.process_record(
            {"subject": entity("P1"), "predicate": "affects", "object": entity("P2")},
            raw,
            str(i),
            i,
        )
    sizes = []
    original = writer._input

    def capture(name, rows, schema):
        if name == "input_payloads":
            assert all("payload_json" not in r for r in rows)
            sizes.append(len(rows))
        return original(name, rows, schema)

    monkeypatch.setattr(writer, "_input", capture)
    try:
        writer.append_observations(extractor, resolver)
        writer.close()
    finally:
        resolver.close()
    assert sizes == [40]
    assert writer.payload_path.stat().st_size < 500000
    assert (
        "RLE_DICTIONARY"
        in pq.ParquetFile(writer.payload_path).metadata.row_group(0).column(4).encodings
    )
    rows = pq.read_table(writer.payload_path).to_pylist()
    assert len(rows) == 40
    assert {r["row_id"] for r in rows} == {str(i) for i in range(40)}
    assert all(json.loads(r["payload_json"]) == raw for r in rows)


def test_payload_batches_bound_unique_bytes_and_preserve_nulls():
    from omnipath_build.writer import _payload_batches, _payload_dictionary

    rows = [{"payload_json": value} for value in ("é" * 100, "é" * 100, None, "x" * 1000, "last")]
    batches = list(_payload_batches(rows, max_bytes=1500))
    assert [row for batch in batches for row in batch] == rows
    assert len(batches) == 3
    for batch in batches:
        encoded = _payload_dictionary(batch)
        assert encoded.to_pylist() == [row["payload_json"] for row in batch]
    assert len(_payload_dictionary(rows).dictionary) == 3


def test_oversized_payload_references_share_dictionary(tmp_path):
    import random
    import pyarrow as pa
    from omnipath_build.writer import _payload_batches, _payload_dictionary, PAYLOAD_STORAGE_SCHEMA

    # Tiny limit exercises the same condition as a >16 MiB production JSON.
    value = random.Random(3).randbytes(5000).hex()
    rows = [
        dict(relation_key=str(i), entity_key=None, source="x", row_id=str(i), payload_json=value)
        for i in range(20)
    ]
    batches = list(_payload_batches(rows, max_bytes=8000))
    assert len(batches) == 1
    path = tmp_path / "oversized.parquet"
    with pq.ParquetWriter(path, PAYLOAD_STORAGE_SCHEMA, store_schema=False) as writer:
        for batch in batches:
            table = pa.Table.from_pylist(batch, schema=PAYLOAD_SCHEMA.remove(4))
            writer.write_table(table.append_column("payload_json", _payload_dictionary(batch)))
    assert pq.read_table(path).to_pylist() == rows
    assert path.stat().st_size < len(value) * 3
    assert len(list(_payload_batches(rows, max_bytes=8000, max_rows=7))) == 3
    extra = dict(rows[0], payload_json="different oversized" * 1000)
    assert len(list(_payload_batches(rows + [extra], max_bytes=8000))) == 2


def test_private_shards_match_single_writer_with_overlapping_event_ids(tmp_path):
    source = [
        (
            "1",
            {
                "subject": entity("P1"),
                "predicate": "affects",
                "object": entity("P2"),
                "annotations": [{"term": "description", "value": "first"}],
            },
        ),
        (
            "2",
            {
                "subject": entity("P1"),
                "predicate": "affects",
                "object": entity("P2"),
                "annotations": [{"term": "description", "value": "second"}],
            },
        ),
        (
            "1",
            {
                "subject": entity("P1"),
                "predicate": "affects",
                "object": entity("P2"),
                "annotations": [{"term": "description", "value": "first"}],
            },
        ),
    ]
    expected = run_rows(tmp_path / "single", source, 3).comparable()
    shards = []
    for i, (row_id, record) in enumerate(source):
        resolver = EntityResolver()
        writer = ParquetWriter(tmp_path / str(i))
        extractor = SilverExtractor("fixture", "interactions")
        extractor.process_record(record, record, row_id, i)
        writer.append_observations(extractor, resolver)
        resolver.close()
        shards.append(writer.seal_observation_shard())
    writer = ParquetWriter(tmp_path / "merged")
    for shard in shards[::-1]:
        writer.import_observation_shard(shard)
    # Shards hold identical work in another order: every table matches the single writer.
    assert read_tables(writer.close()["files"]).comparable() == expected


def test_chemical_names_stay_identifiers_and_group_by_structure(tmp_path):
    chemical = {
        "type": "chemical_entity",
        "identifiers": [
            {"type": "inchikey", "value": "QNAYBMKLOCPYGJ-UQEXSWPGSA-N"},
            {"type": "name", "value": "alanine-d7"},
            {"type": "name", "value": "Alanine"},
        ],
    }
    tables = run_rows(tmp_path, [("one", chemical)], 1)
    (row,) = tables["entity"]
    names = {i["id"] for i in tables.children("entity_identifier", row)}
    assert {"alanine-d7", "Alanine"} <= names
    # The structure group is the first InChIKey block of the reference.
    assert row["reference_entity_key"] == "inchikey:QNAYBMKLOCPYGJ-UQEXSWPGSA-N"
    assert row["group_connectivity"] == "QNAYBMKLOCPYGJ"
    groups = {(g["kind"], g["group_key"]) for g in tables["entity_group"]}
    assert ("connectivity", "QNAYBMKLOCPYGJ") in groups
    terms = {t["term"] for t in tables["entity_term"]}
    assert {"alanine", "alanine-d7", "qnaybmklocpygj-uqexswpgsa-n"} <= terms


def test_outputs_are_sorted_by_key_across_buckets(tmp_path):
    # Enough distinct keys to span many key-prefix buckets.
    source = [
        (str(i), {"subject": entity(f"P{i}"), "predicate": "affects", "object": entity(f"Q{i}")})
        for i in range(120)
    ]
    tables = run_rows(tmp_path, source, 7)
    entity_keys = [row["entity_key"] for row in tables["entity"]]
    relation_keys = [row["relation_key"] for row in tables["relation"]]
    assert len(entity_keys) == 240 and len(relation_keys) == 120
    assert len({key[:2] for key in entity_keys}) > 10
    assert entity_keys == sorted(entity_keys)
    assert relation_keys == sorted(relation_keys)
    # Ids are row positions, continuing across buckets.
    assert [row["entity_id"] for row in tables["entity"]] == list(range(240))
    assert [row["relation_id"] for row in tables["relation"]] == list(range(120))
    # Each relation under both endpoints' entity keys and reference keys.
    endpoints = [(row["key"], row["relation_id"]) for row in tables["relation_endpoint"]]
    assert endpoints == sorted(endpoints) and len(endpoints) == 480
    assert sum(row["key_kind"] == "entity" for row in tables["relation_endpoint"]) == 240
