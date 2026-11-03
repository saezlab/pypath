"""Hub and decision kv writers: key order, value codec, sharding, parity with the Parquet rows."""

from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, str(Path(__file__).parents[3] / "resolver" / "tests"))
from identity_snapshot import downgrade_to_v2, write_hub_index  # noqa: E402

from omnipath_build.identity import build_hub_kv_dir, build_identity_kv  # noqa: E402
from omnipath_resolver import identity_kv  # noqa: E402
from omnipath_resolver.identity_kv import HubKv, KvMissing  # noqa: E402

TRICKY = [
    "a", "A", "a b", "a-b", "a:b", "a_b", "ab", "é", "ü:1", "日本", "Z9", "z9", "\x01x", "~", "CHEBI:1",
    "CHEBI:10", "CHEBI:2", "x" * 300,
]  # fmt: skip


def test_value_codec_prefix_and_compression():
    small, large = ["a", "b"], [["chebi", "x" * 40]] * 50
    raw, packed = identity_kv.pack(small), identity_kv.pack(large)
    assert raw[:1] == identity_kv.RAW and packed[:1] == identity_kv.ZSTD
    assert identity_kv.unpack(raw) == small and identity_kv.unpack(packed) == large
    with pytest.raises(ValueError):
        identity_kv.unpack(b"\x07abc")


def test_duckdb_varchar_order_is_lmdb_byte_order():
    c = duckdb.connect()
    c.execute("CREATE TABLE t(ns VARCHAR, identifier VARCHAR)")
    c.executemany(
        "INSERT INTO t VALUES (?, ?)", [(ns, i) for ns in ("a", "a-b", "ab", "é") for i in TRICKY]
    )
    rows = c.execute("SELECT ns, identifier FROM t ORDER BY ns, identifier").fetchall()
    keys = [ns.encode() + b"\0" + i.encode() for ns, i in rows]
    assert keys == sorted(keys) and len(set(keys)) == len(keys)
    values = [
        r[0]
        for r in c.execute("SELECT x FROM (SELECT unnest(?) x) ORDER BY x", [TRICKY]).fetchall()
    ]
    assert [v.encode() for v in values] == sorted(v.encode() for v in TRICKY)


def hub(directory: Path, records):
    write_hub_index(directory, "chebi", records)
    return directory


def tricky_records():
    records = {}
    for n, local in enumerate(TRICKY):
        rows = [
            ("chebi", local),
            ("name", f"name {n}"),
            ("synonym", "x" * 400 if n % 3 == 0 else "s"),
        ]
        rows += [("cas", f"{n}-0-0"), ("hmdb", TRICKY[(n + 1) % len(TRICKY)])]
        if n == 2:
            rows.append(("kegg", "y" * 600))  # too long for an LMDB key: counted, not stored
        records[local] = (None if n % 2 else "9606", bool(n % 2), rows)
    records["no-names"] = (None, False, [])  # a record without any by_record row
    return records


@pytest.mark.parametrize("shards,workers", [(1, 1), (16, 2)])
def test_hub_kv_holds_every_parquet_row(tmp_path, shards, workers):
    directory = hub(tmp_path / "chebi" / "0123456789ab", tricky_records())
    manifest = build_hub_kv_dir(
        directory,
        id_shards=shards,
        rec_shards=shards,
        memory="256MB",
        threads=1,
        workers=workers,
        min_free_gib=0.01,
    )
    assert manifest["totals"]["rec"]["keys"] == len(tricky_records())
    assert manifest["totals"]["id"]["skipped_long_keys"] > 0  # the 300-character identifier
    assert build_hub_kv_dir(directory, min_free_gib=0.01) == manifest  # idempotent
    assert not (directory / "kv.building").exists()
    kv = HubKv("chebi", directory)
    try:
        by_id = pq.read_table(
            list(directory.glob("by_id/part=*/data.parquet"))[0].parent.parent
        ).to_pylist()
        expected = {}
        for row in by_id:
            expected.setdefault((row["ns"], row["identifier"]), []).append(
                (row["local_id"], row["tag"])
            )
        pairs = [p for p in expected if identity_kv.id_key(*p) is not None]
        found = kv.ids(pairs + [("chebi", "nope")])
        assert set(found) == set(pairs)
        assert all(sorted(found[p]) == sorted(expected[p]) for p in pairs)
        records = {
            r["local_id"]: r for r in pq.read_table(directory / "records.parquet").to_pylist()
        }
        by_record = {}
        for path in directory.glob("by_record/part=*/data.parquet"):
            for r in pq.read_table(path).to_pylist():
                by_record.setdefault(r["local_id"], []).append((r["source_type"], r["value"]))
        usable = [i for i in records if identity_kv.rec_key("chebi", i)]
        heads, rows = kv.heads(usable), kv.rows(usable)
        for local in usable:
            r = records[local]
            assert heads[local] == (r["taxon"], r["anchor"], r["anchor_count"], r["reviewed"])
            assert sorted(rows[local]) == sorted(by_record.get(local, []))
        assert kv.heads(["missing"]) == {} and kv.ids([("chebi", "x" * 600)]) == {}
        if shards == 16:  # a key reads only its own shard
            kv.close()
            kv = HubKv("chebi", directory)
            kv.ids([("chebi", "CHEBI:1")])
            assert {s for _, s in kv.opened()} == {identity_kv.shard_of("CHEBI:1", 16)}
    finally:
        kv.close()


SYMBOLS = ["TP53", "TP53-AS1", "TP5", "TP530", "Tp53", "é", "A", "A1", "A-1", "a b", "x" * 495]
TAXA = ["9606", "10090", None, "0", " 9606 ", "None", ""]  # the last four are no taxon or 9606


def symbol_records():
    """Every symbol in each namespace under records of every kind of taxon."""
    records = {}
    for s, symbol in enumerate(SYMBOLS):
        for t, taxon in enumerate(TAXA):
            rows = [("genesymbol", symbol), ("genesymbol-syn", symbol), ("entrez", f"{s}-{t}")]
            if t % 2 and s:
                rows.append(("genesymbol-syn", SYMBOLS[0]))
            records[f"{s}-{t}"] = (taxon, False, rows)
    return records


def symbol_hub(tmp_path, shards=1):
    directory = tmp_path / "entrez" / "0123456789ab"
    write_hub_index(directory, "entrez", symbol_records())
    build_hub_kv_dir(
        directory, id_shards=shards, rec_shards=shards, memory="256MB", threads=1, min_free_gib=0.01
    )
    return directory


def clean(taxon):
    """The taxon of a record as the runtime reads it."""
    text = (taxon or "").strip()
    return "" if text in ("", "0", "None") else text


@pytest.mark.parametrize("shards", [1, 16])
def test_symbol_rows_are_keyed_by_the_taxon_of_their_record(tmp_path, shards):
    directory = symbol_hub(tmp_path, shards)
    manifest = identity_kv.json.loads((directory / "kv" / "manifest.json").read_text())
    assert manifest["format"] == identity_kv.HUB_KV_FORMAT
    assert manifest["taxon_keyed_namespaces"] == ["genesymbol", "genesymbol-syn"]
    by_id = pq.read_table(directory / "by_id").to_pylist()
    expected = {}
    for row in by_id:
        if row["ns"] in identity_kv.SYMBOL_NAMESPACES:
            key = (row["ns"], row["identifier"], clean(row["taxon"]))
            expected.setdefault(key, set()).add((row["local_id"], row["tag"]))
    kv = HubKv("entrez", directory)
    try:
        assert kv.by_taxon
        taxa = {"9606", "10090", "7227"}
        found = kv.taxon_ids((ns, i, t) for ns, i, _ in expected for t in taxa)
        for (ns, identifier, _), _rows in expected.items():
            for taxon in taxa:
                if identity_kv.id_key(ns, identifier, taxon) is None:
                    continue  # too long for an LMDB key: not stored
                # the rows of that taxon and those of the records without one, nothing else
                want = expected.get((ns, identifier, taxon), set()) | expected.get(
                    (ns, identifier, ""), set()
                )
                got = found.get((ns, identifier, taxon), [])
                assert sorted(got) == sorted(want) and len(got) == len(set(got)), (ns, identifier)
        # the odd taxa ("0", "None", blanks) are no taxon, and " 9606 " is 9606
        tp53 = lambda taxon: {local for local, _ in found[("genesymbol", "TP53", taxon)]}  # noqa: E731
        assert {"0-2", "0-3", "0-5", "0-6"} <= tp53("7227") and "0-0" not in tp53("7227")
        assert "0-4" in tp53("9606") and "0-4" not in tp53("10090")
        # other namespaces keep their keys; a symbol is no longer found without a taxon
        assert kv.ids([("entrez", "0-0")]) == {("entrez", "0-0"): [("0-0", "native")]}
        with pytest.raises(ValueError, match="keyed by taxon"):
            kv.ids([("genesymbol", "TP53")])
        with pytest.raises(ValueError, match="not keyed by taxon"):
            kv.taxon_ids([("entrez", "0-0", "9606")])
        long = "x" * 495  # fits a key of format v2, not one with the taxon: skipped, counted
        assert identity_kv.id_key("genesymbol", long) is not None
        assert identity_kv.id_key("genesymbol", long, "10090") is None
        assert manifest["totals"]["id"]["skipped_long_keys"] > 0
    finally:
        kv.close()


def test_a_symbol_lookup_reads_two_keys_per_taxon(tmp_path):
    directory = symbol_hub(tmp_path)
    kv = HubKv("entrez", directory)
    try:
        env = kv._env("id", 0)
        gets = []
        original = env.get_many
        env.get_many = lambda name, keys: gets.append(sorted(keys)) or original(name, keys)
        kv.taxon_ids([("genesymbol", "TP53", "9606")])
        kv.taxon_ids([("genesymbol", "TP53", "9606"), ("genesymbol", "TP53", "10090")])
        assert [len(g) for g in gets] == [2, 3]  # the key without a taxon is shared
        assert gets[0][0].endswith(b"\0TP53\0") or gets[0][1].endswith(b"\0TP53\0")
    finally:
        kv.close()


def test_format_v2_is_read_and_replaced_only_when_the_new_store_is_complete(tmp_path):
    directory = symbol_hub(tmp_path)
    downgrade_to_v2(directory)
    kv = HubKv("entrez", directory)
    try:
        assert not kv.by_taxon and kv.manifest["format"] == "omnipath-hub-kv-v2"
        symbol = ("genesymbol", "TP53")
        assert len(kv.ids([symbol])[symbol]) > len(SYMBOLS)  # every taxon
        with pytest.raises(ValueError, match="no taxon keys"):
            kv.taxon_ids([("genesymbol", "TP53", "9606")])
    finally:
        kv.close()
    # a partial build of the older format is discarded, and the readable v2 store stays until
    # the new one is complete
    building = directory / "kv.building"
    (building / "format").parent.mkdir()
    (building / "format").write_text("omnipath-hub-kv-v2")
    (building / "id.0").mkdir()
    (building / "id.0" / "stats.json").write_text("{}")
    build_hub_kv_dir(directory, memory="256MB", threads=1, min_free_gib=0.01)
    assert not building.exists()
    kv = HubKv("entrez", directory)
    try:
        assert kv.by_taxon and kv.taxon_ids([("genesymbol", "TP53", "10090")])
    finally:
        kv.close()


def test_hub_kv_refuses_a_changed_index(tmp_path):
    directory = hub(tmp_path / "chebi" / "0123456789ab", tricky_records())
    build_hub_kv_dir(directory, memory="256MB", threads=1, min_free_gib=0.01)
    HubKv("chebi", directory).close()
    (directory / "manifest.json").write_text('{"hub": "chebi", "changed": true}')
    with pytest.raises(KvMissing, match="different index"):
        HubKv("chebi", directory)
    with pytest.raises(KvMissing, match="build-hub-kv"):
        HubKv("chebi", tmp_path)


def test_writer_rejects_keys_out_of_order(tmp_path):
    from omnipath_build.identity.hubkv import Writer

    import pyarrow as pa

    writer = Writer(tmp_path / "w", "id", 0)
    writer.add(pa.array(["b"]), pa.array(["x"]))
    with pytest.raises(AssertionError, match="not ascending"):
        writer.add(pa.array(["a"]), pa.array(["x"]))
    with pytest.raises(AssertionError, match="order"):  # out of order inside one batch
        writer.add(pa.array(["d", "c"]), pa.array(["x", "y"]))
    writer.env.close()


def test_values_round_trip_separators():
    from omnipath_resolver.identity_kv import decode_head, decode_ids, decode_rows

    assert decode_ids("a\x1fnative\x1eb\x1fclaim".encode()) == [("a", "native"), ("b", "claim")]
    assert decode_head("\x1f\x1f0\x1f0".encode()) == (None, None, 0, False)
    value = "9606\x1funiprot:P1\x1f1\x1f1\x1ename\x1fp53\x1egenesymbol\x1fTP53".encode()
    assert decode_head(value) == ("9606", "uniprot:P1", 1, True)
    assert decode_rows(value) == [("name", "p53"), ("genesymbol", "TP53")]
    assert decode_rows("\x1f\x1f0\x1f0".encode()) == []


def test_identity_kv_is_complete_and_optional_candidates(tmp_path):
    from identity_snapshot import build_snapshot

    identity = build_snapshot(tmp_path)
    manifest = build_identity_kv(identity, min_free_gib=0.001)
    assert build_identity_kv(identity, min_free_gib=0.001) == manifest  # idempotent
    assert set(manifest["counts"]) == set(identity_kv.DECISION_DBS)
    assert manifest["counts"]["exc"] == pq.read_table(identity / "exceptions.parquet").num_rows
    assert manifest["counts"]["cand"] > 0
    (identity / "record_candidates.parquet").unlink()  # an older snapshot: empty db
    assert build_identity_kv(identity, force=True, min_free_gib=0.001)["counts"]["cand"] == 0
    kv = identity_kv.DecisionsKv(identity, manifest["fingerprint"])
    try:
        assert kv.get("exc", ["chebi:CHEBI:99999", "nope"]) == {
            "chebi:CHEBI:99999": ["chebi:CHEBI:99999", "quarantined", True]
        }
        assert dict(kv.items("lipid"))
        with pytest.raises(KvMissing, match="another fingerprint"):
            identity_kv.DecisionsKv(identity, "other")
    finally:
        kv.close()
