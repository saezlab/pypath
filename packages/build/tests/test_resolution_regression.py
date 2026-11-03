"""Resolution regression harness: fingerprints, persistence, resolve and diff."""

import json
import sys
import types

import pyarrow.parquet as pq
import pytest

from omnipath_build.regression import diff as regression_diff
from omnipath_build.regression import extract as regression_extract
from omnipath_build.regression import resolve as regression_resolve
from omnipath_build.regression import store
from omnipath_resolver import RawEntityObservation
from omnipath_resolver.observations import key

INCHIKEY = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


def vote(ns, identifier, **extra):
    row = dict(
        ns=ns,
        identifier=identifier,
        scope="",
        anchor="",
        target=1,
        route=1,
        ordinal=0,
        primary=False,
    )
    row.update(extra)
    row["lookup_key"] = key(row["target"], row["route"], ns, row["scope"], identifier)
    return row


def query(target=1, entity_type="small_molecule", namespace="chebi", identifier="CHEBI:1"):
    return dict(
        target=target,
        entity_type=entity_type,
        namespace=namespace,
        identifier=identifier,
        taxon="",
    )


def test_fingerprint_ignores_order_and_input_id_but_not_gene_only():
    a = [vote("chebi", "CHEBI:1", ordinal=0), vote("hmdb", "HMDB0000001", ordinal=1)]
    b = [dict(a[1], ordinal=0), dict(a[0], ordinal=1, input_id="whatever")]
    assert store.observation_fingerprint(1, a) == store.observation_fingerprint(1, b)
    assert store.observation_fingerprint(1, a) != store.observation_fingerprint(2, a)
    flagged = [dict(a[0], gene_only=True), a[1]]
    assert store.observation_fingerprint(1, a) != store.observation_fingerprint(1, flagged)
    assert len(store.observation_fingerprint(1, a)) == 32


def write_observations(root, resource, observations, occurrences=None):
    """observations: list of (query, votes); returns the fingerprints in input order."""
    writer = store.ObservationWriter(root / resource, flush_rows=2)
    ids = []
    for q, rows in observations:
        fingerprint = store.observation_fingerprint(q["target"], rows)
        writer.add(fingerprint, q, rows)
        ids.append(fingerprint)
    for dataset, picks in (occurrences or {}).items():
        writer.add_occurrences(dataset, {ids[i]: n for i, n in picks.items()})
    writer.finalize()
    return ids


def test_persistence_round_trip_streams_matching_votes(tmp_path):
    def votes(i):
        rows = [vote("chebi", f"CHEBI:{i}", ordinal=0)]
        rows += [vote("hmdb", f"HMDB{i:04d}{j}", ordinal=j + 1) for j in range(i % 3)]
        return rows

    # One observation without any vote (they all collapse into one), then 1-3 votes each.
    observations = [(query(identifier="none"), [])] + [
        (query(identifier=f"CHEBI:{i}"), votes(i)) for i in range(1, 10)
    ]
    ids = write_observations(tmp_path, "r", observations)
    directory = tmp_path / "r"
    for name in ("queries", "votes", "occurrences"):
        assert (directory / f"{name}.parquet").is_file()
    assert not (directory / ".parts").exists()
    queries = pq.read_table(directory / "queries.parquet").to_pylist()
    assert [q["input_id"] for q in queries] == sorted(ids)
    assert {q["library"] for q in queries} == {"chemical"}

    seen_queries, seen_votes = [], []
    for batch_queries, batch_votes in store.iter_resolve_batches(
        directory, batch_size=3, vote_chunk=4
    ):
        assert {v["input_id"] for v in batch_votes} <= {q["input_id"] for q in batch_queries}
        seen_queries += [q["input_id"] for q in batch_queries]
        seen_votes += batch_votes
    assert seen_queries == sorted(ids)
    expected = sum(len(rows) for _, rows in observations)
    assert len(seen_votes) == expected
    sample = seen_votes[0]
    assert set(sample) == {f.name for f in store.VOTE_SCHEMA}
    assert isinstance(sample["lookup_key"], bytes) and sample["gene_only"] is False


def test_gene_only_and_binary_key_survive_the_round_trip(tmp_path):
    rows = [vote("hgnc", "HGNC:5", target=2, gene_only=True)]
    write_observations(tmp_path, "g", [(query(target=2, entity_type="transcript"), rows)])
    ((_, votes),) = list(store.iter_resolve_batches(tmp_path / "g", 10))
    assert votes[0]["gene_only"] is True
    assert votes[0]["lookup_key"] == rows[0]["lookup_key"]
    assert votes[0]["target"] == 2


def test_occurrences_aggregate_and_empty_resource(tmp_path):
    obs = [(query(identifier="CHEBI:1"), [vote("chebi", "CHEBI:1")])]
    ids = write_observations(tmp_path, "r", obs, occurrences={"a": {0: 2}, "b": {0: 1}})
    rows = pq.read_table(tmp_path / "r" / "occurrences.parquet").to_pylist()
    assert {(r["dataset"], r["n"]) for r in rows} == {("a", 2), ("b", 1)}
    assert {r["input_id"] for r in rows} == set(ids)
    write_observations(tmp_path, "empty", [])
    assert list(store.iter_resolve_batches(tmp_path / "empty")) == []
    assert store.list_resources(tmp_path) == ["empty", "r"]


def test_collector_dedupes_across_datasets(tmp_path):
    writer = store.ObservationWriter(tmp_path / "r")
    collector = regression_extract._Collector(writer)
    rows = [vote("chebi", "CHEBI:1")]
    fingerprint = store.observation_fingerprint(1, rows)
    item = (fingerprint, query(), rows)
    for dataset in ("one", "two"):
        collector.start(dataset)
        collector.absorb(dict(rows=3, skipped=1, observations=[item, item]))
        collector.flush_occurrences()
    assert collector.dataset["one"]["new_observations"] == 1
    assert collector.dataset["two"]["new_observations"] == 0
    assert collector.dataset["two"]["observations"] == 2
    counts = writer.finalize()
    assert counts == {"queries": 1, "votes": 1, "occurrences": 2}


def test_entity_observations_emit_all_vote_fields_including_gene_only():
    entities = {
        "k1": RawEntityObservation(
            "k1",
            "small_molecule",
            "chebi",
            "CHEBI:15377",
            identifiers=[{"ns": "inchikey", "id": INCHIKEY}],
        ),
        "k2": RawEntityObservation(
            "k2",
            "small_molecule",
            "chebi",
            "CHEBI:15377",
            identifiers=[{"ns": "inchikey", "id": INCHIKEY}],
        ),
        "k3": RawEntityObservation("k3", "transcript", "ensembl", "ENST00000269305", taxon="9606"),
        "k4": RawEntityObservation("k4", "protein", "uniprot", "P04637", taxon="9606"),
        "k5": RawEntityObservation("k5", "macromolecular_complex", "complexportal", "CPX-1"),
    }
    observations, skipped = regression_extract.entity_observations(entities)
    assert skipped == 1  # the complex has no reference library
    by_key = {o[1]["entity_key"]: o for o in observations}
    assert by_key["k1"][0] == by_key["k2"][0]  # identical observations share a fingerprint
    chemical = by_key["k1"][2]
    assert {r["ns"] for r in chemical} == {"chebi", "inchikey"}
    assert [r["anchor"] for r in chemical if r["ns"] == "inchikey"] == [INCHIKEY]
    assert not any(r["gene_only"] for r in chemical)
    transcript = by_key["k3"]
    assert transcript[1]["target"] == 2 and transcript[1]["taxon"] == "9606"
    assert transcript[2] and all(r["gene_only"] for r in transcript[2])
    protein = by_key["k4"][2]
    assert protein and not any(r["gene_only"] for r in protein)


def test_process_records_runs_the_silver_extractor():
    def mapper(row):
        if row["id"] is None:
            return None
        return {
            "type": "protein",
            "identifiers": [{"type": "uniprot", "value": row["id"]}],
        }

    rows = [(0, {"id": "P04637"}), (1, {"id": None}), (2, {"id": "P04637"})]
    observations, skipped = regression_extract.process_records("src", "ds", rows, mapper)
    assert skipped == 0 and len(observations) == 1
    fingerprint, q, votes = observations[0]
    assert [(v["ns"], v["identifier"]) for v in votes] == [("uniprot", "P04637")]
    assert q["entity_type"] == "protein"


class FakeRuntime:
    """Resolve by looking at the first vote's identifier; mimics the identity runtime's shapes."""

    def __init__(self, table):
        self.table = table
        self.calls = []

    def resolve(self, queries, votes):
        self.calls.append((len(queries), len(votes)))
        by_input = {}
        for v in votes:
            by_input.setdefault(v["input_id"], []).append(v)
        results, records = [], {}
        for q in queries:
            first = sorted(by_input.get(q["input_id"], []), key=lambda v: v["ordinal"])
            found = self.table.get(first[0]["identifier"]) if first else None
            entities = sorted(found or [])
            for e in entities:
                records[e] = dict(label=e.upper())
            row = dict(
                input_id=q["input_id"],
                outcome="Resolved" if entities else "NotFound",
                entities=entities,
                candidate_count=len(entities),
            )
            if q["target"] == 2:
                row.update(
                    gene_mapping_status="mapped" if entities else "missing",
                    gene_candidates=[],
                    protein_entity_id=None,
                )
            results.append(row)
        metrics = dict(
            lookup_seconds=0.1,
            decision_seconds=0.2,
            entity_fetch_seconds=0.3,
            total_seconds=0.6,
            unique_keys=len(votes),
        )
        return dict(results=results, records=records), metrics


def synthetic_observations(root):
    obs = [(query(identifier=f"X{i}"), [vote("chebi", f"X{i}")]) for i in range(6)] + [
        (
            query(target=2, entity_type="protein", namespace="uniprot", identifier="P1"),
            [vote("uniprot", "P1", target=2)],
        ),
    ]
    ids = write_observations(root, "res", obs, occurrences={"d1": {0: 4, 1: 1}, "d2": {0: 1}})
    return {f"X{i}": ids[i] for i in range(6)} | {"P1": ids[6]}


def test_resolve_resource_writes_results_and_metrics(tmp_path):
    observations = tmp_path / "obs"
    ids = synthetic_observations(observations)
    runtime = FakeRuntime(
        {"X0": ["inchikey:A"], "X1": ["inchikey:B", "inchikey:C"], "P1": ["uniprot:P1"]}
    )
    metrics = regression_resolve.resolve_all(
        "fake-runtime", observations, tmp_path / "res", batch_size=3, runtime=runtime
    )
    assert [c[0] for c in runtime.calls] == [3, 3, 1]
    table = pq.read_table(tmp_path / "res" / "res" / "results.parquet")
    assert table.schema.names == store.RESULT_SCHEMA.names
    rows = {r["input_id"]: r for r in table.to_pylist()}
    assert len(rows) == 7
    assert rows[ids["X1"]]["entities"] == ["inchikey:B", "inchikey:C"]
    assert rows[ids["X1"]]["entity_labels"] == ["INCHIKEY:B", "INCHIKEY:C"]
    assert rows[ids["X2"]]["outcome"] == "NotFound"
    assert rows[ids["P1"]]["gene_mapping_status"] == "mapped"
    saved = json.loads((tmp_path / "res" / "metrics.json").read_text())
    resource = saved["resources"]["res"]
    assert resource["observations"] == 7 and len(resource["batches"]) == 3
    assert resource["batches"][0]["lookup_seconds"] == pytest.approx(0.1)
    assert resource["decision_seconds"] == pytest.approx(0.6)
    assert resource["record_seconds"] == pytest.approx(0.9)
    assert resource["outcomes"]["chemical"] == {"NotFound": 4, "Resolved": 2}
    assert saved["totals"]["observations"] == 7
    assert metrics["totals"]["observations_per_second"] > 0
    # A second run keeps finished resources unless forced.
    again = FakeRuntime({})
    regression_resolve.resolve_all("fake", observations, tmp_path / "res", runtime=again)
    assert again.calls == []


def test_load_runtime_opens_only_identity_directories(tmp_path, monkeypatch):
    library = tmp_path / "library"
    library.mkdir()
    (library / "manifest.json").write_text(json.dumps({"format": "omnipath-reference-v1"}))
    with pytest.raises(ValueError, match="not an identity directory"):
        regression_resolve.load_runtime(library)

    snapshot = tmp_path / "snapshot"
    snapshot.mkdir()
    (snapshot / "manifest.json").write_text(json.dumps({"format": "omnipath-identity-v2"}))
    created = []

    class Identity:
        def __init__(self, path, cache_dir=None):
            created.append((path, cache_dir))

    module = types.ModuleType("omnipath_resolver.identity_runtime")
    module.IdentityRuntime = Identity
    monkeypatch.setitem(sys.modules, "omnipath_resolver.identity_runtime", module)
    runtime, _ = regression_resolve.load_runtime(snapshot, cache_dir="/cache")
    assert isinstance(runtime, Identity) and created == [(snapshot, "/cache")]


def test_diff_classifies_every_change_and_keeps_votes(tmp_path):
    observations = tmp_path / "obs"
    ids = synthetic_observations(observations)
    a = FakeRuntime(
        {
            "X0": ["inchikey:A"],  # same in both
            "X1": ["inchikey:B"],  # b adds an entity
            "X2": ["inchikey:C", "inchikey:D"],  # b drops an entity
            "X3": ["inchikey:E"],  # resolved -> none
            "X4": ["inchikey:F"],  # different entities
            "P1": ["uniprot:P1"],
        }
    )
    b = FakeRuntime(
        {
            "X0": ["inchikey:A"],
            "X1": ["inchikey:B", "inchikey:B2"],
            "X2": ["inchikey:C"],
            "X5": ["inchikey:G"],  # none -> resolved
            "X4": ["goslin:F"],
            "P1": ["uniprot:P1"],
        }
    )
    regression_resolve.resolve_all("a", observations, tmp_path / "a", runtime=a)
    regression_resolve.resolve_all("b", observations, tmp_path / "b", runtime=b)
    summary = regression_diff.compare(
        observations, tmp_path / "a", tmp_path / "b", tmp_path / "out"
    )

    differences = pq.read_table(tmp_path / "out" / "differences.parquet").to_pylist()
    kinds = {r["input_id"]: r["change_type"] for r in differences}
    assert kinds == {
        ids["X1"]: "b_adds_entities",
        ids["X2"]: "b_drops_entities",
        ids["X3"]: "resolved_to_none",
        ids["X4"]: "different_entities",
        ids["X5"]: "none_to_resolved",
    }
    row = next(r for r in differences if r["input_id"] == ids["X1"])
    assert row["entities_a"] == ["inchikey:B"] and row["entities_b"] == [
        "inchikey:B",
        "inchikey:B2",
    ]
    assert row["votes"] == ["chebi:X1"]
    assert row["occurrences"] == 1 and row["datasets"] == ["d1"]
    assert row["resource"] == "res" and row["library"] == "chemical"

    counts = {(c["library"], c["change_type"]): c["n"] for c in summary["change_types"]}
    assert counts[("chemical", "same")] == 1  # X0
    assert counts[("gene_protein", "same")] == 1  # P1
    assert counts[("chemical", "different_entities")] == 1
    crosstab = {(c["library"], c["outcome_a"], c["outcome_b"]): c["n"] for c in summary["crosstab"]}
    assert crosstab[("chemical", "Resolved", "NotFound")] == 1
    assert crosstab[("chemical", "NotFound", "Resolved")] == 1
    assert crosstab[("chemical", "Resolved", "Resolved")] == 4
    assert ("chemical", "NotFound", "NotFound") not in crosstab

    markdown = (tmp_path / "out" / "summary.md").read_text()
    assert "b_adds_entities" in markdown and "Outcome cross-tab" in markdown
    assert (tmp_path / "out" / "crosstab.parquet").is_file()
    assert (tmp_path / "out" / "change_types.parquet").is_file()


def test_diff_flags_missing_results_and_other_field_changes(tmp_path):
    observations = tmp_path / "obs"
    ids = synthetic_observations(observations)
    a = FakeRuntime({"P1": ["uniprot:P1"]})
    b = FakeRuntime({})
    regression_resolve.resolve_all("a", observations, tmp_path / "a", runtime=a)
    regression_resolve.resolve_all("b", observations, tmp_path / "b", runtime=b)
    # Same (empty) entities but a different gene status only matters when entities match:
    # drop one result from B to make it missing.
    path = tmp_path / "b" / "res" / "results.parquet"
    table = pq.read_table(path)
    keep = [i for i, v in enumerate(table.column("input_id").to_pylist()) if v != ids["X0"]]
    pq.write_table(table.take(keep), path)
    summary = regression_diff.compare(
        observations, tmp_path / "a", tmp_path / "b", tmp_path / "out"
    )
    counts = {(c["library"], c["change_type"]): c["n"] for c in summary["change_types"]}
    assert counts[("chemical", "missing_result")] == 1
    assert counts[("chemical", "same")] == 5
    assert counts[("gene_protein", "resolved_to_none")] == 1
    kinds = {
        r["input_id"]: r["change_type"]
        for r in pq.read_table(tmp_path / "out" / "differences.parquet").to_pylist()
    }
    assert kinds[ids["X0"]] == "missing_result"
    assert kinds[ids["P1"]] == "resolved_to_none"


def test_sample_is_deterministic_stratified_and_shared_by_resolve_and_diff(tmp_path):
    observations = tmp_path / "obs"
    obs = [(query(identifier=f"X{i}"), [vote("chebi", f"X{i}")]) for i in range(40)] + [
        (
            query(target=2, entity_type="protein", namespace="uniprot", identifier=f"P{i}"),
            [vote("uniprot", f"P{i}", target=2)],
        )
        for i in range(25)
    ]
    write_observations(observations, "res", obs)
    first = [
        q["input_id"]
        for qs, _ in store.iter_resolve_batches(observations / "res", 4, sample=10)
        for q in qs
    ]
    again = [
        q["input_id"]
        for qs, _ in store.iter_resolve_batches(observations / "res", 7, sample=10)
        for q in qs
    ]
    assert sorted(first) == sorted(again) and len(first) == 20  # 10 per library
    libraries = {
        q["input_id"]: q["library"]
        for q in pq.read_table(observations / "res" / "queries.parquet").to_pylist()
    }
    assert sorted(libraries[i] for i in first).count("chemical") == 10
    # a sample larger than a library keeps all of it
    everything = [
        q for qs, _ in store.iter_resolve_batches(observations / "res", 100, sample=30) for q in qs
    ]
    assert len(everything) == 30 + 25
    # votes belong to exactly the sampled queries
    for qs, votes in store.iter_resolve_batches(observations / "res", 4, sample=10):
        assert {v["input_id"] for v in votes} == {q["input_id"] for q in qs}

    table = {f"X{i}": [f"inchikey:{i}"] for i in range(40)}
    table.update({f"P{i}": [f"uniprot:P{i}"] for i in range(25)})
    changed = dict(table, **{f"X{i}": ["other"] for i in range(0, 40, 2)})
    for name, mapping in (("a", table), ("b", changed)):
        regression_resolve.resolve_all(
            name, observations, tmp_path / name, runtime=FakeRuntime(mapping), sample=10
        )
    assert pq.read_table(tmp_path / "a" / "res" / "results.parquet").num_rows == 20
    summary = regression_diff.compare(
        observations, tmp_path / "a", tmp_path / "b", tmp_path / "out", sample=10
    )
    counts = {(c["library"], c["change_type"]): c["n"] for c in summary["change_types"]}
    assert sum(counts.values()) == 20 and ("chemical", "missing_result") not in counts
    assert "Sample: 10" in (tmp_path / "out" / "summary.md").read_text()
