"""Bounded annotation preparation retains the exact old four-column link set.

The oracle below is the complete query from commit
5703acf50630e8c1716b0b183c0510087200b66f, before owner-aware preparation. It is
intentionally independent of the replacement helpers and batching strategy.
Synthetic Parquets exercise nullable attribution, scopes and occurrence owners;
no source build, resolver, network or PostgreSQL connection is used.
"""

from collections import Counter
from copy import deepcopy
from types import SimpleNamespace

import duckdb
import pytest

from omnipath_postgres import loader as aligned_loader, projection as aligned_projection
from published_fixture import annotation, fixture_rows, write_fixture


OLD_LINK_SQL = """WITH true_statement_annotation AS (
    SELECT a.* FROM ap_annotation_occurrence a
    WHERE a.owner_kind='relation' AND a.term IS NOT NULL AND NOT EXISTS (
        SELECT 1 FROM ap_annotation_occurrence observed
        WHERE observed.owner_kind='evidence' AND observed.resource=a.resource
          AND observed.version=a.version AND observed.owner_key=a.owner_key
          AND observed.annotation_key=a.annotation_key
          AND observed.source IS NOT DISTINCT FROM a.source
          AND observed.dataset IS NOT DISTINCT FROM a.dataset
          AND (CASE WHEN observed.scope='evidence' THEN 'relation'
                    ELSE coalesce(observed.scope,'relation') END)
            = (CASE WHEN a.scope='evidence' THEN 'relation'
                    ELSE coalesce(a.scope,'relation') END)
    )
), own_annotation AS (
    SELECT * FROM ap_annotation_occurrence WHERE owner_kind='evidence' AND term IS NOT NULL
    UNION ALL SELECT * FROM true_statement_annotation
)
SELECT DISTINCT e.source_id,e.relation_evidence_id,a.annotation_key,sc.id::SMALLINT
FROM own_annotation a JOIN ap_evidence e
ON a.resource=e.resource AND a.version=e.version AND a.owner_key=e.relation_key
AND (a.owner_kind='relation' OR (a.owner_kind='evidence' AND a.evidence_ordinal=e.ordinal))
AND (a.owner_kind='evidence' OR e.synthetic OR
     ((a.source IS NULL OR a.source=e.source) AND (a.dataset IS NULL OR a.dataset=e.dataset)))
JOIN ap_vocab_annotation_scope sc ON sc.name=CASE WHEN a.scope='evidence' THEN 'relation'
                                               ELSE coalesce(a.scope,'relation') END
WHERE a.term IS NOT NULL"""


def _attribute(value, *, source=None, dataset=None, scope="relation", term="description"):
    item = annotation(term, value=value, scope=scope)
    item.update(source=source, dataset=dataset)
    return item


def _statement(key, evidence, annotations):
    statement = deepcopy(fixture_rows()[1][0])
    statement.update(
        relation_key=key,
        evidence=deepcopy(evidence),
        evidence_count=len(evidence or []),
        annotations=deepcopy(annotations),
    )
    return statement


def _evidence(source, dataset, row, attributes=()):
    return {
        "source": source,
        "dataset": dataset,
        "row_id": str(row),
        "upstream_id": "source-record:" + str(row),
        "annotations": deepcopy(list(attributes)),
    }


def _fixtures(root, *, first_version="v1"):
    entities = deepcopy(fixture_rows()[0][:2])
    for item in entities:
        item["annotations"] = []
    union_copy = _attribute("union-copy", source="src-a", dataset="D", scope="evidence")
    null_copy = _attribute("null-vs-empty", source=None, dataset=None)
    owns = [
        union_copy,
        deepcopy(union_copy),
        _attribute("own-conflict", source="other-source", dataset="other-dataset", scope="object"),
        _attribute("multi", scope="object"),
        _attribute("multi", scope="subject"),
        _attribute("empty-own", scope=""),
        _attribute("custom-own", scope="custom_scope"),
        _attribute("untyped", scope="object", term=None),
    ]
    evidence = [
        _evidence("src-a", "D", 10, owns),
        _evidence("src-a", "D", 11),
        _evidence(None, None, 12, [null_copy]),
        _evidence("", "", 13),
    ]
    statement_attributes = [
        dict(union_copy, scope="relation"),
        deepcopy(union_copy),
        _attribute("multi", scope="object"),
        _attribute("multi", scope="subject"),
        dict(null_copy, source="", dataset=""),
        _attribute("wildcard", scope=None),
        _attribute("only-D", dataset="D", scope="subject"),
        _attribute("coalesced-blocked", source="res-a", dataset=aligned_projection.CLAIM_DATASET),
        _attribute("wrong-dataset", source="src-a", dataset="missing-dataset"),
        _attribute("empty-generic", scope=""),
        _attribute("custom-generic", scope="custom_scope"),
        _attribute("untyped", scope="subject", term=None),
    ]
    first = [
        _statement("shared-key", evidence, statement_attributes),
        _statement("qualified-other-2", [_evidence("src-a", "D", 14)], [union_copy]),
        _statement(
            "synthetic",
            None,
            [_attribute("synthetic-only", source="not-owner", dataset="unknown", scope="object")],
        ),
        _statement("empty-owner", [], []),
    ]
    second = [_statement("shared-key", [_evidence("src-a", "D", 20)], [union_copy])]
    selected = []
    for source, version, statements in (("res-a", first_version, first), ("res-b", "v2", second)):
        directory = root / source
        write_fixture(directory, entities=entities, relations=statements, payloads=[])
        selected.append(SimpleNamespace(directory=directory, source=source, version=version))
    assert len(entities) * 2 + len(first) + len(second) == 9
    return tuple(selected)


def _links(connection, plan):
    query = next(item for item in plan.queries if item.table == "relation_evidence_annotation")
    assert query.columns == (
        "source_id",
        "relation_evidence_id",
        "annotation_key",
        "annotation_scope_id",
    )
    return Counter(connection.execute(query.query).fetchall())


def _decoded_links(connection, links):
    ownership = {
        (source, evidence): (resource, version, key, ordinal)
        for resource, version, key, ordinal, source, evidence in connection.execute(
            "SELECT resource,version,relation_key,ordinal,source_id,relation_evidence_id FROM ap_evidence"
        ).fetchall()
    }
    values = dict(connection.execute("SELECT annotation_key,value FROM ap_annotation").fetchall())
    scopes = dict(connection.execute("SELECT id,name FROM ap_vocab_annotation_scope").fetchall())
    decoded = {}
    for (source, evidence, key, scope), multiplicity in links.items():
        decoded.setdefault(ownership[source, evidence], Counter())[
            (values[key], scopes[scope])
        ] += multiplicity
    return decoded


def test_bounded_preparation_matches_frozen_old_query_and_source_owned_semantics(tmp_path):
    selected = _fixtures(tmp_path)
    observed = []

    def observer(event, **fields):
        observed.append((event, fields))

    with duckdb.connect() as connection:
        # Use the loader's actual forwarding boundary. An annotation helper
        # payload named "event" would collide with _emit's positional event.
        plan = aligned_projection.prepare_aligned_release(
            connection,
            selected,
            progress=lambda fields: aligned_loader._emit(observer, "projection_progress", **fields),
        )
        actual = _links(connection, plan)
        expected = Counter(connection.execute(OLD_LINK_SQL).fetchall())
        assert actual == expected
        assert actual and set(actual.values()) == {1}
        common = {
            ("wildcard", "relation"),
            ("empty-generic", ""),
            ("custom-generic", "custom_scope"),
        }
        decoded = _decoded_links(connection, actual)
        assert decoded == {
            ("res-a", "v1", "shared-key", 0): Counter(
                common
                | {
                    ("union-copy", "relation"),
                    ("own-conflict", "object"),
                    ("multi", "subject"),
                    ("multi", "object"),
                    ("empty-own", ""),
                    ("custom-own", "custom_scope"),
                    ("only-D", "subject"),
                }
            ),
            ("res-a", "v1", "shared-key", 1): Counter(common | {("only-D", "subject")}),
            ("res-a", "v1", "shared-key", 2): Counter(common | {("null-vs-empty", "relation")}),
            ("res-a", "v1", "shared-key", 3): Counter(common | {("null-vs-empty", "relation")}),
            ("res-a", "v1", "qualified-other-2", 0): Counter({("union-copy", "relation")}),
            ("res-a", "v1", "synthetic", -1): Counter({("synthetic-only", "object")}),
            ("res-b", "v2", "shared-key", 0): Counter({("union-copy", "relation")}),
        }
        assert all(event == "projection_progress" for event, _fields in observed)
        link_progress = [
            fields for _event, fields in observed if fields["phase"] == "prepare_annotation_links"
        ]
        assert link_progress and all("event" not in fields for fields in link_progress)
        assert [(fields["stage"], fields["state"]) for fields in link_progress[:4]] == [
            ("annotation_input", "start"),
            ("annotation_input", "done"),
            ("evidence_input", "start"),
            ("evidence_input", "done"),
        ]
        assert link_progress[-1]["stage"] == "shard"
        assert link_progress[-1]["state"] == "done"
        assert link_progress[-1]["shard"] == link_progress[-1]["shard_count"] - 1


def test_resource_versions_keep_distinct_evidence_owners_for_same_graph_keys(tmp_path):
    retained = []
    for version in ("version-before", "version-after"):
        selected = _fixtures(tmp_path / version, first_version=version)
        with duckdb.connect() as connection:
            plan = aligned_projection.prepare_aligned_release(connection, selected)
            actual = _links(connection, plan)
            assert actual == Counter(connection.execute(OLD_LINK_SQL).fetchall())
            owned_evidence = set(
                connection.execute(
                    "SELECT source_id,relation_evidence_id FROM ap_evidence WHERE resource='res-a'"
                ).fetchall()
            )
            owned = {
                (source, evidence)
                for source, evidence, _annotation, _scope in actual
                if (source, evidence) in owned_evidence
            }
            canonical = connection.execute("SELECT relation_id FROM ap_relation").fetchall()
            retained.append((owned, canonical))
    assert retained[0][0].isdisjoint(retained[1][0])
    assert retained[0][1] == retained[1][1]


@pytest.mark.parametrize("shard_count", [1, 2, 7, 128])
def test_configured_shards_keep_owner_complete_and_insert_each_link_exactly_once(
    tmp_path, shard_count
):
    selected = _fixtures(tmp_path)
    with duckdb.connect() as connection:
        aligned_projection.prepare_aligned_release(connection, selected)
        expected = Counter(connection.execute(OLD_LINK_SQL).fetchall())
        events, completed, membership = [], [], {}
        owner_membership = {}
        previous = Counter()

        def progress(event):
            nonlocal previous
            events.append(event.copy())
            assert event["shard_count"] == shard_count
            if event["stage"] != "shard":
                return
            shard = event["shard"]
            assert isinstance(shard, int) and 0 <= shard < shard_count
            if event["state"] == "start":
                if not membership:
                    annotation_owners = connection.execute("""SELECT DISTINCT
                        resource,version,owner_key,owner_shard FROM ap_annotation_link_input""").fetchall()
                    evidence_rows = connection.execute("""SELECT
                        resource,version,owner_key,owner_shard,source_id,relation_evidence_id
                        FROM ap_evidence_link_input""").fetchall()
                    for resource, version, key, bucket in [
                        *annotation_owners,
                        *(row[:4] for row in evidence_rows),
                    ]:
                        assert isinstance(bucket, int) and 0 <= bucket < shard_count
                        owner = (resource, version, key)
                        assert owner_membership.setdefault(owner, bucket) == bucket
                    for resource, version, key, bucket, source, evidence in evidence_rows:
                        assert owner_membership[resource, version, key] == bucket
                        membership[source, evidence] = bucket
                assert shard == len(completed)
                return
            assert event["state"] == "done"
            current = Counter(
                connection.execute(
                    "SELECT source_id,relation_evidence_id,annotation_key,annotation_scope_id "
                    "FROM ap_copy_relation_evidence_annotation"
                ).fetchall()
            )
            assert all(current[row] >= count for row, count in previous.items())
            added = current - previous
            # This checks actual appended four-column rows, not an expected
            # hash expression copied from the implementation. No row may be
            # inserted by a bucket other than its evidence owner's bucket.
            assert all(
                membership[source, evidence] == shard
                for source, evidence, _annotation, _scope in added
            )
            assert not set(previous).intersection(added)
            assert all(multiplicity == 1 for multiplicity in current.values())
            previous = current
            completed.append(shard)

        aligned_projection._prepare_relation_evidence_annotation_copy(
            connection,
            shard_count=shard_count,
            on_progress=progress,
        )
        assert previous == expected and previous
        if shard_count > 1:
            assert len(set(owner_membership.values())) > 1
        if shard_count == 2:
            # With five owners, this also exercises unrelated owners sharing
            # a shard, instead of relying solely on one-owner buckets.
            assert max(Counter(owner_membership.values()).values()) > 1
        assert completed == list(range(shard_count))
        assert [(event["stage"], event["state"]) for event in events[:4]] == [
            ("annotation_input", "start"),
            ("annotation_input", "done"),
            ("evidence_input", "start"),
            ("evidence_input", "done"),
        ]
        assert len(events) == 4 + 2 * shard_count
        assert all(event["stage"] == "shard" for event in events[4:])


def test_observed_null_term_retains_old_anti_match_without_creating_an_evidence_link(tmp_path):
    selected = _fixtures(tmp_path)
    with duckdb.connect() as connection:
        aligned_projection.prepare_aligned_release(connection, selected)
        # This contrived staged collision tests the old query boundary: observed
        # evidence rows participate in anti-matching before output term filters.
        # No published Parquet is edited and no NULL term becomes an assertion.
        connection.execute("""INSERT INTO ap_annotation_occurrence
            SELECT * REPLACE ('evidence' AS owner_kind,0::BIGINT AS evidence_ordinal,
                              99::BIGINT AS ordinal,NULL::VARCHAR AS term)
            FROM ap_annotation_occurrence WHERE resource='res-a'
              AND owner_key='qualified-other-2' AND owner_kind='relation' LIMIT 1""")
        expected = Counter(connection.execute(OLD_LINK_SQL).fetchall())
        aligned_projection._prepare_relation_evidence_annotation_copy(connection, shard_count=7)
        actual = Counter(
            connection.execute("SELECT * FROM ap_copy_relation_evidence_annotation").fetchall()
        )
        assert actual == expected
        assert ("res-a", "v1", "qualified-other-2", 0) not in _decoded_links(connection, actual)


@pytest.mark.parametrize("shard_count", [0, -1, True, None, 1.5, "2"])
def test_invalid_shard_bounds_fail_before_any_sql(shard_count):
    class NoSql:
        def execute(self, *_args, **_kwargs):
            pytest.fail("Invalid shard configuration must fail before preparing SQL inputs")

    with pytest.raises(ValueError, match="shard_count"):
        aligned_projection._prepare_relation_evidence_annotation_copy(
            NoSql(), shard_count=shard_count
        )
