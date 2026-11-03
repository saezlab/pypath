from fastapi.testclient import TestClient
from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_core.fixtures import write_resource


def write_source(root, name, rows, relations=()):
    write_resource(root / "resources" / name / "1", rows, relations)


def chemical(
    key,
    identifier,
    *,
    ns="inchikey",
    ref=None,
    label="Chemical",
    taxon="9606",
    kind="chemical_entity",
):
    """A chemical; an InChIKey-identified one is its own reference unless ``ref`` says."""
    if ref is None and ns == "inchikey":
        ref = "inchikey:" + identifier
    return dict(
        entity_key=key,
        entity_type=kind,
        namespace=ns,
        identifier=identifier,
        label=label,
        taxon=taxon,
        reference_entity_key=ref or None,
        identifiers=[{"ns": "name", "id": label}],
    )


def test_connectivity_groups_scope_pagination_and_missing_keys(tmp_path):
    a = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    b = "AAAAAAAAAAAAAA-CCCCCCCCCC-D"
    other = "AAAAAAAAAABBBB-BBBBBBBBBB-C"  # same first ten, different connectivity
    rows = [
        chemical("a", a),
        chemical("b", b),
        chemical("c", "CID:1", ns="pubchem", ref="inchikey:" + a),
        chemical("d", other),
        chemical("e", "missing", ns="pubchem"),
        chemical("f", "bad-key", ref=""),
        chemical("g", "conflict", ns="pubchem", ref="inchikey:" + other, label="Conflict chemical"),
        chemical("protein", a, kind="protein"),
    ]
    write_source(tmp_path, "one", rows)
    write_source(tmp_path, "two", [rows[0]])
    write_source(
        tmp_path,
        "three",
        [chemical("g", "conflict", ns="pubchem", ref="inchikey:" + a, label="Conflict chemical")],
    )
    engine = ParquetQueries(data_root=tmp_path)
    result = engine.search_connectivity_groups(limit=1, member_limit=1)
    group = result["groups"][0]
    assert group["connectivity"] == "AAAAAAAAAAAAAA"
    assert group["member_count"] == 3  # cross-resource copies count once
    assert [m["entityPk"] for m in group["members"]] == ["a"]
    next_members = engine.search_connectivity_groups(
        group_key=group["group_key"], member_cursor=group["nextMemberCursor"], member_limit=1
    )
    assert next_members["groups"][0]["members"][0]["entityPk"] == "b"
    rest = engine.search_connectivity_groups(cursor=result["nextCursor"])
    assert len(rest["groups"]) == 5
    assert rest["groups"][0]["connectivity"] == "AAAAAAAAAABBBB"
    assert all(g["member_count"] == 1 for g in rest["groups"])
    assert rest["nextCursor"] is None
    assert any(
        g["members"][0]["entityPk"] == "protein" and g["connectivity"] is None
        for g in rest["groups"]
    )
    assert engine.search_connectivity_groups(filters={"ncbi_tax_id": ["10090"]})["groups"] == []
    scoped = engine.search_connectivity_groups(filters={"sources": ["two"]})
    assert scoped["groups"][0]["member_count"] == 1
    assert (
        engine.search_connectivity_groups(query="conflict")["groups"][0]["members"][0]["entityPk"]
        == "g"
    )
    client = TestClient(create_app(queries=engine))
    response = client.post("/entities/groups", json={"filters": {"entity_pks": ["b"]}})
    assert response.status_code == 200
    assert response.json()["groups"][0]["member_count"] == 1
    assert response.json()["groups"][0]["members"][0]["entityPk"] == "b"
    assert client.post("/entities/connectivity-groups", json={"limit": 0}).status_code == 422

    assert client.post("/entities/groups", json={"strategy": "unknown"}).status_code == 422


def test_grouping_uses_the_same_text_match_scope_as_entity_search(tmp_path):
    a = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    write_source(
        tmp_path,
        "source",
        [
            chemical("a", a, label="Needle"),
            chemical("b", "AAAAAAAAAAAAAA-CCCCCCCCCC-D", label="Other needle"),
        ],
    )
    engine = ParquetQueries(data_root=tmp_path)
    result = engine.search_entity_groups(query="Needle")
    assert result["groups"][0]["member_count"] == 1
    assert result["groups"][0]["members"][0]["entityPk"] == "a"


def test_group_cards_merge_all_members_and_sort_before_pagination(tmp_path):
    rows = [dict(chemical("small", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", label="Small"), annotations=[])]
    for i in range(7):
        row = chemical(f"large-{i}", f"ZZZZZZZZZZZZZZ-AAAAAAAAA{chr(65 + i)}-A", label=f"Name {i}")
        row["annotations"] = [dict(term="description", value=f"Description {i}", source="one")]
        rows.append(row)
    write_source(tmp_path, "one", rows)
    extra = dict(
        rows[-1], annotations=[dict(term="description", value="Other source", source="two")]
    )
    write_source(tmp_path, "two", [extra])
    engine = ParquetQueries(data_root=tmp_path)
    first = engine.search_entity_groups(limit=1, member_limit=1, include_details=True)
    group = first["groups"][0]
    assert group["member_count"] == 7
    assert len(group["members"]) == 1
    entity = group["entity"]
    assert entity["entityPk"] == "connectivity:ZZZZZZZZZZZZZZ"
    assert entity["groupMemberKeys"] == ["large-0"]
    assert entity["groupMemberCount"] == 7
    assert {f"Name {i}" for i in range(7)} <= {v["identifier"] for v in entity["identifiers"]}
    assert len(entity["entityAttributes"]) == 8
    assert entity["sources"] == ["one", "two"]
    second = engine.search_entity_groups(limit=1, cursor=first["nextCursor"])
    assert second["groups"][0]["member_count"] == 1
    assert second["nextCursor"] is None
    scoped = engine.search_entity_groups(filters={"sources": ["one"]}, include_details=True)[
        "groups"
    ][0]["entity"]
    assert len(scoped["entityAttributes"]) == 7
    assert scoped["groupResources"] == ["one"]


def test_group_cards_skip_hydration_until_opened(tmp_path, monkeypatch):
    rows = [
        dict(
            chemical(str(i), "AAAAAAAAAAAAAA-BBBBBBBBBB-C", label=f"Member {i}"),
            annotations=[dict(term="description", value=f"Details {i}", source="one")],
        )
        for i in range(8)
    ]
    write_source(tmp_path, "one", rows)
    engine = ParquetQueries(data_root=tmp_path)
    queries = []
    original = engine._fetch_dicts

    def record(sql, *args, **kwargs):
        queries.append(sql)
        return original(sql, *args, **kwargs)

    monkeypatch.setattr(engine, "_fetch_dicts", record)
    client = TestClient(create_app(queries=engine))
    summary = client.post("/entities/groups", json={"filters": {"sources": ["one"]}}).json()[
        "groups"
    ][0]["entity"]
    # Grouping, then member rows; no identifiers or annotations.
    assert not any("entity_identifier" in q or "entity_annotation" in q for q in queries)
    assert summary["entityAttributes"] is None
    assert summary["groupDetailsLoaded"] is False
    assert len(summary["groupMemberKeys"]) == 5
    assert summary["groupMemberCount"] == 8
    remaining = client.post(
        "/entities/groups",
        json={
            "group_key": summary["entityPk"],
            "filters": summary["groupFilters"],
            "member_cursor": summary["groupMemberCursor"],
        },
    ).json()["groups"][0]["entity"]
    assert len(remaining["groupMemberKeys"]) == 3
    assert remaining["groupMemberCursor"] is None
    detail = client.post(
        "/entities/groups",
        json={
            "group_key": summary["entityPk"],
            "query": summary["groupQuery"],
            "filters": summary["groupFilters"],
            "resources": summary["groupResources"],
            "include_details": True,
        },
    ).json()["groups"][0]["entity"]
    assert detail["entityPk"] == summary["entityPk"]
    assert detail["groupDetailsLoaded"] is True
    assert len(detail["entityAttributes"]) == 8
    assert detail["groupMemberKeys"] == summary["groupMemberKeys"]


def test_group_cache_reuses_results_and_invalidates_on_source_change(tmp_path, monkeypatch):
    write_source(tmp_path, "one", [chemical("a", "AAAAAAAAAAAAAA-BBBBBBBBBB-C")])
    engine = ParquetQueries(data_root=tmp_path)
    first = engine.search_entity_groups()
    with monkeypatch.context() as patch:
        patch.setattr(
            engine,
            "_fetch_dicts",
            lambda *a, **kw: (_ for _ in ()).throw(AssertionError("repeated scan")),
        )
        cached = engine.search_entity_groups()
        assert cached["groups"] == first["groups"]
    write_source(
        tmp_path,
        "one",
        [
            chemical("a", "AAAAAAAAAAAAAA-BBBBBBBBBB-C"),
            chemical("b", "AAAAAAAAAAAAAA-CCCCCCCCCC-D"),
        ],
    )
    assert engine.search_entity_groups()["groups"][0]["member_count"] == 2


def test_group_detail_collections_are_paged(tmp_path):
    row = chemical("a", "AAAAAAAAAAAAAA-BBBBBBBBBB-C")
    row["annotations"] = [
        dict(term="description", value=f"Note {i}", source="one") for i in range(45)
    ]
    row["identifiers"] += [dict(ns="name", id=f"Name {i}") for i in range(45)]
    write_source(tmp_path, "one", [row])
    client = TestClient(create_app(queries=ParquetQueries(data_root=tmp_path)))
    values = []
    offset = 0
    while True:
        entity = client.post(
            "/entities/groups",
            json={
                "group_key": "connectivity:AAAAAAAAAAAAAA",
                "include_details": True,
                "detail_offset": offset,
            },
        ).json()["groups"][0]["entity"]
        assert len(entity["entityAttributes"]) <= 20
        assert len(entity["identifiers"]) <= 20
        values += [a["value"] for a in entity["entityAttributes"]]
        if entity["detailNextCursor"] is None:
            break
        offset = int(entity["detailNextCursor"])
    assert values == [f"Note {i}" for i in range(45)]


def test_group_relationships_resolve_members_beyond_preview(tmp_path):
    write_source(
        tmp_path,
        "one",
        [chemical(str(i), "AAAAAAAAAAAAAA-BBBBBBBBBB-C") for i in range(7)],
        [
            dict(
                relation_key="r",
                subject_entity_key="6",
                object_entity_key="0",
                predicate="has_part",
                category="membership",
                sources=["one"],
                evidence_count=1,
            )
        ],
    )
    client = TestClient(create_app(queries=ParquetQueries(data_root=tmp_path)))
    response = client.post(
        "/entities/group-relationships", json={"group_key": "connectivity:AAAAAAAAAAAAAA"}
    )
    assert response.status_code == 200
    assert response.json()["relationshipsTotal"] == 1
    assert response.json()["relationships"][0]["groupOutgoing"] is True


def test_group_name_is_the_preferred_member_label_on_every_page(tmp_path):
    first = chemical("a", "QNAYBMKLOCPYGJ-UQEXSWPGSA-N", label="alanine-d7")
    second = chemical("b", "QNAYBMKLOCPYGJ-AAAAAAAAAA-N", label="L-Alanine")
    second["identifiers"] += [
        dict(ns="name", id=f"Long descriptive chemical name {i}") for i in range(40)
    ]
    write_source(tmp_path, "one", [first, second])
    engine = ParquetQueries(data_root=tmp_path)
    group = engine.search_entity_groups(member_limit=1)["groups"][0]
    assert group["entity"]["displayName"] == "L-Alanine"
    # The members' canonical identifiers and labels.
    assert len(group["entity"]["identifiers"]) == 4
    for offset in (0, 20, 40):
        entity = engine.search_entity_groups(
            group_key=group["group_key"], include_details=True, detail_offset=offset
        )["groups"][0]["entity"]
        assert entity["displayName"] == "L-Alanine"
        assert len(entity["identifiers"]) <= 20
        assert entity["identifiersNextCursor"] == (
            str(offset + 20) if offset + 20 < entity["identifiersTotal"] else None
        )
    assert (
        engine.search_entity_groups(filters={"entity_pks": ["a"]})["groups"][0]["entity"][
            "displayName"
        ]
        == "alanine-d7"
    )
