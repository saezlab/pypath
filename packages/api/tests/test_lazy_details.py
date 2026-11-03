from fastapi.testclient import TestClient
from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_api.queries.store.inventory import ReleaseStore
from omnipath_core.fixtures import nested_rows, rewrite_resource, write_resource


def resource(root, version, label):
    entities = [
        dict(
            entity_key=k,
            entity_type="chemical_entity",
            namespace="test",
            identifier=k,
            label=label if k == "a" else "Other",
        )
        for k in ["a", "b"]
    ]
    relations = [
        dict(
            relation_key=f"r{i}",
            subject_entity_key="a",
            object_entity_key="b",
            predicate="has_part",
            category="membership",
            sources=["test"],
            evidence_count=1,
            annotations=[dict(term="note", value=f"{label}-{i}", source="test")],
        )
        for i in range(3)
    ]
    write_resource(root / "resources/test" / version, entities, relations)


def test_lazy_core_does_not_query_relations_and_reuses_cache(tmp_path, monkeypatch):
    resource(tmp_path, "1", "First")
    engine = ParquetQueries(tmp_path)
    client = TestClient(create_app(queries=engine))
    table = engine._table

    def entity_tables_only(name, *args):
        assert not name.startswith("relation"), "relation scan"
        return table(name, *args)

    monkeypatch.setattr(engine, "_table", entity_tables_only)
    first = client.get("/entities/a?includeRelationships=false")
    assert first.status_code == 200 and first.json()["entity"]["label"] == "First"
    monkeypatch.setattr(
        engine,
        "_fetch_entities_by_keys",
        lambda *a, **kw: (_ for _ in ()).throw(AssertionError("uncached entity")),
    )
    assert client.get("/entities/a?includeRelationships=false").json() == first.json()


def test_relationship_pagination_cache_and_release_isolation(tmp_path, monkeypatch):
    resource(tmp_path, "1", "First")
    ReleaseStore(tmp_path).publish(dict(schema_version=1, version="1", resources={"test": "1"}))
    engine = ParquetQueries(tmp_path)
    client = TestClient(create_app(queries=engine))
    first = client.get("/entities/a/relationships?limit=2").json()
    assert first["relationshipsTotal"] == 3 and first["nextCursor"] == "2"
    second = client.get("/entities/a/relationships?limit=2&offset=2").json()
    assert second["nextCursor"] is None
    assert [
        r["relation"]["relationPk"] for r in first["relationships"] + second["relationships"]
    ] == ["r0", "r1", "r2"]
    with monkeypatch.context() as m:
        m.setattr(
            engine,
            "_fetch_dicts",
            lambda *a, **kw: (_ for _ in ()).throw(AssertionError("uncached relationships")),
        )
        assert client.get("/entities/a/relationships?limit=2").json() == first
    resource(tmp_path, "2", "Second")
    assert (
        client.get("/entities/a?includeRelationships=false").json()["entity"]["label"] == "Second"
    )
    latest = client.get("/entities/a/relationships?limit=2").json()
    assert latest["relationships"][0]["annotations"][0]["value"] == "Second-0"
    pinned = client.get(
        "/entities/a/relationships?limit=2", headers={"X-OmniPath-Release": "1"}
    ).json()
    assert pinned == first
    assert client.get("/entities/missing/relationships").status_code == 404
    assert client.get("/entities/a/relationships?limit=0").status_code == 422


def test_group_relationships_include_both_endpoints_once_and_paginate(tmp_path):
    resource(tmp_path, "1", "First")
    client = TestClient(create_app(queries=ParquetQueries(tmp_path)))
    body = {"member_keys": ["a", "b", "a"], "limit": 2}
    response = client.post("/entities/group-relationships", json=body)
    assert response.status_code == 200
    first = response.json()
    assert first["relationshipsTotal"] == 3
    assert first["nextCursor"] == "2"
    second = client.post("/entities/group-relationships", json={**body, "offset": 2}).json()
    assert second["nextCursor"] is None
    assert [
        r["relation"]["relationPk"] for r in first["relationships"] + second["relationships"]
    ] == ["r0", "r1", "r2"]
    assert client.post("/entities/group-relationships", json={"member_keys": []}).status_code == 422
    scoped = client.post("/entities/group-relationships", json={**body, "resources": ["missing"]})
    assert scoped.status_code == 400


def test_normal_entity_details_default_to_bounded_pages(tmp_path):
    resource(tmp_path, "1", "First")
    path = tmp_path / "resources/test/1/entity.parquet"
    rows = nested_rows(path)
    rows[0]["annotations"] = [
        dict(term="description", value=f"Note {i}", source="test") for i in range(45)
    ]
    rewrite_resource(path, entities=rows)
    client = TestClient(create_app(queries=ParquetQueries(tmp_path)))
    first = client.get("/entities/a?includeRelationships=false").json()["entity"]
    assert len(first["entityAttributes"]) == 20
    assert first["entityAttributesTotal"] == 45
    assert first["detailNextCursor"] == "20"
    second = client.get("/entities/a?includeRelationships=false&detail_offset=20").json()["entity"]
    assert len(second["entityAttributes"]) == 20
    third = client.get("/entities/a?includeRelationships=false&detail_offset=40").json()["entity"]
    assert len(third["entityAttributes"]) == 5
    assert third["detailNextCursor"] is None


def test_relationships_and_counts_span_resources(tmp_path):
    resource(tmp_path, "1", "First")
    write_resource(
        tmp_path / "resources/other/1",
        nested_rows(tmp_path / "resources/test/1/entity.parquet"),
        [
            # a self-loop has both endpoints in the entity and counts once
            dict(relation_key=f"s{i}", subject_entity_key="a", object_entity_key="a",
                 predicate="has_member", category="association", sources=["other"],
                 evidence_count=1)
            for i in range(2)
        ],
    )  # fmt: skip
    engine = ParquetQueries(tmp_path)

    def read():
        engine._detail_cache.clear()
        engine._relationship_cache.clear()
        pages = [engine.get_entity_relationships("a", limit=2, offset=o) for o in (0, 2, 4)]
        return engine.get_entity_details("a"), pages

    expected = read()
    details, pages = expected
    assert details["summary"]["interactionCount"] == 5
    assert {a["relationPk"] for a in details["annotations"]} == {"s0", "s1"}
    assert pages[0]["relationshipsTotal"] == 5
    assert [r["relation"]["relationPk"] for p in pages for r in p["relationships"]] == [
        "s0", "s1", "r0", "r1", "r2",
    ]  # fmt: skip


def test_descriptions_lead_the_first_attribute_page():
    from omnipath_api.queries.entity_details import page_details

    entity = dict(
        entityAttributes=[dict(term="in_taxon", value=str(i)) for i in range(30)]
        + [dict(term="up:Function_Annotation", value="Regulates the pump.")]
    )
    first = page_details(entity)["entityAttributes"]
    assert first[0]["term"] == "up:Function_Annotation"
    assert len(first) == 20


def test_one_collection_pages_on_its_own():
    from omnipath_api.queries.entity_details import page_details

    entity = dict(
        entityPk="a",
        identifiers=[dict(id=str(i), type="x") for i in range(5)],
        entityAttributes=[dict(term="in_taxon", value=str(i)) for i in range(30)],
        sources=["s"],
    )
    page = page_details(entity, limit=20, offset=20, field="entityAttributes")
    assert [a["value"] for a in page["entityAttributes"]] == [str(i) for i in range(20, 30)]
    assert page["entityAttributesTotal"] == 30 and page["entityAttributesNextCursor"] is None
    # other collections are left out, so a client keeps what it has
    assert "identifiers" not in page and "sources" not in page and "detailNextCursor" not in page


def test_attributes_beyond_descriptions_are_grouped_by_term():
    from omnipath_api.queries.entity_details import page_details

    terms = ["in_taxon", "has_role", "in_taxon", "description", "has_role"]
    entity = dict(entityAttributes=[dict(term=t, value=str(i)) for i, t in enumerate(terms)])
    assert [(a["term"], a["value"]) for a in page_details(entity)["entityAttributes"]] == [
        ("description", "3"), ("has_role", "1"), ("has_role", "4"),
        ("in_taxon", "0"), ("in_taxon", "2"),
    ]  # fmt: skip
