"""Automatic grouping over mixed resource tables."""

import json

import pytest
from fastapi.testclient import TestClient

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_core.keys import entity_key
from omnipath_core.fixtures import write_resource

CONNECTIVITY = "AAAAAAAAAAAAAA"
INCHI_A = CONNECTIVITY + "-BBBBBBBBBB-C"
INCHI_B = CONNECTIVITY + "-CCCCCCCCCC-D"
INCHI_OTHER = "ZZZZZZZZZZZZZZ-BBBBBBBBBB-C"


def record(kind, ns, identifier, *, ref=None, label=None, aliases=(), taxon="9606"):
    return dict(
        entity_key=entity_key(kind, ns, identifier),
        entity_type=kind,
        namespace=ns,
        identifier=identifier,
        reference_entity_key=ref,
        gene_reference_keys=[ref] if ref and ref.startswith("entrez:") else [],
        label=label or identifier,
        taxon=taxon,
        has_hierarchy=False,
        parent_count=0,
        child_count=0,
        identifiers=[dict(ns=ns_, id=value) for ns_, value in aliases],
        annotations=[dict(term="description", value=f"Detail {identifier}", source="a")],
        evidence=[],
    )


@pytest.fixture
def mixed_engine(tmp_path):
    rows = {
        "gene": record("gene", "entrez", "1", ref="entrez:1", label="Needle gene"),
        "protein": record("protein", "uniprot", "P00001", ref="entrez:1", label="Gene product"),
        "rna": record("rna_product", "entrez", "1", ref="entrez:1", label="Gene RNA"),
        "target": record("gene", "entrez", "2", ref="entrez:2", label="Target"),
        "target_protein": record("protein", "uniprot", "P00002", ref="entrez:2"),
        "conflicting_gene": record("protein", "uniprot", "P00003", ref="entrez:1"),
        "missing_gene": record("protein", "uniprot", "P00004", ref=None),
        "native_gene": record("protein", "uniprot", "P00005", ref="uniprot:P00005", label="Native"),
        # Chemicals group by the connectivity of their InChIKey reference.
        "chemical_a": record(
            "small_molecule",
            "inchikey",
            INCHI_A,
            ref="inchikey:" + INCHI_A,
            label="Needle chemical",
        ),
        "chemical_b": record(
            "small_molecule", "inchikey", INCHI_B, ref="inchikey:" + INCHI_B, label="Stereo B"
        ),
        "chemical_c": record(
            "small_molecule",
            "pubchem",
            "1",
            ref="inchikey:" + INCHI_A,
            label="Alias chemical",
            aliases=[("name", "Alanine")],
        ),
        # Even an accidental gene reference must never attach a chemical to a gene group.
        "chemical_gene_ref": record("small_molecule", "chebi", "1", ref="entrez:1"),
        "conflicting_chemical": record("small_molecule", "pubchem", "2", ref="inchikey:" + INCHI_A),
        "missing_chemical": record("small_molecule", "pubchem", "4"),
        "protein_with_inchikey": record("protein", "inchikey", INCHI_A, ref="inchikey:" + INCHI_A),
    }
    copies = [
        dict(rows["conflicting_gene"], reference_entity_key="entrez:2"),
        dict(rows["missing_gene"], reference_entity_key="entrez:1"),
        dict(rows["conflicting_chemical"], reference_entity_key="inchikey:" + INCHI_OTHER),
        dict(
            rows["chemical_c"],
            annotations=[
                dict(term="description", value="Other source chemical detail", source="b")
            ],
        ),
    ]
    relations = [
        ("chemical", "chemical_c", "target"),
        ("protein", "protein", "chemical_a"),
        ("outside", "native_gene", "target"),
    ]
    for source, entities in [("a", list(rows.values())), ("b", copies)]:
        relation_rows = [
            dict(
                relation_key=name,
                subject_entity_key=rows[subject]["entity_key"],
                object_entity_key=rows[object_]["entity_key"],
                subject_type=rows[subject]["entity_type"],
                object_type=rows[object_]["entity_type"],
                subject_label=rows[subject]["label"],
                object_label=rows[object_]["label"],
                predicate="has_part",
                category="membership",
                taxon="9606",
                sources=[source],
            )
            for name, subject, object_ in (relations if source == "a" else [])
        ]
        write_resource(tmp_path / "resources" / source / "1", entities, relation_rows)
    engine = ParquetQueries(tmp_path)
    return engine, rows


def test_auto_combines_groups_and_preserves_all_singletons(mixed_engine):
    engine, records = mixed_engine
    groups = engine.search_entity_groups(strategy="auto", include_member_keys=True)["groups"]
    # Most studied first (members' relations), then by key: the structure group and
    # gene 2 have two relations, gene 1 one.
    assert [(g["group_key"], g["member_count"]) for g in groups[:2]] == [
        ("connectivity:" + CONNECTIVITY, 3),
        ("gene:entrez:2", 2),
    ]
    by_key = {g["group_key"]: g for g in groups}
    chemical, gene = by_key["connectivity:" + CONNECTIVITY], by_key["gene:entrez:1"]
    assert gene["member_count"] == 3
    assert chemical["entity"]["groupStrategy"] == "chemical_connectivity"
    assert chemical["entity"]["canonicalIdentifierType"] == "connectivity"
    assert chemical["entity"]["entityType"] == "small_molecule"
    # The most connected well-named member: only "Alias chemical" has a relation.
    assert chemical["entity"]["displayName"] == "Alias chemical"
    assert chemical["entity"]["entityAttributes"] is None
    assert gene["entity"]["groupStrategy"] == "gene_reference"
    assert gene["entity"]["memberEntityTypes"] == ["gene", "protein", "rna_product"]
    # named by its gene, not by a product
    assert gene["entity"]["displayName"] == "Needle gene"
    # a group card counts its members' relations
    assert chemical["entity"]["relationCount"] == 2 and gene["entity"]["relationCount"] == 1
    expected_singletons = {
        r["entity_key"]
        for name, r in records.items()
        if name
        not in {
            "gene",
            "protein",
            "rna",
            "target",
            "target_protein",
            "chemical_a",
            "chemical_b",
            "chemical_c",
        }
    }
    assert {g["entity"]["entityPk"] for g in groups if not g["is_group"]} == expected_singletons
    member_keys = [
        key for g in groups for key in g["entity"].get("groupMemberKeys", [g["entity"]["entityPk"]])
    ]
    assert len(member_keys) == len(set(member_keys)) == len(records)
    assert set(member_keys) == {r["entity_key"] for r in records.values()}
    # API compatibility: omitting strategy still groups only chemicals.
    default = engine.search_entity_groups()["groups"]
    assert not any(g["group_key"].startswith("gene:") for g in default)
    assert any(g["group_key"] == "connectivity:" + CONNECTIVITY for g in default)
    gene_only = engine.search_entity_groups(strategy="gene_reference")["groups"]
    assert not any(g["group_key"].startswith("connectivity:") for g in gene_only)
    # Concrete gene followups exclude even malformed chemicals with a gene ref.
    concrete_gene = next(g for g in gene_only if g["group_key"] == "gene:entrez:1")
    assert concrete_gene["member_count"] == 3
    malformed_chemical = next(
        g
        for g in gene_only
        if g["entity"]["entityPk"] == records["chemical_gene_ref"]["entity_key"]
    )
    assert not malformed_chemical["is_group"]


def test_auto_group_and_member_cursors_cover_each_key_once(mixed_engine):
    engine, _ = mixed_engine
    expected = engine.search_entity_groups(strategy="auto", limit=100)["groups"]
    collected, cursor = [], None
    while True:
        result = engine.search_entity_groups(strategy="auto", limit=1, cursor=cursor)
        collected.extend(result["groups"])
        cursor = result["nextCursor"]
        if not cursor:
            break
    assert [g["group_key"] for g in collected] == [g["group_key"] for g in expected]
    for group in [g for g in expected if g["is_group"]][:3]:
        members, cursor = [], None
        while True:
            result = engine.search_entity_groups(
                strategy="auto",
                group_key=group["group_key"],
                member_limit=1,
                member_cursor=cursor,
            )["groups"][0]
            members.extend(result["entity"]["groupMemberKeys"])
            cursor = result["nextMemberCursor"]
            if not cursor:
                break
        assert members == sorted({m["entityPk"] for m in group["members"]})
    with pytest.raises(ValueError, match="Invalid group cursor"):
        engine.search_entity_groups(strategy="auto", cursor=json.dumps([-1, 1, "bad"]))


def test_auto_scope_text_resources_and_filters(mixed_engine):
    engine, records = mixed_engine
    result = engine.search_entity_groups(strategy="auto", query="Needle")["groups"]
    assert {g["group_key"] for g in result} == {"gene:entrez:1", "connectivity:" + CONNECTIVITY}
    assert all(g["member_count"] == 1 for g in result)
    assert {m["entityPk"] for g in result for m in g["members"]} == {
        e["entityPk"] for e in engine.search_entities_api(query="Needle")["entities"]
    }
    alias = engine.search_entity_groups(strategy="auto", query="Alanine")["groups"]
    assert alias[0]["member_count"] == 1
    assert alias[0]["members"][0]["entityPk"] == records["chemical_c"]["entity_key"]
    proteins = engine.search_entity_groups(strategy="auto", filters={"entity_types": ["protein"]})[
        "groups"
    ]
    assert all(not g["connectivity"] for g in proteins)
    assert next(g for g in proteins if g["group_key"] == "gene:entrez:1")["member_count"] == 1
    assert (
        engine.search_entity_groups(strategy="auto", filters={"ncbi_tax_id": ["10090"]})["groups"]
        == []
    )
    scoped = engine.search_entity_groups(strategy="auto", resources=["a"])["groups"]
    assert next(g for g in scoped if g["group_key"] == "gene:entrez:1")["member_count"] == 4
    chemistry = next(g for g in scoped if g["connectivity"] == CONNECTIVITY)
    assert chemistry["member_count"] == 4
    exact = engine.search_entity_groups(
        strategy="auto", filters={"entity_pks": [records["chemical_b"]["entity_key"]]}
    )["groups"]
    assert exact[0]["member_count"] == 1
    assert exact[0]["members"][0]["entityPk"] == records["chemical_b"]["entity_key"]


def test_auto_api_hydration_concrete_strategies_and_relationships(mixed_engine):
    engine, records = mixed_engine
    with TestClient(create_app(queries=engine)) as client:
        response = client.post("/entities/groups", json={"strategy": "auto", "member_limit": 1})
        assert response.status_code == 200
        groups = response.json()["groups"]
        chemical = groups[0]["entity"]
        assert chemical["groupMemberCount"] == 3 and len(chemical["groupMemberKeys"]) == 1
        details = client.post(
            "/entities/groups",
            json={
                "strategy": "auto",
                "group_key": chemical["entityPk"],
                "include_details": True,
                "detail_limit": 100,
            },
        ).json()["groups"][0]["entity"]
        assert details["groupDetailsLoaded"]
        assert len(details["entityAttributes"]) == 4
        assert details["sources"] == ["a", "b"]
        assert {INCHI_A, INCHI_B, "1", "Needle chemical", "Stereo B", "Alias chemical"} <= {
            i["identifier"] for i in details["identifiers"]
        }
        scoped = client.post(
            "/entities/groups",
            json={
                "strategy": "auto",
                "group_key": chemical["entityPk"],
                "include_details": True,
                "filters": {"sources": ["a"]},
                "detail_limit": 100,
            },
        ).json()["groups"][0]["entity"]
        assert scoped["sources"] == ["a"] and len(scoped["entityAttributes"]) == 4
        assert all(a["source"] == "a" for a in scoped["entityAttributes"])
        by_key = {g["entity"]["entityPk"]: g for g in groups}
        for group in (by_key[chemical["entityPk"]], by_key["gene:entrez:1"]):
            entity = group["entity"]
            # Card followups use concrete strategies, including every member.
            hydrated = client.post(
                "/entities/groups",
                json={
                    "strategy": entity["groupStrategy"],
                    "group_key": entity["entityPk"],
                    "include_details": True,
                },
            )
            assert hydrated.status_code == 200
            assert hydrated.json()["groups"][0]["member_count"] == 3
            for strategy in ("auto", entity["groupStrategy"]):
                relationships = client.post(
                    "/entities/group-relationships",
                    json={
                        "strategy": strategy,
                        "group_key": entity["entityPk"],
                        "member_keys": [],
                    },
                )
                assert relationships.status_code == 200
                assert {
                    r["relation"]["relationPk"] for r in relationships.json()["relationships"]
                } == ({"chemical", "protein"} if group["connectivity"] else {"protein"})
        singleton = client.post(
            "/entities/groups",
            json={
                "strategy": "auto",
                "group_key": "entity:" + records["native_gene"]["entity_key"],
                "include_details": True,
            },
        ).json()["groups"][0]
        assert not singleton["is_group"]
        assert singleton["entity"]["entityPk"] == records["native_gene"]["entity_key"]
        assert len(singleton["entity"]["entityAttributes"]) == 1


def test_warming_fills_the_cache_the_explorer_reads(mixed_engine):
    from omnipath_api.warm import warm

    engine, _ = mixed_engine
    client = TestClient(create_app(queries=engine))
    warm(engine)
    engine._fetch_dicts = lambda *a, **kw: (_ for _ in ()).throw(AssertionError("not warmed"))
    # The explorer's landing page reads the curated examples.
    response = client.get("/entities/examples")
    assert response.status_code == 200 and response.json()["kind"] == "examples"
    # The body InteractionFilterSidebar.svelte sends without a scope.
    response = client.post(
        "/relations/scoped-facets", json={"entityIds": [], "taxonomyLimit": 16, "taxonomyQuery": ""}
    )
    assert response.status_code == 200 and response.json()


def test_gene_group_details_combine_their_members_annotations(mixed_engine):
    # A gene's function text sits on the protein record resolved to it, not the first member.
    engine, _ = mixed_engine
    (group,) = engine.search_entity_groups(
        strategy="auto", group_key="gene:entrez:1", include_details=True
    )["groups"]
    values = {a.get("value") for a in group["entity"]["entityAttributes"]}
    members = {m["primaryIdentifier"] for m in group["members"]}
    assert {f"Detail {identifier}" for identifier in members} <= values
