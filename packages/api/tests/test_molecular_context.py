import json
import io
import pyarrow.parquet as pq
import pytest
from omnipath_api.queries.parquet import ParquetQueries
from omnipath_core.fixtures import nested_rows, rewrite_resource, write_resource


@pytest.fixture
def molecular_engine(tmp_path):
    gene, anchor, product, target = [str(i) * 64 for i in range(1, 5)]
    for source in ["a", "b"]:
        folder = tmp_path / "resources" / source / "1"
        entities = [
            dict(
                entity_key=key,
                namespace=ns,
                identifier=identifier,
                entity_type=type_,
                label=label,
                taxon="9606",
                reference_entity_key=ref,
                gene_reference_keys=[ref] if ref else [],
                identifiers=[],
                annotations=[],
            )
            for key, ns, identifier, type_, label, ref in [
                (gene, "entrez", "1", "gene", "ABC", "entrez:1"),
                (anchor, "entrez", "1", "protein", "ABC", "entrez:1"),
                (product, "uniprot", "P00001", "protein", "ABC product", "entrez:1"),
                (target, "entrez", "2", "gene", "TARGET", "entrez:2"),
            ]
        ]
        forms = [
            None,
            dict(protein_entity_key=product, isoform_identifier=dict(ns="uniprot", id="P00001-2")),
            dict(protein_entity_key=product, isoform_identifier=dict(ns="uniprot", id="P00001-3")),
        ]
        entities[1]["evidence"] = [
            dict(
                source=source,
                dataset="standalone",
                row_id=str(index),
                annotations=[],
                molecular_form=form,
            )
            for index, form in enumerate(forms)
        ]
        relation = dict(
            relation_key="relation",
            subject_entity_key=anchor,
            object_entity_key=target,
            subject_reference_entity_key="entrez:1",
            object_reference_entity_key="entrez:2",
            subject_type="protein",
            object_type="gene",
            subject_label="ABC",
            object_label="TARGET",
            predicate="affects",
            category="interaction",
            taxon="9606",
            sources=[source],
            evidence_count=3,
            annotations=[],
            evidence=[
                dict(
                    source=source,
                    row_id="1",
                    dataset="test",
                    subject_molecular_form=form,
                    object_molecular_form=None,
                    annotations=[],
                )
                for form in forms
            ],
        )
        write_resource(folder, entities, [relation])
    return ParquetQueries(tmp_path)


def test_typed_identity_and_gene_group(molecular_engine):
    engine = molecular_engine
    result = engine.search_entities_api(query="ABC")
    assert len(result["entities"]) == 3
    assert {e["entityType"] for e in result["entities"]} == {"protein", "gene"}
    group = engine.search_entity_groups(strategy="gene_reference", query="ABC", member_limit=10)[
        "groups"
    ][0]
    assert group["reference_entity_key"] == "entrez:1"
    assert group["member_count"] == 3
    assert group["entity"]["memberEntityTypes"] == ["gene", "protein"]
    assert (
        len(
            engine.search_entity_groups(
                strategy="gene_reference",
                query="ABC",
                filters={"entity_types": ["protein"]},
                member_limit=10,
            )["groups"][0]["members"]
        )
        == 2
    )


def test_exact_occurrence_filters_and_projection(molecular_engine):
    engine = molecular_engine
    filters = {"protein_entity_keys": ["3" * 64], "isoform_identifiers": ["uniprot:P00001-2"]}
    assert engine.search_relations(filters)["total"] == 2
    rows = engine.search_relations(filters, include_details=True)["rows"]
    assert all(len(r["evidence"]) == 1 for r in rows)
    assert engine.search_relations(dict(filters, molecular_endpoint_mode="target"))["total"] == 0
    assert (
        engine.search_relations(
            {"protein_entity_keys": ["3" * 64], "isoform_identifiers": ["uniprot:absent"]}
        )["total"]
        == 0
    )
    evidence = engine.get_relation_evidence("relation")["evidence"]
    assert len(evidence) == 6
    assert len({e["relationEvidencePk"] for e in evidence}) == 6
    filtered = engine.get_relation_evidence("relation", filters=filters)["evidence"]
    assert len(filtered) == 2
    assert {e["relationEvidencePk"] for e in filtered} <= {
        e["relationEvidencePk"] for e in evidence
    }
    context = engine.get_molecular_context(
        "3" * 64, view="product", isoform_identifier="uniprot:P00001-2"
    )
    assert context["relationsTotal"] == 2
    assert len(context["standaloneEvidence"]) == 2
    # The entity carries the count; its items are paged from their own endpoint.
    core = engine.get_entity_core("2" * 64)["entity"]
    assert core["molecularEvidence"] == [] and core["molecularEvidenceTotal"] == 6
    first = engine.get_entity_evidence("2" * 64, limit=4)
    rest = engine.get_entity_evidence("2" * 64, limit=4, offset=4)
    assert first["evidenceTotal"] == 6 and first["nextCursor"] == "4"
    assert len(first["evidence"]) == 4 and len(rest["evidence"]) == 2 and rest["nextCursor"] is None
    assert len({json.dumps(e, sort_keys=True) for e in first["evidence"] + rest["evidence"]}) == 6
    assert all(len(r["evidence"]) == 1 for r in context["relations"])
    gene_context = engine.get_molecular_context("gene:entrez:1")
    assert len(gene_context["catalogueProducts"]) == 1
    assert gene_context["filters"] == {"reference_entity_keys": ["entrez:1"]}
    assert gene_context["relationsTotal"] == 2
    assert len(gene_context["standaloneEvidence"]) == 6
    assert engine.get_molecular_context("1" * 64)["relationsTotal"] == 2
    data, _ = engine.export_slice(filters)
    exported = pq.read_table(io.BytesIO(data)).to_pylist()
    assert len(exported) == 2
    assert all(len(r["evidence"]) == 1 for r in exported)
    assert {p["entity_key"] for p in exported[0]["referenced_product_records"]} == {"3" * 64}


def test_opposite_endpoint_and_different_occurrences_do_not_combine(molecular_engine):
    engine = molecular_engine
    assert (
        engine.search_relations(
            {
                "protein_entity_keys": ["3" * 64],
                "isoform_identifiers": ["uniprot:P00001-2"],
                "molecular_endpoint_mode": "both",
            }
        )["total"]
        == 0
    )


def test_http_molecular_contract(molecular_engine):
    from fastapi.testclient import TestClient
    from omnipath_api.server import create_app

    client = TestClient(create_app(queries=molecular_engine))
    response = client.get("/relations/relation/evidence")
    assert response.status_code == 200
    assert len(response.json()["evidence"]) == 6
    context = client.get(
        "/entities/" + ("3" * 64) + "/molecular-context",
        params={"view": "product", "isoform_identifier": "uniprot:P00001-2"},
    )
    assert context.status_code == 200 and context.json()["relationsTotal"] == 2
    groups = client.post("/entities/groups", json={"strategy": "gene_reference", "query": "ABC"})
    assert groups.status_code == 200 and groups.json()["groups"][0]["entity"][
        "memberEntityTypes"
    ] == ["gene", "protein"]


def test_filters_do_not_combine_endpoints_or_occurrences(molecular_engine):
    engine = molecular_engine
    for path in engine.data_root.glob("resources/*/*/relation.parquet"):
        original = nested_rows(path)[0]
        split = dict(
            original,
            relation_key="split",
            evidence=[
                dict(
                    source="test",
                    row_id="x",
                    subject_molecular_form=dict(
                        protein_entity_key="3" * 64,
                        isoform_identifier=dict(ns="uniprot", id="P00001-2"),
                    ),
                    object_molecular_form=None,
                    annotations=[],
                ),
                dict(
                    source="test",
                    row_id="y",
                    subject_molecular_form=dict(
                        protein_entity_key="5" * 64,
                        isoform_identifier=dict(ns="uniprot", id="P00001-3"),
                    ),
                    object_molecular_form=None,
                    annotations=[],
                ),
                dict(
                    source="test",
                    row_id="z",
                    subject_molecular_form=dict(protein_entity_key="3" * 64),
                    object_molecular_form=dict(
                        isoform_identifier=dict(ns="uniprot", id="P00001-3")
                    ),
                    annotations=[],
                ),
            ],
        )
        rewrite_resource(path, relations=[original, split])
    engine.reload_resources()
    filters = {"protein_entity_keys": ["3" * 64], "isoform_identifiers": ["uniprot:P00001-3"]}
    assert {r["relation_key"] for r in engine.search_relations(filters)["rows"]} == {"relation"}
    assert (
        engine.search_relations(
            {"reference_entity_keys": ["entrez:1"], "entity_types": ["rna_product"]}
        )["total"]
        == 0
    )
    assert engine.search_relations({"entity_ids": ["entrez:1"]})["total"] == 4


def test_paging_preserves_exact_evidence(molecular_engine):
    engine = molecular_engine
    filters = {"protein_entity_keys": ["3" * 64], "isoform_identifiers": ["uniprot:P00001-2"]}
    result = engine.search_relations(filters, limit=1, include_details=True)
    assert result["total"] == 2
    assert len(result["rows"][0]["evidence"]) == 1
    api = engine.search_relations_api(filters, limit=1)
    assert api["rows"][0]["relation"]["evidenceCount"] == 1
    assert engine.search_relations(filters, limit=1, offset=50)["total"] == 2


def test_fallback_references_do_not_become_gene_groups(molecular_engine):
    engine = molecular_engine
    records = [
        dict(
            entity_key="5" * 64,
            entity_type="small_molecule",
            namespace="chebi",
            identifier="1",
            label="chemical",
            reference_entity_key="chebi:1",
            gene_reference_keys=[],
        ),
        dict(
            entity_key="6" * 64,
            entity_type="protein",
            namespace="uniprot",
            identifier="P99999",
            label="unresolved protein",
            reference_entity_key="uniprot:P99999",
            gene_reference_keys=[],
        ),
        dict(
            entity_key="7" * 64,
            entity_type="protein",
            namespace="uniprot",
            identifier="P88888",
            label="ambiguous protein",
            reference_entity_key="uniprot:P88888",
            gene_reference_keys=["entrez:8", "entrez:9"],
        ),
    ]
    for path in engine.data_root.glob("resources/*/*/entity.parquet"):
        rewrite_resource(path, entities=nested_rows(path) + records)
    engine.reload_resources()
    groups = engine.search_entity_groups(strategy="gene_reference", limit=20)["groups"]
    by_key = {g["entity"]["entityPk"]: g for g in groups}
    assert "gene:entrez:1" in by_key
    for record in records:
        singleton = by_key[record["entity_key"]]
        assert singleton["is_group"] is False
        assert singleton["entity"]["entityType"] == record["entity_type"]
        assert singleton["entity"]["referenceEntityKey"] == record["reference_entity_key"]
        assert singleton["entity"]["canonicalIdentifierType"] == record["namespace"]


def test_product_context_keeps_native_fallback_without_losing_exact_observations(molecular_engine):
    engine = molecular_engine
    # One source keeps the native fallback while another has a supported gene
    # reference for the same reusable product. Both have real anchor evidence.
    path = next(iter(engine._table_paths("entity", ["a"]).values()))
    records = nested_rows(path)
    product = next(row for row in records if row["entity_key"] == "3" * 64)
    product["reference_entity_key"] = "uniprot:P00001"
    product["label"] = "ZZZ fallback"  # Reference choice must not follow display rank.
    rewrite_resource(path, entities=records)
    engine.reload_resources()
    summary = engine.get_entity_core("3" * 64)["entity"]
    assert summary["referenceEntityKey"] == "uniprot:P00001"
    assert summary["geneReferenceKeys"] == ["entrez:1"]
    context = engine.get_molecular_context(
        "3" * 64, view="product", isoform_identifier="uniprot:P00001-2"
    )
    assert context["referenceEntityKey"] == "uniprot:P00001"
    assert context["relationsTotal"] == 2
    assert len(context["standaloneEvidence"]) == 2
    assert len(context["observedForms"]) == 1
    groups = engine.search_entity_groups(strategy="gene_reference", limit=20, member_limit=10)[
        "groups"
    ]
    singleton = next(group for group in groups if group["entity"]["entityPk"] == "3" * 64)
    assert singleton["is_group"] is False
    assert singleton["entity"]["entityType"] == "protein"
    assert singleton["entity"]["referenceEntityKey"] == "uniprot:P00001"


@pytest.mark.parametrize(
    "entity_type,namespace,identifier",
    [
        ("physical_entity", "uniprot", "P04637"),
        ("gene", "uniprot", "P04637"),
        ("gene", "enst", "ENST000001"),
        ("protein", "enst", "ENST000001"),
        ("rna_product", "uniprot", "P04637"),
    ],
)
def test_non_product_source_types_keep_native_context(
    molecular_engine, entity_type, namespace, identifier
):
    from fastapi.testclient import TestClient
    from omnipath_api.server import create_app

    engine = molecular_engine
    key, reference = "5" * 64, f"{namespace}:{identifier}"
    other_key = "6" * 64
    other_type = "gene" if entity_type != "gene" else "physical_entity"
    for path in engine.data_root.glob("resources/*/*/entity.parquet"):
        source = path.parent.parent.name
        record = dict(
            entity_key=key,
            entity_type=entity_type,
            namespace=namespace,
            identifier=identifier,
            label="source-typed record",
            reference_entity_key=reference,
            gene_reference_keys=[],
            evidence=[
                dict(
                    source=source,
                    dataset="standalone",
                    row_id="native",
                    annotations=[],
                    molecular_form=None,
                )
            ],
        )
        other = dict(
            record,
            entity_key=other_key,
            entity_type=other_type,
            label="different typed singleton",
            evidence=[
                dict(
                    source=source,
                    dataset="standalone",
                    row_id="other",
                    annotations=[],
                    molecular_form=None,
                )
            ],
        )
        rewrite_resource(path, entities=nested_rows(path) + [record, other])
    for path in engine.data_root.glob("resources/*/*/relation.parquet"):
        rows = nested_rows(path)
        source = path.parent.parent.name
        native = dict(
            rows[0],
            relation_key="native",
            subject_entity_key=key,
            subject_type=entity_type,
            subject_reference_entity_key=reference,
            subject_label="source-typed record",
            evidence_count=1,
            evidence=[
                dict(
                    source=source,
                    row_id="native",
                    annotations=[],
                    subject_molecular_form=None,
                    object_molecular_form=None,
                )
            ],
        )
        other = dict(
            native,
            relation_key="other_native",
            subject_entity_key=other_key,
            subject_type=other_type,
            subject_label="different typed singleton",
            evidence=[
                dict(
                    source=source,
                    row_id="other",
                    annotations=[],
                    subject_molecular_form=None,
                    object_molecular_form=None,
                )
            ],
        )
        rewrite_resource(path, relations=rows + [native, other])
    engine.reload_resources()
    context = engine.get_molecular_context(key)
    assert context["filters"] == {"entity_pks": [key]}
    assert context["relationsTotal"] == 2
    assert context["observedForms"] == []
    assert len(context["standaloneEvidence"]) == 2
    assert all(record["entityPk"] == key for record in context["standaloneEvidence"])
    assert all(row["subjectEntity"]["entityType"] == entity_type for row in context["relations"])
    assert all(row["relation"]["subjectEntityPk"] == key for row in context["relations"])
    assert all(
        len(row["evidence"]) == 1 and row["evidence"][0]["subject_molecular_form"] is None
        for row in context["relations"]
    )
    with TestClient(create_app(queries=engine)) as client:
        response = client.get(f"/entities/{key}/molecular-context")
        assert response.status_code == 200
        assert response.json()["relationsTotal"] == 2
        assert response.json()["filters"] == {"entity_pks": [key]}


@pytest.mark.parametrize(
    "entity_type,namespace,identifier,other_identifier,product_kind",
    [
        ("protein", "uniprot", "P04637", "P00533", "protein"),
        ("protein", "ensp", "ENSP000001", "ENSP000002", "protein"),
        ("protein", "refseq", "NP_000001", "NP_000002", "protein"),
        ("protein", "refseq_protein", "NP_000001.2", "NP_000002.3", "protein"),
        ("transcript", "enst", "ENST000001", "ENST000002", "transcript"),
        ("transcript", "refseq", "NM_000001", "NM_000002", "transcript"),
    ],
)
def test_compatible_product_type_and_namespace_keep_exact_context(
    molecular_engine, entity_type, namespace, identifier, other_identifier, product_kind
):
    engine = molecular_engine
    product_key, anchor_key, other_key = [str(i) * 64 for i in (5, 6, 7)]
    forms = [
        None,
        {product_kind + "_entity_key": product_key},
        {product_kind + "_entity_key": other_key},
    ]
    for path in engine.data_root.glob("resources/*/*/entity.parquet"):
        source = path.parent.parent.name
        records = [
            dict(
                entity_key=key,
                entity_type=entity_type,
                namespace=ns,
                identifier=id_,
                label=id_,
                reference_entity_key="entrez:7",
                gene_reference_keys=["entrez:7"],
                evidence=[],
            )
            for key, ns, id_ in [
                (product_key, namespace, identifier),
                (anchor_key, "entrez", "7"),
                (other_key, namespace, other_identifier),
            ]
        ]
        records[1]["evidence"] = [
            dict(source=source, row_id=str(index), annotations=[], molecular_form=form)
            for index, form in enumerate(forms)
        ]
        rewrite_resource(path, entities=nested_rows(path) + records)
    for path in engine.data_root.glob("resources/*/*/relation.parquet"):
        rows = nested_rows(path)
        source = path.parent.parent.name
        product_relation = dict(
            rows[0],
            relation_key="product",
            subject_entity_key=anchor_key,
            subject_type=entity_type,
            subject_reference_entity_key="entrez:7",
            subject_label=identifier,
            evidence_count=3,
            evidence=[
                dict(
                    source=source,
                    row_id=str(index),
                    annotations=[],
                    subject_molecular_form=form,
                    object_molecular_form=None,
                )
                for index, form in enumerate(forms)
            ],
        )
        rewrite_resource(path, relations=rows + [product_relation])
    engine.reload_resources()
    context = engine.get_molecular_context(product_key, view="product")
    assert context["filters"] == {product_kind + "_entity_keys": [product_key]}
    assert context["relationsTotal"] == 2
    assert len(context["standaloneEvidence"]) == 2
    assert len(context["observedForms"]) == 1
    assert all(
        len(row["evidence"]) == 1
        and row["evidence"][0]["subject_molecular_form"][product_kind + "_entity_key"]
        == product_key
        for row in context["relations"]
    )
    assert all(
        record["occurrence"]["molecular_form"][product_kind + "_entity_key"] == product_key
        for record in context["standaloneEvidence"]
    )


@pytest.mark.parametrize("relation_count", [0, 2, 31])
def test_shared_cursor_pages_both_streams_to_completion(molecular_engine, relation_count):
    from fastapi.testclient import TestClient
    from omnipath_api.server import create_app

    engine = molecular_engine
    folder = engine.data_root / "resources" / "a" / "1"
    entity_path, relation_path = folder / "entity.parquet", folder / "relation.parquet"
    entities = nested_rows(entity_path)
    # Nullable/repeated occurrence metadata must not make pagination unstable.
    entities[1]["evidence"] = [
        dict(source="a", row_id=str(index), molecular_form=None, annotations=[])
        for index in range(27)
    ]
    rewrite_resource(entity_path, entities=entities)
    relation = nested_rows(relation_path)[0]
    rewrite_resource(
        relation_path,
        relations=[
            dict(relation, relation_key=f"relation-{index}") for index in range(relation_count)
        ],
    )
    other_folder = engine.data_root / "resources" / "b" / "1"
    other_entities = nested_rows(other_folder / "entity.parquet")
    for record in other_entities:
        record["evidence"] = []
    rewrite_resource(other_folder / "entity.parquet", entities=other_entities)
    rewrite_resource(other_folder / "relation.parquet", relations=[])
    engine.reload_resources()
    with TestClient(create_app(queries=engine)) as client:
        endpoint = "/entities/gene:entrez:1/molecular-context"
        first = client.get(endpoint).json()
        assert len(first["standaloneEvidence"]) == 20
        assert first["nextCursor"] == "20"
        assert client.get(endpoint).json()["standaloneEvidence"] == first["standaloneEvidence"]
        second = client.get(endpoint, params={"offset": 20}).json()
        assert len(second["standaloneEvidence"]) == 7
        assert second["nextCursor"] is None
        assert len(first["relations"]) + len(second["relations"]) == relation_count
        observations = first["standaloneEvidence"] + second["standaloneEvidence"]
        assert [record["occurrence"]["row_id"] for record in observations] == [
            str(index) for index in range(27)
        ]


def test_shared_cursor_continues_standalone_after_cross_resource_relations_finish(molecular_engine):
    engine = molecular_engine
    pages = [engine.get_molecular_context("gene:entrez:1", limit=2, offset=i) for i in (0, 2, 4)]
    assert [len(page["relations"]) for page in pages] == [2, 0, 0]
    assert [len(page["standaloneEvidence"]) for page in pages] == [2, 2, 2]
    assert [page["nextCursor"] for page in pages] == ["2", "4", None]
    assert [
        (record["occurrence"]["source"], record["occurrence"]["row_id"])
        for page in pages
        for record in page["standaloneEvidence"]
    ] == [(source, str(index)) for source in ("a", "b") for index in range(3)]


@pytest.mark.parametrize("namespace,identifier", [("uniprot", "UNKNOWN"), ("ensp", "ENSP000001.2")])
def test_native_product_without_catalogue_links_finds_gene_standalone_occurrences(
    molecular_engine, namespace, identifier
):
    engine = molecular_engine
    product, main_gene, other_product = [str(i) * 64 for i in (5, 6, 7)]
    isoform = {"ns": namespace, "id": identifier + "-2"}
    forms = [
        {"protein_entity_key": product, "isoform_identifier": isoform},
        {"protein_entity_key": product, "isoform_identifier": {"ns": namespace, "id": "other"}},
        {"protein_entity_key": other_product, "isoform_identifier": isoform},
        None,
    ]
    for path in engine.data_root.glob("resources/*/*/entity.parquet"):
        source = path.parent.parent.name
        records = [
            dict(
                entity_key=key,
                namespace=namespace,
                identifier=id_,
                entity_type="protein",
                reference_entity_key=f"{namespace}:{id_}",
                gene_reference_keys=[],
                evidence=[],
            )
            for key, id_ in [(product, identifier), (other_product, "other")]
        ]
        records.append(
            dict(
                entity_key=main_gene,
                namespace="entrez",
                identifier="7",
                entity_type="protein",
                reference_entity_key="entrez:7",
                gene_reference_keys=["entrez:7"],
                evidence=[
                    dict(source=source, row_id=str(index), molecular_form=form, annotations=[])
                    for index, form in enumerate(forms)
                ],
            )
        )
        rewrite_resource(path, entities=nested_rows(path) + records)
    engine.reload_resources()
    context = engine.get_molecular_context(
        product, view="product", isoform_identifier=f"{namespace}:{isoform['id']}"
    )
    assert context["filters"] == {
        "protein_entity_keys": [product],
        "isoform_identifiers": [f"{namespace}:{isoform['id']}"],
    }
    assert context["referenceEntityKey"] == f"{namespace}:{identifier}"
    assert context["catalogueProducts"] == []
    assert context["relationsTotal"] == 0
    assert len(context["standaloneEvidence"]) == 2
    assert context["nextCursor"] is None
    assert {record["entityPk"] for record in context["standaloneEvidence"]} == {main_gene}
    assert all(record["occurrence"]["row_id"] == "0" for record in context["standaloneEvidence"])
    assert len(engine.get_molecular_context(product, view="product")["standaloneEvidence"]) == 4
    assert len(engine.get_molecular_context("gene:entrez:7")["standaloneEvidence"]) == 8


@pytest.mark.parametrize(
    "source_type,namespace,identifier,product_type",
    [
        ("protein", "ensembl", "ENSP999999999999.73", "protein"),
        ("rna_product", "ensembl", "ENST999999999999.73", "transcript"),
        ("protein", "uniprot", "Q9Y6K9", "protein"),
        ("protein", "refseq_protein", "NP_999999999999.73", "protein"),
    ],
)
def test_actual_writer_native_source_and_exact_product_views(
    tmp_path, source_type, namespace, identifier, product_type
):
    from omnipath_resolver.resolver import EntityResolver
    from omnipath_build.silver import SilverExtractor
    from omnipath_build.writer import ParquetWriter
    from fastapi.testclient import TestClient
    from omnipath_api.server import create_app

    item = {"type": source_type, "identifiers": [{"type": namespace, "value": identifier}]}
    extractor = SilverExtractor("fixture", "native-product")
    extractor.process_record(item, item, "standalone", 0)
    relation = {
        "subject": item,
        "predicate": "affects",
        "object": {"type": "gene", "identifiers": [{"type": "entrez", "value": "55"}]},
    }
    extractor.process_record(relation, relation, "relation", 1)
    resolver = EntityResolver(library_dir=tmp_path / "absent")
    writer = ParquetWriter(tmp_path / "resources" / "fixture" / "1")
    try:
        writer.append_observations(extractor, resolver)
        files = writer.close()["files"]
        entity_path, relation_path = files["entity"], files["relation"]
    finally:
        resolver.close()
    (row,) = nested_rows(relation_path)
    main_key = row["subject_entity_key"]
    field = product_type + "_entity_key"
    product_key = row["evidence"][0]["subject_molecular_form"][field]
    entities = {entity["entity_key"]: entity for entity in nested_rows(entity_path)}
    assert entities[main_key]["entity_type"] == source_type
    assert entities[product_key]["entity_type"] == product_type
    assert entities[product_key]["identifier"] == identifier
    if namespace == "ensembl":
        assert main_key != product_key
        assert entities[main_key]["identifier"] == identifier.split(".")[0]
    else:
        assert main_key == product_key
    engine = ParquetQueries(tmp_path)
    source_context = engine.get_molecular_context(main_key)
    assert source_context["filters"] == {"entity_pks": [main_key]}
    assert source_context["relationsTotal"] == 1
    assert source_context["relations"][0]["subjectEntity"]["entityType"] == source_type
    assert len(source_context["standaloneEvidence"]) == 1
    assert source_context["standaloneEvidence"][0]["entityPk"] == main_key
    product_context = engine.get_molecular_context(product_key, view="product")
    assert product_context["filters"] == {product_type + "_entity_keys": [product_key]}
    assert product_context["relationsTotal"] == 1
    assert len(product_context["standaloneEvidence"]) == 1
    assert product_context["standaloneEvidence"][0]["entityPk"] == main_key
    with TestClient(create_app(queries=engine)) as client:
        assert client.get(f"/entities/{main_key}/molecular-context").json()["filters"] == {
            "entity_pks": [main_key]
        }
        response = client.get(
            f"/entities/{product_key}/molecular-context", params={"view": "product"}
        )
        assert response.status_code == 200
        assert response.json()["relationsTotal"] == 1
        if main_key != product_key:
            assert (
                client.get(f"/entities/{product_key}/molecular-context").json()["relationsTotal"]
                == 0
            )
        if source_type == "rna_product":
            assert (
                client.get(
                    f"/entities/{main_key}/molecular-context", params={"view": "product"}
                ).status_code
                == 400
            )


def test_shared_native_product_key_does_not_turn_legacy_evidence_into_exact_product(
    molecular_engine,
):
    engine = molecular_engine
    key = "3" * 64
    for path in engine.data_root.glob("resources/*/*/entity.parquet"):
        records = nested_rows(path)
        source = path.parent.parent.name
        for record in records:
            record["evidence"] = []
        product = next(record for record in records if record["entity_key"] == key)
        product["reference_entity_key"] = "uniprot:P00001"
        product["gene_reference_keys"] = []
        product["evidence"] = [
            dict(
                source=source,
                row_id="native",
                molecular_form={"protein_entity_key": key} if source == "a" else None,
            )
        ]
        rewrite_resource(path, entities=records)
    for path in engine.data_root.glob("resources/*/*/relation.parquet"):
        records = nested_rows(path)
        source = path.parent.parent.name
        records[0]["subject_entity_key"] = key
        records[0]["subject_reference_entity_key"] = "uniprot:P00001"
        records[0]["evidence"] = [
            dict(
                source=source,
                row_id="native",
                subject_molecular_form={"protein_entity_key": key} if source == "a" else None,
            )
        ]
        records[0]["evidence_count"] = 1
        rewrite_resource(path, relations=records)
    engine.reload_resources()
    reference = engine.get_molecular_context(key)
    assert reference["filters"] == {"entity_pks": [key]}
    assert reference["relationsTotal"] == 2
    assert len(reference["standaloneEvidence"]) == 2
    product = engine.get_molecular_context(key, view="product")
    assert product["relationsTotal"] == 1
    assert len(product["standaloneEvidence"]) == 1
    assert product["standaloneEvidence"][0]["occurrence"]["source"] == "a"
    assert product["relations"][0]["evidence"][0]["source"] == "a"


def test_isoform_selection_requires_explicit_product_view(molecular_engine):
    from fastapi.testclient import TestClient
    from omnipath_api.server import create_app

    engine = molecular_engine
    with pytest.raises(ValueError, match="Isoform selection requires product view"):
        engine.get_molecular_context("gene:entrez:1", isoform_identifier="uniprot:P00001-2")
    with TestClient(create_app(queries=engine)) as client:
        endpoint = f"/entities/{'3' * 64}/molecular-context"
        assert (
            client.get(endpoint, params={"isoform_identifier": "uniprot:P00001-2"}).status_code
            == 400
        )
        assert (
            client.get(
                endpoint, params={"view": "reference", "isoform_identifier": "uniprot:P00001-2"}
            ).status_code
            == 400
        )
        assert client.get(endpoint, params={"view": "other"}).status_code == 422
        assert (
            client.get(
                "/entities/gene:entrez:1/molecular-context", params={"view": "product"}
            ).status_code
            == 400
        )
        exact = client.get(
            endpoint, params={"view": "product", "isoform_identifier": "uniprot:P00001-2"}
        ).json()
        assert exact["relationsTotal"] == 2
        assert len(exact["standaloneEvidence"]) == 2


@pytest.mark.parametrize("include_relation", [False, True])
def test_actual_writer_native_product_with_gene_is_labeled_before_navigation(
    tmp_path, include_relation
):
    import sys
    from pathlib import Path

    # Import under the bare module name the build tests use, so both share one
    # process-wide reference template (the "packages.build.tests." spelling would
    # be a second module object with its own cache, i.e. a second library build).
    build_tests = str(Path(__file__).resolve().parents[2] / "build" / "tests")
    if build_tests not in sys.path:
        sys.path.insert(0, build_tests)
    from library_fixture import build_fixture_library
    from omnipath_resolver.resolver import EntityResolver
    from omnipath_build.silver import SilverExtractor
    from omnipath_build.writer import ParquetWriter

    item = {
        "type": "protein",
        "identifiers": [
            {"type": "uniprot", "value": "Q9Y6K9"},
            {"type": "entrez", "value": "7157"},
        ],
    }
    extractor = SilverExtractor("fixture", "native-label")
    extractor.process_record(item, item, "standalone", 0)
    if include_relation:
        relation = {
            "subject": item,
            "predicate": "affects",
            "object": {
                "type": "gene",
                "identifiers": [{"type": "entrez", "value": "55"}],
            },
        }
        extractor.process_record(relation, relation, "relation", 1)
    library = build_fixture_library(tmp_path / "reference")
    resolver = EntityResolver(library_dir=library)
    writer = ParquetWriter(tmp_path / "resources" / "fixture" / "1", library_dir=library)
    try:
        writer.append_observations(extractor, resolver)
        writer.close()
    finally:
        resolver.close()
    engine = ParquetQueries(tmp_path)
    context = engine.get_molecular_context("gene:entrez:7157")
    assert context["relationsTotal"] == int(include_relation)
    assert context["catalogueProducts"] == []
    assert len(context["standaloneEvidence"]) == 1
    assert len(context["referencedProducts"]) == 1
    assert len(context["observedForms"]) == 1
    product = context["referencedProducts"][0]
    assert context["observedForms"][0]["protein_entity_key"] == product["entityPk"]
    assert product["canonicalIdentifier"] == product["label"] == "Q9Y6K9"
    assert product["entityType"] == "protein"
    assert product["geneReferenceKeys"] == []
    assert (
        product["entityPk"]
        == context["standaloneEvidence"][0]["occurrence"]["molecular_form"]["protein_entity_key"]
    )


def test_referenced_product_hydration_is_one_bounded_deduplicated_scalar_batch(
    molecular_engine, monkeypatch
):
    from omnipath_core.keys import entity_key

    engine = molecular_engine
    records = [
        dict(
            entity_key=entity_key("protein", "uniprot", f"P9{index:04d}"),
            entity_type="protein",
            namespace="uniprot",
            identifier=f"P9{index:04d}",
            reference_entity_key=f"uniprot:P9{index:04d}",
            gene_reference_keys=[],
            evidence=[],
        )
        for index in range(101)
    ]
    transcript = dict(
        entity_key=entity_key("transcript", "enst", "ENST999999999999.2"),
        entity_type="transcript",
        namespace="enst",
        identifier="ENST999999999999.2",
        reference_entity_key="enst:ENST999999999999.2",
        gene_reference_keys=[],
        evidence=[],
    )
    for path in engine.data_root.glob("resources/*/*/entity.parquet"):
        entities = nested_rows(path)
        entities[1]["evidence"] = [
            dict(
                source=path.parent.parent.name,
                annotations=[],
                molecular_form={
                    "protein_entity_key": records[0]["entity_key"],
                    "transcript_entity_key": transcript["entity_key"],
                },
            )
        ] * 2
        rewrite_resource(path, entities=entities + records + [transcript])
    for path in engine.data_root.glob("resources/*/*/relation.parquet"):
        rows = nested_rows(path)
        rows[0]["evidence"] = [
            dict(
                source=path.parent.parent.name,
                subject_molecular_form={"protein_entity_key": record["entity_key"]},
                object_molecular_form={"transcript_entity_key": transcript["entity_key"]},
            )
            for record in records
        ]
        rows[0]["evidence_count"] = len(records)
        rewrite_resource(path, relations=rows)
    engine.reload_resources()
    original = engine._fetch_entities_by_keys
    calls = []

    def lookup(keys, resources=None, *, slim=False):
        calls.append((keys, resources, slim))
        return original(keys, resources, slim=slim)

    monkeypatch.setattr(engine, "_fetch_entities_by_keys", lookup)
    context = engine.get_molecular_context("gene:entrez:1", resources=["a"], limit=1)
    assert len(calls) == 1  # Endpoints come with their relations; one product batch.
    keys, resources, slim = calls[0]
    assert resources == ["a"] and slim
    assert len(keys) == len(set(keys)) == 100
    assert keys[:2] == [records[0]["entity_key"], transcript["entity_key"]]
    assert len(context["referencedProducts"]) == 100
    assert {product["entityPk"] for product in context["referencedProducts"]} == set(keys)
    assert not set(keys).intersection(
        product["entityPk"] for product in context["catalogueProducts"]
    )


def test_referenced_product_metadata_only_follows_returned_trimmed_pairs_and_standalone_page(
    molecular_engine,
):
    from omnipath_core.keys import entity_key

    engine = molecular_engine
    keys = [entity_key("protein", "uniprot", f"Q9Y6K{index}") for index in range(4, 8)]
    native = [
        dict(
            entity_key=key,
            entity_type="protein",
            namespace="uniprot",
            identifier=f"Q9Y6K{index}",
            reference_entity_key=f"uniprot:Q9Y6K{index}",
            gene_reference_keys=[],
            evidence=[],
        )
        for index, key in zip(range(4, 8), keys)
    ]

    def form(key, iso="Q9Y6K4-2"):
        return dict(protein_entity_key=key, isoform_identifier={"ns": "uniprot", "id": iso})

    for path in engine.data_root.glob("resources/*/*/entity.parquet"):
        entities = nested_rows(path)
        for entity in entities:
            entity["evidence"] = []
        entities[1]["evidence"] = [
            dict(source=path.parent.parent.name, molecular_form=form(key)) for key in keys[:2]
        ]
        rewrite_resource(path, entities=entities + native)
    for path in engine.data_root.glob("resources/*/*/relation.parquet"):
        original = nested_rows(path)[0]
        source = path.parent.parent.name
        rows = [
            dict(
                original,
                relation_key="first",
                evidence=[
                    dict(
                        source=source,
                        subject_molecular_form=form(keys[0]),
                        object_molecular_form=form(keys[2]),
                    ),
                    dict(
                        source=source,
                        subject_molecular_form=form(keys[1]),
                        object_molecular_form=form(keys[3]),
                    ),
                    dict(
                        source=source,
                        subject_molecular_form=form(keys[0], "Q9Y6K4-3"),
                        object_molecular_form=form(keys[3]),
                    ),
                ],
                evidence_count=3,
            ),
            dict(
                original,
                relation_key="second",
                evidence=[
                    dict(
                        source=source,
                        subject_molecular_form=form(keys[0]),
                        object_molecular_form=form(keys[3]),
                    ),
                ],
                evidence_count=1,
            ),
        ]
        rewrite_resource(path, relations=rows)
    engine.reload_resources()
    standalone = engine.get_molecular_context("gene:entrez:1", resources=["a"], limit=1, offset=2)
    assert standalone["standaloneEvidence"] == []
    assert standalone["relations"] == []
    assert standalone["referencedProducts"] == []
    standalone_first = engine.get_molecular_context("gene:entrez:1", resources=["a"], limit=1)
    assert len(standalone_first["standaloneEvidence"]) == 1
    context = engine.get_molecular_context(
        keys[0], resources=["a"], view="product", isoform_identifier="uniprot:Q9Y6K4-2", limit=1
    )
    assert context["relationsTotal"] == 2 and len(context["relations"]) == 1
    assert len(context["standaloneEvidence"]) == 1
    assert len(context["relations"][0]["evidence"]) == 1
    expected = {
        form["protein_entity_key"]
        for row in context["relations"]
        for ev in row["evidence"]
        for form in (ev["subject_molecular_form"], ev["object_molecular_form"])
    }
    assert {product["entityPk"] for product in context["referencedProducts"]} == expected
    assert keys[1] not in expected
    # A lookahead observation must not hydrate products from the next page.
    for path in engine.data_root.glob("resources/*/*/relation.parquet"):
        rewrite_resource(path, relations=[])
    engine.reload_resources()
    for offset, key in enumerate(keys[:2]):
        page = engine.get_molecular_context(
            "gene:entrez:1", resources=["a"], limit=1, offset=offset
        )
        assert page["relations"] == []
        assert {product["entityPk"] for product in page["referencedProducts"]} == {key}
        assert {form["protein_entity_key"] for form in page["observedForms"]} == {key}
