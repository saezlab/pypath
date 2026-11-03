"""Exercise the real-output validator with valid and deliberately damaged Parquets."""

import copy
import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from omnipath_core.keys import entity_key
from scripts.validate_molecular_outputs import Validator, validate_outputs
from omnipath_core.fixtures import nested_rows, rewrite_resource, write_resource


def entity(rows, entity_type, identifier, namespace=None):
    """The fixture entity row of a type and identifier (tables are sorted by key)."""
    return next(
        row
        for row in rows
        if row["entity_type"] == entity_type
        and row["identifier"] == identifier
        and namespace in (None, row["namespace"])
    )


@pytest.fixture
def outputs(tmp_path):
    folder = tmp_path / "resources" / "signor" / "1"
    rows = []
    for type_, ns, identifier, ref in [
        ("gene", "entrez", "7157", "entrez:7157"),
        ("protein", "entrez", "7157", "entrez:7157"),
        ("protein", "uniprot", "P04637", "entrez:7157"),
        ("protein", "entrez", "1956", "entrez:1956"),
        ("protein", "uniprot", "P00533", "entrez:1956"),
        # Catalogue-only links must not turn an ambiguous fallback into a gene group.
        ("protein", "uniprot", "P00001", "uniprot:P00001"),
        ("chemical", "chebi", "1234", "chebi:1234"),
    ]:
        rows.append(
            dict(
                entity_key=entity_key(type_, ns, identifier),
                entity_type=type_,
                namespace=ns,
                identifier=identifier,
                label=identifier,
                taxon="9606",
                reference_entity_key=ref,
                gene_reference_keys=[ref]
                if ref.startswith("entrez:")
                else ["entrez:7157", "entrez:1956"]
                if ns == "uniprot"
                else [],
                identifiers=[],
                annotations=[],
                evidence=[],
            )
        )
    product, counterpart = rows[2]["entity_key"], rows[4]["entity_key"]
    form = dict(
        protein_entity_key=product,
        isoform_identifier=dict(ns="uniprot", id="P04637-2"),
        sequence_identifiers=[dict(ns="uniprot", id="P04637-2")],
        modifications=[
            dict(
                term="MOD:00696",
                residue="S",
                position=15,
                coordinate_reference=dict(
                    identifier=dict(ns="uniprot", id="P04637-2"),
                    coordinate_system="protein",
                    position_base=1,
                ),
            )
        ],
    )
    other = dict(
        protein_entity_key=counterpart, isoform_identifier=dict(ns="uniprot", id="P00533-3")
    )
    rows[1]["evidence"] = [
        dict(
            source="signor",
            dataset="annotations",
            row_id="standalone",
            annotations=[],
            molecular_form=form,
        )
    ]
    evidence = [
        dict(
            source="signor",
            dataset="interactions",
            row_id="row-1",
            upstream_id="source-1",
            annotations=[],
            subject_molecular_form=form,
            object_molecular_form=other,
        ),
        dict(
            source="signor",
            dataset="interactions",
            row_id="row-2",
            annotations=[],
            subject_molecular_form=dict(
                protein_entity_key=product, isoform_identifier=dict(ns="uniprot", id="P04637-4")
            ),
            object_molecular_form=None,
        ),
        dict(
            source="signor",
            dataset="interactions",
            row_id="row-3",
            annotations=[],
            subject_molecular_form=None,
            object_molecular_form=None,
        ),
    ]
    relation = dict(
        relation_key="relation",
        subject_entity_key=rows[1]["entity_key"],
        subject_type="protein",
        subject_label="TP53",
        object_entity_key=rows[3]["entity_key"],
        object_type="protein",
        object_label="EGFR",
        predicate="biolink:affects",
        subject_reference_entity_key="entrez:7157",
        object_reference_entity_key="entrez:1956",
        category="interaction",
        statement_kind="relation",
        sources=["signor"],
        taxon="9606",
        evidence=evidence,
        evidence_count=len(evidence),
        annotations=[],
    )
    write_resource(folder, rows, [relation])
    return tmp_path, folder


def failures(report):
    return [c for r in report["resources"] for c in r["checks"] if c["status"] == "failed"]


def test_progress_logging_omits_query_results_and_parameters(outputs, monkeypatch, capsys):
    monkeypatch.setenv("OMNIPATH_VALIDATION_PROGRESS", "1")
    validator = Validator(outputs[1], examples=1, max_export_relations=100)
    try:
        assert validator.sql("fixture.operation", "SELECT ?", ["private-value"]) == [
            ("private-value",)
        ]
    finally:
        validator.db.close()
    output = capsys.readouterr().err
    assert "private-value" not in output
    assert set(json.loads(output)) == {"resource", "operation", "kind", "elapsed_ms"}


def test_validator_paired_forms_gene_navigation_export_and_unavailable_cases(outputs):
    root, _ = outputs
    report = validate_outputs(root, ["signor"], examples=2)
    assert not failures(report), json.dumps(failures(report), indent=2)
    resource = report["resources"][0]
    assert report["status"] == "passed"
    from scripts import validate_molecular_outputs as validator_module

    assert (
        report["validator_sha256"]
        == hashlib.sha256(Path(validator_module.__file__).read_bytes()).hexdigest()
    )
    assert resource["coverage"]["paired_occurrences"] == 1
    assert resource["coverage"]["modifications"] == 2
    assert any(
        c["name"] == "observed.variants" and c["status"] == "unavailable"
        for c in resource["checks"]
    )
    assert any(c["name"].endswith("product_closure") for c in resource["checks"])
    assert any(c["name"] == "api.paired_occurrences.1" for c in resource["checks"])
    mixed = next(
        example
        for example in resource["examples"]
        if example["name"] == "opposite_endpoint_constraints_1"
    )
    assert mixed["expected_relations"] == 0
    assert mixed["filters"]["isoform_identifiers"] == ["uniprot:P00533-3"]
    assert any(c["name"] == "api.gene_type_filter.protein" for c in resource["checks"])
    assert any(c["name"] == "api.gene_type_filter.gene" for c in resource["checks"])
    assert any(t["kind"] == "api" for t in resource["timings"])
    assert resource["tables"]["entity"]["bytes"] > 0


def test_validator_detects_dangling_product_reference(outputs):
    root, folder = outputs
    rows = nested_rows(folder / "relation.parquet")
    rows[0]["evidence"][0]["subject_molecular_form"]["protein_entity_key"] = "missing-product"
    rewrite_resource(folder, relations=rows)
    report = validate_outputs(root, ["signor/1"], examples=1)
    assert any(c["name"] == "molecular.product_closure_and_type" for c in failures(report))


def test_validator_reports_schema_error_without_unstructured_api_fallback(outputs):
    root, folder = outputs
    table = pq.read_table(folder / "entity.parquet").drop(["reference_entity_key"])
    pq.write_table(table, folder / "entity.parquet")
    report = validate_outputs(root, ["signor"], examples=1)
    assert failures(report)[0]["name"] == "schema.entity"
    assert any(
        c["name"] == "data_and_api_checks" and c["status"] == "unavailable"
        for c in report["resources"][0]["checks"]
    )


def test_validator_detects_form_specific_product_hash(outputs):
    root, folder = outputs
    rows = nested_rows(folder / "entity.parquet")
    entity(rows, "protein", "P04637")["entity_key"] = "form-enumerated-key"
    rewrite_resource(folder, entities=rows)
    report = validate_outputs(root, ["signor"], examples=1)
    assert any(c["name"] == "molecular.no_form_specific_product_keys" for c in failures(report))


def test_validator_checks_physical_arrow_type(outputs):
    root, folder = outputs
    path = folder / "entity.parquet"
    table = pq.read_table(path)
    table = table.set_column(
        table.schema.get_field_index("reference_entity_key"),
        "reference_entity_key",
        pa.array([7157] * table.num_rows, type=pa.int64()),
    )
    pq.write_table(table, path)
    report = validate_outputs(root, ["signor"], examples=1)
    error = next(c for c in failures(report) if c["name"] == "schema.entity")
    assert any(
        "reference_entity_key: actual int64; expected string" in message
        for message in error["errors"]
    )


def add_matching_relation(outputs):
    _, folder = outputs
    rows = nested_rows(folder / "relation.parquet")
    second = copy.deepcopy(rows[0])
    second["relation_key"] = "another-relation"
    rows.append(second)
    rewrite_resource(folder, relations=rows)


@pytest.mark.parametrize("damage", ["invented_empty", "duplicate_key"])
def test_export_audit_rejects_substituted_relation_key_multiset(outputs, monkeypatch, damage):
    add_matching_relation(outputs)
    original = Validator.request

    def corrupted(self, client, method, path, **kwargs):
        response = original(self, client, method, path, **kwargs)
        if path != "/export":
            return response
        table = pq.read_table(io.BytesIO(response.content))
        rows = table.to_pylist()
        if damage == "duplicate_key":
            rows = [rows[0]] * len(rows)
        else:
            for index, row in enumerate(rows):
                row.update(
                    relation_key=f"invented-{index}",
                    evidence=[],
                    evidence_count=0,
                    referenced_product_records=[],
                )
        buffer = io.BytesIO()
        pq.write_table(pa.Table.from_pylist(rows, schema=table.schema), buffer)
        return SimpleNamespace(content=buffer.getvalue())

    monkeypatch.setattr(Validator, "request", corrupted)
    report = validate_outputs(outputs[0], ["signor"], examples=1)
    assert any(c["name"].endswith(".relation_keys") for c in failures(report))
    assert not any(c["name"].endswith(".count") for c in failures(report))
    assert not any(c["name"].startswith("api.") for c in failures(report))
    if damage == "invented_empty":
        assert any(c["name"].endswith(".paired_occurrences") for c in failures(report))


@pytest.mark.parametrize("damage", ["empty_page", "duplicate_page"])
def test_api_audit_requires_populated_unique_search_pages(outputs, monkeypatch, damage):
    add_matching_relation(outputs)
    original = Validator.request

    def corrupted(self, client, method, path, **kwargs):
        response = original(self, client, method, path, **kwargs)
        if path != "/relations/search":
            return response
        data = response.json()
        data["relations"] = (
            [] if damage == "empty_page" else [data["relations"][0]] * len(data["relations"])
        )
        return SimpleNamespace(json=lambda: data)

    monkeypatch.setattr(Validator, "request", corrupted)
    report = validate_outputs(outputs[0], ["signor"], examples=1)
    assert any(
        c["name"].endswith(".page") and ".product_context." not in c["name"]
        for c in failures(report)
    )
    assert not any(c["name"].endswith(".count") for c in failures(report))
    assert not any(c["name"].startswith("export.") for c in failures(report))


@pytest.mark.parametrize("damage", ["empty_page", "duplicate_page", "untrimmed", "opposite_form"])
def test_product_context_audit_checks_page_and_paired_occurrences(outputs, monkeypatch, damage):
    add_matching_relation(outputs)
    raw_evidence = nested_rows(outputs[1] / "relation.parquet")[0]["evidence"]
    original = Validator.request

    def corrupted(self, client, method, path, **kwargs):
        response = original(self, client, method, path, **kwargs)
        if (
            not path.endswith("/molecular-context")
            or kwargs.get("params", {}).get("view") != "product"
        ):
            return response
        data = response.json()
        if damage == "empty_page":
            data["relations"] = []
        elif damage == "duplicate_page":
            data["relations"] = [data["relations"][0]] * len(data["relations"])
        else:
            for row in data["relations"]:
                if damage == "untrimmed":
                    row["evidence"] = copy.deepcopy(raw_evidence)
                else:
                    for ev in row["evidence"]:
                        if ev.get("object_molecular_form"):
                            ev["object_molecular_form"]["isoform_identifier"] = {
                                "ns": "uniprot",
                                "id": "P00533-99",
                            }
        return SimpleNamespace(json=lambda: data)

    monkeypatch.setattr(Validator, "request", corrupted)
    report = validate_outputs(outputs[0], ["signor"], examples=1)
    assert any(".product_context." in c["name"] for c in failures(report))
    assert not any(c["name"].endswith(".product_context") for c in failures(report))
    assert not any(c["name"].startswith("export.") for c in failures(report))


def test_standalone_context_audit_rejects_repeated_observation_with_same_form(outputs, monkeypatch):
    path = outputs[1] / "entity.parquet"
    rows = nested_rows(path)
    source = next(row for row in rows if row["evidence"])
    second = copy.deepcopy(source["evidence"][0])
    second["row_id"] = "another-standalone"
    source["evidence"].append(second)
    rewrite_resource(path, entities=rows)
    original = Validator.request

    def corrupted(self, client, method, path, **kwargs):
        response = original(self, client, method, path, **kwargs)
        if (
            not path.endswith("/molecular-context")
            or kwargs.get("params", {}).get("view") != "product"
        ):
            return response
        data = response.json()
        if len(data["standaloneEvidence"]) > 1:
            data["standaloneEvidence"][1] = copy.deepcopy(data["standaloneEvidence"][0])
        return SimpleNamespace(json=lambda: data)

    monkeypatch.setattr(Validator, "request", corrupted)
    report = validate_outputs(outputs[0], ["signor"], examples=1)
    assert any(c["name"].endswith(".standalone_forms") for c in failures(report))
    assert not any(c["name"].startswith("export.") for c in failures(report))


def test_api_page_proof_allows_valid_sort_ties(outputs, monkeypatch):
    add_matching_relation(outputs)
    original = Validator.request

    def reordered(self, client, method, path, **kwargs):
        response = original(self, client, method, path, **kwargs)
        if path == "/relations/search" or path.endswith("/molecular-context"):
            data = response.json()
            data["relations"].reverse()
            return SimpleNamespace(json=lambda: data)
        return response

    monkeypatch.setattr(Validator, "request", reordered)
    report = validate_outputs(outputs[0], ["signor"], examples=1)
    assert not failures(report), failures(report)


def plan_nodes(plan):
    for node in plan:
        yield node
        yield from plan_nodes(node.get("children", []))


def test_validator_base_expansions_and_example_queries_have_no_delimiter_joins(
    outputs, monkeypatch
):
    original_sql = Validator.sql
    checked = set()
    base_queries = {
        "molecular.product_references",
        "molecular.coverage",
        "molecular.paired_occurrences",
        "baseline.product_examples",
        "baseline.isoform_examples",
        "baseline.paired_examples",
    }

    def checked_sql(self, name, query, params=None):
        if name in base_queries or name == "molecular.reusable_uniprot_products":
            plan = json.loads(
                self.db.execute("EXPLAIN (FORMAT JSON) " + query, params or []).fetchone()[1]
            )
            delimiters = [node for node in plan_nodes(plan) if "DELIM_JOIN" in node["name"]]
            if name in base_queries:
                assert not delimiters, (name, delimiters)
            else:
                # The two short-list provenance EXISTS checks remain MARK
                # joins; evidence-array expansion must never add INNER ones.
                assert len(delimiters) <= 2
                assert all(node["extra_info"]["Join Type"] == "MARK" for node in delimiters)
            checked.add(name)
        return original_sql(self, name, query, params)

    monkeypatch.setattr(Validator, "sql", checked_sql)
    report = validate_outputs(outputs[0], ["signor"], examples=1)
    assert not failures(report), failures(report)
    assert checked == base_queries | {"molecular.reusable_uniprot_products"}


def structural_checks(folder):
    validator = Validator(folder, examples=0, max_export_relations=0)
    try:
        assert validator.schemas()
        validator.structural()
        return {check["name"]: check for check in validator.report["checks"]}
    finally:
        validator.db.close()


def reported_native_product(outputs, identifier="Q9Y6K9-2-PRO_000001"):
    """Use the writer's precise provenance format on every pointed occurrence."""
    _, folder = outputs
    path = folder / "entity.parquet"
    rows = nested_rows(path)
    product = entity(rows, "protein", "P04637")
    old_key = product["entity_key"]
    product["identifier"] = identifier
    product["entity_key"] = entity_key("protein", "uniprot", identifier)
    product["reference_entity_key"] = "uniprot:" + identifier
    product["gene_reference_keys"] = []

    def update(form, annotations, scope):
        if not form or form.get("protein_entity_key") != old_key:
            return
        form["protein_entity_key"] = product["entity_key"]
        form["isoform_identifier"] = {"ns": "uniprot", "id": identifier.split("-PRO_")[0]}
        form["sequence_identifiers"] = [{"ns": "uniprot", "id": identifier}]
        annotations.append(
            dict(
                term="omnipath:protein_mapping_status",
                value="reported",
                source="resolver",
                dataset="molecular_reference",
                scope=scope,
            )
        )

    for row in rows:
        for occurrence in row["evidence"] or []:
            update(occurrence["molecular_form"], occurrence["annotations"], None)
    rewrite_resource(path, entities=rows)
    path = folder / "relation.parquet"
    rows = nested_rows(path)
    for row in rows:
        for occurrence in row["evidence"]:
            for side in ("subject", "object"):
                update(occurrence[side + "_molecular_form"], occurrence["annotations"], side)
    rewrite_resource(path, relations=rows)
    return folder


@pytest.mark.parametrize("identifier", ["Q9Y6K9-2", "Q9Y6K9-PRO_000001", "Q9Y6K9-2-PRO_000001"])
def test_reported_native_uniprot_suffixes_require_every_occurrence_provenance(outputs, identifier):
    checks = structural_checks(reported_native_product(outputs, identifier))
    assert checks["molecular.reusable_uniprot_products"]["status"] == "passed"
    assert checks["molecular.no_form_specific_product_keys"]["status"] == "passed"
    assert all(check["status"] != "failed" for check in checks.values())


@pytest.mark.parametrize("owner", ["relation", "standalone"])
@pytest.mark.parametrize(
    "damage",
    [
        "missing",
        "wrong_scope",
        "wrong_source",
        "wrong_dataset",
        "resolved",
        "missing_assertion",
        "wrong_assertion",
    ],
)
def test_native_uniprot_suffix_cannot_borrow_other_occurrence_or_endpoint_provenance(
    outputs, owner, damage
):
    folder = reported_native_product(outputs)
    path = folder / ("relation.parquet" if owner == "relation" else "entity.parquet")
    rows = nested_rows(path)
    occurrence = (
        rows[0]["evidence"][0]
        if owner == "relation"
        else entity(rows, "protein", "7157", "entrez")["evidence"][0]
    )
    annotations = occurrence["annotations"]
    status = next(
        annotation
        for annotation in annotations
        if annotation["term"] == "omnipath:protein_mapping_status"
    )
    form = occurrence["subject_molecular_form" if owner == "relation" else "molecular_form"]
    if damage == "missing":
        annotations.remove(status)
    elif damage == "wrong_scope":
        if owner == "relation":
            status["scope"] = "object"
        else:
            # Standalone annotations have no scope field in their schema.
            # Provenance on another entity occurrence cannot authorize this one.
            annotations.remove(status)
            entity(rows, "protein", "1956", "entrez")["evidence"] = [
                dict(source="signor", row_id="other", annotations=[status], molecular_form=None)
            ]
    elif damage == "wrong_source":
        status["source"] = "raw"
    elif damage == "wrong_dataset":
        status["dataset"] = "gene_reference"
    elif damage == "resolved":
        status["value"] = "resolved"
    elif damage == "missing_assertion":
        form["sequence_identifiers"] = []
        form["isoform_identifier"] = None
    elif damage == "wrong_assertion":
        # An enclosing isoform is insufficient proof of an exact chain assertion.
        form["sequence_identifiers"] = [{"ns": "uniprot", "id": "Q9Y6K9-2"}]
    rewrite_resource(path, **({"relations": rows} if owner == "relation" else {"entities": rows}))
    check = structural_checks(folder)["molecular.reusable_uniprot_products"]
    assert check["status"] == "failed"
    assert check["invalid_count"] == 1


@pytest.mark.parametrize(
    "identifier", ["NOT_A_VALID_ACCESSION-2", "Q9Y6K9-BAD", "Q9Y6K9-2-PRO_BAD"]
)
def test_reported_provenance_does_not_waive_invalid_uniprot_identifier(outputs, identifier):
    check = structural_checks(reported_native_product(outputs, identifier))[
        "molecular.reusable_uniprot_products"
    ]
    assert check["status"] == "failed"
    assert check["invalid_count"] == 3


def test_both_product_fields_are_checked_without_extending_protein_provenance_to_transcript(
    outputs,
):
    folder = reported_native_product(outputs)
    entity_path = folder / "entity.parquet"
    entities = nested_rows(entity_path)
    transcript_key = entity_key("transcript", "uniprot", "Q9Y6K9-2")
    entities.append(
        dict(
            entity_key=transcript_key,
            entity_type="transcript",
            namespace="uniprot",
            identifier="Q9Y6K9-2",
            reference_entity_key="uniprot:Q9Y6K9-2",
            gene_reference_keys=[],
            evidence=[],
        )
    )
    rewrite_resource(entity_path, entities=entities)
    relation_path = folder / "relation.parquet"
    relations = nested_rows(relation_path)
    relations[0]["evidence"][0]["subject_molecular_form"]["transcript_entity_key"] = transcript_key
    rewrite_resource(relation_path, relations=relations)
    check = structural_checks(folder)["molecular.reusable_uniprot_products"]
    assert check["status"] == "failed"
    assert check["invalid_count"] == 1


@pytest.mark.parametrize("identifier", ["Q9Y6K9-2", "Q9Y6K9-PRO_000001", "Q9Y6K9-2-PRO_000001"])
def test_actual_writer_reported_native_uniprot_provenance_passes_validator(tmp_path, identifier):
    from omnipath_resolver.resolver import EntityResolver
    from omnipath_build.silver import SilverExtractor
    from omnipath_build.writer import ParquetWriter

    folder = tmp_path / "resources" / "fixture" / "1"
    item = {"type": "protein", "identifiers": [{"type": "uniprot", "value": identifier}]}
    extractor = SilverExtractor("fixture", "reported-native")
    extractor.process_record(item, item, "standalone", 0)
    relation = {
        "subject": item,
        "predicate": "affects",
        "object": {"type": "gene", "identifiers": [{"type": "entrez", "value": "55"}]},
    }
    extractor.process_record(relation, relation, "relation", 1)
    resolver = EntityResolver(library_dir=tmp_path / "absent")
    writer = ParquetWriter(folder)
    try:
        writer.append_observations(extractor, resolver)
        writer.close()
    finally:
        resolver.close()
    checks = structural_checks(folder)
    assert checks["molecular.reusable_uniprot_products"]["status"] == "passed"
    assert all(check["status"] != "failed" for check in checks.values())
