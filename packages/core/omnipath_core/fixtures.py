"""Test fixtures: write a resource version's tables from readable nested rows.

Fixtures describe entities with their identifiers, annotations and evidence, and
relations with their annotations and evidence, as the API returns them; this writes
them as the normalized tables the build publishes (``omnipath_core.schema``).
"""

from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from omnipath_core.measurements import QUANTITY_STRUCT
from omnipath_core.schema import PUBLISHED_TABLES, RELATION_QUALIFIERS, SERVING_TABLES
from omnipath_core.serving_tables import write_serving_tables

_QUANTITY = [field.name for field in QUANTITY_STRUCT]


def _flat_quantity(annotation):
    quantity = annotation.get("quantity") or {}
    return {f"quantity_{name}": quantity.get(name) for name in _QUANTITY}


def _connectivity(entity):
    reference = entity.get("reference_entity_key") or ""
    return reference[9:23] if reference.startswith("inchikey:") else None


def write_resource(folder, entities=(), relations=(), payloads=None):
    """Write every table of one resource version into ``folder``; returns the folder."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    entities = sorted(entities, key=lambda e: e["entity_key"])
    relations = sorted(relations, key=lambda r: r["relation_key"])
    relation_count = {}
    for relation in relations:
        for side in ("subject", "object"):
            key = relation.get(f"{side}_entity_key")
            relation_count.setdefault(key, set()).add(relation["relation_key"])
    tables = {name: [] for name in (*PUBLISHED_TABLES, *SERVING_TABLES)}
    for entity_id, entity in enumerate(entities):
        identifiers = entity.get("identifiers") or []
        annotations = entity.get("annotations") or []
        evidence = entity.get("evidence") or []
        tables["entity"].append(
            dict(
                entity_id=entity_id,
                entity_key=entity["entity_key"],
                entity_type=entity.get("entity_type"),
                namespace=entity.get("namespace"),
                identifier=entity.get("identifier"),
                taxon=entity.get("taxon") or "",
                label=entity.get("label"),
                has_hierarchy=bool(entity.get("has_hierarchy")),
                parent_count=entity.get("parent_count") or 0,
                child_count=entity.get("child_count") or 0,
                reference_entity_key=entity.get("reference_entity_key"),
                gene_reference_keys=entity.get("gene_reference_keys"),
                group_connectivity=_connectivity(entity),
                identifier_count=len(identifiers),
                annotation_count=len(annotations),
                evidence_count=len(evidence),
                relation_count=len(relation_count.get(entity["entity_key"], ())),
            )
        )
        for ordinal, item in enumerate(identifiers):
            tables["entity_identifier"].append(
                dict(
                    entity_id=entity_id,
                    ordinal=ordinal,
                    ns=item.get("ns"),
                    id=item.get("id"),
                    is_canonical=bool(item.get("is_canonical")),
                    source=item.get("source") or "",
                )
            )
        for ordinal, item in enumerate(annotations):
            tables["entity_annotation"].append(
                dict(
                    entity_id=entity_id,
                    ordinal=ordinal,
                    term=item.get("term"),
                    value=item.get("value"),
                    **_flat_quantity(item),
                    source=item.get("source"),
                    dataset=item.get("dataset"),
                )
            )
        for ordinal, item in enumerate(evidence):
            tables["entity_evidence"].append(
                dict(
                    item,
                    entity_id=entity_id,
                    ordinal=ordinal,
                    annotations=item.get("annotations") or [],
                )
            )
        reference = entity.get("reference_entity_key")
        groups = [("reference", reference)] if reference else []
        groups += [("gene", key) for key in entity.get("gene_reference_keys") or []]
        if _connectivity(entity):
            groups.append(("connectivity", _connectivity(entity)))
        tables["entity_group"] += [
            dict(group_key=key, kind=kind, entity_id=entity_id) for kind, key in groups
        ]
        terms = {("label", entity.get("label")), ("identifier", entity.get("identifier"))}
        terms |= {("identifier", item.get("id")) for item in identifiers}
        tables["entity_term"] += [
            dict(term=value.lower(), kind=kind, entity_id=entity_id)
            for kind, value in terms
            if value
        ]
    for relation_id, relation in enumerate(relations):
        annotations = relation.get("annotations") or []
        evidence = relation.get("evidence") or []
        qualifiers = {
            name: sorted(
                {
                    a.get("value")
                    for a in annotations
                    if a.get("term") == name and a.get("scope", "relation") == "relation"
                }
            )
            for name in RELATION_QUALIFIERS
        }
        tables["relation"].append(
            dict(
                {name: relation.get(name) for name in PUBLISHED_TABLES["relation"].names},
                relation_id=relation_id,
                statement_kind=relation.get("statement_kind") or "relation",
                taxon=relation.get("taxon") or "",
                is_directed=bool(relation.get("is_directed")),
                sign=relation.get("sign") or 0,
                sources=relation.get("sources") or [],
                evidence_count=relation.get("evidence_count", len(evidence)),
                annotation_count=len(annotations),
                **qualifiers,
            )
        )
        for ordinal, item in enumerate(annotations):
            tables["relation_annotation"].append(
                dict(
                    relation_id=relation_id,
                    ordinal=ordinal,
                    term=item.get("term"),
                    value=item.get("value"),
                    **_flat_quantity(item),
                    source=item.get("source"),
                    dataset=item.get("dataset"),
                    scope=item.get("scope", "relation"),
                )
            )
        for ordinal, item in enumerate(evidence):
            tables["relation_evidence"].append(
                dict(
                    item,
                    relation_id=relation_id,
                    ordinal=ordinal,
                    annotations=item.get("annotations") or [],
                )
            )
        for side in ("subject", "object"):
            for kind, key in (
                ("entity", relation.get(f"{side}_entity_key")),
                ("reference", relation.get(f"{side}_reference_entity_key")),
            ):
                if key:
                    tables["relation_endpoint"].append(
                        dict(
                            key=key,
                            key_kind=kind,
                            side=side,
                            relation_id=relation_id,
                            predicate=relation.get("predicate"),
                            category=relation.get("category"),
                        )
                    )
    tables["evidence_payloads"] = list(payloads or [])
    order = {
        "relation_endpoint": ("key", "relation_id"),
        "entity_group": ("group_key", "entity_id"),
        "entity_term": ("term", "entity_id"),
    }
    for name, schema in {**PUBLISHED_TABLES, **SERVING_TABLES}.items():
        if name == "relation_facet":
            continue
        rows = tables[name]
        if name in order:
            rows = sorted(rows, key=lambda r: tuple(r[c] for c in order[name]))
        pq.write_table(pa.Table.from_pylist(rows, schema=schema), folder / f"{name}.parquet")
    # The relation combinations are aggregates: written by the build's own query.
    import duckdb

    with tempfile.TemporaryDirectory() as work, duckdb.connect() as db:
        write_serving_tables(db, folder, work, tables=("relation_facet",))
    return folder


def write_manifest(folder, resource, version, **extra):
    """A valid build manifest for the tables written in ``folder``; returns it."""
    import hashlib

    from omnipath_core.versioning import BUILD_SCHEMA_VERSION, SERVING_SCHEMA_VERSION

    folder = Path(folder)
    files = {
        path.name: dict(
            rows=pq.read_metadata(path).num_rows,
            size_bytes=path.stat().st_size,
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        )
        for path in sorted(folder.glob("*.parquet"))
    }
    manifest = dict(
        schema_version=BUILD_SCHEMA_VERSION,
        serving_schema_version=SERVING_SCHEMA_VERSION,
        resource=resource,
        version=version,
        files=files,
        **extra,
    )
    (folder / "build_manifest.json").write_text(json.dumps(manifest))
    return manifest


def _folder(path):
    path = Path(path)
    return path.parent if path.suffix == ".parquet" else path


def read_resource(folder):
    """The nested (entities, relations) a resource version was written from; ``folder``
    may also be one of its table files."""
    folder = _folder(folder)

    def rows(name):
        return pq.read_table(folder / f"{name}.parquet").to_pylist()

    def children(name, parent):
        found = {}
        for row in sorted(rows(name), key=lambda r: (r[parent], r["ordinal"])):
            owner = row.pop(parent)
            row.pop("ordinal")
            if "quantity_" + _QUANTITY[0] in row:
                quantity = {q: row.pop("quantity_" + q) for q in _QUANTITY}
                row["quantity"] = (
                    quantity if any(v is not None for v in quantity.values()) else None
                )
            found.setdefault(owner, []).append(row)
        return found

    nested = {
        name: children(name, "entity_id")
        for name in ("entity_identifier", "entity_annotation", "entity_evidence")
    }
    entities = []
    for row in rows("entity"):
        entity_id = row.pop("entity_id")
        for derived in (
            "group_connectivity",
            "identifier_count",
            "annotation_count",
            "evidence_count",
            "relation_count",
        ):
            row.pop(derived)
        row["identifiers"] = nested["entity_identifier"].get(entity_id, [])
        row["annotations"] = nested["entity_annotation"].get(entity_id, [])
        row["evidence"] = nested["entity_evidence"].get(entity_id, [])
        entities.append(row)
    annotations = children("relation_annotation", "relation_id")
    evidence = children("relation_evidence", "relation_id")
    relations = []
    for row in rows("relation"):
        relation_id = row.pop("relation_id")
        for derived in ("annotation_count", *RELATION_QUALIFIERS):
            row.pop(derived)
        row["annotations"] = annotations.get(relation_id, [])
        row["evidence"] = evidence.get(relation_id, [])
        relations.append(row)
    return entities, relations


def rewrite_resource(folder, entities=None, relations=None):
    """Write a resource version again with new entities and/or relations (a new one
    when nothing is written yet)."""
    folder = _folder(folder)
    old_entities, old_relations, payloads = [], [], []
    if (folder / "entity.parquet").exists():
        old_entities, old_relations = read_resource(folder)
        payloads = pq.read_table(folder / "evidence_payloads.parquet").to_pylist()
    write_resource(
        folder,
        old_entities if entities is None else entities,
        old_relations if relations is None else relations,
        payloads,
    )


def nested_rows(path):
    """The nested entities (``entity.parquet``) or relations (``relation.parquet``)."""
    entities, relations = read_resource(path)
    return relations if Path(path).name.startswith("relation") else entities


def fixture_key(value):
    """A deterministic 64-character key for a readable fixture name."""
    return hashlib.sha256(value.encode()).hexdigest()


def write_dataset(root):
    """Small deterministic serving dataset: three resources under ``root/resources``."""
    entities = [
        dict(
            entity_key=fixture_key(ident),
            identifier=ident,
            label=label,
            namespace=ns,
            entity_type=kind,
            taxon="9606",
            has_hierarchy=kind == "ontology_class",
            parent_count=int(kind == "ontology_class"),
            child_count=0,
            identifiers=[dict(ns="uniprot", id="P00813", source="signor")]
            if label == "ADA"
            else [],
            annotations=[
                dict(
                    term="description", value="Example annotation", source="signor", dataset="test"
                )
            ],
        )
        for ident, label, ns, kind in [
            ("P04637", "TP53", "uniprot", "protein"),
            ("P38398", "BRCA1", "uniprot", "protein"),
            ("100", "ADA", "entrez", "protein"),
            ("KW-0001", "KW-0001", "keyword", "ontology_class"),
        ]
    ]
    relations = [
        dict(
            relation_key=fixture_key("relation"),
            subject_entity_key=fixture_key("P04637"),
            subject_label="TP53",
            object_entity_key=fixture_key("P38398"),
            object_label="BRCA1",
            subject_type="protein",
            object_type="protein",
            predicate="affects",
            category="interaction",
            interaction_class="protein_protein",
            taxon="9606",
            is_directed=True,
            sign=1,
            sources=["signor"],
            evidence_count=1,
            evidence=[],
            annotations=[],
        )
    ]
    for source in ["signor", "uniprot", "chebi"]:
        write_resource(
            Path(root) / "resources" / source / "1",
            entities if source != "chebi" else [],
            relations if source == "signor" else [],
            [
                dict(
                    relation_key=fixture_key("relation"),
                    source="signor",
                    row_id="1",
                    payload_json="{}",
                )
            ],
        )
    return Path(root)


def write_resource_version(root, resource, version, *, n_relations=1):
    """One protein and its self-relations as ``root/resources/<resource>/<version>``."""
    root = Path(root)
    target = root / "resources" / resource / version
    entity = dict(entity_key=fixture_key("P1"), identifier="P1", label="Prot", namespace="uniprot")
    relations = [
        dict(
            relation_key=fixture_key(f"r{i}"),
            subject_entity_key=entity["entity_key"],
            object_entity_key=entity["entity_key"],
            predicate="associated_with",
        )
        for i in range(n_relations)
    ]
    write_resource(target, [entity], relations)
    (target / "resolution_stats.json").write_text("{}\n", encoding="utf-8")
    return target
