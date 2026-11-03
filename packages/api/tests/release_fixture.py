"""Tiny resolved releases for product validation, without source builds."""

import hashlib
import json

import pyarrow.parquet as pq

from omnipath_core import SERVING_SCHEMA_VERSION
from omnipath_core.keys import entity_key, relation_key
from omnipath_core.fixtures import write_resource as write_tables


def annotation(
    term, value=None, *, source="fixture", dataset="fixture", scope="relation", quantity=None
):
    return dict(
        term=term, value=value, quantity=quantity, source=source, dataset=dataset, scope=scope
    )


def entity(
    identifier,
    entity_type="small_molecule",
    namespace="chebi",
    *,
    taxon="",
    label=None,
    aliases=(),
    annotations=(),
):
    return dict(
        entity_key=entity_key(entity_type, namespace, identifier),
        entity_type=entity_type,
        namespace=namespace,
        identifier=identifier,
        taxon=taxon,
        label=label or identifier,
        has_hierarchy=False,
        parent_count=0,
        child_count=0,
        identifiers=[
            dict(ns=namespace, id=identifier, is_canonical=True, source="canonical"),
            *[dict(ns=ns, id=value, is_canonical=False, source="raw") for ns, value in aliases],
        ],
        annotations=[
            {key: value for key, value in item.items() if key != "scope"} for item in annotations
        ],
    )


def relation(
    subject,
    predicate,
    obj,
    *,
    source="fixture",
    dataset="fixture",
    row_id="1",
    annotations=(),
    statement_kind="relation",
    upstream_id="",
    sign=0,
):
    key = relation_key(
        subject["entity_key"],
        predicate,
        obj["entity_key"],
        annotations,
        statement_kind=statement_kind,
    )
    return dict(
        relation_key=key,
        statement_kind=statement_kind,
        subject_entity_key=subject["entity_key"],
        subject_label=subject["label"],
        subject_type=subject["entity_type"],
        predicate=predicate,
        object_entity_key=obj["entity_key"],
        object_label=obj["label"],
        object_type=obj["entity_type"],
        taxon="",
        is_directed=predicate != "interacts_with",
        sign=sign,
        category="association",
        interaction_class="directed",
        sources=[source],
        evidence_count=1,
        evidence=[
            dict(
                source=source,
                dataset=dataset,
                row_id=row_id,
                upstream_id=upstream_id,
                annotations=list(annotations),
            )
        ],
        annotations=list(annotations),
    )


def payload(owner, raw, *, source="fixture", row_id="1"):
    return dict(
        relation_key=owner.get("relation_key"),
        entity_key=owner.get("entity_key"),
        source=source,
        row_id=row_id,
        payload_json=json.dumps(raw),
    )


def write_resource(root, source, entities, relations=(), payloads=(), version="1.0.0"):
    # Record caps count original source rows; one source record may yield several spokes.
    assert len({row["row_id"] for rel in relations for row in rel["evidence"]}) <= 20
    directory = root / "resources" / source / version
    directory.mkdir(parents=True)
    write_tables(directory, entities, relations, payloads)
    files = {}
    for path in sorted(directory.glob("*.parquet")):
        files[path.name] = dict(
            rows=pq.read_metadata(path).num_rows,
            size_bytes=path.stat().st_size,
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        )
    (directory / "build_manifest.json").write_text(
        json.dumps(
            dict(
                schema_version=1,
                serving_schema_version=SERVING_SCHEMA_VERSION,
                resource=source,
                version=version,
                max_records=None,
                files=files,
            )
        )
    )
    return directory


def write_release(root, sources, version="2026.09"):
    path = root / "release.json"
    selections = (
        dict(sources) if isinstance(sources, dict) else {source: "1.0.0" for source in sources}
    )
    path.write_text(json.dumps(dict(schema_version=1, version=version, resources=selections)))
    return path
