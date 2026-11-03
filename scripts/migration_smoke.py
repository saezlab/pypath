#!/usr/bin/env python3
"""Build and inspect a bounded, entirely offline SIGNOR migration fixture.

This uses the real inputs_v2 SIGNOR parser and RelationBuilder, native entity
resolver, Parquet writer and immutable resource publication. The small reference
and source rows are synthetic fixtures; no biological assertions are implied.
"""

from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
import csv
import hashlib
from io import StringIO
import json
import os
from pathlib import Path
import shutil
from types import SimpleNamespace
from unittest.mock import patch


FIXTURE = Path(__file__).parent / "fixtures" / "signor.json"
SOURCE_ROWS = 24  # More than the permitted cap, so every run exercises truncation.
RECORD_CAP = 20


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def record_limit(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("max-records must be an integer from 1 to 20") from exc
    if not 1 <= number <= RECORD_CAP:
        raise argparse.ArgumentTypeError("max-records must be between 1 and 20")
    return number


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def source_fixture(path: Path, fixture: dict) -> list[dict]:
    rows = []
    for index in range(SOURCE_ROWS):
        row = dict(fixture["rows"][index % len(fixture["rows"])])
        row["Interaction identifier(s)"] = f"signor:SIGNOR-SMOKE-{index + 1}"
        row["Publication Identifier(s)"] = f"pubmed:{123456 + index}"
        row["Interaction annotation(s)"] += f" (fixture record {index + 1})"
        rows.append(row)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def identity_reference(path: Path, fixture: dict) -> Path:
    """Build an identity library for two synthetic proteins with the production builders.

    This intentionally avoids source hub downloads: a one-hub UniProt export is written
    directly, then indexed, decided and stored by the same code as a full reference build.
    Returns the identity directory the resource build resolves against.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    from omnipath_build.hubs.schema import HUB_SCHEMA, hub_row
    from omnipath_build.identity import (
        build_hub_index,
        build_hub_kv,
        build_identity,
        build_identity_kv,
    )

    if path.exists():
        require(path.is_dir() and not path.is_symlink(), f"Unsafe fixture reference path: {path}")
        marker = path / "fixture.json"
        require(marker.is_file(), f"Refusing to replace a non-fixture reference: {path}")
        shutil.rmtree(path)
    hubs = path / "hubs"
    hubs.mkdir(parents=True)
    (path / "fixture.json").write_text(json.dumps({"description": fixture["description"]}))
    rows = []
    for item in fixture["reference"]:
        accession = item["identifier"]
        rows.append(hub_row("name", item["label"], accession, "9606", "uniprot"))
        for namespace, identifier in item["identifiers"]:
            if namespace == "uniprot" and identifier != accession:
                namespace = "uniprot-sec"
            rows.append(hub_row(namespace, identifier, accession, "9606", "uniprot"))
    require(len(rows) <= RECORD_CAP, "Reference hub fixture exceeds the cap")
    pq.write_table(pa.Table.from_pylist(rows, schema=HUB_SCHEMA), hubs / "uniprot.parquet")
    options = dict(memory="256MB", threads=1, min_free_gib=0)
    build_hub_index("uniprot", hubs, path / "hub-index", **options)
    build_hub_kv(
        "uniprot", path / "hub-index", memory="256MB", threads=1, workers=1, min_free_gib=0
    )
    identity = build_identity(path / "hub-index", path / "identity", **options)
    library = path / "identity" / identity["fingerprint"]
    build_identity_kv(library, min_free_gib=0)
    return library


class LocalDownload:
    """Only a fixture opener is substituted; parsing and mapping stay real."""

    def __init__(self, path: Path):
        self.path = path
        self.url = path.as_uri()

    def open(self, **_kwargs):
        return SimpleNamespace(result=StringIO(self.path.read_text(encoding="utf-8")))


def blocked_download(*_args, **_kwargs):
    raise RuntimeError("Network downloads are disabled for the migration smoke build")


def inspect_resource(directory: Path) -> dict:
    """Check every manifest table's schema and checksum, then query a small DuckDB sample."""
    import duckdb
    import pyarrow.parquet as pq
    from omnipath_core.schema import PUBLISHED_TABLES, SERVING_TABLES

    manifest = json.loads((directory / "build_manifest.json").read_text())
    tables = {**PUBLISHED_TABLES, **SERVING_TABLES}
    require(
        set(manifest["files"]) == {f"{table}.parquet" for table in tables},
        "The manifest does not list exactly the published and serving tables",
    )
    files = {}
    for table, schema in tables.items():
        name = f"{table}.parquet"
        path = directory / name
        parquet = pq.ParquetFile(path)
        expected = manifest["files"][name]
        require(parquet.schema_arrow.equals(schema), f"Unexpected schema in {name}")
        require(parquet.metadata.num_rows == expected["rows"], f"Row count differs in {name}")
        require(path.stat().st_size == expected["size_bytes"], f"Size differs in {name}")
        require(digest(path) == expected["sha256"], f"Checksum differs in {name}")
        files[name] = dict(expected, schema_verified=True, checksum_verified=True)
    with duckdb.connect(config={"threads": "1", "memory_limit": "256MB"}) as connection:
        relations = connection.execute(
            "SELECT s.identifier AS subject,r.subject_label,r.predicate,"
            "o.identifier AS object,r.object_label,r.sign,r.evidence_count "
            "FROM read_parquet(?) r "
            "JOIN read_parquet(?) s ON s.entity_key=r.subject_entity_key "
            "JOIN read_parquet(?) o ON o.entity_key=r.object_entity_key "
            "ORDER BY r.subject_label,r.object_label LIMIT 20",
            [str(directory / "relation.parquet"), *[str(directory / "entity.parquet")] * 2],
        )
        sample = [
            dict(zip((column[0] for column in relations.description), row))
            for row in relations.fetchall()
        ]
    return {
        "resource": manifest["resource"],
        "version": manifest["version"],
        "directory": str(directory),
        "max_records": manifest["max_records"],
        "files": files,
        "relations": sample,
    }


def verify_fixture(result: dict, rows: list[dict]) -> dict:
    import pyarrow.parquet as pq

    def table(name: str) -> list[dict]:
        return pq.read_table(result["files"][name]).to_pylist()

    entities = {entity["entity_id"]: entity for entity in table("entity")}
    identifiers = table("entity_identifier")
    relations = {relation["relation_id"]: relation for relation in table("relation")}
    evidence = table("relation_evidence")
    relation_annotations = table("relation_annotation")
    payloads = table("evidence_payloads")
    canonical = {entity["identifier"]: entity for entity in entities.values()}
    require(
        canonical.keys() == {"P04637", "P0DP23"}
        and all(entity["namespace"] == "uniprot" for entity in entities.values()),
        "The native resolver did not produce both canonical protein identifiers",
    )
    require(all(entity["taxon"] == "9606" for entity in entities.values()), "Entity taxa were lost")
    require(
        all(
            any(
                identifier["entity_id"] == entity_id
                and identifier["is_canonical"]
                and (identifier["ns"], identifier["id"]) == ("uniprot", entity["identifier"])
                for identifier in identifiers
            )
            for entity_id, entity in entities.items()
        ),
        "An entity lacks its canonical identifier row",
    )
    require(result["resolution_stats"]["unresolved_entities"] == 0, "An entity was unresolved")
    require(result["resolution_stats"]["resolved_entities"] >= 2, "Resolution was not exercised")
    require(
        any(
            identifier["id"] == "Q15086"
            and not identifier["is_canonical"]
            and identifier["entity_id"] == canonical["P04637"]["entity_id"]
            for identifier in identifiers
        ),
        "The observed secondary accession was lost",
    )
    entity_keys = {entity["entity_key"] for entity in entities.values()}
    require(
        all(
            {relation["subject_entity_key"], relation["object_entity_key"]} <= entity_keys
            for relation in relations.values()
        ),
        "A relation endpoint is missing from the entity table",
    )
    require(len(payloads) == len(rows), "Not every consumed source record retained a payload")
    require(
        Counter(
            json.dumps(json.loads(payload["payload_json"]), sort_keys=True) for payload in payloads
        )
        == Counter(json.dumps(row, sort_keys=True) for row in rows),
        "Raw source payloads changed",
    )
    relation_keys = {relation["relation_key"] for relation in relations.values()}
    require(
        all(payload["relation_key"] in relation_keys for payload in payloads),
        "A raw payload is not linked to a published relation",
    )
    require(len(evidence) == len(rows), "Source evidence was dropped during consolidation")
    require(
        {item["row_id"] for item in evidence} == {payload["row_id"] for payload in payloads},
        "Relation evidence and raw payloads disagree on the consumed source rows",
    )
    per_relation = Counter(item["relation_id"] for item in evidence)
    require(
        per_relation.keys() == relations.keys()
        and all(
            relation["evidence_count"] == per_relation[relation_id]
            for relation_id, relation in relations.items()
        ),
        "Consolidated evidence counts do not agree",
    )
    require(
        all(relation["predicate"] == "affects" for relation in relations.values()),
        "The SIGNOR Biolink predicate was lost",
    )
    publications = {
        annotation["value"]
        for item in evidence
        for annotation in item["annotations"]
        if annotation["term"] == "publications"
    }
    require(
        publications == {"PMID:" + row["Publication Identifier(s)"].split(":")[1] for row in rows},
        "Publication evidence was lost",
    )
    require(
        all(relation["object_direction_qualifier"] for relation in relations.values()),
        "Biolink relation qualifiers were lost",
    )
    activity = {
        annotation["relation_id"]
        for annotation in relation_annotations
        if annotation["term"] == "object_aspect_qualifier" and annotation["value"] == "activity"
    }
    require(
        all(
            relation["object_aspect_qualifier"] == ["activity"] and relation_id in activity
            for relation_id, relation in relations.items()
        ),
        "The Biolink activity qualifier was lost",
    )
    return {
        "canonical_identifiers": sorted(canonical),
        "secondary_accession_resolved": "Q15086 -> P04637",
        "payloads_and_evidence_preserved": len(rows),
        "resolution": result["resolution_stats"],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("data/migration-smoke"))
    parser.add_argument("--version", default="0.1.0", help="New immutable numeric resource version")
    parser.add_argument("--max-records", type=record_limit, default=RECORD_CAP)
    parser.add_argument(
        "--inspect",
        type=Path,
        metavar="RESOURCE_VERSION_DIR",
        help="Verify and inspect an existing version without building",
    )
    args = parser.parse_args(argv)
    if args.inspect:
        print(json.dumps(inspect_resource(args.inspect.resolve()), indent=2))
        return 0

    from omnipath_build.pipeline import build_resource
    from omnipath_build.versioning import validate_version
    from pypath.inputs_v2 import signor

    validate_version(args.version)
    root = args.output_dir.resolve()
    target = root / "resources" / "signor" / args.version
    if target.exists():
        raise FileExistsError(
            f"Resource signor/{args.version} already exists; choose a new version"
        )
    fixture = json.loads(FIXTURE.read_text())
    fixture_dir = root / "fixtures" / args.version
    fixture_dir.mkdir(parents=True, exist_ok=True)
    rows = source_fixture(fixture_dir / "signor.tsv", fixture)
    reference = identity_reference(fixture_dir / "reference", fixture)
    dataset = signor.resource.datasets()["interactions"]
    raw_parser = dataset._raw_parser
    consumed = 0

    def bounded_parser(opener, **kwargs):
        nonlocal consumed
        for row in raw_parser(opener, **kwargs):
            consumed += 1
            require(consumed <= args.max_records, "The parser consumed a record beyond the cap")
            yield row

    # The isolated worker imports the real mapper from pypath.inputs_v2.signor.
    # Only this parent process opens the synthetic input file.
    os.environ["OMNIPATH_BUILD_DUCKDB_THREADS"] = "1"
    with ExitStack() as stack:
        stack.enter_context(
            patch.object(dataset, "download", LocalDownload(fixture_dir / "signor.tsv"))
        )
        stack.enter_context(patch.object(dataset, "_raw_parser", bounded_parser))
        stack.enter_context(patch("pypath.inputs_v2.base.download_and_open", blocked_download))
        result = build_resource(
            "signor",
            version=args.version,
            output_dir=root,
            datasets=["interactions"],
            max_records=args.max_records,
            batch_workers=1,
            batch_size=20,
            max_batch_records=20,
            max_batch_relations=20,
            duckdb_memory_limit="256MB",
            min_free_disk_bytes=0,
            library_dir=reference,
            cache_dir=root / "pypath-data",
            progress=False,
        )
    require(consumed == args.max_records, "The source cap was not exercised")
    checks = verify_fixture(result, rows[:consumed])
    report = inspect_resource(target)
    before = digest(target / "build_manifest.json")
    try:
        build_resource(
            "signor", version=args.version, output_dir=root, max_records=args.max_records
        )
    except FileExistsError:
        pass
    else:
        raise RuntimeError("A duplicate resource version was accepted")
    require(
        digest(target / "build_manifest.json") == before, "Duplicate build changed the manifest"
    )
    inspect_resource(target)  # Published files must also remain unchanged after rejection.
    report.update(
        fixture=True,
        source_rows_available=len(rows),
        source_rows_consumed=consumed,
        duplicate_version_rejected=True,
        checks=checks,
    )
    report_path = root / f"inspection-{args.version}.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
