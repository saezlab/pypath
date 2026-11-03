"""Replay up to 20 published raw records per resource through current inputs_v2.

Use a private output directory and an existing immutable resolution reference.
This checks the native mapper/resolver/writer; it never downloads source data.
Input JSON files contain dataset, original_row_id, original_payload_sha256 and
payload_json, selected from the original published evidence_payloads.parquet.
These capped artifacts are validation samples, not production replacements.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
from time import perf_counter
from unittest.mock import patch


CAP = 20
SOURCES = ("rhea", "recon3d", "metatlas", "kegg", "reactome", "macdb", "connectomedb2025")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def blocked_download(*_args, **_kwargs):
    raise RuntimeError("Source downloads are disabled for the annotation migration sample")


def inspect_resource(result):
    import pyarrow.parquet as pq
    from omnipath_core.source_attributes import (
        PARTICIPANT_ROLE,
        SOURCE_RECORD_REFERENCE,
        SOURCE_RECORD_SHA256_PREFIX,
        SOURCE_RECORD_TYPE,
        TRAIT_TYPE,
    )

    entities = pq.read_table(result["entities_path"]).to_pylist()
    relations = pq.read_table(result["relations_path"]).to_pylist()
    payloads = pq.read_table(result["payloads_path"]).to_pylist()
    payload_hashes = {
        (p["relation_key"], p["row_id"]): hashlib.sha256(p["payload_json"].encode()).hexdigest()
        for p in payloads
        if p["relation_key"]
    }
    terms = Counter()
    reaction_evidence = 0
    for relation in relations:
        for ev in relation["evidence"]:
            annotations = ev["annotations"] or []
            terms.update(a["term"] for a in annotations)
            if relation["subject_type"] == "molecular_activity" and relation["predicate"] in {
                "has_input",
                "has_output",
                "enabled_by",
            }:
                values = {a["value"] for a in annotations if a["term"] == SOURCE_RECORD_REFERENCE}
                expected = (
                    SOURCE_RECORD_SHA256_PREFIX
                    + payload_hashes[relation["relation_key"], ev["row_id"]]
                )
                require(values == {expected}, "Reaction evidence lost its exact source hash")
                shapes = {a["value"] for a in annotations if a["term"] == SOURCE_RECORD_TYPE}
                require(shapes == {"object"}, "Reaction source shape was not preserved")
                reaction_evidence += 1
            if result["resource"] == "connectomedb2025":
                roles = {
                    (a["scope"], a["value"]) for a in annotations if a["term"] == PARTICIPANT_ROLE
                }
                require(
                    roles
                    in (
                        {("subject", "ligand"), ("object", "receptor")},
                        {("subject", "receptor"), ("object", "ligand")},
                    ),
                    "Ligand/receptor evidence roles do not identify both canonical endpoints",
                )
    for entity in entities:
        terms.update(a["term"] for a in entity["annotations"] or [])
    if result["resource"] == "macdb":
        require(terms[TRAIT_TYPE] > 0, "MACdb trait type was not published")
    if result["resource"] in SOURCES[:5]:
        require(reaction_evidence > 0, "No reaction evidence in the bounded sample")
    return {
        "entities": len(entities),
        "relations": len(relations),
        "payloads": len(payloads),
        "reaction_evidence": reaction_evidence,
        "annotation_counts": dict(sorted(terms.items())),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--library-dir", required=True, type=Path)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()
    root = args.output_dir.resolve()
    library = args.library_dir.resolve(strict=True)
    require(
        root != library and not library.is_relative_to(root),
        "Output must not contain the reference",
    )
    require(not root.exists(), "Use a new, private output directory for the capped samples")
    root.mkdir(parents=True)
    os.environ["OMNIPATH_BUILD_DUCKDB_THREADS"] = "1"
    from omnipath_build.discovery import discover_datasets
    from omnipath_build.pipeline import build_resource
    from omnipath_postgres.releases import load_release

    report = {
        "capped_validation_sample": True,
        "max_records_per_resource": CAP,
        "version": args.version,
        "reference": str(library),
        "resources": {},
    }
    try:
        for source in SOURCES:
            entries = json.loads((args.inputs / f"{source}.json").read_text())
            require(0 < len(entries) <= CAP, f"Source cap violated: {source}")
            rows = defaultdict(list)
            origins = []
            for entry in entries:
                payload = entry["payload_json"]
                require(
                    hashlib.sha256(payload.encode()).hexdigest()
                    == entry["original_payload_sha256"],
                    f"Original input changed: {source}/{entry['original_row_id']}",
                )
                raw = json.loads(payload)
                require(isinstance(raw, dict), "Replay requires a parsed source object")
                dataset = entry["dataset"]
                origins.append(
                    {
                        "row_id": f"{dataset}:{len(rows[dataset])}",
                        "original_row_id": entry["original_row_id"],
                        "original_payload_sha256": entry["original_payload_sha256"],
                    }
                )
                rows[dataset].append(raw)
            cache = root / "pypath-data"
            _, discovered, _ = discover_datasets(source, datasets=list(rows), cache_dir=cache)
            consumed = Counter()

            def replay(dataset):
                def iterator(**_kwargs):
                    for row in rows[dataset]:
                        consumed[dataset] += 1
                        require(
                            sum(consumed.values()) <= CAP, "Resource exceeded the total record cap"
                        )
                        yield row

                return iterator

            started = perf_counter()
            with ExitStack() as stack:
                for dataset in discovered:
                    require(dataset.raw_dataset is not None, "Replay requires a raw dataset mapper")
                    stack.enter_context(
                        patch.object(dataset.raw_dataset, "raw", replay(dataset.dataset_name))
                    )
                stack.enter_context(
                    patch("pypath.inputs_v2.base.download_and_open", blocked_download)
                )
                result = build_resource(
                    source,
                    version=args.version,
                    output_dir=root,
                    library_dir=library,
                    cache_dir=cache,
                    datasets=list(rows),
                    max_records=CAP,
                    batch_workers=1,
                    batch_size=CAP,
                    max_batch_records=CAP,
                    max_batch_relations=CAP,
                    max_batch_bytes=64 * 1024**2,
                    resource_ram_bytes=1024**3,
                    min_free_disk_bytes=1024**3,
                    progress=False,
                )
            require(
                sum(consumed.values()) == len(entries),
                "The replay did not consume exactly the selected records",
            )
            resource = {
                "source_records": len(entries),
                "datasets": dict(consumed),
                "seconds": round(perf_counter() - started, 6),
                "original_records": origins,
                **inspect_resource(result),
            }
            report["resources"][source] = resource
            print(
                json.dumps(
                    {
                        "resource": source,
                        **{k: v for k, v in resource.items() if k != "original_records"},
                    }
                ),
                flush=True,
            )
        manifest = root / "release.json"
        manifest.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "version": args.version,
                    "resources": {s: args.version for s in SOURCES},
                },
                indent=2,
            )
            + "\n"
        )
        pinned = load_release(root, manifest)
        report.update(
            status="success",
            release_manifest=str(manifest),
            release_manifest_sha256=pinned.manifest_sha256,
            total_source_records=sum(r["source_records"] for r in report["resources"].values()),
        )
    except Exception as exc:
        report.update(status="failed", error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        (root / "inspection.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
