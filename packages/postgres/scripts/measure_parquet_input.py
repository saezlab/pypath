#!/usr/bin/env python3
"""Measure a pinned resource's validation traffic and DuckDB projection, without PostgreSQL.

Provide an existing single-resource release; no source/reference build occurs.
Use --max-bytes to explicitly allow a larger scan. Counts and times are printed,
never credentials, paths or raw source data. No persistent Parquet cache is made.
"""

import argparse
from dataclasses import asdict
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import duckdb
from omnipath_core.versioning import validate_build_manifest, validate_release_manifest
from omnipath_postgres.locations import join_location, measure_http_transfer, read_bytes
from omnipath_postgres.projection import prepare_aligned_release
from omnipath_postgres.releases import load_release, verify_release


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest")
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--max-bytes", type=int, default=64 * 1024**2)
    args = parser.parse_args()
    release = validate_release_manifest(json.loads(read_bytes(args.manifest)), strict_fields=True)
    if len(release.resources) != 1 or release.references:
        parser.error("Use a single-resource release without taxonomy for this bounded measurement")
    resource, version = next(iter(release.resources.items()))
    build = validate_build_manifest(
        json.loads(
            read_bytes(
                join_location(args.data_root, "resources", resource, version, "build_manifest.json")
            )
        )
    )
    artifact_bytes = sum(item.size_bytes for item in build.files.values())
    if args.max_bytes < 1 or artifact_bytes > args.max_bytes:
        parser.error(
            "Pinned artifacts exceed --max-bytes; choose a smaller release or explicit budget"
        )
    with measure_http_transfer() as initial:
        pinned = load_release(args.data_root, args.manifest)
    with TemporaryDirectory(prefix="omnipath-input-benchmark-") as folder:
        with duckdb.connect(str(Path(folder) / "projection.duckdb")) as connection:
            connection.execute("SET threads=2")
            connection.execute("SET memory_limit='1GB'")
            started = perf_counter()
            plan = prepare_aligned_release(connection, pinned, dimension_rows={})
            projection_seconds = perf_counter() - started
            counts = dict(plan.counts)
    with measure_http_transfer() as final:
        verify_release(pinned)
    print(
        json.dumps(
            {
                "resource": resource,
                "version": version,
                "artifact_bytes": artifact_bytes,
                "initial_validation": asdict(initial),
                "final_validation": asdict(final),
                "projection_seconds": projection_seconds,
                "projected_counts": counts,
                "note": "Validation bytes include metadata and hashes; DuckDB transfer bytes are not counted.",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
