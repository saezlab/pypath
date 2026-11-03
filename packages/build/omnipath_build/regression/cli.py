"""``omnipath-build resolution-regression {extract,resolve,diff}``."""

from __future__ import annotations

import argparse
import sys
from typing import Any

from .store import DEFAULT_SAMPLE


def add_parser(subparsers: Any) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "resolution-regression",
        help="Extract pre-resolution observations, resolve them with a runtime, diff results",
    )
    sub = parser.add_subparsers(dest="regression_command", required=True)

    extract = sub.add_parser(
        "extract", help="Extract distinct observations (queries and votes) per resource"
    )
    extract.add_argument("--resources", nargs="+", required=True, help="Resource names")
    extract.add_argument("--output", required=True, help="Observation root directory")
    extract.add_argument("--max-records", type=int, help="Rows per dataset (quick runs)")
    extract.add_argument("--workers", type=int, default=2, help="Extraction processes")
    extract.add_argument("--batch-records", type=int, default=2000, help="Rows per worker batch")
    extract.add_argument("--datasets", nargs="+", help="Restrict to these datasets")
    extract.add_argument("--cache-dir", help="pypath download cache (default: env)")
    extract.add_argument("--memory-limit", default="2GB", help="DuckDB limit for the final merge")
    extract.add_argument(
        "--skip-existing", action="store_true", help="Skip resources already extracted OK"
    )

    resolve = sub.add_parser("resolve", help="Resolve stored observations with a runtime")
    resolve.add_argument(
        "--runtime", required=True, help="Identity directory (omnipath-identity-v2)"
    )
    resolve.add_argument("--observations", required=True, help="Observation root directory")
    resolve.add_argument("--output", required=True, help="Result directory")
    resolve.add_argument("--resources", nargs="+", help="Default: every extracted resource")
    resolve.add_argument("--batch-size", type=int, default=5000)
    resolve.add_argument("--cache-dir", help="Identity runtime entity cache directory")
    resolve.add_argument(
        "--sample-per-resource",
        type=int,
        default=DEFAULT_SAMPLE,
        help="Deterministic sample (md5 of the fingerprint) per resource and library; "
        f"0 = everything (default: {DEFAULT_SAMPLE})",
    )
    resolve.add_argument("--force", action="store_true", help="Redo resources with results")

    diff = sub.add_parser("diff", help="Compare two result sets")
    diff.add_argument("--observations", required=True)
    diff.add_argument("--a", required=True, help="Baseline result directory")
    diff.add_argument("--b", required=True, help="Candidate result directory")
    diff.add_argument("--output", required=True)
    diff.add_argument("--resources", nargs="+")
    diff.add_argument("--examples", type=int, default=3, help="Examples per change type")
    diff.add_argument(
        "--sample-per-resource",
        type=int,
        default=DEFAULT_SAMPLE,
        help="Must match the resolve runs; 0 = everything (default: %d)" % DEFAULT_SAMPLE,
    )
    diff.add_argument("--memory-limit", default="2GB")
    return parser


def run(parsed: argparse.Namespace) -> int:
    command = parsed.regression_command
    if command == "extract":
        from .extract import extract_resources

        results = extract_resources(
            parsed.resources,
            parsed.output,
            skip_existing=parsed.skip_existing,
            max_records=parsed.max_records,
            workers=parsed.workers,
            batch_records=parsed.batch_records,
            cache_dir=parsed.cache_dir,
            datasets=parsed.datasets,
            memory_limit=parsed.memory_limit,
        )
        bad = {r: i["status"] for r, i in results.items() if i["status"] != "ok"}
        for resource, status in bad.items():
            print(f"{resource}: {status}: {results[resource]['error']}", file=sys.stderr)
        return 1 if bad else 0
    if command == "resolve":
        from .resolve import resolve_all

        metrics = resolve_all(
            parsed.runtime,
            parsed.observations,
            parsed.output,
            resources=parsed.resources,
            batch_size=parsed.batch_size,
            cache_dir=parsed.cache_dir,
            force=parsed.force,
            sample=parsed.sample_per_resource or None,
        )
        failed = [r for r, m in metrics["resources"].items() if "error" in m]
        for resource in failed:
            print(f"{resource}: {metrics['resources'][resource]['error']}", file=sys.stderr)
        return 1 if failed else 0
    if command == "diff":
        from .diff import compare

        summary = compare(
            parsed.observations,
            parsed.a,
            parsed.b,
            parsed.output,
            resources=parsed.resources,
            memory_limit=parsed.memory_limit,
            examples=parsed.examples,
            sample=parsed.sample_per_resource or None,
        )
        for resource, why in summary["skipped"].items():
            print(f"{resource}: skipped ({why})", file=sys.stderr)
        print(f"wrote {parsed.output}/summary.md")
        return 0
    raise ValueError(command)
