"""Command-line entry point for a single pinned PostgreSQL release."""

import argparse
from dataclasses import asdict
import json
import os

import duckdb
import psycopg2

from .loader import (
    finish_release as finish_release,
    load_release as load_release,
    PRODUCTS,
)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release_manifest", nargs="?", help="Pinned release JSON file or HTTPS URL")
    parser.add_argument(
        "--data-root", help="Local directory or HTTPS base URL containing resources/"
    )
    parser.add_argument("--database-url", default=os.environ.get("OMNIPATH_DATABASE_URL"))
    parser.add_argument(
        "--schema", default="omnipath", help="New destination schema or checkpoint to finish"
    )
    operation = parser.add_mutually_exclusive_group()
    operation.add_argument(
        "--checkpoint-base",
        action="store_true",
        help="Accepted for compatibility: main layout always commits its base before derivation",
    )
    operation.add_argument(
        "--finish",
        action="store_true",
        help="Finish an unpublished base checkpoint from PostgreSQL alone",
    )
    parser.add_argument(
        "--defer-constraints",
        action="store_true",
        help="Accepted for compatibility: main layout creates and validates keys/FKs after COPY",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="Optional source-record audit batch (capped at 1024)",
    )
    parser.add_argument(
        "--validate-source-records",
        action="store_true",
        help="Audit discarded raw source records and their owner/provenance references",
    )
    parser.add_argument(
        "--base-only",
        action="store_true",
        help="Stop after committing main base tables, constraints and indexes",
    )
    parser.add_argument(
        "--retain-published-provenance",
        action="store_true",
        help="Also load exact published crosswalks and repeated occurrences for inspection",
    )
    parser.add_argument(
        "--products",
        nargs="+",
        choices=PRODUCTS,
        default=PRODUCTS,
        help="Select downstream products; completed products are retained on --finish",
    )
    parser.add_argument("--duckdb-threads", type=int)
    parser.add_argument("--memory-limit", help="DuckDB working-memory limit")
    parser.add_argument(
        "--temp-directory",
        help="Local working filesystem; defaults next to a local manifest or in OS temporary storage for HTTPS",
    )
    args = parser.parse_args(argv)
    if not args.database_url:
        parser.error("Set OMNIPATH_DATABASE_URL or pass --database-url")
    if args.finish:
        if (
            any(
                value is not None
                for value in (
                    args.release_manifest,
                    args.data_root,
                    args.batch_size,
                    args.duckdb_threads,
                    args.memory_limit,
                    args.temp_directory,
                )
            )
            or args.validate_source_records
            or args.defer_constraints
            or args.base_only
            or args.retain_published_provenance
        ):
            parser.error(
                "--finish accepts only database/schema selection; loading and source-audit options do not apply"
            )
    elif args.release_manifest is None or args.data_root is None:
        parser.error("Loading requires release_manifest and --data-root")
    try:
        if args.finish:
            result = finish_release(
                args.database_url, schema=args.schema, products=tuple(args.products)
            )
        else:
            if args.validate_source_records or args.batch_size is not None:
                parser.error(
                    "Main layout consumes resolved relational Parquets; raw record auditing belongs to artifact preparation"
                )
            result = load_release(
                args.data_root,
                args.release_manifest,
                args.database_url,
                schema=args.schema,
                duckdb_threads=1 if args.duckdb_threads is None else args.duckdb_threads,
                memory_limit="1GB" if args.memory_limit is None else args.memory_limit,
                temp_directory=args.temp_directory,
                base_only=args.base_only,
                products=tuple(args.products),
                defer_constraints=True,
                retain_published_provenance=args.retain_published_provenance,
            )
    except (ValueError, OSError, psycopg2.Error, duckdb.Error) as exc:
        parser.exit(1, f"Release load failed: {exc}\n")
    print(json.dumps(asdict(result), indent=2))
    return 0
