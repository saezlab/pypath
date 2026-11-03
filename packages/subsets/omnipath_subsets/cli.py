"""Build release-pinned product tables without parsing or resolving resources."""

import argparse
from dataclasses import asdict
import json
import os

import psycopg2

from .build import PRODUCTS, build_subsets


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build", help="Rebuild selected products from a loaded release")
    build.add_argument("--database-url", default=os.environ.get("OMNIPATH_DATABASE_URL"))
    build.add_argument("--schema", default="omnipath")
    build.add_argument("--products", choices=PRODUCTS, nargs="+", default=PRODUCTS)
    build.add_argument(
        "--checkpoint-products",
        action="store_true",
        help="Commit each complete product separately; retry remaining products explicitly",
    )
    args = parser.parse_args(argv)
    if not args.database_url:
        parser.error("Set OMNIPATH_DATABASE_URL or pass --database-url")
    try:
        result = build_subsets(
            args.database_url,
            args.schema,
            products=args.products,
            checkpoint_products=args.checkpoint_products,
        )
    except (ValueError, OSError, psycopg2.Error) as exc:
        parser.exit(1, f"Subset build failed: {exc}\n")
    print(json.dumps(asdict(result), indent=2))
    return 0
