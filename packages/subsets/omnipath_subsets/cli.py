"""Build the subset products of a release: Parquet from the published release, and the
PostgreSQL products that read a loaded release."""

import argparse
from dataclasses import asdict
import json
import os

import psycopg2

from omnipath_core import subset_layout as layout
from .build import PRODUCTS, build_subsets
from .parquet_products import build_release_subsets


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    build = commands.add_parser(
        "build",
        help="Write the Parquet subsets of a published release under <data-root>/subsets/",
    )
    build.add_argument("--release", required=True, help="A published release version")
    build.add_argument("--data-root", required=True)
    build.add_argument("--products", choices=layout.PRODUCTS, nargs="+", default=None)
    build.add_argument("--seed-directory", default=None, help="Seed files; the packaged ones by default")

    check = commands.add_parser(
        "check",
        help="Check the built subsets of a release against their contract and the release",
    )
    check.add_argument("--release", required=True)
    check.add_argument("--data-root", required=True)
    check.add_argument("--products", choices=layout.PRODUCTS, nargs="+", default=None)

    database = commands.add_parser(
        "build-postgres", help="Rebuild selected products from a loaded release"
    )
    database.add_argument("--database-url", default=os.environ.get("OMNIPATH_DATABASE_URL"))
    database.add_argument("--schema", default="omnipath")
    database.add_argument("--products", choices=PRODUCTS, nargs="+", default=PRODUCTS)
    database.add_argument(
        "--checkpoint-products",
        action="store_true",
        help="Commit each complete product separately; retry remaining products explicitly",
    )
    args = parser.parse_args(argv)

    if args.command == "build":
        try:
            manifests = build_release_subsets(
                args.data_root, args.release, args.products, seed_directory=args.seed_directory
            )
        except (layout.SubsetError, OSError) as exc:
            parser.exit(1, f"Subset build failed: {exc}\n")
        print(json.dumps(manifests, indent=2))
        return 0

    if args.command == "check":
        from .checks import check_release

        try:
            found = check_release(args.data_root, args.release, args.products)
        except (layout.SubsetError, OSError) as exc:
            parser.exit(1, f"Subset check failed: {exc}\n")
        print(
            json.dumps(
                {product: [f.as_dict() for f in findings] for product, findings in found.items()},
                indent=2,
            )
        )
        broken = [f for findings in found.values() for f in findings if f.level == "error"]
        return 1 if broken or not found else 0

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
