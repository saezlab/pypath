"""Publish an explicit resource snapshot from already built local artifacts.

With ``--build-database`` it then writes the release's Parquet subsets
(``subsets/<version>/``), which need the ``omnipath-subsets`` package, and then the DuckDB
database the API serves the release from (``serving/<version>.duckdb``), which copies them
in. It keeps the newest ``--keep`` databases. ``--skip-subsets`` serves the release without
subsets.
"""

import argparse
import json
import sys
from pathlib import Path

from omnipath_api.queries.store.inventory import ReleaseStore


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--build-database", action="store_true")
    parser.add_argument(
        "--skip-subsets",
        action="store_true",
        help="with --build-database: do not build the release's subsets first",
    )
    parser.add_argument("--keep", type=int, default=3)
    args = parser.parse_args(argv)
    result = ReleaseStore(args.data_root).publish(json.loads(args.manifest.read_text()))
    print(json.dumps({"version": result["version"], "resources": len(result["resources"])}))
    if args.build_database:
        from omnipath_api.release_database import build_served_database

        if not args.skip_subsets:
            from omnipath_subsets.parquet_products import build_release_subsets

            build_release_subsets(args.data_root, result["version"])

        record = build_served_database(
            args.data_root,
            result["version"],
            keep=args.keep or None,
            progress=lambda message: print(message, file=sys.stderr, flush=True),
        )
        print(json.dumps(record))


if __name__ == "__main__":
    main()
