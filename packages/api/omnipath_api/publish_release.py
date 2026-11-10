"""Publish an explicit resource snapshot from already built local artifacts.

With ``--build-database`` it then writes the DuckDB database the API serves the release
from (``serving/<version>.duckdb``) and keeps the newest ``--keep`` databases.
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
    parser.add_argument("--keep", type=int, default=3)
    args = parser.parse_args(argv)
    result = ReleaseStore(args.data_root).publish(json.loads(args.manifest.read_text()))
    print(json.dumps({"version": result["version"], "resources": len(result["resources"])}))
    if args.build_database:
        from omnipath_api.release_database import build_served_database

        record = build_served_database(
            args.data_root,
            result["version"],
            keep=args.keep or None,
            progress=lambda message: print(message, file=sys.stderr, flush=True),
        )
        print(json.dumps(record))


if __name__ == "__main__":
    main()
