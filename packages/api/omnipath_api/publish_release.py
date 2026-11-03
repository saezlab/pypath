"""Publish an explicit resource snapshot from already built local artifacts."""

import argparse
import json
from pathlib import Path

from omnipath_api.queries.store.inventory import ReleaseStore


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    args = parser.parse_args(argv)
    result = ReleaseStore(args.data_root).publish(json.loads(args.manifest.read_text()))
    print(json.dumps({"version": result["version"], "resources": len(result["resources"])}))


if __name__ == "__main__":
    main()
