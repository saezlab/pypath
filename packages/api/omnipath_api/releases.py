"""Immutable OmniPath release manifests referencing existing resource versions.

No graph files are copied when a release is published. This module also provides
an offline CLI: python -m releases --data-root data publish release.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from omnipath_api.queries.store.inventory import (
    FILES,
    ReleaseStore,
    safe_name,
    validate_manifest,
)

__all__ = ["FILES", "ReleaseStore", "safe_name", "validate_manifest"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default="data")
    commands = parser.add_subparsers(dest="command", required=True)
    publish = commands.add_parser("publish")
    publish.add_argument("manifest", type=Path)
    args = parser.parse_args()
    store = ReleaseStore(args.data_root)
    try:
        result = store.publish(json.loads(args.manifest.read_text()))
        print(f"Published OmniPath {result['version']}")
    except (ValueError, OSError) as exc:
        parser.exit(1, f"{exc}\n")


if __name__ == "__main__":
    main()
