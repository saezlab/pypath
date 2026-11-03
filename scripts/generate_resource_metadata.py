"""Generate core's resource metadata from the inputs_v2 configs, or verify it with --check."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pkgutil

from omnipath_build.discovery import resolve_source_name
from omnipath_build.source_metadata import config_metadata
import pypath.inputs_v2

TARGET = (
    Path(__file__).resolve().parents[1] / "packages/core/omnipath_core/vocab/resource_metadata.yaml"
)


def generate() -> str:
    sources = {}
    for module in pkgutil.iter_modules(pypath.inputs_v2.__path__):
        if module.ispkg or module.name.startswith("_"):
            continue
        metadata = config_metadata(f"pypath.inputs_v2.{module.name}")
        if "website" in metadata:  # the module has a config
            sources[resolve_source_name(module.name)] = metadata
    return json.dumps(dict(sorted(sources.items())), indent=2, ensure_ascii=False) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    content = generate()
    if args.check:
        if TARGET.read_text() != content:
            raise SystemExit(
                "Resource metadata is stale; run scripts/generate_resource_metadata.py"
            )
    else:
        TARGET.write_text(content)


if __name__ == "__main__":
    main()
