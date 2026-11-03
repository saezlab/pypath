"""Generate the checked-in API contract, or verify it with --check."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile

from omnipath_api.server import create_app


def generate() -> str:
    with tempfile.TemporaryDirectory(prefix="omnipath-openapi-") as directory:
        app = create_app(data_root=directory)
        try:
            return json.dumps(app.openapi(), indent=2, sort_keys=True) + "\n"
        finally:
            app.state.queries.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    target = Path(__file__).resolve().parents[1] / "openapi.json"
    content = generate()
    if args.check:
        if target.read_text() != content:
            raise SystemExit(
                "OpenAPI is stale; run scripts/generate_openapi.py and pnpm generate:types"
            )
    else:
        target.write_text(content)


if __name__ == "__main__":
    main()
