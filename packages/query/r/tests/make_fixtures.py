"""Write the R tests' fixtures: the Python client tests' resources, as files.

    uv run python packages/query/r/tests/make_fixtures.py

Writes ``tests/testthat/fixtures/``: ``server/<resource>/<version>/*.parquet`` (the
published and serving tables of three small resources from the core fixtures, and a
second version of one), ``catalog-1.json`` and ``catalog-2.json`` (the catalog before
and after that version, as ``/resources`` returns it) and ``files-<resource>.json``
(``/resources/<resource>/files`` with raw evidence). The tests replace the package's
two HTTP functions with readers of these files.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "packages" / "query" / "python" / "tests"))

from test_client import Service  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "testthat" / "fixtures"


def main() -> None:
    if FIXTURES.exists():
        shutil.rmtree(FIXTURES)
    FIXTURES.mkdir(parents=True)
    service = Service(FIXTURES / "server")

    def catalog():
        return {"resources": [service.record(name) for name in service.versions]}

    (FIXTURES / "catalog-1.json").write_text(json.dumps(catalog(), indent=1) + "\n")
    for name in service.versions:
        (FIXTURES / f"files-{name}-1.json").write_text(
            json.dumps(service.record(name, evidence=True), indent=1) + "\n"
        )
    service.publish("alpha", "2.0")
    (FIXTURES / "catalog-2.json").write_text(json.dumps(catalog(), indent=1) + "\n")
    for name in service.versions:
        (FIXTURES / f"files-{name}-2.json").write_text(
            json.dumps(service.record(name, evidence=True), indent=1) + "\n"
        )
    print(f"wrote {sum(1 for _ in FIXTURES.rglob('*') if _.is_file())} files in {FIXTURES}")


if __name__ == "__main__":
    main()
