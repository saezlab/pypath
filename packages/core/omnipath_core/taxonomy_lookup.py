"""Offline NCBI name lookup for the moving latest view (not pinned releases).

Prepare explicitly: python -m omnipath_core.taxonomy_lookup --data-root /data
The source is the already cached taxdump.tar.gz; no request downloads taxonomy.
"""

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import sqlite3
import tarfile
import tempfile
import threading
from functools import lru_cache


def prepare(root):
    directory = Path(root) / "references/taxonomy"
    archive = directory / "taxdump.tar.gz"
    with archive.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    target = directory / f"names-{digest}.sqlite"
    if not target.exists():
        fd, staging = tempfile.mkstemp(dir=directory, suffix=".sqlite")
        os.close(fd)
        try:
            with sqlite3.connect(staging) as db, tarfile.open(archive, "r:gz") as dump:
                db.execute(
                    "CREATE TABLE names (id TEXT PRIMARY KEY, scientific TEXT, common TEXT) WITHOUT ROWID"
                )
                db.execute("CREATE TABLE merged (id TEXT PRIMARY KEY, current TEXT) WITHOUT ROWID")
                with dump.extractfile("names.dmp") as stream:
                    batch = []
                    for line in io.TextIOWrapper(stream, encoding="utf-8"):
                        taxon, name, _, kind, *_ = [v.strip() for v in line.split("|")]
                        if kind not in ("scientific name", "genbank common name"):
                            continue
                        batch.append(
                            (
                                taxon,
                                name if kind == "scientific name" else None,
                                name if kind == "genbank common name" else None,
                            )
                        )
                        if len(batch) >= 5000:
                            db.executemany(
                                "INSERT INTO names VALUES (?,?,?) ON CONFLICT(id) DO UPDATE SET scientific=coalesce(excluded.scientific,scientific), common=coalesce(excluded.common,common)",
                                batch,
                            )
                            batch.clear()
                    db.executemany(
                        "INSERT INTO names VALUES (?,?,?) ON CONFLICT(id) DO UPDATE SET scientific=coalesce(excluded.scientific,scientific), common=coalesce(excluded.common,common)",
                        batch,
                    )
                with dump.extractfile("merged.dmp") as stream:
                    db.executemany(
                        "INSERT INTO merged VALUES (?,?)",
                        (
                            [v.strip() for v in line.split("|")][:2]
                            for line in io.TextIOWrapper(stream, encoding="utf-8")
                        ),
                    )
            os.replace(staging, target)
        finally:
            Path(staging).unlink(missing_ok=True)
    fd, staging = tempfile.mkstemp(dir=directory)
    with os.fdopen(fd, "w") as stream:
        json.dump({"file": target.name, "source_sha256": digest}, stream)
    os.replace(staging, directory / "latest-names.json")
    return target


def lookup_path(root):
    directory = Path(root) / "references/taxonomy"
    try:
        name = json.loads((directory / "latest-names.json").read_text())["file"]
        if Path(name).name != name:
            return None
        path = directory / name
        return str(path) if path.is_file() else None
    except (OSError, ValueError, KeyError):
        return None


class Names:
    def __init__(self, path):
        self.path = path
        self.local = threading.local()

    @lru_cache(maxsize=16384)
    def get(self, taxon, default=None):
        if not taxon:
            return default
        if not hasattr(self.local, "db"):
            self.local.db = sqlite3.connect(
                Path(self.path).as_uri() + "?mode=ro&immutable=1", uri=True
            )
        db = self.local.db
        seen = set()
        while taxon not in seen:
            seen.add(taxon)
            row = db.execute(
                "SELECT coalesce(common,scientific) FROM names WHERE id=?", (taxon,)
            ).fetchone()
            if row and row[0]:
                return row[0]
            row = db.execute("SELECT current FROM merged WHERE id=?", (taxon,)).fetchone()
            if not row:
                break
            taxon = row[0]
        return default


@lru_cache(maxsize=4)
def read_lookup(path):
    return Names(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default="data")
    print(prepare(parser.parse_args().data_root))
