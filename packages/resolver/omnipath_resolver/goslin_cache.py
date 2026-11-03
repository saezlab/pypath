"""Persistent exact-name results, isolated by normalization code and grammar content."""

from __future__ import annotations

import ast
import hashlib
import inspect
from pathlib import Path
import sqlite3
import textwrap

FIELDS = ("name", "goslin", "level", "status", "specificity")


def cache_fingerprint(normalizer):
    import pygoslin

    digest = hashlib.sha256(b"goslin-name-cache-v1\0")
    digest.update(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(normalizer)))).encode())
    root = Path(pygoslin.__file__).parent
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.suffix in (".py", ".g4", ".G4", ".csv"):
            digest.update(str(path.relative_to(root)).encode() + b"\0")
            digest.update(path.read_bytes())
    return digest.hexdigest()


class NormalizationCache:
    def __init__(self, directory, normalizer):
        self.fingerprint = cache_fingerprint(normalizer)
        self.path = Path(directory) / self.fingerprint / "results.sqlite"
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.connection = sqlite3.connect(self.path, timeout=60)
        self.connection.row_factory = sqlite3.Row
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.execute("""CREATE TABLE IF NOT EXISTS results (
            name TEXT PRIMARY KEY, goslin TEXT, level TEXT,
            status TEXT NOT NULL, specificity INTEGER NOT NULL) WITHOUT ROWID""")
        self.connection.commit()

    def lookup(self, names):
        found = []
        for offset in range(0, len(names), 500):
            batch = names[offset : offset + 500]
            placeholders = ",".join("?" for _ in batch)
            found.extend(
                dict(row)
                for row in self.connection.execute(
                    f"SELECT name,goslin,level,status,specificity FROM results WHERE name IN ({placeholders})",
                    batch,
                )
            )
        return found

    def store(self, results):
        with self.connection:
            self.connection.executemany(
                "INSERT OR IGNORE INTO results VALUES (?,?,?,?,?)",
                [tuple(row[field] for field in FIELDS) for row in results],
            )

    def close(self):
        self.connection.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
