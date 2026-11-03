"""Shared plumbing for the hub-index and identity builds: stages, DuckDB, SQL fragments."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import time
from pathlib import Path

import duckdb

from omnipath_resolver.canonical.identifiers import normalize_id_sql

CHEMICAL = (
    "chebi",
    "pubchem",
    "chembl",
    "hmdb",
    "lipidmaps",
    "swisslipids",
    "bigg",
    "metanetx",
    "refmet",
    "ramp",
)
PROTEIN = ("uniprot", "entrez", "ramp_gene")
# The InChIKeys of an empty structure identify nothing.
EMPTY_KEYS = ("MOSFIJXAXDLOML-UHFFFAOYSA-N", "MOSFIJXAXDLOML-UHFFFAOYNA-N")
PARTS = tuple(f"{n:02x}" for n in range(256))
# Rhea master reactions are the anchors; MetaNetX reactions attach to them (proposal, Reactions).
REACTION = ("rhea", "metanetx_reaction")
HUBS = CHEMICAL + PROTEIN + REACTION
# Names only: ontology terms keep their accession as identity, nothing is matched against it.
ONTOLOGY = ("ontology",)
# Rule 4: preferred hub for the id of a grouped entity (earlier wins).
HUB_PREFERENCE = (
    "chebi",
    "lipidmaps",
    "swisslipids",
    "hmdb",
    "chembl",
    "pubchem",
    "kegg",
    "metanetx",
    "bigg",
    "refmet",
    "ramp",
    "ramp_gene",
    "metanetx_reaction",
)
GOSLIN_HUBS = ("swisslipids", "lipidmaps", "hmdb", "chebi", "refmet")
INCHIKEY_RE = "^[A-Z]{14}-[A-Z]{10}-[A-Z]$"
STRUCTURE_LEVELS = ("full_structure", "complete_structure")


def quote(value) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def scan(path) -> str:
    return f"read_parquet({quote(path)})"


def norm(ns, ident) -> str:
    """SQL normalizing identifier `ident` of namespace `ns`, as the matcher normalizes it."""
    return normalize_id_sql(ns, ident)


def hub_rank(column="hub") -> str:
    return (
        f"CASE {column} "
        + " ".join(f"WHEN {quote(h)} THEN {i}" for i, h in enumerate(HUB_PREFERENCE))
        + f" ELSE {len(HUB_PREFERENCE)} END"
    )


def sql_list(values) -> str:
    return "(" + ",".join(quote(v) for v in values) + ")"


def part_of(value: str) -> str:
    return hashlib.md5(value.encode(), usedforsecurity=False).hexdigest()[:2]


def code_sha256(*names) -> str:
    """Hash of the named modules of this package (build code identity)."""
    digest = hashlib.sha256()
    for name in sorted(names):
        path = Path(__file__).with_name(name)
        digest.update(name.encode() + b"\0" + path.read_bytes())
    return digest.hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def files(*patterns) -> str:
    """DuckDB read_parquet over glob patterns (hive columns off: partition dirs are not columns)."""
    return "read_parquet([" + ",".join(quote(p) for p in patterns) + "],hive_partitioning=false)"


def part_dirs(root: Path, parts) -> list[str]:
    """Glob patterns for the existing `part=XX` directories of the requested parts."""
    return [str(root / f"part={p}" / "*.parquet") for p in parts if (root / f"part={p}").is_dir()]


def log(event, **fields):
    print(
        json.dumps(dict(time=time.strftime("%FT%TZ", time.gmtime()), event=event, **fields)),
        flush=True,
    )


class Context:
    """DuckDB settings and checkpointed stages of one build directory (`out`)."""

    def __init__(self, out, memory="7GB", threads=6, goslin_cache=None, min_free_gib=50):
        self.out = Path(out)
        self.work = self.out / "work"
        self.work.mkdir(parents=True, exist_ok=True)
        self.memory, self.threads = memory, int(threads)
        self.min_free = min_free_gib * 1024**3
        self.goslin_cache = Path(goslin_cache) if goslin_cache else self.work / "goslin-cache"
        self.timings: dict[str, float] = {}
        self.info: dict[str, dict] = {}

    def guard(self):
        if shutil.disk_usage(self.out).free < self.min_free:
            raise RuntimeError("Identity build stopped at the free-disk reserve")

    def connect(self, name: str):
        self.guard()
        c = duckdb.connect()
        c.execute(f"SET memory_limit={quote(self.memory)}")
        c.execute(f"SET threads={self.threads}")
        c.execute("SET preserve_insertion_order=false")
        spill = self.work / "spill" / name
        spill.mkdir(parents=True, exist_ok=True)
        c.execute(f"SET temp_directory={quote(spill)}")
        c.execute("SET max_temp_directory_size='64GiB'")
        c.execute("SET partitioned_write_max_open_files=1024")
        c.execute("SET partitioned_write_flush_threshold=131072")
        return c

    def copy(self, c, query: str, path: Path, options: str = "") -> int:
        """COPY a query to one Parquet file atomically (temp file, then rename)."""
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".tmp")
        n = c.execute(
            f"COPY ({query}) TO {quote(tmp)} (FORMAT PARQUET, COMPRESSION ZSTD{options})"
        ).fetchone()[0]
        os.replace(tmp, path)
        return n

    def copy_partitioned(self, c, query: str, directory: Path, columns: str, pattern: str):
        """COPY a query into hive-partitioned directories; `pattern` keeps file names unique."""
        directory.mkdir(parents=True, exist_ok=True)
        return c.execute(
            f"COPY ({query}) TO {quote(directory)} (FORMAT PARQUET, COMPRESSION ZSTD, "
            f"PARTITION_BY ({columns}), FILENAME_PATTERN {quote(pattern + '{i}')}, "
            "OVERWRITE_OR_IGNORE true, ROW_GROUP_SIZE 131072)"
        ).fetchone()[0]

    def stage(self, name: str, outputs, fn):
        """Run `fn(self, directory)` unless its _SUCCESS.json and listed outputs are present."""
        directory = self.work / name
        directory.mkdir(parents=True, exist_ok=True)
        marker = directory / "_SUCCESS.json"
        outputs = [Path(o) for o in outputs]
        if marker.exists() and all(o.exists() for o in outputs):
            done = json.loads(marker.read_text())
            self.timings[name] = done["seconds"]
            self.info[name] = done["info"]
            log("stage_reused", stage=name, seconds=done["seconds"])
            return done["info"]
        self.guard()
        marker.unlink(missing_ok=True)
        log("stage_start", stage=name)
        start = time.monotonic()
        info = fn(self, directory) or {}
        seconds = round(time.monotonic() - start, 2)
        marker.write_text(json.dumps(dict(seconds=seconds, info=info), indent=2))
        shutil.rmtree(self.work / "spill" / name, ignore_errors=True)
        self.timings[name] = seconds
        self.info[name] = info
        log("stage_done", stage=name, seconds=seconds, **info)
        return info


def part_done(directory: Path, name: str) -> bool:
    return (directory / "done" / name).exists()


def finish_part(directory: Path, name: str):
    (directory / "done").mkdir(exist_ok=True)
    (directory / "done" / name).write_text("")
