"""Release-pinned, offline NCBI names. Entity identities always retain their taxon ID."""

from __future__ import annotations

import hashlib
import io
import os
import tarfile
import tempfile
import warnings
from functools import lru_cache
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

SOURCE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"


def taxon_names(archive: Path, taxon_ids: set[str]) -> list[dict]:
    """Read only requested names; merged IDs inherit names without changing identity."""
    with tarfile.open(archive, "r:gz") as dump:
        merged = {}
        with dump.extractfile("merged.dmp") as stream:
            for line in io.TextIOWrapper(stream, encoding="utf-8"):
                fields = [v.strip() for v in line.split("|")]
                merged[fields[0]] = fields[1]
        canonical = {}
        for taxon in taxon_ids:
            current, seen = taxon, set()
            while current in merged and current not in seen:
                seen.add(current)
                current = merged[current]
            canonical[taxon] = current
        wanted = set(canonical.values())
        names: dict[str, dict[str, str]] = {}
        with dump.extractfile("names.dmp") as stream:
            for line in io.TextIOWrapper(stream, encoding="utf-8"):
                taxon, name, _, kind, *_ = [v.strip() for v in line.split("|")]
                if taxon in wanted and kind in ("scientific name", "genbank common name"):
                    names.setdefault(taxon, {})[kind] = name
    return [
        dict(
            taxon_id=taxon,
            current_taxon_id=current,
            scientific_name=names.get(current, {}).get("scientific name"),
            common_name=names.get(current, {}).get("genbank common name"),
        )
        for taxon, current in sorted(canonical.items())
    ]


def prepare_reference(root: Path, manifest: dict) -> dict:
    """Create an immutable subset at release publication, using a cached NCBI dump."""
    archive = root / "references/taxonomy/taxdump.tar.gz"
    if not archive.is_file():
        raise FileNotFoundError(
            f"Cache {SOURCE_URL} at {archive} before publishing a taxonomy reference"
        )
    ids: set[str] = set()
    for source, version in manifest["resources"].items():
        for filename in ("entity.parquet", "relation.parquet"):
            path = root / "resources" / source / version / filename
            for batch in pq.ParquetFile(path).iter_batches(columns=["taxon"]):
                ids.update(
                    str(value)
                    for value in batch.column(0).to_pylist()
                    if value and str(value) != "0"
                )
    rows = taxon_names(archive, ids)
    missing = [row["taxon_id"] for row in rows if not row["scientific_name"]]
    if missing:
        warnings.warn(
            f"NCBI taxonomy names missing for {len(missing)} IDs: {', '.join(missing)}",
            stacklevel=2,
        )
    source_hash = hashlib.sha256()
    with archive.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            source_hash.update(block)
    schema = pa.schema(
        [
            (name, pa.string())
            for name in ("taxon_id", "current_taxon_id", "scientific_name", "common_name")
        ]
    )
    table = pa.Table.from_pylist(rows, schema=schema)
    buffer = pa.BufferOutputStream()
    pq.write_table(table, buffer, compression="zstd")
    content = buffer.getvalue().to_pybytes()
    digest = hashlib.sha256(content).hexdigest()
    directory = root / "references/taxonomy" / digest
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / "taxonomy.parquet"
    fd, filename = tempfile.mkstemp(dir=directory)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(content)
        os.chmod(filename, 0o644)
        try:
            os.link(filename, target)
        except FileExistsError:
            if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError("Existing taxonomy reference checksum mismatch")
    finally:
        Path(filename).unlink(missing_ok=True)
    return dict(
        version=digest,
        source_url=SOURCE_URL,
        source_sha256=source_hash.hexdigest(),
        taxon_count=len(rows),
        missing_taxon_ids=missing,
    )


@lru_cache(maxsize=16)
def read_reference(path: str) -> dict[str, str]:
    """Cache the small release subset, never the full NCBI taxonomy in API memory."""
    return {
        row["taxon_id"]: row["common_name"] or row["scientific_name"]
        for row in pq.read_table(path).to_pylist()
        if row["common_name"] or row["scientific_name"]
    }
