"""ChEMBL hub: directly queried from the ChEMBL SQLite database."""

from __future__ import annotations

import logging
import sqlite3
import subprocess
from pathlib import Path
from typing import Iterator

from ..common import explode_record
from ..schema import CHEMICAL_TAXON
from ..writer import HubParquetWriter

logger = logging.getLogger(__name__)

CHEMBL_VERSION = 36
CHEMBL_TAR_URL = (
    f"https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/"
    f"chembl_{CHEMBL_VERSION}/chembl_{CHEMBL_VERSION}_sqlite.tar.gz"
)


def _ensure_sqlite_db() -> Path:
    from ...discovery import setup_pypath_cache

    cache = setup_pypath_cache()
    target_dir = cache / "chembl"
    candidates = [target_dir, cache]
    target_dir.mkdir(parents=True, exist_ok=True)

    # 1. Check for existing complete .db file
    for directory in candidates:
        if not directory.is_dir():
            continue
        db_files = [
            directory / f"chembl_{CHEMBL_VERSION}.db",
            directory / f"ChEMBL_SQLite_{CHEMBL_VERSION}.sqlite",
        ]
        for db in db_files:
            if (
                db.is_file() and db.stat().st_size > 10_000_000_000
            ):  # Valid uncompressed ChEMBL is ~29 GB
                return db

    # 2. Check for existing .tar.gz archive and extract
    tar_path: Path | None = None
    for directory in candidates:
        if not directory.is_dir():
            continue
        tar_files = sorted(directory.glob("*.tar.gz"))
        if tar_files:
            tar_path = tar_files[-1]
            break

    if tar_path is None or not tar_path.exists():
        tar_path = target_dir / f"chembl_{CHEMBL_VERSION}_sqlite.tar.gz"
        logger.info("Downloading ChEMBL %s SQLite archive from %s", CHEMBL_VERSION, CHEMBL_TAR_URL)
        subprocess.run(
            [
                "curl",
                "-C",
                "-",
                "-L",
                "--retry",
                "10",
                "--retry-delay",
                "3",
                "-o",
                str(tar_path),
                CHEMBL_TAR_URL,
            ],
            check=True,
        )

    logger.info("Extracting %s into %s", tar_path, target_dir)
    subprocess.run(
        ["tar", "-xzf", str(tar_path), "-C", str(target_dir), "--strip-components=2"],
        check=True,
    )
    db_files = sorted(target_dir.glob("*.db")) or sorted(target_dir.glob("*.sqlite"))
    if not db_files:
        raise FileNotFoundError(f"No .db file found after extracting {tar_path}")
    return db_files[-1]


def _iter_sqlite_molecules(db_path: Path) -> Iterator[dict[str, str | None]]:
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    query = """
        SELECT
            m.chembl_id,
            m.pref_name,
            s.canonical_smiles,
            s.standard_inchi,
            s.standard_inchi_key
        FROM molecule_dictionary m
        LEFT JOIN compound_structures s ON m.molregno = s.molregno
    """
    cursor.execute(query)
    for row in cursor:
        yield {
            "chembl_id": row[0],
            "pref_name": row[1],
            "smiles": row[2],
            "inchi": row[3],
            "inchikey": row[4],
        }
    conn.close()


def emit(writer: HubParquetWriter) -> None:
    db_path = _ensure_sqlite_db()
    for mol in _iter_sqlite_molecules(db_path):
        if writer.full:
            return
        hub_id = str(mol.get("chembl_id") or "").strip()
        if not hub_id:
            continue
        if not explode_record(
            writer,
            hub_id=hub_id,
            hub_type="chembl",
            backend="chembl",
            taxonomy_id=CHEMICAL_TAXON,
            fields={
                "name": mol.get("pref_name"),
                "smiles": mol.get("smiles"),
                "inchi": mol.get("inchi"),
                "inchikey": mol.get("inchikey"),
            },
        ):
            return
