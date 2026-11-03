"""Cached Standard InChI derivation from explicitly supplied molecular SMILES."""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdinchi

POLICY_VERSION = "standard-inchi-smiles-v1"


def fingerprint():
    return hashlib.sha256(
        json.dumps(
            [
                POLICY_VERSION,
                rdBase.rdkitVersion,
                rdinchi.GetInchiVersion(),
            ]
        ).encode()
    ).hexdigest()


def derive(smiles: str) -> dict:
    """Keep all components, isotopes and specified stereochemistry; no cleanup."""
    result = dict(
        smiles=smiles,
        inchikey=None,
        inchi=None,
        status="invalid_smiles",
        message="",
        rdkit_version=rdBase.rdkitVersion,
        inchi_version=rdinchi.GetInchiVersion(),
        policy=POLICY_VERSION,
    )
    params = Chem.SmilesParserParams()
    params.parseName = False
    params.allowCXSMILES = False
    with rdBase.BlockLogs():
        try:
            mol = Chem.MolFromSmiles(smiles, params)
        except (ValueError, RuntimeError) as exc:
            result["message"] = str(exc)
            return result
        if mol is None or not mol.GetNumAtoms():
            return result
        if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
            result["status"] = "generic_structure"
            return result
        try:
            inchi, code, message, _, _ = rdinchi.MolToInchi(mol, "")
        except (ValueError, RuntimeError) as exc:
            result.update(status="inchi_error", message=str(exc))
            return result
    result["message"] = message
    if code not in (0, 1) or not inchi.startswith("InChI=1S/"):
        result["status"] = "inchi_error"
        return result
    key = rdinchi.InchiToInchiKey(inchi)
    if len(key) != 27 or key[23:25] != "SA":
        result["status"] = "nonstandard_inchikey"
        return result
    result.update(inchi=inchi, inchikey=key, status="derived")
    return result


@lru_cache(maxsize=4)
def _connection(directory: str, version: str, process: int):
    path = Path(directory) / version / "results.sqlite"
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path, timeout=60)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.execute(
        "CREATE TABLE IF NOT EXISTS results (smiles TEXT PRIMARY KEY, result TEXT NOT NULL) WITHOUT ROWID"
    )
    con.commit()
    return con


@lru_cache(maxsize=16384)
def _cached(smiles: str, directory: str, version: str):
    con = _connection(directory, version, os.getpid())
    row = con.execute("SELECT result FROM results WHERE smiles=?", [smiles]).fetchone()
    if row:
        return json.loads(row[0])
    result = derive(smiles)
    with con:
        con.execute("INSERT OR IGNORE INTO results VALUES (?,?)", [smiles, json.dumps(result)])
    return result


def cached_derivation(smiles: str):
    directory = os.environ.get("OMNIPATH_STRUCTURE_CACHE", "data/reference/.structure-cache")
    return dict(_cached(smiles.strip(), str(Path(directory).resolve()), fingerprint()))
