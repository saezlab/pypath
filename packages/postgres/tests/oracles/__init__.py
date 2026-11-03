"""Immutable tracked scientific oracles; candidate code never supplies fixtures."""

import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
LEGACY_ROOT = REPO / "legacy/postgres/omnipath_build"
CV_ROOT = Path(__file__).resolve().parent.parent / "oracle_cv_terms"
SOURCE_HASHES = json.loads((Path(__file__).parent / "source_sha256.json").read_text())


def read_checked(path, sha256):
    data = Path(path).read_bytes()
    actual = hashlib.sha256(data).hexdigest()
    if actual != sha256:
        raise AssertionError(f"Immutable oracle changed: {path} ({actual})")
    return data.decode("utf-8")


def read_legacy(relative, sha256=None):
    digest = SOURCE_HASHES["legacy"][str(relative)]
    if sha256 is not None and digest != sha256:
        raise AssertionError(f"Oracle digest declaration disagrees: {relative}")
    return read_checked(LEGACY_ROOT / relative, digest)


def read_frozen(path):
    path = Path(path).resolve()
    if path.is_relative_to(LEGACY_ROOT):
        return read_legacy(path.relative_to(LEGACY_ROOT))
    return read_checked(path, SOURCE_HASHES["cv"][str(path.relative_to(CV_ROOT))])


def read_deployed_cosmos():
    return read_checked(
        Path(__file__).parent / "cosmos_project_edges_951d500.sql",
        "987c6870938d52e1dd5ac6af1837f30ca273e3a1b0fd4f7d6032bf2d8b734071",
    )
