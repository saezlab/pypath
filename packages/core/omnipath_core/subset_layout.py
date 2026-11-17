"""Subset products on disk: ``<data_root>/subsets/<release>/<product>/``.

A product is a few Parquet tables and a ``subset_manifest.json`` that pins the
release it was built from. It follows the rules of resources and releases: it is
written in a temporary directory, validated, then moved into place, and it is
never replaced. A subset whose release manifest has changed since the build is
refused when it is read.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Callable

import pyarrow.parquet as pq

SCHEMA_VERSION = 1
MANIFEST_NAME = "subset_manifest.json"
SUBSETS_DIR = "subsets"

# Products in the order they are built. A product reads the ones it depends on.
DEPENDS_ON = {
    "exclusion_metabolites": (),
    "metsigdb": ("exclusion_metabolites",),
}
PRODUCTS = tuple(DEPENDS_ON)

# The version of each product's row contract: the columns of its tables and what they mean.
# A reader that serves or loads a product checks it, so a change of contract is refused
# instead of misread. Raise it with any change a reader could not absorb.
ROW_CONTRACTS = {
    "exclusion_metabolites": 1,
    "metsigdb": 2,
}


class SubsetError(ValueError):
    """A subset cannot be built, read or trusted."""


class SubsetExists(SubsetError):
    """The release already has this product; subsets are never replaced."""


@dataclass(frozen=True)
class BuildContext:
    """What a product writer needs to know about the build it takes part in."""

    data_root: Path
    release: str
    release_manifest_sha256: str
    partial: bool
    # Directory of each product this one depends on.
    dependencies: dict[str, Path] = field(default_factory=dict)


@dataclass(frozen=True)
class Subset:
    """A subset product found on disk and checked against its release."""

    directory: Path
    manifest: dict
    manifest_sha256: str

    @property
    def tables(self) -> dict[str, Path]:
        return {name: self.directory / name for name in self.manifest["tables"]}


def subset_directory(data_root: str | Path, release: str, product: str) -> Path:
    return Path(data_root) / SUBSETS_DIR / release / product


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(16 * 1024**2), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _release(data_root: Path, release: str) -> tuple[str, bool]:
    """The digest of the published manifest and whether the release pins samples."""
    path = data_root / "releases" / f"{release}.json"
    if not path.is_file():
        raise SubsetError(f"release {release} is not published under {data_root / 'releases'}")
    content = path.read_bytes()
    manifest = json.loads(content)
    return hashlib.sha256(content).hexdigest(), manifest.get("partial_resources") is True


def _code_version() -> str:
    try:
        commit = subprocess.run(
            ["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = ""
    return commit or "unknown"


def _described_tables(directory: Path) -> dict[str, dict]:
    tables = {}
    for path in sorted(directory.glob("*.parquet")):
        try:
            rows = pq.read_metadata(path).num_rows
        except Exception as error:
            raise SubsetError(f"{path.name} is not readable Parquet: {error}") from error
        tables[path.name] = {"rows": rows, "sha256": _sha256(path)}
    if not tables:
        raise SubsetError("the product wrote no table")
    return tables


def read_subset(data_root: str | Path, release: str, product: str) -> Subset | None:
    """The product of a release, checked against the release and its own files.

    Returns None when the release has no such subset. A subset that does not
    match its release manifest or its recorded tables raises ``SubsetError``.
    """
    root = Path(data_root)
    directory = subset_directory(root, release, product)
    manifest_path = directory / MANIFEST_NAME
    if not manifest_path.is_file():
        return None
    manifest = json.loads(manifest_path.read_text())
    release_digest, _ = _release(root, release)
    if manifest["release_manifest_sha256"] != release_digest:
        raise SubsetError(
            f"{product} of release {release} was built from another release manifest "
            f"({manifest['release_manifest_sha256'][:12]}, now {release_digest[:12]})"
        )
    for name, recorded in manifest["tables"].items():
        path = directory / name
        if not path.is_file() or _sha256(path) != recorded["sha256"]:
            raise SubsetError(f"{product} of release {release}: {name} differs from its manifest")
    return Subset(directory, manifest, _sha256(manifest_path))


def present_subsets(data_root: str | Path, release: str) -> list[Subset]:
    """The subsets a release has, in build order."""
    found = (read_subset(data_root, release, product) for product in PRODUCTS)
    return [subset for subset in found if subset is not None]


def build_product(
    data_root: str | Path,
    release: str,
    product: str,
    writer: Callable[[Path, BuildContext], dict | None],
) -> dict:
    """Run ``writer`` in a temporary directory and publish what it wrote.

    ``writer(directory, context)`` writes the product's Parquet tables into
    ``directory`` and may return counts to keep in the manifest. The tables are
    checked, the manifest is added and the directory is moved to
    ``subsets/<release>/<product>/``. Returns the manifest.
    """
    root = Path(data_root)
    if product not in DEPENDS_ON:
        raise SubsetError(f"unknown product {product}; the products are {', '.join(PRODUCTS)}")
    release_digest, partial = _release(root, release)
    final = subset_directory(root, release, product)
    if final.exists():
        raise SubsetExists(f"{product} of release {release} exists; subsets are not replaced")

    dependencies, depends_on = {}, {}
    for name in DEPENDS_ON[product]:
        found = read_subset(root, release, name)
        if found is None:
            raise SubsetError(f"{product} needs {name} of release {release}; build it first")
        dependencies[name] = found.directory
        depends_on[name] = found.manifest_sha256

    context = BuildContext(root, release, release_digest, partial, dependencies)
    parent = final.parent
    parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=f".{product}-", dir=parent))
    try:
        stats = writer(work, context)
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "product": product,
            "release": release,
            "release_manifest_sha256": release_digest,
            "partial": partial,
            "depends_on": depends_on,
            "code_version": _code_version(),
            "row_contract": ROW_CONTRACTS[product],
            "stats": stats if isinstance(stats, dict) else {},
            "tables": _described_tables(work),
            "built_at": datetime.now(timezone.utc).isoformat(),
        }
        (work / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2) + "\n")
        try:
            work.rename(final)
        except OSError as error:
            if final.exists():
                raise SubsetExists(
                    f"{product} of release {release} appeared while it was built"
                ) from error
            raise
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return manifest
