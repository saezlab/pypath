"""Fingerprint the files actually importable in checkouts and wheel installations.

No build hook, Git checkout, resolver call or package import is required. Relative
file names and content hashes make the digest independent of installation paths.
File hashes are cached by filesystem identity/change metadata, not package version.
"""

from __future__ import annotations

from functools import lru_cache
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import subprocess
from typing import Any

_SOFTWARE = {
    "omnipath-build": "omnipath_build",
    "omnipath-core": "omnipath_core",
    "omnipath-resolver": "omnipath_resolver",
    "pypath-omnipath": "pypath",
    "biolink-model": "biolink_model",
}
_NATIVE_SUFFIXES = {".so", ".pyd", ".dll", ".dylib"}
_IGNORED_PARTS = {"__pycache__", "tests", ".git", ".pytest_cache", ".ruff_cache"}


def _package_roots(module: str) -> tuple[Path, ...]:
    spec = importlib.util.find_spec(module)
    if spec is None:
        raise RuntimeError(f"Cannot fingerprint required runtime package {module!r}")
    if spec.submodule_search_locations:
        return tuple(Path(path).resolve() for path in spec.submodule_search_locations)
    if spec.origin and spec.origin not in {"built-in", "frozen"}:
        return (Path(spec.origin).resolve().parent,)
    raise RuntimeError(f"Runtime package {module!r} has no inspectable installed files")


@lru_cache(maxsize=8192)
def _content_hash(path: str, size: int, mtime_ns: int, ctime_ns: int, inode: int) -> str:
    # Change metadata is deliberately part of the key: same-version editable
    # package edits must not keep a stale digest, including restored mtime values.
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fingerprint_files(roots: tuple[Path, ...]) -> dict[str, dict[str, Any]]:
    result = {}
    for root in roots:
        for path in sorted(root.rglob("*")):
            relative = path.relative_to(root)
            if (
                not path.is_file()
                or _IGNORED_PARTS.intersection(relative.parts)
                or path.suffix in {".pyc", ".pyo"}
            ):
                continue
            stat = path.stat()
            fingerprint = {
                "sha256": _content_hash(
                    str(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino
                ),
                "size_bytes": stat.st_size,
            }
            key = relative.as_posix()
            if key in result and result[key] != fingerprint:
                raise RuntimeError(f"Conflicting installed namespace-package file {key!r}")
            result[key] = fingerprint
    if not result:
        raise RuntimeError("Required runtime package contains no inspectable files")
    return result


def _aggregate(files: dict[str, dict[str, Any]]) -> str:
    encoded = json.dumps(files, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _git_provenance(root: Path) -> dict[str, Any] | None:
    """Optional checkout metadata, never inferred from an unrelated ancestor repo."""
    try:
        repo = Path(
            subprocess.run(
                ["git", "rev-parse", "--show-toplevel"],
                cwd=root,
                capture_output=True,
                text=True,
                check=True,
                timeout=5,
            ).stdout.strip()
        ).resolve()
        relative = root.relative_to(repo).as_posix()
        tracked = subprocess.run(
            ["git", "ls-files", "--", relative],
            cwd=repo,
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        ).stdout
        if not any(line.startswith(relative + "/") for line in tracked.splitlines()):
            return None
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo,
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain", "--", relative],
                cwd=repo,
                capture_output=True,
                text=True,
                check=True,
                timeout=5,
            ).stdout.strip()
        )
        return {"commit": revision, "dirty": dirty}
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


def software_provenance() -> dict[str, Any]:
    """Exact installed code/native/vocabulary content, including same-version edits.

    The content digest is authoritative. Optional Git metadata describes the
    checkout and its dirty state; it never substitutes for an installed-file hash.
    """
    packages = {}
    for distribution, module in _SOFTWARE.items():
        roots = _package_roots(module)
        files = _fingerprint_files(roots)
        try:
            version = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            version = None  # source-tree imports still have an exact content digest
        package = {
            "version": version,
            "content_sha256": _aggregate(files),
            "file_count": len(files),
            "native_files": {
                name: value
                for name, value in files.items()
                if Path(name).suffix in _NATIVE_SUFFIXES
            },
            "git": _git_provenance(roots[0]) if len(roots) == 1 else None,
        }
        if module == "omnipath_resolver" and not package["native_files"]:
            raise RuntimeError("Cannot fingerprint the required native resolver extension")
        if module == "omnipath_core":
            vocabulary = {name: value for name, value in files.items() if name.startswith("vocab/")}
        elif module == "pypath":
            vocabulary = {
                name: value
                for name, value in files.items()
                if name.startswith("internals/cv_terms/")
            }
        elif module == "biolink_model":
            vocabulary = {
                name: value for name, value in files.items() if name.startswith("schema/")
            }
        else:
            vocabulary = {}
        if vocabulary:
            package["vocabulary"] = {
                "content_sha256": _aggregate(vocabulary),
                "file_count": len(vocabulary),
            }
        packages[distribution] = package
    identity = {
        name: {"version": package["version"], "content_sha256": package["content_sha256"]}
        for name, package in packages.items()
    }
    return {
        "format": "installed-package-content-v1",
        "sha256": hashlib.sha256(
            json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "packages": packages,
    }
