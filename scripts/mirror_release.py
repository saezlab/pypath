"""Copy one published release into a local data root, checking every checksum.

A PostgreSQL load can read a release over HTTPS, but a large one (171 resources)
failed twice on beauty with DuckDB "Could not resolve hostname" in the projection
staging, after thousands of good requests. Loading from a local copy takes the
network out of the load:

    uv run python scripts/mirror_release.py 2026.10.14.2 /tmp/release \\
        --base-url https://data.omnipath-metabo-dev.schaul.click
    make load-postgres RELEASE_MANIFEST=/tmp/release/releases/2026.10.14.2.json \\
        DATA_ROOT=/tmp/release POSTGRES_SCHEMA=...

The copy has the layout of the data root: the release manifest, each resource's
build manifest and files (checked against the build manifest's SHA-256), and the
taxonomy reference. Files the host does not publish (the serving layer's tables)
answer 404 and are skipped; the loader does not read them. A file already present
with the right checksum is kept, so an interrupted copy can be resumed.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import ssl
import time
import urllib.error
import urllib.request

import certifi

# Standalone Python builds (uv) do not find the system's certificates on every host.
CONTEXT = ssl.create_default_context(cafile=None if os.environ.get("SSL_CERT_FILE") else certifi.where())


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def fetch(base_url: str, target: Path, path: str, sha256: str | None = None) -> int:
    """Copy one file; the number of bytes, 0 for a file the host does not publish."""
    local = target / path
    if local.is_file() and (sha256 is None or _digest(local) == sha256):
        return local.stat().st_size
    local.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, 11):
        partial = local.with_suffix(local.suffix + ".part")
        try:
            with urllib.request.urlopen(f"{base_url}/{path}", timeout=120, context=CONTEXT) as response:
                with partial.open("wb") as out:
                    while block := response.read(1 << 20):
                        out.write(block)
            if sha256 and _digest(partial) != sha256:
                raise ValueError(f"checksum mismatch: {path}")
            partial.replace(local)
            return local.stat().st_size
        except urllib.error.HTTPError as error:
            if error.code == 404:
                print(f"not published: {path}", flush=True)
                return 0
            print(f"retry {attempt} {path}: {error}", flush=True)
        except (OSError, ValueError) as error:
            print(f"retry {attempt} {path}: {error}", flush=True)
        time.sleep(min(60, 2**attempt))
    raise RuntimeError(f"failed: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("release", help="release version, e.g. 2026.10.14.2")
    parser.add_argument("target", type=Path, help="local data root to fill")
    parser.add_argument("--base-url", required=True, help="the data host's root URL")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    base_url = args.base_url.rstrip("/")

    manifest_path = f"releases/{args.release}.json"
    fetch(base_url, args.target, manifest_path)
    manifest = json.loads((args.target / manifest_path).read_text())
    pinned = manifest.get("resource_manifests", {})
    jobs = []
    for source, version in sorted(manifest["resources"].items()):
        directory = f"resources/{source}/{version}"
        fetch(base_url, args.target, f"{directory}/build_manifest.json", pinned.get(source))
        build = json.loads((args.target / directory / "build_manifest.json").read_text())
        jobs += [(f"{directory}/{name}", meta.get("sha256")) for name, meta in build["files"].items()]
    taxonomy = manifest.get("references", {}).get("taxonomy")
    if taxonomy:
        jobs.append((f"references/taxonomy/{taxonomy['version']}/taxonomy.parquet", taxonomy["version"]))

    started, total = time.monotonic(), 0
    with ThreadPoolExecutor(args.workers) as pool:
        sizes = pool.map(lambda job: fetch(base_url, args.target, *job), jobs)
        for number, size in enumerate(sizes, 1):
            total += size
            if number % 200 == 0:
                print(f"{number}/{len(jobs)} files, {total / 1e9:.1f} GB, {time.monotonic() - started:.0f} s", flush=True)
    print(json.dumps({"files": len(jobs), "gigabytes": round(total / 1e9, 2), "seconds": round(time.monotonic() - started)}))


if __name__ == "__main__":
    main()
