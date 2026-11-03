"""Pin an immutable identity library for runtime matching."""

from pathlib import Path


def pin_library(library_dir):
    """The resolved library path, so a later symlink change cannot switch it mid-build."""
    if library_dir is None:
        return None
    return Path(library_dir).resolve()
