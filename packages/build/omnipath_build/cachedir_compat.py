"""Work around a cachedir UnboundLocalError on files with no explicit extension.

``Opener.set_type`` does ``ext = self.ext or ext(self.path)``, which shadows the
imported ``ext`` helper. SIGNOR (and similar) downloads omit ``Download.ext``
and then crash before any rows are read.
"""

from __future__ import annotations


def patch_cachedir_opener() -> None:
    try:
        from cachedir._open import ARCHIVES, COMPRESSED, Opener
        from pkg_infra.utils._misc import ext as path_ext
    except Exception:
        return

    if getattr(Opener.set_type, "_omnipath_patched", False):
        return

    def set_type(self) -> None:
        suffix = self.ext or path_ext(self.path) or "txt"
        suffix = str(suffix).strip(".")
        suffix = "tar.gz" if suffix == "tgz" else suffix
        self.ext = suffix
        self.type = suffix if suffix in COMPRESSED | ARCHIVES else "plain"
        self.type = "tar" if self.type.startswith("tar") else self.type

    set_type._omnipath_patched = True
    Opener.set_type = set_type
