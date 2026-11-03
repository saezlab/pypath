"""Build operations interface for background jobs."""

from __future__ import annotations

import importlib
from typing import Any


def _import_build() -> Any:
    try:
        return importlib.import_module("omnipath_build")
    except ModuleNotFoundError as exc:
        if exc.name != "omnipath_build":
            raise
        return None


class BuildUnavailable(RuntimeError):
    """Raised when omnipath-build is not importable."""


class BuildOps:
    """Lazy wrappers around omnipath-build so serving tests do not import pypath."""

    def hubs_available(self) -> bool:
        if _import_build() is None:
            return False
        try:
            from omnipath_build.hubs.export import export_hubs  # noqa: F401

            return True
        except ModuleNotFoundError as exc:
            if exc.name and (exc.name == "pypath" or exc.name.startswith("pypath.")):
                return False
            raise

    def pipeline_available(self) -> bool:
        if _import_build() is None:
            return False
        try:
            from omnipath_build.pipeline import build_resource  # noqa: F401
            from omnipath_build.silver import SilverExtractor  # noqa: F401

            return True
        except ModuleNotFoundError as exc:
            if exc.name and (exc.name == "pypath" or exc.name.startswith("pypath.")):
                return False
            raise

    def available(self) -> bool:
        return self.hubs_available() or self.pipeline_available()

    def list_sources(self) -> list[dict[str, Any]]:
        if _import_build() is None:
            raise BuildUnavailable("omnipath-build is not installed")
        from omnipath_build.discovery import list_sources

        return list_sources()

    def export_hubs(self, **kwargs: Any) -> dict[str, int]:
        if _import_build() is None:
            raise BuildUnavailable("omnipath-build is not installed")
        from omnipath_build.hubs.export import export_hubs

        return export_hubs(**kwargs)

    def build_all(self, **kwargs: Any) -> dict[str, Any]:
        if _import_build() is None:
            raise BuildUnavailable("omnipath-build is not installed")
        from omnipath_build.cachedir_compat import patch_cachedir_opener
        from omnipath_build.pipeline import build_all

        patch_cachedir_opener()
        return build_all(**kwargs)
