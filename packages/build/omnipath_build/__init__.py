"""OmniPath Full-Parquet Build Pipeline.

This package provides a small, source-neutral pipeline that ingests datasets
from inputs_v2, resolves canonical entity identities and writes each resource version
as normalized Parquet tables (``omnipath_core.schema.PUBLISHED_TABLES``) plus derived
serving lookup tables (``SERVING_TABLES``).
"""

from __future__ import annotations

from typing import Any

__all__ = [
    "build_resource",
    "build_all",
    "orchestrate",
    "Budget",
    "ParquetWriter",
    "EntityResolver",
    "SilverExtractor",
    "discover_datasets",
    "list_sources",
]


def __getattr__(name: str) -> Any:
    if name in {"build_resource", "build_all"}:
        from .pipeline import build_all, build_resource

        return build_resource if name == "build_resource" else build_all
    if name in {"orchestrate", "Budget"}:
        from .orchestrator import orchestrate, Budget

        return orchestrate if name == "orchestrate" else Budget
    if name == "ParquetWriter":
        from .writer import ParquetWriter

        return ParquetWriter
    if name == "EntityResolver":
        from omnipath_resolver import EntityResolver

        return EntityResolver
    if name == "SilverExtractor":
        from .silver import SilverExtractor

        return SilverExtractor
    if name in {"discover_datasets", "list_sources"}:
        from .discovery import discover_datasets, list_sources

        return discover_datasets if name == "discover_datasets" else list_sources
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
