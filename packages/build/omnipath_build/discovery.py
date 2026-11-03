"""Source and dataset discovery for the inputs_v2 ecosystem.

Walks available inputs_v2 packages, extracts ResourceConfig metadata, and identifies
callable dataset objects that produce Silver Entity and Relation observations.
"""

from __future__ import annotations

import importlib
import inspect
import json
import logging
import os
import pkgutil
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any


logger = logging.getLogger(__name__)


@dataclass(slots=True)
class DiscoveredDataset:
    """A discovered dataset callable with metadata."""

    source: str
    dataset_name: str
    qualified_module: str
    call: Callable[..., Iterable[Any]]
    output_kind: str = "entity"
    raw_dataset: Any = None
    mapper: Callable[[Any], Any] | None = None


class DiscoveryError(RuntimeError):
    """Raised when resource or dataset discovery fails."""


class UnknownResourceError(DiscoveryError):
    """Raised when the requested source or dataset is not found."""


SOURCE_OVERRIDES: dict[str, str] = {
    "connectomedb": "connectomedb2025",
}


def canonical_json(value: Any) -> str:
    """Format any nested value as deterministic canonical JSON."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def resolve_source_name(module_name: str) -> str:
    """Map inputs_v2 module name to standard source slug."""
    return SOURCE_OVERRIDES.get(module_name.lower(), module_name.lower())


def setup_pypath_cache(cache_dir: str | Path | None = None) -> Path:
    """Ensure pypath uses project-local or specified cache directory."""
    if cache_dir:
        path = Path(cache_dir).resolve()
    elif "PYPATH_DOWNLOAD_DATADIR" in os.environ:
        path = Path(os.environ["PYPATH_DOWNLOAD_DATADIR"]).resolve()
    else:
        path = Path.cwd() / "data" / "pypath-data"
    os.environ["PYPATH_DOWNLOAD_DATADIR"] = str(path)
    path.mkdir(parents=True, exist_ok=True)
    from pypath.share import settings

    settings.setup(progressbars=False, cachedir=str(path))
    return path


def extract_resource_config(module_or_name: Any) -> dict[str, Any]:
    """Serialize the discovered ResourceConfig into a clean dictionary."""
    try:
        if isinstance(module_or_name, str):
            module = importlib.import_module(module_or_name)
        else:
            module = module_or_name
        from pypath.inputs_v2.base import Resource

        resources = [v for v in vars(module).values() if isinstance(v, Resource)]
        if not resources:
            return {}
        config = resources[0].config
        cfg_dict = {k: v for k, v in vars(config).items() if not k.startswith("_")}
        return json.loads(canonical_json(cfg_dict))
    except (ImportError, AttributeError, TypeError, ValueError) as exc:
        raise DiscoveryError(
            f"Unable to load resource metadata from {module_or_name}: {exc}"
        ) from exc


def _collect_sources(
    *,
    inputs_package: str,
    cache_dir: str | Path | None,
    source_filter: str | None = None,
) -> tuple[dict[str, list[DiscoveredDataset]], dict[str, dict[str, Any]]]:
    setup_pypath_cache(cache_dir)

    try:
        root_module = importlib.import_module(inputs_package)
    except ImportError as exc:
        raise DiscoveryError(f"Unable to import inputs package '{inputs_package}': {exc}") from exc

    package_paths = getattr(root_module, "__path__", None)
    if package_paths is None:
        raise DiscoveryError(f"Inputs package '{inputs_package}' is not a namespace package")

    prefix = f"{inputs_package}."
    all_sources: dict[str, list[DiscoveredDataset]] = {}
    config_by_source: dict[str, dict[str, Any]] = {}

    from pypath.inputs_v2.base import ArtifactDataset, Dataset, Resource

    dataset_types = (Dataset, ArtifactDataset)

    for module_info in pkgutil.walk_packages(package_paths, prefix):
        module_name = module_info.name
        rel_name = module_name[len(prefix) :]
        if not rel_name or rel_name.split(".")[-1].startswith("_"):
            continue
        slug = resolve_source_name(rel_name)
        if source_filter is not None and slug != source_filter:
            continue

        try:
            mod = importlib.import_module(module_name)
        except Exception as exc:
            if source_filter is not None:
                raise DiscoveryError(
                    f"Unable to load requested source {module_name}: {exc}"
                ) from exc
            logger.warning("Skipping broken source module %s: %s", module_name, exc, exc_info=True)
            continue

        resource_objs = [obj for _, obj in inspect.getmembers(mod) if isinstance(obj, Resource)]
        if resource_objs and slug not in config_by_source:
            config_by_source[slug] = extract_resource_config(mod)

        found_datasets: list[tuple[str, Any]] = [
            (name, obj)
            for name, obj in inspect.getmembers(mod)
            if isinstance(obj, dataset_types) and getattr(obj, "kind", None) != "id_translation"
        ]

        seen = {name for name, _ in found_datasets}
        for res in resource_objs:
            for ds_name, ds_obj in res.datasets().items():
                if getattr(ds_obj, "kind", None) != "id_translation" and ds_name not in seen:
                    found_datasets.append((ds_name, ds_obj))
                    seen.add(ds_name)

        callables: list[DiscoveredDataset] = []
        for ds_name, ds_obj in found_datasets:
            raw_ds = (
                ds_obj if isinstance(ds_obj, Dataset) else getattr(ds_obj, "_raw_dataset", None)
            )
            mapper = getattr(raw_ds, "mapper", None) if raw_ds is not None else None
            callables.append(
                DiscoveredDataset(
                    source=slug,
                    dataset_name=ds_name,
                    qualified_module=module_name,
                    call=ds_obj,
                    raw_dataset=raw_ds,
                    mapper=mapper,
                )
            )

        if callables:
            all_sources.setdefault(slug, []).extend(callables)

    return all_sources, config_by_source


def list_sources(
    *,
    inputs_package: str = "pypath.inputs_v2",
    cache_dir: str | Path | None = None,
) -> list[dict[str, Any]]:
    """List discoverable inputs_v2 sources and their dataset names."""
    all_sources, _config = _collect_sources(inputs_package=inputs_package, cache_dir=cache_dir)
    return [
        {
            "source": slug,
            "datasets": [item.dataset_name for item in datasets],
        }
        for slug, datasets in sorted(all_sources.items())
    ]


def discover_datasets(
    source: str,
    *,
    inputs_package: str = "pypath.inputs_v2",
    datasets: Sequence[str] | None = None,
    cache_dir: str | Path | None = None,
) -> tuple[str, list[DiscoveredDataset], dict[str, Any]]:
    """Discover all executable datasets and configuration for a given source."""
    source_lower = source.lower()
    all_sources, config_by_source = _collect_sources(
        inputs_package=inputs_package,
        cache_dir=cache_dir,
        source_filter=source_lower,
    )

    if source_lower not in all_sources:
        # Keep the complete available-source diagnostic for invalid requests.
        all_sources, config_by_source = _collect_sources(
            inputs_package=inputs_package,
            cache_dir=cache_dir,
        )
        avail = ", ".join(sorted(all_sources.keys()))
        raise UnknownResourceError(f"Unknown inputs_v2 source '{source}'; available: {avail}")

    selected = all_sources[source_lower]
    if datasets:
        requested = set(datasets)
        selected = [d for d in selected if d.dataset_name in requested]
        found = {d.dataset_name for d in selected}
        missing = sorted(requested - found)
        if missing:
            raise UnknownResourceError(f"Unknown dataset(s) for {source}: {', '.join(missing)}")

    if not selected:
        raise UnknownResourceError(f"No executable datasets found for source '{source}'")

    config = config_by_source.get(source_lower, {})
    return source_lower, selected, config
