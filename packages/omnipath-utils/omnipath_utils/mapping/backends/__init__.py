"""Mapping backend registry."""

from __future__ import annotations

from typing import TYPE_CHECKING
import logging

if TYPE_CHECKING:
    from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)

_BACKENDS: dict[str, type[MappingBackend]] = {}


def register(name: str, backend_cls: type[MappingBackend]) -> None:
    """Register a mapping backend."""

    _BACKENDS[name] = backend_cls


def get_backend(name: str) -> MappingBackend | None:
    """Get a backend instance by name. Lazy-loads backend modules."""

    if name not in _BACKENDS:
        _try_import(name)

    cls = _BACKENDS.get(name)

    return cls() if cls else None


def _try_import(name: str) -> None:
    """Try to import a backend module to trigger registration."""

    try:
        import importlib

        importlib.import_module(
            f'omnipath_utils.mapping.backends._{name}',
        )
    except ImportError:
        _log.debug('Backend %s not available', name)
