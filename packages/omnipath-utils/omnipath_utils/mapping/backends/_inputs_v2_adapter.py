"""Generic ``inputs_v2`` adapter backend (Milestone I).

Consumes a pypath ``inputs_v2`` resource's raw structure-bearing rows
(``resource.<dataset>.raw()`` → dict rows like ``chebi_id``/``inchikey``/
``inchi``/``smiles``) into the ``dict[str, set]`` mapping contract — without
building Entity records (no ``Dataset.__call__`` overhead) and without RDKit
(structure strings come straight from the source; Principle II keeps
cheminformatics out of utils).

Subclass :class:`InputsV2Backend`, set ``resource_module`` (the
``pypath.inputs_v2`` module) and ``dataset`` (e.g. ``molecules``/``lipids``), and
``register(...)``. The canonical id_type → raw-column mapping comes from
``id_types.yaml`` (the base class resolves ``src_col``/``tgt_col``), so a backend
is a few lines and ``id_types.yaml`` wiring.
"""

from __future__ import annotations

import importlib
import logging
from itertools import islice
from typing import Any, Iterable

from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)

# Parsed raw rows are cached per (resource, dataset) for the process — a build
# reads several id-type pairs from the same resource and the parse is expensive.
_RAW_CACHE: dict[tuple[str, str], list[dict[str, Any]]] = {}


def _as_values(value: Any) -> list[str]:
    """Normalise a raw cell to a list of non-empty strings (cells may be lists)."""
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(v).strip() for v in value if v is not None and str(v).strip()]
    text = str(value).strip()
    return [text] if text else []


def build_mapping(
    rows: Iterable[dict[str, Any]],
    src_col: str,
    tgt_col: str,
) -> dict[str, set[str]]:
    """Project raw rows into ``{src_value: {tgt_value, ...}}`` (pure, testable)."""
    result: dict[str, set[str]] = {}
    for row in rows:
        sources = _as_values(row.get(src_col))
        targets = _as_values(row.get(tgt_col))
        if not sources or not targets:
            continue
        for source in sources:
            result.setdefault(source, set()).update(targets)
    return result


def raw_rows(
    resource_module: str,
    dataset: str,
    limit: int | None = None,
) -> list[dict[str, Any]]:
    """Materialise (and cache) an inputs_v2 dataset's raw rows. Needs pypath."""
    cache_key = (resource_module, dataset)
    if limit is None and cache_key in _RAW_CACHE:
        return _RAW_CACHE[cache_key]
    module = importlib.import_module(f'pypath.inputs_v2.{resource_module}')
    dataset_obj = getattr(module.resource, dataset)
    iterator = dataset_obj.raw()
    if limit is not None:
        return list(islice(iterator, limit))
    rows = list(iterator)
    _RAW_CACHE[cache_key] = rows
    return rows


class InputsV2Backend(MappingBackend):
    """Map between columns of an ``inputs_v2`` resource's raw rows."""

    #: The ``pypath.inputs_v2`` module name (e.g. ``chebi``).
    resource_module: str = ''
    #: The dataset attribute on ``resource`` whose ``.raw()`` is read.
    dataset: str = ''

    def _read_via_pypath(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        *,
        src_col: str,
        tgt_col: str,
        limit: int | None = None,
        **kwargs: object,
    ) -> dict[str, set[str]]:
        importlib.import_module('pypath')  # ImportError → base falls to _read_direct
        rows = raw_rows(self.resource_module, self.dataset, limit)
        return build_mapping(rows, src_col, tgt_col)

    def _read_direct(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        *,
        src_col: str,
        tgt_col: str,
        **kwargs: object,
    ) -> dict[str, set[str]]:
        raise ImportError(f'{self.name} requires pypath inputs_v2')
