"""
Parser for NicheNet v2 ligand-receptor RDS files.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any
import warnings

import rdata


def _first_handle(opener) -> Any | None:
    if not opener or not opener.result:
        return None
    if isinstance(opener.result, dict):
        return next(iter(opener.result.values()), None)
    if hasattr(opener.result, '__next__'):
        return opener.result
    return opener.result


def _read_bytes(opener) -> bytes:
    handle = _first_handle(opener)
    if handle is None:
        return b''
    if hasattr(handle, 'read'):
        data = handle.read()
        return data if isinstance(data, bytes) else data.encode('utf-8')
    return b''.join(
        chunk if isinstance(chunk, bytes) else str(chunk).encode('utf-8')
        for chunk in handle
    )


def _clean(value: Any) -> str:
    if value is None:
        return ''
    if hasattr(value, 'item'):
        value = value.item()
    text = str(value).strip()
    return '' if text.lower() in {'nan', 'none', 'na'} else text


def iter_lr_network(
    opener,
    *,
    taxon_id: str,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield NicheNet v2 ligand-receptor records."""
    data = _read_bytes(opener)
    if not data:
        return
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message='Unknown file type: assumed RDS',
            category=UserWarning,
        )
        parsed = rdata.parser.parse_data(data)
        warnings.filterwarnings(
            'ignore',
            message='Missing constructor for R class "tbl.*',
            category=UserWarning,
        )
        network = rdata.conversion.convert(parsed)
    if network is None:
        return

    for _, row in network.iterrows():
        ligand = _clean(row.get('from'))
        receptor = _clean(row.get('to'))
        if not ligand or not receptor:
            continue
        yield {
            'ligand': ligand,
            'receptor': receptor,
            'database': _clean(row.get('database')),
            'source': _clean(row.get('source')),
            'taxon_id': str(taxon_id),
        }
