"""
Parser for NeuronChat .rda files.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

import rdata


def _first_handle(opener) -> Any | None:
    """Extract the first file handle from an opener result."""
    if not opener or not opener.result:
        return None
    # Opener.result can be a handle, a dict, or a generator/iterator
    if isinstance(opener.result, dict):
        return next(iter(opener.result.values()), None)
    
    # If it's a generator/iterator (which it seems to be from the error)
    if hasattr(opener.result, '__next__'):
        return opener.result

    return opener.result


def iter_neuronchat(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    Parse NeuronChat .rda files and yield interaction records.

    Args:
        opener: The opener object for handling binary R data files.

    Yields:
        Dictionary representing a single interaction record.
    """
    handle = _first_handle(opener)
    if not handle:
        return

    # rdata expects a file-like object or bytes.
    # If handle is an iterator (common in Opener), we may need to read it all.
    # Opener usually wraps the file in a way that iterating yields lines.
    # For binary files, we want the whole blob.
    
    try:
        if hasattr(handle, 'read'):
            data_bytes = handle.read()
        else:
            # Fallback if it's an iterator of chunks/lines
            data_bytes = b''.join(handle)
            
        parsed = rdata.parser.parse_data(data_bytes)
        converted = rdata.conversion.convert(parsed)
    except Exception as e:
        print(f"Error parsing .rda data: {e}")
        return

    if not converted:
        return

    db_key = list(converted.keys())[0]
    db = converted[db_key]

    if not isinstance(db, dict):
        return

    for interaction_name, data in db.items():
        record = {
            'interaction_name': _to_str(data.get('interaction_name')),
            'lig_contributor': _to_list(data.get('lig_contributor')),
            'receptor_subunit': _to_list(data.get('receptor_subunit')),
            'ligand_type': _to_str(data.get('ligand_type')),
            'interaction_type': _to_str(data.get('interaction_type')),
        }
        yield record


def _to_str(value: Any) -> str:
    """Helper to convert NumPy array or single value to string."""
    if value is None:
        return ''
    if hasattr(value, 'tolist'): # NumPy array
        value = value.tolist()
    if isinstance(value, list) and len(value) > 0:
        return str(value[0])
    return str(value)


def _to_list(value: Any) -> list[str]:
    """Helper to convert NumPy array or single value to list of strings."""
    if value is None:
        return []
    if hasattr(value, 'tolist'): # NumPy array
        value = value.tolist()
    if not isinstance(value, list):
        value = [value]
    return [str(v) for v in value if v]
