"""
Read the edge attribute files of the Laudanna Lab.

The files are Cytoscape edge attribute files (``.EA``): a title line, then one
attribute of one edge per line, the edge named by its two nodes and the
interaction type ``pp`` (protein-protein):

    Signaling_flow
    F2R (pp) PAWR = state_change

    SignalingDirection
    EDNRA (pp) NBN = activation

The nodes are HGNC gene symbols. Four lines of the signaling flow file lack the
second node (``PKMYT1 (pp)  = state_change``); they are yielded with a target
of None. An edge can have several lines.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import read_text

# The title line of each file and the values its edges can have.
SIGNALING_FLOW = 'Signaling_flow'
SIGNALING_DIRECTION = 'SignalingDirection'
FLOW_VALUES = frozenset({'state_change', 'complex_formation'})
DIRECTION_VALUES = frozenset({'activation', 'inhibition', 'docking'})

_LINE = re.compile(r'(\S+) \(pp\) (\S*) = (\S+)')


def iter_edges(
    opener,
    title: str,
    values: frozenset[str],
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield the lines of an edge attribute file.

    ``title`` is the first line the file must have and ``values`` the
    attribute values it may have; a file with another title, a line that is not
    an edge attribute, a value outside ``values`` and a file without edges
    raise ``ValueError``.

    Keys: ``source`` and ``target`` (gene symbols, ``target`` None where the
    line has none) and ``value`` (the attribute, as written).
    """
    lines = read_text(opener).splitlines()
    if not lines or lines[0].strip() != title:
        raise ValueError(
            f'Laudanna file does not start with the title {title!r}: {lines[:1]}'
        )
    count = 0
    for line in lines[1:]:
        if not line.strip():
            continue
        match = _LINE.fullmatch(line.strip())
        if not match:
            raise ValueError(f'Laudanna file has a malformed line: {line!r}')
        source, target, value = match.groups()
        if value not in values:
            raise ValueError(f'Unknown value in the {title} file: {value!r}')
        count += 1
        yield {'source': source, 'target': target or None, 'value': value}
    if not count:
        raise ValueError(f'Laudanna file {title} has no edges')
