"""
Parser of the Human Signaling Network file of the Wang Lab.

The file is comma separated and holds three tables one after the other, each
announced by a line of its own and separated by empty lines:

    Notes:,,,,,
    A-A,Anti-Apoptic,,,,
    Cy,Intracellular,,,
    ...
    +,activation,,,
    _,inhibition,,,
    0,physical interaction,,,

    Node Number in the  network,Gene_Name,Entrez Gene ID,Molecular Function,Cellular Location
    10,CIAPIN1,57019,NA,NA
    11,CREBBP,1387,Adapter,Nu
    11,CREBBP,1387,Kinase,Cm

    From (Node Number),To (Node Number),Link_Type,,
    1453,1226,+

The notes are the key of the abbreviations and of the link types. A node has
one line for each of its function and location pairs. A link names its two
nodes by number. A physical interaction (``0``), which has no direction, is
often listed in both orientations (``A,B,0`` and ``B,A,0``); a link with a
sign is a statement of its own and is never repeated this way.
"""

from __future__ import annotations

from collections.abc import Generator
import csv
import io
from typing import Any

from pypath.inputs_v2.parsers.base import read_bytes

NOTES_HEADER = 'Notes:'
NODE_HEADER = (
    'Node Number in the  network',
    'Gene_Name',
    'Entrez Gene ID',
    'Molecular Function',
    'Cellular Location',
)
LINK_HEADER = ('From (Node Number)', 'To (Node Number)', 'Link_Type')
# The link types as the notes of the file define them. The links themselves
# write the inhibition as ``-``.
# The link type of an undirected link.
UNDIRECTED = '0'
LINK_TYPE_KEY = {
    '+': 'activation',
    '_': 'inhibition',
    '0': 'physical interaction',
}


def _append(values: list[str], value: str) -> None:
    if value and value not in values:
        values.append(value)


def iter_network(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the links of the network, each with the records of its two nodes.

    Output fields, for ``source`` and ``target``:

    - ``<role>_node``: the node number
    - ``<role>_symbol``, ``<role>_entrez``: ``Gene_Name`` and ``Entrez Gene ID``
    - ``<role>_functions``, ``<role>_locations``: the distinct ``Molecular
      Function`` and ``Cellular Location`` values of the node's lines, in
      source order and as written (abbreviations are not expanded)

    and ``link_type``, the ``Link_Type`` as written, and ``mirror``: True on
    a physical interaction (``0``) whose reverse, the same two nodes the other
    way round, is an earlier link of the file. The two carry the same
    information: the file gives a link nothing but its two nodes and its type,
    and the node records are the same in both.

    Raises ValueError if a table is missing, the notes do not define the link
    types as expected, a node number stands for two genes, or a link names a
    node the node table does not have.
    """
    text = read_bytes(opener).decode('utf-8')
    notes: dict[str, str] = {}
    nodes: dict[str, dict[str, Any]] = {}
    section = None
    count = 0
    undirected: set[tuple[str, str]] = set()
    for fields in csv.reader(io.StringIO(text)):
        fields = [field.strip() for field in fields]
        if not any(fields):
            continue
        if fields[0] == NOTES_HEADER:
            section = 'notes'
        elif tuple(fields[:5]) == NODE_HEADER:
            section = 'nodes'
        elif tuple(fields[:3]) == LINK_HEADER:
            if not nodes:
                raise ValueError('Wang network file has no node table')
            for code, meaning in LINK_TYPE_KEY.items():
                if notes.get(code) != meaning:
                    raise ValueError(
                        f'Wang network file defines link type {code!r} as '
                        f'{notes.get(code)!r}, expected {meaning!r}'
                    )
            section = 'links'
        elif section == 'notes' and len(fields) > 1 and fields[1]:
            notes[fields[0]] = fields[1]
        elif section == 'nodes':
            if len(fields) < 5:
                raise ValueError(f'Malformed Wang network node: {fields}')
            number, symbol, entrez, function, location = fields[:5]
            if not number.isdigit():
                raise ValueError(f'Malformed Wang network node: {fields}')
            node = nodes.setdefault(
                number,
                {'symbol': symbol, 'entrez': entrez, 'functions': [], 'locations': []},
            )
            if (node['symbol'], node['entrez']) != (symbol, entrez):
                raise ValueError(
                    f'Wang network node {number} stands for two genes: '
                    f'{node["symbol"]} and {symbol}'
                )
            _append(node['functions'], function)
            _append(node['locations'], location)
        elif section == 'links':
            if len(fields) < 3:
                raise ValueError(f'Malformed Wang network link: {fields}')
            record: dict[str, Any] = {}
            for role, number in zip(('source', 'target'), fields[:2], strict=True):
                if number not in nodes:
                    raise ValueError(f'Wang network link to unknown node {number}')
                record[f'{role}_node'] = number
                for field, value in nodes[number].items():
                    record[f'{role}_{field}'] = (
                        list(value) if isinstance(value, list) else value
                    )
            record['link_type'] = fields[2]
            pair = (record['source_node'], record['target_node'])
            record['mirror'] = (
                fields[2] == UNDIRECTED
                and pair[0] != pair[1]
                and pair[::-1] in undirected
            )
            if fields[2] == UNDIRECTED:
                undirected.add(pair)
            count += 1
            yield record
    if not count:
        raise ValueError('Wang network file has no links')
