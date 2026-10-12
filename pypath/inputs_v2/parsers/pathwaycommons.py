"""
Parser of the Pathway Commons extended binary SIF file.

``pc-hgnc.txt.gz`` (the "TXT or EXTENDED_BINARY_SIF" format of Pathway Commons)
is a tab-separated text file of two sections, each with its own header line and
separated by one blank line.

The first section has the binary relations (edges), one relation of a type
between two participants in each line:

- ``PARTICIPANT_A``, ``PARTICIPANT_B``: the participants, identified by HGNC
  symbol for genes and gene products, ``chebi:<number>`` (or
  ``CHEBI:<number>``) for small molecules;
- ``INTERACTION_TYPE``: the relation type, for example
  ``controls-state-change-of``;
- ``INTERACTION_DATA_SOURCE``: the databases the relation is derived from,
  semicolon separated;
- ``INTERACTION_PUBMED_ID``: PubMed IDs, semicolon separated;
- ``PATHWAY_NAMES``: names of the pathways of the relation, semicolon
  separated;
- ``MEDIATOR_IDS``: URIs of the BioPAX interactions the relation is derived
  from, semicolon separated.

The second section, after the edges, describes each participant once:

- ``PARTICIPANT``: the participant, as in the edges;
- ``PARTICIPANT_TYPE``: the BioPAX reference types that merged into the
  participant, semicolon separated (``ProteinReference``,
  ``SmallMoleculeReference``, ``RnaReference``, ...);
- ``PARTICIPANT_NAME``: names of the BioPAX objects, semicolon separated;
- ``UNIFICATION_XREF``, ``RELATIONSHIP_XREF``: cross-references as
  ``database:accession``, semicolon separated.

A parsed row is one edge, with the description of both participants joined to
it. The participant section is at the end of the file, so the parser reads the
file twice: the second section first, then the edges. The second reading
opens the downloaded file again with the opener's own class and arguments.

A parsed row has:

- ``interaction_type``: the relation type, as written;
- ``participant_a``, ``participant_b``: dicts with ``id`` (as in the edge),
  ``types`` (list), ``names`` (list), ``unification_xrefs`` and
  ``relationship_xrefs`` (lists of ``(database, accession)`` pairs);
- ``data_sources``, ``pubmed``: lists of the semicolon-separated values, as
  written;
- ``pathway_names``, ``mediator_ids``: the cells as written.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable, Iterator
import csv
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

EDGE_COLUMNS = (
    'PARTICIPANT_A',
    'INTERACTION_TYPE',
    'PARTICIPANT_B',
    'INTERACTION_DATA_SOURCE',
    'INTERACTION_PUBMED_ID',
    'PATHWAY_NAMES',
    'MEDIATOR_IDS',
)
NODE_COLUMNS = (
    'PARTICIPANT',
    'PARTICIPANT_TYPE',
    'PARTICIPANT_NAME',
    'UNIFICATION_XREF',
    'RELATIONSHIP_XREF',
)


def _split(value: str) -> list[str]:
    return value.split(';') if value else []


def _xrefs(value: str) -> list[tuple[str, str]]:
    """``database:accession`` items; the accession can have colons itself."""
    pairs = []
    for item in _split(value):
        database, separator, accession = item.partition(':')
        if not separator or not database or not accession:
            raise ValueError(f'Pathway Commons cross-reference {item!r} is malformed')
        pairs.append((database, accession))
    return pairs


def _edge_lines(lines: Iterable[str]) -> Iterator[str]:
    """The lines before the blank line that ends the edge section."""
    for line in lines:
        if not line.strip():
            return
        yield line


def _node_lines(lines: Iterable[str]) -> Iterator[str]:
    """The lines after the blank line that ends the edge section."""
    lines = iter(lines)
    for line in lines:
        if not line.strip():
            break
    else:
        raise ValueError('Pathway Commons file has no participant section')
    yield from lines


# A cell of the file can be longer than the reader's default limit of 131072
# characters: ``MEDIATOR_IDS`` of the longest relation has 13297 URIs.
_CELL_LIMIT = 16 * 1024 * 1024


def _table(lines: Iterable[str], columns: tuple[str, ...]):
    if csv.field_size_limit() < _CELL_LIMIT:
        csv.field_size_limit(_CELL_LIMIT)
    return iter_tsv(
        SimpleNamespace(result=lines),
        quotechar=None,
        strict=True,
        required=columns,
    )


def _reopen(opener):
    """The same file, opened again from its start, to read the other section."""
    path = getattr(opener, 'path', None)
    if path is None:
        raise ValueError('Pathway Commons parser needs the path of the downloaded file')
    return type(opener)(
        path=path,
        ext=getattr(opener, 'ext', None),
        large=True,
        default_mode='r',
        encoding='utf-8',
    )


def iter_participants(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One dict for each participant of the second section."""
    for row in _table(_node_lines(opener.result), NODE_COLUMNS):
        yield {
            'id': row['PARTICIPANT'],
            'types': _split(row['PARTICIPANT_TYPE']),
            'names': _split(row['PARTICIPANT_NAME']),
            'unification_xrefs': _xrefs(row['UNIFICATION_XREF']),
            'relationship_xrefs': _xrefs(row['RELATIONSHIP_XREF']),
        }


def iter_interactions(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One dict for each edge, with its participants described.

    A header that lacks a column, a line with another number of fields, an
    edge with a participant that the second section does not describe and a
    file without edges or without a second section raise ``ValueError``.
    """
    participants = {
        participant['id']: participant
        for participant in iter_participants(_reopen(opener))
    }
    if not participants:
        raise ValueError('Pathway Commons file has no participants')
    for row in _table(_edge_lines(opener.result), EDGE_COLUMNS):
        try:
            participant_a = participants[row['PARTICIPANT_A']]
            participant_b = participants[row['PARTICIPANT_B']]
        except KeyError as error:
            raise ValueError(
                f'Pathway Commons participant {error.args[0]!r} is not described'
            ) from None
        yield {
            'interaction_type': row['INTERACTION_TYPE'],
            'participant_a': participant_a,
            'participant_b': participant_b,
            'data_sources': _split(row['INTERACTION_DATA_SOURCE']),
            'pubmed': _split(row['INTERACTION_PUBMED_ID']),
            'pathway_names': row['PATHWAY_NAMES'],
            'mediator_ids': row['MEDIATOR_IDS'],
        }
