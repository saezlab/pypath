"""Read the Cancer Cell Map binary network of Pathway Commons 1.

The last release of Pathway Commons 1 (October 2011) distributes the Cancer
Cell Map as a "tab delimited network" of two tables: the edges, which name two
Pathway Commons records, a binary interaction rule and PubMed IDs, and the
nodes, which give each record's gene symbol, UniProt accessions, Entrez Gene
IDs, type and organism.

    CPATH_RECORD_ID_A  INTERACTION_TYPE  CPATH_RECORD_ID_B  GENE_SYMBOL_A  GENE_SYMBOL_B  INTERACTION_DATA_SOURCE  INTERACTION_PUBMED_ID
    101318             INTERACTS_WITH    140665             MAP3K4         AXIN1          CELL_MAP                 PUBMED:15262978;

The reader yields the edge rows with the attributes of both nodes added.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

# Columns of the edge table, in source order.
EDGE_COLUMNS = (
    'CPATH_RECORD_ID_A',
    'INTERACTION_TYPE',
    'CPATH_RECORD_ID_B',
    'GENE_SYMBOL_A',
    'GENE_SYMBOL_B',
    'INTERACTION_DATA_SOURCE',
    'INTERACTION_PUBMED_ID',
)
# Columns of the node table, in source order.
NODE_COLUMNS = (
    'CPATH_RECORD_ID',
    'GENE_SYMBOL',
    'UNIPROT_ACCESSION',
    'ENTREZ_GENE_ID',
    'CHEBI_ID',
    'NODE_TYPE',
    'NCBI_TAX_ID',
)
# Node attributes the reader adds to each edge, once per participant
# (``UNIPROT_ACCESSION_A``, ``UNIPROT_ACCESSION_B``, ...). The gene symbol is
# already in the edge table.
NODE_ATTRIBUTES = ('UNIPROT_ACCESSION', 'ENTREZ_GENE_ID', 'NODE_TYPE', 'NCBI_TAX_ID')
# Pathway Commons writes this where it has no value.
MISSING = 'NOT_SPECIFIED'
# Undirected rules are written once in each orientation.
UNDIRECTED = frozenset({'INTERACTS_WITH', 'IN_SAME_COMPONENT', 'REACTS_WITH'})


def _check_columns(row: dict[str, Any], columns: tuple[str, ...], table: str) -> None:
    missing = [column for column in columns if column not in row]
    if missing:
        raise ValueError(
            f'Cancer Cell Map {table} table lacks expected columns: {", ".join(missing)}'
        )


def join_nodes(
    edges: Iterable[dict[str, Any]],
    nodes: Iterable[dict[str, Any]],
) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of the edge table with the attributes of their two nodes.

    ``NOT_SPECIFIED`` becomes None. ``MIRROR`` is True on a row of an undirected
    rule that repeats, with the participants exchanged, a row already read.
    """
    by_id: dict[str, dict[str, Any]] = {}
    for node in nodes:
        _check_columns(node, NODE_COLUMNS, 'node')
        by_id[node['CPATH_RECORD_ID']] = node
    if not by_id:
        raise ValueError('Cancer Cell Map node table has no rows')
    read: set[tuple[str, str, str, str]] = set()
    for edge in edges:
        _check_columns(edge, EDGE_COLUMNS, 'edge')
        row = dict(edge)
        for side in 'AB':
            record_id = row[f'CPATH_RECORD_ID_{side}']
            if record_id not in by_id:
                raise ValueError(
                    f'Cancer Cell Map edge names a record without a node: {record_id!r}'
                )
            for attribute in NODE_ATTRIBUTES:
                row[f'{attribute}_{side}'] = by_id[record_id][attribute]
        row = {
            column: None if value == MISSING else value
            for column, value in row.items()
        }
        first, second = row['CPATH_RECORD_ID_A'], row['CPATH_RECORD_ID_B']
        statement = (row['INTERACTION_TYPE'], row['INTERACTION_PUBMED_ID'])
        row['MIRROR'] = (
            row['INTERACTION_TYPE'] in UNDIRECTED
            and first != second
            and (second, first, *statement) in read
        )
        if not row['MIRROR']:
            read.add((first, second, *statement))
        yield row
    if not read:
        raise ValueError('Cancer Cell Map edge table has no rows')


def iter_interactions(
    opener,
    force_refresh: bool = False,
    nodes_download=None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield the edges of the Cancer Cell Map, each with its nodes' attributes.

    ``opener`` is the edge table; ``nodes_download`` is the ``Download`` of the
    node table.
    """
    nodes = nodes_download.open(force_refresh=force_refresh)
    yield from join_nodes(iter_tsv(opener), iter_tsv(nodes))
