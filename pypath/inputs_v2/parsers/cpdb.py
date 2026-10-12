"""
Parser of the ConsensusPathDB human protein interaction file.

``ConsensusPathDB_human_PPI.gz`` is a tab-separated text file. Its first line
is a comment with the release, its second line is the header (also written as
a comment, ``#  source_databases``), and each following line is one
interaction of the integrated network. The columns:

- ``source_databases``: the databases that report the interaction, comma
  separated;
- ``interaction_publications``: PubMed IDs, comma separated. ConsensusPathDB
  writes ``unassigned<number>`` where a source gave none, and a few other
  tokens; they are kept as written, for the mapper to judge;
- six columns that identify the participants, comma separated and aligned
  (the n-th item of each column belongs to the n-th participant):
  UniProt entry name, UniProt accession, gene name, HGNC ID, Entrez gene ID
  and Ensembl gene ID. An item is empty where the participant has no such
  identifier;
- ``interaction_confidence``: a number, or ``NA``.

ConsensusPathDB's data access page states how to read the participants: a
row with one protein is a self-interaction, a row with more than two proteins
is a complex interaction, and where a participant is known only by genomic
identifiers that map to several UniProt entries, the entry names are
concatenated with dots (``RL40_HUMAN.RS27A_HUMAN.UBIQ_HUMAN``), because it is
unclear which of the gene products interacts. The other identifier columns
use the dot in the same way for several genes.

A parsed row has:

- ``source_databases``: list of database names, as written;
- ``publications``: list of the publication tokens, as written;
- ``participants``: list with one dict for each participant; the keys are
  ``uniprot_entry``, ``uniprot_id``, ``gene_name``, ``hgnc_id``,
  ``entrez_gene`` and ``ensembl_gene``, each a list of the dot-separated
  values (empty for an empty item);
- ``confidence``: the confidence as written, or None for ``NA``.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

# The header line, written as a comment: the first column name has the hash.
_FIRST_COLUMN = '#  source_databases'
_PARTICIPANT_COLUMNS = {
    'uniprot_entry': 'interaction_participants__uniprot_entry',
    'uniprot_id': 'interaction_participants__uniprot_id',
    'gene_name': 'interaction_participants__genename',
    'hgnc_id': 'interaction_participants__hgnc_id',
    'entrez_gene': 'interaction_participants__entrez_gene',
    'ensembl_gene': 'interaction_participants__ensembl_gene',
}
COLUMNS = (
    _FIRST_COLUMN,
    'interaction_publications',
    *_PARTICIPANT_COLUMNS.values(),
    'interaction_confidence',
)
_NO_CONFIDENCE = 'NA'


def _split(value: str, separator: str) -> list[str]:
    return value.split(separator) if value else []


def _participants(row: dict[str, Any]) -> list[dict[str, list[str]]]:
    """The participants of a row, from the aligned identifier columns."""
    columns = {
        key: row[column].split(',') for key, column in _PARTICIPANT_COLUMNS.items()
    }
    counts = {len(items) for items in columns.values()}
    if len(counts) != 1:
        raise ValueError(
            f'ConsensusPathDB participant columns are not aligned: '
            f'{ {key: len(items) for key, items in columns.items()} }'
        )
    return [
        {
            key: _split(columns[key][index], '.')
            for key in _PARTICIPANT_COLUMNS
        }
        for index in range(counts.pop())
    ]


def iter_interactions(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One dict for each interaction line, with the fields described above.

    A header that lacks a column, a line with another number of fields,
    participant columns that are not aligned and an empty file raise
    ``ValueError``.
    """
    rows = iter_tsv(
        opener,
        skip_rows=1,
        quotechar=None,
        strict=True,
        required=COLUMNS,
    )
    for row in rows:
        confidence = row['interaction_confidence']
        yield {
            'source_databases': _split(row[_FIRST_COLUMN], ','),
            'publications': _split(row['interaction_publications'], ','),
            'participants': _participants(row),
            'confidence': None if confidence == _NO_CONFIDENCE else confidence,
        }

