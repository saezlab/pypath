"""Read the complex files of the three versions of hu.MAP.

hu.MAP (Human Protein Complex Map) publishes the complexes of each version in
a file of its own, with other identifiers and another way to write the
members. Each reader yields one dictionary per complex with these fields:

- ``complex``: the complex identifier of the version as a string
  (``HuMAP2_00000``, ``huMAP3_00000.1``; the first file has none, see below)
- ``confidence``: the confidence level of the complex as the file writes it
  (``'1'``); None for version 1, which has no confidence
- ``members``: the proteins of the complex, in file order, each as a
  dictionary with the keys ``uniprot``, ``entrez``, ``ensembl``, ``entry`` (a
  UniProt entry name) and ``genename``; the key that the file gives is the
  string, the others are None

Version 1 (``clusters.txt``) has one complex per line, written as the
Entrez Gene IDs of its members separated by spaces; some members are Ensembl
gene IDs and one is a UniProt entry name (``GFP_AEQVI``, the green fluorescent
protein of Aequorea victoria, which is not a human protein). The file has no identifier. The site's own node table
(``nodeTable.txt``) numbers the complexes from 0 (``clustid``) in the order of
this file (checked: all 4659 complexes have the same members in both), so the
identifier of a complex is its position in the file, from 0.

Versions 2 and 3 are CSV files with the columns complex identifier, confidence,
UniProt accessions (separated by single spaces) and gene names. In version 2
the gene names are in the same positions as the accessions, a protein without
accession having an empty position, and a protein listed twice in a complex is
one member. In version 3 the gene names are not in the positions of the
accessions (some entries hold several names separated by ``;``, some proteins
have none), so they are not read and the members are the accessions. A member
written as a number is an Entrez Gene ID, as the file writes it.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_csv, read_text

_UNIPROT = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?'
)
_ENSEMBL_GENE = re.compile(r'ENSG\d{11}')
_UNIPROT_ENTRY = re.compile(r'[A-Z0-9]{1,5}_[A-Z0-9]{3,5}')

# Version -> the columns of its table: identifier, confidence.
TABLE_COLUMNS = {
    2: ('HuMAP2_ID', 'Confidence', 'Uniprot_ACCs', 'genenames'),
    3: ('HuMAP3_ID', 'ComplexConfidence', 'Uniprot_ACCs', 'genenames'),
}


def _protein(
    token: str | None,
    genename: str | None = None,
    *,
    complex_id: str,
) -> dict[str, str | None]:
    """A member from its token: a UniProt accession or entry name, an Entrez or Ensembl gene ID."""
    member = {
        'uniprot': None,
        'entrez': None,
        'ensembl': None,
        'entry': None,
        'genename': genename,
    }
    if token is None:
        return member
    if _UNIPROT.fullmatch(token):
        member['uniprot'] = token
    elif token.isdigit():
        member['entrez'] = token
    elif _ENSEMBL_GENE.fullmatch(token):
        member['ensembl'] = token
    elif _UNIPROT_ENTRY.fullmatch(token):
        member['entry'] = token
    else:
        raise ValueError(f'hu.MAP complex {complex_id}: unknown member {token!r}')
    return member


def iter_humap1(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The complexes of hu.MAP 1: one line per complex, the members separated by spaces."""
    count = 0
    for index, line in enumerate(read_text(opener).splitlines()):
        tokens = line.split()
        if not tokens:
            raise ValueError(f'hu.MAP 1 line {index + 1} has no members')
        count += 1
        complex_id = str(index)
        yield {
            'complex': complex_id,
            'confidence': None,
            'members': [_protein(token, complex_id=complex_id) for token in tokens],
        }
    if not count:
        raise ValueError('hu.MAP 1 file has no complexes')


def iter_humap_table(
    opener,
    version: int,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """The complexes of the CSV file of hu.MAP 2 or 3; a changed header is an error."""
    id_column, confidence_column, accession_column, name_column = TABLE_COLUMNS[version]
    count = 0
    for row in iter_csv(opener, strict=True, quotechar=None):
        missing = [
            column
            for column in (id_column, confidence_column, accession_column, name_column)
            if column not in row
        ]
        if missing:
            raise ValueError(f'hu.MAP {version} table lacks columns {missing}')
        count += 1
        complex_id = row[id_column]
        # Single spaces separate the positions; an empty position has no accession.
        accessions = [item or None for item in row[accession_column].split(' ')]
        if version == 2:
            names = row[name_column].split(' ')
            if len(names) != len(accessions):
                raise ValueError(
                    f'hu.MAP 2 complex {complex_id}: {len(accessions)} accessions '
                    f'and {len(names)} gene names'
                )
            pairs = list(dict.fromkeys(zip(accessions, names, strict=True)))
        else:
            pairs = [(accession, None) for accession in dict.fromkeys(accessions)]
        yield {
            'complex': complex_id,
            'confidence': row[confidence_column],
            'members': [
                _protein(accession, name, complex_id=complex_id)
                for accession, name in pairs
                if accession or name
            ],
        }
    if not count:
        raise ValueError(f'hu.MAP {version} table has no complexes')
