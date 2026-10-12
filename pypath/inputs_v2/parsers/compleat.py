"""
Reader for the COMPLEAT complex download (``ComplexDownload?org=Human``, older export).

The file has 12 tab-separated columns, no header, CRLF line ends and a last
line ``null``. The column names are those of the web interface where it has
one (source, database, name, species, method, reference).
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['COLUMNS', 'SOURCES', 'iter_compleat_complexes']

COLUMNS = (
    'compleat_id',
    'size',
    'source',  # Literature | Predicted
    'top_term',  # representative term: the name of a curated complex, a GO term of a predicted one
    'terms',  # ';'-separated names or GO terms
    'extra',  # empty, except in one malformed row
    'method',  # prediction method, or detection methods of the curated complex
    'name',  # name in the source database; empty for predicted complexes
    'database',
    'species',  # species the complex was obtained from before the orthology mapping
    'pubmeds',  # ','-separated
    'members',  # ' '-separated Entrez Gene IDs
)
SOURCES = ('Literature', 'Predicted')

_entrez = re.compile(r'[1-9]\d*')
_accession = re.compile(r'HC\d+')
# ``MI:0019- coimmunoprecipitation | MI:0096- pull down``
_mi_term = re.compile(r'(MI:\d{4})- ?\S.*')


def iter_compleat_complexes(
    opener, source: str | None = None, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Complexes of the file, those of one ``source`` (Literature | Predicted) if given.

    The columns of :data:`COLUMNS` are yielded as they are, except:

    - ``members``: the list of Entrez Gene IDs;
    - ``publications``: the PubMed IDs of ``pubmeds`` that are well formed
      (one to eight digits), in order;
    - ``detection_methods``: the PSI-MI accessions of the terms ``method``
      writes as ``MI:0019- coimmunoprecipitation``; its other text (BioGRID
      experimental systems, free text, prediction methods) is not read.

    The last line of the file, ``null``, is not a complex. Any other row that
    does not have the 12 columns, an accession, a value of ``source`` and
    members that match ``size`` is an error, and so is a file without
    complexes.
    """
    if source is not None and source not in SOURCES:
        raise ValueError(f'Unknown COMPLEAT source: {source!r}')
    rows = 0
    for record in iter_tsv(opener, fieldnames=COLUMNS, quotechar=None):
        if record['source'] is None:
            if record['compleat_id'] != 'null':
                raise ValueError(f'COMPLEAT row with too few fields: {record}')
            continue
        if None in record:
            raise ValueError(f'COMPLEAT row with too many fields: {record}')
        rows += 1
        if not _accession.fullmatch(record['compleat_id']):
            raise ValueError(f'Unknown COMPLEAT accession: {record["compleat_id"]!r}')
        if record['source'] not in SOURCES:
            raise ValueError(f'Unknown COMPLEAT source: {record["source"]!r}')
        members = record['members'].split()
        if not members or not all(_entrez.fullmatch(member) for member in members):
            raise ValueError(f'COMPLEAT {record["compleat_id"]} has malformed members')
        if str(len(members)) != record['size']:
            raise ValueError(f'COMPLEAT {record["compleat_id"]} size does not match its members')
        if source is not None and record['source'] != source:
            continue
        publications = []
        for token in record['pubmeds'].split(','):
            token = token.strip()
            if is_pmid(token) and token not in publications:
                publications.append(token)
        methods = []
        for token in record['method'].split('|'):
            match = _mi_term.fullmatch(token.strip())
            if match and match.group(1) not in methods:
                methods.append(match.group(1))
        yield {
            **record,
            'members': members,
            'publications': publications,
            'detection_methods': methods,
        }
    if not rows:
        raise ValueError('COMPLEAT file has no complexes')
