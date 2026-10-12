"""Read the KEA phosphosite-level kinase-substrate file.

The file is tab separated, has no header row and holds one
kinase-site-reference-database combination per line:

    ABL1    ABL1_Y226    10964922    phosphosite
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = ('enzyme', 'substrate_site', 'reference', 'source_database')
# "<gene symbol>_<residue letter><position>", e.g. "HLA-A_Y344".
_SITE_REGEX = re.compile(r'^(\S+)_([A-Z])([1-9]\d*)$')
# Stricter than the shared PubMed identifier rule (1 to 8 digits), because this
# file has a known defect: it holds "0", 6-digit and 9-digit numbers in the
# reference column, and the 6-digit ones are truncated PubMed IDs, not papers
# on the pair. Looked up at NCBI: 898671 is a 1977 adenovirus paper while the
# pair LCK-VAV1 Y174 is from PMID 8986718; 602896 is "Cystitis emphysematosa";
# 934118 is a paper on childhood constipation. So a reference counts as a
# PubMed ID only with 7 or 8 digits and no leading zero.
_PMID_REGEX = re.compile(r'^[1-9]\d{6,7}$')
# KEA fills a missing kinase with "NA".
_MISSING = {'', 'NA'}


def iter_kea(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per line of the KEA phosphosite database.

    The file has no header row: the shared reader names the four columns and
    rejects a line with another number of fields.

    Output fields:

    - ``enzyme``: first column, the kinase's gene symbol; None if "NA".
    - ``substrate_site``: second column as written, e.g. ``ABL1_Y226``.
    - ``substrate``, ``residue``, ``position``: gene symbol, residue letter and
      integer position read from ``substrate_site``; all None if it does not
      have the ``SYMBOL_X123`` form.
    - ``reference``: third column as written (a PubMed ID, another reference
      label, a number that is not a PubMed ID, or empty).
    - ``pmid``: ``reference`` if it is a PubMed ID of 7 or 8 digits, else None
      (another label, "0" and numbers of other lengths are not PubMed IDs).
    - ``source_database``: fourth column.
    """
    for line in iter_tsv(opener, fieldnames=COLUMNS, strict=True, quotechar=None):
        row: dict[str, object] = {column: value.strip() for column, value in line.items()}
        site = _SITE_REGEX.match(row['substrate_site'])
        symbol, residue, position = site.groups() if site else (None,) * 3
        row['substrate'] = symbol
        row['residue'] = residue
        row['position'] = int(position) if position else None
        row['pmid'] = row['reference'] if _PMID_REGEX.match(row['reference']) else None
        if row['enzyme'] in _MISSING:
            row['enzyme'] = None
        yield row
