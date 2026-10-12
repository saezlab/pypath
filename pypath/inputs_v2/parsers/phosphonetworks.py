"""
Reader for the PhosphoNetworks high-resolution network.

The tab separated lines are read by the shared TSV reader;
:func:`iter_network` is the PhosphoNetworks interpretation of the file, in which
a ``>SYMBOL`` line names the substrate of the kinase-site lines below it.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['iter_network']


_COLUMNS = ('window', 'site', 'enzyme', 'score')
# "S177": residue letter and its position in the substrate's sequence.
_SITE_REGEX = re.compile(r'^([A-Z])([1-9]\d*)$')
# Fifteen residues around the site, the site in the middle.
_WINDOW_REGEX = re.compile(r'^[A-Z-]{15}$')
_WINDOW_CENTER = 7


def iter_network(opener, **_kwargs: object) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per kinase-site line of the high-resolution network.

    The file has no header row. Substrate lines (``>SYMBOL``) are not rows of
    their own: the symbol is added to the kinase-site lines that follow.

    Output fields:

    - ``substrate``: gene symbol of the closest substrate line above.
    - ``window``: first column, the sequence window of the site as written.
    - ``site``: second column as written, e.g. ``S177``.
    - ``residue``, ``position``: residue letter and integer position read from
      ``site``.
    - ``enzyme``: third column, the kinase's gene symbol.
    - ``score``: fourth column as written.

    A line that is neither a substrate line nor a kinase-site line of four
    fields, a site or window of another form, a window whose middle residue
    is not the site's, a kinase-site line before the first substrate line and
    an empty file raise ``ValueError``.
    """
    substrate = None
    found = False
    for number, line in enumerate(
        iter_tsv(opener, fieldnames=_COLUMNS, quotechar=None, allow_repeated=True),
        start=1,
    ):
        cells = {
            key: value.strip() if value is not None else None
            for key, value in line.items()
            if key is not None
        }
        extra = line.get(None)
        first = cells['window'] or ''
        if first.startswith('>'):
            if any(cells[key] for key in _COLUMNS[1:]) or extra:
                raise ValueError(f'PhosphoNetworks line {number} is not a substrate line: {first!r}')
            substrate = first[1:].strip()
            if not substrate:
                raise ValueError(f'PhosphoNetworks line {number} has an empty substrate')
            continue
        if extra or not all(cells[key] for key in _COLUMNS):
            raise ValueError(f'PhosphoNetworks line {number} is not a kinase-site line')
        if substrate is None:
            raise ValueError(f'PhosphoNetworks line {number} comes before any substrate')
        site = _SITE_REGEX.match(cells['site'])
        if not site or not _WINDOW_REGEX.match(first):
            raise ValueError(f'PhosphoNetworks line {number} has an unreadable site: {cells!r}')
        if first[_WINDOW_CENTER] != site[1]:
            raise ValueError(f'PhosphoNetworks line {number}: window does not centre on {site[0]}')
        found = True
        yield {
            'substrate': substrate,
            'window': first,
            'site': cells['site'],
            'residue': site[1],
            'position': int(site[2]),
            'enzyme': cells['enzyme'],
            'score': cells['score'],
        }
    if not found:
        raise ValueError('PhosphoNetworks file holds no kinase-site lines')
