"""
Reader for PepCyber:P~Pep, from OmniPath's rescued copies of its web pages.

PepCyber:P~Pep (http://www.pepcyber.org/PPEP/) was a database of interactions
between proteins with a phosphoprotein-binding domain (PPBD: SH2, 14-3-3, PTB,
FHA, BRCT, WW, ...) and their phosphorylated substrates. The website is gone
(the domain is parked). OmniPath keeps two kinds of page:

- the result page of the search for all interactions, one table row per
  interaction (9,530 of them): :func:`iter_pepcyber`;
- the detail page of some interactions (``pepcyber_<n>.html``, where ``n`` is
  the interaction number of the search page), with the UniProt and RefSeq
  accessions of the two proteins and the publications of the records. Only 704
  of the 9,530 pages were rescued; the directory listing of the rescue server
  says which. :func:`iter_pepcyber` reads the pages the listing names, and
  none of the others.

The search page is written in UTF-8 with a few stray Latin-1 bytes, so it is
read as Latin-1.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import (
    iter_html_links,
    iter_html_tables,
    read_text,
)

__all__ = ['SEARCH_COLUMNS', 'iter_detail_numbers', 'iter_pepcyber', 'read_detail']

ENCODING = 'latin-1'

# Cells of an interaction row of the search page, after the checkbox cell, with
# the field each is read into.
SEARCH_COLUMNS = (
    'ppbd_class',
    'ppbd_symbol',
    'substrate_symbol',
    'binding_sequence',
    'index_site',
    'evidence',
    'n_records',
    'category',
)
# The header cells of the search table that follow the checkbox cell.
_SEARCH_HEADER = (
    'PPBD class',
    'PPBD protein',
    'Substrate protein',
    'Binding sequence',
    'Index site',
    'Evidence',
    '#Records',
    'Category',
)
_DETAIL_FILE = re.compile(r'pepcyber_(\d+)\.html')
_DETAIL_LINK = re.compile(r'idetail\.php\?iid=(\d+)')
_ENTREZ_LINK = re.compile(r'[?&]list_uids=(\d+)')
_PUBMED_LINK = re.compile(r'[?&]list_uids=(\d+)')
_GENE_NAME = re.compile(r'(?P<symbol>[^(]*?)\s*(?:\(\s*Description:\s*(?P<description>.*?)\s*\))?')
_NOT_AVAILABLE = 'N/A'
# A lower case ``p`` before an upper case residue letter marks the phosphorylated residue.
_PHOSPHO_MARK = re.compile(r'p([A-Z])')


def iter_detail_numbers(opener, **_kwargs: Any) -> Generator[int, None, None]:
    """The interaction numbers of the rescued detail pages named by the directory listing."""
    html = read_text(opener, encoding=ENCODING)
    found = False
    for link in iter_html_links(html):
        match = _DETAIL_FILE.fullmatch(link['href'])
        if match:
            found = True
            yield int(match.group(1))
    if not found:
        raise ValueError('PepCyber detail directory lists no pepcyber_<n>.html page')


def _protein(table: dict[str, Any]) -> dict[str, Any]:
    cells = {row[0]['text']: row[1] for row in table['rows'] if len(row) >= 2}
    missing = [k for k in ('Gene name', 'RefSeq accession', 'SwissProt accession') if k not in cells]
    if missing:
        raise ValueError('PepCyber protein table lacks: ' + ', '.join(missing))
    name = cells['Gene name']
    match = _GENE_NAME.fullmatch(name['text'])
    entrez = next(
        (m.group(1) for link in name['links'] if (m := _ENTREZ_LINK.search(link))),
        None,
    )
    return {
        'symbol': match.group('symbol') if match else name['text'],
        'description': (match.group('description') or None) if match else None,
        'entrez': entrez,
        'refseq': cells['RefSeq accession']['text'] or None,
        'swissprot': cells['SwissProt accession']['text'] or None,
    }


def read_detail(html: str) -> dict[str, Any]:
    """The content of one detail page.

    Keys:

    - ``ppbd``, ``substrate``: for the two proteins, ``symbol`` (the gene name),
      ``description`` (the gene's description, if the page has one), ``entrez``
      (the Gene ID in the link of the gene name), ``refseq`` and ``swissprot``
      (the accessions the page gives, which are not always of the same
      species as the other protein's);
    - ``records``: the records behind the interaction, each ``category``,
      ``evidence`` (as written), ``pmid`` (the PubMed ID in the reference's
      link) and ``reference`` (the citation text).

    A page without the two protein tables or the table of records raises
    ``ValueError``.
    """
    proteins = []
    records = None
    for table in iter_html_tables(html):
        first = table['rows'][0] if table['rows'] else []
        if first and first[0]['text'] == 'Gene name':
            proteins.append(_protein(table))
        elif tuple(cell['text'] for cell in first) == ('Category', 'Evidence', 'Reference'):
            records = []
            for row in table['rows'][1:]:
                if len(row) != 3:
                    raise ValueError(f'PepCyber record row has {len(row)} cells, expected 3')
                pmid = next(
                    (m.group(1) for link in row[2]['links'] if (m := _PUBMED_LINK.search(link))),
                    None,
                )
                records.append(
                    {
                        'category': row[0]['text'],
                        'evidence': row[1]['text'],
                        'pmid': pmid,
                        'reference': row[2]['text'],
                    }
                )
    if len(proteins) != 2 or records is None:
        raise ValueError('PepCyber detail page lacks the two protein tables or the records table')
    return {'ppbd': proteins[0], 'substrate': proteins[1], 'records': records}


def _phospho_residue(sequence: str | None, site: int | None) -> str | None:
    """The residue after the one ``p`` of the peptide, if the page gives a site.

    A peptide with two phosphorylated residues (one row has) or none, and a
    peptide without an index site, give no residue: which residue the site
    is would be a guess.
    """
    if not sequence or site is None:
        return None
    marked = _PHOSPHO_MARK.findall(sequence)
    return marked[0] if len(marked) == 1 else None


def _search_rows(html: str) -> Generator[list[dict[str, Any]], None, None]:
    """The interaction rows of the search result table, as lists of cells."""
    header = None
    for table in iter_html_tables(html):
        for row in table['rows']:
            texts = tuple(cell['text'] for cell in row)
            if texts[1 : 1 + len(_SEARCH_HEADER)] == _SEARCH_HEADER:
                header = row
            elif header is not None and len(row) == len(header):
                yield row
    if header is None:
        raise ValueError(
            'PepCyber page has no table with the columns: ' + ', '.join(_SEARCH_HEADER)
        )


def iter_pepcyber(
    opener,
    details=None,
    detail_numbers=None,
    force_refresh: bool = False,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield the interactions of the PepCyber search result page.

    ``details`` is the ``Download`` of a detail page (its URL and file name
    take the keyword ``number``) and ``detail_numbers`` the ``Download`` of the
    directory listing that names the rescued pages; with both, the pages of the
    interactions that have one are read too, each from the cache or downloaded
    once. Without them only the search page is read.

    Fields of a yielded record, all as written on the page unless said:

    - ``number``: PepCyber's number of the interaction (int);
    - ``ppbd_class``: the binding domain family (``SH2``, ``14-3-3``, ``WW``...);
    - ``ppbd_symbol``, ``substrate_symbol``: the Entrez gene names of the two
      proteins;
    - ``substrate_entrez``: the Gene ID in the link of the substrate's name
      (None for 84 rows without a link);
    - ``binding_sequence``: the peptide of the substrate that binds, a lower
      case ``p`` before the phosphorylated residue, or None;
    - ``index_site``: the position of the phosphosite in the substrate (int),
      None where the page says ``N/A``;
    - ``phospho_residue``: the one letter code of the phosphorylated residue,
      the letter after the ``p`` of the peptide, when the peptide has exactly
      one ``p`` and the row has an index site; else None;
    - ``evidence``: the evidence types as written (list of strings; the source
      writes ``,Inference`` in one row, whose empty entry is dropped);
    - ``n_records``: the number of records behind the interaction (int);
    - ``category``: ``A`` (known to be PPBD-mediated, peptide identified),
      ``B`` (PPBD-mediated, peptide not identified) or ``C`` (not known whether
      a PPBD mediates it), after Gong et al. 2008;
    - ``detail``: the output of :func:`read_detail` for the interaction, or
      None where no page was rescued.

    A page without the table, a row with a field that is not as expected and
    a duplicated interaction number raise ``ValueError``.
    """
    html = read_text(opener, encoding=ENCODING)
    available: set[int] = set()
    if details is not None and detail_numbers is not None:
        available = set(iter_detail_numbers(detail_numbers.open(force_refresh=force_refresh)))
    seen: set[int] = set()
    for row in _search_rows(html):
        cells = row[1:]
        number = next(
            (
                int(m.group(1))
                for cell in row[-1:]
                for link in cell['links']
                if (m := _DETAIL_LINK.search(link))
            ),
            None,
        )
        if number is None:
            raise ValueError('PepCyber interaction row has no detail link')
        if number in seen:
            raise ValueError(f'PepCyber interaction {number} is listed twice')
        seen.add(number)
        record = {
            field: cell['text'].strip() or None
            for field, cell in zip(SEARCH_COLUMNS, cells, strict=False)
        }
        record['number'] = number
        record['substrate_entrez'] = next(
            (m.group(1) for link in cells[2]['links'] if (m := _ENTREZ_LINK.search(link))),
            None,
        )
        site = record['index_site']
        if site is not None and site != _NOT_AVAILABLE and not site.isdigit():
            raise ValueError(f'PepCyber interaction {number}: index site {site!r}')
        record['index_site'] = int(site) if site is not None and site.isdigit() else None
        record['phospho_residue'] = _phospho_residue(
            record['binding_sequence'], record['index_site']
        )
        count = record['n_records']
        if count is None or not count.isdigit():
            raise ValueError(f'PepCyber interaction {number}: record count {count!r}')
        record['n_records'] = int(count)
        # One row writes ``,Inference``: the empty first entry is dropped.
        record['evidence'] = [
            part.strip() for part in (record['evidence'] or '').split(',') if part.strip()
        ]
        record['detail'] = None
        if number in available:
            record['detail'] = read_detail(
                read_text(
                    details.open(force_refresh=force_refresh, number=number),
                    encoding=ENCODING,
                )
            )
        yield record
