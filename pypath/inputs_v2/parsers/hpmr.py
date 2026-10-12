"""
Read the receptor pages of the Human Plasma Membrane Receptome (HPMR).

:func:`parse_receptor_page` and :func:`iter_interactions` interpret the HPMR
receptor pages; the tables of a page are read by the shared HTML reader.

HPMR has one page for each receptor protein. Under the heading ``Receptor X
interacts with:`` the page has a table with one cell for each interaction
entry of the database. An entry lists one or more partners, each marked as a
``Ligand``, a ``Receptor`` or a ``non-gene ligand``, followed by a comment
with links to the PubMed records that support the entry.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

import bs4

from pypath.inputs_v2._publications import PMID_PATTERN
from pypath.inputs_v2.parsers.base import content_text, iter_html_tables

__all__ = [
    'NON_GENE',
    'iter_interactions',
    'parse_receptor_page',
    'pubmed_ids',
]


def _text(value: str) -> str:
    return ' '.join(value.split())


# The heading of the interaction table names the page's receptor.
_heading_re = re.compile(r'Receptor (\S+) interacts with:')
# The page header: gene symbol and the title of the protein record.
_header_re = re.compile(r'(\S+)\s+-\s+(.*\S)', re.S)
# A partner links to its own page; the address holds the site's key for it.
_partner_href_re = re.compile(r'Param=([^&]*)&ProtId=(\d+)&ProtType=(\w+)')
_partner_text_re = re.compile(r'(\S+) - (.*)', re.S)
_list_uids_re = re.compile(r'list_uids=(.*)$', re.S)
_gene_roles = ('Ligand', 'Receptor')
NON_GENE = 'non-gene ligand'
_non_gene_label = NON_GENE + ':'
# A PubMed ID as the reference links write it; some carry a `PMID:` prefix
# or a stray closing bracket or full stop.
_pmid_re = re.compile(rf'(?:PMID:\s*)?({PMID_PATTERN})[).]?')
# When the site was copied its server answered some requests with an error
# report instead of the page.
_error_title_re = re.compile(r'\w+ at /SearchDB/getGenePage\.asp')


def _entry(cell: dict[str, Any]) -> dict[str, Any] | None:
    """One interaction entry from a cell of the interaction table."""
    partners: list[dict[str, Any]] = []
    pending: dict[str, Any] | None = None
    interaction_id = None
    in_comment = False
    comment: list[str] = []
    references: list[str] = []
    for item in cell['content']:
        if in_comment:
            if isinstance(item, str):
                comment.append(item)
                continue
            comment.append(item['text'])
            if uids := _list_uids_re.search(item['href']):
                references.append(uids.group(1))
            continue
        if not isinstance(item, str):
            if pending is not None:
                # A non-gene ligand without a name.
                partners.append(pending)
                pending = None
            if item['attrs'].get('id', '').startswith('OpenComm_'):
                interaction_id = item['attrs']['id'].removeprefix('OpenComm_')
                in_comment = True
                continue
            target = _partner_href_re.search(item['href'])
            if not target or item['text'] not in _gene_roles:
                raise ValueError(
                    f'Unexpected link in an HPMR interaction: {item["href"]!r}'
                )
            pending = {
                'role': item['text'],
                'symbol': None,
                'title': None,
                'page_id': target.group(1),
            }
            continue
        text = _text(item)
        if text == 'and':
            continue
        if text == _non_gene_label:
            if pending is not None:
                partners.append(pending)
            pending = {
                'role': NON_GENE,
                'symbol': None,
                'title': None,
                'page_id': None,
            }
        elif pending is None:
            raise ValueError(f'Unexpected text in an HPMR interaction: {text!r}')
        elif pending['role'] == NON_GENE:
            pending['symbol'] = text
            partners.append(pending)
            pending = None
        else:
            named = _partner_text_re.fullmatch(text)
            if not named:
                raise ValueError(f'Unexpected HPMR partner: {text!r}')
            pending['symbol'], pending['title'] = named.groups()
            partners.append(pending)
            pending = None
    if pending is not None:
        partners.append(pending)
    if not partners:
        return None
    return {
        'interaction_id': interaction_id,
        'partners': partners,
        'comment': _text(' '.join(comment)) or None,
        'references': references,
    }


def parse_receptor_page(html: str) -> list[dict[str, Any]] | None:
    """Interaction records of one HPMR receptor page.

    Returns None for a server error report saved in place of a page, and an
    empty list for a receptor without interactions. There is one record for
    each partner of each interaction entry:

    - ``receptor``: gene symbol of the page's receptor, from the heading of
      the interaction table.
    - ``receptor_title``: title of the page's protein record.
    - ``receptor_entrez``: NCBI Gene ID from the page's ``Entrez Gene`` link.
    - ``interaction_id``: the site's number of the entry's comment, shared by
      the pages that show the same entry; None for an entry without comment.
    - ``partner_role``: ``Ligand``, ``Receptor`` or ``non-gene ligand``.
    - ``partner``: gene symbol of the partner, or the name of a non-gene
      ligand; None if the page gives none.
    - ``partner_title``: the text after the partner's symbol, which the site
      cuts at about 45 characters.
    - ``partner_page_id``: the site's key of the partner's own page (a mix of
      NCBI Gene IDs, GI numbers and protein accessions).
    - ``entry``: every partner of the entry as ``role symbol``, in page order.
    - ``comment``: the entry's comment.
    - ``references``: the PubMed IDs of the comment's links, as written.
    """
    soup = bs4.BeautifulSoup(html, 'html.parser')
    interactions = soup.find('div', id='GeneInts')
    if interactions is None:
        title = soup.find('title')
        if title and _error_title_re.search(title.get_text()):
            return None
        raise ValueError('HPMR page has no interaction section')
    heading = interactions.find_previous_sibling('span')
    receptor = _heading_re.search(_text(heading.get_text()) if heading else '')
    if not receptor:
        raise ValueError('HPMR page does not name its receptor')
    header = soup.find('span', class_='supertitle')
    titled = _header_re.fullmatch(header.get_text().strip()) if header else None
    entrez = {
        match.group(1)
        for link in soup.find_all('a')
        if link.get_text().strip() == 'Entrez Gene'
        and (match := re.search(r'list_uids=(\d+)', link.get('href') or ''))
    }
    if len(entrez) > 1:
        raise ValueError(
            f'HPMR page of {receptor.group(1)} links several genes: {sorted(entrez)}'
        )
    page = {
        'receptor': receptor.group(1),
        'receptor_title': _text(titled.group(2)) if titled else None,
        'receptor_entrez': next(iter(entrez), None),
    }
    records = []
    # The shared reader takes markup: the interaction section, written out.
    for table in iter_html_tables(str(interactions)):
        for row in table['rows']:
            for cell in row:
                entry = _entry(cell)
                if entry is None:
                    continue
                listed = [
                    f'{partner["role"]} {partner["symbol"] or ""}'.strip()
                    for partner in entry['partners']
                ]
                for partner in entry['partners']:
                    records.append(
                        {
                            **page,
                            'interaction_id': entry['interaction_id'],
                            'partner_role': partner['role'],
                            'partner': partner['symbol'],
                            'partner_title': partner['title'],
                            'partner_page_id': partner['page_id'],
                            'entry': listed,
                            'comment': entry['comment'],
                            'references': entry['references'],
                        }
                    )
    return records


def pubmed_ids(references: Any) -> list[str]:
    """PubMed IDs of a record's ``references`` in page order, each once; an entry
    that is not a PubMed ID is left out."""
    out: list[str] = []
    for reference in references or ():
        match = _pmid_re.fullmatch(str(reference).strip())
        if match and match.group(1) not in out:
            out.append(match.group(1))
    return out


def iter_interactions(
    opener,
    roles: tuple[str, ...] | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield the interaction records of every receptor page in the archive.

    The records are those of :func:`parse_receptor_page`, with ``page``, the
    name of the archive member, added. ``roles`` keeps only partners with
    these roles.

    HPMR has a page for each protein record of a gene and shows an entry on
    the page of every receptor it involves. A record is yielded once: for the
    first page, in name order, on which an entry lists a partner for a
    receptor gene, and once for the two views of an entry between two
    receptors.
    """
    result = getattr(opener, 'result', None) if opener else None
    if not result:
        return
    if not isinstance(result, dict):
        raise ValueError('Expected the archive of HPMR receptor pages')
    seen: set[tuple[Any, ...]] = set()
    pages = 0
    count = 0
    for name in sorted(result):
        if not name.endswith('.html'):
            continue
        records = parse_receptor_page(content_text(result[name]))
        if records is None:
            continue
        pages += 1
        for record in records:
            if roles is not None and record['partner_role'] not in roles:
                continue
            pair: Any = (record['receptor'], record['partner'])
            if record['interaction_id'] is None:
                key = (None, record['partner_role'], pair, *record['references'])
            elif record['partner_role'] == 'Receptor':
                key = (record['interaction_id'], 'Receptor', frozenset(pair))
            else:
                key = (record['interaction_id'], record['partner_role'], pair)
            if key in seen:
                continue
            seen.add(key)
            count += 1
            yield {**record, 'page': name.rsplit('/', 1)[-1]}
    if not pages:
        raise ValueError('The HPMR archive holds no receptor page')
    if not count:
        raise ValueError('The HPMR receptor pages hold no interaction')
