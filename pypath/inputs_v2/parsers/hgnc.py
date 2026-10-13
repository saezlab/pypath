"""
Parser of the HGNC gene groups.

HGNC (the HUGO Gene Nomenclature Committee) groups genes into gene groups, also
called gene families, by shared structure, function or origin. A group has members
and can have subgroups. The groups are published as tables of the database
(``genefamily_db_tables``), which this parser joins:

- ``family.csv``: ``id`` (the group id), ``abbreviation``, ``name``,
  ``pubmed_ids``, ``desc_comment`` (the description), ``desc_go`` (a GO term that
  describes the group, as a UniProt URL) and other columns the parser leaves
  unread. ``NULL`` and the empty string both stand for no value.
- ``family_alias.csv``: ``family_id``, ``alias``.
- ``hierarchy.csv``: ``parent_fam_id``, ``child_fam_id``. A group can have several
  parents.
- ``gene_has_family.csv``: ``hgnc_id`` (digits), ``family_id``. A gene is listed with
  the groups it belongs to directly; most genes of a subgroup are not repeated in its
  parents, some (120 pairs of a parent and a subgroup) are.

The tables name genes by HGNC id only. The approved symbols come from the web
download interface of HGNC (``custom`` download with the columns HGNC ID and
Approved symbol). Genes that are no longer approved have no symbol there.
"""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_csv, iter_tsv

FAMILY_COLUMNS = ('id', 'abbreviation', 'name', 'pubmed_ids', 'desc_comment', 'desc_go')
_GO_URL = re.compile(r'http://purl\.uniprot\.org/go/(\d{7})')
_NULL = ('', 'NULL')


def _value(row: dict[str, Any], column: str) -> str | None:
    text = (row.get(column) or '').strip()
    return None if text in _NULL else text


def _table(download, columns, force_refresh, **kwargs):
    return iter_csv(
        download.open(force_refresh=force_refresh, **kwargs), required=columns
    )


def symbols(
    symbols_download, force_refresh: bool = False, **kwargs: Any
) -> dict[str, str]:
    """HGNC id (digits) -> approved symbol, from the two-column custom download."""
    result = {}
    opener = symbols_download.open(force_refresh=force_refresh, **kwargs)
    for row in iter_tsv(opener, quotechar=None, required=('HGNC ID', 'Approved symbol')):
        match = re.fullmatch(r'HGNC:(\d+)', row['HGNC ID'] or '')
        symbol = (row['Approved symbol'] or '').strip()
        if not match or not symbol or match.group(1) in result:
            raise ValueError(f'HGNC symbol table row not understood: {row!r}')
        result[match.group(1)] = symbol
    return result


def iter_groups(
    opener,
    *,
    alias_download,
    hierarchy_download,
    members_download,
    symbols_download,
    force_refresh: bool = False,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    One dictionary per gene group.

    - ``id``: the group id, digits.
    - ``name``, ``abbreviation`` (or None), ``aliases``: the names of the group.
    - ``description``: the description text, or None.
    - ``pmids``: the PubMed ids of the group's references, as written.
    - ``go``: the GO term of the description, as ``GO:`` and seven digits, or None.
    - ``children``: the ids of the subgroups, in the order of the hierarchy table.
    - ``members``: the genes of the group as ``{'hgnc_id': ..., 'symbol': ...}``;
      the symbol is None for a gene without an approved symbol.

    A group that has no name, a group id or an HGNC id twice in a table, a reference
    to a group that does not exist, or a table without rows is an error.
    """
    aliases = defaultdict(list)
    for row in _table(alias_download, ('family_id', 'alias'), force_refresh, **kwargs):
        if _value(row, 'alias'):
            aliases[row['family_id']].append(row['alias'].strip())
    children = defaultdict(list)
    for row in _table(
        hierarchy_download, ('parent_fam_id', 'child_fam_id'), force_refresh, **kwargs
    ):
        children[row['parent_fam_id']].append(row['child_fam_id'])
    approved = symbols(symbols_download, force_refresh=force_refresh, **kwargs)
    members = defaultdict(list)
    seen_members = set()
    for row in _table(members_download, ('hgnc_id', 'family_id'), force_refresh, **kwargs):
        pair = (row['hgnc_id'], row['family_id'])
        if not row['hgnc_id'].isdigit() or pair in seen_members:
            raise ValueError(f'HGNC gene_has_family row not understood: {row!r}')
        seen_members.add(pair)
        members[row['family_id']].append(
            {'hgnc_id': row['hgnc_id'], 'symbol': approved.get(row['hgnc_id'])}
        )

    ids = set()
    for row in iter_csv(opener, required=FAMILY_COLUMNS):
        group = (row['id'] or '').strip()
        name = _value(row, 'name')
        if not group.isdigit() or not name or group in ids:
            raise ValueError(f'HGNC gene group row not understood: {row!r}')
        ids.add(group)
        go = _value(row, 'desc_go')
        if go:
            match = _GO_URL.fullmatch(go)
            if not match:
                raise ValueError(f'HGNC gene group {group}: GO link not understood: {go!r}')
            go = f'GO:{match.group(1)}'
        pmids = [
            token
            for token in re.split(r'[,\s]+', _value(row, 'pubmed_ids') or '')
            if token
        ]
        if not all(is_pmid(token) for token in pmids):
            raise ValueError(f'HGNC gene group {group}: PubMed ids not understood: {pmids}')
        yield {
            'id': group,
            'name': name,
            'abbreviation': _value(row, 'abbreviation'),
            'aliases': list(dict.fromkeys(aliases.get(group, []))),
            'description': _value(row, 'desc_comment'),
            'pmids': list(dict.fromkeys(pmids)),
            'go': go,
            'children': children.get(group, []),
            'members': members.get(group, []),
        }
    if not ids:
        raise ValueError('HGNC gene group table has no rows')
    unknown = (
        {c for cs in children.values() for c in cs}
        | set(children)
        | set(members)
        | set(aliases)
    ) - ids
    if unknown:
        raise ValueError(f'HGNC tables refer to gene groups that do not exist: {sorted(unknown)[:5]}')
