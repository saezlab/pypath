"""Read the ACSN 2.0 binary relations between proteins.

ACSN publishes a projection of its master map onto pairs of proteins, a table
with one row per unordered pair of map entity names::

    INTERACTOR1  INTERACTION_TYPE                           INTERACTOR2  PMID
    VEGFA        HETERODIMER_ASSOCIATION;PROTEIN_INTERACTION  FLT1       10022831;11166270;8837607

``INTERACTION_TYPE`` is a ``;`` separated set of CellDesigner reaction and
modification types, ``PMID`` a ``;`` separated list (or ``N/A``). A second file,
``ACSN2_HUGO_Correspondence.gmt``, gives the HUGO names of every map entity
(one gene set per entity name). The reader yields, for each row of the table
and each statement its types give, the row with:

- ``HUGO1`` and ``HUGO2``: the HUGO names of the two entities, from the
  correspondence file (an empty list for an entity it does not list);
- ``STATEMENT``: the statement that the types of the row give, or None when
  none of its types gives one;
- ``STATEMENT_TYPES``: the interaction types of the row that give it.

The statement of an interaction type is a choice of the mapping: the reader
takes the types' statements from the caller (``interaction_types``) and the
order in which a row's statements are yielded (``statements``).
"""

from __future__ import annotations

from collections.abc import Collection, Generator, Iterable, Mapping
from typing import Any

from pypath.inputs_v2.parsers.base import iter_gmt, iter_tsv

# Columns of the relations table, in source order.
COLUMNS = ('INTERACTOR1', 'INTERACTION_TYPE', 'INTERACTOR2', 'PMID')


def hugo_names(sets: Iterable[dict[str, Any]]) -> dict[str, list[str]]:
    """The HUGO names of every map entity, from the correspondence file.

    A few entity names are written with a space before or after; the names
    are compared without it.
    """
    names: dict[str, list[str]] = {}
    for entity in sets:
        members = names.setdefault(entity['name'].strip(), [])
        for member in entity['members']:
            member = member.strip()
            if member not in members:
                members.append(member)
    if not names:
        raise ValueError('ACSN HUGO correspondence file has no entities')
    return names


def join_names(
    rows: Iterable[dict[str, Any]],
    names: dict[str, list[str]],
    interaction_types: Mapping[str, str | None],
    statements: Collection[str],
) -> Generator[dict[str, Any], None, None]:
    """Yield one row per statement of each row of the relations table.

    Adds ``HUGO1`` and ``HUGO2`` (the HUGO names of the two entities),
    ``STATEMENT`` and ``STATEMENT_TYPES`` (the interaction types of the row
    that give the statement). A row whose types give no statement is yielded
    once, with ``STATEMENT`` None. ``interaction_types`` gives the statement
    of every known interaction type (None for one that gives none; an unknown
    type is an error) and ``statements`` lists the statements in the order a
    row's statements are yielded.
    """
    count = 0
    for row in rows:
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(
                f'ACSN relations table lacks expected columns: {", ".join(missing)}'
            )
        count += 1
        row = dict(row)
        for side in '12':
            row[f'INTERACTOR{side}'] = (row[f'INTERACTOR{side}'] or '').strip()
            row[f'HUGO{side}'] = list(names.get(row[f'INTERACTOR{side}'], ()))
        by_statement: dict[str | None, list[str]] = {}
        for interaction_type in (row['INTERACTION_TYPE'] or '').split(';'):
            if interaction_type not in interaction_types:
                raise ValueError(
                    f'Unknown ACSN interaction type: {interaction_type!r}'
                )
            by_statement.setdefault(
                interaction_types[interaction_type], []
            ).append(interaction_type)
        stated = [statement for statement in statements if statement in by_statement]
        for statement in stated or [None]:
            yield {
                **row,
                'STATEMENT': statement,
                'STATEMENT_TYPES': by_statement[statement],
            }
    if not count:
        raise ValueError('ACSN relations table has no rows')


def iter_relations(
    opener,
    force_refresh: bool = False,
    names_download=None,
    interaction_types: Mapping[str, str | None] | None = None,
    statements: Collection[str] = (),
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield the statements of the relations table with the HUGO names of the pair.

    ``opener`` is the relations table; ``names_download`` is the ``Download`` of
    the HUGO correspondence file; ``interaction_types`` and ``statements`` are
    as for :func:`join_names`.
    """
    names = hugo_names(iter_gmt(names_download.open(force_refresh=force_refresh)))
    yield from join_names(iter_tsv(opener), names, interaction_types, statements)
