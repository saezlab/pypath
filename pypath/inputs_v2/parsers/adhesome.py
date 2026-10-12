"""
Parser of the integrin adhesome tables.

Two kinds of file are read.

The supplementary tables of Zaidel-Bar et al. 2007 (spreadsheets, read with
:func:`iter_excel`): Table S2 has one interaction per row, from a source
component to a target component, and Table S1 describes the components.
:func:`join_components` adds the Table S1 record of both components to each
row of Table S2.

The adhesome.org website's ``interactions.csv`` (read with :func:`iter_csv`)
lists every interaction twice, once in each orientation.
:func:`iter_website_pairs` merges the two rows of a pair and leaves out the
rows of the 2007 network, which are read from the publication.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_csv, iter_excel

# The interaction column of Table S2, whose name is its own key.
INTERACTION = 'Interaction (1=activating; 2=inhibiting; 3=binding)'
INTERACTION_COLUMNS = ('Source', INTERACTION, 'target', 'PMID')
# Columns of Table S1 in source order (``Clasification`` is the source's
# spelling). The parsed row carries them for both components as
# ``source_<field>`` and ``target_<field>``.
COMPONENT_FIELDS = {
    'Protein name': 'name',
    'Synonyms': 'synonyms',
    'Gene name': 'gene',
    'swiss prot ID': 'uniprot',
    'Function': 'function',
    'Category': 'category',
    'Clasification': 'classification',
}

WEBSITE_COLUMNS = ('Source', 'Target', 'Effect', 'Type', 'PMID', 'data_source')
# ``data_source`` of the website's rows of the 2007 network.
CURATED_SOURCE = 'adhesome.org'

_uniprot_re = re.compile(
    r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}'
)


def _text(value: Any) -> str | None:
    """A cell as text; a spreadsheet number without its ``.0``."""
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    return str(value).strip() or None if value is not None else None


def _require(row: dict[str, Any], columns: Iterable[str], table: str) -> None:
    missing = [column for column in columns if column not in row]
    if missing:
        raise ValueError(
            f'Adhesome {table} lacks expected columns: {", ".join(missing)}'
        )


def uniprot_accessions(cell: str | None) -> list[str]:
    """The well-formed UniProt accessions of a ``swiss prot ID`` cell.

    A cell holds one accession or a comma separated list of them; a token that
    is not an accession is left out.
    """
    return [
        token.strip()
        for token in (cell or '').split(',')
        if _uniprot_re.fullmatch(token.strip())
    ]


def join_components(
    interactions: Iterable[dict[str, Any]],
    components: Iterable[dict[str, Any]],
) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of Table S2 with the Table S1 record of both components.

    Table S2 writes the component names in capitals (``PAXILLIN``) and Table S1
    as in the text (``paxillin``); they are matched ignoring case. A name
    Table S1 does not have leaves the component's fields None.
    """
    by_name: dict[str, dict[str, Any]] = {}
    for component in components:
        _require(component, COMPONENT_FIELDS, 'component table')
        name = _text(component['Protein name'])
        if not name:
            continue
        if name.lower() in by_name:
            raise ValueError(f'Adhesome component listed twice: {name!r}')
        by_name[name.lower()] = component
    if not by_name:
        raise ValueError('Adhesome component table has no rows')
    count = 0
    for row in interactions:
        _require(row, INTERACTION_COLUMNS, 'interaction table')
        count += 1
        out = {
            'source': _text(row['Source']),
            'target': _text(row['target']),
            'interaction': _text(row[INTERACTION]),
            'pmid': _text(row['PMID']),
        }
        for role in ('source', 'target'):
            component = by_name.get((out[role] or '').lower(), {})
            for column, field in COMPONENT_FIELDS.items():
                out[f'{role}_{field}'] = _text(component.get(column))
        yield out
    if not count:
        raise ValueError('Adhesome interaction table has no rows')


def iter_interactions(
    opener, force_refresh: bool = False, components_download=None, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield the interactions of Table S2, each with its components of Table S1.

    ``components_download`` is the download of Table S1.
    """
    if components_download is None:
        raise ValueError('Adhesome: the download of the component table is needed')
    components = iter_excel(components_download.open(force_refresh=force_refresh))
    yield from join_components(iter_excel(opener), components)


def iter_website_pairs(
    opener, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield the website's interactions that are not of the 2007 network.

    The two mirrored rows of an interaction (same components, same PubMed ID)
    become one row: ``component_a`` and ``component_b`` in alphabetical order,
    ``PMID``, and the distinct ``Type``, ``Effect`` and ``data_source`` values
    of the rows as sorted lists.
    """
    pairs: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    count = 0
    for row in iter_csv(opener):
        _require(row, WEBSITE_COLUMNS, 'website table')
        count += 1
        if row['data_source'] == CURATED_SOURCE:
            continue
        key = (*sorted((row['Source'].strip(), row['Target'].strip())), row['PMID'].strip())
        pairs.setdefault(key, []).append(row)
    if not count:
        raise ValueError('Adhesome website table has no rows')
    for (component_a, component_b, pmid), rows in pairs.items():
        yield {
            'component_a': component_a,
            'component_b': component_b,
            'PMID': pmid,
            'types': sorted({row['Type'] for row in rows}),
            'effects': sorted({row['Effect'] for row in rows}),
            'data_sources': sorted({row['data_source'] for row in rows}),
        }
