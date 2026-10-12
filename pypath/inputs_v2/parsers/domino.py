"""
Reader for the DOMINO interaction file.

DOMINO was distributed as a tab-separated file that starts like PSI-MI TAB:
its first 15 columns are those of MITAB 2.5, with PSI-MI terms written
without the ``psi-mi`` database name (``MI:0407(direct interaction)``). The
25 columns that follow are DOMINO's own and are not the columns of MITAB 2.6
or later. The header row names two columns ``interactor type B``, the first
of which is the type of interactor A, so the rows cannot be keyed by it:
:func:`iter_domino` reads the file with the shared table reader and gives
every column a field name of its own.

The cells of the columns that DOMINO added have their own syntax, read by
:func:`parse_binding_site` and :func:`parse_ptms`.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['FIELDS', 'iter_domino', 'parse_binding_site', 'parse_ptms']

# Header of the file and the field each column is read into, in file order.
COLUMNS = (
    ('ID interactor A', 'id_a'),
    ('ID interactor B', 'id_b'),
    ('Alt. ID interactor A', 'alt_ids_a'),
    ('Alt. ID interactor B', 'alt_ids_b'),
    ('Alias(es) interactor A', 'aliases_a'),
    ('Alias(es) interactor B', 'aliases_b'),
    ('Interaction detection method(s)', 'detection_methods'),
    ('Publication 1st author(s)', 'first_authors'),
    ('Publication Identifier(s)', 'publications'),
    ('Taxid interactor A', 'taxid_a'),
    ('Taxid interactor B', 'taxid_b'),
    ('Interaction type(s)', 'interaction_types'),
    ('Source database(s)', 'source_databases'),
    ('Interaction identifier(s)', 'interaction_ids'),
    ('Confidence value(s)', 'confidence'),
    ('expansion', 'expansion'),
    ('biological role A', 'biological_role_a'),
    ('biological role B', 'biological_role_b'),
    ('experimental role A', 'experimental_role_a'),
    ('experimental role B', 'experimental_role_b'),
    # The header names this column ``interactor type B`` too.
    ('interactor type B', 'interactor_type_a'),
    ('interactor type B', 'interactor_type_b'),
    ('xrefs A', 'xrefs_a'),
    ('xrefs B', 'xrefs_b'),
    ('xrefs Interaction', 'xrefs_interaction'),
    ('Annotations A', 'annotations_a'),
    ('Annotations B', 'annotations_b'),
    ('Interaction Annotations', 'annotations_interaction'),
    ('Host organism taxid', 'host_organism'),
    ('parameters Interaction', 'parameters'),
    ('dataset', 'dataset'),
    ('Caution Interaction', 'caution'),
    ('binding site A', 'binding_site_a'),
    ('binding site B', 'binding_site_b'),
    ('ptm A', 'ptm_a'),
    ('ptm B', 'ptm_b'),
    ('mutations A', 'mutations_a'),
    ('mutations B', 'mutations_b'),
    ('negative', 'negative'),
    ('inference', 'inference'),
)
_HEADER = tuple(header for header, _ in COLUMNS)
FIELDS = tuple(field for _, field in COLUMNS)


def iter_domino(opener, **_kwargs: Any) -> Generator[dict[str, str], None, None]:
    """Yield the rows of the DOMINO file as dictionaries keyed by :data:`FIELDS`.

    Values are the cells as written; ``-`` and the empty string both occur for
    a missing value. The header must be the 40 column names of the DOMINO
    release and every row must have 40 cells (one trailing empty cell is
    ignored): anything else raises ``ValueError``, as does a file without
    rows.
    """
    count = 0
    # The header cannot name the columns, so it is read as the first row.
    # Quotes are part of the cells.
    rows = iter_tsv(opener, fieldnames=FIELDS, quotechar=None)
    for number, row in enumerate(rows, 1):
        cells = [cell for cell in row.values() if isinstance(cell, str)]
        cells += row.pop(None, [])
        if not any(cell.strip() for cell in cells):
            continue
        if len(cells) == len(COLUMNS) + 1 and not cells[-1].strip():
            cells.pop()
        if number == 1:
            if tuple(cells) != _HEADER:
                raise ValueError(
                    'DOMINO file does not start with the expected header row'
                )
            continue
        if len(cells) != len(COLUMNS):
            raise ValueError(
                f'DOMINO row {number} has {len(cells)} cells, '
                f'expected {len(COLUMNS)}'
            )
        count += 1
        yield row
    if not count:
        raise ValueError('DOMINO file has no rows')


_RANGE_RE = re.compile(r'(\d+)-(\d+)')
_TERM_RE = re.compile(r'(MI:\d{4})\((.*)\)')
_INTERPRO_RE = re.compile(r'interpro:(IPR\d{6})(?:\((.*)\))?')
# The note of a modification that has a window of the protein sequence around
# it: ``range MINT-1497979: original sequence: 998-998 : HLPSPTDSNFYRALMDEEDM``.
# Where the modification is a prerequisite of the interaction, the source
# writes ``prerequisite-ptm`` straight after the window.
_WINDOW_RE = re.compile(
    r'range MINT-\d+: original sequence: (\d+)-(\d+) : ([A-Z]+)(?:prerequisite-ptm)?'
)


def _term(text: str) -> tuple[str | None, str | None]:
    """Accession and name of ``MI:0033(deletion analysis)``."""
    match = _TERM_RE.fullmatch(text.strip())
    if not match:
        return None, None
    return match.group(1), match.group(2).strip().strip('"') or None


def _positions(text: str) -> tuple[int | None, int | None]:
    """Start and end of ``12-45`` or of the single position ``12``.

    DOMINO writes 0 for a position that is not known; that and a range whose
    end precedes its start give ``(None, None)``.
    """
    text = text.strip()
    match = _RANGE_RE.fullmatch(text)
    if match:
        start, end = int(match.group(1)), int(match.group(2))
    elif text.isdigit():
        start = end = int(text)
    else:
        return None, None
    return (start, end) if 0 < start <= end else (None, None)


def parse_binding_site(cell: object) -> dict[str, Any] | None:
    """Read a ``binding site`` cell: ``ranges;method;domains``.

    - ``ranges``: residue ranges separated by commas (``15-92,168-244``),
      ``0`` when the position is not known, or ``-``;
    - ``method``: the PSI-MI term of the method that located the site
      (``MI:0033(deletion analysis)``), or ``-``;
    - ``domains``: InterPro entries separated by commas, with or without the
      entry's short name (``interpro:IPR001452(SH3)``), or ``-``.

    Returns None for ``-`` and for a cell that does not have three parts.
    Otherwise a dictionary with

    - ``ranges``: ``(start, end)`` pairs of the ranges that are valid
      one-based positions, in file order;
    - ``unknown_ranges``: the range entries that are not (``0``, ``351-349``);
    - ``method``: the method's PSI-MI accession, or None;
    - ``domains``: ``(accession, name)`` pairs of the InterPro entries;
    - ``unknown_domains``: domain entries that are not InterPro accessions
      (``interpro:CSD(CSD)``).
    """
    text = str(cell if cell is not None else '').strip()
    if not text or text == '-':
        return None
    parts = text.split(';')
    if len(parts) != 3:
        return None
    ranges, method, domains = (part.strip() for part in parts)
    site: dict[str, Any] = {
        'ranges': [],
        'unknown_ranges': [],
        'method': _term(method)[0],
        'domains': [],
        'unknown_domains': [],
    }
    for entry in ranges.split(',') if ranges != '-' else ():
        start, end = _positions(entry)
        if start is None:
            site['unknown_ranges'].append(entry.strip())
        else:
            site['ranges'].append((start, end))
    for entry in domains.split(',') if domains != '-' else ():
        match = _INTERPRO_RE.fullmatch(entry.strip())
        if match:
            site['domains'].append((match.group(1), match.group(2) or None))
        else:
            site['unknown_domains'].append(entry.strip())
    return site


def parse_ptms(cell: object) -> list[dict[str, Any]]:
    """Read a ``ptm`` cell: entries ``position;type;note;method`` joined by ``|``.

    - ``position``: the modified residue (``137``), a range, or ``0`` when the
      position is not known;
    - ``type``: the PSI-MI term of the modification
      (``MI:0176(o-phospho-serine)``);
    - ``note``: ``prerequisite-ptm``, a note with the sequence around the
      residue, or ``-``;
    - ``method``: the PSI-MI term of the method that showed it, or ``-``.

    Returns one dictionary per entry that has these four parts and a PSI-MI
    term as its type, with ``start`` and ``end`` (None when the position is
    not known), ``term`` and ``name`` of the modification, ``note``, the
    ``method`` accession and ``window``: None, or for a note with a window of
    the sequence the dictionary ``start`` and ``end`` (the range the note
    gives) and ``sequence`` (the window's residues, as written).
    """
    text = str(cell if cell is not None else '').strip()
    if not text or text == '-':
        return []
    modifications = []
    for entry in text.split('|'):
        parts = entry.split(';')
        if len(parts) != 4:
            continue
        position, modification, note, method = (part.strip() for part in parts)
        term, name = _term(modification)
        if term is None:
            continue
        start, end = _positions(position)
        window = _WINDOW_RE.fullmatch(note)
        modifications.append(
            {
                'start': start,
                'end': end,
                'term': term,
                'name': name,
                'note': None if note == '-' else note,
                'method': _term(method)[0],
                'window': {
                    'start': int(window.group(1)),
                    'end': int(window.group(2)),
                    'sequence': window.group(3),
                }
                if window
                else None,
            }
        )
    return modifications
