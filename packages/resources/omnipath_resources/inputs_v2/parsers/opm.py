"""Read the Orientations of Proteins in Membranes (OPM) database exports.

OPM serves its tables as CSV files from its backend
(``https://opm-back.cc.lehigh.edu/opm-backend/<table>?fileFormat=csv``), the
same files the "download" buttons of its browsing pages give. Two groups of
tables are read:

* the vocabularies: the membranes (24, with numeric ids) and the four levels
  of the classification of the structures, each with its own numeric ids:
  types (3), classes (11, the ``classtypes`` table), superfamilies (530) and
  families (1090). A class names its type, a superfamily its class and a family
  its superfamily;
* the primary structures (the ``primary_structures`` table, 8915 rows): one
  row per PDB entry, with its membrane, family, the UniProt codes of the
  proteins of the entry and the calculated orientation in the membrane.

:func:`iter_terms` yields the vocabularies as rows of one kind each,
:func:`iter_structures` the structures. Neither knows the schema.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.base import Download
from pypath.inputs_v2.parsers.base import iter_csv

__all__ = ['iter_structures', 'iter_terms']

# The columns each table is read with; a missing one is an error.
MEMBRANE_COLUMNS = ('id', 'name', 'short_name', 'topology_in', 'topology_out')
TYPE_COLUMNS = ('id', 'name')
CLASS_COLUMNS = ('id', 'name', 'type_id')
SUPERFAMILY_COLUMNS = ('id', 'name', 'comment', 'tcdb', 'pfam', 'classtype_id')
FAMILY_COLUMNS = (
    'id', 'name', 'tcdb', 'pfam', 'interpro', 'superfamily_id',
)
STRUCTURE_COLUMNS = (
    'id', 'pdbid', 'name', 'uniprotcode', 'membrane_id', 'family_id',
    'thickness', 'thicknesserror', 'tilt', 'tilterror', 'gibbs', 'tau',
)

# The cross-references of the classification, by the syntax of the accession.
# A Pfam entry (``PF00556``) or clan (``CL0328``), an InterPro entry, a TCDB
# class (``1.G``), subclass, family or subfamily number (``5.B.4``, ``1.A.1.19``).
_XREF_PATTERNS = {
    'pfam': re.compile(r'PF\d{5}|CL\d{4}'),
    'interpro': re.compile(r'IPR\d{6}'),
    'tcdb': re.compile(r'\d\.[A-Z](?:\.\d+){0,2}'),
}
# The export is written for a spreadsheet: a PDB code is a formula
# (``="1qjp"``) so that it is not read as a number.
_PDB_FORMULA = re.compile(r'="([0-9][0-9a-z]{3})"')


def _clean(value: Any) -> str | None:
    """A cell without the white space and the byte order mark around it."""
    text = str(value if value is not None else '').strip('﻿ \t\r\n')
    return text or None


def _xrefs(row: dict[str, Any]) -> dict[str, str]:
    """The cross-references of a row, by database.

    Cells of the classification tables have a trailing space (a TCDB number) or
    a byte order mark (a Pfam entry), which are removed, and two TCDB numbers
    end with a dot (``1.B.42.``), which is removed. One accession has the wrong
    syntax (the Pfam entry ``PF0036`` has one digit too few) and is left out:
    it is not repaired.
    """
    refs = {}
    for database, pattern in _XREF_PATTERNS.items():
        value = _clean(row.get(database))
        if value and database == 'tcdb':
            value = value.rstrip('.')
        if value and pattern.fullmatch(value):
            refs[database] = value
    return refs


def _rows(download: Download, columns: tuple[str, ...], **kwargs: Any):
    opener = download.open(**kwargs)
    return iter_csv(opener, required=columns)


def iter_terms(
    opener,
    *,
    types: Download,
    classes: Download,
    superfamilies: Download,
    families: Download,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield the membranes, then the types, classes, superfamilies and families.

    ``opener`` is the download of the membranes table; the other tables are the
    downloads given as arguments. Every row has ``kind`` (``membrane``,
    ``type``, ``class``, ``superfamily`` or ``family``), ``id`` and ``name``;
    the others are:

    * membrane: ``short_name``, ``inside`` and ``outside``, the sides of the
      membrane as OPM names them (``topology_in``, ``topology_out``);
    * class, superfamily, family: ``parent_kind`` and ``parent_id``;
    * superfamily: ``comment``, the note of the three superfamilies that have
      one; superfamily and family: ``xrefs``, a dictionary of the
      cross-references the row has (``pfam``, ``interpro``, ``tcdb``).

    A class, superfamily or family whose parent is not in the tables is an
    error.
    """
    forced = {key: kwargs[key] for key in ('force_refresh',) if key in kwargs}
    seen: dict[str, set[str]] = {}

    def check(kind: str, identifier: str, name: str | None) -> None:
        if not identifier or not name:
            raise ValueError(f'OPM {kind} without an id or a name: {identifier!r}')
        if identifier in seen.setdefault(kind, set()):
            raise ValueError(f'OPM {kind} id {identifier} is repeated')
        seen[kind].add(identifier)

    for row in iter_csv(opener, required=MEMBRANE_COLUMNS):
        check('membrane', row['id'], _clean(row['name']))
        yield {
            'kind': 'membrane',
            'id': row['id'],
            'name': _clean(row['name']),
            'short_name': _clean(row['short_name']),
            'inside': _clean(row['topology_in']),
            'outside': _clean(row['topology_out']),
        }
    for row in _rows(types, TYPE_COLUMNS, **forced):
        check('type', row['id'], _clean(row['name']))
        yield {'kind': 'type', 'id': row['id'], 'name': _clean(row['name'])}

    levels = (
        ('class', classes, CLASS_COLUMNS, 'type', 'type_id'),
        ('superfamily', superfamilies, SUPERFAMILY_COLUMNS, 'class', 'classtype_id'),
        ('family', families, FAMILY_COLUMNS, 'superfamily', 'superfamily_id'),
    )
    for kind, download, columns, parent_kind, parent_column in levels:
        for row in _rows(download, columns, **forced):
            check(kind, row['id'], _clean(row['name']))
            if row[parent_column] not in seen[parent_kind]:
                raise ValueError(
                    f'OPM {kind} {row["id"]} is in {parent_kind} '
                    f'{row[parent_column]!r}, which is not in the {parent_kind} table'
                )
            term = {
                'kind': kind,
                'id': row['id'],
                'name': _clean(row['name']),
                'parent_kind': parent_kind,
                'parent_id': row[parent_column],
            }
            if kind != 'class':
                term['xrefs'] = _xrefs(row)
            if kind == 'superfamily' and _clean(row['comment']):
                term['comment'] = _clean(row['comment'])
            yield term


def iter_structures(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the primary structures, one per row of the table.

    The fields:

    * ``id``: the OPM record number (two records can have the same PDB code);
    * ``pdb_id``: the PDB code, from the spreadsheet formula the export writes;
    * ``name``: the name of the structure;
    * ``entries``: the UniProt codes of the proteins of the entry (entry names
      such as ``OMPA_ECOLI``, a few accessions), in the order of the export
      and without repeats, empty where OPM states none;
    * ``membrane_id`` and ``family_id``: ids of the membrane and the family
      tables;
    * ``thickness``, ``thickness_error``, ``tilt``, ``tilt_error``,
      ``transfer_energy``, ``curvature_radius``: the calculated orientation of
      the entry in the membrane, as written (``None`` where empty).
    """
    for row in iter_csv(opener, required=STRUCTURE_COLUMNS):
        formula = _PDB_FORMULA.fullmatch(row['pdbid'] or '')
        if not formula:
            raise ValueError(
                f'OPM structure {row["id"]}: unexpected PDB code {row["pdbid"]!r}'
            )
        if not row['membrane_id'] or not row['family_id']:
            raise ValueError(f'OPM structure {row["id"]} has no membrane or family')
        yield {
            'id': row['id'],
            'pdb_id': formula.group(1),
            'name': _clean(row['name']),
            'entries': list(dict.fromkeys((row['uniprotcode'] or '').split())),
            'membrane_id': row['membrane_id'],
            'family_id': row['family_id'],
            'thickness': _clean(row['thickness']),
            'thickness_error': _clean(row['thicknesserror']),
            'tilt': _clean(row['tilt']),
            'tilt_error': _clean(row['tilterror']),
            'transfer_energy': _clean(row['gibbs']),
            'curvature_radius': _clean(row['tau']),
        }

