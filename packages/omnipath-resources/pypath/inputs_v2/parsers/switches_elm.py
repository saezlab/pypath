"""
Read the switches.ELM table of motif-based molecular switches.

Source file (``switches.ELM-v1.txt``, tab separated, with a header; lines end
with a bare carriage return; one line is one interaction of one switch):

    Switch ID  Status  Interaction ID  Intramolecular  ID A  Bindingsite A ID
    Bindingsite A Start  Bindingsite A End  ID B  Bindingsite B ID
    Bindingsite B Start  Bindingsite B End  Affected interactor  Switch type
    Switch subtype  Switch mechanism  Switch direction  Switch outcome direction
    Switch outcome  Modification  Modification sites  Modifying enzymes
    Effector  Cell cycle phase  Localisation  Pathway  PMID
    SWTI000001  Active  INTI000001    UNIPROT:O43561  ELM:LIG_SH2_STAT5  161  164
    UNIPROT:P19174  PFAM:PF00017  550  639  ID A  Binary
    Physicochemical compatibility  PTM-dependent  Reversible  Positive  Induction
    MOD:00696  Y161    ...  PMID: 20610546;PMID: 15065860

An interaction of a switch is between two participants. The participant A has
the motif (an ELM class) and B the domain that binds it (a Pfam or InterPro
entry), but the participants can also be a small molecule, a nucleic acid, a
membrane or a complex that an earlier interaction formed (``INTI000326``).
``Affected interactor`` says which participant the modification is added to or
the effector binds.
"""

from __future__ import annotations

from collections.abc import Generator
import io
import re
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv, read_text

COLUMNS = (
    'Switch ID',
    'Status',
    'Interaction ID',
    'Intramolecular',
    'ID A',
    'Bindingsite A ID',
    'Bindingsite A Start',
    'Bindingsite A End',
    'ID B',
    'Bindingsite B ID',
    'Bindingsite B Start',
    'Bindingsite B End',
    'Affected interactor',
    'Switch type',
    'Switch subtype',
    'Switch mechanism',
    'Switch direction',
    'Switch outcome direction',
    'Switch outcome',
    'Modification',
    'Modification sites',
    'Modifying enzymes',
    'Effector',
    'Cell cycle phase',
    'Localisation',
    'Pathway',
    'PMID',
)
# A UniProt accession, with the isoform number some participants have.
_UNIPROT = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(?:-\d+)?'
)
_INTERACTION = re.compile(r'INTI\d+')
_SITE = re.compile(r'([A-Z])([1-9]\d*)')


def _cells(value: object) -> list[str]:
    """The items of a cell with items joined by semicolons, without blanks."""
    return [item.strip() for item in str(value or '').split(';') if item.strip()]


def _participant(cell: object) -> dict[str, str | None]:
    """The participant of ``ID A`` or ``ID B``.

    ``kind`` is ``uniprot`` (``id`` is the accession, with its isoform
    number), ``interaction`` (a complex formed by that interaction of the
    table), or the lowercase prefix of any other identifier (``chebi``,
    ``ena``, ``go``), with ``id`` the rest.
    """
    text = str(cell or '').strip()
    prefix, _, rest = text.partition(':')
    if prefix == 'UNIPROT' and _UNIPROT.fullmatch(rest):
        return {'kind': 'uniprot', 'id': rest}
    if _INTERACTION.fullmatch(text):
        return {'kind': 'interaction', 'id': text}
    if rest:
        return {'kind': prefix.lower(), 'id': rest}
    return {'kind': None, 'id': None}


def _binding_sites(row: dict[str, Any], side: str) -> list[dict[str, Any]]:
    """The binding sites of one participant: ``id`` (as written, e.g.
    ``ELM:LIG_SH2_STAT5``, None without one), ``start`` and ``end``.

    The cells hold the identifiers and the positions joined by semicolons. A
    single identifier with several positions (a domain that occurs several
    times) is repeated for each; otherwise the items pair by position.
    """
    ids = _cells(row[f'Bindingsite {side} ID'])
    starts = [int(x) for x in _cells(row[f'Bindingsite {side} Start'])]
    ends = [int(x) for x in _cells(row[f'Bindingsite {side} End'])]
    if len(starts) != len(ends):
        raise ValueError(
            f'Switch {row["Switch ID"]}: {len(starts)} start and {len(ends)} '
            f'end positions of binding site {side}'
        )
    if len(ids) == 1 and len(starts) > 1:
        ids = ids * len(starts)
    if starts and len(ids) != len(starts):
        raise ValueError(
            f'Switch {row["Switch ID"]}: {len(ids)} binding site '
            f'identifiers and {len(starts)} positions on {side}'
        )
    if not starts:
        return [{'id': item, 'start': None, 'end': None} for item in ids]
    return [
        {'id': item, 'start': start, 'end': end}
        for item, start, end in zip(ids, starts, ends)
    ]


def _modification_sites(cell: object) -> list[dict[str, Any]]:
    """The modified residues (``Y161``) as ``residue`` and ``position``."""
    sites = []
    for item in _cells(cell):
        match = _SITE.fullmatch(item)
        if not match:
            raise ValueError(f'Unreadable modification site: {item!r}')
        sites.append({'residue': match[1], 'position': int(match[2])})
    return sites


def _row(row: dict[str, Any]) -> dict[str, Any]:
    """The decoded fields of one line.

    Keys: the columns of the table that are text as written
    (``switch_id``, ``interaction_id``, ``status``, ``switch_type``,
    ``switch_subtype``, ``switch_mechanism``, ``switch_direction``,
    ``outcome_direction``, ``outcome``, ``modification``, ``effector``,
    ``modifying_enzymes``, ``cell_cycle_phase``, ``localisation``,
    ``pathway``, with empty cells as None), ``intramolecular`` (bool),
    ``a`` and ``b`` (the participants, see ``_participant``, each with its
    ``binding_sites``), ``affected`` (``A``, ``B`` or None),
    ``modification_sites`` and ``pmids`` (the items of the reference cell,
    without the prefix ``PMID:``).
    """
    def text(column: str) -> str | None:
        return str(row.get(column) or '').strip() or None

    affected = text('Affected interactor')
    if affected not in (None, 'ID A', 'ID B'):
        raise ValueError(f'Unknown affected interactor: {affected!r}')
    intramolecular = text('Intramolecular')
    if intramolecular not in (None, 'TRUE'):
        raise ValueError(f'Unknown value of Intramolecular: {intramolecular!r}')
    return {
        'switch_id': text('Switch ID'),
        'interaction_id': text('Interaction ID'),
        'status': text('Status'),
        'intramolecular': intramolecular == 'TRUE',
        'a': {**_participant(row['ID A']), 'binding_sites': _binding_sites(row, 'A')},
        'b': {**_participant(row['ID B']), 'binding_sites': _binding_sites(row, 'B')},
        'affected': affected[-1] if affected else None,
        'switch_type': text('Switch type'),
        'switch_subtype': text('Switch subtype'),
        'switch_mechanism': text('Switch mechanism'),
        'switch_direction': text('Switch direction'),
        'outcome_direction': text('Switch outcome direction'),
        'outcome': text('Switch outcome'),
        'modification': text('Modification'),
        'modification_sites': _modification_sites(row['Modification sites']),
        'modifying_enzymes': text('Modifying enzymes'),
        'effector': text('Effector'),
        'cell_cycle_phase': text('Cell cycle phase'),
        'localisation': text('Localisation'),
        'pathway': text('Pathway'),
        'pmids': [
            item.removeprefix('PMID:').strip() for item in _cells(row.get('PMID'))
        ],
    }


def iter_switches(opener, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    Yield the interactions of the switches.ELM table, decoded (see ``_row``).

    The file ends its lines with a bare carriage return, which is changed to
    a newline before the table is read. A header without one of the columns
    of ``COLUMNS``, a line with another number of fields, a binding site
    whose identifiers and positions do not pair, an unreadable modification
    site and a file without lines are errors.
    """
    text = read_text(opener).replace('\r\n', '\n').replace('\r', '\n')
    lines = iter_tsv(
        SimpleNamespace(result=io.StringIO(text)),
        strict=True,
        quotechar=None,
        required=COLUMNS,
        **kwargs,
    )
    for line in lines:
        yield _row(line)
