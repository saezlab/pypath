"""Read the flat files of the Human Protein Reference Database (HPRD).

HPRD is a database of human proteins (its name; the flat files carry no species
column). The release read here is the archive ``HPRD_FLAT_FILES_041310.tar.gz``,
whose members are in the directory ``FLAT_FILES_072010/``. Each table is a
tab-separated file with no header row; the column names are those of the
``README`` of the archive. A missing value is written ``-``.

Three tables are read:

- ``POST_TRANSLATIONAL_MODIFICATIONS.txt``: one row per modification site of a
  substrate isoform, with the enzyme where HPRD names one.
  :func:`iter_modification_sites` yields every row, :func:`iter_enzyme_sites`
  those that name an enzyme.
- ``BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt``: one row per interacting pair
  (:func:`iter_interactions`).
- ``PROTEIN_COMPLEXES.txt``: one row per complex member, the complex named
  ``COM_<number>`` (:func:`iter_complexes`, which groups them).

Two further tables of the archive are not read: ``BINARY_PROTEIN_NONPROTEIN_
INTERACTIONS.txt``, whose partners are free text, and ``HPRD_ID_MAPPINGS.txt``.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

ARCHIVE_DIRECTORY = 'FLAT_FILES_072010/'
PTM_MEMBER = ARCHIVE_DIRECTORY + 'POST_TRANSLATIONAL_MODIFICATIONS.txt'
INTERACTION_MEMBER = ARCHIVE_DIRECTORY + 'BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt'
COMPLEX_MEMBER = ARCHIVE_DIRECTORY + 'PROTEIN_COMPLEXES.txt'

# Column names of the archive's README, in file order.
PTM_COLUMNS = (
    'substrate_hprd_id',
    'substrate_gene_symbol',
    'substrate_isoform_id',
    'substrate_refseq_id',
    'site',
    'residue',
    'enzyme_name',
    'enzyme_hprd_id',
    'modification_type',
    'experiment_type',
    'reference_id',
)
INTERACTION_COLUMNS = (
    'interactor_1_geneSymbol',
    'interactor_1_hprd_id',
    'interactor_1_refseq_id',
    'interactor_2_geneSymbol',
    'interactor_2_hprd_id',
    'interactor_2_refseq_id',
    'experiment_type',
    'reference_id',
)
COMPLEX_COLUMNS = (
    'hprd_interaction_id',
    'interactor_hprd_id',
    'interactor_geneSymbol',
    'interactor_refseq_id',
    'experiment_type',
    'reference_id',
)

_HPRD_ID_REGEX = re.compile(r'^\d{5}$')
_ISOFORM_ID_REGEX = re.compile(r'^\d{5}_\d+$')
_COMPLEX_ID_REGEX = re.compile(r'^COM_[1-9]\d*$')
# "12": one position; "12;-": a position and an empty second one; "416;462":
# two positions (the cysteines of a disulfide bridge). A few sites are "0",
# which is no position of a sequence.
_SITE_REGEX = re.compile(r'^(\d+)(?:;(-|\d+))?$')
_RESIDUE_REGEX = re.compile(r'^[A-Z]$')
_UNIPROT = (
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
)
# The "refseq" columns hold the sequence of the isoform in one of four kinds
# of accession, not only RefSeq's.
_ACCESSION_KINDS = (
    ('refseq_protein', re.compile(r'^(?:AP|NP|XP|YP|WP)_\d+\.\d+$')),
    ('genbank', re.compile(r'^[A-Z]{3}\d{5,7}\.\d+$')),
    ('uniprot', re.compile(rf'^{_UNIPROT}$')),
    # Three rows of one protein (13685, CCDC56) have a GenBank nucleotide record.
    ('genbank_nucleotide', re.compile(r'^[A-Z]{1,2}\d{5,6}\.\d+$')),
)
# Spellings of the experiment column that state no experiment type.
_NO_EXPERIMENT_TYPE = frozenset({'', 'None'})


def _rows(
    opener, member: str, columns: tuple[str, ...]
) -> Generator[dict[str, Any], None, None]:
    """Rows of one member of the archive, cells stripped, with their line number.

    The member is picked by its exact name; a download without it is an error,
    as is a line with another number of fields than the table has columns.
    """
    result = getattr(opener, 'result', None)
    if not isinstance(result, dict) or member not in result:
        raise ValueError(f'HPRD archive has no member {member}')
    table = SimpleNamespace(result=result[member])
    count = 0
    # A quote character in the files would be data, but none occurs.
    for count, row in enumerate(
        iter_tsv(table, fieldnames=columns, strict=True, quotechar=None), 1
    ):
        cells = {name: value.strip() for name, value in row.items()}
        cells['line'] = count
        yield cells
    if not count:
        raise ValueError(f'HPRD table {member} has no rows')


def _missing(value: str) -> str | None:
    """The source writes a missing value as a hyphen."""
    return None if value in ('', '-') else value


def _hprd_id(value: str, where: str) -> str:
    if not _HPRD_ID_REGEX.match(value):
        raise ValueError(f'HPRD identifier not understood in {where}: {value!r}')
    return value


def accession_kind(value: str, where: str) -> str:
    """The kind of the accession in a "refseq" column: ``refseq_protein`` for
    ``NP_000680.2``, ``genbank`` for ``AAH11772.1``, ``uniprot`` for ``Q14643``,
    ``genbank_nucleotide`` for ``BC002698.2`` (a nucleotide record).

    Any other form is an error.
    """
    for kind, regex in _ACCESSION_KINDS:
        if regex.match(value):
            return kind
    raise ValueError(f'HPRD sequence accession not understood in {where}: {value!r}')


def experiment_types(value: str) -> list[str]:
    """The experiment types of a cell, in order and without repeats.

    The cell lists types separated by ``;`` (``in vitro;in vivo``), at times
    with a space after the separator or a type repeated. An empty cell and
    the word ``None`` state no experiment type. The words are not interpreted.
    """
    types = [part.strip() for part in value.split(';')]
    return list(
        dict.fromkeys(part for part in types if part not in _NO_EXPERIMENT_TYPE)
    )


def references(value: str) -> list[str]:
    """The reference IDs of a cell, comma separated, as written.

    The source says nothing else of them; they are PubMed IDs where the cell
    holds numbers.
    """
    return [part.strip() for part in value.split(',') if part.strip()]


def iter_modification_sites(
    opener, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per row of the modification table.

    Fields:

    - ``line``: line number in the file.
    - ``substrate_hprd_id``, ``substrate_isoform_id``: HPRD's five-digit
      protein identifier and the identifier of the isoform (``00001_1``).
    - ``substrate_gene_symbol``: None where the source has ``-``.
    - ``substrate_accession``, ``substrate_accession_kind``: the accession of
      the isoform's sequence and its kind (see :func:`accession_kind`).
    - ``site``: the site cell as written; ``position`` its first number and
      ``second_position`` the second one, None where the cell has ``-`` in
      its place or no second part.
    - ``residue``: one-letter code of the modified residue.
    - ``enzyme_name``, ``enzyme_hprd_id``: the enzyme as HPRD names it (a gene
      symbol) and its HPRD identifier, None where the source has ``-``;
      ``has_enzyme`` is true if either is given.
    - ``modification_type``: HPRD's modification type as written.
    - ``experiment_type``: the cell as written; ``experiment_types`` the list
      of :func:`experiment_types`.
    - ``references``: the list of :func:`references`.
    """
    for row in _rows(opener, PTM_MEMBER, PTM_COLUMNS):
        where = f'{PTM_MEMBER} line {row["line"]}'
        site = _SITE_REGEX.match(row['site'])
        if not site:
            raise ValueError(f'HPRD site not understood in {where}: {row["site"]!r}')
        if not _RESIDUE_REGEX.match(row['residue']):
            raise ValueError(f'HPRD residue not understood in {where}: {row["residue"]!r}')
        if not row['modification_type']:
            raise ValueError(f'HPRD modification type missing in {where}')
        isoform = row['substrate_isoform_id']
        if not _ISOFORM_ID_REGEX.match(isoform):
            raise ValueError(f'HPRD isoform identifier not understood in {where}: {isoform!r}')
        enzyme_name = _missing(row['enzyme_name'])
        enzyme_id = _missing(row['enzyme_hprd_id'])
        if enzyme_id is not None:
            _hprd_id(enzyme_id, where)
        accession = row['substrate_refseq_id']
        yield {
            'line': row['line'],
            'substrate_hprd_id': _hprd_id(row['substrate_hprd_id'], where),
            'substrate_gene_symbol': _missing(row['substrate_gene_symbol']),
            'substrate_isoform_id': isoform,
            'substrate_accession': accession,
            'substrate_accession_kind': accession_kind(accession, where),
            'site': row['site'],
            'position': int(site[1]),
            'second_position': (
                int(site[2]) if site[2] and site[2] != '-' else None
            ),
            'residue': row['residue'],
            'enzyme_name': enzyme_name,
            'enzyme_hprd_id': enzyme_id,
            'has_enzyme': enzyme_name is not None or enzyme_id is not None,
            'modification_type': row['modification_type'],
            'experiment_type': row['experiment_type'],
            'experiment_types': experiment_types(row['experiment_type']),
            'references': references(row['reference_id']),
        }


def iter_enzyme_sites(
    opener, **kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of the modification table that name an enzyme.

    The other rows are modification sites without an enzyme: they are not
    enzyme-substrate statements. The ``experiment_type`` is not used to select
    rows (the legacy reader kept only ``in vivo`` ones).
    """
    for row in iter_modification_sites(opener, **kwargs):
        if row['has_enzyme']:
            yield row


def _interactor(row: dict[str, str], number: int, where: str) -> dict[str, Any]:
    accession = row[f'interactor_{number}_refseq_id']
    return {
        f'interactor_{number}_gene_symbol': _missing(
            row[f'interactor_{number}_geneSymbol']
        ),
        f'interactor_{number}_hprd_id': _hprd_id(
            row[f'interactor_{number}_hprd_id'], where
        ),
        f'interactor_{number}_accession': accession,
        f'interactor_{number}_accession_kind': accession_kind(accession, where),
    }


def iter_interactions(
    opener, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per row of the binary protein-protein table.

    Fields: ``line``; for each of the two interactors, ``interactor_<n>_
    gene_symbol`` (None where the source has ``-``), ``interactor_<n>_hprd_id``,
    ``interactor_<n>_accession`` and ``interactor_<n>_accession_kind`` (as for
    the substrate of :func:`iter_modification_sites`); ``experiment_type`` as
    written and ``experiment_types`` (:func:`experiment_types`); ``references``.

    A pair can be listed in both orientations and a line can repeat; every
    line is yielded.
    """
    for row in _rows(opener, INTERACTION_MEMBER, INTERACTION_COLUMNS):
        where = f'{INTERACTION_MEMBER} line {row["line"]}'
        yield {
            'line': row['line'],
            **_interactor(row, 1, where),
            **_interactor(row, 2, where),
            'experiment_type': row['experiment_type'],
            'experiment_types': experiment_types(row['experiment_type']),
            'references': references(row['reference_id']),
        }


def iter_complexes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per complex of the complex table.

    The table has one row per member; the rows of a complex are together, and
    its experiment type and references are the same on each. A complex whose
    rows are not together, or differ in these cells, is an error.

    Fields:

    - ``complex_id``: ``COM_<number>``.
    - ``members``: the identified members, each a dictionary with ``hprd_id``,
      ``gene_symbol`` (None where the source has ``-``), ``accession`` and
      ``accession_kind`` (None where the source has ``-``; see
      :func:`accession_kind`). A member listed on several rows is given once;
      two rows of one HPRD identifier that differ are an error.
    - ``unidentified_rows``: the number of rows whose member has no HPRD
      identifier (the cell holds ``None``, the symbol and accession ``-``).
      They name no protein and are not in ``members``.
    - ``duplicate_rows``: the number of rows repeating the row of an identified
      member given before.
    - ``first_line``: line number of the complex's first row.
    - ``experiment_type`` as written and ``experiment_types``; ``references``.
    """
    current: dict[str, Any] | None = None
    seen: set[str] = set()

    def finish(complex_: dict[str, Any]) -> dict[str, Any]:
        complex_['members'] = list(complex_['members'].values())
        return complex_

    for row in _rows(opener, COMPLEX_MEMBER, COMPLEX_COLUMNS):
        where = f'{COMPLEX_MEMBER} line {row["line"]}'
        complex_id = row['hprd_interaction_id']
        if not _COMPLEX_ID_REGEX.match(complex_id):
            raise ValueError(f'HPRD complex identifier not understood in {where}: {complex_id!r}')
        if current is None or current['complex_id'] != complex_id:
            if current is not None:
                yield finish(current)
            if complex_id in seen:
                raise ValueError(f'HPRD complex {complex_id} is not in one block ({where})')
            seen.add(complex_id)
            current = {
                'complex_id': complex_id,
                'first_line': row['line'],
                'members': {},
                'unidentified_rows': 0,
                'duplicate_rows': 0,
                'experiment_type': row['experiment_type'],
                'experiment_types': experiment_types(row['experiment_type']),
                'references': references(row['reference_id']),
            }
        if (
            row['experiment_type'] != current['experiment_type']
            or references(row['reference_id']) != current['references']
        ):
            raise ValueError(f'HPRD complex {complex_id} differs in experiment or references ({where})')
        hprd_id = row['interactor_hprd_id']
        if hprd_id == 'None':
            if (row['interactor_geneSymbol'], row['interactor_refseq_id']) != ('-', '-'):
                raise ValueError(f'HPRD member without identifier has a symbol or accession ({where})')
            current['unidentified_rows'] += 1
            continue
        member = {
            'hprd_id': _hprd_id(hprd_id, where),
            'gene_symbol': _missing(row['interactor_geneSymbol']),
            'accession': _missing(row['interactor_refseq_id']),
        }
        member['accession_kind'] = (
            accession_kind(member['accession'], where) if member['accession'] else None
        )
        earlier = current['members'].get(hprd_id)
        if earlier is None:
            current['members'][hprd_id] = member
        elif earlier == member:
            current['duplicate_rows'] += 1
        else:
            raise ValueError(f'HPRD complex {complex_id} lists {hprd_id} twice, differently ({where})')
    if current is None:
        raise ValueError('HPRD complex table has no rows')
    yield finish(current)
