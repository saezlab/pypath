"""
Reader for the ELM interaction exports.

ELM publishes its motif-domain interactions twice, as two exports of the same
records in the same order:

- ``interactions/as_tsv``, a table with the motif class, the domain, the two
  proteins, the residues of the motif instance and of the domain, the measured
  affinity, the publications and the species. It has no record number;
- ``interactions/as_mitab``, a PSI-MI TAB 2.5 file with the record number
  (``ELMINT000120``), the detection methods and the gene names, but not the
  motif or the domain.

:func:`iter_interactions` reads both and joins them by position, as
:func:`join_evidence` explains. Many records share both proteins, both species
and all publications and differ only in the motif or the domain, so a join by
position could pair such records wrongly without any column disagreeing. A
third export, ``interactions/as_psimi`` (PSI-MI XML), lists the same records by
number with the motif and the domain of each; :func:`read_features` reads it
and :func:`join_evidence` checks every pair against it.
:func:`parse_position` and :func:`split_methods` read two cells in ELM's own
way.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Generator, Iterable, Mapping
from itertools import zip_longest
import re
from typing import Any

from pypath.inputs_v2 import _mitab
from pypath.inputs_v2.parsers.base import iter_tsv, read_xml_root

__all__ = [
    'COLUMNS',
    'EVIDENCE_COLUMNS',
    'iter_interactions',
    'join_evidence',
    'parse_position',
    'read_features',
    'split_methods',
]

# Columns of the interaction table, in file order. The header row ends with a
# tab, which the reader gives back as a column without a name.
COLUMNS = (
    'Elm',
    'Domain',
    'interactorElm',
    'interactorDomain',
    'StartElm',
    'StopElm',
    'StartDomain',
    'StopDomain',
    'AffinityMin',
    'AffinityMax',
    'PMID',
    'taxonomyElm',
    'taxonomyDomain',
)
# Columns of the PSI-MI TAB export that the joined rows keep. ELM names the
# columns in its own way; interactor A is the protein with the motif.
EVIDENCE_COLUMNS = (
    'Aliases for A',
    'Aliases for B',
    'Methods',
    'Interaction_types',
    'Interaction_identifier',
)
_EVIDENCE_KEY_COLUMNS = (
    'Unique_identifier_A',
    'Unique_identifier_B',
    'publication',
)


_TAXON_RE = re.compile(r'(?:taxid:)?"?(\d+)"?\(')
_RECORD_RE = re.compile(r'ELMINT\d+')


def _taxon(cell: object) -> str:
    """The NCBI taxon number of ``"559292"(Saccharomyces ...)`` or of
    ``taxid:559292(Saccharomyces ...)``, the first one of a list."""
    match = _TAXON_RE.match(str(cell or '').split('|')[0].strip())
    return match.group(1) if match else ''


def _key(evidence_row: dict[str, Any]) -> tuple[str, str, frozenset[str], str, str]:
    """Interactors, publications and species of a row of the PSI-MI TAB export."""
    first, second = (
        next(
            (
                item.accession
                for item in _mitab.parse_mitab_cell(evidence_row.get(column))
                if item.database == 'uniprotkb'
            ),
            '',
        )
        for column in _EVIDENCE_KEY_COLUMNS[:2]
    )
    publications = frozenset(
        item.accession
        for item in _mitab.parse_mitab_cell(evidence_row.get('publication'))
        if item.database == 'pubmed'
    )
    return (
        first,
        second,
        publications,
        _taxon(evidence_row.get('Taxon_A')),
        _taxon(evidence_row.get('Taxon_B')),
    )


Feature = tuple[str, str, int | None, int | None]


def read_features(opener) -> dict[str, Counter]:
    """The motif and the domain of each record of the PSI-MI XML export.

    Returns, for each record number (``ELMINT000007``), the features of its
    participants as a count of ``(UniProt accession, feature name, start,
    end)``. The feature name is the motif class or the domain accession, and
    a position the export leaves out, or writes 0, is None. A document that
    is not the export, a record without a number, a number that occurs twice
    and a feature of a participant without an interactor raise ``ValueError``.
    """
    root = read_xml_root(opener)
    if root is None or not root.tag.endswith('entrySet'):
        raise ValueError('The ELM PSI-MI XML export is missing or is not PSI-MI XML')
    features: dict[str, Counter] = {}
    for entry in root.iterfind('{*}entry'):
        accessions = {}
        for interactor in entry.iterfind('{*}interactorList/{*}interactor'):
            primary = interactor.find('{*}xref/{*}primaryRef')
            accessions[interactor.get('id')] = (
                primary.get('id') if primary is not None else None
            )
        for interaction in entry.iterfind('{*}interactionList/{*}interaction'):
            reference = interaction.find('{*}xref/{*}primaryRef')
            record = reference.get('id') if reference is not None else None
            if not _RECORD_RE.fullmatch(record or ''):
                raise ValueError(f'ELM PSI-MI XML interaction without a record number: {record!r}')
            if record in features:
                raise ValueError(f'ELM PSI-MI XML lists {record} twice')
            found: Counter = Counter()
            for participant in interaction.iterfind('{*}participantList/{*}participant'):
                interactor = participant.findtext('{*}interactorRef')
                accession = accessions.get((interactor or '').strip())
                for feature in participant.iterfind('{*}featureList/{*}feature'):
                    if not accession:
                        raise ValueError(f'{record}: a feature of a participant without an interactor')
                    begin = feature.find('{*}featureRangeList/{*}featureRange/{*}begin')
                    end = feature.find('{*}featureRangeList/{*}featureRange/{*}end')
                    found[
                        (
                            accession,
                            (feature.findtext('{*}names/{*}shortLabel') or '').strip(),
                            parse_position(begin.get('position')) if begin is not None else None,
                            parse_position(end.get('position')) if end is not None else None,
                        )
                    ] += 1
            features[record] = found
    if not features:
        raise ValueError('The ELM PSI-MI XML export has no interactions')
    return features


def join_evidence(
    interactions: Iterable[dict[str, Any]],
    evidence: Iterable[dict[str, Any]],
    features: Mapping[str, Counter],
) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of the interaction table with the fields of the PSI-MI
    TAB export of the same records.

    The two exports list the records in the same order and the table has no
    record number, so the rows are paired by position. A pair must agree on
    the two interactors, the two species and the publications. Records that
    share all of these differ in the motif or the domain, which the PSI-MI
    TAB export does not have, so each pair is also checked against
    ``features`` (see :func:`read_features`): the record of the pair, by its
    number, must have exactly the motif and the domain that the table row
    states, each with its protein and its positions. A disagreement, a record
    that ``features`` lacks, a table of another length and a missing column
    raise ``ValueError``.

    Fields of a yielded row, all strings as written in the exports:

    - the columns of :data:`COLUMNS` of the interaction table: ``Elm`` (the
      motif class), ``Domain`` (a Pfam, InterPro or SMART accession),
      ``interactorElm`` and ``interactorDomain`` (UniProt accessions of the
      protein with the motif and of the one with the domain, some with an
      isoform number), the residues of the motif instance (``StartElm``,
      ``StopElm``) and of the domain (``StartDomain``, ``StopDomain``;
      ``None`` or ``0`` when not given), the affinity (``AffinityMin``,
      ``AffinityMax``; ``None`` when not measured), ``PMID`` (PubMed IDs
      joined by commas) and the species (``taxonomyElm``,
      ``taxonomyDomain``, ``9606(Homo sapiens)``);
    - the columns of :data:`EVIDENCE_COLUMNS` of the PSI-MI TAB export: the
      gene names of the two proteins (``Aliases for A``, ``Aliases for B``),
      ``Methods`` (PSI-MI terms joined by semicolons), ``Interaction_types``
      and ``Interaction_identifier`` (ELM's record number).
    """
    count = 0
    for number, (row, evidence_row) in enumerate(
        zip_longest(interactions, evidence), 1
    ):
        if row is None or evidence_row is None:
            raise ValueError(
                'ELM interaction table and PSI-MI TAB export differ in length '
                f'(from row {number})'
            )
        missing = [column for column in COLUMNS if column not in row] + [
            column
            for column in EVIDENCE_COLUMNS + _EVIDENCE_KEY_COLUMNS
            if column not in evidence_row
        ]
        if missing:
            raise ValueError(
                f'ELM exports lack expected columns: {", ".join(missing)}'
            )
        key = (
            (row['interactorElm'] or '').strip(),
            (row['interactorDomain'] or '').strip(),
            frozenset(filter(None, (row['PMID'] or '').split(','))),
        )
        key = (*key, _taxon(row['taxonomyElm']), _taxon(row['taxonomyDomain']))
        if key != _key(evidence_row):
            raise ValueError(
                f'Row {number} of the ELM interaction table and of the PSI-MI '
                'TAB export are not the same record'
            )
        record = (evidence_row['Interaction_identifier'] or '').strip()
        expected = Counter(
            [
                (
                    key[0],
                    (row['Elm'] or '').strip(),
                    parse_position(row['StartElm']),
                    parse_position(row['StopElm']),
                ),
                (
                    key[1],
                    (row['Domain'] or '').strip(),
                    parse_position(row['StartDomain']),
                    parse_position(row['StopDomain']),
                ),
            ]
        )
        if features.get(record) != expected:
            raise ValueError(
                f'Row {number} of the ELM interaction table is paired with '
                f'{record!r} of the PSI-MI TAB export, but the PSI-MI XML export '
                f'gives that record the features {dict(features.get(record) or {})}, '
                f'not the motif and the domain of the row: {dict(expected)}'
            )
        count += 1
        yield {
            **{column: row[column] for column in COLUMNS},
            **{column: evidence_row[column] for column in EVIDENCE_COLUMNS},
        }
    if not count:
        raise ValueError('ELM interaction table has no rows')


def iter_interactions(
    opener,
    force_refresh: bool = False,
    evidence=None,
    features=None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield the ELM interactions, each with its record number and methods.

    ``opener`` is the interaction table; ``evidence`` and ``features`` are the
    downloads of the PSI-MI TAB export and of the PSI-MI XML export of the
    same records (anything with an ``open`` method that takes
    ``force_refresh``). The rows are those of :func:`join_evidence`.
    """
    if evidence is None or features is None:
        raise ValueError(
            'The PSI-MI TAB and PSI-MI XML exports of the ELM interactions are required'
        )
    yield from join_evidence(
        iter_tsv(opener),
        iter_tsv(evidence.open(force_refresh=force_refresh)),
        read_features(features.open(force_refresh=force_refresh)),
    )


def parse_position(value: object) -> int | None:
    """A residue number of the interaction table; None for ``None`` and for
    0, which ELM writes for a domain whose position it does not give."""
    text = str(value if value is not None else '').strip()
    return int(text) if text.isdigit() and int(text) > 0 else None


def split_methods(cell: object) -> list[str]:
    """The detection methods of a ``Methods`` cell, one PSI-MI TAB item each.

    ELM joins the methods with semicolons, not with the ``|`` of PSI-MI TAB.
    """
    return str(cell or '').split(';')
