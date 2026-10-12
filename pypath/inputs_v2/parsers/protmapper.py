"""Join the two tables of the ProtMapper corpus of phosphosite annotations.

The supplementary archive holds ``export.csv``, one row per site annotation
(a controller, a target protein and a residue), and ``evidences.csv``, one row
per piece of evidence for an annotation. The two share the ``ID`` column.

An evidence row has a site of its own: ``MAPPED_*`` (protein, residue,
position after ProtMapper's mapping, filled when ``DESCRIPTION`` is not
"VALID") or else ``ORIG_*``. The annotation's ``TARGET_*`` is "the final
(valid or mapped) site" (README.txt of the archive). In version 4 of the corpus 5,301
of the 115,667 evidence rows have an own site that differs from the
annotation's, each a repeat of an evidence row of another annotation of the
same controller; the join leaves those rows out (``join_evidences``).
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Generator, Iterable
import logging
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2.parsers.base import iter_csv

_log = logging.getLogger(__name__)

ANNOTATIONS_FILE = 'data/export.csv'
EVIDENCES_FILE = 'data/evidences.csv'

ANNOTATION_COLUMNS = (
    'ID',
    'CTRL_NS',
    'CTRL_ID',
    'CTRL_GENE_NAME',
    'CTRL_IS_KINASE',
    'TARGET_UP_ID',
    'TARGET_GENE_NAME',
    'TARGET_RES',
    'TARGET_POS',
    'SOURCES',
    'BELIEF',
)
EVIDENCE_COLUMNS = (
    'ID',
    'SOURCE',
    'PMID',
    'DBID',
    'TEXT',
    'DESCRIPTION',
    'ORIG_UP_ID',
    'ORIG_RES',
    'ORIG_POS',
    'MAPPED_UP_ID',
    'MAPPED_RES',
    'MAPPED_POS',
)
ORIG_COLUMNS = ('ORIG_UP_ID', 'ORIG_RES', 'ORIG_POS')
MAPPED_COLUMNS = ('MAPPED_UP_ID', 'MAPPED_RES', 'MAPPED_POS')


def _checked(
    rows: Iterable[dict[str, Any]], columns: tuple[str, ...], table: str
) -> Generator[dict[str, Any], None, None]:
    """Rows of a table; a changed header is an error, not an empty dataset."""
    count = 0
    for row in rows:
        missing = [column for column in columns if column not in row]
        if missing:
            raise ValueError(f'ProtMapper table {table} lacks columns {missing}')
        count += 1
        yield row
    if not count:
        raise ValueError(f'ProtMapper table {table} has no rows')


def evidence_site(row: dict[str, Any]) -> tuple[str, str, str]:
    """The protein, residue and position an evidence row states for itself.

    ``MAPPED_*`` is the site after ProtMapper's mapping; it is filled for the
    rows whose ``DESCRIPTION`` is not "VALID". Otherwise the site is the
    original one (``ORIG_*``), which was valid as it stood. The three
    ``MAPPED_*`` columns are filled together or all empty.
    """
    mapped = tuple(row[column] for column in MAPPED_COLUMNS)
    if any(mapped) and not all(mapped):
        raise ValueError(f'ProtMapper evidence has a partly empty mapped site: {row!r}')
    site = mapped if all(mapped) else tuple(row[column] for column in ORIG_COLUMNS)
    if not all(site):
        raise ValueError(f'ProtMapper evidence has no site: {row!r}')
    return site


def site_mismatch(annotation: dict[str, Any], evidence: dict[str, Any]) -> str | None:
    """Why an evidence row does not support the site of its annotation, or None.

    The evidence row's own site (``evidence_site``) is compared with the
    annotation's ``TARGET_UP_ID``, ``TARGET_RES`` and ``TARGET_POS``, all as
    written. The reason, the first that applies:

    - ``protein``: another UniProt entry;
    - ``isoform``: the same entry, another accession (canonical sequence or
      isoform), so the positions count another sequence;
    - ``position``: the same protein, another position;
    - ``residue``: the same protein and position, another residue letter.
    """
    protein, residue, position = evidence_site(evidence)
    target = annotation['TARGET_UP_ID']
    if protein != target:
        same_entry = protein.split('-')[0] == target.split('-')[0]
        return 'isoform' if same_entry else 'protein'
    if position != annotation['TARGET_POS']:
        return 'position'
    if residue != annotation['TARGET_RES']:
        return 'residue'
    return None


def join_evidences(
    annotations: Iterable[dict[str, Any]],
    evidences: Iterable[dict[str, Any]],
    stats: Counter | None = None,
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per evidence row, with the columns of its annotation.

    The annotation's columns (``CTRL_*``, ``TARGET_*``, ``SOURCES``,
    ``BELIEF``) are repeated on each of its evidence rows; the evidence's own
    columns (``SOURCE``, ``PMID``, ``DBID``, ``TEXT``, ``DESCRIPTION``,
    ``ORIG_*``, ``MAPPED_*``) are as written. The rows come in the order of
    the evidence table.

    An evidence row whose own site differs from the site of its annotation
    (``site_mismatch``) is not yielded: its publication and source belong to
    another site. In version 4 of the corpus each of these rows is repeated, with the
    same source, PMID, database entry and text, under the annotation of the
    same controller and the row's own site. An annotation left with no
    evidence row is yielded once, at the end, with the evidence columns set
    to None. An evidence row that names an unknown annotation is an error.

    ``stats``, if given, is filled with the counts: ``evidence_rows``,
    ``evidence_attached``, ``evidence_dropped``, ``evidence_dropped_<reason>``,
    ``annotations``, ``annotations_with_dropped_evidence`` and
    ``annotations_without_evidence``, which is split in
    ``annotations_without_evidence_rows`` (the table has none for them) and
    ``annotations_with_only_dropped_evidence``.
    """
    stats = Counter() if stats is None else stats
    by_id: dict[str, dict[str, Any]] = {}
    for row in _checked(annotations, ANNOTATION_COLUMNS, 'export.csv'):
        if row['ID'] in by_id:
            raise ValueError(f'ProtMapper annotation ID is not unique: {row["ID"]!r}')
        by_id[row['ID']] = row
    stats['annotations'] = len(by_id)
    attached: set[str] = set()
    dropped: set[str] = set()
    for row in _checked(evidences, EVIDENCE_COLUMNS, 'evidences.csv'):
        annotation = by_id.get(row['ID'])
        if annotation is None:
            raise ValueError(
                f'ProtMapper evidence names an unknown annotation: {row["ID"]!r}'
            )
        stats['evidence_rows'] += 1
        reason = site_mismatch(annotation, row)
        if reason:
            stats['evidence_dropped'] += 1
            stats[f'evidence_dropped_{reason}'] += 1
            dropped.add(row['ID'])
            continue
        stats['evidence_attached'] += 1
        attached.add(row['ID'])
        yield {**annotation, **row}
    without_evidence = dict.fromkeys(EVIDENCE_COLUMNS[1:])
    stats['annotations_with_dropped_evidence'] = len(dropped)
    for identifier, annotation in by_id.items():
        if identifier not in attached:
            stats['annotations_without_evidence'] += 1
            if identifier in dropped:
                stats['annotations_with_only_dropped_evidence'] += 1
            else:
                stats['annotations_without_evidence_rows'] += 1
            yield {**annotation, **without_evidence}


def iter_evidences(
    opener, stats: Counter | None = None, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Read both tables from the opened archive and join them (``join_evidences``).

    The counts of ``join_evidences`` go to ``stats`` and, when the join is
    complete, to the log.
    """
    if not opener or not opener.result:
        return
    members = opener.result
    missing = [
        name
        for name in (ANNOTATIONS_FILE, EVIDENCES_FILE)
        if not isinstance(members, dict) or name not in members
    ]
    if missing:
        raise ValueError(f'ProtMapper archive lacks {missing}')
    stats = Counter() if stats is None else stats
    yield from join_evidences(
        iter_csv(SimpleNamespace(result=members[ANNOTATIONS_FILE])),
        iter_csv(SimpleNamespace(result=members[EVIDENCES_FILE])),
        stats,
    )
    _log.info('ProtMapper evidence join: %s', dict(sorted(stats.items())))
