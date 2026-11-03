"""Columns of the iPTMnet download tables and the join of scores to PTM rows.

``ptm.txt`` and ``score.txt`` are tab separated and have no header row; the
column names are those of the download's ``readme.txt``. The tables are read
with the shared delimited reader.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Generator, Iterable, Sequence
import logging
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

_log = logging.getLogger(__name__)

PTM_COLUMNS = (
    'ptm_type',
    'source',
    'substrate_UniProtAC',
    'substrate_genename',
    'organism',
    'site',
    'enzyme_UniProtAC',
    'enzyme_genename',
    'note',
    'pmid',
)
SCORE_COLUMNS = (
    'substrate_UniProtAC',
    'site',
    'enzyme_UniProtAC',
    'ptm_type',
    'score',
)
# A UniProt accession, with the isoform number where there is one.
UNIPROT = (
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|'
    r'[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?'
)
UNIPROT_REGEX = re.compile(rf'^{UNIPROT}$')
# The Protein Ontology (PRO) names the protein of a UniProt entry, and an
# isoform, by the accession: ``PR:Q92793``, ``PR:Q92793-1``. Each such term has
# the same accession as its UniProtKB cross-reference (checked for the 63
# accessions of ``ptm.txt`` in PRO release 73.1, 2026-09-15).
PRO_PROTEIN_REGEX = re.compile(rf'^PR:({UNIPROT})$')
# The other PRO terms, ``PR:000027108``, are protein classes, modified forms,
# families or complexes; the table does not say which.
PRO_TERM_REGEX = re.compile(r'^PR:\d{9}$')
# PhosphoGRID names a yeast kinase by the systematic name of its gene.
YEAST_LOCUS_REGEX = re.compile(r'^Y[A-P][LR]\d{3}[WC](?:-[A-Z])?$')


def enzyme_kind(value: str) -> str | None:
    """The identifier system of an ``enzyme_UniProtAC`` cell, which holds four.

    ``uniprot``, ``pro_protein`` (PRO term of a UniProt entry),
    ``pro_term`` (any other PRO term) or ``yeast_locus``; None for a value of
    no known form, such as a Gene Ontology term.
    """
    for kind, regex in (
        ('uniprot', UNIPROT_REGEX),
        ('pro_protein', PRO_PROTEIN_REGEX),
        ('pro_term', PRO_TERM_REGEX),
        ('yeast_locus', YEAST_LOCUS_REGEX),
    ):
        if regex.match(value):
            return kind
    return None


def _score_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    # The two tables write the PTM type in different case.
    return (
        row['substrate_UniProtAC'],
        row['site'],
        row['enzyme_UniProtAC'],
        row['ptm_type'].upper(),
    )


def join_scores(
    ptm_rows: Iterable[dict[str, Any]],
    score_rows: Iterable[dict[str, Any]],
    stats: Counter | None = None,
) -> Generator[dict[str, Any], None, None]:
    """Yield the PTM rows that name an enzyme, each with its ``score``.

    A row names an enzyme if ``enzyme_UniProtAC`` or ``enzyme_genename`` is
    filled; the other rows are sites without an enzyme and are skipped.
    ``score`` is the score of the row with the same substrate accession,
    site, enzyme accession and PTM type in the score table, all as written,
    or None if the score table has no such row. Nothing is matched on a
    shortened or translated accession. A row that names its enzyme only by
    ``enzyme_genename`` (``enzyme_UniProtAC`` empty) has no score: the score
    table does not name the enzyme of its rows without an accession.

    A row with an enzyme is malformed, and skipped, if its
    ``enzyme_UniProtAC`` is of no known form (``enzyme_kind``), or its
    substrate accession or site is empty. ``stats``, if given, is filled with
    the counts: ``rows``, ``sites_without_enzyme``, ``enzyme_rows``,
    ``enzyme_label_only`` (the enzyme only by label, yielded),
    ``malformed_enzyme_accession``, ``malformed_substrate_accession`` and
    ``malformed_site``.
    """
    stats = Counter() if stats is None else stats
    scores: dict[tuple[str, str, str, str], str] = {}
    for row in score_rows:
        if not row['enzyme_UniProtAC']:
            continue
        key = _score_key(row)
        if scores.setdefault(key, row['score']) != row['score']:
            raise ValueError(f'iPTMnet score table has two scores for {key}')
    for row in ptm_rows:
        stats['rows'] += 1
        accession = row['enzyme_UniProtAC']
        if not (accession or row['enzyme_genename']):
            stats['sites_without_enzyme'] += 1
            continue
        stats['enzyme_rows'] += 1
        if accession and enzyme_kind(accession) is None:
            stats['malformed_enzyme_accession'] += 1
            _log.warning('iPTMnet row with enzyme %r of no known form skipped', accession)
            continue
        if not row['substrate_UniProtAC']:
            stats['malformed_substrate_accession'] += 1
            _log.warning('iPTMnet row without a substrate accession skipped: %r', row)
            continue
        if not row['site']:
            stats['malformed_site'] += 1
            _log.warning('iPTMnet row without a site skipped: %r', row)
            continue
        if not accession:
            stats['enzyme_label_only'] += 1
        yield {**row, 'score': scores.get(_score_key(row))}
    if not stats['rows']:
        raise ValueError('iPTMnet PTM table has no rows')


def iter_table(opener, columns: Sequence[str]) -> Generator[dict[str, str], None, None]:
    """Rows of an iPTMnet table, which has no header row, with cells stripped.

    Quotes are not interpreted, and a line with another number of fields
    than ``columns`` is an error.
    """
    for row in iter_tsv(opener, fieldnames=columns, strict=True, quotechar=None):
        yield {column: value.strip() for column, value in row.items()}


def iter_enzyme_sites(
    opener,
    scores_download,
    force_refresh: bool = False,
    stats: Counter | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of ``ptm.txt`` that name an enzyme, with their ``score``.

    ``opener`` is the opened ``ptm.txt``; ``scores_download`` is the
    ``Download`` of ``score.txt``, opened here with the same ``force_refresh``
    (``join_scores`` describes the fields).
    """
    stats = Counter() if stats is None else stats
    scores = scores_download.open(force_refresh=force_refresh)
    yield from join_scores(
        iter_table(opener, PTM_COLUMNS),
        iter_table(scores, SCORE_COLUMNS),
        stats,
    )
    _log.info('iPTMnet ptm.txt: %s', dict(sorted(stats.items())))
