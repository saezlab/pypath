"""
Parser of the BioGRID all-organism archive (BioGRID Tab 3.0 format).

The archive ``BIOGRID-ALL-<release>.tab3.zip`` has one text member,
``BIOGRID-ALL-<release>.tab3.txt``, with a header line. One row is one
interaction between two genes reported by one publication with one
experimental system. See
https://wiki.thebiogrid.org/doku.php/biogrid_tab_version_3.0

Besides the columns, every row gets the key ``Ontology Terms`` (``TERMS``):
the six ontology columns decoded into one dictionary per term. A caller that
passes its table of term kinds (``term_kinds``) also gets the terms grouped:
``TERMS_BY_KIND`` and ``QUALIFIERS_BY_KIND`` (accessions per kind of term) and
``PHENOTYPES`` (the phenotype terms of a genetic interaction per vocabulary).
BioGRID's category of a term is not reliable, so the kind of a term is the
caller's table entry for the category and the vocabulary of the accession, and
the table of cell lines (see ``cell_line_accessions``) repairs the one case the
caller names.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Collection, Generator, Mapping
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2 import _mitab
from pypath.inputs_v2.parsers.base import iter_tsv

# The values of the column "Experimental System Type".
SYSTEM_TYPES = ('physical', 'genetic')
# The values of the column "Throughput"; a study of both has both, joined by
# ``|``.
THROUGHPUTS = ('High Throughput', 'Low Throughput')
# Swiss-Prot and TrEMBL accessions, with the isoform number of an accession.
UNIPROT_ACCESSION = (
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|'
    r'[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?'
)
COLUMNS = (
    '#BioGRID Interaction ID',
    'Entrez Gene Interactor A',
    'Entrez Gene Interactor B',
    'BioGRID ID Interactor A',
    'BioGRID ID Interactor B',
    'Systematic Name Interactor A',
    'Systematic Name Interactor B',
    'Official Symbol Interactor A',
    'Official Symbol Interactor B',
    'Synonyms Interactor A',
    'Synonyms Interactor B',
    'Experimental System',
    'Experimental System Type',
    'Publication Source',
    'Organism ID Interactor A',
    'Organism ID Interactor B',
    'Throughput',
    'Score',
    'SWISS-PROT Accessions Interactor A',
    'TREMBL Accessions Interactor A',
    'REFSEQ Accessions Interactor A',
    'SWISS-PROT Accessions Interactor B',
    'TREMBL Accessions Interactor B',
    'REFSEQ Accessions Interactor B',
)

# The columns of the modification and of the ontology terms of an interaction
# (the experimental context). The six ontology columns are lists, one entry per
# term, separated by ``|``; the qualifiers of one term are joined by ``^``.
MODIFICATION = 'Modification'
TERM_LISTS = (
    'Ontology Term IDs',
    'Ontology Term Names',
    'Ontology Term Categories',
    'Ontology Term Qualifier IDs',
    'Ontology Term Qualifier Names',
    'Ontology Term Types',
)
TERM_COLUMNS = (MODIFICATION, *TERM_LISTS)
# The key of the decoded terms in each row, and of their grouping (see
# ``group_terms``).
TERMS = 'Ontology Terms'
TERMS_BY_KIND = 'Ontology Terms By Kind'
QUALIFIERS_BY_KIND = 'Ontology Qualifiers By Kind'
PHENOTYPES = 'Phenotype Terms By Vocabulary'
# BioGRID's category of a term that is a cell line.
CELL_LINE = 'cell line'
# An accession of an ontology term: ``<PREFIX>:<local identifier>``. All 733,990
# rows with terms of release 5.0.262 write the local identifier as digits, in
# every term and every qualifier.
ACCESSION = re.compile(r'([A-Za-z][A-Za-z0-9]*):([0-9]+)')

# The forms of the column "Score" that are a number: plain (``2.5219``),
# negative (``-5.6431``) and with an exponent (``4.75e-06``).
SCORE = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?')
_NUMBER = re.compile(r'\d+')
_UNIPROT = re.compile(UNIPROT_ACCESSION)


def _items(
    cell: str, valid: Callable[[str], bool], none: tuple[str, ...] = ('-',)
) -> bool:
    """Whether a cell is ``-`` or has only valid items, separated by ``|``."""
    return cell in none or all(map(valid, cell.split('|')))


def _checks() -> dict[str, Callable[[str], bool]]:
    """Check of each column of identifiers and categories (not the score)."""
    return {
        'Organism ID Interactor A': lambda c: _NUMBER.fullmatch(c) and int(c) > 0,
        'Organism ID Interactor B': lambda c: _NUMBER.fullmatch(c) and int(c) > 0,
        'BioGRID ID Interactor A': _NUMBER.fullmatch,
        'BioGRID ID Interactor B': _NUMBER.fullmatch,
        'Entrez Gene Interactor A': lambda c: c == '-' or _NUMBER.fullmatch(c),
        'Entrez Gene Interactor B': lambda c: c == '-' or _NUMBER.fullmatch(c),
        'Throughput': lambda c: _items(c, THROUGHPUTS.__contains__) and c != '-',
        **{
            f'{database} Accessions Interactor {side}': lambda c: _items(
                c, _UNIPROT.fullmatch
            )
            for database in ('SWISS-PROT', 'TREMBL')
            for side in 'AB'
        },
        **{
            # The file writes no RefSeq accession as ``-`` and, in 9,105 rows
            # of Interactor A and 9,819 of Interactor B (release 5.0.262),
            # as an empty cell.
            f'REFSEQ Accessions Interactor {side}': lambda c: _items(
                c, _mitab.refseq_namespace, ('-', '')
            )
            for side in 'AB'
        },
    }


def decode_terms(
    row: Mapping[str, str],
    number: str,
    types: Collection[str] | None = None,
) -> list[dict[str, Any]]:
    """
    The ontology terms of a row, one dictionary per term, in file order.

    The six ontology columns must be ``-`` together (no term) or lists of the
    same length. A dictionary has ``accession`` (``BTO:0000007``), ``category``
    (BioGRID's category as written), ``qualifiers`` (the accessions of the
    qualifiers, which may be none) and ``type`` (as written, None for ``-``).
    The term names are not carried. An accession of another form than
    ``<PREFIX>:<digits>``, a qualifier entry that is not ``-`` or accessions
    joined by ``^`` with as many names as accessions, and, if the caller says
    which types it understands (``types``), another type are errors.
    """
    cells = [row[column] for column in TERM_LISTS]
    if all(cell == '-' for cell in cells):
        return []
    lists = [cell.split('|') for cell in cells]
    if len({len(items) for items in lists}) != 1:
        raise ValueError(
            f'Ontology term columns of interaction {number} have lists of '
            f'unequal length: {dict(zip(TERM_LISTS, map(len, lists)))}'
        )
    terms = []
    for accession, _name, category, q_ids, q_names, type_ in zip(*lists):
        match = ACCESSION.fullmatch(accession)
        if not match:
            raise ValueError(
                f'Unknown ontology term accession {accession!r} '
                f'in interaction {number}'
            )
        if types is not None and type_ != '-' and type_ not in types:
            raise ValueError(
                f'Unknown ontology term type {type_!r} in interaction {number}'
            )
        qualifiers = [] if q_ids == '-' else q_ids.split('^')
        if len(qualifiers) != (0 if q_names == '-' else len(q_names.split('^'))):
            raise ValueError(
                f'Qualifier names {q_names!r} of ontology term {accession} '
                f'do not match qualifier identifiers {q_ids!r} '
                f'in interaction {number}'
            )
        for qualifier in qualifiers:
            if not ACCESSION.fullmatch(qualifier):
                raise ValueError(
                    f'Unknown ontology term qualifier {qualifier!r} '
                    f'in interaction {number}'
                )
        terms.append(
            {
                'accession': accession,
                'category': category,
                'qualifiers': qualifiers,
                'type': None if type_ == '-' else type_,
            }
        )
    return terms


def cell_line_accessions(rows) -> frozenset[str]:
    """
    The accessions that BioGRID's category column calls a cell line in at
    least one row of the table.

    BioGRID labels the same cell line ``cell line`` in some rows and
    ``phenotype`` in others (HEK-293 cell, HeLa cell), so a label is only
    taken from the file as a whole. ``rows`` are the rows of the table; the
    terms are read by the position of the entries in the two columns, and a
    row whose lists do not match raises later, when its terms are decoded.
    """
    found = set()
    for row in rows:
        accessions = row.get(TERM_LISTS[0])
        categories = row.get(TERM_LISTS[2])
        if accessions is None or categories is None or accessions == '-':
            continue
        found.update(
            accession
            for accession, category in zip(
                accessions.split('|'), categories.split('|')
            )
            if category == CELL_LINE
        )
    return frozenset(found)


def group_terms(
    terms: list[dict[str, Any]],
    number: str,
    row_type: str,
    term_kinds: Mapping[tuple[str, str], str],
    cell_lines: Collection[str] = (),
    cell_line_repairs: Collection[tuple[str, str]] = (),
    phenotype_kind: str | None = None,
    phenotype_vocabularies: Collection[str] = (),
) -> tuple[dict[str, list[str]], dict[str, list[str]], dict[str, list[str]]]:
    """
    The decoded terms of a row grouped by kind.

    ``term_kinds`` maps (BioGRID category, prefix of the accession) to the kind
    of the term; a combination it lacks is an error. A term whose combination
    is in ``cell_line_repairs`` and whose accession is in ``cell_lines`` has
    the kind of the same vocabulary's ``cell line`` category instead. Each
    qualifier has the kind of its term. The result is the accessions of the
    terms by kind, the accessions of the qualifiers by kind and, for a genetic
    interaction, the accessions of the terms of kind ``phenotype_kind`` by
    prefix (these are not among the terms by kind): the prefix must be in
    ``phenotype_vocabularies``. Every list is free of repeated accessions and
    in file order.
    """
    by_kind: dict[str, list[str]] = {}
    qualifiers: dict[str, list[str]] = {}
    phenotypes: dict[str, list[str]] = {}
    for term in terms:
        accession = term['accession']
        prefix = accession.partition(':')[0]
        combination = (term['category'], prefix)
        kind = term_kinds.get(combination)
        if kind is None:
            raise ValueError(
                f'Unknown BioGRID ontology term category {term["category"]!r} '
                f'for the vocabulary {prefix!r} (term {accession}) '
                f'in interaction {number}'
            )
        if combination in cell_line_repairs and accession in cell_lines:
            kind = term_kinds[(CELL_LINE, prefix)]
        if row_type == 'genetic' and kind == phenotype_kind:
            if prefix not in phenotype_vocabularies:
                raise ValueError(
                    f'Unknown BioGRID phenotype vocabulary {prefix!r} '
                    f'(term {accession}) in genetic interaction {number}'
                )
            phenotypes.setdefault(prefix, []).append(accession)
        else:
            by_kind.setdefault(kind, []).append(accession)
        if term['qualifiers']:
            qualifiers.setdefault(kind, []).extend(term['qualifiers'])
    return tuple(  # type: ignore[return-value]
        {key: list(dict.fromkeys(values)) for key, values in group.items()}
        for group in (by_kind, qualifiers, phenotypes)
    )


def iter_interactions(
    opener,
    system_type: str = 'physical',
    release: str | None = None,
    systems: Mapping[str, Collection[str]] | None = None,
    term_types: Collection[str] | None = None,
    modifications: Collection[str] | None = None,
    term_kinds: Mapping[tuple[str, str], str] | None = None,
    cell_line_repairs: Collection[tuple[str, str]] = (),
    phenotype_kind: str | None = None,
    phenotype_vocabularies: Collection[str] = (),
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Rows of one experimental system type from the all-organism archive.

    ``release`` is the release the archive must be of (its member is
    ``BIOGRID-ALL-<release>.tab3.txt``); with none, any release. ``systems``
    maps each system type to the names of the experimental systems the caller
    understands. A missing or additional member, a missing column, a row with
    another number of fields and a value the caller does not understand are
    errors: they mean the download or its layout changed. They are checked in
    every row of the file, not only in the rows of the system type asked for.
    The values checked are the experimental system and its type, the
    throughput, the organisms, the BioGRID, Entrez Gene, UniProt and RefSeq
    identifiers, the ontology term columns (see ``decode_terms``); the score
    is ``-`` or a number. ``term_types`` and ``modifications`` are the term
    types and modifications the caller understands; with none, any.

    ``term_kinds``, ``cell_line_repairs``, ``phenotype_kind`` and
    ``phenotype_vocabularies`` are the caller's table of the kinds of term
    (see ``group_terms``); a combination of category and vocabulary that
    ``term_kinds`` lacks, in any row, is an error. With ``cell_line_repairs``
    the archive member is read twice (it must be seekable), first for the
    accessions that some row calls a cell line (see ``cell_line_accessions``).

    Each row is a dictionary keyed by the 37 header names, all values text as
    written; ``-`` stands for no value. It has four more keys: ``TERMS``
    (``Ontology Terms``), the list of terms that ``decode_terms`` returns, and
    the three groupings of ``group_terms`` (``TERMS_BY_KIND``,
    ``QUALIFIERS_BY_KIND``, ``PHENOTYPES``), which are empty dictionaries
    without ``term_kinds``.
    """
    if system_type not in SYSTEM_TYPES:
        raise ValueError(
            f'BioGRID interactions are {SYSTEM_TYPES}, not {system_type!r}'
        )
    if not opener or not opener.result:
        return
    pattern = re.compile(
        r'BIOGRID-ALL-%s\.tab3\.txt'
        % (re.escape(release) if release else r'[^/]+')
    )
    members = [name for name in opener.result if pattern.fullmatch(name)]
    if len(members) != 1:
        raise ValueError(
            f'BioGRID archive has {len(members)} members named {pattern.pattern}'
        )
    member = members[0]
    handle = opener.result[member]
    source = SimpleNamespace(result=handle)
    cell_lines: frozenset[str] = frozenset()
    if term_kinds is not None and cell_line_repairs:
        if not (hasattr(handle, 'seekable') and handle.seekable()):
            raise ValueError(f'{member} cannot be read a second time')
        cell_lines = cell_line_accessions(
            _mitab.checked_rows(
                iter_tsv(source, quotechar=None, **kwargs), COLUMNS, member
            )
        )
        handle.seek(0)
    rows = _mitab.checked_rows(
        iter_tsv(source, quotechar=None, **kwargs), COLUMNS, member
    )
    checks = _checks()
    valid: dict[str, set[str]] = {column: set() for column in checks}
    for row in rows:
        number = row['#BioGRID Interaction ID']
        row_type = row['Experimental System Type']
        if row_type not in SYSTEM_TYPES:
            raise ValueError(
                f'{member}: unknown experimental system type {row_type!r} '
                f'in interaction {number}'
            )
        system = row['Experimental System']
        if systems is not None and system not in systems.get(row_type, ()):
            raise ValueError(
                f'{member}: unknown {row_type} experimental system '
                f'{system!r} in interaction {number}'
            )
        score = row['Score']
        if score != '-' and not SCORE.fullmatch(score):
            raise ValueError(
                f'{member}: unknown Score {score!r} in interaction {number}'
            )
        # The columns of the experimental context, checked after the layout
        # above so that a header without them is reported by name.
        missing = [c for c in TERM_COLUMNS if row.get(c) is None]
        if missing:
            raise ValueError(
                f'{member}: interaction {number} lacks the columns {missing}'
            )
        modification = row[MODIFICATION]
        if (
            modifications is not None
            and modification != '-'
            and modification not in modifications
        ):
            raise ValueError(
                f'{member}: unknown {MODIFICATION} {modification!r} '
                f'in interaction {number}'
            )
        row[TERMS] = decode_terms(row, number, term_types)
        if term_kinds is None:
            row[TERMS_BY_KIND] = row[QUALIFIERS_BY_KIND] = row[PHENOTYPES] = {}
        else:
            (
                row[TERMS_BY_KIND],
                row[QUALIFIERS_BY_KIND],
                row[PHENOTYPES],
            ) = group_terms(
                row[TERMS], number, row_type, term_kinds, cell_lines,
                cell_line_repairs, phenotype_kind, phenotype_vocabularies,
            )
        for column, check in checks.items():
            cell = row[column]
            if cell not in valid[column]:
                if not check(cell):
                    raise ValueError(
                        f'{member}: unknown {column} {cell!r} in '
                        f'interaction {number}'
                    )
                valid[column].add(cell)
        if row_type == system_type:
            yield row
