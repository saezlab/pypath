"""Read the significant preys of the Human Cell Map (HEK293 BioID proximity map).

``preys-v1.txt`` (the version 1 of ``preys-latest.txt``) lists the proteins
detected as significant interactors of a BioID bait, with the localization the
map predicts for them by two methods, one row per prey::

    symbol  NMF rank  MMF localization  SAFE domain  SAFE localization
    Entrez  BioGRID  HGNC  RefSeq  UniProt

``MMF localization`` is the column of the non-negative matrix factorization
(NMF); the file's header spells the abbreviation with an M. ``NMF rank`` is
the rank of the factorization the prey is assigned to. ``SAFE domain`` is the
domain of the spatial analysis of functional enrichment (SAFE) network, ``-``
for a prey outside every domain. Both localization columns hold one or more
labels of the map's compartments separated by commas, ``no prediction`` or
``-`` for none, and, for NMF, ``miscellaneous`` for the rank of the preys
that the factorization did not place in a compartment. ``UniProt`` is empty
for 50 of the preys.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = (
    'symbol',
    'NMF rank',
    'MMF localization',
    'SAFE domain',
    'SAFE localization',
    'Entrez',
    'HGNC',
    'UniProt',
)
NMF = 'NMF'
SAFE = 'SAFE'
# What the localization columns write where the map has no localization to
# give, and the label of the NMF rank that has no compartment.
NO_LABELS = frozenset({'-', 'no prediction'})
MISCELLANEOUS = 'miscellaneous'

# The localization column, and the column of the rank or domain, of each method.
_METHODS = ((NMF, 'MMF localization', 'NMF rank'), (SAFE, 'SAFE localization', 'SAFE domain'))


def _labels(value: str) -> list[str]:
    return [label.strip() for label in value.split(',') if label.strip()]


def iter_localizations(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One row for each label of each method of each prey.

    Fields: ``uniprot`` (None where the file has none), ``symbol``, ``entrez``,
    ``hgnc``, ``method`` (NMF or SAFE), ``label`` and ``group`` (the NMF rank or
    the SAFE domain). A prey with no label of a method has no row for it, and
    neither has the NMF label ``miscellaneous``, which is no compartment.
    """
    for row in iter_tsv(opener, quotechar=None, required=COLUMNS):
        for method, column, group_column in _METHODS:
            for label in _labels(row[column]):
                if label in NO_LABELS or label == MISCELLANEOUS:
                    continue
                yield {
                    'uniprot': row['UniProt'] or None,
                    'symbol': row['symbol'],
                    'entrez': row['Entrez'] or None,
                    'hgnc': row['HGNC'] or None,
                    'method': method,
                    'label': label,
                    'group': row[group_column],
                }


def iter_terms(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The labels of the map, each once: ``label`` and the ``methods`` that use it."""
    methods: dict[str, set[str]] = {}
    for row in iter_localizations(opener):
        methods.setdefault(row['label'], set()).add(row['method'])
    for label, used in methods.items():
        yield {'label': label, 'methods': sorted(used)}
