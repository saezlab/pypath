"""Read the protein/RNA and experiment tables of ExoCarta and Vesiclepedia.

Both databases publish two tab-separated tables with the same layout (the
legacy reader served both with one function).

``*_PROTEIN_MRNA_DETAILS`` lists every molecule found in an experiment, one
row per molecule and experiment::

    CONTENT ID  CONTENT TYPE  ENTREZ GENE ID  GENE SYMBOL  SPECIES
    EXPERIMENT ID  METHODS

``CONTENT TYPE`` is ``protein``, ``mRNA`` or, in ExoCarta, the type of another
RNA, spelled with varying case (``Protein``, ``mrna``, ``protein `` with a
trailing space). ``METHODS`` are the methods that detected the molecule,
separated by ``|``. ``SPECIES`` is the species of the molecule.

``*_EXPERIMENT_DETAILS`` lists the studies::

    EXPERIMENT ID  PUBMED ID  SPECIES  EXPERIMENT DESCRIPTION  SAMPLE
    SAMPLE SOURCE  SAMPLE NAME  IDENTIFICATIONS  METHODS  YEAR
    ISOLATION METHOD  [VESICLE TYPE]

``PUBMED ID`` is 0 for a study with none. ``SAMPLE`` is the kind of cell or
fluid the vesicles are from (``Urine``, ``Melanoma cells``), ``SAMPLE NAME``
the cell line or tissue, ``SAMPLE SOURCE`` ``cl``, ``ts`` or ``fl`` (or both
of the first two, separated by ``|``) and 0 where it has none. Vesiclepedia's
table has the ``VESICLE TYPE`` column, with one or several types separated by
``/``; ExoCarta's has none, and every ExoCarta study is of exosomes.

Whitespace around values is removed: the tables write species, vesicle types
and content types with trailing spaces.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.base import Download
from pypath.inputs_v2.parsers.base import iter_tsv

EXOCARTA = 'exocarta'
VESICLEPEDIA = 'vesiclepedia'
# ExoCarta is the database of exosomes; its experiment table has no vesicle type.
EXOCARTA_VESICLE_TYPES = ('Exosomes',)

DETECTION_COLUMNS = (
    'CONTENT TYPE', 'ENTREZ GENE ID', 'GENE SYMBOL', 'SPECIES', 'EXPERIMENT ID', 'METHODS',
)
EXPERIMENT_COLUMNS = (
    'EXPERIMENT ID', 'PUBMED ID', 'SPECIES', 'SAMPLE', 'SAMPLE SOURCE', 'SAMPLE NAME',
    'ISOLATION METHOD',
)

# The content types of the tables, by their lower-case spelling, with the
# spelling this module uses. A type that is not here raises.
CONTENT_TYPES = {
    'protein': 'protein',
    'mrna': 'mRNA',
    'snrna': 'snRNA',
    'snorna': 'snoRNA',
    'trna': 'tRNA',
    'rrna': 'rRNA',
    'scarna': 'scaRNA',
    'lncrna': 'lncRNA',
    'lincrna': 'lincRNA',
    'yrna': 'yRNA',
    'vtrna': 'vtRNA',
    'ncrna': 'ncRNA',
    'pseudogene': 'pseudogene',
}
# The values of SAMPLE SOURCE, as the tables write them (cell line, tissue and
# fluid, by the look of the samples); any other value raises.
SAMPLE_SOURCES = frozenset({'cl', 'ts', 'fl'})
# What the tables write where a field has no value.
_NO_VALUE = frozenset({'', '0', '-'})


def _list(value: str | None, separator: str) -> list[str]:
    """The parts of a cell, stripped and without repeats or empty parts."""
    parts: list[str] = []
    for part in (value or '').split(separator):
        part = part.strip()
        if part and part not in parts:
            parts.append(part)
    return parts


def vesicle_types(experiment: dict[str, str], database: str) -> list[str]:
    """The vesicle types of an experiment row."""
    if database == EXOCARTA:
        return list(EXOCARTA_VESICLE_TYPES)
    if database != VESICLEPEDIA:
        raise ValueError(f'Unknown database: {database!r}')
    return _list(experiment['VESICLE TYPE'], '/')


def _sample_sources(row: dict[str, str]) -> list[str]:
    sources = [s for s in _list(row['SAMPLE SOURCE'], '|') if s not in _NO_VALUE]
    unknown = [s for s in sources if s not in SAMPLE_SOURCES]
    if unknown:
        raise ValueError(
            f'Unknown sample source {unknown} in experiment {row["EXPERIMENT ID"].strip()}'
        )
    return sources


def _experiment(row: dict[str, str], database: str) -> dict[str, Any]:
    pmid = row['PUBMED ID'].strip()
    return {
        'experiment': row['EXPERIMENT ID'].strip(),
        'pmid': None if pmid in _NO_VALUE else pmid,
        'sample': row['SAMPLE'].strip() or None,
        'sample_name': row['SAMPLE NAME'].strip() or None,
        'sample_source': _sample_sources(row),
        'isolation_methods': [s for s in _list(row['ISOLATION METHOD'], '|') if s not in _NO_VALUE],
        'vesicle_types': vesicle_types(row, database),
    }


def _experiments(
    experiment_download: Download, database: str, force_refresh: bool
) -> dict[str, dict[str, Any]]:
    required = EXPERIMENT_COLUMNS + (('VESICLE TYPE',) if database == VESICLEPEDIA else ())
    opener = experiment_download.open(force_refresh=force_refresh)
    experiments = {}
    for row in iter_tsv(opener, quotechar=None, required=required):
        experiment = _experiment(row, database)
        if experiment['experiment'] in experiments:
            raise ValueError(f'Experiment {experiment["experiment"]} is listed twice')
        experiments[experiment['experiment']] = experiment
    return experiments


def iter_detections(
    opener,
    *,
    experiment_download: Download,
    database: str,
    force_refresh: bool = False,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """One row for each vesicle type of each molecule found in an experiment.

    Fields: ``entrez``, ``symbol``, ``species`` (stripped), ``content_type``
    (see ``CONTENT_TYPES``), ``methods`` (list), ``vesicle_type`` and the
    fields of the experiment: ``experiment``, ``pmid`` (None for 0),
    ``sample``, ``sample_name``, ``sample_source`` (list) and
    ``isolation_methods`` (list). A row of an experiment that is not in the
    experiment table, or of an unknown content type, raises.
    """
    experiments = _experiments(experiment_download, database, force_refresh)
    for row in iter_tsv(opener, quotechar=None, required=DETECTION_COLUMNS):
        experiment = experiments.get(row['EXPERIMENT ID'].strip())
        if experiment is None:
            raise ValueError(f'Unknown experiment {row["EXPERIMENT ID"]!r}')
        content_type = CONTENT_TYPES.get(row['CONTENT TYPE'].strip().lower())
        if content_type is None:
            raise ValueError(f'Unknown content type {row["CONTENT TYPE"]!r}')
        fields = {key: value for key, value in experiment.items() if key != 'vesicle_types'}
        for vesicle_type in experiment['vesicle_types']:
            yield {
                'entrez': row['ENTREZ GENE ID'].strip(),
                'symbol': row['GENE SYMBOL'].strip(),
                'species': row['SPECIES'].strip(),
                'content_type': content_type,
                'methods': _list(row['METHODS'], '|'),
                'vesicle_type': vesicle_type,
                **fields,
            }


def iter_vesicle_types(opener, *, database: str, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The vesicle types of the experiment table, each once: ``name``."""
    required = EXPERIMENT_COLUMNS + (('VESICLE TYPE',) if database == VESICLEPEDIA else ())
    seen: set[str] = set()
    for row in iter_tsv(opener, quotechar=None, required=required):
        for name in vesicle_types(row, database):
            if name not in seen:
                seen.add(name)
                yield {'name': name}
