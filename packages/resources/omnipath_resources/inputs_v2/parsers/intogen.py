"""
Readers for the IntOGen release archives.

An IntOGen release is two zip archives, each with a folder named after the
release date: ``IntOGen-Drivers`` (``Compendium_Cancer_Genes.tsv``, the driver
genes per cohort, and ``Unfiltered_drivers.tsv``, the candidates before the
filters) and ``IntOGen-Cohorts`` (``cohorts.tsv``). Both tables are tab
separated, with a header row. The shared TSV reader reads them;
the functions here select the table in the archive, check the header and give
the cells their shape.
"""

from __future__ import annotations

from collections.abc import Generator
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['iter_cancer_types', 'iter_drivers']

DRIVER_COLUMNS = (
    'SYMBOL',
    'TRANSCRIPT',
    'COHORT',
    'CANCER_TYPE',
    'METHODS',
    'MUTATIONS',
    'SAMPLES',
    'QVALUE_COMBINATION',
    'ROLE',
    'IS_DRIVER',
)
COHORT_COLUMNS = ('COHORT', 'CANCER', 'CANCER_NAME')
# What the file writes for a missing value.
EMPTY = ('', 'NA', 'nan', 'NULL')


def _member(opener, name: str):
    """The table ``name`` of the archive, as an opener of its own."""
    if not opener or not isinstance(opener.result, dict):
        raise ValueError('IntOGen archive not available')
    for key, handle in opener.result.items():
        if key.rsplit('/', 1)[-1] == name and handle is not None:
            return SimpleNamespace(result=handle)
    raise ValueError(f'IntOGen archive has no {name}: {sorted(opener.result)}')


def _cell(value: str | None) -> str | None:
    value = None if value is None else value.strip()
    return None if value in EMPTY else value


def iter_drivers(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one record per driver gene of a cohort, from ``Compendium_Cancer_Genes.tsv``.

    The compendium holds the genes that pass the filters of IntOGen (its column
    ``IS_DRIVER`` is ``True`` on each line; another value raises). Output
    fields, the columns of the file with empty cells as None:

    - ``symbol``, ``transcript``: the gene's HUGO symbol and the Ensembl
      transcript the analysis used (``SYMBOL``, ``TRANSCRIPT``).
    - ``cohort``, ``cancer_type``: the cohort and the tumor type of the cohort.
    - ``methods``: the list of the methods that find the gene a driver, as
      written (``oncodriveclustl``, ``combination``).
    - ``mutations``, ``samples``: the number of mutations of the gene in the
      cohort and of samples that have them, as written.
    - ``qvalue``: the significance of the combined output.
    - ``role``: the mode of action of the gene, ``Act``, ``LoF`` or
      ``ambiguous``.
    """
    for row in iter_tsv(_member(opener, 'Compendium_Cancer_Genes.tsv'), required=DRIVER_COLUMNS):
        if row['IS_DRIVER'] != 'True':
            raise ValueError(f'IntOGen compendium: IS_DRIVER is {row["IS_DRIVER"]!r} for {row["SYMBOL"]}')
        methods = _cell(row['METHODS'])
        yield {
            'symbol': _cell(row['SYMBOL']),
            'transcript': _cell(row['TRANSCRIPT']),
            'cohort': _cell(row['COHORT']),
            'cancer_type': _cell(row['CANCER_TYPE']),
            'methods': methods.split(',') if methods else [],
            'mutations': _cell(row['MUTATIONS']),
            'samples': _cell(row['SAMPLES']),
            'qvalue': _cell(row['QVALUE_COMBINATION']),
            'role': _cell(row['ROLE']),
        }


def iter_cancer_types(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield each cancer type of ``cohorts.tsv`` once, in order of first appearance.

    Output fields: ``cancer_type``, the acronym of the type (``CANCER``), and
    ``name``, its long name (``CANCER_NAME``). A type with two names raises.
    """
    names: dict[str, str] = {}
    for row in iter_tsv(_member(opener, 'cohorts.tsv'), required=COHORT_COLUMNS):
        code, name = _cell(row['CANCER']), _cell(row['CANCER_NAME'])
        if code is None or name is None:
            raise ValueError(f'IntOGen cohort {row["COHORT"]} has no cancer type or name')
        if code not in names:
            names[code] = name
            yield {'cancer_type': code, 'name': name}
        elif names[code] != name:
            raise ValueError(f'IntOGen cancer type {code} has two names: {names[code]!r}, {name!r}')
