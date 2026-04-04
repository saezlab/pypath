"""BindingDB TSV parser (2026-03+ schema)."""

from __future__ import annotations

from collections.abc import Generator
import csv
import re


_CHEMBL_RUN_RE = re.compile(r'(CHEMBL\d+)(?=CHEMBL\d+)')


def _normalize_row(row: dict[str, str]) -> dict[str, str]:
    value = (row.get('ChEMBL ID of Ligand') or '').strip()
    if value and value.count('CHEMBL') > 1 and '::' not in value and ';' not in value and '|' not in value:
        row['ChEMBL ID of Ligand'] = _CHEMBL_RUN_RE.sub(r'\1::', value)
    return row


def _raw(opener, max_lines: int | None = None, **_kwargs: object) -> Generator[dict[str, str], None, None]:
    """Parse BindingDB TSV rows.

    Expects the current BindingDB schema with chain-indexed target columns
    (e.g. ``UniProt ... Target Chain 1``).
    """
    if not opener or not opener.result:
        return

    for file_handle in opener.result.values():
        reader = csv.DictReader(file_handle, delimiter='\t')

        if max_lines:
            for i, row in enumerate(reader):
                if i >= max_lines:
                    break
                yield _normalize_row(row)
        else:
            for row in reader:
                yield _normalize_row(row)

        break
