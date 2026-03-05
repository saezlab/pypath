"""BindingDB TSV parser (2026-03+ schema)."""

from __future__ import annotations

from collections.abc import Generator
import csv


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
                yield row
        else:
            yield from reader

        break
