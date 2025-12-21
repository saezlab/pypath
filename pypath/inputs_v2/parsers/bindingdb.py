"""
BindingDB TSV parser.

Parses BindingDB's tab-separated data files with column filtering.
"""

from __future__ import annotations

from collections.abc import Generator
import csv


def _raw(opener, max_lines: int | None = None, **_kwargs: object) -> Generator[dict[str, str], None, None]:
    """
    Parse BindingDB TSV file, keeping only the first 49 columns.

    Args:
        opener: File opener from download_and_open
        max_lines: Optional limit on number of records to process

    Yields:
        Dictionary for each binding record
    """
    if not opener or not opener.result:
        return
    for file_handle in opener.result.values():
        header_line = file_handle.readline().strip()
        header = header_line.split('\t')

        columns_to_keep = min(49, len(header))
        filtered_header = header[:columns_to_keep]

        def filtered_rows():
            for line in file_handle:
                columns = line.strip().split('\t')
                yield '\t'.join(columns[:columns_to_keep])

        reader = csv.DictReader(filtered_rows(), fieldnames=filtered_header, delimiter='\t')

        if max_lines:
            for i, row in enumerate(reader):
                if i >= max_lines:
                    break
                yield row
        else:
            yield from reader
        break
