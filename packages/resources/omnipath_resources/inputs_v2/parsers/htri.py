"""
Read the HTRIdb file.

Source file (semicolon separated, with a header and a ``;`` ending every line;
one experimental observation per line):

    OID;GENEID_TF;SYMBOL_TF;GENEID_TG;SYMBOL_TG;TECHNIQUE;PUBMED_ID;
    32399;142;PARP1;675;BRCA2;Electrophoretic Mobility Shift Assay;18990703;
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_semicolon

COLUMNS = (
    'OID',
    'GENEID_TF',
    'SYMBOL_TF',
    'GENEID_TG',
    'SYMBOL_TG',
    'TECHNIQUE',
    'PUBMED_ID',
)


def iter_observations(
    opener,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield the lines of the HTRIdb file keyed by the header.

    Keys: the seven columns of ``COLUMNS``, as the file has them. A header
    without one of them, a line with another number of fields than the header
    and a file without lines are errors.
    """
    count = 0
    for row in iter_semicolon(opener, strict=True, quotechar=None, **kwargs):
        # The empty name of the column after the last semicolon.
        row.pop('', None)
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(f'HTRIdb file lacks the columns {missing}')
        count += 1
        yield row
    if not count:
        raise ValueError('HTRIdb file has no rows')
