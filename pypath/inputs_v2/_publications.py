"""
The form of a PubMed identifier, shared by the inputs.

A PubMed identifier is a whole number of one to eight digits with no leading
zero. Parsers use the rule to tell an identifier from other content of a
reference field; mappers use it to build the ``PMID:`` CURIE. A number of
this form that a source has cut or run together is still not the paper the
source means: a parser that knows such a defect of its file excludes those
values itself and says why.
"""

from __future__ import annotations

import re

__all__ = ['PMID_PATTERN', 'PMID_EXTRACT', 'is_pmid', 'pmid_curie']

# Without anchors or a group, for use inside a larger expression.
PMID_PATTERN = r'[1-9]\d{0,7}'
# A whole value that is one identifier, as a ``FieldConfig`` extract pattern.
PMID_EXTRACT = rf'^({PMID_PATTERN})$'

_PMID = re.compile(PMID_PATTERN)


def is_pmid(value: object) -> bool:
    """Whether the value, as written, is one PubMed identifier."""
    return isinstance(value, (str, int)) and not isinstance(value, bool) and bool(
        _PMID.fullmatch(str(value).strip())
    )


def pmid_curie(value: object) -> str | None:
    """``PMID:123`` for a PubMed identifier; None for anything else."""
    return f'PMID:{str(value).strip()}' if is_pmid(value) else None
