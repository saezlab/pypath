"""
Parser of the InnateDB interactions in PSI-MI TAB 2.7.

The file has a header line, whose first name starts with a hash
(``#unique_identifier_A``), and one row per experimental observation of an
interaction. Rows are dictionaries keyed by the header names; every value is
text as written, with ``-`` for no value. Cells that hold several items
separate them with ``|`` (see ``pypath.inputs_v2._mitab``).
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2 import _mitab
from pypath.inputs_v2.parsers.base import iter_tsv

# The columns the mapper reads, of the 42 of PSI-MI TAB 2.7. The header line
# starts with a hash.
COLUMNS = (
    '#unique_identifier_A',
    'unique_identifier_B',
    'alt_identifier_A',
    'alt_identifier_B',
    'alias_A',
    'alias_B',
    'ncbi_taxid_A',
    'ncbi_taxid_B',
    'biological_role_A',
    'biological_role_B',
    'exp_role_A',
    'exp_role_B',
    'interactor_type_A',
    'interactor_type_B',
    'participant_identification_method_A',
    'participant_identification_method_B',
    'interaction_detection_method',
    'pmid',
    'interaction_type',
    'source_database',
    'idinteraction_in_source_db',
    'ncbi_taxid_host_organism',
    'negative',
)


def iter_interactions(
    opener, **kwargs: Any
) -> Generator[dict[str, str], None, None]:
    """Rows of the file; another layout is an error, not an empty dataset."""
    yield from _mitab.checked_rows(
        iter_tsv(opener, **kwargs), COLUMNS, 'InnateDB interactions'
    )
