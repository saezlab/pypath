"""
Reader for the DIP interaction file in PSI-MI TAB 2.5.

DIP (Database of Interacting Proteins) distributes its interactions as a
tab separated file with a header row. The 16 columns are those of MITAB 2.5
plus the DIP columns ``Confidence value(s)`` (the quality status) and
``Processing Status``. The header ends with an empty name, a seventeenth
column that holds nothing, and every row has an eighteenth field after it,
``-``, that the header does not name.

    ID interactor A     DIP-617N|refseq:NP_000607|uniprotkb:P01730
    Interaction type(s) MI:0407(direct interaction)|MI:0407(direct interaction)
    Confidence value(s) dip-quality-status:core
    Processing Status   dip:0004(small scale)|dip:0004(small scale)

A row is one DIP interaction (``DIP-42E``) between two DIP nodes (``DIP-617N``).
Four of its columns list the evidence items of the interaction, one item each
and in the same order: ``Interaction detection method(s)``, ``Interaction
type(s)``, ``Processing Status`` and ``Publication Identifier(s)``, which has
two entries for an item, the PubMed ID and a DIP record of the article
(``pubmed:DIP-17838S``). Terms of the PSI-MI ontology are written without the
database name ``psi-mi`` (``MI:0114(x-ray crystallography)``), and the DIP
identifiers without a database name. The cells of other PSI-MI TAB sources are
read by ``pypath.inputs_v2._mitab``.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2 import _mitab
from pypath.inputs_v2.parsers.base import iter_tsv

# Header of the file and the field each column is read into, in file order.
COLUMNS = (
    ('ID interactor A', 'id_a'),
    ('ID interactor B', 'id_b'),
    ('Alt. ID interactor A', 'alt_ids_a'),
    ('Alt. ID interactor B', 'alt_ids_b'),
    ('Alias(es) interactor A', 'aliases_a'),
    ('Alias(es) interactor B', 'aliases_b'),
    ('Interaction detection method(s)', 'detection_methods'),
    ('Publication 1st author(s)', 'first_authors'),
    ('Publication Identifier(s)', 'publications'),
    ('Taxid interactor A', 'taxid_a'),
    ('Taxid interactor B', 'taxid_b'),
    ('Interaction type(s)', 'interaction_types'),
    ('Source database(s)', 'source_databases'),
    ('Interaction identifier(s)', 'interaction_ids'),
    ('Confidence value(s)', 'confidence'),
    ('Processing Status', 'processing_status'),
    # The header names a seventeenth column with an empty string.
    ('', 'unnamed'),
)
# The one quality status the file has a name for: the other rows have the
# database name and no value.
_QUALITY_PREFIX = 'dip-quality-status:'
QUALITY_STATUSES = frozenset({'core'})


def iter_dip(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per row of the file.

    Every cell is text, as written but for the whitespace around it (the
    column ``Alt. ID interactor A`` starts with a space); ``-`` stands for no
    value. A file with other columns, a row with another number of fields or
    without the empty seventeenth column, evidence columns that do not list the
    same number of items, a quality status other than ``core`` or none, or no
    row at all is an error.

    The dictionary has the fields named in ``COLUMNS`` (the last one is left
    out), and:

    - ``quality_status``: ``core`` for ``dip-quality-status:core``; None for
      ``dip-quality-status:`` without a value.
    - ``evidence_count``: the number of evidence items of the row.
    """
    count = 0
    for row in iter_tsv(opener, quotechar=None):
        if tuple(name for name in row if name is not None) != tuple(
            name for name, _ in COLUMNS
        ):
            raise ValueError(
                f'DIP file has columns {[name for name in row if name is not None]}'
            )
        # The eighteenth field of every row; the shared reader keeps the fields
        # beyond the header in a list under the key None.
        if row.get(None) != ['-'] or row[''] != '':
            raise ValueError(
                f'DIP row {count + 2} has other fields beyond the header: '
                f'{row.get(None)!r}'
            )
        record = {
            field: (row[name] or '').strip() for name, field in COLUMNS[:-1]
        }
        confidence = record['confidence']
        if not confidence.startswith(_QUALITY_PREFIX):
            raise ValueError(f'DIP quality status not understood: {confidence!r}')
        status = confidence[len(_QUALITY_PREFIX) :]
        if status and status not in QUALITY_STATUSES:
            raise ValueError(f'Unknown DIP quality status: {confidence!r}')
        record['quality_status'] = status or None
        # One entry for each item, and two for each item in the publications.
        items = [
            len(_mitab.parse_mitab_cell(record[field]))
            for field in ('detection_methods', 'interaction_types', 'processing_status')
        ]
        publications = len(_mitab.parse_mitab_cell(record['publications']))
        if len(set(items)) != 1 or publications != 2 * items[0] or not items[0]:
            raise ValueError(
                f'DIP {record["interaction_ids"]}: evidence columns list '
                f'{items} items and {publications} publications'
            )
        record['evidence_count'] = items[0]
        count += 1
        yield record
    if not count:
        raise ValueError('DIP file has no rows')
