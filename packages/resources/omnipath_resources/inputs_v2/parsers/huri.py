"""
Reader for the PSI-MI TAB files of the CCSB human binary interactome datasets.

The portal (interactome-atlas.org) offers every systematic screen of the Center
for Cancer Systems Biology as a PSI-MI TAB 2.7 file with 42 columns and no
header line, one row per ORF pair and assay. The last-but-one comment of a row,
in the column of the interaction annotations, tells in which screens the pair
was found: ``comment:Found in screens 4,7.`` The HuRI file also holds, besides
the pairs of the HuRI map, pairs the portal marks as found in the "test space
screens" (``comment:Found in test space screens 8,3,12.``): those belong to the
separate test space dataset (``Test_space_screens-19.tsv``) and not to the
published map. :func:`iter_screen` can leave them out and counts them.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Generator
import logging
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['MITAB_COLUMNS', 'iter_screen', 'parse_screens']

_log = logging.getLogger(__name__)

# PSI-MI TAB 2.7 columns, in file order. The files have no header line.
MITAB_COLUMNS = (
    'id_a',
    'id_b',
    'alt_ids_a',
    'alt_ids_b',
    'aliases_a',
    'aliases_b',
    'detection_method',
    'first_author',
    'publication',
    'taxon_a',
    'taxon_b',
    'interaction_type',
    'source_database',
    'interaction_id',
    'confidence',
    'expansion',
    'biological_role_a',
    'biological_role_b',
    'experimental_role_a',
    'experimental_role_b',
    'interactor_type_a',
    'interactor_type_b',
    'xrefs_a',
    'xrefs_b',
    'xrefs_interaction',
    'annotations_a',
    'annotations_b',
    'annotations_interaction',
    'host_organism',
    'parameters',
    'creation_date',
    'update_date',
    'checksum_a',
    'checksum_b',
    'checksum_interaction',
    'negative',
    'features_a',
    'features_b',
    'stoichiometry_a',
    'stoichiometry_b',
    'identification_method_a',
    'identification_method_b',
)
_SCREENS_RE = re.compile(r'comment:Found in (test space )?screens ([\d,]*)\.')


def parse_screens(cell: object) -> tuple[list[str], bool]:
    """Read the interaction annotation ``comment:Found in screens 4,7.``.

    Returns the screen numbers as written (an empty list where the comment has
    none) and whether the pair was found in the test space screens. Any other
    cell is an error: the layout of the file changed.
    """
    match = _SCREENS_RE.fullmatch(str(cell if cell is not None else '').strip())
    if match is None:
        raise ValueError(f'Unexpected interaction annotation in a CCSB file: {cell!r}')
    numbers = [number for number in match.group(2).split(',') if number]
    return numbers, bool(match.group(1))


def iter_screen(
    opener,
    exclude_test_space: bool = False,
    stats: Counter | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, str], None, None]:
    """Yield the rows of a PSI-MI TAB file of the portal, keyed by :data:`MITAB_COLUMNS`.

    A line with another number of fields is an error. With
    ``exclude_test_space`` the rows of pairs found in the test space screens
    are not yielded; ``stats``, if given, is filled with the counts ``rows``
    and ``test_space_screens_excluded``, which are also logged.
    """
    stats = Counter() if stats is None else stats
    for row in iter_tsv(opener, fieldnames=MITAB_COLUMNS, strict=True):
        stats['rows'] += 1
        if exclude_test_space and parse_screens(row['annotations_interaction'])[1]:
            stats['test_space_screens_excluded'] += 1
            continue
        yield row
    _log.info('CCSB PSI-MI TAB file: %s', dict(sorted(stats.items())))
