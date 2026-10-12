"""
Parser of the ENCODE-Nets network files of the Gerstein Lab.

The networks that accompany "Architecture of the human regulatory network
derived from ENCODE data" (Gerstein et al., Nature 2012) are published on
encodenets.gersteinlab.org as text files without a header, one edge per line
with three columns:

    CTCF proximal_filtered MYBPC3          (enets2.Proximal_filtered.txt)
    ATF3 distal ICD7                       (enets3.Distal.txt)
    BDP1 (TF-miRNA) hsa-miR-19b-1*         (enets10.TF-miRNA.txt)

 1. ``source``: the regulator, a transcription factor, or in the lines of
    ``(miRNA-TF)`` a microRNA
 2. ``edge_type``: the network (``proximal_filtered``, ``distal``) or, in the
    microRNA file, the kind of edge
 3. ``target``: the regulated gene, microRNA or transcription factor

The columns are separated by a space in the proximal and the microRNA file and
by a tab in the distal file (which has Windows line ends). The microRNA file
holds two kinds of edges: ``(TF-miRNA)``, a transcription factor with a peak
near a microRNA, and ``(miRNA-TF)``, a microRNA that targets the transcription
factor, which the paper takes from TargetScan.
"""

from __future__ import annotations

from collections.abc import Generator
import functools
from typing import Any

from pypath.inputs_v2.parsers.base import iter_csv

COLUMNS = ('source', 'edge_type', 'target')
PROXIMAL = 'proximal_filtered'
DISTAL = 'distal'
TF_MIRNA = '(TF-miRNA)'
MIRNA_TF = '(miRNA-TF)'


def iter_edges(
    opener,
    *,
    delimiter: str,
    edge_type: str,
    other_edge_types: tuple[str, ...] = (),
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    The lines of a network file whose edge type is ``edge_type``.

    The keys are ``COLUMNS``. A line of an edge type in ``other_edge_types`` is
    left out; any other edge type, a line that does not have three columns, or
    a file with no line of the edge type is an error.
    """
    count = 0
    for row in iter_csv(
        opener,
        delimiter=delimiter,
        fieldnames=COLUMNS,
        strict=True,
        quotechar=None,
    ):
        record = {key: value.strip() for key, value in row.items()}
        if record['edge_type'] in other_edge_types:
            continue
        if record['edge_type'] != edge_type:
            raise ValueError(
                f"Unexpected ENCODE-Nets edge type {record['edge_type']!r}, "
                f'expected {edge_type!r}'
            )
        count += 1
        yield record
    if not count:
        raise ValueError(f'ENCODE-Nets file has no {edge_type!r} edge')


iter_proximal_filtered = functools.partial(
    iter_edges, delimiter=' ', edge_type=PROXIMAL
)
iter_distal = functools.partial(iter_edges, delimiter='\t', edge_type=DISTAL)
iter_tf_mirna = functools.partial(
    iter_edges,
    delimiter=' ',
    edge_type=TF_MIRNA,
    other_edge_types=(MIRNA_TF,),
)
