"""
Parser of the LMPID XML export with literature references.

The document is a ``RECORDSET`` of ``RECORD`` elements, one per linear motif
instance and interaction partner::

    <RECORD>
    <Bait_UniProt_ID>P54253</Bait_UniProt_ID>
    <Motif_Instance>RRWSAP</Motif_Instance>
    <Sequence_Position>772-777</Sequence_Position>
    <Interacting_Domain>14-3-3 domain</Interacting_Domain>
    <Prey_UniProt_ID>P92177</Prey_UniProt_ID>
    <References>20037628, 12757707</References>
    </RECORD>

The bait is the protein that carries the motif; the prey is the protein
whose domain binds it.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any
import xml.etree.ElementTree as ET

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import read_bytes

_FIELDS = {
    'Bait_UniProt_ID': 'bait_uniprot',
    'Prey_UniProt_ID': 'prey_uniprot',
    'Motif_Instance': 'motif_instance',
    'Sequence_Position': 'sequence_position',
    'Interacting_Domain': 'interacting_domain',
    'References': 'references',
}
# The source writes this where the binding domain is not known.
_UNKNOWN = 'Unknown'
_POSITION_RE = re.compile(r'(\d+)-(\d+)')


def iter_lmpid(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    Yield one dictionary per ``RECORD``.

    Fields, all from the element of the corresponding name:

    - ``bait_uniprot``, ``prey_uniprot``: UniProt accessions, whitespace removed.
    - ``motif_instance``: the motif's residues in the bait.
    - ``sequence_position``: the source's ``start-end`` text;
      ``motif_start`` and ``motif_end`` are its two numbers, or None when
      the text is not such a range. The coordinates are the source's own.
    - ``interacting_domain``: name of the prey's domain; None for ``Unknown``.
    - ``references``: the source text; ``pmids`` are its comma-separated
      PubMed IDs, in source order, without entries that are not PubMed IDs.

    An element the record lacks, or leaves empty, is None.
    """
    data = read_bytes(opener)
    if not data.strip():
        return
    root = ET.fromstring(data)
    records = root.findall('RECORD')
    if not records:
        raise ValueError(
            f'LMPID: no RECORD elements under <{root.tag}>; '
            'the download is not the expected XML export'
        )
    for record in records:
        row: dict[str, Any] = dict.fromkeys(_FIELDS.values())
        for element in record:
            field = _FIELDS.get(element.tag)
            if field:
                row[field] = (element.text or '').strip() or None
        if row['interacting_domain'] == _UNKNOWN:
            row['interacting_domain'] = None
        position = _POSITION_RE.fullmatch(row['sequence_position'] or '')
        row['motif_start'] = int(position.group(1)) if position else None
        row['motif_end'] = int(position.group(2)) if position else None
        references = (row['references'] or '').split(',')
        row['pmids'] = [
            reference.strip()
            for reference in references
            if is_pmid(reference)
        ]
        yield row
