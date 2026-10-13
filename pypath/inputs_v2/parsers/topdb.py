"""Read the TOPDB transmembrane protein entries.

TOPDB (the Topology Data Bank of Transmembrane Proteins, now served by UniTmp)
publishes its entries as one XML document (``topdb_all.xml``, release 3.3, all
species). An entry is a ``TOPDB`` element:

    <TOPDB ID="ADA2C_HUMAN" type="Alpha_Polytopic">
      <Name>Alpha-2C adrenergic receptor</Name>
      <CrossRef><PDB ID="6kuw">...</PDB></CrossRef>
      <Sequence Len="462" Checksum="...">...</Sequence>
      <Membranes><Membrane>Cell membrane;Cytoplasm;Extracellular</Membrane></Membranes>
      <Topology>
        <Reliability>70.6</Reliability><Numtm Count="6"/>
        <Regions><Region Begin="1" End="53" Loc="Outside"/>...</Regions>
      </Topology>
      <Experiments>...</Experiments>
    </TOPDB>

The ``ID`` is the UniProt entry name, or for entries that come from a PDB
structure without a UniProt entry the PDB code and the chain (``2k3c_A``).
Some entries are in the file twice, identically.

:func:`iter_entries` yields the fields the module uses, one dictionary per
distinct entry, and :func:`iter_terms` the membranes and types the entries
have. Of the topology the number of transmembrane segments (``Numtm``) and the
reliability are read, and the number is checked against the ``Membrane`` regions
of the entry; the regions themselves, the sequence, the cross-references and the
experiments are not read.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable
import math
import re
from typing import Any

from lxml import etree

__all__ = ['iter_entries', 'iter_terms']

_NS = '{https://topdb.unitmp.org}'
_TOPDB = f'{_NS}TOPDB'
_ENTRY_NAME = re.compile(r'[A-Z0-9]+_[A-Z0-9]+')
# A PDB structure without a UniProt entry: the PDB code (lower case) and a chain.
_PDB_CHAIN = re.compile(r'[0-9][a-z0-9]{3}_[A-Za-z0-9]+')


def _chunks(result: Any) -> Iterable[bytes | str]:
    if hasattr(result, 'read'):
        while chunk := result.read(1 << 20):
            yield chunk
    else:
        yield from result


def _elements(opener) -> Generator[etree._Element, None, None]:
    """The ``TOPDB`` elements of the document, each released after it is used."""
    result = getattr(opener, 'result', None) if opener else None
    if result is None:
        return
    if isinstance(result, dict):
        result = next(iter(result.values()))
    parser = etree.XMLPullParser(events=('end',), tag=_TOPDB)
    for chunk in _chunks(result):
        parser.feed(chunk)
        for _, element in parser.read_events():
            yield element
            element.clear()
            while element.getprevious() is not None:
                del element.getparent()[0]
    parser.close()


def _membrane(text: str | None) -> dict[str, str]:
    """A membrane with the two compartments on its sides, as the cell writes it:
    ``Cell membrane;Cytoplasm;Extracellular``."""
    parts = [part.strip() for part in (text or '').split(';')]
    if len(parts) != 3 or not all(parts):
        raise ValueError(f'TOPDB membrane is not "name;side;side": {text!r}')
    return {'name': parts[0], 'side_1': parts[1], 'side_2': parts[2]}


def _topology(element: etree._Element, identifier: str) -> tuple[int, str]:
    """The number of transmembrane segments and the reliability of an entry.

    ``Numtm`` must equal the number of ``Membrane`` regions of the topology, and
    the reliability must be a number from 0 to 100 (the average of the posterior
    probabilities along the predicted path, in percent, as the site defines it).
    """
    topology = element.find(f'{_NS}Topology')
    if topology is None:
        raise ValueError(f'TOPDB entry {identifier} has no topology')
    numtm = topology.find(f'{_NS}Numtm')
    reliability = topology.find(f'{_NS}Reliability')
    if numtm is None or reliability is None:
        raise ValueError(f'TOPDB entry {identifier} has no Numtm or no Reliability')
    try:
        segments = int(numtm.attrib['Count'])
        value = float((reliability.text or '').strip())
    except (KeyError, ValueError):
        raise ValueError(
            f'TOPDB entry {identifier} has no valid Numtm Count or Reliability'
        ) from None
    if segments < 0 or not math.isfinite(value) or not 0 <= value <= 100:
        raise ValueError(f'TOPDB entry {identifier} has Numtm {segments} and Reliability {value}')
    regions = topology.findall(f'{_NS}Regions/{_NS}Region')
    membrane = sum(region.attrib.get('Loc') == 'Membrane' for region in regions)
    if segments != membrane:
        raise ValueError(
            f'TOPDB entry {identifier} has Numtm {segments} but {membrane} Membrane regions'
        )
    return segments, reliability.text.strip()


def _entry(element: etree._Element) -> dict[str, Any]:
    identifier = element.attrib.get('ID')
    kind = element.attrib.get('type')
    if not identifier or not kind:
        raise ValueError('TOPDB entry without an ID or a type')
    segments, reliability = _topology(element, identifier)
    name = element.find(f'{_NS}Name')
    return {
        'id': identifier,
        'type': kind,
        'name': (name.text or '').strip() or None if name is not None else None,
        'membranes': [
            _membrane(membrane.text)
            for membrane in element.findall(f'{_NS}Membranes/{_NS}Membrane')
        ],
        'segments': segments,
        'reliability': reliability,
    }


def iter_entries(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the distinct entries of the document.

    The fields: ``id`` (see above), ``kind`` as ``uniprot_entry`` (an entry
    name) or ``pdb_chain``, ``type`` (the ``type`` attribute: ``Alpha_Polytopic``,
    ``Alpha_Bitopic`` or ``Beta_Polytopic``), ``name`` (the protein name) and
    ``membranes``, a list of dictionaries with ``name``, ``side_1`` and
    ``side_2`` (empty for the entries without a membrane), ``segments`` (the
    ``Numtm`` count, an integer) and ``reliability`` (the topology's
    ``Reliability``, as written).

    An entry that is in the file twice is yielded once; two entries with one ID
    and different content, an ID that is neither an entry name nor a PDB chain,
    an entry without a topology, a ``Numtm`` that is not the number of
    ``Membrane`` regions of the entry, a reliability that is not a number from
    0 to 100, and a document without entries are errors.
    """
    seen: dict[str, str] = {}
    count = 0
    for element in _elements(opener):
        entry = _entry(element)
        count += 1
        signature = repr(entry)
        if entry['id'] in seen:
            if seen[entry['id']] != signature:
                raise ValueError(f'TOPDB entry {entry["id"]} is repeated with other content')
            continue
        seen[entry['id']] = signature
        if _ENTRY_NAME.fullmatch(entry['id']):
            entry['kind'] = 'uniprot_entry'
        elif _PDB_CHAIN.fullmatch(entry['id']):
            entry['kind'] = 'pdb_chain'
        else:
            raise ValueError(f'TOPDB entry ID is not an entry name or a PDB chain: {entry["id"]!r}')
        yield entry
    if not count:
        raise ValueError('TOPDB document has no entries')


def iter_terms(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the distinct membranes and types of the entries, membranes first.

    The membrane rows have ``kind`` as ``membrane``, ``id`` and ``name`` (the
    membrane's name), ``side_1`` and ``side_2``; the type rows ``kind`` as
    ``type``, ``id`` and ``name`` (the attribute's value). A membrane name that
    comes with different sides is an error.
    """
    membranes: dict[str, dict[str, str]] = {}
    types: dict[str, dict[str, str]] = {}
    for entry in iter_entries(opener):
        types.setdefault(
            entry['type'], {'kind': 'type', 'id': entry['type'], 'name': entry['type']}
        )
        for membrane in entry['membranes']:
            previous = membranes.setdefault(
                membrane['name'],
                {'kind': 'membrane', 'id': membrane['name'], **membrane},
            )
            if (previous['side_1'], previous['side_2']) != (
                membrane['side_1'], membrane['side_2'],
            ):
                raise ValueError(f'TOPDB membrane {membrane["name"]!r} has two side definitions')
    yield from sorted(membranes.values(), key=lambda row: row['id'])
    yield from sorted(types.values(), key=lambda row: row['id'])
