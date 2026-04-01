"""
Parser for STITCH chemical-protein interaction data.

Merges the actions file (directionality, mode) with the links file
(confidence sub-scores) into unified interaction records.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Generator, Iterable
from typing import Any, NamedTuple

from pypath.inputs_v2.base import _first_handle


RE_STITCH_ID = re.compile(r'((?:\d+)?)\.?(CID|ENS[A-Z]*P)([ms]?)(0*)(\d+)')


class StitchEntity(NamedTuple):
    id: str
    type: str
    ncbi_tax_id: int | None
    stereospecific: bool


def _parse_entity(raw_id: str) -> StitchEntity | None:
    """
    Parse a raw STITCH entity ID into a StitchEntity.

    Args:
        raw_id: Raw STITCH identifier string, e.g. '9606.CIDm10461'
            or 'ENSHSAp00000004761'.

    Returns:
        A StitchEntity if the ID matches the expected pattern, otherwise None.
    """
    m = RE_STITCH_ID.match(raw_id)
    if not m:
        return None
    tax, prefix, stereo, zeros, _id = m.groups()
    if prefix == 'CID':
        return StitchEntity(
            id=_id,
            type='small_molecule',
            ncbi_tax_id=None,
            stereospecific=stereo == 's',
        )
    return StitchEntity(
        id=f'{prefix}{zeros}{_id}',
        type='protein',
        ncbi_tax_id=int(tax) if tax else None,
        stereospecific=False,
    )


def _parse_links(handle: Any) -> dict[tuple, dict[str, int]]:
    """
    Parse the STITCH links file into a confidence score lookup.

    Args:
        handle: File handle for the tab-separated links file.

    Returns:
        Dict mapping (chemical_id, protein_id, ncbi_tax_id, stereospecific)
        to a dict of integer confidence sub-scores and combined score.
    """
    lookup = {}
    for row in csv.DictReader(handle, delimiter='\t'):
        chem = _parse_entity(row['chemical'])
        prot = _parse_entity(row['protein'])
        if not chem or not prot:
            continue
        lookup[(chem.id, prot.id, prot.ncbi_tax_id, chem.stereospecific)] = {
            'experimental': int(row['experimental']),
            'prediction': int(row['prediction']),
            'database': int(row['database']),
            'textmining': int(row['textmining']),
            'combined_score': int(row['combined_score']),
        }
    return lookup


def _parse_action_row(row: dict[str, str]) -> dict | None:
    """
    Parse a single row from the STITCH actions file.

    Args:
        row: Dict of raw field values from the actions TSV.

    Returns:
        Parsed action entry dict, or None if either entity ID is
        unrecognised by the regex.
    """
    a = _parse_entity(row['item_id_a'])
    b = _parse_entity(row['item_id_b'])
    if not a or not b:
        return None
    return {
        'a': a,
        'b': b,
        'mode': row['mode'],
        'a_is_acting': row['a_is_acting'].lower() == 't',
        'action': row['action'],
    }


def _iter_deduplicated(
        entries: Iterable[dict],
) -> Generator[dict, None, None]:
    """
    Remove consecutive bidirectional duplicate action entries.

    STITCH lists the same interaction twice for both directions
    (a→b and b→a). For each consecutive pair sharing the same entity
    IDs and mode, retains the directed entry where a_is_acting is True.
    If neither entry is directed, retains the first.

    Args:
        entries: Stream of parsed action entry dicts.

    Yields:
        Deduplicated action entry dicts.
    """
    buf = []
    for entry in entries:
        buf.append(entry)
        if len(buf) == 2:
            pair_0 = (frozenset((buf[0]['a'].id, buf[0]['b'].id)), buf[0]['mode'])
            pair_1 = (frozenset((buf[1]['a'].id, buf[1]['b'].id)), buf[1]['mode'])
            if pair_0 == pair_1:
                if not any(r['a_is_acting'] for r in buf):
                    to_emit = [buf[0]]
                else:
                    to_emit = [r for r in buf if r['a_is_acting']]
                buf = []
            else:
                to_emit = [buf.pop(0)]
            yield from to_emit
    yield from buf


def _iter_unique(
        entries: Iterable[dict],
) -> Generator[dict, None, None]:
    """
    Remove non-consecutive duplicate action entries by entity pair and mode.

    Catches residual duplicates that are not adjacent in the stream and
    therefore not handled by _iter_deduplicated.

    Args:
        entries: Stream of parsed action entry dicts.

    Yields:
        Action entry dicts with all duplicates removed.
    """
    seen = set()
    for entry in entries:
        key = (frozenset((entry['a'].id, entry['b'].id)), entry['mode'])
        if key not in seen:
            seen.add(key)
            yield entry


def _merge_entry(
        entry: dict,
        link_lookup: dict[tuple, dict[str, int]],
        score_threshold: int,
) -> dict[str, str] | None:
    """
    Merge an action entry with its corresponding links confidence scores.

    Determines directionality (chem_is_acting, prot_is_acting) from the
    a_is_acting flag. Both are False when directionality is unknown.

    Args:
        entry: Parsed action entry dict from _parse_action_row.
        link_lookup: Confidence score lookup returned by _parse_links.
        score_threshold: Minimum combined_score required to retain an
            interaction.

    Returns:
        Flat dict of string values for EntityBuilder mapping, or None if
        no matching link exists or the combined score is below the threshold.
    """
    if entry['a'].type == 'small_molecule':
        chem = entry['a']
        prot = entry['b']
    else:
        chem = entry['b']
        prot = entry['a']

    if entry['a_is_acting']:
        chem_is_acting = entry['a'].type == 'small_molecule'
        prot_is_acting = entry['a'].type == 'protein'
    else:
        chem_is_acting = False
        prot_is_acting = False

    link = link_lookup.get((chem.id, prot.id, prot.ncbi_tax_id, chem.stereospecific))
    if link is None or link['combined_score'] <= score_threshold:
        return None

    return {
        'interaction_name': f'{chem.id}_{prot.id}',
        'chemical_id': chem.id,
        'stereospecific': str(chem.stereospecific),
        'protein_id': prot.id,
        'ncbi_tax_id': str(prot.ncbi_tax_id),
        'mode': entry['mode'],
        'action': entry['action'],
        'chem_is_acting': str(chem_is_acting),
        'prot_is_acting': str(prot_is_acting),
        'combined_score': str(link['combined_score']),
        'experimental': f'experimental:{link["experimental"]}',
        'prediction': f'prediction:{link["prediction"]}',
        'database': f'database:{link["database"]}',
        'textmining': f'textmining:{link["textmining"]}',
    }


def iter_stitch_interactions(
        opener: Any,
        score_threshold: int = 700,
        **kwargs: Any,
) -> Generator[dict[str, str], None, None]:
    """
    Parse STITCH interaction data and yield flat interaction records.

    Merges the actions file (directionality, interaction mode) with the
    links file (confidence sub-scores), deduplicates bidirectional entries,
    and filters by combined_score.

    Args:
        opener: Opener object with .links and .actions attributes, each
            wrapping the corresponding STITCH download file handle.
        score_threshold: Minimum combined_score to retain an interaction.
            Defaults to 700.
        **kwargs: Accepted for compatibility with the Dataset interface.

    Yields:
        Flat dict of string values representing one interaction record.
    """
    links_handle = _first_handle(opener.links)
    actions_handle = _first_handle(opener.actions)
    if links_handle is None or actions_handle is None:
        return
    link_lookup = _parse_links(links_handle)
    if not link_lookup:
        return
    parsed_rows = (
        _parse_action_row(row)
        for row in csv.DictReader(actions_handle, delimiter='\t')
    )
    for entry in _iter_unique(_iter_deduplicated(
        r for r in parsed_rows if r is not None
    )):
        rec = _merge_entry(entry, link_lookup, score_threshold)
        if rec is not None:
            yield rec
