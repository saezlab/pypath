"""
Raw parsers for TCDB (Transporter Classification Database).

Two parsers are provided:

- transporters: parses acc2tcid.py (UniProt → TCID) as the primary download
  and fetches families.py (family prefix → family name) as a secondary
  download injected via functools.partial.

- substrates: parses getSubstrates.py (TCID → CHEBI substrates) as the
  primary download and fetches acc2tcid.py (UniProt → TCID) as a secondary
  download injected via functools.partial to map TCIDs to UniProt accessions.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from collections.abc import Generator
from typing import Any


_FAMILY_PREFIX_RE = re.compile(r'(\d+\.[A-Z]+\.\d+)')


def _parse_families(handle) -> dict[str, str]:
    """Parse families.py into a family prefix to family name mapping.

    Each line in families.py has the form::

        family_prefix<TAB>family name

    Args:
        handle: An iterable of text lines from the families.py file.

    Returns:
        A dict mapping each family prefix (e.g. ``'1.A.1'``) to its cleaned
        family name string.
    """
    families: dict[str, str] = {}
    for row in csv.reader(handle, delimiter='\t'):
        if len(row) != 2:
            continue
        prefix, name = row
        families[prefix] = name
    return families


def _build_tcid_to_uniprots(handle) -> dict[str, list[str]]:
    """Build a TCID-to-UniProt index from an acc2tcid.py file handle.

    Each line in acc2tcid.py has the form::

        uniprot_accession<TAB>full_tcid

    Args:
        handle: An iterable of text lines from the acc2tcid.py file.

    Returns:
        A dict mapping each TC number to the list of UniProt accessions
        that carry that classification.
    """
    result: dict[str, list[str]] = defaultdict(list)
    for row in csv.reader(handle, delimiter='\t'):
        if len(row) != 2:
            continue
        result[row[1]].append(row[0])
    return dict(result)


def _family_prefix(tcid: str) -> str:
    """Extract the three-component family prefix from a full TCID.

    For example ``'1.A.1.1.1'`` → ``'1.A.1'``.

    Args:
        tcid: A full TC number string.

    Returns:
        The family prefix string, or an empty string if the pattern is not
        found.
    """
    m = _FAMILY_PREFIX_RE.search(tcid)
    return m.group(1) if m else ''


def transporters(
    opener,
    *,
    families_download,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Parse acc2tcid.py and yield one record per (UniProt, TCID) pair.

    Each record includes the full TC number and the transporter family name
    looked up from families.py via the injected ``families_download``.

    Args:
        opener: Download opener for acc2tcid.py (UniProt → TCID mapping).
        families_download: A :class:`~pypath.inputs_v2.base.Download` instance
            for families.py, injected via ``functools.partial`` in the main
            ``tcdb`` module. Used to resolve family prefixes to family names.
        **kwargs: Passed through; ``force_refresh`` is forwarded to the
            secondary families.py download.

    Yields:
        dict with keys ``uniprot``, ``tcid``, and ``family_name``.
    """
    if not opener or not opener.result:
        return

    force_refresh = kwargs.get('force_refresh', False)
    families_opener = families_download.open(force_refresh=force_refresh)
    if not families_opener or not families_opener.result:
        return

    families = _parse_families(families_opener.result)

    for tcid, uniprots in _build_tcid_to_uniprots(opener.result).items():
        prefix = _family_prefix(tcid)
        family_name = families.get(prefix, '')
        for uniprot in uniprots:
            yield {
                'uniprot': uniprot,
                'tcid': tcid,
                'family_name': family_name,
            }


def substrates(
    opener,
    *,
    acc2tc_download,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Parse getSubstrates.py and yield one record per transporter-substrate pair.

    The substrates file maps TC numbers to CHEBI-identified small molecules.
    A secondary fetch of acc2tcid.py (via ``acc2tc_download``) is performed
    to resolve TC numbers to UniProt accessions.

    Args:
        opener: Download opener for getSubstrates.py. Each line has the
            form ``TCID<TAB>CHEBI:id;name|CHEBI:id;name|...``.
        acc2tc_download: A :class:`~pypath.inputs_v2.base.Download` instance
            for acc2tcid.py, injected via ``functools.partial`` in the main
            ``tcdb`` module. Used to build the TCID → UniProt mapping.
        **kwargs: Passed through; ``force_refresh`` is forwarded to the
            secondary acc2tcid.py download.

    Yields:
        dict with keys ``tcid``, ``transporter_uniprot``, ``substrate_id``
        (e.g. ``'CHEBI:24870'``), and ``substrate_name``.
    """
    if not opener or not opener.result:
        return

    force_refresh = kwargs.get('force_refresh', False)
    acc2tc_opener = acc2tc_download.open(force_refresh=force_refresh)
    if not acc2tc_opener or not acc2tc_opener.result:
        return

    tcid_to_uniprots = _build_tcid_to_uniprots(acc2tc_opener.result)

    for row in csv.reader(opener.result, delimiter='\t'):
        if len(row) != 2:
            continue

        tcid, substrates_raw = row
        uniprots = tcid_to_uniprots.get(tcid, [])

        for entry in substrates_raw.split('|'):
            entry = entry.strip()
            if not entry:
                continue

            substrate_id, substrate_name = entry.split(';', maxsplit=1)

            for uniprot in uniprots:
                yield {
                    'tcid': tcid,
                    'transporter_uniprot': uniprot,
                    'substrate_id': substrate_id,
                    'substrate_name': substrate_name,
                }
