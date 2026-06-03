#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

from __future__ import annotations

import json
import urllib.parse
from typing import Iterable, Iterator, NamedTuple

import pypath.resources.urls as urls
import pypath.share.session as session
from pypath.share.downloads import dm, download_and_open

_logger = session.Logger(name = 'pubchem_input')
_log = _logger._log


class PubchemMapping(NamedTuple):
    """A single ``CID -> value`` row from a PubChem Compound ``Extras`` table."""

    cid: str
    value: str


# Our id-type name -> (PubChem ``Extras`` file token, value column index).
# The source is always the Compound CID (column 0). The combined
# ``CID-InChI-Key`` table carries the InChI (column 1) and the standard
# InChIKey (column 2) side by side, so both map to the same file.
_PUBCHEM_COMPOUND_TARGETS: dict[str, tuple[str, int]] = {
    'inchikey': ('InChI-Key', 2),
    'inchi': ('InChI-Key', 1),
    'smiles': ('SMILES', 1),
    'iupac': ('IUPAC-Name', 1),
    'synonym': ('Synonym-unfiltered', 1),
    'parent-cid': ('Parent', 1),
    'component-cid': ('Component', 1),
    'preferred-cid': ('Preferred', 1),
}


def pubchem_mapping(
        target: str,
        source: str = 'cid',
    ) -> Iterator[PubchemMapping]:
    """
    Stream identifier/structure translation data from PubChem.

    Yields ``(cid, value)`` rows from the PubChem FTP "Extras" tables
    (``Compound/Extras/CID-<TYPE>.gz``) -- e.g. CID -> standard InChIKey,
    InChI or canonical SMILES. The table is fetched and cached via the
    shared dlmachine download manager, then read line by line (no in-memory
    table and no SQLite side-cache), so a consumer can ``itertools.islice``
    for a capped subset.

    Args
        target (str): The target ID/structure type. One of the keys of
            ``_PUBCHEM_COMPOUND_TARGETS`` (``inchikey``, ``inchi``, ``smiles``,
            ``iupac``, ``synonym``, ``parent-cid``, ``component-cid``,
            ``preferred-cid``) or a raw PubChem ``Extras`` file token.
        source (str): The source ID type. Only the PubChem Compound CID
            (``cid``) is supported.

    Yields
        PubchemMapping: ``(cid, value)`` named tuples.
    """

    if source not in {'cid', 'pubchem-cid'}:

        msg = (
            'Only the PubChem Compound (CID) source is supported, not `%s`.'
            % source
        )
        _log(msg)
        raise ValueError(msg)

    file_token, value_col = _PUBCHEM_COMPOUND_TARGETS.get(target, (target, 1))

    url = urls.urls['pubchem']['ftp'] % ('Compound', 'CID', file_token)
    opener = download_and_open(
        url,
        filename = 'CID-%s.gz' % file_token,
        subfolder = 'pubchem',
        large = True,
        ext = 'gz',
    )

    for line in opener.result:

        if isinstance(line, bytes):
            line = line.decode('utf-8', errors = 'replace')

        fields = line.rstrip('\n').split('\t')

        if len(fields) <= value_col:
            continue

        cid = fields[0].strip()
        value = fields[value_col].strip()

        if cid and value:

            yield PubchemMapping(cid = cid, value = value)


def pubchem_name_cids(name: str) -> set[str]:
    """
    PubChem CIDs for a compound name via the PUG REST API.

    Each name–URL response is fetched and cached on disk by the shared
    dlmachine download manager, so repeated calls across sessions incur no
    HTTP requests.  Unknown names (PubChem 404) and network errors return an
    empty set.

    Args:
        name: Compound name or synonym (e.g. ``'ATP'``, ``'NAD+'``).

    Returns:
        Set of PubChem CID strings matching *name*.  Empty set if the name
        is not found or the request fails.
    """

    url = urls.urls['pubchem']['name_cids'] % urllib.parse.quote(name)

    try:
        path = dm.download(url)
    except Exception:
        return set()

    if not path:
        return set()

    try:
        with open(path, encoding = 'utf-8') as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        return set()

    # PubChem returns {"Fault": {...}} with HTTP 4xx for unknown names.
    if not isinstance(data, dict) or 'Fault' in data:
        return set()

    return {
        str(cid)
        for cid in data.get('IdentifierList', {}).get('CID', [])
    }


def pubchem_names_cids(
        names: Iterable[str],
        show_progress: bool = True,
) -> dict[str, set[str]]:
    """
    PubChem CIDs for a collection of compound names.

    Calls :func:`pubchem_name_cids` once per unique name.  Each result is
    cached on disk by pypath's curl infrastructure; subsequent calls for the
    same names incur no HTTP requests.

    Args:
        names: Iterable of compound name strings.
        show_progress: Show a tqdm progress bar (default ``True``).
            Silently ignored if tqdm is not installed.

    Returns:
        Dict mapping each input name to its set of PubChem CID strings.
        Names not found in PubChem map to an empty set.
    """

    names = list(names)

    try:
        from tqdm import tqdm
        iterator = tqdm(
            names,
            desc = 'PubChem name lookup',
            unit = 'name',
            disable = not show_progress,
        )
    except ImportError:
        iterator = names

    return {name: pubchem_name_cids(name) for name in iterator}
