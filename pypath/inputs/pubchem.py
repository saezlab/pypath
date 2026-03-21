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
import os
import urllib.parse
from typing import Iterable

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.cache as cache
import pypath.share.lookup as lookup
import pypath.share.session as session

_logger = session.Logger(name = 'pubchem_input')
_log = _logger._log


def pubchem_mapping(target, source = 'cid'):
    """
    Identifier translation data from PubChem.

    Args
        target (str): The target ID type, either as it is used in the file
            names in the PubChem FTP service or as simpler, all lowercase
            strings used in this module. Possible values are parent-cid,
            component-cid, inchi, iupac, preferred-cid, sid, smiles,
            synonym.
        source (str): The source ID type. Either sid or cid.

    Returns
        (dict): A dict of sets with the source identifiers as keys and sets
            of target identifiers as values.
    """

    id_types = {
        'parent-cid': 'Parent',
        'component-cid': 'Component',
        'inchi': 'InChi',
        'iupac': 'IUPAC',
        'preferred-cid': 'Preferred',
        'pubchem-sid': 'SID',
        'sid': 'SID',
        'smiles': 'SMILES',
        'synonym': 'Synonym-unfiltered',
        'cid': 'CID',
        'pubchem-cid': 'CID',
    }

    _target = id_types.get(target, target)
    _source = id_types.get(source, source)

    if _source not in {'CID', 'SID'}:

        msg = (
            'The source identifier type must be either CID or SID, '
            'not `%s`.' % source
        )
        _log(msg)
        raise ValueError(msg)

    ftp_dir = (
        {
            'SID': 'Substance',
            'CID': 'Compound',
        }[_source]
    )

    url = urls.urls['pubchem']['ftp'] % (ftp_dir, _source, _target)
    c = curl.Curl(url, large = True, silent = False)

    db_path = os.path.join(
        cache.get_cachedir(),
        'pubchem_%s_%s.sqlite' % (_source, _target)
    )

    with lookup.ManyToMany(db_path) as result:

        result.populate(fileobj = c._gzfile_mode_r)


def pubchem_name_cids(name: str) -> set[str]:
    """
    PubChem CIDs for a compound name via the PUG REST API.

    Each name–URL is cached on disk by pypath's curl infrastructure, so
    repeated calls across sessions incur no HTTP requests.  Unknown names
    (PubChem 404) and network errors return an empty set; the 404 response
    is also cached so the lookup is not retried on every run.

    Args:
        name: Compound name or synonym (e.g. ``'ATP'``, ``'NAD+'``).

    Returns:
        Set of PubChem CID strings matching *name*.  Empty set if the name
        is not found or the request fails.
    """

    url = urls.urls['pubchem']['name_cids'] % urllib.parse.quote(name)
    c = curl.Curl(url, silent = True, large = False)

    if not c.result:
        return set()

    try:
        data = json.loads(c.result)
    except (json.JSONDecodeError, TypeError):
        return set()

    # PubChem returns {"Fault": {...}} with HTTP 4xx for unknown names.
    if 'Fault' in data:
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
