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

from typing import Iterable
import json
import re

import pypath.share.progress as progress
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.session as session
import pypath.share.common as common

"""
Clients to the NCBI Entrez Utils (E-utils) API.
"""

_logger = session.Logger(name = 'eutils_input')
_REPREFIX = re.compile(r'[A-z]*(\d+)')


def esummary(
        ids: str | Iterable[str],
        db: str,
        cache_small: bool | int = 10,
    ) -> list | dict:
    """
    Record metainformation from NCBI databases.

    A client to the eSummary API endpoint:
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESummary

    Args:
        ids:
            One or more record identifiers for the database of interst (`db`).
            The API understands the numeric IDs, letter prefixes if present
            will be stripped.
        db:
            The database, such as "pubmed" or "geoprofiles".
        cache_small:
            Small requests querying less than 10 IDs by default are not cached,
            except if this parameter is True or is set to a lower number.
    """

    ids = list(map(str, common.to_list(ids)))
    ids = sorted(_REPREFIX.sub('\\1', i) for i in ids)
    url = urls.urls['eutils']['esummary']
    cache = len(ids) < cache_small
    result = {}
    _logger._log(
        f'Retrieving data from NCBI E-utils eSummary: querying {len(ids)} '
        f'identifiers ({",".join(ids[:5])}...).'
    )
    prg = progress.Progress(
        len(ids) / 100 + 1,
        'Retrieving data from NCBI E-utils',
        1,
        percent = False,
    )

    for chunk in common.paginate(ids, 100):

        prg.step()
        post = {
            'id': ','.join(chunk),
            'retmode': 'json',
            'db': db,
        }

        for i in range(3):

            try:
                c = curl.Curl(
                    url,
                    silent = False,
                    cache = cache,
                    post = post,
                    override_post = True,
                )
                res = json.loads(c.result)
                if 'error' in res:
                    _logger._log(f'Error from NCBI E-utils: `{res["error"]}`.')
                result.update(res['result'].items())
                break

            except ValueError:
                _logger._log('Failed to process JSON from NCBI E-utils:')
                _logger._log_tracebakc()

    prg.terminate()
    _logger._log(
        'Finished retrieving data from NCBI E-utils eSummary, '
        f'result contains {len(result)} items.'
    )

    return result
