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

import re
import time
import datetime

import timeloop

import pypath.share.session as session_mod
import pypath.share.settings as settings
import pypath.inputs.uniprot as uniprot

_logger = session_mod.Logger(name = 'uniprot_db')

db = {}
_cleanup_period = settings.get('mapper_cleanup_interval')
_lifetime = 300
_last_used = {}


def all_uniprots(organism = 9606, swissprot = None):

    return get_db(organism = organism, swissprot = swissprot)


def all_swissprots(organism = 9606):

    return get_db(organism = organism, swissprot = True)


def all_trembls(organism = 9606):

    return get_db(organism = organism, swissprot = False)


def init_db(organism = 9606, swissprot = None):

    _swissprot = uniprot._swissprot_param(swissprot)
    _logger._log(
        'Loading list of all UniProt IDs for '
        'organism `%u` (only SwissProt: %s).' % (
            organism,
            str(swissprot),
        )
    )

    key = (organism, _swissprot)

    globals()['db'][key] = uniprot._all_uniprots(
        organism = organism,
        swissprot = swissprot,
    )
    globals()['_last_used'][key] = time.time()


def get_db(organism = 9606, swissprot = None):

    _swissprot = uniprot._swissprot_param(swissprot)
    key = (organism, _swissprot)

    if key not in globals()['db']:

        init_db(organism = organism, swissprot = swissprot)

    globals()['_last_used'][key] = time.time()

    return globals()['db'][key]


def is_uniprot(name, organism = 9606, swissprot = None):
    """
    Tells if ``name`` is a UniProt ID of ``organism``.
    If ``swissprot`` is None then both SwissProt and TrEMBL IDs will be
    considered.
    """

    return name in get_db(organism = organism, swissprot = swissprot)


def is_swissprot(name, organism = 9606):
    """
    Tells if ``name`` is a SwissProt ID of ``organism``.
    For TrEMBL IDs returns False.
    """

    return is_uniprot(name, organism = organism, swissprot = True)


def is_trembl(name, organism = 9606):
    """
    Tells if ``name`` is a TrEMBL ID of ``organism``.
    For SwissProt IDs returns False.
    """

    return is_uniprot(name, organism = organism, swissprot = False)


_cleanup_timeloop = timeloop.Timeloop()
_cleanup_timeloop.logger.setLevel(9999)

@_cleanup_timeloop.job(
    interval = datetime.timedelta(
        seconds = _cleanup_period
    )
)
def _cleanup():

    keys = list(globals()['db'].keys())

    for key in keys:

        if time.time() - globals()['_last_used'][key] > _lifetime:

            _remove(key)

_cleanup_timeloop.start(block = False)


def _remove(key):

    if key in globals()['db']:

        _logger._log(
            'Removing UniProt ID list for '
            'organism `%u` (SwissProt: %s)' % (
                key[0],
                str(key[1]),
            )
        )
        del globals()['db'][key]

    if key in globals()['_last_used']:

        del globals()['_last_used'][key]
