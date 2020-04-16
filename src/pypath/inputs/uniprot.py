#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from future.utils import iteritems

import re
import time
import datetime
import timeloop
timeloop.app.logging.disable(level = 9999)

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session_mod
import pypath.share.settings as settings

_logger = session_mod.Logger(name = 'uniprot_input')

db = {}
_cleanup_period = settings.get('mapper_cleanup_interval')
_lifetime = 300
_last_used = {}


_redatasheet = re.compile(r'([A-Z\s]{2})\s*([^\n\r]+)[\n\r]+')


def _all_uniprots(organism = 9606, swissprot = None):

    swissprot = 'yes' if swissprot == True else swissprot
    rev = '' if not swissprot else ' AND reviewed: %s' % swissprot
    url = urls.urls['uniprot_basic']['url']
    get = {
        'query': 'organism:%s%s' % (str(organism), rev),
        'format': 'tab',
        'columns': 'id',
    }
    c = curl.Curl(url, get = get, silent = False)
    data = c.result

    return [
        l.strip() for l in data.split('\n')[1:] if l.strip()
    ]


def all_uniprots(organism = 9606, swissprot = None):

    return get_db(organism = organism, swissprot = swissprot)


def init_db(organism = 9606, swissprot = None):

    _logger._log(
        'Loading list of all UniProt IDs for '
        'organism `%u` (only SwissProt: %s).' % (
            organism,
            str(swissprot == True),
        )
    )

    key = (organism, swissprot == True)

    globals()['db'][key] = _all_uniprots(
        organism = organism,
        swissprot = swissprot,
    )
    globals()['_last_used'][key] = time.time()


def get_db(organism = 9606, swissprot = None):

    key = (organism, swissprot == True)

    if key not in globals()['db']:

        init_db(organism = organism, swissprot = swissprot)

    globals()['_last_used'][key] = time.time()

    return globals()['db'][key]


def is_uniprot(name, organism = 9606, swissprot = None):
    """
    Tells if ``name`` is a UniProt ID of ``organism``.
    """

    return name in get_db(organism = organism, swissprot = swissprot)


_cleanup_timeloop = timeloop.Timeloop()

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
            'organism `%u` (only SwissProt: %s)' % (
                key[0],
                str(key[1]),
            )
        )
        del globals()['db'][key]

    if key in globals()['_last_used']:

        del globals()['_last_used'][key]


def protein_datasheet(identifier):

    url = urls.urls['uniprot_basic']['datasheet'] % identifier
    c = curl.Curl(url, silent = True, large = False)

    return _redatasheet.findall(c.result) if c.result else []


def get_uniprot_sec(organism = 9606):
    """
    Downloads and processes the mapping between secondary and
    primary UniProt IDs.

    Yields pairs of secondary and primary UniProt IDs.

    :param int organism:
        NCBI Taxonomy ID of the organism.
    """

    if organism is not None:
        proteome = all_uniprots(organism=organism)
        proteome = set(proteome)

    sec_pri = []
    url = urls.urls['uniprot_sec']['url']
    c = curl.Curl(url, silent = False, large = True, timeout = 2400)

    for line in filter(
        lambda line:
            len(line) == 2 and (organism is None or line[1] in proteome),
            map(
                lambda i:
                    i[1].split(),
                filter(
                    lambda i: i[0] >= 30,
                    enumerate(c.result)
                )
            )
        ):

        yield line


_uniprot_fields = {
    'function': 'comment(FUNCTION)',
    'activity_regulation': 'comment(ACTIVITY REGULATION)',
    'tissue_specificity': 'comment(TISSUE SPECIFICITY)',
    'developmental_stage': 'comment(DEVELOPMENTAL STAGE)',
    'induction': 'comment(INDUCTION)',
    'intramembrane': 'feature(INTRAMEMBRANE)',
    'signal_peptide': 'feature(SIGNAL)',
    'subcellular_location': 'comment(SUBCELLULAR LOCATION)',
    'transmembrane': 'feature(TRANSMEMBRANE)',
    'comment': 'comment(MISCELLANEOUS)',
    'topological_domain': 'feature(TOPOLOGICAL DOMAIN)'
}


def uniprot_data(field, organism = 9606, reviewed = True):
    """
    Retrieves a field from UniProt for all proteins of one organism, by
    default only the reviewed (SwissProt) proteins.
    For the available fields refer to the ``_uniprot_fields`` attribute of
    this module or the UniProt website.
    """
    
    rev = (
        ' AND reviewed: yes'
            if reviewed == True or reviewed == 'yes' else
        ' AND reviewed: no'
        if reviewed == False or reviewed == 'no' else
        ''
    )
    _field = _uniprot_fields[field] if field in _uniprot_fields else field
    url = urls.urls['uniprot_basic']['url']
    get = {
        'query': 'organism:%s%s' % (str(organism), rev),
        'format': 'tab',
        'columns': 'id,%s' % _field,
        'compress': 'yes',
    }
    
    c = curl.Curl(url, get = get, silent = False, large = True, compr = 'gz')
    _ = next(c.result)

    return dict(
        id_value
        for id_value in
        (
            line.strip('\n\r').split('\t')
            for line in c.result if line.strip('\n\r')
        )
        if id_value[1]
    )
