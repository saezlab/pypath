#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Olga Ivanova
#           Sebastian Lobentanzer
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import collections

import pypath.share.session as session
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy
import pypath.inputs.rdata as rdata

_logger = session.Logger(name = 'progeny_input')
_log = _logger._log


def progeny_raw(organism = 9606):
    """
    Args:
        organism (int,str): Name or NCBI Taxonomy ID of the organism. Human
            and mouse are supported.
    """

    _organism = taxonomy.ensure_common_name(organism)

    if _organism not in ('Human', 'Mouse'):

        msg = (
            'Wrong organism: `%s`; '
            'only human and mouse are available.' % organism
        )
        _log(msg)
        raise ValueError(msg)

    _organism = _organism.lower()

    url = urls.urls['progeny']['url'] % _organism
    c = curl.Curl(url, large = True, silent = False)

    rdata_path = c.fileobj.name
    c.fileobj.close()

    rdata_parsed = rdata.rdata.parser.parse_file(rdata_path)
    rdata_converted = rdata.rdata.conversion.convert(rdata_parsed)

    key = 'model_%s_full' % _organism

    return rdata_converted[key]
