#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  (Planned for) centrally handling cache for all databases/resources.
#
#  Copyright
#  2014-2022
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

"""
File cache management using beaker
(https://beaker.readthedocs.io/en/latest/index.html).

Notes:
    "_dbm.error: cannot add item to database": sometimes
    deleting cache helps. The dbm.db file goes to 4.4 TB
    on one specific biomodel download and can't be repaired.
"""

from beaker.cache import CacheManager, cache_regions
from beaker.util import parse_cache_config_options

import pypath.share.settings as settings
import pypath.share.session as session

_logger = session.Logger(name = 'beaker')
_log = _logger._log

# def get_cache_regions():
#     cache_regions.update({
#         'short_term':{
#             'expire':60,
#             'type':'memory'
#         },
#         'long_term':{
#             'expire':'false',
#             'type':'dbm',
#         }
#     })

#     return cache_regions

def get_cache_manager():
    """
    Get beaker cache manager for two use cases ("regions"): short term 
    buffering and long term storage (the latter being the standard use
    case for pypath downloads). Long term storage does not expire and 
    can be manually purged.

    Returns:
        cache manager instance

    Todo:
        explicit purge functionality

        do we want a "daily" type (for updating lists such as all 
        biomodels)?
    """
    _log("Returning cache manager.")

    cache_opts = {
        'cache.data_dir': settings.get('cachedir') + '/data',
        'cache.lock_dir': settings.get('cachedir') + '/lock',
        'cache.regions': 'short_term, long_term',
        'cache.short_term.type': 'memory',
        'cache.short_term.expire': '60',
        'cache.long_term.type': 'dbm',
    }

    cache = CacheManager(**parse_cache_config_options(cache_opts))

    return cache