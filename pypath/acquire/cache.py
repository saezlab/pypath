#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
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
This file contains a custom cache system to collect data from primary
sources. Here are the (preliminary) requirements:


- For each download the URL, GET and POST parameters, the time of access
  and path to the local file should be stored in an inventory

- One item might have multiple versions

- The users should be able to discard or disable versions based on
  custom criteria, e.g. "I want fresh download for everything from
  UniProt, but use the cached contents for anything else" or "Use only
  cache content that is not older than 6 months" or "Delete all cache
  contents which is older than 1 year"

- Sebastian: I think that the "regions" of beaker can be very useful, in
  terms of determining the life time of a particular input generally.
  This way, we don't require explicit decisions from the user (eg
  "discard everything older than 6 months" may be overkill for resources
  that do not change often). Regions could be defined in settings YAML
  as a "sensible default", and overwritten at run-time if the user
  desires a particular manual change.

- This inventory should be stored in an SQlite database, so we don't
  need to deal with concurrency and data integrity

- OmnipathR has already this kind of facility, although it works by a
  JSON file, we can copy many ideas from there:
  https://github.com/saezlab/OmnipathR/blob/master/R/cache.R

- The logic of the implementation looks like this:

    - The requesting procedure provides an URL with parameters

    - The cache manager checks if it already exists in the cache

    - If exists, does it have any version which meets all the current
      conditions for cache use?

    - If the item or version does not exist, create a new one, and write
      in the database that "download started"

    - Start the download

    - Once the download finished, write into the database the timestamp
      and local path

    - Return the local path or open file pointer to the requesting
      procedure
"""

import functools

class CacheManager():
    def region(func):
        @functools.wraps(func)
        def wrapper_decorator(*args, **kwargs):
            # check if data (url + params) already in db
            # if so, get date last modified
            # compare to (individual) shelf life of the input
            # if shelf life exceeded, re-download
            #   lock-file, "download started at" ... ?
            #   once download finished, update metadata
            # return stored data
            
            value = func(*args, **kwargs)
            return value
        return wrapper_decorator