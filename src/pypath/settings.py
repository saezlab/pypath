#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module.
#  Settings for PyPath
#
#  Copyright
#  2014-2019
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

#TODO move to yaml file

from future.utils import iteritems

import os
import copy
import collections

import pypath.common as common


_defaults = {
    # name of the module
    'module_name': 'pypath',
    # The absolute root directory.
    # This should not be necessary, why is it here?
    'path_root': '/',
    # The basedir for every files and directories in the followings.
    'basedir': os.getcwd(),
    
    'progressbars': False,
    # verbosity for messages printed to console
    'console_verbosity': -1,
    # verbosity for messages written to log
    'log_verbosity': 0,
    # log flush time interval in seconds
    'log_flush_interval': 2,
    # check for expired mapping tables and delete them
    # (period in seconds)
    'mapper_cleanup_interval': 60,
    # If None will be the same as ``basedir``.
    'data_basedir': None,
    'acsn_names': 'acsn_names.gmt',
    'alzpw_ppi': 'alzpw-ppi.csv',
    'goose_annot_sql': 'goose_annotations.sql',
    'webpage_main': 'main.html',
    'nrf2ome': 'nrf2ome.csv',
    'ppoint': 'phosphopoint.csv',
    'slk3_nodes': 'signalink3_nodes.tsv',
    'acsn': 'acsn_ppi.txt',
    'arn': 'arn_curated.csv',
    'goose_ancest_sql': 'goose_ancestors.sql',
    'goose_terms_sql': 'goose_terms.sql',
    'lmpid': 'LMPID_DATA_pubmed_ref.xml',
    'nci_pid': 'nci-pid-strict.csv',
    'old_dbptm': 'old_dbptm.tab',
    'slk3_edges': 'signalink3_edges.tsv',
    'slk01human': 'slk01human.csv',
    'cachedir': None,
    'pubmed_cache': 'pubmed.pickle',
    'mapping_use_cache': True,
    'use_intermediate_cache': True,
    'default_organism': 9606,
    'default_name_types': {
        'protein': 'uniprot',
        'mirna': 'mirbase',
        'drug': 'chembl',
        'lncrna': 'lncrna-genesymbol',
    },
    'trip_preprocessed': 'trip_preprocessed.pickle',
    'deathdomain': 'deathdomain.tsv',
    'hpmr_preprocessed': 'hpmr_preprocessed.pickle',
    'network_expand_complexes': False,
    'network_keep_original_names': True,
    'network_pickle_cache': True,
    'go_pickle_cache': True,
    'go_pickle_cache_fname': 'goa__%u.pickle',
    'network_extra_directions': {
        'Wang',
        'KEGG',
        'STRING',
        'ACSN',
        'PhosphoSite',
        'PhosphoPoint',
        'CancerCellMap',
        'PhosphoSite_dir',
        'PhosphoSite_noref',
        'PhosphoNetworks',
        'MIMP',
        'HPRD-phos',
    },
    'keep_noref': False,
}

in_datadir = {
    'acsn_names',
    'alzpw_ppi',
    'goose_annot_sql',
    'webpage_main',
    'nrf2ome',
    'ppoint',
    'slk3_nodes',
    'acsn',
    'arn',
    'goose_ancest_sql',
    'goose_terms_sql',
    'lmpid',
    'nci_pid',
    'old_dbptm',
    'slk3_edges',
    'slk01human',
    'deathdomain',
}


in_cachedir = {
    'pubmed_cache',
    'trip_preprocessed',
    'hpmr_preprocessed',
}


def reset_all():
    
    settings = collections.namedtuple('Settings', list(_defaults.keys()))
    
    for k in _defaults.keys():
        
        val = getattr(defaults, k)
        
        if k in in_datadir:
            val = os.path.join(common.ROOT, 'data', val)
        
        setattr(settings, k, val)
    
    if settings.cachedir is None:
        
        settings.cachedir = os.path.join(
            os.path.expanduser('~'),
            '.pypath',
            'cache',
        )
    
    for k in in_cachedir:
        
        setattr(settings, k, os.path.join(settings.cachedir, _defaults[k]))
    
    globals()['settings'] = settings


def setup(**kwargs):
    
    for param, value in iteritems(kwargs):
        
        setattr(settings, param, value)


def get(param):
    
    if hasattr(settings, param):
        
        return getattr(settings, param)


def get_default(param):
    
    if hasattr(defaults, param):
        
        return getattr(defaults, param)


def reset(param):
    
    setup(param, get_default(param))


defaults = common._const()

for k, v in iteritems(_defaults):
    
    setattr(defaults, k, v)

reset_all()
