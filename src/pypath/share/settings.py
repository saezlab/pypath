#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module.
#  Settings for PyPath
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

#TODO move to yaml file

from future.utils import iteritems

import os
import copy
import collections

ROOT = os.path.join(
    *os.path.split(
        os.path.abspath(
            os.path.dirname(__file__)
        )
    )[:-1]
)


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
    'mapper_translate_deleted_uniprot': False,
    'mapper_keep_invalid_uniprot': False,
    'mapper_trembl_swissprot_by_genesymbol': True,
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
    # directory in datadir with licenses
    'license_dir': 'licenses',
    # password file for within company license-free redistribution
    'secrets_dir': None,
    'license_secret': 'license_secret',
    'server_default_license': 'academic',
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
    'default_label_types': {
        'protein': 'genesymbol',
        'mirna': 'mir-mat-name',
        'lncrna': 'lncrna-genesymbol',
    },
    'uniprot_uploadlists_chunk_size': 10000,
    'trip_preprocessed': 'trip_preprocessed.pickle',
    'deathdomain': 'deathdomain.tsv',
    'hpmr_preprocessed': 'hpmr_preprocessed.pickle',
    'network_expand_complexes': False,
    'network_allow_loops': False,
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
    'msigdb_email': 'omnipathdb@gmail.com',

    # the annotation classes should infer complex annotations
    # from protein annotations
    'annot_infer_complexes': True,
    # the resource name for annotation categories
    # combined from multiple original resources
    'annot_composite_database_name': 'OmniPath',

    # load small, specific categories from CellPhoneDB
    # in the intercell database
    'intercell_cellphonedb_categories': True,
    # same for Baccin2019 and some others
    'intercell_baccin_categories': True,
    'intercell_hpmr_categories': True,
    'intercell_surfaceome_categories': True,
    'intercell_gpcrdb_categories': True,
    'intercell_icellnet_categories': True,

    # parameters for pypath.omnipath
    'timestamp_format': '%Y%m%d',

    # tfregulons levels
    'tfregulons_levels': {'A', 'B', 'C', 'D'},

    # datasets
    'datasets': [
       'omnipath',
       'curated',
       'complex',
       'annotations',
       'intercell',
       'tf_target',
       'tf_mirna',
       'mirna_mrna',
       'lncrna_mrna',
       'enz_sub',
    ],

    'omnipath_mod': 'network',
    'curated_mod': 'network',
    'complex_mod': 'complex',
    'annotations_mod': 'annot',
    'intercell_mod': 'intercell',
    'enz_sub_mod': 'enz_sub',
    'tf_target_mod': 'network',
    'tf_mirna_mod': 'network',
    'mirna_mrna_mod': 'network',
    'lncrna_mrna_mod': 'network',

    'omnipath_args': {
        'use_omnipath': True,
        'kinase_substrate_extra': True,
        'ligand_receptor_extra': True,
        'pathway_extra': True,
    },

    # only for pypath.omnipath.app and pypath.core.network
    'dorothea_expand_levels': True,

    'dependencies': {
        'intercell': ('annotations',),
        'annotations': ('complex',),
    },

    'omnipath_pickle': 'network_omnipath.pickle',
    'curated_pickle': 'network_curated.pickle',
    'complex_pickle': 'complexes.pickle',
    'annotations_pickle': 'annotations.pickle',
    'intercell_pickle': 'intercell.pickle',
    'enz_sub_pickle': 'enz_sub_%u.pickle',
    'tf_target_pickle': 'tftarget.pickle',
    'tf_mirna_pickle': 'tfmirna.pickle',
    'mirna_mrna_pickle': 'mirna_mrna.pickle',
    'lncrna_mrna_pickle': 'lncrna_mrna.pickle',

    'pickle_dir': None,

    # directory for exported tables
    'tables_dir': 'omnipath_tables',

    # directory for figures
    'figures_dir': 'omnipath_figures',

    # directory for LaTeX
    'latex_dir': 'omnipath_latex',

    # include a timestamp in directory names
    'timestamp_dirs': True,

    # maximum lenght of the strings in UniProt info printed tables
    'uniprot_info_maxlen': 500,
    # at downloading UniProt datasheets the default very long timeouts,
    # what we use in the curl module especially because many of our downloads
    # are huge, are too long and better to start the next attempt sooner if
    # the first fails to respond
    'uniprot_datasheet_connect_timeout': 10,
    'uniprot_datasheet_timeout': 20,
    'genecards_datasheet_connect_timeout': 10,
    'genecards_datasheet_timeout': 20,

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
    'license_dir',
}


in_cachedir = {
    'pubmed_cache',
    'trip_preprocessed',
    'hpmr_preprocessed',
}


in_secrets_dir = {
    'license_secret',
}


class Settings(object):


    def __init__(self, **kwargs):

        self.__dict__.update(kwargs)


Defaults = collections.namedtuple(
    'Defaults',
    sorted(_defaults.keys()),
)


def reset_all():

    settings = Settings()

    for k in _defaults.keys():

        val = getattr(defaults, k)

        if k in in_datadir:

            val = os.path.join(ROOT, 'data', val)

        setattr(settings, k, val)

    # special directories with built in default at user level
    pypath_dirs = (
        ('cachedir', 'cache'),
        ('pickle_dir', 'pickles'),
        ('secrets_dir', 'secrets'),
    )

    for _key, _dir in pypath_dirs:

        if getattr(settings, _key) is None:

            setattr(
                settings,
                _key,
                os.path.join(
                    os.path.expanduser('~'),
                    '.pypath',
                    _dir,
                )
            )

    for k in in_cachedir:

        setattr(settings, k, os.path.join(settings.cachedir, _defaults[k]))

    for k in in_secrets_dir:

        setattr(settings, k, os.path.join(settings.secrets_dir, _defaults[k]))

    globals()['settings'] = settings


def setup(**kwargs):

    for param, value in iteritems(kwargs):

        setattr(settings, param, value)


def get(param, value = None):
    """
    Retrieves the current value of a parameter.

    :param str param:
        The key for the parameter.
    :param object,NoneType value:
        If this value is not None it will be returned instead of the settings
        value. It is useful if the parameter provided at the class or method
        level should override the one in settings.
    """

    if value is not None:

        return value

    if hasattr(settings, param):

        return getattr(settings, param)


def get_default(param):

    if hasattr(defaults, param):

        return getattr(defaults, param)


def reset(param):

    setup(**{param: get_default(param)})


defaults = Defaults(**_defaults)


reset_all()
