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

import os
import csv
import collections
import base64
import json

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session_mod
import pypath.share.settings as settings
import pypath.utils.mapping as mapping
import pypath.inputs.credentials as credentials

_logger = session_mod.Logger(name = 'cosmic_input')
_log = _logger._log


def cancer_gene_census_annotations(
        user = None,
        passwd = None,
        credentials_fname = 'cosmic_credentials',
    ):
    """
    Retrieves a list of cancer driver genes (Cancer Gene Census) from
    the Sanger COSMIC (Catalogue of Somatic Mutations in Cancer) database.
    Returns dict of annotations.
    """

    try:

        cosmic_cred = credentials.credentials(
            user = user,
            passwd = passwd,
            resource = 'COSMIC',
            from_file = credentials_fname,
        )

    except RuntimeError:

        _log(
            'No credentials available for the COSMIC website. '
            'Either set the `cosmic_credentials` key in the `settings` '
            'module (e.g. `{\'user\': \'myuser\', '
            '\'passwd\': \'mypassword\'}`), or pass them directly to the '
            '`pypath.inputs.cosmic.cancer_gene_census_annotations` '
            'method.'
        )

        return {}

    CancerGeneCensusAnnotation = collections.namedtuple(
        'CancerGeneCensusAnnotation',
        (
            'tier',
            'hallmark',
            'somatic',
            'germline',
            'tumour_types_somatic',
            'tumour_types_germline',
            'cancer_syndrome',
            'tissue_type',
            'genetics',
            'role',
            'mutation_type',
        ),
    )


    def multi_field(content):

        return (
            tuple(sorted(i.strip() for i in content.split(',')))
                if content.strip() else
            ()
        )


    url = urls.urls['cgc']['url_new']

    auth_str = base64.b64encode(
        ('%s:%s\n' % (cosmic_cred['user'], cosmic_cred['passwd'])).encode()
    )

    req_hdrs = ['Authorization: Basic %s' % auth_str.decode()]

    c = curl.Curl(
        url,
        large = False,
        silent = False,
        req_headers = req_hdrs,
        cache = False,
    )

    access_url = json.loads(c.result)

    if 'url' not in access_url:

        _log(
            'Could not retrieve COSMIC access URL. '
            'Most likely the authentication failed. '
            'The reply was: `%s`' % c.result
        )

        return None

    c = curl.Curl(
        access_url['url'],
        large = True,
        silent = False,
        bypass_url_encoding = True,
    )

    data = csv.DictReader(c.fileobj, delimiter = ',')
    result = collections.defaultdict(set)

    for rec in data:

        uniprots = mapping.map_name(
            rec['Gene Symbol'],
            'genesymbol',
            'uniprot',
        )

        for uniprot in uniprots:

            result[uniprot].add(
                CancerGeneCensusAnnotation(
                    tier = int(rec['Tier']),
                    hallmark = rec['Hallmark'].strip().lower() == 'yes',
                    somatic = rec['Somatic'].strip().lower() == 'yes',
                    germline = rec['Germline'].strip().lower() == 'yes',
                    tumour_types_somatic = (
                        multi_field(rec['Tumour Types(Somatic)'])
                    ),
                    tumour_types_germline = (
                        multi_field(rec['Tumour Types(Germline)'])
                    ),
                    cancer_syndrome = (
                        multi_field(rec['Cancer Syndrome'])
                    ),
                    tissue_type = (
                        multi_field(rec['Tissue Type'].replace(' ', ''))
                    ),
                    genetics = rec['Molecular Genetics'].strip() or None,
                    role = (
                        multi_field(rec['Role in Cancer'])
                    ),
                    mutation_type = (
                        multi_field(rec['Mutation Types'])
                    ),
                )
            )

    return dict(result)
