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
import urllib

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session_mod
import pypath.share.settings as settings
import pypath.utils.mapping as mapping
import pypath.inputs.credentials as credentials

_logger = session_mod.Logger(name = 'cosmic_input')
_log = _logger._log


def _is_html(path: str) -> bool:

    with open(path) as f:

        return f.read(10).strip()[:2] == '<!'


def _cosmic_cgc_download(
        user = None,
        passwd = None,
        credentials_fname = 'cosmic_credentials',
        force_download = False,
    ):

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

        return


    url = urls.urls['cgc']['cgc_tsv']
    cookies = {}

    c_nocall = curl.Curl(
        url,
        call = False,
        process = False,
        bypass_url_encoding = True,
    )

    if (
        not os.path.exists(c_nocall.cache_file_name) or
        os.path.getsize(c_nocall.cache_file_name) == 0 or
        force_download or
        _is_html(c_nocall.cache_file_name)
    ):

        req_headers = [settings.get('user_agent')]
        login_url = urls.urls['cgc']['login']
        pgsmith_cookie = r'%7B%22z%22%3A%22n%22%2C%22a%22%3A%22e%22%7D'
        login_payload = {
            'email': cosmic_cred['user'],
            'pass': cosmic_cred['passwd'],
            'r_url': '',
            'd': '0',
        }
        login_payload = urllib.parse.urlencode(login_payload)

        c_login_1 = curl.Curl(
            login_url,
            cache = False,
            req_headers = req_headers,
            write_cache = False,
            process = False,
            large = False,
            silent = True,
            empty_attempt_again = False,
        )

        cookies = c_login_1.cookies()
        cookies['Pagesmith'] = pgsmith_cookie

        c_login_2 = curl.Curl(
            login_url,
            post = login_payload,
            cookies = cookies,
            req_headers = req_headers,
            cache = False,
            large = False,
            silent = True,
            empty_attempt_again = False,
        )

        cookies = c_login_2.cookies()
        cookies['redirect_to'] = (
            '%2Fcensus%2Fall%3Fhome%3Dy%26name%3Dall%26tier%3D%26'
            'export%3Djson%26sEcho%3D1'
        )
        cookies['Pagesmith'] = pgsmith_cookie
        cookies['DNT'] = '1'

    c = curl.Curl(
        url,
        large = True,
        silent = False,
        cookies = cookies,
        bypass_url_encoding = True,
    )

    return getattr(c, 'fileobj', None)


def _cosmic_cgc_rescued():

    url = urls.urls['cgc']['tsv_rescued']
    c = curl.Curl(url, large = True, silent = False)

    return getattr(c, 'fileobj', None)


def cancer_gene_census_raw(
        user = None,
        passwd = None,
        credentials_fname = 'cosmic_credentials',
        force_download = False,
    ):
    """
    Retrieves a list of cancer driver genes (Cancer Gene Census) from
    the Sanger COSMIC (Catalogue of Somatic Mutations in Cancer) database.
    Returns dict of annotations.
    """

    if True:

        fileobj = _cosmic_cgc_rescued()

    else:

        fileobj = _cosmic_cgc_download(
            user = user,
            passwd = passwd,
            credentials_fname = credentials_fname,
            force_download = force_download,
        )

    if fileobj is None:

        _log('Failed to download COSMIC Cancer Gene Census.')
        result = []

    else:

        result = csv.DictReader(fileobj, delimiter = '\t')

    yield from result


def cancer_gene_census_annotations(
        user = None,
        passwd = None,
        credentials_fname = 'cosmic_credentials',
        force_download = False,
    ):
    """
    Retrieves a list of cancer driver genes (Cancer Gene Census) from
    the Sanger COSMIC (Catalogue of Somatic Mutations in Cancer) database.
    Returns dict of annotations.
    """

    CancerGeneCensusAnnotation = collections.namedtuple(
        'CancerGeneCensusAnnotation',
        (
            'tier',
            'hallmark',
            'somatic',
            'germline',
            'chr_band',
            'tumour_types_somatic',
            'tumour_types_germline',
            'cancer_syndrome',
            'tissue_type',
            'genetics',
            'role',
            'mutation_type',
            'translocation_partner',
        ),
    )


    def multi_field(content):

        return (
            tuple(sorted(i.strip() for i in content.split(',')))
                if content.strip() else
            ()
        )


    raw = cancer_gene_census_raw(
        user = user,
        passwd = passwd,
        credentials_fname = credentials_fname,
        force_download = force_download,
    )

    result = collections.defaultdict(set)

    for rec in raw:

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
                    chr_band = rec['Chr Band'],
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
                    translocation_partner = (
                        multi_field(rec['Translocation Partner']),
                    ),
                )
            )

    return dict(result)
