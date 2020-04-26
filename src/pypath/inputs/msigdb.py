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

import os
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.settings as settings
import pypath.share.session as session
import pypath.utils.mapping as mapping

_logger = session.Logger(name = 'inputs.msigdb')
_log = _logger._log


def msigdb_download(
        registered_email = None,
        collection = 'msigdb',
        id_type = 'symbols',
        force_download = False,
    ):
    """
    Downloads and preprocesses a collection of gmt format gene sets from
    MSigDB. Returns dict of sets with gene set names as keys and molecular
    identifiers as values.

    :arg str,NoneType registered_email:
        An email address registered at MSigDB. If `None` the `msigdb_email`
        from ``pypath.settings`` will be used.
    :arg str collection:
        The name of the gene set collection. For available collections (e.g.
        `h.all` or `c2.cpg`) refer to the MSigDB website:
        http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
        The default value `msigdb` contains all the genesets however you
        won't be able to distinguish which geneset comes from which
        collection. For this you need to download the collections one by one.
    :arg str id_type:
        MSigDB provides Gene Symbols (`symbols`) and Entrez Gene IDs
        (`entrez`).
    :arg bool force_download:
        Download even if cache content is available.
    """

    registered_email = registered_email or settings.get('msigdb_email')

    if not registered_email:
        _log(
            'To download MSigDB you must provide an email address '
            'you have previously registered at '
            '`http://software.broadinstitute.org/gsea/register.jsp`. '
            'Could not proceed, returning empty dict.'
        )

        return {}

    url = urls.urls['msigdb']['url'] % (
        collection,
        id_type,
    )

    req_headers_1 = []

    c_nocall = curl.Curl(
        url,
        call = False,
        process = False,
        bypass_url_encoding = True,
    )

    if (
        not os.path.exists(c_nocall.cache_file_name) or
        os.path.getsize(c_nocall.cache_file_name) == 0 or
        force_download
    ):
        c_login_1 = curl.Curl(
            urls.urls['msigdb']['login1'],
            cache = False,
            write_cache = False,
            follow = False,
            large = False,
            silent = True,
        )

        jsessionid = ''

        if hasattr(c_login_1, 'resp_headers'):
            for hdr in c_login_1.resp_headers:
                if hdr.lower().startswith(b'set-cookie'):
                    jsessionid = hdr.split(b':')[1].split(b';')[0].strip()
                    jsessionid = jsessionid.decode('ascii')
                    _log('msigdb cookie obtained: `%s`.' % jsessionid)

                    break

        if not jsessionid:
            _log('msigdb: could not get cookie, returning empty list.')

            return {}

        req_headers = ['Cookie: %s' % jsessionid]

        c_login_2 = curl.Curl(
            urls.urls['msigdb']['login2'],
            cache = False,
            write_cache = False,
            large = False,
            silent = True,
            req_headers = req_headers,
            post = {
                'j_username': registered_email,
                'j_password': 'password',
            },
            follow = False,
            empty_attempt_again = False,
        )

        jsessionid_1 = ''

        if hasattr(c_login_2, 'resp_headers'):
            for hdr in c_login_2.resp_headers:
                if hdr.lower().startswith(b'set-cookie'):

                    jsessionid_1 = hdr.split(b':')[1].split(b';')[0].strip()
                    jsessionid_1 = jsessionid_1.decode('ascii')

            _log(
                'msigdb: logged in with email `%s`, '
                'new cookie obtained: `%s`.' % (
                    registered_email,
                    jsessionid_1
                )
            )

        if not jsessionid_1:
            _log(
                'msigdb: could not log in with email `%s`, '
                'returning empty dict.' % registered_email
            )

            return {}

        req_headers_1 = ['Cookie: %s' % jsessionid_1]

    c = curl.Curl(
        url,
        req_headers = req_headers_1,
        silent = False,
        large = True,
        bypass_url_encoding = True,
        cache = not force_download,
    )

    result = {}

    for gset in c.result:
        gset = gset.strip().split('\t')

        result[gset[0]] = set(gset[2:])

    return result


def msigdb_download_collections(
        registered_email = None,
        only_collections = None,
        exclude = ('c5',),
        id_type = 'symbols',
    ):
    """
    Downloads all or some MSigDB gene set collections.
    Returns a dict of dicts where upper level keys are collections while
    lower level keys are geneset names and values are molecular identifiers.

    :arg str,NoneType registered_email:
        An email address registered at MSigDB. If `None` the `msigdb_email`
        from ``pypath.settings`` will be used.
    :arg set,NoneType only_collections:
        Limit the annotations only to these collections. For available
        collections e.g. ``{'h.all', 'c2cgp'}`` refer to the MSigDB webpage:
        http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
    :arg tuple exclude:
        Exclude the collections having their name starting with any of the
        strings in this tuple. By default `c5` (Gene Ontology) is excluded.
    """

    MsigdbAnnotation = collections.namedtuple(
        'MsigdbAnnotation',
        [
            'collection',
            'geneset',
        ],
    )

    all_collections = {
        'hallmark': 'h.all',
        'positional': 'c1.all',
        'chemical_and_genetic_perturbations': 'c2.cgp',
        'biocarta_pathways': 'c2.cp.biocarta',
        'kegg_pathways': 'c2.cp.kegg',
        'pid_pathways': 'c2.cp.pid',
        'reactome_pathways': 'c2.cp.reactome',
        'mirna_targets': 'c3.mir',
        'tf_targets': 'c3.tft',
        'cancer_gene_neighborhoods': 'c4.cgn',
        'cancer_modules': 'c4.cm',
        'go_biological_process': 'c5.bp',
        'go_molecular_function': 'c5.mf',
        'go_cellular_component': 'c5.cc',
        'oncogenic_signatures': 'c6.all',
        'immunologic_signatures': 'c7.all',
    }

    collection_data = {}

    for collection, label in iteritems(all_collections):
        if (
            (
                only_collections and
                label not in only_collections
            ) or
            any(label.startswith(ex) for ex in exclude)
        ):

            continue

        _log(
            'MSigDB: downloading collection `%s` (%s).' % (collection, label)
        )

        collection_data[(collection, label)] = (
            msigdb_download(
                registered_email = registered_email,
                collection = label,
                id_type = id_type,
            )
        )

    return collection_data


def msigdb_annotations(
        registered_email = None,
        only_collections = None,
        exclude = ('c5',),
    ):
    """
    Downloads all or some MSigDB gene set collections and processes them
    to an annotation type dictionary.

    :arg str,NoneType registered_email:
        An email address registered at MSigDB. If `None` the `msigdb_email`
        from ``pypath.settings`` will be used.
    :arg set,NoneType only_collections:
        Limit the annotations only to these collections. For available
        collections e.g. ``{'h.all', 'c2cgp'}`` refer to the MSigDB webpage:
        http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
    :arg tuple exclude:
        Exclude the collections having their name starting with any of the
        strings in this tuple. By default `c5` (Gene Ontology) is excluded.

    """

    MsigdbAnnotation = collections.namedtuple(
        'MsigdbAnnotation',
        [
            'collection',
            'geneset',
        ],
    )


    annotations = collections.defaultdict(set)

    collection_data = msigdb_download_collections(
        registered_email = registered_email,
        only_collections = only_collections,
        exclude = exclude,
    )

    for (collection, label), genesets in iteritems(collection_data):

        for geneset, genesymbols in iteritems(genesets):

            this_annot = MsigdbAnnotation(
                collection = collection,
                geneset = geneset,
            )

            for uniprot in mapping.map_names(
                genesymbols,
                'genesymbol',
                'uniprot',
            ):
                annotations[uniprot].add(this_annot)

    return dict(annotations)
