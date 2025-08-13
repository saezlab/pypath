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

from future.utils import iteritems

import os
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.settings as settings
import pypath.share.session as session
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy

_logger = session.Logger(name = 'msigdb_input')
_log = _logger._log


ALL_COLLECTIONS = {
    'hallmark': ('h.all', 'mh.all'),
    'positional': ('c1.all', 'm1.all'),
    'chemical_and_genetic_perturbations': ('c2.cgp', 'm2.cgp'),
    'biocarta_pathways': ('c2.cp.biocarta', 'm2.cp.biocarta'),
    'kegg_pathways': ('c2.cp.kegg_legacy', None),
    'kegg_medicus_pathways': ('c2.cp.kegg_medicus', None),
    'pid_pathways': ('c2.cp.pid', None),
    'reactome_pathways': ('c2.cp.reactome', 'm2.cp.reactome'),
    'wikipathways': ('c2.cp.wikipathways', 'm2.cp.wikipathways'),
    'mirna_targets_mirdb': ('c3.mir.mirdb', 'm3.mirdb'),
    'mirna_targets_legacy': ('c3.mir.mir_legacy', None),
    'tf_targets_gtrf': ('c3.tft.gtrd', 'm3.gtrd'),
    'tf_targets_legacy': ('c3.tft.tft_legacy', None),
    'cancer_cell_atlas': ('c4.3ca', None),
    'cancer_gene_neighborhoods': ('c4.cgn', None),
    'cancer_modules': ('c4.cm', None),
    'go_biological_process': ('c5.go.bp', 'm5.go.bp'),
    'go_molecular_function': ('c5.go.mf', 'm5.go.mf'),
    'go_cellular_component': ('c5.go.cc', 'm5.go.cc'),
    'human_phenotype_ontology': ('c5.hpo', None),
    'mouse_phenotype_ontology': (None, 'm5.mpt'),
    'oncogenic_signatures': ('c6.all', None),
    'immunesigdb': ('c7.immunesigdb', None),
    'vaccine_response': ('c7.vax', None),
    'cell_type_signatures': ('c8.all', 'm8.all'),
}


def _msigdb_login(registered_email: str | None = None) -> dict[str, str]:
    """
    Logs in to MSigDB and returns a dict of cookies.
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

    c_login_1 = curl.Curl(
        urls.urls['msigdb']['login1'],
        cache = False,
        write_cache = False,
        process = False,
        large = True,
        silent = True,
        empty_attempt_again = False,
    )

    cookies = c_login_1.cookies()

    if not cookies:

        _log('msigdb: could not obtain cookie, returning empty list.')

        return {}

    c_login_2 = curl.Curl(
        urls.urls['msigdb']['login2'],
        cache = False,
        write_cache = False,
        large = False,
        silent = True,
        cookies = cookies,
        post = {
            'username': registered_email,
            'password': 'password',
        },
        process = False,
        empty_attempt_again = False,
    )

    if hasattr(c_login_2, 'resp_headers'):

        cookies = c_login_2.cookies()

        _log(
            'msigdb: logged in with email `%s`, '
            'new cookie obtained: `JSESSIONID=%s`.' % (
                registered_email,
                cookies['JSESSIONID']
            )
        )

    globals()['COOKIES'] = cookies

    return cookies


def _is_html(path: str) -> bool:

    with open(path) as f:

        return f.read(10).strip()[:2] == '<!'


def msigdb_download(
        registered_email = None,
        collection = 'msigdb',
        id_type = 'symbols',
        force_download = False,
        organism = 'human',
        version = '2025.1',
        cookies = None,
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


    organisms = {9606: 'Hs', 10090: 'Mm'}
    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)
    msigdb_org = organisms.get(ncbi_tax_id, None)

    if not msigdb_org:

        _log(f'Could not recognize organism: `{organism}`.')

        return {}

    version = version or settings.get('msigdb_version')

    #http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2022.1.Mm/mh.all.v2022.1.Mm.symbols.gmt
    #http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2022.1.Hs/h.all.v2022.1.Hs.symbols.gmt


    url = urls.urls['msigdb']['url'] % (
        version,
        msigdb_org,
        collection,
        version,
        msigdb_org,
        id_type,
    )

    c_nocall = curl.Curl(
        url,
        call = False,
        process = False,
        bypass_url_encoding = True,
    )

    if (
        (
            not os.path.exists(c_nocall.cache_file_name) or
            os.path.getsize(c_nocall.cache_file_name) == 0 or
            force_download or
            _is_html(c_nocall.cache_file_name)
        ) and
        not cookies
    ):

        cookies = globals().get('COOKIES')

        if not cookies:

            _log('msigdb: no cookies available, logging in...')
            cookies = _msigdb_login(registered_email)

            if not cookies:

                _log('msigdb: could not log in, returning empty list.')
                return {}

        else:

            _log('msigdb: using cookies from previous login.')

    c = curl.Curl(
        url,
        silent = False,
        large = True,
        cookies = cookies,
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
        exclude = ('c5', 'm5'),
        id_type = 'symbols',
        organism = 'human',
        version = '2025.1',
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
        collections e.g. ``{'h.all', 'c2.cgp'}`` refer to the MSigDB webpage:
        http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
    :arg tuple exclude:
        Exclude the collections having their name starting with any of the
        strings in this tuple. By default `c5` and `m5` (Gene Ontology and
        Human/Mouse Phenotype Ontology) is excluded.
    """

    collection_data = {}

    organisms = {9606: 0, 10090: 1}
    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)
    idx = organisms.get(ncbi_tax_id, None)

    for collection, labels in iteritems(ALL_COLLECTIONS):

        label = labels[idx]

        if (
            not label or
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
                organism = organism,
                version = version,
            )
        )

    return collection_data


def msigdb_annotations(
        registered_email = None,
        only_collections = None,
        exclude = ('c5', 'm5'),
        organism = 'human',
        version = '2025.1',
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
        organism = organism,
        version = version,
    )

    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)

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
                ncbi_tax_id = ncbi_tax_id,
            ):

                annotations[uniprot].add(this_annot)

    return dict(annotations)
