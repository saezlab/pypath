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
import pickle
import hashlib
import collections

import pypath.resources.urls as urls
import pypath.share.cache as cache
import pypath.share.session as session

_logger = session.Logger(name = 'pathophenodb_input')
_log = _logger._log


DiseasePathogen = collections.namedtuple(
    'DiseasePathogen',
    (
        'disease_id',
        'disease',
        'pathogen_taxid',
        'pathogen',
        'evidence'
    ),
)


def disease_pathogen_interactions():
    """
    Retrieves disease pathogen relationships from PathoPhenoDb.

    Returns:
        Disease-pathogen relationships as a list of tuples.
    """

    query = """#EX3:List all diseases which caused by pathogens
    PREFIX SIO: <http://semanticscience.org/resource/SIO_>
    PREFIX RO: <http://purl.obolibrary.org/obo/RO_>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    SELECT distinct ?Disease_ID ?Disease ?Pathogen_ID ?Pathogen ?evidence_Code
    FROM <http://patho.phenomebrowser.net>
    WHERE
    {
        ?Disease_ID SIO:000255 ?o .
        ?o RO:0002558 ?o1 .
        ?o RO:0002556  ?Pathogen_ID .
        ?Disease_ID rdfs:label ?Disease .
        ?Pathogen_ID rdfs:label ?Pathogen .
        ?o1 rdfs:label ?evidence_Code .
    }
    """

    url = urls.urls['pathophenodb']['url']
    urlmd5 = hashlib.md5(f'{url}{query}'.encode()).hexdigest()
    cache_path = os.path.join(cache.get_cachedir(), urlmd5)

    if os.path.exists(cache_path):

        try:

            _log(f'Loading from cache: `{cache_path}`.')

            with open(cache_path, 'rb') as fp:

                return pickle.load(fp)

        except:

            _log(
                f'Failed to load from `{cache_path}`, '
                'falling back to download.'
            )

    try:

        from SPARQLWrapper import SPARQLWrapper, JSON

    except ModuleNotFoundError:

        _logger._console(
            'No module `SPARQLWrapper` is available. '
            'Please install it to access PathoPhenoDB: '
            'pip install sparqlwrapper'
        )
        _log('Returning empty result!')
        return []

    sparql = SPARQLWrapper(url)
    sparql.setReturnFormat(JSON)
    sparql.setQuery(query)
    response = sparql.queryAndConvert()

    result = set()

    for r in response['results']['bindings']:

        pair = DiseasePathogen(
            disease_id = (
                r['Disease_ID']['value'].split('/')[-1].replace('_',':')
            ),
            disease = r['Disease']['value'],
            pathogen_taxid = (
                r['Pathogen_ID']['value'].split('/')[-1].split('_')[1]
            ),
            pathogen = r['Pathogen']['value'],
            evidence = r['evidence_Code']['value'],
        )

        result.add(pair)

    result = list(result)

    with open(cache_path, 'wb') as fp:

        _log(f'Saving to cache: `{cache_path}`.')

        pickle.dump(
            obj = result,
            file = fp,
        )

    return result
