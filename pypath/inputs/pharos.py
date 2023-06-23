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

"""
Retrieve data from the NIH Pharos database.
"""

from __future__ import annotations

import json

from pypath.share.curl import Curl
from pypath.resources.urls import urls
import pypath.share.session as session

_logger = session.Logger(name = 'pharos_input')
_log = _logger._log


QUERY_TYPES = (
    'expression',
    'gtex',
    'orthologs',
    'ligands',
    'xrefs',
    'diseases',
)


PHAROS_QUERY = """
    query targetDetails(
        $chunk_size: Int!,
        $step: Int!,
        $getExpressions: Boolean!,
        $getGtex: Boolean!,
        $getOrthologs: Boolean!,
        $getLigands: Boolean!,
        $getXrefs: Boolean!,
        $getDiseases: Boolean!,
    ) {

        targets {

            targets(top: $chunk_size skip: $step) {

                name
                sym
                uniprot

                expressions(top: 10000) @include(if: $getExpressions) {

                    expid
                    type
                    tissue
                    value

                    uberon {
                        name
                        uid
                    }

                    pub {
                        pmid
                    }
                }

                gtex @include(if: $getGtex) {

                    tissue
                    tpm
                    log2foldchange

                    uberon {
                        name
                        uid
                    }
                }

                orthologs(top: 10000) @include(if: $getOrthologs) {
                    name
                    species
                    orid
                    dbid
                    geneid
                    source
                }

                ligands(top: 10000 isdrug: true) @include(if: $getLigands) {

                    ligid
                    name

                    synonyms {
                        name
                        value
                    }

                    activities(all: true) {
                        actid
                        type
                        moa
                        value

                        pubs {
                            pmid
                            __typename
                            }
                    }
                }

                xrefs(source: "Ensembl") @include(if: $getXrefs) {
                    name
                }

                diseases(top:10000) @include(if: $getDiseases) {

                    name
                    mondoID

                    dids {
                        id
                        dataSources
                        doName
                    }
                }
            }
        }
    }
    """


def pharos_general(
        query: str,
        variables: dict[str, bool] | None = None,
    ) -> dict:
    """
    Query the NIH Pharos database by GraphQL.

    Read about Pharos here: https://pharos.nih.gov/about

    Args:
        query:
            A GraphQL query.
        variables:
            Variables to retrieve. A dict of variable names and boolean values.

    Return:
        The JSON response parsed into a dict.
    """

    url = urls['pharos_api']['url']

    req_headers = {
        'Accept-Encoding': 'gzip, deflate, br',
        'Content-Type': 'application/json',
        'Connection': 'keep-alive',
        'DNT': '1',
        'Origin': 'https://pharos-api.ncats.io',
    }

    query_param = {'query': query}

    if variables:

        _log(
            'Querying Pharos, variables: '
             f'{", ".join(k for k, v in variables.items() if v)}'
        )
        query_param['variables'] = variables

    binary_data = json.dumps(query_param).encode('utf-8')

    c = Curl(
        url=url,
        req_headers=req_headers,
        binary_data=binary_data,
        compressed=True,
        compr='gz',
    )
    result = json.loads(c.result)

    result = result['data']

    return result


def pharos_targets(
        chunk_size: int = 100,
        expression: bool = False,
        gtex: bool = False,
        orthologs: bool = False,
        ligands: bool = False,
        xrefs: bool = False,
        diseases: bool = False,
    ) -> list:
    """
    Query the NIH Pharos database by GraphQL.

    The queried data is fetched in chunks, by default 100 records each. The
    complete data consists of thousands of chunks, the retrieval takes
    about half hour.

    Args:
        chunk_size:
            Records in one batch. Better stay 100 because higher numbers
            likely to cause timeout errors.

    Return:
        Records as a list of dicts.
    """

    variables = {
        'chunk_size': chunk_size,
        'step': 0,
        'getExpressions': expression,
        'getGtex': gtex,
        'getOrthologs': orthologs,
        'getLigands': ligands,
        'getXrefs': xrefs,
        'getDiseases': diseases,
    }
    result = []

    while True:

        _log(f'Pharos query, chunk #{variables["step"]}')
        response = pharos_general(PHAROS_QUERY, variables)
        response = response['targets']['targets']

        if not response:

            break

        result.extend(response)
        variables['step'] += chunk_size

    return result


def _create_query_functions():

    for qtype in QUERY_TYPES:

        args = {qtype: True}
        name = f'pharos_{qtype}'

        doc = f"""
            Retrieve `{qtype}` records from Pharos.

            Note: data retrieval might take about half an hour.

            Args:
                chunk_size:
                    Records in one batch. Better stay 100 because higher
                    numbers likely to cause timeout errors.

            Return:
                Records as a list of dicts.
            """

        def query_func(chunk_size: int = 100) -> list:

            return pharos_targets(chunk_size = chunk_size, **args)


        query_func.__name__ = name
        query_func.__doc__ = doc

        globals()[name] = query_func


_create_query_functions()
