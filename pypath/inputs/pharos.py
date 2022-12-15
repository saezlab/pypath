
import json

from pypath.share.curl import Curl
from pypath.resources.urls import urls


def pharos_general(query, variables=None):

    url = urls['pharos_api']['url']

    req_headers = {
        "Accept-Encoding": "gzip, deflate, br",
        "Content-Type": "application/json",
        "Connection": "keep-alive",
        "DNT": "1",
        "Origin": "https://pharos-api.ncats.io",
    }

    if variables:

        binary_data = json.dumps(
            {
                'query': query,
                'variables': variables
            }
        )

    else:

        binary_data = json.dumps(
            {
                'query': query,
            }
        )

    binary_data = binary_data.encode('utf-8')

    c = Curl(
        url=url,
        req_headers=req_headers,
        binary_data=binary_data,
        compressed=True,
        compr='gz',
    )
    result = json.loads(c.result)['data']

    return result


def get_targets(
    # It better stay 100 because higher numbers likely to cause timeout errors
    chunk_size=100,
    getExpressions=False,
    getGtex=False,
    getOrthologs=False,
    getLigands=False,
    getXrefs=False,
    getDiseases=False,
):

    step = 0
    result = list()
    variables = {
        "chunk_size": chunk_size,
        "step": step,
        "getExpressions": getExpressions,
        "getGtex": getGtex,
        "getOrthologs": getOrthologs,
        "getLigands": getLigands,
        "getXrefs": getXrefs,
        "getDiseases": getDiseases,
    }

    query = """
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

                        activities(all: true) {
                            actid
                            type
                            moa
                            value
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

    while True:

        response = pharos_general(query, variables)["targets"]["targets"]

        if not response:
            break

        result.extend(response)
        variables["step"] += chunk_size

    return result


def getExpressions(chunk_size=100):
    return get_targets(chunk_size, getExpressions=True)

def getGtex(chunk_size=100):
    return get_targets(chunk_size, getGtex=True)

def getOrthologs(chunk_size=100):
    return get_targets(chunk_size, getOrthologs=True)

def getLigands(chunk_size=100):
    return get_targets(chunk_size, getLigands=True)

def getXrefs(chunk_size=100):
    return get_targets(chunk_size, getXrefs=True)

def getDiseases(chunk_size=100):
    return get_targets(chunk_size, getDiseases=True)
