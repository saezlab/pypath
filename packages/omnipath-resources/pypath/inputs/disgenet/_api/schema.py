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

from __future__ import annotations

import pypath.share.common as common
from .. import _valid
from .._records import (
    DiseaseDiseaseAssociation,
    GeneDiseaseAssociation,
    VariantDiseaseAssociation,
)

from typing import Any, Callable, Generator, Literal


ARG_DOCS = {
    'email': 'Email address of an account registered at DisGeNet.',
    'password': 'Password for the DisGeNet account.',
    'query_param': 'Query string (HTTP GET) parameters for the DisGeNet API. '
        'See details below, or in the official API docs: '
        'https://www.disgenet.org/api/',
    'diseases': 'One or more disease identifiers.',
    'id_type': 'Disease ID type.',
    'source': 'Restrict query to this resource.',
    'gene': 'Restrict query to this gene (Entrez or HGNC ID).',
    'variant': 'Restrict the query to this variant.',
    'by': 'Group the assciations by this kind of entity.',
    'fail_on_invalid': 'Raise error if the query validation fails.',
    'raw': 'Yield raw records, i.e. the JSON returned by the server converted '
        'to dict, instead of processed named tuples.',
}


QUERY_PARAM_DOCS = {
    'min_score': 'Min value of the {entity_type}-disease score range.',
    'max_score': 'Max value of the {entity_type}-disease score range.',
    'min_ei': 'Min value of the evidence index range.',
    'max_ei': 'Max value of the evidence index range.',
    'disease_type': 'DisGeNET disease type.',
    'disease_class': 'MeSH disease classes.',
    'min_dsi': 'Min value of the DSI range for the {entity_type}.',
    'max_dsi': 'Max value of the DSI range for the {entity_type}.',
    'min_dpi': 'Min value of the DPI range for the {entity_type}.',
    'max_dpi': 'Max value of the DPI range for the {entity_type}.',
    'min_pli': 'Min value of the pLI range.',
    'max_pli': 'Max value of the pLI range.',
    'min_year': 'Ignore evidences from before this year.',
    'max_year': 'Ignore evidences from after this year.',
    'offset': 'Start from this record.',
    'limit': 'Number of {record_type}s to retrieve.',
    'pvalue': 'P-value of the Jaccard index based on the shared genes.',
    'source': 'Restrict query to this resource.',
}


def ARG_TYPES(
    email: str | None = None,
    password: str | None = None,
    query_type: Literal['gda', 'vda'] | None = None,
    query_param: dict[str, Any] | None = None,
    diseases: str | list[str] | None = None,
    id_type: str | None = None,
    gene: str | None = None,
    variant: str | None = None,
    by: str | None = None,
    fail_on_invalid: bool = True,
    raw: bool = False,
): pass


def QUERY_PARAM_TYPES(
    source: str | None = None,
    min_score: float | None = None,
    max_score: float | None = None,
    min_ei: float | None = None,
    max_ei: float | None = None,
    disease_type: Literal['disease', 'phenotype', 'group'] | None = None,
    disease_class: str | list[str] | None = None,
    min_dsi: float | None = None,
    max_dsi: float | None = None,
    min_dpi: float | None = None,
    max_dpi: float | None = None,
    min_pli: float | None = None,
    max_pli: float | None = None,
    limit: int | None = None,
    pvalue: float | None = None,
    min_year: int | None = None,
    max_year: int | None = None,
    offset: int | None = None,
): pass


def RETURN_TYPES(
    dda: Generator[DiseaseDiseaseAssociation],
    gda: Generator[GeneDiseaseAssociation],
    vda: Generator[VariantDiseaseAssociation],
    evs: Generator[GeneDiseaseAssociation | VariantDiseaseAssociation],
): pass


def _argnames(fun: Callable) -> list[str]:

    return fun.__code__.co_varnames[:fun.__code__.co_argcount]


def _argnames_except(
        fun: Callable,
        exc: str | set[str] | None = None,
    ) -> set[str]:

    return set(_argnames(fun)) - common.to_set(exc)


def _query_param(exc: str | set[str] | None = None) -> set[str]:

    return _argnames_except(QUERY_PARAM_TYPES, exc)


def _args(exc: str | set[str] | None = None) -> set[str]:

    return _argnames_except(ARG_TYPES, exc)


def _return(key) -> str:

    return RETURN_TYPES.__annotations__[key]


QUERY_PARAM_EXCEPT = {
    'offset',
    'min_year',
    'max_year',
    'pvalue',
}


RECORD_LABELS = {
    'dda': 'disease-disease association',
    'gda': 'gene-disease association',
    'vda': 'variant-disease association',
}


API_FUNCTIONS = {
    'gene_disease': {
        'args': {'query_type': 'gda'},
        'args_user': _args({'variant'}),
        'query_param_user': _query_param(QUERY_PARAM_EXCEPT),
        'return': _return('gda'),
    },
    'variant_disease': {
        'args': {'query_type': 'vda'},
        'args_user': _args(),
        'query_param_user': _query_param(
            QUERY_PARAM_EXCEPT |
            {'min_pli', 'max_pli'}
        ),
        'return': _return('vda'),
    },
    'disease_disease': {
        'args': {'query_type': 'dda'},
        'args_user': _args({'gene', 'variant'}),
        'query_param_user': {'limit', 'pvalue'},
        'return': _return('dda'),
    },
    'evidences': {
        'args': {'evidences': True},
        'args_user': _args(),
        'query_param_user': {
            'offset', 'min_year', 'max_year',
            'limit', 'min_score', 'max_score',
        },
        'return': _return('evs'),
    },
}
