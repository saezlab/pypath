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
Methods working on objects from `pypath.core` which might be useful to access
without importing the main modules.
"""

from future.utils import iteritems

import re

import pandas as pd

import pypath.share.common as common
import pypath_common._constants as _const


def filter_network_df(
        df,
        resource = None,
        entity_type = None,
        data_model = None,
        interaction_type = None,
        only_directed = None,
        only_undirected = None,
        only_signed = None,
        only_proteins = None,
        effect = None,
        entities = None,
        source_entities = None,
        target_entities = None,
        swap_undirected = True,
        remove_loops = True,
        entities_or = False,
        **kwargs
    ):
    """
    Filters a network data frame.
    """

    args = locals().copy()
    args.update(kwargs)

    query_elements = {}

    _node_attrs = {'id', 'entity_type'}

    _synonyms = {
        'resource': 'sources',
        'data_model': 'dmodel',
        'interaction_type': 'type',
        'entities': 'id',
        'entity': 'id',
    }

    if swap_undirected and not only_directed and not only_signed:

        undirected = df.query('not directed', inplace = False)

        with pd.option_context('mode.chained_assignment', None):

            undirected.rename(
                {
                    'id_a': 'id_b',
                    'id_b': 'id_a',
                    'type_a': 'type_b',
                    'type_b': 'type_a',
                },
                inplace = True,
                axis = 'columns',
            )

        df = df.append(undirected, ignore_index = True, sort = False)

    if only_directed:

        args['directed'] = True

    if only_signed:

        args['effect'] = {1, -1}

    if only_undirected:

        args['directed'] = False

    if only_proteins:

        args['entity_type'] = 'protein'

    for var, val in iteritems(args):

        if val is None:

            continue

        node_postfix = None

        if var.startswith('source_') or var.endswith('_a'):

            node_postfix = ['_a']
            var = re.sub(r'(?:source_)?(\w+)(?:_a)?', r'\1', var)

        if var.startswith('target_') or var.endswith('_b'):

            node_postfix = ['_b']
            var = re.sub(r'(?:target_)?(\w+)(?:_b)?', r'\1', var)

        if var == 'type':

            raise ValueError(
                '`type` is ambiguous, use either '
                '`interaction_type` or `entity_type`.'
            )


        var = _synonyms[var] if var in _synonyms else var

        node_postfix = (
            ['_a', '_b']
                if var in _node_attrs and not node_postfix else
            node_postfix or ['']
        )

        for pf in node_postfix:

            var_pf = '%s%s' % (var, pf)
            var_pf = var_pf.replace('entity_type', 'type')

            if var_pf in df.columns:

                query_elements[var_pf] = val

    query = []
    entity_query = []

    for var, val in iteritems(query_elements):

        first = var
        second = '@query_elements["%s"]' % var

        if df.dtypes[var].name == 'object':

            first, second = second, first
            op = 'in' if isinstance(val, _const.SIMPLE_TYPES) else '&'

        else:

            op = '==' if isinstance(val, _const.SIMPLE_TYPES) else 'in'

        q = '%s %s %s' % (first, op, second)

        if var.startswith('id_'):

            entity_query.append(q)

        else:

            query.append(q)

    if entities_or and len(entity_query) > 1:

        query.append('(%s)' % ' or '.join(entity_query))

    else:

        query.extend(entity_query)

    if remove_loops:

        query.append('id_a.astype("str") != id_b.astype("str")')

    query = ' and '.join(query)

    result = df.query(query) if query else df

    return result
