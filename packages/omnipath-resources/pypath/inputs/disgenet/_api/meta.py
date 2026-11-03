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

from typing import Callable, Generator, NamedTuple

import functools
import textwrap

import pypath.inputs.disgenet._client as _client
import pypath.inputs.disgenet._api.schema as _schema


class ApiMeta:


    def __init__(self, key: str) -> None:

        self._key = key
        self._schema = _schema.API_FUNCTIONS[key]


    def main(self) -> Callable:

        self._prepare()
        self._function()
        self._name()
        self._signature()
        self._docs()
        self._finalize()

        return self._api_function_final


    def _prepare(self) -> None:

        pass


    def _function(self) -> None:

        inner_args = self._schema['args']
        argnames = self._schema['args_user'] - set(inner_args.keys())
        argnames = [
            a for a in _schema._argnames(_schema.ARG_TYPES)
            if a in argnames
        ]
        inner_query_param = self._schema.get('query_param', {})

        def api_function(*args, **kwargs):

            maybe_positional = [a for a in argnames if a not in kwargs]

            if len(maybe_positional) < len(args):

                raise TypeError('Too many positional arguments')

            # arguments
            kw_argnames = kwargs.keys() & set(argnames)
            valid_kwargs = {a: kwargs.pop(a) for a in kw_argnames}
            args = {**dict(zip(maybe_positional, args)), **valid_kwargs}
            args.update(inner_args)
            id_key = 'gene' if args['id_type'] == 'uniprot' else 'diseases'
            args['identifiers'] = args.pop(id_key)
            raw = args.pop('raw', False)

            # query_param
            query_param = {**(args.pop('query_param') or {}), **kwargs}
            query_param = {
                k: v for k, v in query_param.items()
                if v in self._schema.get('query_param_user', {})
            }
            query_param.update(self._schema.get('query_param', {}))
            query_param = {
                k: str(v) for k, v in query_param.items()
                if v is not None
            }
            forbidden = set(inner_query_param.keys()) & set(query_param.keys())

            if forbidden:

                raise TypeError(
                    'Invalid parameters for this query: ' +
                    ', '.join(forbidden)
                )

            query_param.update(inner_query_param)

            client = _client.DisgenetClient(query_param = query_param, **args)

            return client.iterraw() if raw else client.__iter__()

        self._api_function = api_function
        self._argnames = argnames


    def _name(self) -> None:

        self._api_function.__name__ = f'disgenet_{self._key}'


    def _signature(self) -> None:

        atypes = _schema.ARG_TYPES
        rtypes = _schema.RETURN_TYPES

        sig = ', '.join(
            f'{a}: {atypes.__annotations__[a]} = {default}'
            for a, default in self._arg_defaults(atypes, self._argnames)
        )
        sig += ', **kwargs'
        ret = rtypes.__annotations__.get(self._qtype, 'Generator[NamedTuple]')
        args = ', '.join(f'{a} = {a}' for a in self._argnames)

        self._fun_sig_str = (
            f'def _api_func_sig({sig}) -> {ret}: '
            f'return api_func({args}, **kwargs)'
        )


    def _docs(self) -> None:

        space = lambda n: ' ' * n
        record_label = _schema.RECORD_LABELS.get(self._qtype, 'data')
        docs = f'Retrieve {record_label}s from the DisGeNet API.\n\nArgs:\n'

        for a in self._argnames:

            a_doc = _schema.ARG_DOCS.get(a, 'Not documented.')
            a_doc = textwrap.wrap(a_doc, 66)
            a_doc = f'\n{space(8)}'.join(a_doc)
            docs += f'{space(4)}{a}:\n'
            docs += f'{space(8)}{a_doc}\n'

        docs += '\n'
        docs += f'Details:\n{space(4)}'
        details_head = (
            'The following arguments can be provided in the '
            '`query_param` dict or as extra keyword arguments.'
        )
        docs += f'\n{space(4)}'.join(textwrap.wrap(details_head, 70))
        docs += '\n'

        qtypes = _schema.QUERY_PARAM_TYPES.__annotations__
        qdefaults = dict(self._arg_defaults(_schema.QUERY_PARAM_TYPES))

        for q in self._schema.get('query_param_user', ()):

            q_doc = _schema.QUERY_PARAM_DOCS[q].format(
                record_type = record_label,
                entity_type = record_label.split('-', maxsplit = 1)[0],
            )
            q_doc = f'* {q}: {qtypes[q]} = {qdefaults[q]} -- {q_doc}'
            q_doc = textwrap.wrap(q_doc, 70)
            q_doc = f'\n{space(6)}'.join(q_doc)
            docs += f'{space(4)}{q_doc}\n'

        docs += '\nYields:\n'
        docs += (
            f'{space(4)}Named tuples or dicts '
            f'representing {record_label}s.\n'
        )
        docs = textwrap.indent(docs, space(4))
        docs = f'\n{docs}'

        self._api_function.__doc__ = docs


    def _finalize(self) -> None:

        fun_dict = {}
        exec(self._fun_sig_str, {'api_func': self._api_function}, fun_dict)
        with_signature = fun_dict['_api_func_sig']
        with_signature.__doc__ = self._api_function.__doc__
        with_signature.__name__ = self._api_function.__name__

        self._api_function_final = with_signature


    @staticmethod
    def _arg_defaults(
            fun: Callable,
            include: set[str] | None = None,
        ) -> Generator[tuple[str, Any]]:

        for a, default in zip(_schema._argnames(fun), fun.__defaults__):

            if include is None or a in include:

                yield a, default


    @property
    def _qtype(self) -> str | None:

        return self._schema.get('args', {}).get('query_type', None)


for key in _schema.API_FUNCTIONS.keys():

    _fun = ApiMeta(key).main()
    globals()[_fun.__name__] = _fun
    del _fun
