#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Sebastian Lobentanzer
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from __future__ import annotations

"""
Utilities for XML processing.

Here you find a few functions that might make it easier to process XML
with `lxml.etree`. The alternatives of these are be xpath and the plain
etree API. You can combine all these three, as you can embed `xpath`
expressions and any API call or custom external procedures in the functions
below.
"""

from typing import Iterable, Iterator

import functools
import collections

from lxml import etree

import pypath.share.common as common

EtreeMethod = collections.namedtuple(
    'EtreeMethod',
    (
        'value',
        'subject',
        'direct',
        'tag_arg',
    ),
    defaults = (Iterator, 'descendants', True, None),
)

METHODS = {
    'find': EtreeMethod(etree.Element, tag_arg = True),
    'findall': EtreeMethod(list, tag_arg = True),
    'getchildren': EtreeMethod(list, tag_arg = False),
    'iterchildren': EtreeMethod(),
    'iterfind': EtreeMethod(tag_arg = True),
    'iterdescendants': EtreeMethod(direct = False),
    'itersiblings': EtreeMethod(subject = 'siblings'),
    'iterancestors': EtreeMethod(subject = 'ancestors'),
    'iter': EtreeMethod(direct = False),
    'getnext': EtreeMethod(etree.Element, tag_arg = False),
    'findtext': EtreeMethod(str, tag_arg = True),
    'itertext': EtreeMethod(Iterable, 'descendants', False, None),
    'xpath': EtreeMethod(common.to_list, direct = False),
}


def fetch(
        elem: etree.ElementBase,
        path: Iterable,
        namespaces: str | dict | None = None,
        clear: bool = True,
    ) -> str | list | dict | None:
    """
    Fetch data from an `etree.Element`.

    Args:
        elem:
            The `etree.Element` to fetch data from.
        path:
            Path specification of the data to fetch.
            It is a list of lookup steps that ultimately lead to the desired
            data. One step can be defined as simple as a string, in this case
            the string is a tag name and the `find` method is called to access
            it. It other method, such as `findall` or `iterchildren` should be
            used, the step might consist of the tag name and the method name,
            such as `('tag', 'iterchildren')`. Methods yielding multiple tags
            will create arrays in the output. If the step is a dict, multiple
            methods can be applied on the same element, and an associative
            array will be created in the output. Values of such a dict follow
            the same format as the whole path definition. The last step of
            the path actually extracts the data. If it is None, the `text`
            attribute will be accessed. Alternatively, an array of attribute
            names van be provided, creating an associative array of attributes
            in the output. Namespaces can be controlled by the third elements
            of the tuples specifying each step.
        namespaces:
            One namespace or a dict of namespaces. If a string is provided, it
            will be used as namespace for all tags. For a dict of namespaces,
            the third elements in each path step tuple point to namespaces by
            keys of the dict.

    Returns:
        The data extracted from the `etree.Element`. For content of a single
        tag the value is a string; for multiple tags it is a list, for
        attributes and dict step specifications dicts are returned.
    """

    namespaces = namespaces or ''
    namespaces = (
        {None: namespaces}
            if isinstance(namespaces, str) else
        namespaces
    )

    if isinstance(path, dict):

        fetch_ns = functools.partial(
            fetch,
            elem,
            namespaces = namespaces,
            clear = False,
        )
        return common.compr(path, apply = fetch_ns)

    print(f'processing path: {path}')
    step_ns = functools.partial(step, namespaces = namespaces)
    result = functools.reduce(step_ns, path, elem)

    if clear:

        elem.clear(keep_tail = True)

    return result


def step(
        elem: etree._Element,
        spec: str | tuple | dict,
        namespaces: str | dict | None = None,
    ) -> str | list | dict:
    """
    Dispatch the execution of one processing step.
    """

    print(f'calling step, spec: {spec}')

    if isinstance(spec, set):

        spec = _simple_fields(spec)

    if isinstance(elem, etree._Element):

        if isinstance(spec, dict):

            print('applying array of steps on a single element')
            result = {
                key: fetch(elem, path, namespaces, clear = False)
                for key, path in spec.items()
            }
            return result

        else:

            print(f'calling simple step, spec = {spec}')
            result = simple_step(elem, spec, namespaces)
            print(f'result: {result}')
            return result

    elif isinstance(elem, common.simple_types):

        return elem

    else:

        print(f'applying step on array, spec = {spec}')
        print(f'the array: {elem}')
        this_step = functools.partial(
            step,
            spec = spec,
            namespaces = namespaces,
        )
        result = common.compr(elem, apply = this_step)
        print(f'result: {result}')
        return result


def simple_step(
        elem: etree._Element | None,
        spec: str | tuple | dict,
        namespaces: str | dict | None = None,
    ) -> str | list | dict:
    """
    Execute one processing step.
    """

    if not isinstance(elem, etree._Element):

        return elem

    if spec is None:

        return elem.text

    if callable(spec):

        return spec(elem)

    tag, method, ns, *_ = common.to_tuple(spec) + (None,) * 3
    method = method or 'find'
    ns = ns or namespaces.get(None, '')
    param = METHODS.get(method)
    args = (f'{ns}{tag}',) if tag and param.tag_arg is not False else ()
    expand = list if param.value == Iterable else common.identity

    return expand(getattr(elem, method)(*args))


def _simple_fields(fields: Iterable[str]) -> dict[str, tuple]:
    """
    Expand a list of fields names to make possible a more concise notation.
    """

    return {f: (f, None) for f in fields}


def _array_fields(fields: Iterable[tuple[str, str]]) -> dict[str, tuple]:
    """
    Expand a list of flat and homogenous array type fields.
    """

    return {f[0]: (f[0], (f[1], 'findall'), None) for f in fields}
