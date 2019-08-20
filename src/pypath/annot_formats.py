#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Meta-classes used for annotation functionality
#
#  Copyright
#  2014-2019
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
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

import collections

"""
Annotations are defined by a ``name``, a ``source`` and an ``args``
parameter. If the former is a string it will be first looked up among the
annotation resources in ``pypath.annot.db``. Otherwise among the keys of the
classes in the ``CustomAnnotation`` object or in the dictionary of the
class_definitions. If ``source`` is a `set`, it will be used as a category
without further processing. If it is callable it will be called and should
return a `set`. If ``bool(args)`` is `False`, in case of annotations in
``pypath.annot.db`` the ``to_set`` method will be called. Otherwise ``args``
will be passed to the ``get_subset`` method. If ``source`` is callable,
``args`` will be passed if available.
"""
AnnotDef = collections.namedtuple(
    'AnnotDef',
    ['name', 'source', 'args'],
)
AnnotDef.__new__.__defaults__ = (None,)

"""
Annotation operations consist of list of annotation definitions or names as
they can be looked up in the ``class_definitions`` of the
``CustomAnnotation`` object and an operator to be called on the sets
(union, intersection or difference).
"""
AnnotOp = collections.namedtuple(
    'AnnotOp',
    ['annots', 'op'],
)
