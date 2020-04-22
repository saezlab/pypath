#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Meta-classes used for annotation functionality
#
#  Copyright
#  2014-2020
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

try:
    import collections.abc as collections_abc
except:
    import collections as collections_abc

import pypath.share.settings as settings

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
    ['name', 'resource', 'source', 'aspect', 'scope', 'args', 'exclude'],
)
AnnotDef.__new__.__defaults__ = (None,) * 6

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


class AnnotationGroup(collections_abc.Set):


    def __init__(
            self,
            members,
            name = None,
            aspect = 'functional',
            source = 'resource_specific',
            scope = 'specific',
            resource = None,
        ):

        collections_abc.Set.__init__(self)
        self.members = set(members)
        self.name = name or 'unnamed'
        self.aspect = aspect
        self.source = source
        self.scope = scope
        self.resource = (
            resource or
            settings.get('annot_composite_database_name') or
            'Unknown'
        )


    def __iter__(self):

        return self.members.__iter__()


    def __contains__(self, other):

        return other in self.members


    def __len__(self):

        return len(self.members)


    @classmethod
    def _from_iterable(cls, iterable):

        return set(iterable)


    def __repr__(self):

        return '<AnnotationGroup `%s` from %s, %u elements>' % (
            self.name,
            self.resource,
            len(self),
        )
