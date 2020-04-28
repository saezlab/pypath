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
import pypath.share.common as common
import pypath.core.entity as entity


class AnnotDef(
        collections.namedtuple(
            'AnnotDefBase',
            [
                'name',
                'resource',
                'parent',
                'aspect',
                'scope',
                'args',
                'source',
                'exclude',
                'transmitter',
                'receiver',
                'resource_name',
            ]
        )
    ):
    """
    Annotations are defined by a ``name``, a ``resource`` and an ``args``
    parameter. If the former is a string it will be first looked up among the
    annotation resources in ``pypath.annot.db``. Otherwise among the keys of
    the classes in the ``CustomAnnotation`` object or in the dictionary of
    the class_definitions. If ``source`` is a `set`, it will be used as a
    category without further processing. If it is callable it will be called
    and should return a `set`. If ``bool(args)`` is `False`, in case of
    annotations in ``pypath.annot.db`` the ``to_set`` method will be called.
    Otherwise ``args`` will be passed to the ``get_subset`` method.
    If ``resource`` is callable, ``args`` will be passed if available.
    """


    def __new__(
            cls,
            name,
            resource,
            parent = None,
            aspect  = 'functional',
            scope = 'specific',
            source = 'resource_specific',
            args = None,
            exclude = None,
            transmitter = None,
            receiver = None,
            resource_name = None,
        ):

        return super().__new__(
            cls,
            name = name,
            resource = resource,
            parent = parent,
            aspect = aspect,
            scope = scope,
            args = args,
            source = source,
            exclude = exclude,
            transmitter = transmitter,
            receiver = receiver,
            resource_name = resource_name,
        )


    @property
    def key(self):

        return (
            self.name,
            self.parent,
            self.resource_str,
        )


    @property
    def resource_str(self):

        return (
            self.resource
                if isinstance(self.resource, common.basestring) else
            self.resource_name
        )


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
    """
    Represents a set of molecular entities sharing a custom defined
    annotation. This class behaves like a ``set`` and set operations on it
    result set objects. Normally this class is instantiated by
    ``pypath.core.annot.CustomAnnotation`` in the process of populating
    categories and the contents of the groups defined in
    ``pypath.core.intercell_annot`` in case of annotations of the
    intercellular communication roles.
    For detailed definitions of the parameter values see the Supplementary
    Table S10 in Turei et al. 2020 (in prep).

    :param list,set,tuple members:
        The identifiers of the entities in the category.
    :param str name:
        The name of the category.
    :param str parent:
        The name of the parent category; might be the same as ``name`` in
        case of high level (generic) categories.
    :param str aspect:
        Either *functional* or *locational*.
    :param str source:
        Either *resource_specific* or *composite*.
    :param str scope:
        Either *specific* or *generic*.
    :param str resource:
        The resource (database) name; in case of composite categories it
        should be the name of the database you are actually building, this
        by default is `OmniPath` and you can change by the
        ``pypath.share.settings`` module using the
        ``annot_composite_database_name`` key.
    :param bool transmitter:
        Whether the category contains transmitters of signaling information
        from the cell expressing the molecular entities in direction of other
        cells.
    :param bool receiver:
        Whether the category contains receivers of signaling information
        from other cells in direction of the cells expressing the molecular
        entites in the category.
    """

    def __init__(
            self,
            members,
            name = None,
            parent = None,
            aspect = 'functional',
            source = 'resource_specific',
            scope = 'specific',
            resource = None,
            transmitter = None,
            receiver = None,
        ):

        collections_abc.Set.__init__(self)
        self.members = set(members)
        self.name = name or 'unnamed'
        self.parent = parent or self.name
        self.aspect = aspect
        self.source = source
        self.scope = scope
        self.resource = (
            resource or
            settings.get('annot_composite_database_name') or
            'Unknown'
        )
        self.transmitter = transmitter
        self.receiver = receiver


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


    @property
    def label(self):

        return (
            '%s%s__%s' % (
                self.name,
                '__%s' % self.parent if self.parent != self.name else '',
                self.resource,
            )
        )


    @property
    def key(self):

        return (
            self.name,
            self.parent,
            self.resource,
        )


    def filter_entity_type(self, entity_type = None):
        """
        Returns a copy of the group with only the selected entity types.
        If ``entity_type`` is None returns the object itself.
        """

        if entity_type is None:

            return self

        else:

            members = entity.Entity.filter_entity_type(
                self.members,
                entity_type = entity_type,
            )

            return AnnotationGroup(
                members = members,
                name = self.name,
                parent = self.parent,
                aspect = self.aspect,
                source = self.source,
                scope = self.scope,
                resource = self.resource,
                transmitter = self.transmitter,
                receiver = self.transmitter,
            )


    def count_entity_type(self, entity_type = None):

        return len(
            entity.Entity.count_entity_type(
                self.members,
                entity_type = entity_type,
            )
        )
