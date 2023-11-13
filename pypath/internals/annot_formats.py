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

import collections

try:
    import collections.abc as collections_abc
except:
    import collections as collections_abc

import pypath.share.settings as settings
import pypath.share.common as common
import pypath.core.entity as entity


AnnotDefKey = collections.namedtuple(
    'AnnotDefKey',
    [
        'name',
        'parent',
        'resource',
    ],
)


class AnnotDef(
        collections.namedtuple(
            'AnnotDefBase',
            [
                'name',
                'resource',
                'parent',
                'aspect',
                'scope',
                'source',
                'args',
                'exclude',
                'transmitter',
                'receiver',
                'resource_name',
                'limit',
                'avoid',
                'enabled',
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
            limit = None,
            avoid = None,
            enabled = True,
        ):

        resource_name = (
            resource
                if cls._is_resource_name(resource) else
            (
                resource_name or
                settings.get('annot_composite_database_name') or
                'Unknown'
            )
        )

        return super().__new__(
            cls,
            name = name,
            resource = resource,
            parent = parent or name,
            aspect = aspect,
            scope = scope,
            source = source,
            args = args,
            exclude = exclude,
            transmitter = transmitter,
            receiver = receiver,
            resource_name = resource_name,
            limit = cls._zero_one_or_more(limit),
            avoid = cls._zero_one_or_more(avoid),
            enabled = enabled,
        )


    @property
    def key(self):

        return AnnotDefKey(
            name = self.name,
            parent = self.parent,
            resource = self.resource_name,
        )


    @staticmethod
    def _is_resource_name(name):

        return (
            isinstance(name, str) and
            not (
                name.startswith('~') or
                name.startswith('#')
            )
        )


    @staticmethod
    def _zero_one_or_more(arg):

        return (
            ()
                if not arg else
            (arg,)
                if isinstance(arg, (str, _annot_type)) else
            arg
        )



class AnnotOp(
        collections.namedtuple(
            'AnnotOpBase',
            ['annots', 'op'],
        )
    ):
    """
    Annotation operations consist of list of annotation definitions or names as
    they can be looked up in the ``class_definitions`` of the
    ``CustomAnnotation`` object and an operator to be called on the sets
    (union, intersection or difference).
    """

    def __new__(cls, annots, op = set.union):

        if op in AnnotationGroup._set_methods:

            op = getattr(AnnotationGroup, op.__name__)

        return super().__new__(cls, annots, op)


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
    :param str,set limit:
        Limit to this or these categories. E.g. if it's 'extracellular'
        the result will be the intersection of this category and
        'extracellular' i.e. the category will be limited to extracellular
        proteins.
    :param str,set avoid:
        Avoid elements of this or these categories. E.g. if it's
        'cell_surface' then all cell_surface proteins will be removed from
        this category.
    """


    _set_methods = {
        set.union,
        set.intersection,
        set.difference,
        set.symmetric_difference,
    }


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
            limit = None,
            avoid = None,
            enabled = True,
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
        self.limit = common.to_set(limit)
        self.avoid = common.to_set(avoid)
        self.enabled = enabled


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

        return '%s%s@%s' % (
            self.parent,
            '::%s' % self.name if self.name != self.parent else '',
            self.resource
        )


    @property
    def name_label(self):

        return common.upper0(self.name).replace('_', ' ')


    @property
    def key(self):

        return AnnotDefKey(
            name = self.name,
            parent = self.parent,
            resource = self.resource,
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

        return (
            entity.Entity.count_entity_type(
                self.members,
                entity_type = entity_type,
            )
        )


    @property
    def args(self):

        return dict(**self)


    def keys(self):

        return ['name', 'parent', 'source', 'scope', 'aspect']


    def __getitem__(self, key):

        return getattr(self, key)


    @property
    def n_proteins(self):

        return self.count_entity_type(entity_type = 'protein')


    @property
    def n_mirnas(self):

        return self.count_entity_type(entity_type = 'mirna')


    @property
    def n_complexes(self):

        return self.count_entity_type(entity_type = 'complex')


    @property
    def proteins(self):

        return self.filter_entity_type(entity_type = 'protein')


    @property
    def complexes(self):

        return self.filter_entity_type(entity_type = 'complex')


    @property
    def mirnas(self):

        return self.filter_entity_type(entity_type = 'mirna')


    @staticmethod
    def sets(*args):

        return (
            (
                a
                    if isinstance(a, set) else
                a.members
                    if hasattr(a, 'members') else
                common.to_set(a)
            )
            for a in args
        )


    @classmethod
    def union(cls, *args):

        return set.union(*cls.sets(*args))


    @classmethod
    def intersection(cls, *args):

        return set.intersection(*cls.sets(*args))


    @classmethod
    def difference(cls, *args):

        return set.difference(*cls.sets(*args))


    @classmethod
    def symmetric_difference(cls, *args):

        return set.symmetric_difference(*cls.sets(*args))


    @classmethod
    def isdisjoint(cls, *args):

        return set.isdisjoint(*cls.sets(*args))


_set_type = (set, AnnotationGroup)
_annot_type = _set_type + (AnnotDef, AnnotOp)
