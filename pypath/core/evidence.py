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
Provides classes for representing and processing evidences supporting
relationships. The evidences hold information about the databases and
literature references, they can be organized into collections. A number
of operations are available on evidences and their collections, for
example they can be combined or filtered.
"""

from __future__ import annotations

from future.utils import iteritems

import importlib as imp
import copy

import pypath.internals.refs as refs
import pypath.share.common as common
import pypath.share.session as session_mod
import pypath.share.constants as _const
import pypath.core.entity as entity
import pypath.core.attrs as attrs_mod
import pypath.resources.network as netres
import pypath_common._constants as _const

_logger = session_mod.Logger(name = 'evidence')
_log = _logger._log


class Evidence(attrs_mod.AttributeHandler):
    """
    Represents an evidence supporting a relationship such as molecular
    interaction, molecular complex, enzyme-PTM interaction, annotation, etc.

    The evidence consists of two main parts: the database and the literature
    references. If a relationship is supported by multiple databases, for
    each one `Evidence` object should be created and

    :arg pypath.resource.ResourceAttributes resource:
        An object derived from :py:class:`pypath.resource.ResourceAttributes`.
    :arg str,list,set,NoneType references:
        Optional, one or more literature references (preferably PubMed IDs).
    """

    __slots__ = [
        'resource',
        'references',
        'dataset',
    ]


    def __init__(self, resource, references = None, attrs = None):

        self.resource = resource
        self.dataset = getattr(resource, 'dataset', None)
        self.references = self._process_references(references)
        attrs_mod.AttributeHandler.__init__(self, attrs)


    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    @staticmethod
    def _process_references(references):

        references = common.to_set(references)

        return (
            set(
                (
                    refs.Reference(ref)
                        if not isinstance(ref, refs.Reference) else
                    ref
                )
                for ref in references
            )
        )


    def __hash__(self):

        return self.resource.__hash__()


    def __eq__(self, other):

        return (
            self.resource == other or
            (
                hasattr(other, 'resource') and
                self.resource == other.resource and
                (
                    self.resource.interaction_type ==
                    self.resource.interaction_type
                )
            )
        )


    def __iadd__(self, other):
        """
        This will ignore if the other evidence is from different resource:
        still better than attributing wrong references to a resource.
        """

        if self == other:

            self.dataset = netres.choose_dataset(self.dataset, other.dataset)
            self.references.update(other.references)
            attrs_mod.AttributeHandler.__iadd__(other)

        else:

            _log(
                'Warning: attempt to merge evidences from different '
                'resources. Ignoring the second evidence.'
            )

        return self


    def __add__(self, other):

        dataset = netres.choose_dataset(self.dataset, other.dataset)

        new = self.__class__(
            resource = self.resource,
            references = self.references | other.references,
        )
        new.dataset = dataset
        new.update_attrs(self.attrs.copy())
        new.update_attrs(other.attrs.copy())

        return new


    @property
    def key(self):

        return self.resource.key


    def merge(self, other):
        """
        Merges two evidences. Returns set of either one or two evidences
        depending on whether the two evidences are from the same resource.
        """

        if self == other:

            self += other
            return {self}

        else:

            return {self, other}


    def __repr__(self):

        return '<Evidence %s (%s%u references)>' % (
            self.resource.name,
            'via %s, ' % self.resource.via if self.resource.via else '',
            len(self.references),
        )


    def __str__(self):

        return self.resource.name


    def __copy__(self):

        return self.__class__(
            resource = self.resource,
            references = copy.copy(self.references),
            attrs = self.attrs.copy(),
        )


    def __contains__(self, other):
        """
        :arg str,tuple,Reference other:
            Either a reference or a database name, or a tuple of a database
            name and an interaction type or a tuple of a database, interaction
            type and a primary database (or None if the query should be
            limited only to primary databases).
        """

        return (
            self._contains(self, other) or
            attrs_mod.AttributeHandler.__contains__(self, other)
        )


    def contains_database(self, database):

        return self.resource.name == database


    def contains_reference(self, reference):

        return reference in self.references


    def has_database_via(self, database, via):

        return (
            self.resource.name == database and
            self.resource.via == via
        )


    def has_interaction_type(
            self,
            interaction_type,
            database = None,
            via = False,
        ):
        """
        If ``via`` is ``False`` then it will be ignored, otherwise if ``None``
        only primary resources are considered.
        """

        return (
            self.resource.interaction_type == interaction_type and
            (
                not database or
                self.resource.name == database
            ) and
            (
                via == False or
                self.resource.via == via
            )
        )


    @staticmethod
    def _contains(obj, other):

        if isinstance(other, int):

            other = '%u' % other

        if isinstance(other, str) and other.isdigit():

            other = refs.Reference(other)

        if isinstance(other, refs.Reference):

            return obj.contains_reference(other)

        # this makes possible to accept a NetworkResource or a
        # NetworkResourceKey:
        if (
            hasattr(other, 'name') and
            hasattr(other, 'interaction_type') and
            hasattr(other, 'via')
        ):

            other = (other.name, other.interaction_type, other.via)

        other = other if isinstance(other, tuple) else (other,)

        return (
            obj.contains_database(other[0]) and
            (
                len(other) == 1 or
                obj.has_interaction_type(other[1], other[0])
            ) and
            (
                len(other) <= 2 or
                obj.has_database_via(other[0], other[2])
            )
        )


    def has_data_model(self, data_model):

        return self.resource.data_model == data_model


    def match(
            self,
            resource = None,
            data_model = None,
            interaction_type = None,
            via = False,
            references = None,
            datasets = None,
        ):

        def _match(attr, value):

            return (
                getattr(self.resource, attr) in value
                    if isinstance(value, _const.LIST_LIKE) else
                getattr(self.resource, attr) == value
            )


        resource = (
            resource.resource
                if isinstance(resource, Evidence) else
            resource
        )

        interaction_type = (
            resource.interaction_type
                if (
                    interaction_type is None and
                    hasattr(resource, 'interaction_type')
                ) else
            interaction_type
        )

        via = (
            resource.via
                if (
                    via is None and
                    hasattr(resource, 'via')
                ) else
            via
        )

        data_model = (
            resource.data_model
                if hasattr(resource, 'data_model') else
            data_model
        )

        references = common.to_set(references)

        return (
            (
                resource is None or (
                    self.resource.name in resource
                        if isinstance(resource, set) else
                    self.resource == resource
                )
            ) and
            (
                interaction_type is None or
                _match('interaction_type', interaction_type)
            ) and
            (
                via is None or
                (via == False and not self.resource.via) or
                (via == True and self.resource.via) or
                _match('via', via)
            ) and
            (
                not references or
                self.references & references
            ) and
            (
                not data_model or
                _match('data_model', data_model)
            ) and
            (
                datasets is None or
                _match('dataset', datasets)
            )
        )


    def __str__(self):

        return self.resource.name


    @property
    def pubmeds(self) -> list[str]:
        """
        PubMed IDs of the references supporting this evidence.
        """

        return [r.pmid for r in self.references]


    def asdict(self) -> dict:
        """
        Dictionary representation of the evidence.
        """

        return {
            'resource': self.resource.name,
            'references': self.pubmeds,
            'dataset': self.dataset,
            'via': self.resource.via,
            'attrs': self.attrs,
        }


class Evidences(object):
    """
    A collection of evidences. All evidences supporting a relationship such
    as molecular interaction, molecular complex, enzyme-PTM interaction,
    annotation, etc should be collected in one `Evidences` object. This way
    the set of evidences can be queried a comprehensive way.

    :arg tuple,list,set,Evidences evidences:
        An iterable providing :py:class:`Evidence` instances. It is possible
        to create an empty evidence collection and populate it later or to
        show this way that certain relationship has no supporting evidences.
    """

    __slots__ = [
        'evidences',
    ]


    def __init__(self, evidences = ()):

        self.evidences = {}
        self.__iadd__(evidences)


    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

        new_ev_class = getattr(mod, 'Evidence')

        for ev in self:

            ev.__class__ = new_ev_class


    def __iadd__(self, other):

        if isinstance(other, str):

            other = Evidence(other)

        other = (
            other
                if (
                    hasattr(other, '__iter__') and
                    not isinstance(other, Evidence)
                ) else
            (other,)
                if isinstance(other, self.__class__) else
            ()
        )

        for ev in other:

            if ev.key in self.evidences:

                self.evidences[ev.key] = self.evidences[ev.key] + ev

            else:

                self.evidences[ev.key] = ev.__copy__()

        return self


    def __add__(self, other):

        if not isinstance(other, self.__class__):

            return self.__copy__()

        return Evidences(
            (
                self.evidences[key].__copy__()
                    if key not in other.evidences else
                other.evidences[key].__copy__()
                    if key not in self.evidences else
                self.evidences[key] + other.evidences[key]
            )
            for key in
            set(self.evidences.keys()) | set(other.evidences.keys())
        )


    def __radd__(self, other):

        return self.__add__(other)


    def __sub__(self, other):

        return Evidences(
            ev
            for ev in self
            if ev not in other
        )


    def intersection(self, other):

        return Evidences(
            self.evidences[key] + other.evidences[key]
            for key in
            set(self.evidences.keys()) & set(other.evidences.keys())
        )


    def __iter__(self):

        for ev in self.evidences.values():

            yield ev


    def __repr__(self):

        return '<Evidences: %s (%u references)>' % (
            (
                ', '.join(sorted(set(ev.resource.name for ev in self)))
                    if self else
                'None'
            ),
            (
                len(set.union(*(ev.references for ev in self)))
                    if self else
                0
            ),
        )


    def __copy__(self):

        return Evidences((ev.__copy__() for ev in self))


    def __bool__(self):

        return bool(len(self.evidences))


    def __contains__(self, other):
        """
        :arg str,tuple,Reference other:
            Either a reference or a database name, or a tuple of a database
            name and an interaction type or a tuple of a database, interaction
            type and a primary database (or None if the query should be
            limited only to primary databases).
        """

        return Evidence._contains(self, other)


    def __and__(self, other):

        other = self._foreign_resources_set(other)
        this = self._resident_resources_set(other)

        return this & other


    def __or__(self, other):

        other = self._foreign_resources_set(other)
        this = self._resident_resources_set(other)

        return this | other


    @staticmethod
    def _foreign_resources_set(resources):

        other = common.to_set(resources)

        return {
            (
                res.resource
                    if hasattr(res, 'resource') else
                res
            )
            for res in resources
        }


    def _resident_resources_set(self, other = None):

        return (
            {ev.resource.name for ev in self}
                if (
                    hasattr(other, '__iter__') and
                    all(isinstance(res, str) for res in other)
                ) else
            {ev.resource for ev in self}
        )


    def __eq__(self, other):

        return {ev.resource for ev in self} == {ev.resource for ev in other}


    def __len__(self):

        return self.count_resources()


    def count_resources(self, **kwargs):

        return len(list(self.filter(**kwargs)))


    def get_resources(self, **kwargs):

        return {
            ev.resource
            for ev in self.filter(**kwargs)
        }


    def get_resources_via(self, **kwargs):

        return {
            (ev.resource, ev.resource.via)
            for ev in self.filter(**kwargs)
        }


    def get_resource_names_via(self, **kwargs):

        return {
            (ev.resource.name, ev.resource.via)
            for ev in self.filter(**kwargs)
        }


    def count_references(self, **kwargs):

        return len(self.get_references(**kwargs))


    def get_references(self, **kwargs):

        evidences = self.filter(**kwargs)

        return {
            ref
            for ev in evidences
            for ref in ev.references
        }


    def count_curation_effort(self, **kwargs):

        return len(self.get_curation_effort(**kwargs))


    def get_curation_effort(self, **kwargs):

        evidences = self.filter(**kwargs)

        return {
            (ev.resource, ref)
            for ev in evidences
            for ref in ev.references
        }


    def contains_database(self, database, **kwargs):

        return any(
            ev.resource.name == database
            for ev in self.filter(**kwargs)
        )


    def contains_reference(self, reference, **kwargs):

        return any(reference in ev.references for ev in self.filter(**kwargs))


    def has_database_via(self, database, via, **kwargs):

        return any(
            ev.has_database_via(database, via)
            for ev in self.filter(**kwargs)
        )


    def has_interaction_type(
            self,
            interaction_type,
            database = None,
            via = False,
        ):
        """
        If ``via`` is ``False`` then it will be ignored, otherwise if ``None``
        only primary resources are considered.
        """

        return any(
            ev.has_interaction_type(interaction_type, database, via)
            for ev in self
        )


    def has_data_model(self, data_model, **kwargs):

        return any(
            ev.has_data_model(data_model)
            for ev in self.filter(**kwargs)
        )


    def get_resources(self, **kwargs):

        return {ev.resource for ev in self}


    def get_resource_names(self, **kwargs):

        return {ev.resource.name for ev in self.filter(**kwargs)}


    def get_interaction_types(self, **kwargs):

        return {ev.resource.interaction_type for ev in self.filter(**kwargs)}


    def get_data_models(self, **kwargs):

        return {ev.resource.data_model for ev in self.filter(**kwargs)}


    def get_datasets(self, **kwargs):

        return {
            ds for ds in
            (ev.dataset for ev in self.filter(**kwargs))
            if ds
        }


    def has_dataset(self, dataset: str, **kwargs) -> bool:
        """
        Contains evidence(s) from a given dataset meeting the criteria.

        Args:
            dataset:
                Name of the dataset.
            kwargs:
                Filtering criteria for evidences.
        """

        return dataset in self.get_datasets(**kwargs)


    def __isub__(self, other):

        if isinstance(other, self.__class__):

            self.evidences = dict(
                (key, ev)
                for key, ev in iteritems(self.evidences)
                if key not in other.evidences or other.evidences[key] != ev
            )

        else:

            self.remove(other)

        return self


    def remove(self, resource = None, interaction_type = None, via = False):

        self.evidences = dict(
            (key, ev)
            for key, ev in iteritems(self.evidences)
            if not ev.match(
                resource = resource,
                interaction_type = interaction_type,
                via = via,
            )
        )


    def filter(
            self,
            resource = None,
            data_model = None,
            interaction_type = None,
            via = False,
            references = None,
            datasets = None,
        ):

        return (
            ev for ev in self
            if ev.match(
                resource = resource,
                data_model = data_model,
                interaction_type = interaction_type,
                via = via,
                references = references,
                datasets = datasets,
            )
        )


    def match(
            self,
            resource = None,
            data_model = None,
            interaction_type = None,
            via = False,
            references = None,
            datasets = None,
        ):

        return bool(
            tuple(
                self.filter(
                    resource = resource,
                    data_model = data_model,
                    interaction_type = interaction_type,
                    via = via,
                    references = references,
                    datasets = datasets,
                )
            )
        )


    def __getitem__(self, key):
        """
        Key is a :py:class:`pypath.internals.resource.NetworkResourceKey` or
        an equivalent tuple.
        """

        return self.evidences.get(key, None) or self.simple_dict[key]


    def keys(self):
        """
        Returns
            (dict_keys): The keys of this dictionary are
                :py:class:`pypath.internals.resource.NetworkResourceKey`
                objects.
        """

        return self.evidences.keys()


    def items(self):
        """
        Returns
            (dict_items): The evidences as a mapping, with
                :py:class:`pypath.internals.resource.NetworkResourceKey`
                objects as keys and :py:class:`pypath.core.evidence.Evidence`
                objects as values.
        """

        return self.evidences.items()


    @property
    def simple_dict(self) -> dict[str, evidence.Evidence]:
        """
        Returns
            Keys are resource labels, values are ``Evidence`` objects.
        """

        return dict(
            (res.last, ev)
            for res, ev in iteritems(self)
        )


    def serialize_attrs(self, top_key_prefix: bool = True) -> str:
        """
        Serialize the extra attributes of the evidences as a JSON string.

        Returns
            A JSON serialized string with the evidences from each resource.
        """

        return attrs_mod.AttributeHandler._serialize(
            self.simple_dict,
            top_key_prefix = top_key_prefix,
            default = lambda obj: obj.serialize(),
        )


    @property
    def datasets(self) -> set:
        """
        Datasets in this evidence set.
        """

        return {ev.dataset for ev in self}


    @property
    def attrs(self) -> dict:
        """
        Combines the custom attributes from all evidences within this set.
        """

        return common.combine_attrs([ev.attrs for ev in self])


    def asdict(self) -> list[dict]:
        """
        Evidence set as a list of dictionaries.
        """

        return [ev.asdict() for ev in self]
