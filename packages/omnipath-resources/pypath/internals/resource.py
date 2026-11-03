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

"""
Generic objects for representing resources.
"""

from future.utils import iteritems

from typing import Iterable, Mapping, TYPE_CHECKING

if TYPE_CHECKING:

    import pypath.internals.license as License

import os
import collections
import copy

try:
    import cPickle as pickle
except:
    import pickle


import pypath.inputs as inputs
import pypath.share.common as common
import pypath.share.session as session_mod


class AbstractResource(session_mod.Logger):
    """
    Generic class for downloading, processing and serving
    data from a resource.
    """


    def __init__(
            self,
            name,
            ncbi_tax_id = 9606,
            input_method = None,
            input_args = None,
            dump = None,
            data_attr_name = None,
            **kwargs
        ):
        """
        name : str
            Custom name for the resource.
        input_method : callable
            Method providing the input data.
        """

        if not hasattr(self, '_log_name'):

            session_mod.Logger.__init__(self, name = 'resource')

        self.dump = dump
        self.name = name
        self._data_attr_name = data_attr_name or 'data'
        self._input_method = input_method
        self.input_args = input_args or {}
        self.ncbi_tax_id = ncbi_tax_id


    def load(self):

        self.set_method()
        from_dump = self.from_dump()

        if not from_dump:

            self.load_data()
            self.process()

        if hasattr(self, 'data'):

            delattr(self, 'data')


    def set_method(self):
        """
        Sets the data input method by looking up in ``inputs`` module if
        necessary.
        """

        if callable(self._input_method):

            self.input_method = self._input_method

        elif self._input_method:

            self.input_method = inputs.get_method(self._input_method)


    def load_data(self):
        """
        Loads the data by calling ``input_method``.
        """

        self._log('Loading data from `%s`.' % self.name)

        self.set_method()

        if hasattr(self, 'input_method'):

            self.data = self.input_method(**self.input_args)


    def process(self):
        """
        Calls the ``_process_method``.
        """

        self._log('Processing data from `%s`.' % self.name)
        self._process_method()


    def _process_method(self):

        pass


    def from_dump(self):

        if self.dump is not None:

            if (
                isinstance(self.dump, str) and
                os.path.exists(self.dump)
            ):

                with open(self.dump, 'rb') as fp:

                    self._from_dump = pickle.load(fp)

            else:

                self._from_dump = self.dump

            self._from_dump_callback()

            return True

        return False


    def _from_dump_callback(self):

        if hasattr(self, '_from_dump'):

            setattr(self, self._data_attr_name, self._from_dump)
            delattr(self, '_from_dump')
            delattr(self, 'dump')


    def save_to_pickle(self, pickle_file):

        with open(pickle_file, 'wb') as fp:

            pickle.dump(
                obj = getattr(self, self._data_attr_name),
                file = fp,
            )


class ResourceAttributes(object):


    def __init__(
            self,
            name,
            data_type,
            evidence_types = None,
            dataset = None,
            **kwargs
        ):

        self.name = name
        self.data_type = data_type
        self.evidence_types = evidence_types or set()
        self.resource_attrs = {}

        for attr, value in iteritems(kwargs):

            setattr(self, attr, value)

        self.dataset = dataset


    def __eq__(self, other):

        return (
            self.name == other.name and self.data_type == other.data_type
                if isinstance(other, self.__class__) else
            self.name == other
        )


    def __str__(self):

        return self.name


    @property
    def dataset(self):

        return self._dataset


    @dataset.setter
    def dataset(self, dataset):

        self._dataset = dataset

        networkinput = getattr(self, 'networkinput', None)

        if hasattr(self, 'networkinput'):

            netinput_new = copy.deepcopy(networkinput)
            netinput_new.dataset = dataset
            self.networkinput = netinput_new


class NetworkResourceKey(
        collections.namedtuple(
            'NetworkResourceKeyBase',
            [
                'name',
                'data_type',
                'interaction_type',
                'data_model',
                'via',
            ]
        )
    ):


    def __new__(cls, *args, **kwargs):

        return super(NetworkResourceKey, cls).__new__(cls, *args, **kwargs)


    @property
    def label(self):
        """
        Returns
            (str): A label containing the resource name, and if it's a
                secondary resource, the name of the primary resource
                separated by an underscore.
        """

        return '%s_%s' % (self.name, self.via) if self.via else self.name


    @property
    def last(self):
        """
        Returns
            (str): The name of the resource where the data directly came from
                ignoring the primary resource.
        """

        return self.via or self.name


class NetworkResource(ResourceAttributes):


    _key = NetworkResourceKey


    def __init__(
            self,
            name,
            interaction_type = 'PPI',
            data_model = None,
            evidence_types = None,
            via = None,
            dataset = None,
            **kwargs
        ):

        if not dataset and 'networkinput' in kwargs:

            dataset = kwargs['networkinput'].dataset

        ResourceAttributes.__init__(
            self,
            name = name,
            data_type = 'network',
            interaction_type = interaction_type,
            evidence_types = evidence_types,
            data_model = data_model,
            via = via,
            dataset = dataset,
            **kwargs
        )


    def __hash__(self):

        return hash(self.key)


    @property
    def key(self):

        return self._key(
            name = self.name,
            data_type = self.data_type,
            interaction_type = self.interaction_type,
            data_model = self.data_model,
            via = self.via,
        )


    def __eq__(self, other):

        return (
            self.name == other
                if isinstance(other, str) else
            self.__hash__() == other.__hash__()
        )


    def __repr__(self):

        return '<NetworkResource: %s (%s, %s)>' % (
            self.name,
            self.interaction_type,
            self.data_model,
        )


    def is_primary(self):

        return self.via is None


    @property
    def data_model_label(self):

        return (
            self.data_model.capitalize().replace('_', ' ')
                if self.data_model else
            'Unknown'
        )

    @property
    def license(self) -> license.License | None:

        return self.resource_attrs.get('license', None)


class NetworkDataset(collections.abc.MutableMapping):


    def __init__(
            self,
            name: str,
            resources: dict | list | None = None,
        ):
        """
        A set of network resources.

        Formerly the network datasets were represented by dicts. This is
        only a thin wrapper around that solution to better organise metadata
        of the datasets and resources within.
        """

        self._name = name
        self._resources = {}
        self.add(resources)


    def __repr__(self):

        it = ', '.join(self.interaction_types)

        return f'<NetworkDataset: {self.name} ({len(self)} resources; {it})>'


    def __iter__(self):

        return (r for r in self._resources.values())


    def __len__(self):

        return len(self._resources)


    @property
    def interaction_types(self):

        return sorted({r.interaction_type for r in self})


    def __setitem__(self, key, value):

        self.add(value, key)


    def __getitem__(self, key):

        return self._resources[key]


    def __delitem__(self, key):

        del self._resources[key]


    def __contains__(self, key):

        return (
            key in self._resources or
            any(r.name == key for r in self._resources.values())
        )


    def __eq__(self, other):

        return (
            self._name == other or
            self._name == getattr(other, '_name', None)
        )


    def items(self):

        return self._resources.items()


    def values(self):

        return self._resources.values()


    def keys(self):

        return self._resources.keys()


    def add(self, value, key = None):

        if isinstance(value, Mapping):

            for label, resource in value.items():

                self.add(resource, label)

        elif isinstance(value, Iterable):

            for resource in value:

                self.add(resource)

        elif isinstance(value, NetworkResource):

            resource = copy.deepcopy(value)
            resource.dataset = self.name
            self._resources[key or resource.name] = resource


    update = add


    @property
    def name(self):

        return self._name


    @name.setter
    def name(self, name):

        for resource in self.values():

            resource.networkinput.dataset = name

        self._name = name


    def __copy__(self):

        return self.__class__(
            name = self.name,
            resources = self._resources,
        )


    def rename(self, name: str):

        new = self.__class__(name = name)
        new.add(self)

        return new


    def remove(self, remove: str | set | None):

        remove = common.to_set(remove)

        self._resources = {
            k: v for k, v in self.items()
            if k not in remove and v.name not in remove
        }


    def without(self, exclude: str | set | None):

        new = copy.copy(self)
        new.remove(exclude)

        return new


EnzymeSubstrateResourceKey = collections.namedtuple(
    'EnzymeSubstrateResourceKey',
    [
        'name',
        'data_type',
        'via',
    ]
)


class EnzymeSubstrateResource(ResourceAttributes):


    _key = EnzymeSubstrateResourceKey


    def __init__(
            self,
            name,
            input_method,
            evidence_types = None,
            via = None,
            id_type_enzyme = 'uniprot',
            id_type_substrate = 'uniprot',
            organisms_supported = False,
            organisms = None,
            resource_attrs = None,
            extra_attrs = None,
            **kwargs
        ):

        ResourceAttributes.__init__(
            self,
            name = name,
            input_method = input_method,
            data_type = 'enzyme_substrate',
            evidence_types = evidence_types,
            via = via,
            id_type_enzyme = id_type_enzyme,
            id_type_substrate = id_type_substrate,
            organisms_supported = organisms_supported,
            organisms = organisms,
            resource_attrs = resource_attrs or {},
            extra_attrs = extra_attrs or {},
            **kwargs
        )


    def __hash__(self):

        return hash(self.key)


    @property
    def key(self):

        return self._key(
            name = self.name,
            data_type = self.data_type,
            via = self.via,
        )


    def __eq__(self, other):

        return (
            self.name == other
                if isinstance(other, str) else
            self.__hash__() == other.__hash__()
        )


    def __repr__(self):

        return '<EnzymeSubstrateResource: %s>' % (
            self.name,
        )


    def is_primary(self):

        return self.via is None


    def get_via(self, name):
        """
        Returns a copy of the same resource attributes but the ``name`` set
        to ``name`` and the ``via`` set to the original name. This means
        the data comes from the resource ``name`` via the resource ``via``.
        """

        args = dict(
            (k, getattr(self, k))
            for k in self.__dir__()
            if (
                not k.startswith('__') and
                not callable(getattr(self, k))
            )
        )
        args['via'] = self.name
        args['name'] = name
        _ = args.pop('data_type', None)
        _ = args.pop('key', None)

        return EnzymeSubstrateResource(**args)
