#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from __future__ import annotations

from future.utils import iteritems
from typing import Any, Iterable, Generator

import os
import collections
import hashlib

try:
    import cPickle as pickle
except:
    import pickle


import pypath.inputs as inputs
import pypath.share.common as common
import pypath.share.constants as constants
import pypath.share.cache as cache
import pypath.share.session as session_mod
import pypath.utils.taxonomy as taxonomy


class ResourceBase(session_mod.Logger):
    """
    One resource with all metainformation and data access methods.
    """


    def __init__(self, name: str, **metadata):
        """
        Create a resource description.

        This object holds all information about one resource as described
        in `resources.json`. It is able to create and operate the loader
        classes, hence provide data of the desired type and the desired
        level of processing (raw or processed).

        Args:
            name:
                Name of the resource.
            metadata:
                Arbitrary metadata. Includes the data type specific
                input definitions, license information, etc.
        """

        self.name = name

        for k, v in metadata.items():

            setattr(self, k, v)


    def __str__(self) -> str:

        return self.name


    def __repr__(self) -> str:

        return f'<Resource: {self.name}>'


    def data_types(self) -> tuple[str]:
        """
        Data types available from this resource.
        """

        return tuple(sorted(getattr(self, 'inputs', ())))


    @staticmethod
    def _data_type_pascal_case(data_type: str) -> str:

        return ''.join(i.capitalize() for i in data_type.split('_'))


    def loader(
            self,
            data_type: str,
            organism: str | int | None = None,
            **kwargs
        ) -> callable:
        """
        Loader for a certain data type.
        """

        class_name = f'{self._data_type_pascal_case(data_type)}Resource'
        cls = globals().get(class_name)
        ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)

        if cls:

            return cls(
                name = self.name,
                ncbi_tax_id = ncbi_tax_id,
                **kwargs
            )


    def raw(
            self,
            data_type: str,
            organism: str | int | None = None,
            **kwargs
        ) -> dict | list | Iterable:
        """
        Data after resource specific preprocessing.

        As provided by the ``pypath.inputs`` module: the relevant fields
        extracted, certain identifiers translated, some wrong or undesired
        records removed and organised into straightforward, simple data
        structure.
        """

        return self._loader_op('raw', locals())


    def records(
            self,
            data_type: str,
            organism: str | int | None = None,
            **kwargs
        ) -> Any:
        """
        Data after full preprocessing.

        Records as stored in the database obejcts in ``pypath.core``.
        """

        return self._loader_op('records', locals())


    def _loader_op(self, op: str, loc: dict) -> dict | list | Iterable:
        """
        Execute an operation on a loader for a certain data type.
        """

        loader = self.loader(
            data_type = loc['data_type'],
            organism = loc['organism'],
            **loc['kwargs']
        )

        return getattr(loader, op)()


class AbstractResource(session_mod.Logger):
    """
    Generic class for downloading, processing and serving data from a resource.
    """


    def __init__(
            self,
            name: str,
            organism: str | int = 9606,
            input_method: callable | None = None,
            input_args: dict | None = None,
            data_attr_name: str | None = None,
            data_type: str | None = None,
            evidence_types: set[str] | None = None,
            key: callable | tuple[str] | None = None,
            raw_pickle: str | None = None,
            records_pickle: str | None = None,
            database_pickle: str | None = None,
            **kwargs
        ):
        """
        Generic data resource.

        Args:
            name:
                Custom name for the resource.
            organism:
                Name or NCBI Taxonomy ID of an organism supported by the
                resource.
            input_method:
                Method providing the input data. If not provided, a standard
                name derived from the resource name will be used.
            input_args:
                Arguments for the input function.
            data_attr_name:
                A synonym for the attribute that strores the processed records.
            data_type:
                The database domain this resource belongs to: network,
                complexes, enz_sub, annotations.
            key:
                A function or class suitable to identify this resource.
            raw_pickle:
                Path to a pickle file with raw, minimally preprocessed data.
            records_pickle:
                Path to a pickle file with processed records.
            database_pickle:
                Path to a pickle file with database object.
            kwargs:
                Ignored.
        """

        if not hasattr(self, '_log_name'):

            session_mod.Logger.__init__(self, name = 'resource_{name}')

        self.dump = dump
        self.name = name
        self._data_attr_name = data_attr_name or 'data'
        self._data_type = data_type
        self._input_method = input_method
        self.input_args = input_args or {}
        self.ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)
        self._raw_pickle = raw_pickle
        self._records_pickle = records_pickle
        self._database_pickle = database_pickle
        self.evidence_types = evidence_types or set()
        self._key = key
        self._data = {}


    def raw(self) -> Any:
        """
        Load the data with resource specific preprocessing.
        """

        return self._get_data('raw')


    def records(self):
        """
        Load and process data
        """

        return self._get_data('records')


    def database(self):
        """
        Load and process data
        """

        return self._get_data('database')


    def _get_data(self, stage: str) -> Any:

        for method in (self.load, self.build):

            if stage not in self._data:

                method(stage)

        return self._data.get(stage)


    def set_method(self):
        """
        Sets the data input method by looking up in ``inputs`` module if
        necessary.
        """

        if callable(self._input_method):

            self.input_method = self._input_method

        else:

            self._input_method = (
                self._input_method or
                f'{self.name}.{self.name}_{self._data_type}'.lower()
            )

            self.input_method = inputs.get_method(self._input_method)


    def build(self, stage: str) -> Any:
        """
        Builds resource data at a certain stage of processing.
        """

        self._log(f'Building `{stage}` stage data from `{self.name}`.')

        getattr(self, f'_build_{stage}')()

        if not getattr(self, f'_{stage}_pickle'):

            self.save(stage, self._data[stage])


    def _build_raw(self):

        self.set_method()

        if hasattr(self, 'input_method'):

            self._data['raw'] = self.input_method(**self.input_args))


    def _build_records(self):
        """
        Calls the ``_process_method``.
        """

        self._log('Processing data from `%s`.' % self.name)
        return self._process_method()


    def _process_method(self):

        self._data['records'] = self._data['raw']


    def _build_database(self):

        self._data['database'] = None


    def load(self, stage: str) -> Any:
        """
        Load data from a pickle dump.
        """

        path = self.pickle_path(stage)

        if isinstance(path, common.basestring) and os.path.exists(path):

            with open(self.dump, 'rb') as fp:

                self._data[stage] = pickle.load(fp)


    def save(self, stage: str, data: Any | None = None):
        """
        Save data into a pickle dump.
        """

        path = self.pickle_path(stage)
        data = data or getattr(self, stage)()

        if isinstance(data, Generator):

            data = list(data)

        with open(path, 'wb') as fp:

            self._log(f'Saving stage `{stage}` data to `{path}`.')
            pickle.dump(obj = data, file = fp)


    def pickle_path(self, stage: str) -> str:
        """
        Path to a pickle dump of data from this resource at a certain stage.
        """

        return getattr(self, f'_{stage}_pickle') or self.cache_path(stage)


    def __eq__(self, other: Any) -> bool:

        return (
            self.name == other.name and self.data_type == other.data_type
                if isinstance(other, self.__class__) else
            self.name == other
        )


    @property
    def data_type(self) -> str:

        return ''.join(i.capitalize() for i in self._data_type.split('_'))


    def __len__(self) -> int:

        return len(getattr(self, self._data_attr_name, ()))


    def __str__(self) -> str:

        return self.name


    def __repr__(self) -> str:

        return f'<{self.data_type}Resource: {self.name} ({len(self)} records)>'


    @property
    def key_class_name(self) -> str:
        """
        Name of the class that creates unique keys to identify a resource.
        """

        return (
            self._key.__name__
                if isinstance(self._key, callable) else
            self._key
                if isinstance(self._key, common.basestring) else
            f'{self.data_type}ResourceKey'
        )


    @property
    def key_class(self) -> callable:
        """
        The class that creates unique keys to identify a resource.
        """

        if isinstance(self._key, callable):

            return self._key

        name = self.key_class_name

        the_class = globals().get(name) or self._create_key_class()

        return the_class


    def _create_key_class(self) -> callable:

        name = self.key_class_name

        fields = (
            self._key
                if isinstance(self._key, (tuple, list)) else
            ('name', 'data_type', 'via')
        )


        class the_class(collections.namedtuple(f'{name}Base', fields)):


            def __new__(cls, *args, **kwargs):

                return super(the_class, cls).__new__(cls, *args, **kwargs)


            @property
            def label(self) -> str:
                """
                Returns:
                    A label containing the resource name, and if it's a
                    secondary resource, the name of the primary resource
                    separated by an underscore.
                """

                return f'{self.name}_{self.via}' if self.via else self.name


            @property
            def last(self) -> str:
                """
                Returns:
                    The name of the resource where the data directly came
                    from ignoring the primary resource.
                """

                return self.via or self.name


        the_class.__name__ = name
        globals()[name] = the_class

        return the_class


    @property
    def key(self) -> callable:
        """
        Key for identifying this resource.
        """

        cls = self.key_class

        args = {k: getattr(self, k) for k in cls._fields}

        return cls(**args)


    def cache_key(self) -> str:
        """
        Unique identifier of the resource with a set of parameters.
        """

        return (
            hashlib.md5(
                (
                    common.dict_str(self.key.asdict(), sep = ',') +
                    f',{self.ncbi_tax_id},' +
                    common.dict_str(self.input_args, sep = ',')
                ).
                encode('utf8')
            ).
            hexdigest()
        )


    def cache_path(self, stage: str) -> str:
        """
        Path to the cached pickle dump of resource data at a certain stage.
        """

        return (
            os.path.join(
                cache.get_cachedir(),
                f'{self.cache_key()}-{stage}.pickle'
            )
        )


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
            **kwargs
        ):

        ResourceAttributes.__init__(
            self,
            name = name,
            data_type = 'network',
            interaction_type = interaction_type,
            evidence_types = evidence_types,
            data_model = data_model,
            via = via,
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
                if isinstance(other, common.basestring) else
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
                if isinstance(other, common.basestring) else
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
