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

import collections

import pypath.share.session as _session

"""
Classes for reading and use serving ID mapping data from custom file,
function, UniProt, UniProt ID Mapping, Ensembl BioMart,
PRO (Protein Ontology), miRBase or pickle file.
"""

MappingTableKey = collections.namedtuple(
    'MappingTableKey',
    [
        'id_type',
        'target_id_type',
        'ncbi_tax_id',
    ],
)
MappingTableKey.__new__.__defaults__ = ('protein', 9606)


class MappingTable(_session.Logger):
    """
    This is the class directly handling ID translation data.
    It does not care about loading it or what kind of IDs these
    only accepts the translation dictionary.

    lifetime : int
        If this table has not been used for longer than this preiod it is
        to be removed at next cleanup. Time in seconds.
    """

    def __init__(
            self,
            data,
            id_type,
            target_id_type,
            ncbi_tax_id,
            lifetime = 300,
        ):
        """
        Wrapper around a dictionary of identifier mapping. The dictionary
        is located in the `data` attribute, keys are the source identifiers,
        values are sets of target identifiers. Most often the mapping is
        unambigous, which means one target identifier for each source
        identifier.

        Args
            data (dict): The identifier translation dictionary.
            id_type (str): The source ID type.
            target_id_type (str): The target ID type.
            ncbi_tax_id (int): NCBI Taxonomy identifier of the organism.
            lifetime (int): Time in seconds to keep the table loaded in
                the memory. If not used, the table will be unloaded after
                this time. Each usage resets the expiry time.
        """

        _session.Logger.__init__(self, name = 'mapping')

        self.id_type = id_type
        self.target_id_type = target_id_type
        self.ncbi_tax_id = ncbi_tax_id
        self.data = data
        self.lifetime = lifetime
        self._used()


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def __getitem__(self, key):

        self._used()

        if key in self.data:

            return self.data[key]

        return set()


    def __contains__(self, key):

        self._used()

        return key in self.data


    def __len__(self):

        return len(self.data)


    def _used(self):

        self._last_used = time.time()


    def _expired(self):

        return time.time() - self._last_used > self.lifetime


    def get_key(self):
        """
        Creates a mapping table key, a tuple with all the defining properties
        of the mapping table.
        """

        return MappingTableKey(
            id_type = self.id_type,
            target_id_type = self.target_id_type,
            ncbi_tax_id = self.ncbi_tax_id,
        )


    @property
    def key(self):

        return MappingTableKey(
            id_type = self.id_type,
            target_id_type = self.target_id_type,
            ncbi_tax_id = self.ncbi_tax_id,
        )


    def __repr__(self):

        return '<MappingTable from=%s, to=%s, taxon=%u (%u IDs)>' % (
            self.key + (len(self),)
        )


    @property
    def items(self):

        return self.data.items


    @property
    def keys(self):

        return self.data.keys


    @property
    def values(self):

        return self.data.values
