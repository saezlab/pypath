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

from future.utils import iteritems
from past.builtins import xrange, range

import os
import sys
import math
import re
import importlib as imp
import collections
import functools
import datetime
import time

import urllib

if not hasattr(urllib, 'urlencode'):

    import urllib.parse
    _urllib = urllib
    urllib = _urllib.parse

import json
try:
    import cPickle as pickle
except:
    import pickle

from typing import Iterable, List, Literal, Optional, Set, Union

import pandas as pd

# from pypath:
import pypath.share.progress as progress
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.share.cache as cache_mod
import pypath.internals.maps as maps
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs as inputs
import pypath.inputs.uniprot as uniprot_input
import pypath.inputs.uniprot_db as uniprot_db
import pypath.inputs.pro as pro_input
import pypath.inputs.biomart as biomart_input
import pypath.inputs.unichem as unichem_input
import pypath.inputs.ramp as ramp_input
import pypath.inputs.hmdb as hmdb_input
import pypath.internals.input_formats as input_formats
import pypath.utils.reflists as reflists
import pypath.utils.taxonomy as taxonomy
import pypath.share.settings as settings
import pypath.share.session as session_mod
_logger = session_mod.log()

import pypath.utils._mapping._cleanup as _cleanup


__all__ = ['MapReader', 'MappingTable', 'Mapper']

_logger = session_mod.Logger(name = 'mapping')
_log = _logger._log

try:
    UNICHEM_NAME_TYPES = set(unichem_input.unichem_sources().values())
except Exception as e:
    exc = sys.exc_info()
    _log('Failed to retrieve UniChem ID types:')
    _logger._log_traceback()
    UNICHEM_NAME_TYPES = ()

RESOURCES_EXPLICIT = ('uniprot', 'basic', 'mirbase', 'ipi')

RESOURCES_IMPLICIT = (
    (
        input_formats.AC_MAPPING,
        'uniprot',
        input_formats.UniprotListMapping,
    ),
    (
        input_formats.PRO_MAPPING,
        'pro',
        input_formats.ProMapping,
    ),
    (
        input_formats.BIOMART_MAPPING,
        'biomart',
        input_formats.BiomartMapping,
    ),
    (
        input_formats.ARRAY_MAPPING,
        'array',
        input_formats.ArrayMapping,
    ),
    (
        {n: n for n in UNICHEM_NAME_TYPES},
        'unichem',
        input_formats.UnichemMapping,
    ),
    (
        dict(
            **{
                it: it
                for it in ramp_input.ramp_id_types_2('compound')
            },
            **input_formats.RAMP_MAPPING,
        ),
        'ramp',
        input_formats.RampMapping,
    ),
    (
        dict(
            **{
                it: it
                for it in hmdb_input.ID_FIELDS
            },
            **input_formats.HMDB_MAPPING,
        ),
        'hmdb',
        input_formats.HmdbMapping,
    ),
)

UNIPROT_ID_TYPES = {
    'uniprot',
    'trembl',
    'swissprot',
    'uniprot-pri',
    'uniprot-sec',
}


class Mapper(session_mod.Logger):

    default_name_types = settings.get('default_name_types')
    default_label_types = settings.get('default_label_types')

    def _get_label_type_to_id_type(default_name_types):

        label_type_to_id_type = dict(
            (
                label_type,
                default_name_types[entity_type],
            )
            for entity_type, label_type in
            iteritems(settings.get('default_label_types'))
        )
        #TODO: some nicer solution
        label_type_to_id_type['mir-name'] = 'mir-pre'

        return label_type_to_id_type

    label_type_to_id_type = _get_label_type_to_id_type(default_name_types)


    def __init__(
            self,
            ncbi_tax_id = None,
            lifetime = 300,
            translate_deleted_uniprot = None,
            keep_invalid_uniprot = None,
            trembl_swissprot_by_genesymbol = None,
        ):
        """
        lifetime : int
            If a table has not been used for longer than this preiod it is
            to be removed at next cleanup.
        translate_deleted_uniprot : bool
            Do an extra attempt to translate deleted or obsolete UniProt IDs
            by retrieving their archived datasheet and use the gene symbol
            to find the corresponding valid UniProt ID?
        keep_invalid_uniprot : bool
            If the target ID is UniProt, keep the results if they fit the
            format for UniProt IDs (we won't check if they are deleted or
            from a different taxon). The alternative is to keep only those
            which are in the list of all UniProt IDs for the given organism.
        trembl_swissprot_by_genesymbol : bool
            Attempt to translate TrEMBL IDs to SwissProt by translating to
            gene symbols and then to SwissProt.
        """

        session_mod.Logger.__init__(self, name = 'mapping')

        self._translate_deleted_uniprot = settings.get(
            'mapper_translate_deleted_uniprot',
            translate_deleted_uniprot,
        )
        self._keep_invalid_uniprot = settings.get(
            'mapper_keep_invalid_uniprot',
            keep_invalid_uniprot,
        )
        self._trembl_swissprot_by_genesymbol = settings.get(
            'mapper_trembl_swissprot_by_genesymbol',
            trembl_swissprot_by_genesymbol,
        )

        # regex for matching UniProt AC format
        self.reuniprot = re.compile(r'^(?:%s)$' % uniprot_input.reac.pattern)
        self.remipreac = re.compile(r'^MI\d{7}$')
        self.remimatac = re.compile(r'^MIMAT\d{7}$')
        self.remipreid = re.compile(
            r'^[a-z]{3}-'
            r'(?:mir|MIR|let|lsy|lin)-?'
            r'\d+-?[A-z\*]*(?:-((?!p)[\w\*\.-])+)?$'
        )
        self.remimatid = re.compile(
            r'^[a-z]{3}-'
            r'(?:miR|let|lsy|lin)-?'
            r'\d+[a-z\*]*(?:-((?!p)[\w\*])+)?(?:-(3|5)p)?$'
        )
        self.cachedir = cache_mod.get_cachedir()
        self.ncbi_tax_id = ncbi_tax_id or settings.get('default_organism')

        self.unmapped = []
        self.tables = {}
        self.uniprot_mapped = []
        self.trace = []
        self.uniprot_static_names = {
            'uniprot_id': 'UniProtKB-ID',
            'embl': 'EMBL-CDS',
            'embl_id': 'EMBL',
            'entrez': 'GeneID',
            'gi': 'GI',
            'refseqp': 'RefSeq',
            'refseqn': 'RefSeq_NT',
            'ensembl': 'Ensembl',
            'ensg': 'ENSEMBL',
            'ensp': 'ENSEMBL_PRO_ID',
            'enst': 'ENSEMBL_TRS',
            'hgnc': 'HGNC',
        }
        self.names_uniprot_static = (
            common.swap_dict_simple(self.uniprot_static_names)
        )

        _cleanup.register(self)


    def reload(self):
        """
        Reload the class from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def get_table_key(
            self,
            id_type,
            target_id_type,
            ncbi_tax_id = None,
        ):
        """
        Returns a tuple unambigously identifying a mapping table.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        return MappingTableKey(
            id_type = id_type,
            target_id_type = target_id_type,
            ncbi_tax_id = ncbi_tax_id,
        )


    def which_table(
            self,
            id_type,
            target_id_type,
            load = True,
            ncbi_tax_id = None,
        ):
        """
        Returns the table which is suitable to convert an ID of
        id_type to target_id_type. If no such table have been loaded
        yet, it attempts to load from UniProt. If all attempts failed
        returns `None`.
        """

        tbl = None
        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id


        def check_loaded():

            return self.which_table(
                id_type = id_type,
                target_id_type = target_id_type,
                load = False,
                ncbi_tax_id = ncbi_tax_id,
            )


        tbl_key = self.get_table_key(
            id_type = id_type,
            target_id_type = target_id_type,
            ncbi_tax_id = ncbi_tax_id,
        )
        tbl_key_noorganism = self.get_table_key(
            *tbl_key[:-1],
            ncbi_tax_id = _const.NOT_ORGANISM_SPECIFIC,
        )

        tbl_key_rev = self.get_table_key(
            id_type = target_id_type,
            target_id_type = id_type,
            ncbi_tax_id = ncbi_tax_id,
        )
        tbl_key_rev_noorganism = self.get_table_key(
            *tbl_key_rev[:-1],
            ncbi_tax_id = _const.NOT_ORGANISM_SPECIFIC,
        )

        if tbl_key in self.tables:

            tbl = self.tables[tbl_key]

        elif tbl_key_noorganism in self.tables:

            tbl = self.tables[tbl_key_noorganism]

        elif tbl_key_rev in self.tables:

            self.create_reverse(tbl_key_rev)
            tbl = self.tables[tbl_key_rev]

        elif tbl_key_rev_noorganism in self.tables:

            self.create_reverse(tbl_key_rev_noorganism)
            tbl = self.tables[tbl_key_rev_noorganism]

        elif load:

            self._log(
                'Requested to load ID translation table from '
                '`%s` to `%s`, organism: %u.' % (
                    id_type,
                    target_id_type,
                    ncbi_tax_id,
                )
            )

            if id_type == 'complex' or target_id_type == 'complex':

                raise ValueError('Can not translate protein complexes.')

            id_types = (id_type, target_id_type)
            id_types_rev = tuple(reversed(id_types))
            resource = None

            for resource_attr in RESOURCES_EXPLICIT:

                resources = getattr(maps, resource_attr)

                if id_types in resources:

                    resource = resources[id_types]
                    load_a_to_b = True
                    load_b_to_a = False

                elif id_types_rev in resources:

                    resource = resources[id_types_rev]
                    load_a_to_b = False
                    load_b_to_a = True

                if resource:

                    self._log(
                        'Chosen built-in defined ID translation table: '
                        'resource=%s, id_type_a=%s, id_type_b=%s' % (
                            resource_attr,
                            resource.id_type_a,
                            resource.id_type_b,
                        )
                    )

                    self.load_mapping(
                        resource = resource,
                        load_a_to_b = load_a_to_b,
                        load_b_to_a = load_b_to_a,
                        ncbi_tax_id = ncbi_tax_id,
                    )

                    tbl = check_loaded()

                    break

                if tbl is not None:

                    break

            if tbl is None:

                basic_services = {'hmdb', 'ramp', 'uniprot', 'unichem'}

                for (service_ids, service_id_type, input_cls) in (
                    RESOURCES_IMPLICIT
                ):

                    if (
                        (
                            input_cls.possible(
                                id_type,
                                target_id_type,
                                ncbi_tax_id,
                            ) and
                            id_type != target_id_type
                        ) or (
                            service_id_type == 'pro' and (
                                (
                                    id_type in service_ids or
                                    target_id_type in service_ids
                                ) and
                                (
                                    id_type == service_id_type or
                                    target_id_type == service_id_type
                                )
                            )
                        ) or (
                            service_id_type == 'biomart' and (
                                (
                                    id_type in service_ids and
                                    target_id_type in service_ids
                                )
                            )
                        ) or (
                            service_id_type == 'array' and (
                                (
                                    id_type in service_ids and
                                    target_id_type in {'ensg', 'enst', 'ensp'}
                                ) or
                                (
                                    target_id_type in service_ids and
                                    id_type in {'ensg', 'enst', 'ensp'}
                                )
                            )
                        )
                    ):

                        if target_id_type == service_id_type:

                            _id_type, _target_id_type = (
                                target_id_type,
                                id_type,
                            )
                            load_a_to_b = False
                            load_b_to_a = True

                        else:

                            _id_type, _target_id_type = (
                                id_type,
                                target_id_type,
                            )
                            load_a_to_b = True
                            load_b_to_a = False

                        self._log(
                            'Chosen ID translation table from service: '
                            'service=%s, id_type_a=%s, id_type_b=%s' % (
                                service_id_type,
                                _id_type,
                                _target_id_type,
                            )
                        )

                        if service_id_type in {'hmdb', 'ramp', 'unichem'}:

                            ncbi_tax_id = _const.NOT_ORGANISM_SPECIFIC
                            tbl_key = tbl_key_noorganism
                            tbl_key_rev = tbl_key_rev_noorganism

                        # for uniprot/idmapping or PRO or array
                        # we create here the mapping params
                        this_param = input_cls(
                            id_type_a = _id_type,
                            id_type_b = _target_id_type,
                            ncbi_tax_id = ncbi_tax_id,
                        )

                        reader = MapReader(
                            param = this_param,
                            ncbi_tax_id = ncbi_tax_id,
                            load_a_to_b = load_a_to_b,
                            load_b_to_a = load_b_to_a,
                            uniprots = None,
                            lifetime = 300,
                            resource_id_types = service_ids,
                        )

                        self.tables[tbl_key] = getattr(
                            reader,
                            'mapping_table_%s_to_%s' % (
                                reader.id_type_side(tbl_key.id_type),
                                reader.id_type_side(tbl_key.target_id_type),
                            )
                        )

                    tbl = check_loaded()

                    if tbl:

                        break

            if tbl is None and id_type == 'genesymbol5':

                self.load_genesymbol5(ncbi_tax_id = ncbi_tax_id)

                tbl = check_loaded()

            if tbl is None:

                if (
                    settings.get('mapping_uniprot_static') and
                    id_type in self.uniprot_static_names and
                    target_id_type == 'uniprot'
                ):

                    self.load_uniprot_static([id_type])

                    tbl = check_loaded()

        if tbl is None:

            self._log(
                'Could not find suitable ID translation table '
                f'between id types `{id_type}` and `{target_id_type}` '
                f'for organism `{ncbi_tax_id}`.'
            )

        if hasattr(tbl, '_used'):

            tbl._used()

        return tbl


    @staticmethod
    def reverse_mapping(mapping_table):
        """
        Creates an opposite direction `MappingTable` by swapping the
        dictionary inside an existing `MappingTable` object.

        Args
            mapping_table (MappingTable): A `MappingTable` object.

        Returns
            A new `MappingTable` object.
        """

        rev_data = common.swap_dict(mapping_table.data)

        return MappingTable(
            data = rev_data,
            id_type = mapping_table.target_id_type,
            target_id_type = mapping_table.id_type,
            ncbi_tax_id = mapping_table.ncbi_tax_id,
            lifetime = mapping_table.lifetime,
        )


    def reverse_key(self, key):
        """
        For a mapping table key returns a new key with the identifiers
        reversed.

        Args
            key (tuple): A mapping table key.

        Returns
            A tuple representing a mapping table key, identifiers swapped.
        """

        return (
            self.get_table_key(
                id_type = key.target_id_type,
                target_id_type = key.id_type,
                ncbi_tax_id = key.ncbi_tax_id,
            )
        )


    def create_reverse(self, key):
        """
        Creates a mapping table with ``id_type`` and ``target_id_type``
        (i.e. direction of the ID translation) swapped.
        """

        table = self.tables[key]
        rev_key = self.reverse_key(key)

        self.tables[rev_key] = self.reverse_mapping(table)


    def map_name0(
            self,
            name,
            id_type = None,
            target_id_type = None,
            ncbi_tax_id = None,
            strict = False,
            expand_complexes = None,
            uniprot_cleanup = None,
        ):
        """
        Translates the name and returns only one of the resulted IDs. It
        means in case of ambiguous ID translation, a random one of them
        will be picked and returned. Recommended to use only if the
        translation between the given ID types is mostly unambigous and
        the loss of information can be ignored. See more details at
        `map_name`.
        """

        names = self.map_name(
            name = name,
            id_type = id_type,
            target_id_type = target_id_type,
            ncbi_tax_id = ncbi_tax_id,
            strict = strict,
            expand_complexes = expand_complexes,
            uniprot_cleanup = uniprot_cleanup,
        )

        return list(names)[0] if names else None


    @common.ignore_unhashable
    @functools.lru_cache(maxsize = int(1e5))
    def map_name(
            self,
            name,
            id_type = None,
            target_id_type = None,
            ncbi_tax_id = None,
            strict = False,
            expand_complexes = True,
            uniprot_cleanup = True,
        ):
        """
        Translates one instance of one ID type to a different one.
        Returns set of the target ID type.

        This function should be used to convert individual IDs.
        It takes care about everything and ideally you don't need to
        think on the details.

        How does it work: looks up dictionaries between the original
        and target ID type, if doesn't find, attempts to load from the
        predefined inputs.
        If the original name is genesymbol, first it looks up among the
        preferred gene names from UniProt, if not found, it takes an
        attempt with the alternative gene names. If the gene symbol
        still couldn't be found, and strict = False, the last attempt
        only the first 5 characters of the gene symbol matched. If the
        target name type is uniprot, then it converts all the ACs to
        primary. Then, for the Trembl IDs it looks up the preferred gene
        names, and find Swissprot IDs with the same preferred gene name.

        Args
            name (str): The original name to be converted.
            id_type (str): The type of the name. Available by default:
                - genesymbol (gene name)
                - entrez (Entrez Gene ID \[#\])
                - refseqp (NCBI RefSeq Protein ID \[NP\_\*|XP\_\*\])
                - ensp (Ensembl protein ID \[ENSP\*\])
                - enst (Ensembl transcript ID \[ENST\*\])
                - ensg (Ensembl genomic DNA ID \[ENSG\*\])
                - hgnc (HGNC ID \[HGNC:#\])
                - gi (GI number \[#\])
                - embl (DDBJ/EMBL/GeneBank CDS accession)
                - embl_id (DDBJ/EMBL/GeneBank accession)
                And many more, see the code of
                ``pypath.internals.input_formats``
            target_id_type (str): The name type to translate to, more or
                less the same values are available as for ``id_type``.
            ncbi_tax_id (int): NCBI Taxonomy ID of the organism.
            strict (bool): In case a Gene Symbol can not be translated,
                try to add number "1" to the end, or try to match only
                its first five characters. This option is rarely used,
                but it makes possible to translate some non-standard
                gene names typically found in old, unmaintained resources.
            expand_complexes (bool): When encountering complexes,
                translated the IDs of its components and return a set
                of IDs. The alternative behaviour is to return the
                `Complex` objects.
            uniprot_cleanup (bool): When the `target_id_type` is UniProt
                ID, call the `uniprot_cleanup` function at the end.
        """

        if not name:

            return set()

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        # we support translating from more name types
        # at the same time
        if isinstance(id_type, (list, set, tuple)):

            return set.union(
                *(
                    self.map_name(
                        name = name,
                        id_type = this_id_type,
                        target_id_type = target_id_type,
                        strict = strict,
                        ncbi_tax_id = ncbi_tax_id,
                    )
                    for this_id_type in id_type
                )
            )

        # complexes
        if hasattr(name, 'components'):

            if expand_complexes:

               return set(name.components.keys())

            else:

                return {name}

        # translating from an ID type to the same ID type?
        elif id_type == target_id_type:

            if target_id_type != 'uniprot' or not uniprot_cleanup:

                # no need for translation
                return {name}

            else:

                # we still try to search the primary UniProt
                mapped_names = {name}

        # actual translation comes here
        elif id_type.startswith('refseq'):

            # RefSeq is special
            mapped_names = self._map_refseq(
                refseq = name,
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
                strict = strict,
            )

        elif id_type == 'ensp':

            mapped_names = self._map_ensp(
                ensp = name,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        elif target_id_type == 'ensp':

            mapped_names = self._map_to_ensp(
                name = name,
                id_type = id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        elif (
            (
                id_type in input_formats.ARRAY_MAPPING and
                not target_id_type.startswith('ens')
            ) or (
                target_id_type in input_formats.ARRAY_MAPPING and
                not id_type.startswith('ens')
            )
        ):

            # microarray probe IDs we are able to directly translate
            # only to and from Ensembl gene, transcript and protein IDs
            # if the other ID is different (such as uniprot), we translate
            # in two steps, via Ensembl peptide ID:
            mapped_names = self.chain_map(
                name = name,
                id_type = id_type,
                by_id_type = 'ensp',
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
                strict = strict,
                expand_complexes = expand_complexes,
                uniprot_cleanup = uniprot_cleanup,
            )

        else:

            # all the other ID types
            mapped_names = self._map_name(
                name = name,
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        # as ID translation tables for PRO IDs are not organism specific
        # we need an extra step to limit the results to the target organism
        if id_type == 'pro' and target_id_type == 'uniprot':

            mapped_names = (
                mapped_names &
                reflists.get_reflist(
                    id_type = 'uniprot',
                    ncbi_tax_id = ncbi_tax_id,
                )
            )

        # by default the uniprot-genesymbol tables contain only SwissProt
        if id_type == 'uniprot' and target_id_type == 'genesymbol':

            mapped_names = self._map_name(
                name = name,
                id_type = 'trembl',
                target_id_type = 'genesymbol',
                ncbi_tax_id = ncbi_tax_id,
            )

            if not mapped_names:

                uniprots = self._map_name(
                    name = name,
                    id_type = 'uniprot-sec',
                    target_id_type = 'uniprot-pri',
                    ncbi_tax_id = ncbi_tax_id,
                )

                if uniprots:

                    mapped_names = self.map_names(
                        names = uniprots,
                        id_type = 'uniprot',
                        target_id_type = 'genesymbol',
                        ncbi_tax_id = ncbi_tax_id,
                    )

        # further attempts to set it right if
        # first attempt was not successful

        # for miRNAs if the translation from mature miRNA name failed
        # we still try if maybe it is a hairpin name
        # or the other way around
        if not mapped_names and id_type in {'mir-mat-name', 'mir-name'}:

            for id_type0, id_type1, target_id_type0, target_id_type1 in (
                ('mir-name', 'mir-mat-name', 'mir-pre', 'mirbase'),
                ('mir-mat-name', 'mir-name', 'mirbase', 'mir-pre'),
            ):

                if id_type == id_type0:

                    mapped_names = self._map_name(
                        name = name,
                        id_type = id_type1,
                        target_id_type = target_id_type1,
                        ncbi_tax_id = ncbi_tax_id,
                    )

                    if mapped_names and target_id_type == target_id_type0:

                        mapped_names = self.map_names(
                            names = mapped_names,
                            id_type = target_id_type1,
                            target_id_type = target_id_type0,
                            ncbi_tax_id = ncbi_tax_id,
                        )

                    if mapped_names:

                        break

        # for genesymbol, we automatically try 2 steps mapping via uniprot
        if (
            not mapped_names and (
                id_type == 'genesymbol' or
                target_id_type == 'genesymbol'
            ) and
            id_type not in UNIPROT_ID_TYPES and
            target_id_type not in UNIPROT_ID_TYPES
        ):

            mapped_names = self.chain_map(
                name = name,
                id_type = id_type,
                by_id_type = 'uniprot',
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        if not mapped_names:

            # maybe it should be all uppercase (e.g. human gene symbols)?
            mapped_names = self._map_name(
                name = name.upper(),
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        if (
            not mapped_names and
            id_type not in {'uniprot', 'trembl', 'uniprot-sec'}
        ):

            # maybe should be capitalized (e.g. rodent gene symbols)?
            mapped_names = self._map_name(
                name = name.capitalize(),
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        if (
            not mapped_names and
            id_type not in {'uniprot', 'trembl', 'uniprot-sec'}
        ):

            # maybe it should be all lowercase?
            mapped_names = self._map_name(
                name = name.lower(),
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        if (
            not mapped_names and
            id_type.startswith('ens') and
            '.' in name
        ):

            # trying to split the part after the dot:
            mapped_names = self._map_name(
                name = name.upper().split('.')[0],
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        if (
            not mapped_names and
            ':' in name
        ):

            # trying to remove the prefix which sometimes
            # shows the ID type, e.g. CHEBI:4956 should become 4956
            mapped_names = self._map_name(
                name = common.remove_prefix(name, ':'),
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        # if a gene symbol could not be translated by the default
        # conversion table, containing only the primary gene symbols
        # in next step we try the secondary (synonym) gene symbols
        if (
            not mapped_names and
            id_type == 'genesymbol'
        ):

            mapped_names = self._map_name(
                name = name,
                id_type = 'genesymbol-syn',
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

            # for gene symbols we might try one more thing,
            # sometimes the source gene symbol missing some isoform
            # information or number because it refers to the first
            # or all isoforms or subtypes; or the opposite: the
            # original resource contains a gene symbol with a number
            # appended which is not part of the official primary
            # gene symbol
            #
            # here we try to translate by adding a number `1` or
            # by matching only the first few letters;
            # obviously we can not exclude mistranslation here
            #
            # by setting `strict = True` this step is disabled
            if not strict and not mapped_names:

                mapped_names = self._map_name(
                    name = '%s1' % name,
                    id_type = 'genesymbol',
                    target_id_type = target_id_type,
                    ncbi_tax_id = ncbi_tax_id,
                )

                if not mapped_names and target_id_type == 'uniprot':

                    mapped_names = self._map_name(
                        name = name,
                        id_type = 'genesymbol5',
                        target_id_type = target_id_type,
                        ncbi_tax_id = ncbi_tax_id,
                    )


        # for UniProt IDs we do a few more steps to
        # try to find out the primary SwissProt ID
        if uniprot_cleanup and target_id_type == 'uniprot':

            mapped_names = self.uniprot_cleanup(
                uniprots = mapped_names,
                ncbi_tax_id = ncbi_tax_id,
            )

        return mapped_names


    def uniprot_cleanup(self, uniprots, ncbi_tax_id = None):
        """
        We use this function as a standard callback when the target ID
        type is UniProt. It checks if the format of the IDs are correct,
        if they are part of the organism proteome, attempts to translate
        secondary and deleted IDs to their primary, recent counterparts.

        Args
            uniprots (str,set): One or more UniProt IDs.
            ncbi_tax_id (int): The NCBI Taxonomy identifier of the organism.

        Returns
            Set of checked and potentially translated UniProt iDs. Elements
            which do not fit the criteria will be discarded.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        uniprots = common.to_set(uniprots)

        # step 1: translate secondary IDs to primary
        uniprots = self.primary_uniprot(uniprots)

        # step 2: translate TrEMBL to SwissProt by gene symbols
        if self._trembl_swissprot_by_genesymbol:

            uniprots = self.trembl_swissprot(
                uniprots,
                ncbi_tax_id = ncbi_tax_id,
            )

        # step 3: translate deleted IDs by gene symbols
        if self._translate_deleted_uniprot:

            uniprots = self.translate_deleted_uniprots_by_genesymbol(
                uniprots
            )

        # step 4: check if the IDs exist in the proteome of the organism
        if not self._keep_invalid_uniprot:

            uniprots = self.only_valid_uniprots(
                uniprots,
                ncbi_tax_id = ncbi_tax_id,
            )

        # step 5: ensure the format validity
        uniprots = self.only_uniprot_ac(uniprots)

        return uniprots


    def map_names(
            self,
            names,
            id_type = None,
            target_id_type = None,
            ncbi_tax_id = None,
            strict = False,
            expand_complexes = True,
            uniprot_cleanup = True,
        ):
        """
        Same as ``map_name`` but translates multiple IDs at once. These two
        functions could be seamlessly implemented as one, still I created
        separate functions to always make it explicit if a set of translated
        IDs come from multiple original IDs.

        Args
            name (str): The original name to be converted.
            id_type (str): The type of the name. Available by default:
                - genesymbol (gene name)
                - entrez (Entrez Gene ID \[#\])
                - refseqp (NCBI RefSeq Protein ID \[NP\_\*|XP\_\*\])
                - ensp (Ensembl protein ID \[ENSP\*\])
                - enst (Ensembl transcript ID \[ENST\*\])
                - ensg (Ensembl genomic DNA ID \[ENSG\*\])
                - hgnc (HGNC ID \[HGNC:#\])
                - gi (GI number \[#\])
                - embl (DDBJ/EMBL/GeneBank CDS accession)
                - embl_id (DDBJ/EMBL/GeneBank accession)
                And many more, see the code of
                ``pypath.internals.input_formats``
            target_id_type (str): The name type to translate to, more or
                less the same values are available as for ``id_type``.
            ncbi_tax_id (int): NCBI Taxonomy ID of the organism.
            strict (bool): In case a Gene Symbol can not be translated,
                try to add number "1" to the end, or try to match only
                its first five characters. This option is rarely used,
                but it makes possible to translate some non-standard
                gene names typically found in old, unmaintained resources.
            expand_complexes (bool): When encountering complexes,
                translated the IDs of its components and return a set
                of IDs. The alternative behaviour is to return the
                `Complex` objects.
            uniprot_cleanup (bool): When the `target_id_type` is UniProt
                ID, call the `uniprot_cleanup` function at the end.
        """

        return set.union(
            *(
                self.map_name(
                    name = name,
                    id_type = id_type,
                    target_id_type = target_id_type,
                    ncbi_tax_id = ncbi_tax_id,
                    strict = strict,
                )
                for name in names
            )
        ) if names else set()


    def chain_map(
            self,
            name,
            id_type,
            by_id_type,
            target_id_type,
            ncbi_tax_id = None,
            **kwargs
        ):
        """
        Translate IDs which can not be directly translated in two steps:
        from `id_type` to `via_id_type` and from there to `target_id_type`.

        Args
            name (str): The original name to be converted.
            id_type (str): The type of the name.
            by_id_type (str): The intermediate name type.
            target_id_type (str): The name type to translate to, more or
                less the same values are available as for ``id_type``.
            ncbi_tax_id (int): The NCBI Taxonomy identifier of the organism.
            kwargs: Passed to `map_name`.

        Returns
            Set of IDs of type `target_id_type`.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        mapped_names = self.map_names(
            names =
                self.map_name(
                    name = name,
                    id_type = id_type,
                    target_id_type = by_id_type,
                    ncbi_tax_id = ncbi_tax_id,
                    **kwargs
                ),
            id_type = by_id_type,
            target_id_type = target_id_type,
            ncbi_tax_id = ncbi_tax_id,
            **kwargs
        )

        return mapped_names


    def _map_refseq(
            self,
            refseq,
            id_type,
            target_id_type,
            ncbi_tax_id = None,
            strict = False,
        ):
        """
        ID translation adapted to the specialities of RefSeq IDs.
        """

        mapped_names = set()
        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        # try first as it is
        mapped_names = self._map_name(
            name = refseq,
            id_type = id_type,
            target_id_type = target_id_type,
            ncbi_tax_id = ncbi_tax_id,
        )

        # then with the number at the end removed
        # this is disabled if `strict = True`
        if not mapped_names and not strict:

            mapped_names = self._map_name(
                name = refseq.split('.')[0],
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        if not mapped_names and not strict:

            rstem = refseq.split('.')[0]

            # try some other numbers
            # this risky and is disabled if `strict = True`
            for n in xrange(49):

                mapped_names.update(
                    self._map_name(
                        name = '%s.%u' % (rstem, n),
                        id_type = id_type,
                        target_id_type = target_id_type,
                        ncbi_tax_id = ncbi_tax_id,
                    )
                )

        return mapped_names


    def _map_ensp(
            self,
            ensp,
            target_id_type,
            ncbi_tax_id = None,
        ):
        """
        Special ID translation from ENSP (Ensembl peptide IDs).
        """

        mapped_names = set()
        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        # try first UniProt ID Mapping
        # then Ensembl BioMart
        for id_type in ('ensp', 'ensp_biomart'):

            if not mapped_names:

                mapped_names = self._map_name(
                    name = ensp,
                    id_type = id_type,
                    target_id_type = target_id_type,
                    ncbi_tax_id = ncbi_tax_id,
                )

        if not mapped_names:

            tax_ensp = '%u.%s' % (ncbi_tax_id, ensp)

            # this uses UniProt ID Mapping with STRING ID type
            mapped_names = self._map_name(
                name = tax_ensp,
                id_type = 'ensp_string',
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        return mapped_names


    def _map_to_ensp(
            self,
            name,
            id_type,
            ncbi_tax_id = None,
        ):
        """
        Special ID translation to ENSP (Ensembl peptide IDs).
        """

        mapped_names = set()
        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        # try first UniProt ID Mapping
        # then Ensembl BioMart
        for target_id_type in ('ensp', 'ensp_biomart'):

            if not mapped_names:

                mapped_names = self._map_name(
                    name = name,
                    id_type = id_type,
                    target_id_type = target_id_type,
                    ncbi_tax_id = ncbi_tax_id,
                )

        if not mapped_names:

            # this uses UniProt ID Mapping with STRING type
            mapped_names = self._map_name(
                name = name,
                id_type = id_type,
                target_id_type = 'ensp_string',
                ncbi_tax_id = ncbi_tax_id,
            )

            mapped_names = {n.split('.')[-1] for n in mapped_names}

        return mapped_names


    def _map_name(
            self,
            name,
            id_type,
            target_id_type,
            ncbi_tax_id = None,
        ):
        """
        Once we have defined the name type and the target name type,
        this function looks it up in the most suitable dictionary.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        tbl = self.which_table(
            id_type,
            target_id_type,
            ncbi_tax_id = ncbi_tax_id,
        )

        return tbl[name] if tbl else set()


    def translation_dict(
            self,
            id_type: str,
            target_id_type: str,
            ncbi_tax_id: int | None = None,
        ) -> MappingTable | None:
        """
        Translation table as a dict.
        """

        return self.which_table(
            id_type,
            target_id_type,
            ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id,
        )


    def translation_df(
            self,
            id_type: str,
            target_id_type: str,
            ncbi_tax_id: int | None = None,
        ) -> pd.DataFrame | None:
        """
        Translation table as a data frame.
        """

        tbl = self.translation_dict(id_type, target_id_type, ncbi_tax_id)

        if tbl:

            return pd.DataFrame(
                (
                    (source_id, target_id)
                    for source_id, target_ids in tbl.data.items() for
                    target_id in target_ids
                ),
                columns = [id_type, target_id_type],
            )


    #
    # ID specific translation methods
    #


    def label(
            self,
            name,
            entity_type = None,
            id_type = None,
            ncbi_tax_id = None,
        ):
        """
        For any kind of entity, either protein, miRNA or protein complex,
        returns the preferred human readable label. For proteins this means
        Gene Symbols, for miRNAs miRNA names, for complexes a series of
        Gene Symbols.
        """

        if isinstance(name, _const.LIST_LIKE):

            return [
                self.label(
                    _name,
                    entity_type = entity_type,
                    id_type = id_type,
                    ncbi_tax_id = ncbi_tax_id,
                )
                for _name in name
            ]

        elif hasattr(name, 'genesymbol_str'):

            return name.genesymbol_str

        elif isinstance(name, str):

            ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

            entity_type = (
                entity_type or
                (
                    'small_molecule'
                        if ncbi_tax_id == _const.NOT_ORGANISM_SPECIFIC else
                    'protein'
                )
            )

            if name.startswith('MIMAT'):

                return map_name0(
                    name,
                    id_type or 'mirbase',
                    'mir-mat-name',
                    ncbi_tax_id = ncbi_tax_id,
                ) or name

            elif name.startswith('MI'):

                return self.map_name0(
                    name,
                    id_type or 'mir-pre',
                    'mir-name',
                    ncbi_tax_id = ncbi_tax_id,
                ) or name

            elif entity_type in self.default_label_types:

                id_type = id_type or self.default_name_types[entity_type]
                target_id_type = self.default_label_types[entity_type]

                return self.map_name0(
                    name,
                    id_type = id_type,
                    target_id_type = target_id_type,
                    ncbi_tax_id = ncbi_tax_id,
                ) or name

            else:

                return self.map_name0(
                    name,
                    id_type or 'uniprot',
                    'genesymbol',
                    ncbi_tax_id = ncbi_tax_id,
                ) or name

        else:

            return str(name)


    def identifier(
            self,
            label: Union[str, Iterable[str]],
            ncbi_tax_id: Optional[int] = None,
            id_type: Optional[str] = None,
            entity_type:
                Optional[
                    Literal[
                        'drug',
                        'lncrna',
                        'mirna',
                        'protein',
                        'small_molecule',
                    ]
                ] = None,
        ) -> Union[Set[str], List[Set[str]]]:
        """
        For a label returns the corresponding primary identifier. The type
        of default identifiers is determined by the settings module. Note,
        this kind of translation is not always unambigous, one gene symbol
        might correspond to multiple UniProt IDs.
        """

        if not common.is_str(label):

            return [
                self.identifier(
                    _label,
                    entity_type = entity_type,
                    id_type = id_type,
                    ncbi_tax_id = ncbi_tax_id,
                )
                for _label in label
            ]

        elif hasattr(label, 'components'):

            return label.__str__()

        elif common.is_str(label):

            ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

            entity_type = (
                entity_type or
                (
                    'small_molecule'
                        if ncbi_tax_id == _const.NOT_ORGANISM_SPECIFIC else
                    'protein'
                )
            )

            id_type = (
                id_type or
                settings.get('default_label_types')[entity_type]
            )

            target_id_type = settings.get('default_name_types')[entity_type]

            return self.map_name(
                label,
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        else:

            return str(name)


    def identifier0(
            self,
            label: Union[str, Iterable[str]],
            ncbi_tax_id: Optional[int] = None,
            id_type: Optional[str] = None,
            entity_type:
                Optional[
                    Literal[
                        'drug',
                        'lncrna',
                        'mirna',
                        'protein',
                        'small_molecule',
                    ]
                ] = None,
        ) -> Union[str, List[str]]:

        args = locals()
        _ = args.pop('self')
        ids = self.identifier(**args)

        return (
            common.first(ids)
                if isinstance(label, str) else
            list(map(common.first, ids))
        )


    def guess_type(self, name, entity_type = None):
        """
        From a string, tries to guess the ID type and optionally the entity
        type. Returns a tuple of strings: ID type and entity type.
        """

        if (
            (
                not entity_type or
                entity_type == 'protein'
            ) and
            self.reuniprot.match(name)
        ):

            return 'uniprot', 'protein'

        if not entity_type or entity_type == 'mirna':

            if self.remipreac.match(name):

                return 'mir-pre', 'mirna'

            if self.remimatac.match(name):

                return 'mirbase', 'mirna'

            if self.remimatid.match(name):

                return 'mir-mat-name', 'mirna'

            if self.remipreid.match(name):

                return 'mir-name', 'mirna'

        return None, entity_type


    def id_from_label(
            self,
            label,
            label_id_type = 'genesymbol',
            ncbi_tax_id = None,
        ):

        if label_id_type in self.label_type_to_id_type:

            ids = self.map_name(
                label,
                label_id_type,
                self.label_type_to_id_type[label_id_type],
                ncbi_tax_id = ncbi_tax_id,
            )

        return ids or {label}

    def id_from_label0(
            self,
            label,
            label_id_type = 'genesymbol',
            ncbi_tax_id = None,
        ):

        return next(
            self.id_from_label(
                label = label,
                label_id_type = label_id_type,
                ncbi_tax_id = ncbi_tax_id
            ).__iter__()
        )


    def primary_uniprot(self, uniprots, ncbi_tax_id = None):
        """
        For an iterable of UniProt IDs returns a set with the secondary IDs
        changed to the corresponding primary IDs. Anything what is not a
        secondary UniProt ID left intact.
        """

        primaries = set()
        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        for uniprot in uniprots:

            primary = self.map_name(
                name = uniprot,
                id_type = 'uniprot-sec',
                target_id_type = 'uniprot-pri',
                ncbi_tax_id = ncbi_tax_id,
            )

            if primary:

                primaries.update(primary)

            else:

                # most probably this UniProt is already primary
                primaries.add(uniprot)

        return primaries


    def trembl_swissprot(self, uniprots, ncbi_tax_id = None):
        """
        For an iterable of TrEMBL and SwissProt IDs, returns a set with
        only SwissProt, mapping from TrEMBL to gene symbols, and
        then back to SwissProt. If this kind of translation is not successful
        for any of the IDs it will be kept in the result, no matter if it's
        not a SwissProt ID. If the
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id
        swissprots = set()

        for uniprot in uniprots:

            swissprot = None
            genesymbols = self.map_name(
                name = uniprot,
                id_type = 'trembl',
                target_id_type = 'genesymbol',
                ncbi_tax_id = ncbi_tax_id,
            )

            this_swissprots = self.map_names(
                names = genesymbols,
                id_type = 'genesymbol',
                target_id_type = 'swissprot',
                ncbi_tax_id = ncbi_tax_id,
            )

            if not this_swissprots:

                swissprots.add(uniprot)

            else:

                swissprots.update(this_swissprots)

        return swissprots


    def translate_deleted_uniprots_by_genesymbol(
            self,
            uniprots,
            ncbi_tax_id = None,
        ):

        if isinstance(uniprots, str):

            return self.translate_deleted_uniprot_by_genesymbol(
                uniprots,
                ncbi_tax_id = ncbi_tax_id,
            )

        else:

            ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

            return set.union(*(
                self.translate_deleted_uniprot_by_genesymbol(
                    uniprot,
                    ncbi_tax_id = ncbi_tax_id,
                )
                for uniprot in uniprots
            )) if uniprots else set()


    def translate_deleted_uniprot_by_genesymbol(
            self,
            uniprot,
            ncbi_tax_id = None,
        ):
        """
        Due to potentially ambiguous translation always returns set.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        if uniprot_db.is_uniprot(uniprot, organism = ncbi_tax_id):

            return {uniprot}

        elif self.other_organism_uniprot(uniprot, ncbi_tax_id = ncbi_tax_id):

            return set()

        else:

            genesymbol, taxid = self.deleted_uniprot_genesymbol(uniprot)

            if genesymbol and taxid == ncbi_tax_id:

                return self.map_name(
                    genesymbol,
                    'genesymbol',
                    'uniprot',
                    ncbi_tax_id = ncbi_tax_id,
                    uniprot_cleanup = False,
                )

        return {uniprot}


    def other_organism_uniprot(self, uniprot, ncbi_tax_id = None):
        """
        Tells if ``uniprot`` is an UniProt ID from some other organism than
        ``ncbi_tax_id``.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        uniprot_taxid = taxonomy.uniprot_taxid(uniprot)

        return uniprot_taxid and uniprot_taxid != ncbi_tax_id


    def deleted_uniprot_genesymbol(self, uniprot):

        return uniprot_input.deleted_uniprot_genesymbol(uniprot)


    def only_valid_uniprots(self, uniprots, ncbi_tax_id = None):

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        if isinstance(uniprots, str):

            return self.valid_uniprot(uniprots, ncbi_tax_id = ncbi_tax_id)

        else:

            return {
                uniprot
                for uniprot in uniprots
                if uniprot_db.is_uniprot(uniprot, organism = ncbi_tax_id)
            }


    def valid_uniprot(self, uniprot, ncbi_tax_id = None):
        """
        If the UniProt ID ``uniprot`` exist in the proteome of the organism
        ``ncbi_tax_id`` returns the ID, otherwise returns None.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        if uniprot_db.is_uniprot(uniprot, organism = ncbi_tax_id):

            return uniprot


    def only_uniprot_ac(self, uniprots):
        """
        For one or more strings returns only those which match the format
        of UniProt accession numbers.
        The format is defined here:
        https://www.uniprot.org/help/accession_numbers

        If string provided, returns string or None.
        If iterable provided, returns set (potentially empty if none of the
        strings are valid).
        """

        if isinstance(uniprots, str):

            return self._only_uniprot_ac(uniprots)

        else:

            return {
                validated
                for validated in
                (
                    self._only_uniprot_ac(uniprot)
                    for uniprot in uniprots
                )
                if validated
            }


    def _only_uniprot_ac(self, uniprot):

        return uniprot if uniprot_input.valid_uniprot(uniprot) else None

    #
    # Mapping table management methods
    #

    @staticmethod
    def mapping_tables():
        """
        List of mapping tables available to load.

        Returns
            (list): A list of tuples, each representing an ID translation
                table, with the ID types, the data source and the loader
                class.
        """

        MappingTableDefinition = collections.namedtuple(
            'MappingTableDefinition',
            (
                'id_type_a',
                'id_type_b',
                'resource',
                'input_class',
                'resource_id_type_a',
                'resource_id_type_b',
            ),
        )

        MappingTableDefinition.__new__.__defaults__ = (None, None)

        result = []

        for resource_attr in RESOURCES_EXPLICIT:

            resources = getattr(maps, resource_attr)

            for (id_type_a, id_type_b), inputdef in iteritems(resources):

                result.append(
                    MappingTableDefinition(
                        id_type_a = id_type_a,
                        id_type_b = id_type_b,
                        resource = resource_attr,
                        input_class = inputdef.__class__.__name__,
                        resource_id_type_a = inputdef._resource_id_type_a,
                        resource_id_type_b = inputdef._resource_id_type_b,
                    )
                )

        for service_ids, service_id_type, input_cls in RESOURCES_IMPLICIT:

            service_ids = (
                iteritems(service_ids)
                    if isinstance(service_ids, dict) else
                zip(*(service_ids,) * 2)
            )

            for id_type, resource_id_type in service_ids:

                id_type_b = 'pro' if service_id_type == 'pro' else None

                result.append(
                    MappingTableDefinition(
                        id_type_a = id_type,
                        id_type_b = id_type_b,
                        resource = service_id_type,
                        input_class = input_cls.__name__,
                        resource_id_type_a = resource_id_type,
                        resource_id_type_b = None,
                    )
                )

        return result


    @classmethod
    def id_types(cls):
        """
        A list of all identifier types that can be handled by any of the
        resources.

        Returns
            (list): A list of tuples with the identifier type labels used
                in pypath and in the original resource. If the latter is
                None, typically the ID type has no name in the original
                resource.
        """

        IdType = collections.namedtuple(
            'IdType',
            (
                'pypath',
                'original',
            ),
        )

        return {
            IdType(
                pypath = getattr(mapdef, 'id_type_%s' % side),
                original = getattr(mapdef, 'resource_id_type_%s' % side),
            )
            for mapdef in cls.mapping_tables()
            for side in ('a', 'b')
            if getattr(mapdef, 'id_type_%s' % side)
        }


    def has_mapping_table(
            self,
            id_type,
            target_id_type,
            ncbi_tax_id = None,
        ):
        """
        Tells if a mapping table is loaded. If it's loaded, it resets the
        expiry timer so the table remains loaded.

        Returns
            (bool): True if the mapping table is loaded.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        key = self.get_table_key(
            id_type = id_type,
            target_id_type = target_id_type,
            ncbi_tax_id = ncbi_tax_id,
        )

        if key in self.tables:

            self.tables[key]._used()

        return key in self.tables


    def load_mapping(
            self,
            resource,
            **kwargs
        ):
        """
        Loads a single mapping table based on input definition in
        ``resource``. ``**kwargs`` passed to ``MapReader``.
        """

        if (
            resource.type in {'file', 'pickle'} and
            not (
                os.path.exists(resource.input) or
                inputs.get_method(resource.input)
            )
        ):

            self._log(
                'Could not load mapping: no such '
                'file or function: `%s`.' % resource.input
            )
            return

        ncbi_tax_id = kwargs.get('ncbi_tax_id', resource.ncbi_tax_id)

        self._log(
            'Loading mapping table for organism `%s` '
            'with identifiers `%s` and `%s`, '
            'input type `%s`' % (
                ncbi_tax_id,
                resource.id_type_a,
                resource.id_type_b,
                resource.type,
            )
        )

        reader = MapReader(param = resource, **kwargs)

        a_to_b = reader.mapping_table_a_to_b
        b_to_a = reader.mapping_table_b_to_a

        for sides in (('a', 'b'), ('b', 'a')):

            table = locals()['%s_to_%s' % sides]

            if table:
                self._log(
                    'Sucessfully loaded mapping table for organism `%s` '
                    'with identifiers `%s` to `%s`.' % (
                        str(ncbi_tax_id),
                        getattr(resource, f'id_type_{sides[0]}'),
                        getattr(resource, f'id_type_{sides[1]}'),
                    )
                )
                self.tables[table.get_key()] = table


    def swissprots(self, uniprots, ncbi_tax_id = None):
        """
        Creates a dict translating a set of potentially secondary and
        non-reviewed UniProt IDs to primary SwissProt IDs (whenever
        is possible).
        """

        swissprots = {}

        for uniprot in uniprots:

            swissprots[uniprot] = self.map_name(
                name = uniprot,
                id_type = 'uniprot',
                target_id_type = 'uniprot',
                ncbi_tax_id = ncbi_tax_id,
            )

        return swissprots


    def load_genesymbol5(self, ncbi_tax_id = None):
        """
        Creates a Gene Symbol to UniProt mapping table with the first
        5 characters of each Gene Symbol.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        genesymbol_table = self.which_table(
            id_type = 'genesymbol',
            target_id_type = 'uniprot',
            ncbi_tax_id = ncbi_tax_id,
        )
        genesymbol_syn_table = self.which_table(
            id_type = 'genesymbol-syn',
            target_id_type = 'uniprot',
            ncbi_tax_id = ncbi_tax_id,
        )

        genesymbol5_data = collections.defaultdict(set)

        for table in (genesymbol_table, genesymbol_syn_table):

            for genesymbol, uniprots in iteritems(table.data):

                if len(genesymbol) >= 5:

                    genesymbol5 = genesymbol[:5]

                    genesymbol5_data[genesymbol5].update(uniprots)

        mapping_table = MappingTable(
            data = genesymbol5_data,
            id_type = 'genesymbol5',
            target_id_type = 'uniprot',
            ncbi_tax_id = ncbi_tax_id,
        )

        self.tables[mapping_table.get_key()] = mapping_table


    def load_uniprot_static(
            self,
            keys,
            ncbi_tax_id = None,
        ):
        """
        Loads mapping tables from the huge static mapping file from UniProt.
        Takes long to download and process, also requires more memory. This
        is the last thing we try if everything else failed.
        """

        cachedir = cache_mod.get_cachedir()
        data = dict((key, collections.defaultdict(set)) for key in keys)
        cache_files = {}
        to_load = set()
        id_type_b = 'uniprot'

        # attempting to load them from Pickle
        for key in keys:

            mapping_id = common.md5(
                json.dumps(
                    (
                        key,
                        'uniprot_static',
                    )
                )
            )

            cachefile = os.path.join(cachedir, mapping_id)
            cache_files[key] = cachefile

            if os.path.exists(cachefile):

                with open(cachefile, 'rb') as fp:

                    data[key] = pickle.load(fp)

            else:

                to_load.add(key)

        # loading the remaining from the big UniProt mapping file:
        if to_load:

            url = urls.urls['uniprot_idmap_ftp']['url']
            c = curl.Curl(url, silent = False, large = True)

            prg = progress.Progress(
                c.size,
                'Processing ID conversion list',
                99,
            )

            for line in c.result:

                prg.step(len(line))

                line = common.decode(line, 'ascii').strip().split('\t')

                if len(line) > 2 and line[1] in self.names_uniprot_static:

                    id_type_a = self.names_uniprot_static[line[1]]

                    key_a_to_b = MappingTableKey(
                        id_type = id_type_a,
                        target_id_type = id_type_b,
                        ncbi_tax_id = ncbi_tax_id,
                    )
                    key_b_to_a = MappingTableKey(
                        id_type = id_type_b,
                        target_id_type = id_type_a,
                        ncbi_tax_id = ncbi_tax_id,
                    )

                    this_uniprot = line[0].split('-')[0]

                    if key_a_to_b in to_load:

                        data[key_a_to_b][line[2]].add(this_uniprot)

                    if key_b_to_a in to_load:

                        data[key_b_to_a][this_uniprot].add(line[2])

            prg.terminate()

            for key, this_data in iteritems(data):

                pickle.dump(this_data, open(cache_files[key], 'wb'))

        for key, this_data in iteritems(data):

            table = MappingTable(
                data = this_data,
                id_type = key,
                target_id_type = id_type_b,
                ncbi_tax_id = ncbi_tax_id,
                lifetime = 600,
            )

            self.tables[key] = table


    def remove_table(self, id_type, target_id_type, ncbi_tax_id):
        """
        Removes the table defined by the ID types and organism.
        """

        key = MappingTableKey(
            id_type = id_type,
            target_id_type = target_id_type,
            ncbi_tax_id = ncbi_tax_id,
        )

        self.remove_key(key)


    def remove_key(self, key):
        """
        Removes the table with key ``key`` if exists.
        """

        if key in self.tables:

            if key and len(key) == 3:

                self._log(
                    'Removing mapping table `%s` '
                    'to `%s` for organism `%u`.' % key
                )

            del self.tables[key]


    def remove_expired(self):
        """
        Removes tables last used a longer time ago than their lifetime.
        """

        to_remove = set()

        for key, table in iteritems(self.tables):

            if not table or table._expired():

                to_remove.add(key)

        for key in to_remove:

            self.remove_key(key)
