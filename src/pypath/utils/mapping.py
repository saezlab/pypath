#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
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

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import os
import sys
import math
import re
import importlib as imp
import copy
import itertools
import collections
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

import timeloop
# we use this for simple little tasks only
# and don't want engage another logger
timeloop.app.logging.disable(level = 9999)

# from pypath:
import pypath.share.progress as progress
import pypath.share.common as common
import pypath.share.cache as cache_mod
import pypath.internals.maps as maps
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs as inputs
import pypath.inputs.uniprot as uniprot_input
import pypath.inputs.pro as pro_input
import pypath.inputs.biomart as biomart_input
import pypath.internals.input_formats as input_formats
import pypath.utils.reflists as reflists
import pypath.utils.taxonomy as taxonomy
import pypath.share.settings as settings
import pypath.share.session as session_mod
_logger = session_mod.get_log()


__all__ = ['MapReader', 'MappingTable', 'Mapper']

"""
Classes for reading and use serving ID mapping data
from UniProt, file, mysql or pickle.
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


class MapReader(session_mod.Logger):
    """
    Reads ID translation data and creates ``MappingTable`` instances.
    When initializing ID conversion tables for the first time
    data is downloaded from UniProt and read into dictionaries.
    It takes a couple of seconds. Data is saved to pickle
    dumps, this way later the tables load much faster.

    :arg source_type str:
        Type of the resource, either `file`, `uniprot` or `unprotlist`.
    """

    def __init__(
            self,
            param,
            ncbi_tax_id = None,
            entity_type = None,
            load_a_to_b = True,
            load_b_to_a = False,
            uniprots = None,
            lifetime = 300,
        ):
        """
        entity_type : str
            An optional, custom string showing the type of the entities,
            e.g. `protein`. This is not mandatory for the identification
            of mapping tables, hence the same name types can't be used
            for different entities. E.g. if both proteins and miRNAs have
            Entrez gene IDs then these should be different ID types (e.g.
            `entrez_protein` and `entrez_mirna`) or both protein and miRNA
            IDs can be loaded into one mapping table and simply called
            `entrez`.
        uniprots : set
            UniProt IDs to query in case the source of the mapping table
            is the UniProt web service.
        lifetime : int
            If this table has not been used for longer than this preiod it is
            to be removed at next cleanup. Time in seconds. Passed to
            ``MappingTable``.
        """

        session_mod.Logger.__init__(self, name = 'mapping')

        self.ncbi_tax_id = (
            ncbi_tax_id or
            param.ncbi_tax_id or
            settings.get('default_organism')
        )

        self._log(
            'Reader created for ID translation table, parameters: '
            '`ncbi_tax_id=%u, id_a=%s, id_b=%s, '
            'load_a_to_b=%u, load_b_to_a=%u, '
            'input_type=%s (%s)`.' % (
                self.ncbi_tax_id,
                param.id_type_a,
                param.id_type_b,
                load_a_to_b,
                load_b_to_a,
                param.type,
                param.__class__.__name__,
            )
        )

        self.cachedir = cache_mod.get_cachedir()

        self.id_type_a = param.id_type_a
        self.id_type_b = param.id_type_b
        self.load_a_to_b = load_a_to_b
        self.load_b_to_a = load_b_to_a
        self.entity_type = entity_type
        self.source_type = param.type
        self.param = param
        self.lifetime = lifetime
        self.a_to_b = None
        self.b_to_a = None
        self.uniprots = uniprots

        self.load()


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def load(self):
        """
        The complete process of loading mapping tables. First sets up the
        paths of the cache files, then loads the tables from the cache files
        or the original sources if necessary. Upon successful loading from an
        original source writes the results to cache files.
        """

        self.use_cache = settings.get('mapping_use_cache')
        self.setup_cache()

        if self.use_cache:

            self.read_cache()

        if not self.tables_loaded():

            # read from the original source
            self.read()

            if self.tables_loaded():

                # write cache only at successful loading
                self.write_cache()


    @property
    def mapping_table_a_to_b(self):
        """
        Returns a ``MappingTable`` instance created from the already
        loaded data.
        """

        return self._get_mapping_table('a', 'b')


    @property
    def mapping_table_b_to_a(self):
        """
        Returns a ``MappingTable`` instance created from the already
        loaded data.
        """

        return self._get_mapping_table('b', 'a')


    def id_type_side(self, id_type):

        return (
            'a'
                if id_type == self.id_type_a else
            'b'
                if id_type == self.id_type_b else
            None
        )


    def _get_mapping_table(self, *args):

        data = getattr(self, '%s_to_%s' % args)
        id_type = getattr(self, 'id_type_%s' % args[0])
        target_id_type = getattr(self, 'id_type_%s' % args[1])

        if isinstance(data, dict):

            return MappingTable(
                data = data,
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = self.ncbi_tax_id,
                lifetime = self.lifetime,
            )


    def tables_loaded(self):
        """
        Tells if the requested tables have been created.
        """

        return (
            (bool(self.a_to_b) or not self.load_a_to_b) and
            (bool(self.b_to_a) or not self.load_b_to_a)
        )


    def write_cache(self):
        """
        Exports the ID translation data into pickle files.
        """

        self._write_cache('a', 'b')
        self._write_cache('b', 'a')


    def _write_cache(self, *args):

        data = getattr(self, '%s_to_%s' % args)

        if self._to_be_loaded(*args) and data:

            cachefile = self._attr('cachefile', *args)

            self._remove_cache_file(*args)

            pickle.dump(data, open(cachefile, 'wb'))


    def read_cache(self):
        """
        Reads the ID translation data from a previously saved pickle file.
        """

        self._read_cache('a', 'b')
        self._read_cache('b', 'a')


    def _read_cache(self, *args):

        if self._to_be_loaded(*args):

            cachefile = self._attr('cachefile', *args)

            if os.path.exists(cachefile):

                setattr(
                    self,
                    '%s_to_%s' % args,
                    pickle.load(open(cachefile, 'rb')),
                )
                self._log(
                    'Loading `%s` to `%s` mapping table '
                    'from pickle file `%s`.' % (
                        self.param.id_type_a,
                        self.param.id_type_b,
                        cachefile,
                    )
                )


    def _to_be_loaded(self, *args):

        return self._attr('load', *args)


    def _attr(self, attr, *args):

        return getattr(self, self._attr_name(attr, *args))


    @staticmethod
    def _attr_name(attr, *args):

        return '%s_%s_to_%s' % ((attr,) + args)


    def read(self):
        """
        Reads the ID translation data from the original source.
        """

        method = 'read_mapping_%s' % self.source_type

        if hasattr(self, method):

            getattr(self, method)()


    def setup_cache(self):
        """
        Constructs the cache file path as md5 hash of the parameters.
        """

        self._setup_cache('a', 'b')
        self._setup_cache('b', 'a')


    def _setup_cache(self, *args):

        mapping_id_attr = self._attr_name('mapping_id', *args)
        cachefile_attr  = self._attr_name('cachefile', *args)

        setattr(
            self,
            mapping_id_attr,
            self._get_mapping_id(*args),
        )

        setattr(
            self,
            cachefile_attr,
            os.path.join(self.cachedir, getattr(self, mapping_id_attr)),
        )


    def _get_mapping_id(self, *args):
        """
        Returns an md5 checksum unambigously identifying the mapping table
        by the identifiers, the direction of translation, the organism
        and other parameters like, for example, the source URL.
        """

        return common.md5(
            json.dumps(
                (
                    getattr(self, 'id_type_%s' % args[0]),
                    getattr(self, 'id_type_%s' % args[1]),
                    self.ncbi_tax_id,
                    sorted(self.param.__dict__.items())
                )
            )
        )


    def _cache_files_exist(self):
        """
        Checks if both cache files are either not necessary or exist.
        """

        return (
            self.cache_file_exists('a', 'b') and
            self.cache_file_exists('b', 'a')
        )


    def _cache_file_exists(self, *args):
        """
        Checks if a cache file is either not necessary or exists.
        """

        return (
            not self._attr('load', *args) or
            os.path.isfile(self._attr('cachefile', *args))
        )


    def _remove_cache_file(self, *args):

        cachefile = self._attr('cachefile', *args)

        if os.path.exists(cachefile):

            self._log('Removing mapping table cache file `%s`.' % cachefile)
            os.remove(cachefile)


    def read_mapping_file(self):

        if not os.path.exists(self.param.input):

            method = inputs.get_method(self.param.input)

            if not method:

                return {}

            else:

                input_args = (
                    self.param.input_args
                        if hasattr(self.param, 'input_args') else
                    {}
                )
                infile = method(**input_args)

        else:

            infile = open(self.param.input, encoding = 'utf-8', mode = 'r')
            total = os.path.getsize(self.param.input)

        a_to_b = collections.defaultdict(set)
        b_to_a = collections.defaultdict(set)

        for i, line in enumerate(infile):

            if self.param.header and i < self.param.header:

                continue

            if hasattr(line, 'decode'):

                line = line.decode('utf-8')

            if hasattr(line, 'rstrip'):

                line = line.rstrip().split(self.param.separator)

            if len(line) < max(self.param.col_a, self.param.col_b):

                continue

            id_a = line[self.param.col_a]
            id_b = line[self.param.col_b]

            if self.load_a_to_b:

                a_to_b[id_a].add(id_b)

            if self.load_b_to_a:

                b_to_a[id_b].add(id_a)

        if hasattr(infile, 'close'):

            infile.close()

        self.a_to_b = a_to_b if self.load_a_to_b else None
        self.b_to_a = b_to_a if self.load_b_to_a else None


    def read_mapping_uniprot_list(self):
        """
        Builds a mapping table by downloading data from UniProt's
        upload lists service.
        """

        a_to_b = collections.defaultdict(set)
        b_to_a = collections.defaultdict(set)

        if not self.uniprots:

            self.set_uniprot_space(swissprot = False)

        # We need a list to query this service, and we have method only for
        # getting a proteome wide list of UniProt IDs. If the translated
        # ID type is not UniProt, then first we need to translate the
        # proteome wide reference list from UniProt to the target ID type.
        if self.param.id_type_a != 'uniprot':

            u_target = self._read_mapping_uniprot_list(
                uniprot_id_type_a = 'ACC',
                uniprot_id_type_b = self.param.uniprot_id_type_a,
            )

            upload_ac_list = [l.split('\t')[1].strip() for l in u_target]

        else:

            upload_ac_list = self.uniprots

        uniprot_data = self._read_mapping_uniprot_list(
            upload_ac_list = upload_ac_list,
        )

        for l in uniprot_data:

            if not l:

                continue

            l = l.strip().split('\t')

            if self.load_a_to_b:

                a_to_b[l[0]].add(l[1])

            if self.load_b_to_a:

                b_to_a[l[1]].add(l[0])

        self.a_to_b = a_to_b if self.load_a_to_b else None
        self.b_to_a = b_to_a if self.load_b_to_a else None


    def set_uniprot_space(self, swissprot = None):
        """
        Sets up a search space of UniProt IDs.
        """

        swissprot = self.param.swissprot if swissprot is None else swissprot

        self.uniprots = uniprot_input.all_uniprots(
            self.ncbi_tax_id,
            swissprot = swissprot,
        )


    def _read_mapping_uniprot_list(
            self,
            uniprot_id_type_a = None,
            uniprot_id_type_b = None,
            upload_ac_list = None,
            chunk_size = None,
        ):
        """
        Reads a mapping table from UniProt "upload lists" service.
        """

        chunk_size = (
            chunk_size or
            settings.get('uniprot_uploadlists_chunk_size')
        )
        uniprot_id_type_a = uniprot_id_type_a or self.param.uniprot_id_type_a
        uniprot_id_type_b = uniprot_id_type_b or self.param.uniprot_id_type_b

        upload_ac_list = upload_ac_list or self.uniprots
        upload_ac_list = sorted(upload_ac_list)

        self._log(
            'Querying the UniProt uploadlists service for ID translation '
            'data. Querying a list of %u IDs.' % len(upload_ac_list)
        )

        url = urls.urls['uniprot_basic']['lists']

        result = []

        # loading data in chunks of 10,000 by default
        for i in range(math.ceil(len(upload_ac_list) / chunk_size)):

            this_chunk = upload_ac_list[i * chunk_size:(i + 1) * chunk_size]

            self._log(
                'Request to UniProt uploadlists, chunk #%u with %u IDs.' % (
                    i,
                    len(this_chunk),
                )
            )

            post = {
                'from': uniprot_id_type_a,
                'format': 'tab',
                'to': uniprot_id_type_b,
                'uploadQuery': ' '.join(sorted(this_chunk)),
            }

            c = curl.Curl(url, post = post, large = True, silent = False)

            # 3 extra attempts
            if c.result is None:

                for i in xrange(3):

                    c = curl.Curl(
                        url,
                        post = post,
                        large = True,
                        silent = False,
                        cache = False,
                    )

                    if c.result is not None:

                        break

            if c.result is None or c.fileobj.read(5) == '<!DOC':

                _logger.console(
                    'Error at downloading ID mapping data from UniProt.', -9
                )

                c.result = ''

            c.fileobj.seek(0)

            # removing the header row
            _ = next(c.result)

            result.extend(list(c.fileobj)[1:])

        return result


    def read_mapping_uniprot(self):
        """
        Downloads ID mappings directly from UniProt.
        See the names of possible identifiers here:
        http://www.uniprot.org/help/programmatic_access
        """

        resep = re.compile(r'[\s;]')
        recolend = re.compile(r'$;')

        a_to_b = collections.defaultdict(set)
        b_to_a = collections.defaultdict(set)


        rev = (
            ''
                if not self.param.swissprot else
            ' AND reviewed:%s' % self.param.swissprot
        )
        query = 'organism:%u%s' % (int(self.ncbi_tax_id), rev)

        url = urls.urls['uniprot_basic']['url']
        post = {
            'query': query,
            'format': 'tab',
            'columns': 'id,%s%s' % (
                self.param.field,
                ''
                    if self.param.subfield is None else
                '(%s)' % self.param.subfield
            ),
        }

        url = '%s?%s' % (url, urllib.urlencode(post))
        c = curl.Curl(url, silent = False)
        # fallback to empty string: Curl returns None in case of
        # empty file but in case of UniProt, especially for under-researched
        # taxons it can happen there is no result for certain queries
        data = c.result or ''

        data = [
            [
                [xx]
                    if self.param.field == 'protein names' else
                [
                    xxx for xxx in resep.split(recolend.sub('', xx.strip()))
                    if len(xxx) > 0
                ]
                for xx in x.split('\t') if len(xx.strip()) > 0
            ]
            for x in data.split('\n') if len(x.strip()) > 0
        ]

        if data:

            del data[0]

            for l in data:

                if len(l) >= 2:

                    l[1] = (
                        self._process_protein_name(l[1][0])
                            if self.param.field == 'protein names' else
                        l[1]
                    )

                    for other in l[1]:

                        if self.load_a_to_b:

                            a_to_b[other].add(l[0][0])

                        if self.load_b_to_a:

                            b_to_a[l[0][0]].add(other)

        self.a_to_b = a_to_b if self.load_a_to_b else None
        self.b_to_a = b_to_a if self.load_b_to_a else None


    def read_mapping_pro(self):

        pro_data = pro_input.pro_mapping(target_id_type = self.param.id_type)

        pro_to_other = collections.defaultdict(set)

        for pro, other in pro_data:

            pro_to_other[pro].add(other)

        self.a_to_b = (
            None
                if not self.load_a_to_b else
            common.swap_dict(pro_to_other)
                if self.param.to_pro else
            dict(pro_to_other)
        )
        self.b_to_a = (
            None
                if not self.load_b_to_a else
            dict(pro_to_other)
                if self.param.to_pro else
            common.swap_dict(pro_to_other)
        )


    def read_mapping_biomart(self):

        biomart_data = biomart_input.biomart_query(attrs = self.param.attrs)

        a_to_b = collections.defaultdict(set)
        b_to_a = collections.defaultdict(set)

        for rec in biomart_data:

            id_a = getattr(rec, self.param.biomart_id_type_a)
            id_b = getattr(rec, self.param.biomart_id_type_b)

            if id_a and id_b:

                if self.load_a_to_b:

                    a_to_b[id_a].add(id_b)

                if self.load_b_to_a:

                    b_to_a[id_b].add(id_a)



        self.a_to_b = dict(a_to_b) if self.load_a_to_b else None
        self.b_to_a = dict(b_to_a) if self.load_b_to_a else None


    @staticmethod
    def _process_protein_name(name):

        rebr = re.compile(r'\(([^\)]{3,})\)')
        resq = re.compile(r'\[([^\]]{3,})\]')

        names = [name.split('(')[0]]
        names += rebr.findall(name)
        others = common.flat_list([x.split(';') for x in resq.findall(name)])
        others = [x.split(':')[1] if ':' in x else x for x in others]
        others = [x.split('(')[1] if '(' in x else x for x in others]
        names += others

        return {x.strip() for x in names}


class MappingTable(session_mod.Logger):
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

        session_mod.Logger.__init__(self, name = 'mapping')

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


    def _used(self):

        self._last_used = time.time()


    def _expired(self):

        return time.time() - self._last_used > self.lifetime


    def get_key(self):

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

        return 'MappingTable from=`%s`, to=`%s`, taxon=`%u`' % self.key


class Mapper(session_mod.Logger):

    default_name_types = settings.get('default_name_types')

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
            cleanup_period = 10,
            lifetime = 300,
            translate_deleted_uniprot = None,
            keep_invalid_uniprot = None,
            trembl_swissprot_by_genesymbol = None,
        ):
        """
        cleanup_period : int
            Periodically check and remove unused mapping data.
            Time in seconds. If `None` tables kept forever.
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

        cleanup_period = settings.get(
            'mapper_cleanup_interval',
            cleanup_period
        )
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

        self._mapper_cleanup_timeloop = timeloop.Timeloop()

        for job in self._mapper_cleanup_timeloop.jobs:

            if job.is_alive():

                job.stop()
                job.stopped.set()

        self._mapper_cleanup_timeloop.jobs = []


        @self._mapper_cleanup_timeloop.job(
            interval = datetime.timedelta(
                seconds = cleanup_period
            )
        )
        def _cleanup():

            self.remove_expired()


        self._mapper_cleanup_timeloop.start(block = False)

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


    def reload(self):

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
        tbl_key_noorganism = tbl_key[:-1] + (-1,)

        tbl_key_rev = self.get_table_key(
            target_id_type = target_id_type,
            id_type = id_type,
            ncbi_tax_id = ncbi_tax_id,
        )
        tbl_key_rev_noorganism = tbl_key_rev[:-1] + (-1,)

        if tbl_key in self.tables:

            tbl = self.tables[tbl_key]

        elif tbl_key_noorganism in self.tables:

            tbl = self.tables[tbl_key_noorganism]

        elif tbl_key_rev in self.tables:

            self.create_reverse(tbl_key_rev)
            tbl = self.tables[tbl_key]

        elif tbl_key_rev_noorganism in self.tables:

            self.create_reverse(tbl_key_rev_noorganism)
            tbl = self.tables[tbl_key_rev_noorganism]

        elif load:

            id_types = (id_type, target_id_type)
            id_types_rev = tuple(reversed(id_types))
            resource = None

            for resource_attr in ['uniprot', 'basic', 'mirbase', 'ipi']:

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

                for service_ids, service_id_type, input_cls in (
                    (
                        input_formats.ac_mapping,
                        'uniprot',
                        input_formats.UniprotListMapping,
                    ),
                    (
                        input_formats.pro_mapping,
                        'pro',
                        input_formats.ProMapping,
                    ),
                    (
                        input_formats.biomart_mapping,
                        'biomart',
                        input_formats.BiomartMapping,
                    )
                ):

                    if (
                        (
                            service_id_type == 'uniprot' and (
                                id_type in service_ids and
                                target_id_type in service_ids and
                                id_type != target_id_type
                            )
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

                        # for uniprot/uploadlists or PRO
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

                if id_type in self.uniprot_static_names:

                    self.load_uniprot_static([id_type])

                    tbl = check_loaded()

        if hasattr(tbl, '_used'):

            tbl._used()

        return tbl


    @staticmethod
    def reverse_mapping(mapping_table):

        rev_data = common.swap_dict(mapping_table.data)

        return MappingTable(
            data = rev_data,
            lifetime = mapping_table.lifetime,
        )


    def reverse_key(self, key):

        self.get_table_key(
            id_type = key.target_id_type,
            target_id_type = key.id_type,
            ncbi_tax_id = key.ncbi_tax_id,
        )


    def create_reverse(self, key):
        """
        Creates a mapping table with ``id_type`` and ``target_id_type``
        (i.e. direction of the ID translation) swapped.
        """

        table = self.mappings[key]
        rev_key = self.reverse_key(key)

        self.tables[rev_key] = self.reverse_mapping(table)


    def map_name0(
            self,
            name,
            id_type = None,
            target_id_type = None,
            ncbi_tax_id = None,
            strict = False,
            silent = True,
            nameType = None,
            targetNameType = None,
        ):

        names = self.map_name(
            name = name,
            id_type = id_type,
            target_id_type = target_id_type,
            ncbi_tax_id = ncbi_tax_id,
            strict = strict,
            silent = silent,
            nameType = nameType,
            targetNameType = targetNameType,
        )

        return list(names)[0] if names else None

    def map_name(
            self,
            name,
            id_type = None,
            target_id_type = None,
            ncbi_tax_id = None,
            strict = False,
            silent = True,
            expand_complexes = True,
            uniprot_cleanup = True,
            nameType = None,
            targetNameType = None,
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

        name : str
            The original name to be converted.
        id_type : str
            The type of the name.
            Available by default:
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
            To use other IDs, you need to define the input method
            and load the table before calling
            :py:func:Mapper.map_name().
        target_id_type : str
            The name type to translate to, more or less the same values
            are available as for ``id_type``.
        nameType : str
            Deprecated. Synonym for ``id_type`` for backwards
            compatibility.
        targetNameType : str
            Deprecated. Synonym for ``target_id_type``
            for backwards compatibility.
        """

        if not name:

            return set()

        id_type = id_type or nameType
        target_id_type = target_id_type or targetNameType

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
                        silent = silent,
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

            if target_id_type != 'uniprot':

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

        else:

            # all the other ID types
            mapped_names = self._map_name(
                name = name,
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = ncbi_tax_id,
            )

        # as ID translation tables for PRO IDs are not organism specific
        # we need an extra step to limit the results for the target organism
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

                if not mapped_names:

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

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        if isinstance(uniprots, common.basestring):

            uniprots = {uniprots}

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
            silent = True,
            nameType = None,
            targetNameType = None,
        ):
        """
        Same as ``map_name`` with multiple IDs.
        """

        return set.union(
            *(
                self.map_name(
                    name = name,
                    id_type = id_type,
                    target_id_type = target_id_type,
                    ncbi_tax_id = ncbi_tax_id,
                    strict = strict,
                    silent = silent,
                    nameType = nameType,
                    targetNameType = targetNameType,
                )
                for name in names
            )
        ) if names else set()


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
            ncbi_tax_id = ncbi_tax_id
        )

        return tbl[name] if tbl else set()

    #
    # ID specific translation methods
    #


    def label(self, name, id_type = None, ncbi_tax_id = None):
        """
        For any kind of entity, either protein, miRNA or protein complex,
        returns the preferred human readable label. For proteins this means
        Gene Symbols, for miRNAs miRNA names, for complexes a series of
        Gene Symbols.
        """

        if isinstance(name, common.list_like):

            return [
                self.label(
                    _name,
                    id_type = id_type,
                    ncbi_tax_id = ncbi_tax_id
                )
                for _name in name
            ]

        elif hasattr(name, 'genesymbol_str'):

            return name.genesymbol_str

        elif isinstance(name, common.basestring):

            ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

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

            else:

                return self.map_name0(
                    name,
                    id_type or 'uniprot',
                    'genesymbol',
                    ncbi_tax_id = ncbi_tax_id,
                ) or name

        else:

            return str(name)


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


    def primary_uniprot(self, uniprots):
        """
        For an iterable of UniProt IDs returns a set with the secondary IDs
        changed to the corresponding primary IDs. Anything what is not a
        secondary UniProt ID left intact.
        """

        primaries = set()

        for uniprot in uniprots:

            primary = self.map_name(
                name = uniprot,
                id_type = 'uniprot-sec',
                target_id_type = 'uniprot-pri',
                ncbi_tax_id = 0,
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

            if not swissprots:

                swissprots.add(uniprot)

            else:

                swissprots.update(this_swissprots)

        return swissprots


    def translate_deleted_uniprots_by_genesymbol(
            self,
            uniprots,
            ncbi_tax_id = None,
        ):

        if isinstance(uniprots, common.basestring):

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

        if uniprot_input.is_uniprot(uniprot, organism = ncbi_tax_id):

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

        if isinstance(uniprots, common.basestring):

            return self.valid_uniprot(uniprots, ncbi_tax_id = ncbi_tax_id)

        else:

            return {
                uniprot
                for uniprot in uniprots
                if uniprot_input.is_uniprot(uniprot, organism = ncbi_tax_id)
            }


    def valid_uniprot(self, uniprot, ncbi_tax_id = None):
        """
        If the UniProt ID ``uniprot`` exist in the proteome of the organism
        ``ncbi_tax_id`` returns the ID, otherwise returns None.
        """

        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id

        if uniprot_input.is_uniprot(uniprot, organism = ncbi_tax_id):

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

        if isinstance(uniprots, common.basestring):

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

    def has_mapping_table(
            self,
            id_type,
            target_id_type,
            ncbi_tax_id = None,
        ):

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

        self._log(
            'Loading mapping table for organism `%s` '
            'with identifiers `%s` and `%s`, '
            'input type `%s`' % (
                kwargs['ncbi_tax_id']
                    if 'ncbi_tax_id' in kwargs else
                resource.ncbi_tax_id,
                resource.id_type_a,
                resource.id_type_b,
                resource.type,
            )
        )

        reader = MapReader(param = resource, **kwargs)

        a_to_b = reader.mapping_table_a_to_b
        b_to_a = reader.mapping_table_b_to_a

        if a_to_b:

            self._log(
                '`%s` to `%s` mapping table for organism `%s` '
                'successfully loaded.' % (
                    resource.id_type_a,
                    resource.id_type_b,
                    str(resource.ncbi_tax_id),
                )
            )
            self.tables[a_to_b.get_key()] = a_to_b

        if b_to_a:

            self._log(
                '`%s` to `%s` mapping table for organism `%s` '
                'successfully loaded.' % (
                    resource.id_type_b,
                    resource.id_type_a,
                    str(resource.ncbi_tax_id),
                )
            )
            self.tables[b_to_a.get_key()] = b_to_a


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
            ncbi_tax_id=None,
        ):
        """
        Loads mapping tables from the huge static mapping file from UniProt.
        Takes long to download and process.
        """

        cachedir = cache_mod.get_cachedir()
        data = dict((key, collections.defaultdict(set)) for key in keys)
        cache_files = {}
        to_load = set()

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

                data[key] = pickle.load(open(cachefile, 'rb'))

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

            id_type_b = 'uniprot'

            for line in c.result:

                prg.step(len(line))

                line = line.decode('ascii').strip().split('\t')

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
                id_type = key.id_type,
                target_id_type = key.target_id_type,
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

            if table._expired():

                to_remove.add(key)

        for key in to_remove:

            self.remove_key(key)


    def __del__(self):

        if hasattr(self._mapper_cleanup_timeloop, 'stop'):

            for job in self._mapper_cleanup_timeloop.jobs:

                if job.is_alive():

                    job.stop()
                    job.stopped.set()


def init(**kwargs):

    if 'mapper' in globals():

        globals()['mapper'].__del__()

    globals()['mapper'] = Mapper(**kwargs)


def get_mapper(**kwargs):

    if 'mapper' not in globals():

        init(**kwargs)

    return globals()['mapper']


def map_name(
        name,
        id_type,
        target_id_type,
        ncbi_tax_id = None,
        strict = False,
        silent = True,
        expand_complexes = True,
        nameType = None,
        targetNameType = None,
    ):

    mapper = get_mapper()

    return mapper.map_name(
        name = name,
        id_type = id_type,
        target_id_type = target_id_type,
        ncbi_tax_id = ncbi_tax_id,
        strict = strict,
        silent = silent,
        expand_complexes = expand_complexes,
        nameType = nameType,
        targetNameType = targetNameType,
    )


def map_name0(
        name,
        id_type,
        target_id_type,
        ncbi_tax_id = None,
        strict = False,
        silent = True,
        nameType = None,
        targetNameType = None,
    ):

    mapper = get_mapper()

    return mapper.map_name0(
        name = name,
        id_type = id_type,
        target_id_type = target_id_type,
        ncbi_tax_id = ncbi_tax_id,
        strict = strict,
        silent = silent,
        nameType = nameType,
        targetNameType = targetNameType,
    )


def map_names(
        names,
        id_type = None,
        target_id_type = None,
        ncbi_tax_id = None,
        strict = False,
        silent = True,
        nameType = None,
        targetNameType = None,
    ):

    mapper = get_mapper()

    return mapper.map_names(
        names = names,
        id_type = id_type,
        target_id_type = target_id_type,
        ncbi_tax_id = ncbi_tax_id,
        strict = strict,
        silent = silent,
        nameType = nameType,
        targetNameType = targetNameType,
    )


def label(name, id_type = None, ncbi_tax_id = 9606):
    """
    For any kind of entity, either protein, miRNA or protein complex,
    returns the preferred human readable label. For proteins this means
    Gene Symbols, for miRNAs miRNA names, for complexes a series of
    Gene Symbols.
    """

    mapper = get_mapper()

    return mapper.label(
        name = name,
        id_type = id_type,
        ncbi_tax_id = ncbi_tax_id,
    )


def guess_type(name, entity_type = None):
    """
    From a string, tries to guess the ID type and optionally the entity
    type. Returns a tuple of strings: ID type and entity type.
    """

    mapper = get_mapper()

    return mapper.guess_type(name = name, entity_type = entity_type)


def id_from_label(label, label_id_type = 'genesymbol', ncbi_tax_id = None):
    """
    For a label (e.g. Gene Symbol) returns the corresponding IDs (e.g.
    UniProt IDs).
    """

    mapper = get_mapper()

    return mapper.id_from_label(
        label = label,
        label_id_type = label_id_type,
        ncbi_tax_id = ncbi_tax_id,
    )


def id_from_label0(label, label_id_type = 'genesymbol', ncbi_tax_id = None):
    """
    For a label (e.g. Gene Symbol) returns a single ID (e.g. UniProt IDs).
    """

    mapper = get_mapper()

    return mapper.id_from_label0(
        label = label,
        label_id_type = label_id_type,
        ncbi_tax_id = ncbi_tax_id,
    )