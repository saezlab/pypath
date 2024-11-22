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
Parse chemical table (SDF) files.
"""

from typing import IO
from collections.abc import Generator, Iterable

import importlib as imp
import sys
import os
import re
import collections
import itertools
import types
import warnings

from pypath_common import _misc
import pypath.share.session as session


resyn = re.compile(
    r'(^[A-Z]{2,})\(([0-9]+:[0-9]+)\(.*\)/([0-9]+:[0-9]+)\(.*\)\)'
)
rehg = re.compile(r'^([A-Z]{2,4})(\(.*\))')
refa = re.compile(r'C([0-9]+:[0-9]+)n?-?[-0-9]*$')
refa2 = re.compile(r'([0-9]{1,2}:[0-9])\(?[0-9EZ]*\)?$')
hgsyn = {
    'TG': 'TAG',
    'DG': 'DAG'
}


class SdfReader(session.Logger):

    names_default = {
        'PUBCHEM_CID': 'pubchem',
        'CHEBI_ID': 'chebi',
        'SYNONYMS': 'synonym',
        'INCHI': 'inchi',
        'INCHI_KEY': 'inchikey',
        'COMMON_NAME': 'commname',
        'SYSTEMATIC_NAME': 'sysname',
        'SMILES': 'smiles',
        'FORMULA': 'formula',
    }


    def __init__(
            self,
            sdf: str | IO,
            names: dict | None = None,
            fields: set | None = None,
            silent: bool = True,
            mol: bool = False,
        ):
        """
        Processes and serves data from an SDF file.

        Builds an index of the file and retrieve the records on demand.
        Note, sdf is not a well defined or well kept standard, this reader
        has been developed to process the LipidMaps database. Once there
        is a need to use with other databases we are happy to adapt to
        their formats.

        Args:
            sdf:
                Path or an open file pointer to the SDF file.
            names:
                These are the names to build indexes for. Once indexing is done
                it's possible to search and retrieve records by these IDs and
                names. By deafult the names in `names_default` are used. Names
                provided here are added to the defaults. Keys of the dict are
                labels as used in the sdf, values of the dict are the attribute
                names of the indexes.
            fields:
                Additional fields to be read. These are the data to be
                retrieved with each record. Works the same way as `names`.
            silent:
                Print number of records at the end of indexing.
            mol:
                Include the structure as chemical tab with the records.
        """

        session.Logger.__init__(self, name = 'sdf')

        self.sdf = sdf
        self.indexed = False
        self.silent = silent
        self.name_fields = names or {}
        self.name_fields.update(self.names_default)
        self.fields = fields or set()
        self.store_mol = mol

        self._empty()
        self._open()
        self._byte_mode()
        self._file_size()


    def reload(self):
        """
        Reload object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def _empty(self):

        self.data = {}
        self.synonym = collections.defaultdict(set)
        self.names = {
            name: collections.defaultdict(set)
            for name in itertools.chain(('id',), self.name_fields.values())
        }


    def _open(self):

        if isinstance(self.sdf, str):

            if os.path.exists(self.sdf):

                self._log(f'Opening SDF file `{self.sdf}`.')
                self.fp = open(self.sdf, 'rb')

            else:

                self._log(f'SDF file `{self.sdf}` does not exist.')
                raise FileNotFoundError(self.sdf)

        else:

            self.fp = self.sdf


    def _byte_mode(self):

        if not isinstance(self.fp.read(1), bytes):

            if 'b' not in self.fp.mode:

                self.fp.close()
                self.fp = open(self.fp.name, 'rb')

        self.fp.seek(0)
        self.fname = self.fp.name
        self._log(f'Opened SDF file `{self.fname}`.')


    def _file_size(self):

        self.fp.seek(-1, 2)
        self.eof = self.fp.tell()


    def parse(
            self,
            go_to: int = 0,
            store_mol: bool | None = None,
        ) -> Generator[tuple[dict, int]]:
        """
        Performs all reading operations on the sdf file.

        This method is able to read the entire file, scan the file and build
        an index of records, and retrieve one record.

        Args:
            go_to:
                Start reading from this byte offset in the file.
            store_mol:
                Include in the records also the structure as chemical tab.
        """

        if not isinstance(store_mol, bool):

            store_mol = self._store_mol(store_mol)

        self.fp.seek(go_to)

        # Parser state
        expect_new = True
        molpart = None
        namepart = False
        name_or_id = False
        this_offset = None
        offset = 0

        # Store parsed record
        _id = None
        source = None
        mol = ''
        name = {}
        annot = {}
        namekey = None


        for line in self.fp:

            llen = len(line)
            line = line.decode('utf-8').strip()

            if not molpart:

                if name_or_id and line and line[0] != '>' and line[0] != '$':

                    expect_new = True

                if namekey and namekey in self.name_fields:

                    name[namekey] = line
                    namekey = None

                if namekey and (not self.fields or namekey in self.fields):

                    annot[namekey] = line
                    namekey = None

                if line[:3] == '> <':

                    name_or_id = False
                    namepart = True
                    namekey = line[3:-1]

                if namepart and not line:

                    name_or_id = True

                if expect_new and llen:

                    # Start of new record
                    _id = line
                    name = {}
                    annot = {}
                    this_offset = offset
                    expect_new = False
                    molpart = 1
                    comment = ''

            elif molpart == 1:

                source = line
                molpart += 1

            elif molpart == 2 and line:

                if not line[0].isdigit():

                    comment += f' {line}'

                else:

                    molpart = 3

            expect_new = line == '$$$$'

            if molpart == 3:

                if store_mol:

                    mol += line

                if line == 'M  END':

                    molpart = None
                    namepart = True
                    name_or_id = True

            offset += llen

            if expect_new or self.fp.tell() == self.eof:

                yield (
                    {
                        'id': _id,
                        'source': source,
                        'comment': comment,
                        'mol': mol,
                        'name': name,
                        'annot': annot,
                    },
                    this_offset,
                )


    def __iter__(self):

        for record, _ in self.parse():

            yield record


    def load(
            self,
            store_mol: bool | None = None,
            reload: bool = False,
        ) -> None:
        """
        Loads the data from the SDF file into the dict under `data`.

        Args:
            store_mol:
                Include in the records also the structure as chemical tab.
            reload:
                Reload even if the `data` dict is already populated.
        """

        self.index(
            store = True,
            store_mol = store_mol,
            reindex = reload or not self.data,
        )


    def index(
            self,
            store: bool = False,
            store_mol: bool | None = None,
            reindex: bool = False,
        ) -> None:
        """
        Builds an index of the SDF file.

        Args:
            store:
                Store the parsed records in the dict under the `data`
                attribute.
            store_mol:
                Store also the structure as chemical tab.
            index_only:
                Do not read the file, only build an index.
        """

        if self.indexed and not reindex:

            return

        store_mol = self._store_mol(store_mol)

        self._empty()

        for rec, offset in self.parse(store_mol = store and store_mol):

            # this is indexing: we build dicts of names
            self.names['id'][rec['id']].add(offset)
            name = rec['name']

            if m := refa2.match(name.get('COMMON_NAME', '')):

                name['SYNONYMS'] = ';'.join(
                    _misc.to_list(name.get('SYNONYMS', None)) +
                    [f'FA({m.groups()[0]})']
                )

            for k, v in self.name_fields.items():

                if k not in name:

                    continue

                if k == 'SYNONYMS':

                    syns = {syn.strip() for syn in name[k].split(';')}

                    for rexp, templ in zip(
                        (rehg, resyn, refa),
                        ('%s%s', '%s(%s/%s)', 'FA(%s)')
                    ):

                        syns |= {
                            templ % ((hgsyn.get(g[0], g[0]),) + g[1:])
                            for syn in syns
                            if (
                                (m := rexp.match(syn)) and
                                (g := m.groups())
                            )
                        }

                    for syn in syns:

                        self.names['synonym'][syn].add(offset)

                    name['SYNONYMS'] = ';'.join(syns)

                else:

                    self.names[v][name[k]].add(offset)

            if store:

                self.data[offset] = rec

        self.index_info()
        self.indexed = True


    def _store_mol(self, store_mol: bool | None = None) -> bool:

        return (
            bool(self.store_mol)
                if not isinstance(store_mol, bool) else
            store_mol
        )


    def get_record(
            self,
            name: str | int,
            typ: str | None = None,
        ) -> list[dict]:
        """
        Retrieves all records matching `name`.

        Returns list of records or empty list if none found.
        Each record is a dict of processed values from the sdf file.

        Args
            name:
                Molecule name or identifier.
            typ:
                Type of name or identifier. These are the attribute names of
                the index dicts which are taken from the values in the `names`
                dict.
        """

        result = []

        if isinstance(name, int):

            result = [self.by_offset(name)]

        if (
            (index := self.names.get(typ or 'id', None)) and
            (offsets := index.get(name, None))
        ):

            result = list(map(self.by_offset, _misc.to_set(offsets)))

        return result


    def by_offset(self, offset: int) -> dict:
        """
        Get a record by its byte offset in the SDF file.
        """

        return self.data.get(offset, next(self.parse(go_to = offset)))


    def get_obmol(self, name, typ, use_mol = False):
        """
        Returns generator yielding `pybel.Molecule` instances for `name`.

        Args
            name:
                Molecule name or ID.
            typ:
                Type of the name or identifier.
            use_mol:
                Process structures from mol format.
                By default structures are processed from InChI.
        """

        records = self.get_record(name, typ)

        for rec in records:

            if use_mol:

                mol = self.record_to_obmol_mol(rec)

            else:

                mol = self.record_to_obmol(rec)

            mol.db_id = rec['id']
            title = []
            if 'COMMON_NAME' in rec['name']:
                title.append(rec['name']['COMMON_NAME'])
            if 'SYNONYMS' in rec['name']:
                title.extend(rec['name']['SYNONYMS'].split(';'))
            if 'SYSTEMATIC_NAME' in rec['name']:
                title.append(rec['name']['SYSTEMATIC_NAME'])
            mol.title = '|'.join(n.strip() for n in title)
            mol.lipidmaps = rec['id']
            if 'INCHI' in rec['name']:
                mol.inchi = rec['name']['INCHI']
            if 'PUBCHEM_CID' in rec['name']:
                mol.pubchem = rec['name']['PUBCHEM_CID']
            if 'CHEBI_ID' in rec['name']:
                mol.chebi = rec['name']['CHEBI_ID']
            if 'COMMON_NAME' in rec['name']:
                mol.name = rec['name']['COMMON_NAME']

            yield mol

    def _ensure_openbabel(self) -> types.ModuleType:

        def _log_and_warn(msg):

            self._log_traceback()
            self._log(msg)
            warnings.warn(msg)


        try:

            import openbabel.pybel as pybel

            if 'ipykernel' not in sys.modules and pybel.tk is None:

                try:

                    import tkinter
                    import PIL
                    import PIL.ImageTk
                    pybel.tk = tkinter
                    pybel.PIL = PIL.Image
                    pybel.piltk = PIL.ImageTk

                except Exception as e:

                    self.log_and_warn(
                        '`PIL` or `tkinter` not available. '
                        '`pybel` won\'t be able to draw molecules.'
                    )

        except Exception as e:
            self.log_and_warn('Module `pybel` not available.')


    def to_obmol(self, record):
        """
        Processes a record to `pybel.Molecule` object.
        """

        pybel = self._ensure_openbabel()

        if 'INCHI' in record['name']:

            return pybel.readstring('inchi', record['name']['INCHI'])

        else:

            sys.stdout.write(
                'No InChI for `%s`!\n' % record['name']['COMMON_NAME']
            )

    def record_to_obmol_mol(self, record):

        pybel = self._ensure_openbabel()

        return pybel.readstring('mol', self.get_mol(record))


    @staticmethod
    def format_mol(record):
        """
        Returns structure as a string in mol format.
        """

        return '%s\n  %s\n%s\n%s' % (
            record['id'],
            record['source'],
            record['comment'],
            record['mol']
        )


    def write_mol(self, name, typ, outf = None, return_data = False):
        """
        Writes a record into file in mol format.
        """

        outf = outf or '%s_%s_%s.mol'

        rr = self.get_record(name, typ)

        if not rr:

            return None

        if type(rr) is not list:

            rr = [rr]

        for r in rr:

            _outf = outf % (
                name.replace('/', '.'),
                r['name']['COMMON_NAME'].replace('/', '.').replace(' ', '..')
                    if 'COMMON_NAME' in r['name']
                    else '',
                r['id']
            )

            r['molfile'] = _outf

            with open(_outf, 'w') as fp:

                _ = fp.write(
                    self.get_mol(r)
                )

        if return_data:

            return rr


    def iter_records(self):
        """
        Iterates over all records in the SDF file.
        """

        for offset in self.names['id'].values():

            yield self.get_record(offset)


    def iter_obmol(self):
        """
        Iterates all structures in the file and yields `pybel.Molecule`
        objects.
        """

        for _id in self.names['id'].keys():

            for mol in self.get_obmol(_id, typ = 'id'):

                yield mol


    def index_info(self):
        """
        Prints number of records indexed and name of the source file.
        """

        msg = f"Indexed {len(self.names['id'])} records from `{self.fname}`."

        if not self.silent:

            sys.stdout.write(f'\t:: {msg}\n')

    def reopen(self):
        """
        Reopens the SDF file.
        """

        self.close()

        if hasattr(self, 'name'):

            self.sdf = self.fname
            self._open()
            self._byte_mode()


    def __del__(self):

        self.close()


    def close(self):
        """
        Close the SDF file.
        """

        if hasattr(self, 'fp') and hasattr(self.fp, 'close'):

            self._log(f'Closing SDF file `{self.fp.name}`.')
            self.fp.close()


    def __len__(self):

        return len(self.names['id'])


    def __repr__(self):

        return f'<SDF file `{self.fname}`: {len(self)} records>'


    def __contains__(self, name: str) -> bool:

        return any(name in n for n in self.names.values())


    def __getitem__(
            self,
            key: str | int | Iterable[str | int],
        ) -> dict | list[dict] | None:

        result = list(itertools.chain(*(

           self.get_record(key, typ)

           for key, typ in itertools.product(
               _misc.to_list(key),
               self.names.keys()
           )

        )))

        return result[0] if len(result) == 1 else result or None
