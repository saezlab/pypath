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
import itertools
import csv
import re
import asyncio
import inspect

from concurrent.futures.thread import ThreadPoolExecutor

from abc import ABC, abstractmethod
from typing import Iterable, Literal

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.common as common

_logger = session.Logger(name = 'kegg_api')
_log = _logger._log

_url = urls.urls['kegg_api']['url']


def _generate_relation_functions():

    _entity_types = ('disease', 'drug', 'gene', 'pathway')

    for etypes in itertools.combinations(_entity_types, 2):

        for args in (etypes, reversed(etypes)):

            args = tuple(args)
            name = f'{args[0]}_to_{args[1]}'
            synopsis = f'{args[0].capitalize()}-{args[1]} relations from KEGG.'

            def _relation_function(organism):

                if 'gene' in args:

                    args = args + (organism,)

                return _kegg_relations(*args)


            _relation_function.__name__ = name
            _relation_function.__doc__ = synopsis

            if 'gene' not in args:

                sig = inspect.signature(_relation_function)
                sig.replace(parameters = ())
                _relation_function.__signature__ = sig

            else:

                _relation_function.__doc__ += (
                    '\n\nArgs\n    organism:\n        Name of the organism. '
                    'Gene relations are organism specific.\n'
                )

            globals()[name] = _relation_function


def drug_to_drug(
    drugs: list | tuple | None = None,
    join: bool = True,
    asynchronous: bool = False
) -> dict[str, tuple]:
    """
    Downloads drug-drug interaction data from KEGG database.

    Args
        drugs:
            Drug IDs as a list or a tuple.
        join:
            If it's True, returns individual interactions of queried list.
            Else, joins them together and returns mutual interactions.
        asynchronous:
            Yet to be implemented.

    Returns
        A dict with disease IDs as keys and drug-drug interactions as values.
    """

    DrugToDrugInteraction = collections.namedtuple(
        'DrugToDrugInteraction',
        (
            'type',
            'name',
            'interactions',
        ),
    )

    Interaction = collections.namedtuple(
        'Interaction',
        (
            'type',
            'id',
            'name',
            'contraindication',
            'precaution',
        )
    )

    entry_types = {'d': 'drug', 'c': 'compound'}
    entry_dbs = {'drug': _Drug(), 'compound': _Compound()}
    interactions = collections.defaultdict(
        lambda: {
            'interactions': collections.defaultdict(list),
        }
    )

    join = join and (len(drugs) > 0)
    asynchronous = not drugs or asynchronous
    drugs = drugs or entry_dbs['drug'].data.keys()
    entries = _kegg_ddi(drugs, join = join, async_=asynchronous)

    for entry in entries:

        partners = dict(
            (
                role,
                {
                    'type': entry_types.get(entry[i][0].lower(), None),
                    'id': entry[i].split(':')[-1],
                    'name': (
                        entry_dbs[
                            entry_types.get(entry[i][0].lower(), None)
                        ].
                        get(entry[i].split(':')[-1], None)
                    ),
                }
            )
            for i, role in enumerate(('source', 'target'))
        )

        labels = entry[2].split(',')
        contraindication = 'CI' in labels
        precaution = 'P' in labels

        interaction = Interaction(
            type = partners['target']['type'],
            id = partners['target']['id'],
            name = partners['target']['name'],
            contraindication = contraindication,
            precaution = precaution,
        )

        disease_id = partners['source']['id']
        try:
            interactions[disease_id]['interactions'].append(interaction)
        except AttributeError:
            interactions[disease_id]['interactions'] = [interaction]
        interactions[disease_id]['type'] = partners['source']['type']
        interactions[disease_id]['name'] = partners['source']['name']

    interactions = dict(
        (
            key,
            DrugToDrugInteraction(
                value['type'],
                value['name'],
                tuple(value['interactions']),
            )
        )
        for key, value in interactions.items()
    )

    return interactions


def _generate_conv_functions():

    _id_types = (
        ('drug', ('chebi',)),
        ('gene', ('ncbi-geneid', 'uniprot')),
    )

    labels = {
        'chebi': 'ChEBI',
        'ncbi-geneid': 'NCBI Gene',
        'uniprot': 'UniProt',
    }

    for entity, id_types in _id_types:

        for id_type in id_types:

            args_ = (entity, id_type)

            for args in (args_, reversed(args_)):

                synopsis = (
                    'Translation dict between ' +
                    ' and '.join(
                        f'{labels.get(a, f"KEGG {a}")} IDs'
                        for a in args
                    ) +
                    '.'
                )

                def _conv_function(organism):

                    splits = [a != 'gene' for a in args]
                    args = [a if s else organism for s, a in zip(splits, args)]

                    return _kegg_conv(*args, *splits)


                name = (
                    '_to_'.join(
                        f'kegg_{a}' if a == entity else a
                        for a in args
                    ).
                    replace('-', '_')
                )
                _conv_function.__name__ = name
                _conv_function.__doc__ = synopsis

                if entity != 'gene':

                    sig = inspect.signature(_conv_function)
                    sig.replace(parameters = ())
                    _conv_function.__signature__ = sig

                else:

                    _conv_function.__doc__ += (
                        '\n\nArgs\n    organism:\n        Name of the '
                        'organism. Gene relations are organism specific.\n'
                    )

                globals()[name] = _conv_function


def _kegg_general(
    operation: str,
    *arguments: str,
) -> list[list[str]]:

    arguments = [arg for arg in arguments if arg is not None]

    url = '/'.join([_url % operation] + list(arguments))
    curl_args = {'url': url, 'silent': True, 'large': False}

    c = curl.Curl(**curl_args)

    lines = getattr(c, 'result', []).split('\n') or []

    return [line.split('\t') for line in lines if line]


async def _kegg_general_async(
    operation: str,
    *arguments: str,
) -> list[list[str]]:

    #TODO Yet to be implemented
    # This function doesn't work but it better
    # stay so we can implement it without
    # changing the structure of the module

    return _kegg_general(operation, *arguments)


def _kegg_list(
    database: str,
    option: str | None = None,
    organism: str | int | None = None,
) -> list[list[str]]:

    args = ['list', database]

    if database == 'brite' and option is not None:
        args += common.to_list(option)
    elif database == 'pathway' and organism is not None:
        args += common.to_list(organism)

    return _kegg_general(*args)


def _kegg_conv(
    source_db: str,
    target_db: str,
    source_split: bool = False,
    target_split: bool = False,
) -> dict[str, set[str]]:

    result = _kegg_general('conv', target_db, source_db)
    conversion_table = collections.defaultdict(set)

    for source, target in result:

        source = source.split(':')[1] if source_split else source
        target = target.split(':')[1] if target_split else target
        conversion_table[source].add(target)

    return dict(conversion_table)


def _kegg_link(source_db: str, target_db: str) -> list[list[str]]:

    return _kegg_general('link', target_db, source_db)


def _kegg_ddi(drug_ids: str | Iterable[str], join=True, async_: bool = False):


    if join and not isinstance(drug_ids, str):

        drug_ids = '+'.join(common.to_list(drug_ids))

    if async_:

        pool = ThreadPoolExecutor()

        return pool.submit(asyncio.run, _kegg_ddi_async(drug_ids)).result()

    return _kegg_ddi_sync(drug_ids)


def _kegg_ddi_sync(drug_ids: str | Iterable[str]):

    return list(itertools.chain(*(
        _kegg_general('ddi', drug_id)
        for drug_id in common.to_list(drug_ids)
    )))


async def _kegg_ddi_async(drug_ids):

    #TODO Yet to be implemented
    # This function doesn't work but it better
    # stay so we can implement it without
    # changing the structure of the module

    result = []

    for response in asyncio.as_completed([
        _kegg_general_async('ddi', drug_id)
        for drug_id in common.to_list(drug_ids)
    ]):
        the_response = await response
        result.extend(common.to_list(the_response))

    return result


def _kegg_relations(
    source_db: Literal['gene', 'pathway', 'disease', 'drug'],
    target_db: Literal['gene', 'pathway', 'disease', 'drug'],
    # should have human as a default, instead of triggering an error:
    organism: str | None = None,
) -> tuple:

    l_organism = common.to_list(organism)
    data = {}

    record = collections.namedtuple(
        'KeggEntry',
        (
            'id',
            'name',
            'type',
            'ncbi_gene_ids',
            'uniprot_ids',
            'chebi_ids',
        )
    )


    def get_data(name, cls_prefix = ''):

        if name not in data:

            cls = f'_{cls_prefix}{name.capitalize()}'
            data[name] = globals()[cls](*l_organism)

        return data[name]

    def db(name):

        return get_data(name)


    def ids(name):

        return get_data(name, cls_prefix = 'KeggTo')


    def process(entry, type_):

        id_ = db(type_).proc_key(entry)
        name = db(type_).get(id_, None)
        ncbi = ids('ncbi').get(id_) if type_ == 'gene' else ()
        uniprot = ids('uniprot').get(id_) if type_ == 'gene' else ()
        chebi = ids('chebi').get(id_) if type_ == 'drug' else ()

        return record(
            id = id_,
            name = name,
            type = type_,
            ncbi_gene_ids = ncbi,
            uniprot_ids = uniprot,
            chebi_ids = chebi,
        )

    args = [organism if db == 'gene' else db for db in (source_db, target_db)]
    entries = _kegg_link(*args)

    interactions = [(process(e[0], source_db), process(e[1], target_db)) for e in entries]

    return interactions


class _KeggDatabase(ABC):

    _data = None
    _query_args = None


    def __init__(self, *args):

        self.load(*args)


    @abstractmethod
    def proc_key(self, entry):

        return entry


    @abstractmethod
    def proc_value(self, entry):

        return entry


    def load(self, *args):

        entries = _kegg_list(*common.to_list(self._query_args), *args)

        self._data = {
            self.proc_key(entry[0]): self.proc_value(entry[1])
            for entry in entries
        }


    def get(self, index, default = None):

        return self._data.get(index, default)


    def __getitem__(self, index):

        return self.get(index)


    @property
    def data(self):

        return self._data


class _Organism(_KeggDatabase):

    _query_args = 'organism'

    def load(self, *args):

        entries = _kegg_list(*common.to_list(self._query_args), *args)

        self._data = {
            self.proc_key(entry[1]): self.proc_value(entry[0], entry[2])
            for entry in entries
        }


    def proc_value(self, entry):

        return self.get(entry)


    def proc_key(self, entry):

        return entry


class _Gene(_KeggDatabase):


    def __init__(self, organism):

        super().__init__(organism)

    def load(self, *args):

        entries = _kegg_list(*common.to_list(self._query_args), *args)

        self._data = {
            self.proc_key(entry[0]): self.proc_value(entry[-1])
            for entry in entries
        }


    def proc_key(self, entry):

        return entry


    def proc_value(self, entry):

        return entry.rsplit(';', maxsplit = 1)[-1].strip(' ')


class _Pathway(_KeggDatabase):

    _re_pathway = re.compile(r'\d+')
    _query_args = 'pathway'


    def proc_value(self, entry):

        return entry


    def proc_key(self, entry):

        pathway_id = self._re_pathway.search(entry)

        # is this correct?
        # there are pathway prefixes in KEGG other than "map"
        return f'map{pathway_id.group()}'


class _SplitDatabase(_KeggDatabase):


    def proc_key(self, entry):

        return entry[0].split(':')[1]


    def proc_value(self, entry):

        return entry[1]


class _Disease(_SplitDatabase):

    _query_args = 'disease'


class _Drug(_SplitDatabase):

    _query_args = 'drug'


class _Compound(_SplitDatabase):

    _query_args = 'compound'


class _ConversionTable:

    _table = {}


    def __init__(
        self,
        *id_types: str,
        source_split: bool = False,
        target_split: bool = False,
    ):

        self._id_types = id_types
        self._splits = {
            'source_split': source_split,
            'target_split': target_split,
        }
        self.load()


    @abstractmethod
    def load(self):

        self._table.update(_kegg_conv(*self._id_types, **self._splits))


    def get(self, index, default = None):

        return self._table.get(index, default)


    def __getitem__(self, index):

        return self._table.get(index, None)


    @property
    def table(self):

        return self._table


class _KeggToNcbi(_ConversionTable):


    def __init__(self, organism):

        super().__init__(organism, 'ncbi-geneid', target_split = True)


class _NcbiToKegg(_ConversionTable):


    def __init__(self, organism):

        super().__init__('ncbi-geneid', organism, source_split = True)


class _KeggToUniprot(_ConversionTable):


    def __init__(self, organism):

        super().__init__(organism, 'uniprot', target_split = True)


class _UniprotToKegg(_ConversionTable):


    def __init__(self, organism):

        super().__init__('uniprot', organism, source_split = True)


class _KeggToChebi(_ConversionTable):


    def __init__(self):

        super().__init__(
            'drug',
            'chebi',
            source_split = True,
            target_split = True,
        )


class _ChebiToKegg(_ConversionTable):


    def __init__(self):

        super().__init__(
            'chebi',
            'drug',
            source_split = True,
            target_split = True,
        )


_generate_relation_functions()
_generate_conv_functions()
