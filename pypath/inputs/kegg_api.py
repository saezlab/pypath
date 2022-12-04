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
#           Melih Darcan
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.organism/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.organism/
#

from __future__ import annotations

import collections
import csv
import re
import asyncio

from concurrent.futures.thread import ThreadPoolExecutor

from abc import ABC, abstractmethod
from collections.abc import Iterable

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.common as common

_logger = session.Logger(name = 'kegg_api')
_log = _logger._log

_url = urls.urls['kegg_api']['url']


def gene_to_pathway(organism):

    return _kegg_from_source_to_target('gene', 'pathway', organism)


def pathway_to_gene(organism):

    return _kegg_from_source_to_target('pathway', 'gene', organism)


def gene_to_drug(organism):

    return _kegg_from_source_to_target('gene', 'drug', organism)


def drug_to_gene(organism):

    return _kegg_from_source_to_target('drug', 'gene', organism)


def gene_to_disease(organism):

    return _kegg_from_source_to_target('gene', 'disease', organism)


def disease_to_gene(organism):

    return _kegg_from_source_to_target('disease', 'gene', organism)


def pathway_to_drug():

    return _kegg_from_source_to_target('pathway', 'drug')


def drug_to_pathway():

    return _kegg_from_source_to_target('drug', 'pathway')


def pathway_to_disease():

    return _kegg_from_source_to_target('pathway', 'disease')


def disease_to_pathway():

    return _kegg_from_source_to_target('disease', 'pathway')


def disease_to_drug():

    return _kegg_from_source_to_target('disease', 'drug')


def drug_to_disease():

    return _kegg_from_source_to_target('drug', 'disease')


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

    join = join and drugs
    asynchronous = not drugs or asynchronous
    drugs = drugs or entry_dbs['drug'].get_data().keys()
    entries = _kegg_ddi(drugs, join = join, asynchronous=asynchronous)

    for entry in entries:

        partners = dict(
            (
                role,
                {
                    'type': entry_types.get(entry[i][0].lower(), None),
                    'id': entry[i].split(':')[-1],
                    'name': entry_dbs[entry_type].get(entry_id, None)
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
        interactions[disease_id]['interactions'].append(interaction)
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


def kegg_gene_id_to_ncbi_gene_id(organism):

    return _kegg_conv(organism, 'ncbi-geneid', target_split=True)


def ncbi_gene_id_to_kegg_gene_id(organism):

    return _kegg_conv('ncbi-geneid', organism, source_split=True)


def kegg_gene_id_to_uniprot_id(organism):

    return _kegg_conv(organism, 'uniprot', target_split=True)


def uniprot_id_to_kegg_gene_id(organism):

    return _kegg_conv('uniprot', organism, source_split=True)


def kegg_drug_id_to_chebi_id():

    return _kegg_conv('drug', 'chebi', source_split=True, target_split=True)


def chebi_id_to_kegg_drug_id():

    return _kegg_conv('chebi', 'drug', source_split=True, target_split=True)


async def _kegg_general(
    operation: str,
    *arguments: str,
    async_: bool = False,
) -> list[list[str]]:

    url = '/'.join([_url % operation] + arguments)
    curl_args = {'url': url, 'silent': True, 'large': False}

    if _async:
        c = await curl.Curl(**curl_args)
    else:
        c = curl.Curl(**curl_args)

    lines = getattr(c, 'result', []) or []

    return [line.split('\t') for line in lines if line]


async def _kegg_general_async(
    operation: str,
    *arguments: str,
) -> list[list[str]]:

    #TODO Yet to be implemented
    # This function doesn't work but it better
    # stay so we can implement it without
    # changing the structure of the module

    return _kegg_general(operation, *arguments, async_ = False)


def _kegg_list(
    database: str,
    option: str | None = None,
    organism: str | int | None = None,
) -> list[list[str]]:

    args = (
        ['list', database] +
        common.to_list(option) if database == 'brite' else [] +
        common.to_list(organism) if database == 'pathway' else []
    )

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


def _kegg_ddi(drug_ids, join=True, asynchronous=False):

    if join and not isinstance(drug_ids, str):

        drug_ids = ['+'.join(drug_ids)]

    if asynchronous:

        pool = ThreadPoolExecutor()

        return pool.submit(asyncio.run, _kegg_ddi_async(drug_ids)).result()

    return _kegg_ddi_sync(drug_ids)


def _kegg_ddi_sync(drug_ids):

    result = list()

    if isinstance(drug_ids, Iterable):

        for drug_id in drug_ids:

            result.extend(_kegg_general('ddi', drug_id))

        return result


async def _kegg_ddi_async(drug_ids):

    #TODO Yet to be implemented
    # This function doesn't work but it better
    # stay so we can implement it without
    # changing the structure of the module

    result = []

    if isinstance(drug_ids, common.LIST_LIKE):

        for response in asyncio.as_completed([
            _kegg_general_async('ddi', drug_id)
            for drug_id in drug_ids
        ]):
            the_response = await response
            result.extend(common.to_list(the_response))

    return result


def _kegg_from_source_to_target(source_db, target_db, organism=None) -> tuple:

    db_name_list = [source_db, target_db]
    db_list = list()

    for db in db_name_list:

        if db == 'gene' and organism is not None:

            db_list.append(_Gene(organism))

            kegg_to_ncbi = _KeggToNcbi(organism)
            kegg_to_uniprot = _KeggToUniprot(organism)

        elif db == 'pathway':

            db_list.append(_Pathway())

        elif db == 'disease':

            db_list.append(_Disease())

        elif db == 'drug':

            db_list.append(_Drug())

            kegg_to_chebi = _KeggToChebi()

        else:
            print('Problem in function call. Check arguments.')
            exit()

    if target_db == 'gene':
        TargetDbEntry = collections.namedtuple(
            f'{target_db.capitalize()}Entry',
            [
                f'{target_db}_id',
                f'{target_db}_name',
                'ncbi_gene_id',
                'uniprot_ids'
            ]
        )

    elif target_db == 'drug':
        TargetDbEntry = collections.namedtuple(
            f'{target_db.capitalize()}Entry',
            [
                f'{target_db}_id',
                f'{target_db}_name',
                'chebi_id',
            ]
        )

    else:
        TargetDbEntry = collections.namedtuple(
            f'{target_db.capitalize()}Entry',
            [
                f'{target_db}_id',
                f'{target_db}_name',
            ]
        )


    if source_db == 'gene':
        Interaction = collections.namedtuple(
            f'{source_db.capitalize()}To{target_db.capitalize()}Interaction',
            [
                f'{source_db}_name',
                f'{target_db.capitalize()}Entries',
                'ncbi_gene_id',
                'uniprot_ids'
            ]
        )

    elif source_db == 'drug':
        Interaction = collections.namedtuple(
            f'{source_db.capitalize()}To{target_db.capitalize()}Interaction',
            [
                f'{source_db}_name',
                f'{target_db.capitalize()}Entries',
                'chebi_id'
            ]
        )

    else:
        Interaction = collections.namedtuple(
            f'{source_db.capitalize()}To{target_db.capitalize()}Interaction',
            [
                f'{source_db}_name',
                f'{target_db.capitalize()}Entries',
            ]
        )

    source = db_list[0]
    target = db_list[1]

    source_url = source_db if source_db != 'gene' else org
    target_url = target_db if target_db != 'gene' else org

    entries = _kegg_link(source_url, target_url)
    interactions = collections.defaultdict(
        lambda: collections.defaultdict(list)
    )

    for entry in entries:

        source_id = source.handle(entry[0])
        target_id = target.handle(entry[1])

        source_name = source.get(source_id, None)
        target_name = target.get(target_id, None)

        if target_db == 'gene':

            ncbi_gene_id = kegg_to_ncbi.get(target_id)
            uniprot_ids = kegg_to_uniprot.get(target_id)

            if isinstance(ncbi_gene_id, list):
                ncbi_gene_id = tuple(ncbi_gene_id)

            if isinstance(uniprot_ids, list):
                uniprot_ids = tuple(uniprot_ids)

            target_db_entry = TargetDbEntry(
                target_id,
                target_name,
                ncbi_gene_id,
                uniprot_ids
            )

        elif target_db == 'drug':

            chebi_id = kegg_to_chebi.get(target_id)

            if isinstance(chebi_id, list):
                chebi_id = tuple(chebi_id)

            target_db_entry = TargetDbEntry(
                target_id,
                target_name,
                chebi_id
            )

        else:

            target_db_entry = TargetDbEntry(
                target_id,
                target_name,
            )

        interactions[source_id][f'{target_db}_entries'].append(target_db_entry)
        interactions[source_id][f'{source_db}_name'] = source_name

        if source_db == 'gene':

            interactions[source_id]['ncbi_gene_id'] = (
                kegg_to_ncbi.get(source_id)
            )
            interactions[source_id]['uniprot_ids'] = (
                kegg_to_uniprot.get(source_id)
            )

        elif source_db == 'drug':

            interactions[source_id]['chebi_id'] = (
                kegg_to_chebi.get(source_id)
            )

    for key, value in interactions.items():

        if source_db == 'gene':
            interaction = Interaction (
                value[f'{source_db}_name'],
                tuple(value[f'{target_db}_entries']),
                value['ncbi_gene_id'],
                value['uniprot_ids']
            )

        elif source_db == 'drug':
            interaction = Interaction (
                value[f'{source_db}_name'],
                tuple(value[f'{target_db}_entries']),
                value['chebi_id'],
            )

        else:
            interaction = Interaction (
                value[f'{source_db}_name'],
                tuple(value[f'{target_db}_entries'])
            )

        interactions[key] = interaction

    if organism is not None:
        organism = _Organism()
        org_id, org_name = organism.get(organism)
        interactions['org_id'] = org_id
        interactions['org_name'] = org_name

    return interactions


class _KeggDatabase(ABC):

    _data = None


    @abstractmethod
    def __init__(self):
        pass


    @abstractmethod
    def handle(self):
        pass


    @abstractmethod
    def download_data(self):
        pass


    def get(self, index):
        return self._data[index]


    def get_data(self):
        return self._data


class _Organism(_KeggDatabase):

    def __init__(self):
        self.download_data()


    def download_data(self):
        entries = _kegg_list('organism')
        self._data = {self.handle(organism) : [org_id, org_name] for (org_id, organism, org_name, _) in entries}


    def handle(self, organism):
        return org


class _Gene(_KeggDatabase):

    def __init__(self, organism):
        self.download_data(organism)


    def download_data(self, organism):

        entries = _kegg_list(organism)

        gene_slice = [row[0] for row in entries]

        name_slice = [row[-1] for row in entries]
        name_slice = [name.split(';')[-1] for name in name_slice]
        name_slice = [name.strip(' ') for name in name_slice]

        entries = zip(gene_slice, name_slice)
        self._data = {self.handle(gene) : gene_name for (gene, gene_name) in entries}


    def handle(self, gene):

        return gene


class _Pathway(_KeggDatabase):

    def __init__(self, organism=None):
        self.download_data()

    def download_data(self, organism=None):

        if organism is not None:

            entries = _kegg_list('pathway', organism)

        else:

            entries = _kegg_list('pathway')

        self._data = {self.handle(pathway) : pathway_name for (pathway, pathway_name) in entries}

    def handle(self, pathway):

        pathway_re = re.compile(r'\d+')
        pathway_id = pathway_re.search(pathway)

        return 'map' + pathway_id.group()


class _SplitDatabase(_KeggDatabase):
    def __init__(self, entry_url):
        self.download_data(entry_url)


    def download_data(self, entry_url):
        entries = _kegg_list(entry_url)
        self._data = {self.handle(entry) : entry_name for (entry, entry_name) in entries}


    def handle(self, entry):
        return entry.split(':')[1]


class _Disease(_SplitDatabase):

    def __init__(self):
        super().__init__('disease')


class _Drug(_SplitDatabase):

    def __init__(self):
        super().__init__('drug')


class _Compound(_SplitDatabase):

    def __init__(self):
        super().__init__('compound')


class _ConversionTable:

    _table = dict()

    def __init__(self):
        self.download_table()


    @abstractmethod
    def download_table(self):
        pass


    def get(self, index):
        try:
            return self._table[index]
        except KeyError:
            return None


    def get_table(self):
        return self._table


class _OrgTable(_ConversionTable):

    def __init__(self, organism=None):
        if organism is not None:
            self.download_table(organism)


class _KeggToNcbi(_OrgTable):

    def download_table(self, organism):
        table = _kegg_conv(organism, 'ncbi-geneid', target_split=True)
        self._table.update(table)


class _NcbiToKegg(_OrgTable):

    def download_table(self, organism):
        table = _kegg_conv('ncbi-geneid', organism, source_split=True)
        self._table.update(table)


class _KeggToUniprot(_OrgTable):

    def download_table(self, organism):
        table = _kegg_conv(organism, 'uniprot', target_split=True)
        self._table.update(table)


class _UniprotToKegg(_OrgTable):

    def download_table(self, organism):
        table = _kegg_conv('uniprot', organism, source_split=True)
        self._table.update(table)


class _KeggToChebi(_ConversionTable):

    def download_table(self):
        table = _kegg_conv('drug', 'chebi', source_split=True, target_split=True)
        self._table = table


class _ChebiToKegg(_ConversionTable):

    def download_table(self):
        table = _kegg_conv('chebi', 'drug', source_split=True, target_split=True)
        self._table = table
