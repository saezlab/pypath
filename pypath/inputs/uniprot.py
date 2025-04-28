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

from typing import Iterable

import re
import json
import collections
import itertools
import functools
import urllib.parse

import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.settings as settings
import pypath.share.session as session_mod
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.utils.taxonomy as taxonomy
from pypath.inputs.uniprot_idmapping import idtypes as idmapping_idtypes

_logger = session_mod.Logger(name = 'uniprot_input')

_redatasheet = re.compile(r'([A-Z\s]{2})\s*([^\n\r]+)[\n\r]+')

# regex for matching UniProt AC format
# from https://www.uniprot.org/help/accession_numbers
reac = re.compile(
    r'[OPQ][0-9][A-Z0-9]{3}[0-9]|'
    r'[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
)
_rename = re.compile(r'Name=([\w\(\)-]+)\W')
_retaxid = re.compile(r'=(\d+)[^\d]')


def _all_uniprots(organism = 9606, swissprot = None):

    swissprot = _swissprot_param(swissprot)
    rev = '' if swissprot is None else ' AND reviewed: %s' % swissprot
    url = urls.urls['uniprot_basic']['url']
    get = {
        'query': 'organism_id:%s%s' % (str(organism), rev),
        'format': 'tsv',
        'fields': 'accession',
    }

    if organism == '*':
        get['query'] = rev.strip(' AND ')

    c = curl.Curl(url, get = get, silent = False, slow = True)
    data = c.result

    return {
        l.strip() for l in data.split('\n')[1:] if l.strip()
    }


def _swissprot_param(swissprot):

    return (
        'true'
            if swissprot in {'true', 'True', 'yes', 'YES', True} else
        'false'
            if swissprot in {'false', 'False', 'no', 'NO', False} else
        None
    )


def valid_uniprot(name):
    """
    Checks if ``name`` fits the format requirements for UniProt accession
    numbers.
    """

    return bool(reac.match(name))


def protein_datasheet(identifier):

    url = urls.urls['uniprot_basic']['datasheet'] % identifier.strip()

    datasheet =  _protein_datasheet(url)

    if not datasheet:

        _logger._log(
            'UniProt ID `%s` returns empty response, it might be and an old '
            'ID which has been deleted from the database. Attempting to '
            'find its history and retrieve either an archived version or '
            'the find the new ID which replaced this one.' % identifier
        )
        return uniprot_history_recent_datasheet(identifier)

    else:

        return datasheet


def deleted_uniprot_genesymbol(identifier):
    """
    Retrieves the archived datasheet for a deleted UniProt ID and returns
    the Gene Symbol and the NCBI Taxonomy ID from the datasheet.
    """

    datasheet = uniprot_history_recent_datasheet(identifier)
    genesymbol = None
    ncbi_tax_id = None

    for tag, line in datasheet:

        if tag == 'GN':

            m = _rename.search(line.strip())

            if m:

                genesymbol = m.groups()[0]

        if tag == 'OX':

            ncbi_tax_id = int(_retaxid.search(line).groups()[0])
            break

    return genesymbol, ncbi_tax_id


def _protein_datasheet(url):

    cache = True

    for a in range(3):

        c = curl.Curl(
            url,
            silent = True,
            large = False,
            cache = cache,
            connect_timeout = (
                settings.get('uniprot_datasheet_connect_timeout')
            ),
            timeout = settings.get('uniprot_datasheet_timeout'),
        )

        if not c.result or c.result.startswith('<!DOCTYPE'):

            cache = False

        else:

            break

    if not c.result:

        _logger._log(
            'Could not retrieve UniProt datasheet by URL `%s`.' % url
        )

    return _redatasheet.findall(c.result) if c.result else []


def uniprot_history_recent_datasheet(identifier):

    recent_version = uniprot_recent_version(identifier)

    if recent_version:

        if recent_version.replaced_by:

            new = recent_version.replaced_by.split(';')[0]
            url = urls.urls['uniprot_basic']['datasheet'] % new
            _logger._log(
                'UniProt ID `%s` is obsolete, has been replaced by '
                '`%s`: `%s`.' % (
                    identifier,
                    new,
                    url,
                )
            )
            return protein_datasheet(new)

        else:

            version = int(recent_version.entry_version)
            url = '%s?version=%u' % (
                urls.urls['uniprot_basic']['datasheet'] % identifier,
                version,
            )
            _logger._log(
                'UniProt ID `%s` is obsolete, downloading archived '
                'version %u: `%s`.' % (
                    identifier,
                    version,
                    url,
                )
            )
            c = curl.Curl(url, silent = True, large = False)
            return _protein_datasheet(url)

    return []


UniprotRecordHistory = collections.namedtuple(
    'UniprotRecordHistory',
    [
        'entry_version',
        'sequence_version',
        'entry_name',
        'database',
        'number',
        'date',
        'replaces',
        'replaced_by',
    ],
)


def uniprot_history(identifier):
    """
    Retrieves the history of a record.
    Returns a generator iterating over the history from most recent to the
    oldest.
    """

    if valid_uniprot(identifier):

        url_history = urls.urls['uniprot_basic']['history'] % identifier
        c_history = curl.Curl(
            url_history,
            silent = True,
            large = True,
        )

        if c_history.result:

            line0 = next(c_history.result)

            if not line0.startswith('<!DOCTYPE'):

                for line in c_history.result:

                    if line:

                        yield UniprotRecordHistory(
                            *(
                                field.strip() for field in line.split('\t')
                            )
                        )


def uniprot_recent_version(identifier):

    for version in uniprot_history(identifier):

        if (
            (
                version.entry_version != '0' and
                version.entry_name != 'null'
            ) or version.replaced_by
        ):

            return version


def uniprot_deleted(confirm = True):

    return swissprot_deleted() | trembl_deleted(confirm = confirm)


def _uniprot_deleted(swissprot = True, confirm = True):

    if not swissprot and confirm:

        resp = input(
            'Loading the list of deleted TrEMBL IDs requires '
            '>5GB memory. Do you want to proceed [y/n] '
        )

        if not resp or resp[0].lower() != 'y':

            return set()

    key = 'deleted_%s' % ('sp' if swissprot else 'tr')
    url = urls.urls['uniprot_basic'][key]
    c = curl.Curl(url, silent = False, large = True)

    result = set()

    for line in c.result:

        m = reac.match(line.strip())

        if m:

            result.add(m.groups()[0])

    return result


def swissprot_deleted():

    return _uniprot_deleted(swissprot = True)


def trembl_deleted(confirm = True):

    return _uniprot_deleted(swissprot = False, confirm = True)


def get_uniprot_sec(organism = 9606):
    """
    Downloads and processes the mapping between secondary and
    primary UniProt IDs.

    Yields pairs of secondary and primary UniProt IDs.

    :param int organism:
        NCBI Taxonomy ID of the organism.
    """

    _organism = organism not in (None, _const.NOT_ORGANISM_SPECIFIC)

    if _organism:

        from pypath.inputs import uniprot_db
        proteome = uniprot_db.all_uniprots(organism=organism)
        proteome = set(proteome)

    sec_pri = []
    url = urls.urls['uniprot_sec']['url']
    c = curl.Curl(url, silent = False, large = True, timeout = 2400)

    for i, line in enumerate(c.result):

        if i < 30:

            continue

        line = line.split()

        if len(line) == 2 and (not _organism or line[1] in proteome):

            yield line


class UniprotQuery:

    _PROCESS = {
        'dict': '_process_dict',
        'list': '_process_list',
    }
    _OP = ('_AND', '_NOT', '_OR')
    _OPSTART = re.compile(r'^(OR|AND)')
    _OPEND = re.compile(r'(OR|AND)$')
    _FIELDSEP = re.compile(r'[\s;]')
    _FIELDEND = re.compile(r';$')
    _SYNONYMS = {
        'organism': 'organism_id',
        'ncbi_tax_id': 'organism_id',
    }
    _FIELD_SYNONYMS = {
        'function': 'cc_function',
        'activity_regulation': 'cc_activity_regulation',
        'tissue_specificity': 'cc_tissue_specificity',
        'developmental_stage': 'cc_developmental_stage',
        'induction': 'cc_induction',
        'intramembrane': 'ft_intramem',
        'signal_peptide': 'ft_signal',
        'subcellular_location': 'cc_subcellular_location',
        'transmembrane': 'ft_transmem',
        'comment': 'cc_miscellaneous',
        'topological_domain': 'ft_topo_dom',
        'family': 'protein_families',
        'interactor': 'cc_interaction',
        'keywords': 'keyword',
    }


    def __init__(
            self,
            *query,
            fields: str | Iterable[str] | None = None,
            **kwargs
        ):
        """
        Constructs a query for the UniProt REST API.

        Args:
            query:
                Query elements: can be a ready query or its components, bypassing
                the processing in this function or performing only simple
                concatenation. Alternatively, it can be a nested structure of lists
                and dicts describing more complex queries. See the examples below.
            kwargs:
                Same as passing a dict to ``query``.

        Details:
            The query can be built in several ways:
            - Simple string or concatenation of strings:
              query_builder('kinase AND organism_id:9606')
              query_builder('kinase', 'organism_id:9606')
              query_builder('kinase', organism_id = 9606)
              The above 3 examples all return the same query:
              `kinase AND organism_id:9606`
            - The default operator within lists is `OR` and within dicts is `AND`:
              query_builder(organism = [9606, 10090, 10116])
              `organism_id:9606 OR organism_id:10090 OR organism_id:10116`
              query_builder({'organism_id': 9606, 'reviewed': True})
              `organism_id:9606 AND reviewed:true`
            - These default operators can be changed by including the `op` key in
              dicts or including the operator with an underscore in lists:
              query_builder({'length': (500,), 'mass': (50000,), 'op': 'OR'})
              `length:[500 TO *] OR mass:[50000 TO *]`
              query_builder(lit_author = ['Huang', 'Kovac', '_AND'])
              `lit_author:Huang AND lit_author:Kovac`
            - The nested structures translate into nested parentheses in the query:
              query_builder({'organism_id': [9606, 10090], 'reviewed': True})
              `(organism_id:9606 OR organism_id:10090) AND reviewed:true`
            - Values are converted to strings, intervals can be provided as tuples:
              query_builder({'length': (100, None), 'organism_id': 9606})
              `length:[100 TO *] AND organism_id:9606`

            For a complete reference of the available parameters, see
            https://www.uniprot.org/help/query-fields and
            https://www.uniprot.org/help/text-search for additional syntax
            elements.

            For the available fields refer to the ``_FIELD_SYNONYMS`` attribute of
            this class or the UniProt website:
            https://www.uniprot.org/help/return_fields

        Methods:
            __iter__:
                Perform the query and iterate over the lines in the results,
                skipping the header and the empty lines, stripping the
                linebreaks and splitting by tab.

                Yields:
                    A list of fields for each line.

        Attributes:
            fail_on_empty:
                If set to True, an error will be raised if the UniProt API
                returns empty response. By default no error is raised.
            name_process:
                If set to True, a different processing will be applied on the
                results. This is appropriate especially for identifier type
                fields.
        """

        self.fields = common.to_list(fields)
        self._args = query, kwargs
        self._process_main()
        # tolerate empty result: Curl returns None in case of
        # empty file but in case of UniProt, especially for under-researched
        # taxons it can happen there is no result for certain queries
        self.fail_on_empty = False
        self.name_process = False


    @classmethod
    def _value(
            cls,
            val: str | int | bool | tuple,
            field: str | None = None,
        ) -> str:

        field = cls._SYNONYMS.get(field, field)

        if field == 'organism_id':

            result = str(taxonomy.ensure_ncbi_tax_id(val) or val)

        elif isinstance(val, tuple):

            aux = tuple(map(cls._value, val))

            if len(aux) > 1 or field != 'organism_id':

                val = (aux + ('*',))[:2]
                result = '[%s TO %s]' % val

            else:

                result = aux[0]

        elif val is None:

            if field in ('reviewed', 'organism_id'):

                result = ''
                field = None

            else:

                result = '*'

        elif isinstance(val, bool):

            result = str(val).lower()

        else:

            result = str(val)

        if field:

            result = f'{field}:{result}'

        return result


    def _process_main(self):

        query, kwargs = self._args
        op = kwargs.pop('_op', 'AND')
        query = list(query)
        query.append(kwargs)
        result = []

        for q in query:

            q = self._process(q).strip()

            if (
                result and
                q and
                not self._OPEND.match(result[-1]) and
                not self._OPSTART.match(q)
            ):

                result.append(op)

            if q:

                result.append(q)

        self.query = ' '.join(result)


    @classmethod
    def _process(
            cls,
            query: str | list | dict,
            field: str | None = None,
        ) -> str:

        method = cls._PROCESS.get(type(query).__name__, '_value')

        return getattr(cls, method)(query, field)


    @classmethod
    def _process_list(cls, query: list, field: str | None = None) -> str:

        op = '_OR'

        for _op in cls._OP:

            if _op in query:

                op = query.pop(query.index(_op))

        op = f' {op[1:]} '

        query = [cls._process(i, field) for i in query]

        return cls._par(op.join(query))


    @classmethod
    def _process_dict(cls, query: dict, field: str | None = None) -> str:

        query = query.copy()
        op = ' %s ' % query.pop('op', ' AND ').strip()
        result = op.join(
            it for k, v in query.items()
            if (it := cls._process(v, k))
        )

        return cls._par(result) if len(query) > 1 else result


    @staticmethod
    def _par(value: str) -> str:

        return f'({value})' if value else ''


    @property
    def _get(self) -> dict[str, str]:

        field_qs = ','.join(
            ['accession'] +
            [self._FIELD_SYNONYMS.get(f, f) for f in self.fields]
        )

        return {
            'query': self.query,
            'format': 'tsv',
            'fields': field_qs,
            'compressed': 'true',
        }


    @property
    def _baseurl(self) -> str:

        return urls.urls['uniprot_basic']['url']


    @property
    def url(self) -> str:
        """
        UniProt REST API URL (urlencoded).

        Returns:
            A valid query suitable for the UniProt REST API.
        """

        return f'{self._baseurl}?{urllib.parse.urlencode(self._get)}'


    @property
    def url_plain(self) -> str:
        """
        UniProt REST API URL (plain).
        """

        return urllib.parse.unquote_plus(self.url)


    def __iter__(self):

        c = curl.Curl(
            self._baseurl,
            get = self._get,
            silent = False,
            large = True,
            compr = 'gz',
            slow = True,
        )
        result = c.result if c.result or self.fail_on_empty else [0].__iter__()
        _ = next(result)
        _proc0 = functools.partial(self._FIELDEND.sub, '')
        _proc1 = self._FIELDSEP.split if self.name_process else common.identity

        for line in result:

            line = line.strip('\n\r')

            if line.strip():

                yield [_proc1(_proc0(f)) for f in line.split('\t')]


    def perform(self) -> list[str] | dict[str, str] | dict[str, dict[str, str]]:
        """
        Perform the query and preprocess the result.

        Returns:
            - A list of UniProt IDs if no fields were provided.
            - A dict of UniProt IDs and corresponding field values if
              exactly one field was provided.
            - A dict with field names as top level keys and dicts of the
              kind described in the previous point as values.
        """

        _id, *variables = zip(*self)
        _id = list(map(common.sfirst, _id))

        if variables:

            result = {
                f: {i: v for i, v in zip(_id, vs) if i}
                for f, vs in zip(self.fields, variables)
            }

            result = (
                common.first(result.values())
                    if len(result) == 1 else
                result
            )

        else:

            result = list(_id)

        return result


def query_builder(*query, **kwargs) -> str:
    """
    Build a query for the UniProt web site and REST API.

    Args:
        query:
            Query elements: can be a ready query or its components, bypassing
            the processing in this function or performing only simple
            concatenation. Alternatively, it can be a nested structure of lists
            and dicts describing more complex queries. See the examples below.
        kwargs:
            Same as passing a dict to ``query``.

    Details:
        The query can be built in several ways:
        - Simple string or concatenation of strings:
          query_builder('kinase AND organism_id:9606')
          query_builder('kinase', 'organism_id:9606')
          query_builder('kinase', organism_id = 9606)
          The above 3 examples all return the same query:
          `kinase AND organism_id:9606`
        - The default operator within lists is `OR` and within dicts is `AND`:
          query_builder(organism = [9606, 10090, 10116])
          `organism_id:9606 OR organism_id:10090 OR organism_id:10116`
          query_builder({'organism_id': 9606, 'reviewed': True})
          `organism_id:9606 AND reviewed:true`
        - These default operators can be changed by including the `op` key in
          dicts or including the operator with an underscore in lists:
          query_builder({'length': (500,), 'mass': (50000,), 'op': 'OR'})
          `length:[500 TO *] OR mass:[50000 TO *]`
          query_builder(lit_author = ['Huang', 'Kovac', '_AND'])
          `lit_author:Huang AND lit_author:Kovac`
        - The nested structures translate into nested parentheses in the query:
          query_builder({'organism_id': [9606, 10090], 'reviewed': True})
          `(organism_id:9606 OR organism_id:10090) AND reviewed:true`
        - Values are converted to strings, intervals can be provided as tuples:
          query_builder({'length': (100, None), 'organism_id': 9606})
          `length:[100 TO *] AND organism_id:9606`

        For a complete reference of the available parameters, see
        https://www.uniprot.org/help/query-fields and
        https://www.uniprot.org/help/text-search for additional syntax
        elements.

    Returns:
        A query that can be inserted into the UniProt search field.
    """

    return UniprotQuery(*query, **kwargs).query


def uniprot_data(
        *query,
        fields: str | Iterable[str] | None = None,
        organism: str | int | None = 9606,
        reviewed: bool | None = True,
        **kwargs
    ) -> dict[str, str] | dict[str, dict[str, str]]:
    """
    Basic client for the UniProt REST API.

    Retrieves one or more fields from UniProt, by default for all reviewed
    (SwissProt) proteins of one organism

    Args:
        query:
            Query elements: can be a ready query or its components, bypassing
            the processing in this function or performing only simple
            concatenation. Alternatively, it can be a nested structure of lists
            and dicts describing more complex queries. See the examples below.
        fields:
            One or more UniProt field name. See details.
        organism:
            Organism name or identifier, e.g. "human", or "Homo sapiens",
            or 9606.
        reviewed:
            Restrict the query to SwissProt (True), to TrEMBL (False), or
            cover both (None).
        kwargs:
            Same as passing a dict to ``query``.

    Details:
        The query can be built in several ways:
        - Simple string or concatenation of strings:
          query_builder('kinase AND organism_id:9606')
          query_builder('kinase', 'organism_id:9606')
          query_builder('kinase', organism_id = 9606)
          The above 3 examples all return the same query:
          `kinase AND organism_id:9606`
        - The default operator within lists is `OR` and within dicts is `AND`:
          query_builder(organism = [9606, 10090, 10116])
          `organism_id:9606 OR organism_id:10090 OR organism_id:10116`
          query_builder({'organism_id': 9606, 'reviewed': True})
          `organism_id:9606 AND reviewed:true`
        - These default operators can be changed by including the `op` key in
          dicts or including the operator with an underscore in lists:
          query_builder({'length': (500,), 'mass': (50000,), 'op': 'OR'})
          `length:[500 TO *] OR mass:[50000 TO *]`
          query_builder(lit_author = ['Huang', 'Kovac', '_AND'])
          `lit_author:Huang AND lit_author:Kovac`
        - The nested structures translate into nested parentheses in the query:
          query_builder({'organism_id': [9606, 10090], 'reviewed': True})
          `(organism_id:9606 OR organism_id:10090) AND reviewed:true`
        - Values are converted to strings, intervals can be provided as tuples:
          query_builder({'length': (100, None), 'organism_id': 9606})
          `length:[100 TO *] AND organism_id:9606`

        For a complete reference of the available parameters, see
        https://www.uniprot.org/help/query-fields and
        https://www.uniprot.org/help/text-search for additional syntax
        elements.

        For the available fields refer to the ``_FIELD_SYNONYMS`` attribute of
        the UniprotQuery class or the UniProt website:
        https://www.uniprot.org/help/return_fields

    Returns:
        - A list of UniProt IDs if no fields were provided.
        - A dict of UniProt IDs and corresponding field values if
          exactly one field was provided.
        - A dict with field names as top level keys and dicts of the
          kind described in the previous point as values.
    """

    for arg in ('organism', 'reviewed'):

        if locals()[arg] is not None:

            kwargs[arg] = locals()[arg]

    return uniprot_query(*query, fields = fields, **kwargs)


def uniprot_query(
        *query,
        fields: str | Iterable[str] | None = None,
        **kwargs
    ) -> dict[str, str] | dict[str, dict[str, str]]:
    """
    Basic client for the UniProt REST API.

    Args:
        query:
            Query elements: can be a ready query or its components, bypassing
            the processing in this function or performing only simple
            concatenation. Alternatively, it can be a nested structure of lists
            and dicts describing more complex queries. See the examples below.
        fields:
            One or more UniProt field name. See details.
        kwargs:
            Same as passing a dict to ``query``.

    Details:
        The query can be built in several ways:
        - Simple string or concatenation of strings:
          query_builder('kinase AND organism_id:9606')
          query_builder('kinase', 'organism_id:9606')
          query_builder('kinase', organism_id = 9606)
          The above 3 examples all return the same query:
          `kinase AND organism_id:9606`
        - The default operator within lists is `OR` and within dicts is `AND`:
          query_builder(organism = [9606, 10090, 10116])
          `organism_id:9606 OR organism_id:10090 OR organism_id:10116`
          query_builder({'organism_id': 9606, 'reviewed': True})
          `organism_id:9606 AND reviewed:true`
        - These default operators can be changed by including the `op` key in
          dicts or including the operator with an underscore in lists:
          query_builder({'length': (500,), 'mass': (50000,), 'op': 'OR'})
          `length:[500 TO *] OR mass:[50000 TO *]`
          query_builder(lit_author = ['Huang', 'Kovac', '_AND'])
          `lit_author:Huang AND lit_author:Kovac`
        - The nested structures translate into nested parentheses in the query:
          query_builder({'organism_id': [9606, 10090], 'reviewed': True})
          `(organism_id:9606 OR organism_id:10090) AND reviewed:true`
        - Values are converted to strings, intervals can be provided as tuples:
          query_builder({'length': (100, None), 'organism_id': 9606})
          `length:[100 TO *] AND organism_id:9606`

        For a complete reference of the available parameters, see
        https://www.uniprot.org/help/query-fields and
        https://www.uniprot.org/help/text-search for additional syntax
        elements.

        For the available fields refer to the ``_FIELD_SYNONYMS`` attribute of
        the UniprotQuery class or the UniProt website:
        https://www.uniprot.org/help/return_fields

    Returns:
        - A list of UniProt IDs if no fields were provided.
        - A dict of UniProt IDs and corresponding field values if
          exactly one field was provided.
        - A dict with field names as top level keys and dicts of the
          kind described in the previous point as values.
    """

    return UniprotQuery(*query, fields = fields, **kwargs).perform()


def uniprot_preprocess(field, organism = 9606, reviewed = True):

    relabel = re.compile(r'[A-Z\s]+:\s')
    reisoform = re.compile(r'\[[-\w\s]+\]:?\s?')
    retermsep = re.compile(r'\s?[\.,]\s?')
    reref = re.compile(r'\{[-\w :\|,\.]*\}')

    result = collections.defaultdict(set)

    data = uniprot_data(
        fields = field,
        organism = organism,
        reviewed = reviewed,
    )

    for uniprot, raw in iteritems(data):

        raw = raw.split('Note=')[0]
        raw = relabel.sub('', raw)
        raw = reref.sub('', raw)
        raw = reisoform.sub('', raw)
        raw = retermsep.split(raw)

        for item in raw:

            if item.startswith('Note'):

                continue

            item = item.split('{')[0]
            elements = tuple(
                it0
                for it0 in
                (
                    common.upper0(it.strip(' .;,'))
                    for it in item.split(';')
                )
                if it0
            )

            if elements:

                result[uniprot].add(elements)

    return result


def uniprot_locations(organism = 9606, reviewed = True):


    UniprotLocation = collections.namedtuple(
        'UniprotLocation',
        [
            'location',
            'features',
        ],
    )


    result = collections.defaultdict(set)

    data = uniprot_preprocess(
        field = 'subcellular_location',
        organism = organism,
        reviewed = reviewed,
    )

    for uniprot, locations in iteritems(data):

        for location in locations:

            result[uniprot].add(
                UniprotLocation(
                    location = location[0],
                    features = location[1:] or None,
                )
            )

    return dict(result)


def uniprot_keywords(organism = 9606, reviewed = True):

    UniprotKeyword = collections.namedtuple(
        'UniprotKeyword',
        [
            'keyword',
        ],
    )


    result = collections.defaultdict(set)

    data = uniprot_data(
        fields = 'keywords',
        organism = organism,
        reviewed = reviewed,
    )

    for uniprot, keywords in iteritems(data):

        for keyword in keywords.split(';'):

            result[uniprot].add(
                UniprotKeyword(
                    keyword = keyword.strip(),
                )
            )

    return dict(result)


def uniprot_families(organism = 9606, reviewed = True):

    refamily = re.compile(r'(.+) (?:super)?family(?:, (.*) subfamily)?')


    UniprotFamily = collections.namedtuple(
        'UniprotFamily',
        [
            'family',
            'subfamily',
        ],
    )


    result = collections.defaultdict(set)

    data = uniprot_data(
        fields = 'family',
        organism = organism,
        reviewed = reviewed,
    )

    for uniprot, family in iteritems(data):

        if not family:

            continue

        family, subfamily = refamily.search(family).groups()

        result[uniprot].add(
            UniprotFamily(
                family = family,
                subfamily = subfamily,
            )
        )

    return dict(result)


def uniprot_topology(organism = 9606, reviewed = True):

    retopo = re.compile(r'TOPO_DOM (\d+)\.\.(\d+);\s+/note="(\w+)"')
    retm = re.compile(r'(TRANSMEM|INTRAMEM) (\d+)\.\.(\d+);')


    UniprotTopology = collections.namedtuple(
        'UniprotTopology',
        [
            'topology',
            'start',
            'end',
        ],
    )


    result = collections.defaultdict(set)

    transmem = uniprot_data(
        fields = 'transmembrane',
        organism = organism,
        reviewed = reviewed,
    )

    intramem = uniprot_data(
        fields = 'intramembrane',
        organism = organism,
        reviewed = reviewed,
    )

    signal = uniprot_data(
        fields = 'signal_peptide',
        organism = organism,
        reviewed = reviewed,
    )

    data = uniprot_data(
        fields = 'topological_domain',
        organism = organism,
        reviewed = reviewed,
    )

    for uniprot, topo in iteritems(data):

        for topo_dom in retopo.findall(topo):

            start, end, topology = topo_dom
            start = int(start)
            end = int(end)

            result[uniprot].add(
                UniprotTopology(
                    topology = topology,
                    start = start,
                    end = end,
                )
            )

    for uniprot, tm in itertools.chain(
        iteritems(transmem),
        iteritems(intramem),
        iteritems(signal),
    ):

        for mem, start, end in retm.findall(tm):

            topology = (
                '%s%s' % (
                    mem.capitalize(),
                    'brane' if mem.endswith('MEM') else ''
                )
            )
            start = int(start)
            end = int(end)

            result[uniprot].add(
                UniprotTopology(
                    topology = topology,
                    start = start,
                    end = end,
                )
            )

    return dict(result)


def uniprot_tissues(organism = 9606, reviewed = True):

    reref = re.compile(r'\s?\{.*\}\s?')
    resep = re.compile(
        r',?(?:'
            r' in almost all |'
            r' but also in |'
            r' but also at |'
            r' within the |'
            r', in |'
            r' in |'
            r' but |'
            r', and |'
            r' and |'
            r' such as |'
            r' \(both |'
            r' as well as |'
            r' as |'
            r' or |'
            r' at the |'
            r' at |'
            r' including |'
            r' during |'
            r' especially |'
            r' to |'
            r' into |'
            r' = |'
            r' > |'
            r'; |'
            r', '
        r')(?=[^\d])'
    )
    relabel = re.compile(r'^TISSUE SPECIFICITY: ')
    repubmed = re.compile(r'\(?PubMed:?\d+\)?')
    respeci = re.compile(r'(\w+)[-\s]specific')
    rethe = re.compile(
        r'\s?(?:'
           r'[Tt]he |'
           r'[Ii]n |'
           r'[Ss]ome|'
           r'[Ii]n the|'
           r'[Ww]ithin the|'
           r'[Ww]ithin|'
           r'[Ii]nto|'
           r'[Ww]ith only|'
           r'[Ww]ith the|'
           r'[Ww]ith an|'
           r'[Ww]ith |'
           r'[Ii]s |'
           r'[Mm]any  |'
           r'[Aa] variety of '
           r'[Aa] |'
           r'[Ii]t |'
           r'[Tt]o |'
           r'[Oo]n |'
           r'[Oo]f |'
           r'[Tt]hose |'
           r'[Ff]rom |'
           r'[Aa]lso|'
           r'[Bb]y |'
           r'[Pp]articularly|'
           r'[Pp]articular|'
           r'[Pp]atients|'
           r'[Aa]n |'
           r'\'|'
           r':|'
           r'/'
        r')?(.*)'
    )
    reand = re.compile(r'(?: and| of| from| or| than)$')
    replevel = re.compile(r'\(at \w+ levels?\)')
    reiso = re.compile(r'[Ii]soform \w+')
    reindef = re.compile(
        r'\w'
        r'(?:'
           r'ifferent parts of |'
           r'ariety of tissues |'
           r' variety of tissues |'
           r' number of |'
           r'everal regions of '
        r')'
    )

    level_kw = (
        ('low', 'low'),
        ('weak', 'low'),
        ('lesser extent', 'low'),
        ('minimal level', 'low'),
        ('decrease', 'low'),
        ('moderate', 'low'),
        ('barely', 'low'),
        ('minor level', 'low'),
        ('reduced', 'low'),
        ('lesser', 'low'),
        ('down-regulated', 'low'),
        ('high', 'high'),
        ('elevated', 'high'),
        ('strong', 'high'),
        ('prominent', 'high'),
        ('greatest level', 'high'),
        ('concentrated', 'high'),
        ('predominant', 'high'),
        ('increase', 'high'),
        ('enrich', 'high'),
        ('abundant', 'high'),
        ('primarily', 'high'),
        ('induced', 'high'),
        ('up-regulated', 'high'),
        ('up regulated', 'high'),
        ('expression is restricted', 'high'),
        ('amplified', 'high'),
        ('basal l', 'basal'),
        ('not detected', 'none'),
        ('absent', 'none'),
        ('expressed', 'undefined'),
        ('detect', 'undefined'),
        ('found', 'undefined'),
        ('present', 'undefined'),
        ('expression', 'undefined'),
        ('localized', 'undefined'),
        ('produced', 'undefined'),
        ('confined', 'undefined'),
        ('transcribed', 'undefined'),
        ('xpressed', 'undefined'),
        ('synthesized', 'undefined'),
        ('secreted', 'undefined'),
        ('seen', 'undefined'),
        ('prevalent', 'undefined'),
        ('released', 'undefined'),
        ('appears', 'undefined'),
        ('varying levels', 'undefined'),
        ('various levels', 'undefined'),
        ('identified', 'undefined'),
        ('observed', 'undefined'),
        ('occurs', 'undefined'),
    )

    wide_kw = (
        ('widely', 'wide'),
        ('wide tissue distribution', 'wide'),
        ('wide range of tissues', 'wide'),
        ('wide range of adult tissues', 'wide'),
        ('wide range of cells', 'wide'),
        ('wide variety of normal adult tissues', 'wide'),
        ('widespread', 'wide'),
        ('ubiquitous', 'ubiquitous'),
        ('variety of tissues', 'wide'),
        ('many tissues', 'wide'),
        ('many organs', 'wide'),
        ('various organs', 'wide'),
        ('various tissues', 'wide'),
    )

    tissue_exclude = {
        'Adult',
        'All',
        'Apparently not',
        'Areas',
        'Are likely',
        'Both',
        'By contrast',
        'Normal cells',
        'Not only',
        'A',
        '[]: Localized',
        'Early',
        'Change from a quiescent',
        'Central',
        'Beta',
        'This layer',
        'With little',
        'Preferential occurrence',
        'Stage III',
        'Take up',
        'Hardly',
        'Only seen',
        'Prevalent',
        'Inner segment',
        'Memory',
        'Many fetal',
        'Tissues',
        '0 kb',
        '9 kb',
        'A 2',
        'A 3',
        'A 5',
        'A 6',
        '1-7',
        '1b-1',
        '2 is widely',
        '8 and 4',
        'Often amplified',
        'Other',
        'Others',
        'Those',
        'Tissues examined',
        'Tissues with',
        'Tissues (e)',
        'Probably shed',
        'Reports that',
        'Primitive',
        'Prolactin',
        'Overlap',
        'A smaller 0',
        'A smaller form',
        'A smaltissues',
        'Different levels',
        'Different amounts',
        'Disappears',
        'Digestion',
        'Very similar',
        'Vivo',
        'Contrary',
        'Contrast',
        'Not',
        'Not all',
        'Has it',
        'Has little',
        'All stages',
        'Soon',
        'Specific',
        'Stage',
        'Stage I',
        'Stage II',
        'Stages II',
        'Ends',
        'A minor degree',
        'A much smaller extent',
        'Lost',
        'Varies',
        'Various',
        'Mostly restricted',
        'Mostly',
        'Most probably',
        'Much more stable',
        'Naive',
        'Neither',
        'Nor',
        'None',
    }

    exclude_startswith = (
        'Were',
        'Where',
        'Which',
        'While',
        'When',
        'There',
        'Their',
        'Then',
        'These',
        'Level',
        'This',
        'Almost',
        'If',
        'Control',
        'Be ',
        'Although',
        'Than',
        'Addition',
    )

    exclude_in = (
        'kb transcript',
        'compared',
        'soform',
        'concentration of'
    )


    UniprotTissue = collections.namedtuple(
        'UniprotTissue',
        [
            'tissue',
            'level',
        ],
    )


    data = uniprot_data(
        fields = 'tissue_specificity',
        organism = organism,
        reviewed = reviewed,
    )

    result = collections.defaultdict(set)

    for uniprot, raw in iteritems(data):

        raw = relabel.sub('', raw)
        raw = reref.sub('', raw)
        raw = replevel.sub('', raw)
        raw = reiso.sub('', raw)
        raw = repubmed.sub('', raw)
        raw = reindef.sub('', raw)
        raw = raw.replace('adult and fetal', '')

        raw = raw.split('.')

        for phrase in raw:

            tokens = tuple(resep.split(phrase))
            level = None

            for token in tokens:

                level_token = False
                wide_token = False
                tissue = None

                token_lower = token.lower()

                for kw, lev in level_kw:

                    if kw in token_lower:

                        level = lev
                        level_token = True
                        break

                if level_token:

                    for kw, wide in wide_kw:

                        if kw in token_lower:

                            tissue = wide
                            wide_token = True
                            break

                if not level_token or wide_token:

                    if not wide_token:

                        specific = respeci.search(token)

                        tissue = (
                            specific.groups()[0].lower()
                                if specific else
                            token
                        )

                        if specific and not level:

                            level = 'high'

                    if tissue.strip():

                        if any(e in tissue for e in exclude_in):

                            continue

                        tissue = rethe.match(tissue).groups()[0]
                        tissue = rethe.match(tissue).groups()[0]
                        tissue = rethe.match(tissue).groups()[0]

                        if tissue.endswith('+'):

                            tissue = '%s cells' % tissue

                        tissue = tissue.strip(')(.,;- ')

                        if '(' in tissue and ')' not in tissue:

                            tissue = '%s)' % tissue

                        tissue = reand.sub('', tissue)
                        tissue = common.upper0(tissue)
                        tissue = tissue.replace('  ', ' ')

                        if any(
                            tissue.startswith(e)
                            for e in exclude_startswith
                        ):

                            continue

                        if tissue in tissue_exclude or len(tissue) < 3:

                            continue

                        result[uniprot].add(
                            UniprotTissue(
                                tissue = tissue,
                                level = level or 'undefined',
                            )
                        )

    return dict(result)


def uniprot_taxonomy(
        ncbi_tax_ids: bool = False,
    ) -> dict[str, set[str]] | dict[str, int]:
    """
    From UniProt IDs to organisms

    Args:
        ncbi_tax_ids:
            Translate the names to NCBI Taxonomy numeric identifiers.

    Returns:
        A dictionary with SwissProt IDs as keys and sets of various taxon
        names as values.
    """

    rename = re.compile(r'\(?(\w[\w\s\',/\.-]+\w)\)?')
    reac = re.compile(r'\s*\w+\s+\(([A-Z\d]+)\)\s*,')

    url = urls.urls['uniprot_basic']['speindex']
    c = curl.Curl(url, large = True, silent = False)

    result = collections.defaultdict(set)

    for line in c.result:

        if line[0] != ' ':

            names = set(rename.findall(line))

        else:

            for ac in reac.findall(line):

                result[ac].update(names)

    if ncbi_tax_ids:

        new_result = {}

        for ac, names in result.items():

            for name in names:

                nti = taxonomy.ensure_ncbi_tax_id(name)

                if nti:

                    new_result[ac] = nti
                    break

        result = new_result

    return dict(result)


Taxon = collections.namedtuple(
    'Taxon',
    [
        'ncbi_id',
        'latin',
        'english',
        'latin_synonym',
    ]
)
Taxon.__new__.__defaults__ = (None, None)


def uniprot_ncbi_taxids():

    url = urls.urls['uniprot_basic']['taxids']

    with settings.context(curl_timeout = 10000):

        c = curl.Curl(
            url,
            large = True,
            silent = False,
            compr = 'gz',
        )

    _ = next(c.result)

    result = {}

    for line in c.result:

        line = line.split('\t')

        if line[0].isdigit() and len(line) > 2:

            taxid = int(line[0])

            result[taxid] = Taxon(
                ncbi_id = taxid,
                latin = line[2],
                english = line[1] or None,
            )

    return result


def uniprot_ncbi_taxids_2():

    reline = re.compile(
        r'(?:([A-Z\d]+)\s+)?' # code
        r'(?:([A-Z]))?\s+' # kingdom
        r'(?:(\d+): )?' # NCBI Taxonomy ID
        r'([A-Z])=' # name type
        r'([ \w\(\),/\.\'-]+)[\n\r\s]*' # the name
    )

    url = urls.urls['uniprot_basic']['speclist']
    c = curl.Curl(url, large = True, silent = False)

    result = {}
    entry = {}

    for line in c.result:

        m = reline.match(line)

        if m:

            _code, _kingdom, _taxid, _name_type, _name = m.groups()

            if _taxid:

                if entry and 'ncbi_id' in entry:

                    result[entry['ncbi_id']] = Taxon(**entry)

                entry = {}
                entry['ncbi_id'] = int(_taxid)

            if _name_type == 'N':

                entry['latin'] = _name

            elif _name_type == 'C':

                entry['english'] = _name

            elif _name_type == 'S':

                entry['latin_synonym'] = _name

    if entry and 'ncbi_id' in entry:

        result[entry['ncbi_id']] = Taxon(**entry)

    return result
