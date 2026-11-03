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
Interface to UniProt protein datasheets.
"""

from future.utils import iteritems

import os
import sys
import re
import shutil
import importlib as imp
import collections
import itertools

import pypath.inputs.uniprot as uniprot_input
import pypath.inputs.genecards as genecards_input
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.share.settings as settings
import pypath.core.entity as entity


class UniprotProtein(object):

    _relength = re.compile(r'([0-9]+) AA')
    _rename = re.compile(r'Name=([\w\(\)-]+)\W')
    _rerecname = re.compile(r'(?:Rec|Sub)Name: Full=([^;\{]+)(?: \{.*\})?;')
    _recc = re.compile(r'-!- ([A-Z ]+):\s?(.*)')
    _remw = re.compile(r'([0-9]+) MW')
    _redb = re.compile(r'([^;]+);\s?(.*)\s?\.\s?(?:\[(.*)\])?')
    _redbsep = re.compile(r'\s?;\s?')
    _retaxid = re.compile(r'=(\d+)[^\d]')
    _rexref = re.compile(r'[\.,]?\s?\{[^\}]+\}')
    _reec = re.compile(r'EC=(\d+(?:\.[-\d]+)+)')

    def __init__(self, uniprot_id):

        self.uniprot_id = uniprot_id.strip()
        self.load()


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def load(self):

        self.raw = uniprot_input.protein_datasheet(self.uniprot_id)


    @property
    def is_reviewed(self):

        return 'Reviewed' in self.raw[0][1]


    @property
    def id(self):

        return self.raw[0][1].split()[0]


    @property
    def ac(self):

        return next(self.itertag('AC')).split(';')[0]


    @property
    def length(self):
        """
        Returns the length (number of residues) of the canonical sequence.
        """

        return int(self._relength.search(self.raw[0][1]).groups()[0])


    @property
    def organism(self):

        return int(self._retaxid.search(next(self.itertag('OX'))).groups()[0])


    @property
    def full_name(self):

        return self._rerecname.search(next(self.itertag('DE'))).groups()[0]


    @property
    def ec(self):

        return set(self._reec.findall(''.join(self.itertag('DE'))))


    @property
    def info(self):

        if not hasattr(self, '_info'):

            self.update_info()

        return self._info


    def update_info(self):

        result = collections.defaultdict(list)
        title = None

        for cc in self.itertag('CC'):

            if cc.startswith('---'):

                break

            m = self._recc.match(cc)

            if m:

                title, cc = m.groups()

            line = cc.strip()

            if line:

                result[title].append(line)

        self._info = dict(
            (
                title,
                ' '.join(line)
            )
            for title, line in iteritems(result)
        )


    @property
    def function_genecards(self):

        summaries = genecards_input.genecards_summaries(self.genesymbol)

        return ' '.join(
            '%s: %s' % (resource, summary)
            for resource, summary in iteritems(summaries)
        )


    @property
    def function_with_xrefs(self):

        return self.info_section('FUNCTION')


    @property
    def function(self):

        return self.remove_xrefs(self.function_with_xrefs)


    @property
    def function_with_genecards(self):

        return '%s %s' % (
            self.function,
            self.function_genecards,
        )


    @property
    def function_or_genecards(self):

        return self.function or self.function_genecards


    @property
    def subcellular_location(self):

        return self.remove_xrefs(self.subcellular_location_with_xrefs)


    @property
    def tissue_specificity(self):

        return self.remove_xrefs(self.tissue_specificity_with_xrefs)


    @property
    def subunit(self):

        return self.remove_xrefs(self.subunit_with_xrefs)


    @property
    def interaction(self):

        return self.remove_xrefs(self.interaction_with_xrefs)


    @property
    def sequence_caution(self):

        return self.remove_xrefs(self.sequence_caution_with_xrefs)


    @property
    def catalytic_activity(self):

        return self.remove_xrefs(self.catalytic_activity_with_xrefs)


    @property
    def activity_regulation(self):

        return self.remove_xrefs(self.activity_regulation_with_xrefs)


    @property
    def alternative_products(self):

        return self.remove_xrefs(self.alternative_products_with_xrefs)


    @property
    def ptm(self):

        return self.remove_xrefs(self.ptm_with_xrefs)


    @property
    def disease(self):

        return self.remove_xrefs(self.disease_with_xrefs)


    @property
    def similarity(self):

        return self.remove_xrefs(self.similarity_with_xrefs)


    @property
    def web_resource(self):

        return self.remove_xrefs(self.web_resource_with_xrefs)


    @property
    def subcellular_location_with_xrefs(self):

        return self.info_section('SUBCELLULAR LOCATION')


    @property
    def tissue_specificity_with_xrefs(self):

        return self.info_section('TISSUE SPECIFICITY')


    @property
    def subunit_with_xrefs(self):

        return self.info_section('SUBUNIT')


    @property
    def interaction_with_xrefs(self):

        return self.info_section('INTERACTION')


    @property
    def sequence_caution_with_xrefs(self):

        return self.info_section('SEQUENCE CAUTION')


    @property
    def catalytic_activity_with_xrefs(self):

        return self.info_section('CATALYTIC ACTIVITY')


    @property
    def activity_regulation_with_xrefs(self):

        return self.info_section('ACTIVITY REGULATION')


    @property
    def alternative_products_with_xrefs(self):

        return self.info_section('ALTERNATIVE PRODUCTS')


    @property
    def ptm_with_xrefs(self):

        return self.info_section('PTM')


    @property
    def disease_with_xrefs(self):

        return self.info_section('DISEASE')


    @property
    def similarity_with_xrefs(self):

        return self.info_section('SIMILARITY')


    @property
    def web_resource_with_xrefs(self):

        return self.info_section('WEB RESOURCE')


    @property
    def lengths(self):
        """
        Returns the length of all isoforms as a list.
        """

        return [
            int(self._relength.search(sq).groups()[0])
            for sq in self.itertag('SQ')
        ]


    @property
    def weight(self):
        """
        Returns the molecular weight of the canonical isoform in Daltons.
        """

        try:

            return int(
                self._remw.search(next(self.itertag('SQ'))).groups()[0]
            )

        except StopIteration:

            return None


    @property
    def weights(self):
        """
        Returns the molecular weights of all isoforms as a list.
        """

        return [
            int(self._remw.search(sq).groups()[0])
            for sq in self.itertag('SQ')
        ]


    @property
    def databases(self):
        """
        Returns the database identifiers (cross-references) as a dict of
        database names and identifiers.
        """

        if not hasattr(self, '_databases'):

            self.update_databases()

        return self._databases


    def update_databases(self):

        result = collections.defaultdict(set)

        for db in self.itertag('DR'):

            m = self._redb.match(db)

            if m:

                dbname, ids, subtype = m.groups()
                ids = self._redbsep.split(ids)
                ids = tuple(_id for _id in ids if _id != '-')

                if subtype:

                    ids += (subtype,)

                ids = ids[0] if len(ids) == 1 else ids
                result[dbname].add(ids)

        self._databases = dict(result)


    def info_section(self, title):
        """
        Retrieves a section from the description. If the section is not
        availeble, returns ``None``.
        """

        info = self.info

        if title in info:

            return info[title]

    @property
    def genesymbol(self):

        try:

            m = self._rename.search(next(self.itertag('GN')))

            return m.groups()[0] if m else self.ac

        except StopIteration:

            return self.ac


    @property
    def keywords_with_xrefs(self):
        """
        Returns the keywords as a list with keeping the cross-references.
        """

        return [
            kw for kw in
            itertools.chain(
                *(
                    self._redbsep.split(kw.strip('.'))
                    for kw in self.itertag('KW')
                )
            )
            if kw
        ]


    @property
    def keywords(self):
        """
        Returns the keywords as a list.
        """

        return (
            self.remove_xrefs(
                '\t'.join(self.keywords_with_xrefs)
            ).split('\t')
        )


    @classmethod
    def remove_xrefs(cls, value):

        return cls._rexref.sub('', value) if value else value


    @property
    def sequence(self):
        """
        Returns the canonical sequence (the first one) as a string of
        standard capital letter residue symbols.
        """

        result = []
        collect = False

        for tag, line in self:

            if not collect and tag == 'SQ':

                collect = True

            elif collect:

                if tag == '  ':

                    result.append(line)

                else:

                    break

        return ''.join(x.replace(' ', '') for x in result)


    def __iter__(self):

        return self.raw.__iter__()


    def itertag(self, tag):

        for _tag, line in self:

            if _tag == tag:

                yield line


    def has_tag(self, tag):

        return any(line[0] == tag for line in self)


    def __repr__(self):

        return '<UniProt datasheet %s (%s)>' % (self.ac, self.genesymbol)


def _update_methods():

    for method_name in UniprotProtein.__dict__.keys():

        if method_name.startswith('_'):

            continue

        def create_method(method_name):

            def method(uniprot_id, *args, **kwargs):

                bound_m = getattr(UniprotProtein(uniprot_id), method_name)

                if isinstance(getattr(UniprotProtein, method_name), property):

                    return bound_m

                else:

                    return bound_m(*args, **kwargs)

            return method

        _method = create_method(method_name)

        common.add_method(
            cls = sys.modules[__name__],
            method_name = method_name,
            method = _method,
            doc = getattr(UniprotProtein, method_name).__doc__,
        )


_update_methods()


def query(*uniprot_ids):
    """
    Queries the datasheet of one or more UniProt IDs.
    Returns a single ``UniprotProtein`` object or a list of those objects.
    """

    if (
        len(uniprot_ids) > 0 and
        isinstance(uniprot_ids[0], _const.LIST_LIKE)
    ):

        uniprot_ids = uniprot_ids[0]

    uniprot_ids = common.to_list(uniprot_ids)
    uniprot_ids = entity.Entity.only_proteins(uniprot_ids)

    single_id = len(uniprot_ids) == 1

    result = [
        UniprotProtein(uniprot_id)
        for uniprot_id in uniprot_ids
    ]
    result = [u for u in result if u.raw]

    return common.first(result) if single_id else result


def collect(uniprot_ids, *features):
    """
    Collects data about one or more UniProt IDs.

    :param str,list uniprot_ids:
        One or more UniProt IDs.
    :param *str,list features:
        Features to query: these must be method (property) names of the
        ``UniprotProtein`` class. E.g. ``['ac', 'genesymbol', 'function']``.

    :return:
        A ``collections.OrderedDict`` object with feature names as keys and
        list of values for each UniProt ID as values.
    """

    uniprot_ids = entity.Entity.only_proteins(uniprot_ids)

    resources = [
        UniprotProtein(uniprot_id)
        for uniprot_id in uniprot_ids
    ]
    # this is mainly for removal of obsolate records
    # where the response from the server is empty
    # most of the times it removes nothing
    resources = [u for u in resources if u.raw]

    features = features or default_features

    if 'ac' not in features:

        features = ['ac'] + list(features)

    table = collections.OrderedDict(
        (
            feature_name,
            [
                getattr(resource, feature_name)
                for resource in resources
            ]
        )
        for feature_name in features
    )

    return table


def features_table(
        uniprot_ids,
        *features,
        width = 40,
        maxlen = None,
        tablefmt = 'fancy_grid',
        **kwargs
    ):
    """
    Returns a table with the requested features of a list of UniProt IDs.
    The underlying table formatting module is ``tabulate``, a versatile
    module to export various ascii tables as well as HTML or LaTeX --
    check the docs for formatting options:
    https://github.com/astanin/python-tabulate

    Args
        kwargs:
            Passed to ``tabulate.tabulate``.

    Returns
        The table as a string.
    """

    maxlen = maxlen or settings.get('uniprot_info_maxlen')

    features = features or default_features

    tbl = collect(uniprot_ids, *features)

    return common.table_format(
        tbl,
        width = width,
        maxlen = maxlen,
        tablefmt = tablefmt,
        **kwargs
    )


def print_features(
        uniprot_ids,
        *features,
        fileobj = None,
        width = None,
        maxlen = None,
        tablefmt = 'fancy_grid',
        **kwargs
    ):
    """
    Prints a table with the requested features of a list of UniProt IDs.
    The underlying table formatting module is ``tabulate``, a versatile
    module to export various ascii tables as well as HTML or LaTeX --
    check the docs for formatting options:
    https://github.com/astanin/python-tabulate

    Args
        kwargs:
            Passed to ``tabulate.tabulate``.
    """

    maxlen = maxlen or settings.get('uniprot_info_maxlen')
    features = features or default_features
    term_width = (shutil.get_terminal_size().columns - 60) * 2 + 40
    width = width or int(term_width / len(features)) if term_width else 40
    fileobj = fileobj or sys.stdout

    fileobj.write(
        features_table(
            uniprot_ids,
            *features,
            width = width,
            maxlen = maxlen,
            tablefmt = tablefmt,
            **kwargs
        )
    )
    fileobj.write(os.linesep)
    fileobj.flush()


def info(
        *uniprot_ids,
        features = None,
        fileobj = None,
        header = None,
        **kwargs
    ):
    """
    Prints a table with the most important (or the requested) features of a
    list of UniProt IDs.
    """

    if (
        len(uniprot_ids) == 1 and
        isinstance(uniprot_ids, _const.LIST_LIKE)
    ):

        uniprot_ids = uniprot_ids[0]

    features = features or default_features

    fileobj = fileobj or sys.stdout

    header = (
        header or
        '=====> [%u proteins] <=====\n' % len(
            list(
                entity.Entity.filter_entity_type(
                    common.to_list(uniprot_ids),
                    entity_type = 'protein',
                )
            )
        )
    )

    fileobj.write(header)

    print_features(
        common.to_list(uniprot_ids),
        *features,
        fileobj = fileobj,
        **kwargs
    )


def browse(groups, start = 0, fileobj = None, **kwargs):
    """
    Browses through a series of protein groups, printing an information table
    for each group. ``kwargs`` passed to ``info`` and then to print_features``.
    Parameters for ``common.table_format`` can be provided.
    """

    labels = sorted(groups.keys())
    n_groups = len(labels)
    stop = False
    maxlen_default = kwargs['maxlen'] if 'maxlen' in kwargs else 500
    fileobj = fileobj or sys.stdout

    for n, label in enumerate(labels):

        if start > n + 1:

            continue

        if stop:

            break

        kwargs['maxlen'] = maxlen_default

        while True:

            uniprots = groups[label]
            uniprots = (
                uniprots.members
                    if hasattr(uniprots, 'members') else
                uniprots
            )

            header = (
                '[%u/%u] =====> %s <===== [%u proteins]\n' % (
                    n + 1,
                    n_groups,
                    label,
                    len(uniprots)
                )
            )

            info(uniprots, fileobj = fileobj, header = header, **kwargs)

            inp = input()

            if inp == 'q':

                stop = True
                break

            elif inp.isdigit():

                kwargs['maxlen'] = int(inp)

            else:

                fileobj.write(os.linesep * 2)
                break

    sys.stdout.write(os.linesep)
    sys.stdout.flush()


default_features = (
    'ac',
    'genesymbol',
    'length',
    'weight',
    'full_name',
    'function_or_genecards',
    'keywords',
    'subcellular_location',
)
