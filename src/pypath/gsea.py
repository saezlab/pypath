#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import re
import bs4
import sys
import os
try:
    import cPickle as pickle
except ImportError:
    import pickle

from collections import OrderedDict

# from this module:
import pypath.dataio as dataio
import pypath.data_formats as data_formats
import pypath.enrich as enrich
import pypath.mapping as mapping
import pypath.progress as progress
import pypath.common as common


class GSEA(object):
    def __init__(self, user=None, mapper=None):
        self.user = user if user is not None \
            else globals()['MSIGDB_USER'] if 'MSIGDB_USER' in globals() \
            else None
        if self.user is not None:
            self.login()
            self.mapper = mapper if mapper is not None else mapping.Mapper()
            self.info = {}
            self.groups = {}
            self.sets = {}
            self.collections = {}
            self.list_collections()
            self.ids = {'entrez': 'entrez', 'symbol': 'genesymbol'}
            self.target_id = 'uniprot'
        else:
            sys.stdout.write('\t:: Please provide an MSigDB username by \n'
                             '``pypath.gsea.GSEA(user = \'\', ...)``, or by \n'
                             'setting ``MSIGDB_USER`` global variable.\n\n')
            sys.stdout.flush()

    def login(self):
        url = data_formats.urls['msigdb']['login1']
        self.pre_session = dataio.curl(
            url, init_url=url, silent=False, cache=False, init_headers=True)
        url = data_formats.urls['msigdb']['login2']
        post = {'j_username': self.user, 'j_password': 'password'}
        self.session = dataio.curl(
            url,
            init_url=url,
            post=post,
            req_headers=self.pre_session,
            silent=False,
            cache=False,
            init_headers=True)

    def list_collections(self):
        renm = re.compile(r'(.+)\([^0-9]*([0-9]*)[^0-9]*\)')
        url = data_formats.urls['msigdb']['coll']
        html = dataio.curl(url, req_headers=self.session, silent=False)
        soup = bs4.BeautifulSoup(html, 'lxml')
        for col in soup.find('table', class_='lists1').find_all('tr'):
            lname, num = renm.findall(col.find('th').text.replace('\n', ''))[0]
            sname = col.find('a').attrs['name']
            urls = dict(
                [(d.attrs['href'].split('.')[-2],
                  data_formats.urls['msigdb']['url_stem'] % d.attrs['href'])
                 for d in col.find_all('a')[-3:]])
            self.collections[sname] = {
                'name': lname,
                'count': int(num),
                'urls': urls
            }

    def show_collections(self):
        s = '\n :: Available gene set collections:\n\n'\
            + '\tID\t\t\t#genes\tDescription\n\t%s\n\t' % ('-' * 75) \
            + '\n\t'.join('%s\t\t\t%u\t%s' % (sname, inf['count'], inf['name'])
                          for sname, inf in iteritems(self.collections)) \
            + '\n'
        sys.stdout.write(s)
        sys.stdout.flush()

    def load_collection(self,
                        collname,
                        id_type='entrez',
                        map_ids=True,
                        cachedir='cache'):
        if os.path.exists(os.path.join(cachedir, 'gsea-%s.pickle' % collname)):
            self.load([collname])
            return None
        url = self.collections[collname]['urls'][id_type]
        data = dataio.curl(
            url,
            req_headers=self.session,
            silent=False,
            cache=False,
            write_cache=True)
        data = data.split('\n')
        names = []
        prg = progress.Progress(len(data), 'Loading gene sets', 1)
        for line in (l.split('\t') for l in data if len(l) > 0):
            prg.step()
            setname = line[0].strip()
            self.write_set(line[2:], setname, id_type, map_ids)
            self.get_desc(setname)
            names.append(setname)
        prg.terminate()
        self.groups[collname] = set(names)
        self.save([collname], cachedir=cachedir)

    def get_desc(self, setname):
        url = data_formats.urls['msigdb']['one_set'] % setname
        txt = dataio.curl(url, req_headers=self.session, silent=True)
        self.info[setname] = txt.split('\n')[1][2:]

    def load_set(self, setname, map_ids=True):
        url = data_formats.urls['msigdb']['one_set'] % setname
        data = dataio.curl(url, req_headers=self.session, silent=True)
        data = data.split('\n')
        self.info[setname] = data[1][2:]
        self.write_set((j for j in (i.strip() for i in data[2:])
                        if len(j) > 0), setname, 'symbol', map_ids)

    def write_set(self, id_list, setname, id_type, map_ids=True):
        self.sets[setname] = set(common.uniqList(common.flatList(
            self.mapper.map_name(n, self.ids[id_type], self.target_id)
            for n in id_list))) if map_ids \
            else set(id_list)

    def save(self, collections, cachedir='cache'):
        if collections is None:
            collections = self.groups.keys()
        if type(collections) is not list:
            collections = [collections]
        collections = set(collections) & set(self.groups.keys())
        for coll in collections:
            pfile = os.path.join(cachedir, 'gsea-%s.pickle' % coll)
            sets = dict([(k, v) for k, v in iteritems(self.sets)
                         if k in self.groups[coll]])
            info = dict([(k, v) for k, v in iteritems(self.info)
                         if k in self.groups[coll]])
            group = self.groups[coll]
            pickle.dump((sets, info, group), open(pfile, 'wb'))

    def load(self, collections, cachedir='cache'):
        if type(collections) is not list:
            collections = [collections]
        for coll in collections:
            pfile = os.path.join(cachedir, 'gsea-%s.pickle' % coll)
            if os.path.exists(pfile):
                sets, info, group = pickle.load(open(pfile, 'rb'))
                self.sets = dict(self.sets.items() + sets.items())
                self.info = dict(self.info.items() + info.items())
                self.groups[coll] = group


class GSEABinaryEnrichmentSet(enrich.EnrichmentSet):
    def __init__(self,
                 basic_set,
                 gsea=None,
                 geneset_ids=None,
                 alpha=0.05,
                 correction_method='hommel',
                 user=None,
                 mapper=None):
        if type(gsea) is not GSEA and geneset_ids is None:
            console('Please give either a `pypath.gsea.GSEA` object'
                    'or a list of geneset names.')
        if geneset_ids is None:
            geneset_ids = gsea.sets.keys()
        self.user = user
        self.mapper = mapper
        if type(gsea) is not GSEA:
            gsea = GSEA(user=self.user, mapper=self.mapper)
            for geneset_id in geneset_ids:
                gsea.load_set(geneset_id)
        self.geneset_ids = geneset_ids
        self.gsea = gsea
        self.alpha = alpha
        self.correction_method = correction_method
        self.basic_set = set(basic_set)
        self.counts_pop = self.count(self.basic_set)
        self.pop_size = len(self.basic_set)
        self.set_size = None
        self.counts_set = None
        self.top_geneset_ids = self.top_ids
        self.top_genesets = self.top_names

    def count(self, this_set):
        return dict((i, len(this_set & s))
                    for i, s in iteritems(self.gsea.sets))

    def new_set(self, set_names):
        if type(set_names) is not set:
            set_names = set(set_names)
        self.set_size = len(set_names)
        self.counts_set = self.count(set_names)
        self.calculate()

    def calculate(self):
        data = dict([(gset_id, (cnt, self.counts_pop[gset_id], self.set_size,
                                self.gsea.info[gset_id]))
                     for gset_id, cnt in iteritems(self.counts_set)])
        enrich.EnrichmentSet.__init__(
            self,
            data,
            self.pop_size,
            alpha=self.alpha,
            correction_method=self.correction_method)

    def toplist(self,
                length=None,
                alpha=None,
                significant=True,
                min_set_size=0,
                groups=None,
                filtr=lambda x: True,
                **kwargs):
        args = get_args(locals(), ['filtr', 'groups'])
        if groups is None:
            groups = self.gsea.groups.keys()  # all by default
        sets = set(
            common.flatList(s for g, s in iteritems(self.gsea.groups)
                            if g in groups))
        return super(GSEABinaryEnrichmentSet, self).toplist(
            filtr=lambda x: x[0] in sets and filtr(x), **args)

    def __str__(self):
        if self.counts_set is None:
            resp = '\n\t:: No calculations performed yet. Please define '\
                'a set of genes with `new_set()`.\n\n'
        else:
            resp = '\n :: Top significantly enriched genesets (max. 10):\n\n\t'\
                + '\n\t'.join([t[0].upper() + t[1:] for t in
                               self.top_genesets(length=10, significant=True)]) + '\n'
        return resp
