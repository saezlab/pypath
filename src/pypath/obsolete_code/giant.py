#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Connection to Giant database. This module was never used.
#
#  Copyright
#  2014-2019
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

# other modules:
import sys
import os
import json
import bs4
import cPickle as pickle
import hashlib
import random

# from this module:
import data_formats
import dataio
import mapping
import progress
from common import *


class Giant(object):
    '''
    This class is a client for GIANT (Genome wide Integrated Analysis of Networks in Tissues)
    webservice. The main input is the `graph` parameter of the `Giant.network()` function, 
    which is an iGraph object having UniProt IDs in its `name` node attribute. The `tissues`
    list is a list of tissue names (strings) or IDs (integers) as those are used in GIANT.
    If `tissues` is not given, data for all tissues will be downloaded. 

    Notice: downloading data for one tissue may take hours, depending on the size of the net-
    work, and downloading data for all the 144 tissues may take one week or more. But all 
    the values are cached and you can save them to a pickle file by Giant.save(), and you
    will be able to continue downloading from the point where it happend to be interrupted 
    or whatever reason.
    '''

    def __init__(self, mapper=None, debug=False):
        self.mapper = None
        self.init_mapper(mapper)
        self.tissnm = self.tissue_ids()
        self.tissid = dict((v, k) for k, v in self.tissnm.iteritems())
        self.cached = {}
        self.debug = debug
        if self.debug:
            self.trace = {
                'entrez_list': None,
                'input_edges': None,
                'url': None,
                'raw': None,
                'result': None,
                'edges_set': None,
                'edges_from_cache': None
            }

    def load(self):
        if os.path.exists(self.cache):
            self.cached = pickle.load(open(self.cache, 'rb'))

    def save(self):
        pickle.dump(self.cached, open(self.cache, 'wb'))

    def init_mapper(self, mapper):
        self.mapper = mapper if type(mapper) is mapping.Mapper \
            else self.mapper if type(self.mapper) is mapping.Mapper \
            else mapping.Mapper()

    def tissue_ids(self):
        result = {}
        url = data_formats.urls['giant']['init_url']
        html = dataio.curl(url, silent=False)
        soup = bs4.BeautifulSoup(html)
        for tis in soup.find('select', {'id': 'tissue'}).findAll('option'):
            result[int(tis.attrs['value'])] = tis.text
        return result

    def query(self, tissue, entrez):
        url = data_formats.urls['giant']['url'] % (
            tissue, 'entrez=%s' % '&entrez='.join([str(e) for e in entrez]))
        data = dataio.curl(
            url,
            silent=True,
            cache=False,
            debug=self.debug,
            override_post=True)
        if self.debug:
            self.trace['url'] = url
            self.trace['raw'] = data
            self.trace['result'] = json.loads(data)
        return data if data is None else json.loads(data)

    def genes_dict(self, data):
        # apparently edge source and target refer to the indexes of
        # genes in the genes list, so we need to translate those to
        # Entrez Gene IDs
        return dict((i, g['id']) for i, g in enumerate(data['genes']))

    def network(self,
                graph,
                tissues=None,
                cache_file='giant_cache.pickle',
                mapper=None,
                size=100):
        self.cache = cache_file
        self.load()
        self.init_mapper(mapper)
        if 'giant' not in graph.es.attributes():
            graph.es['giant'] = [{} for _ in graph.es]
        tissues = self.tissid.keys() if tissues is None else tissues
        tissues = [
            int(t) if type(t) in charTypes and t.isdigit() and
            int(t) in self.tissnm else t for t in tissues
        ]
        tissues = [self.tissid[t] if t in self.tissid else t for t in tissues]
        self.ventrez = [
            None if len(e) == 0 else e[0]
            for e in [
                mapper.map_name(v['name'], 'uniprot', 'entrez')
                for v in graph.vs
            ]
        ]
        self.eentrez = dict(((self.ventrez[e.source], self.ventrez[e.target]),
                             e.index) for e in graph.es)
        if len(set(tissues) - set(self.tissnm.keys())) > 0:
            sys.stdout.write(
                '\t:: The following tissues could not be identified in GIANT tissues list:\n\t\t%s\n'
                % str(list(set(tissues) - set(self.tissnm.keys()))))
            sys.stdout.flush()
        tissues = list(set(tissues) & set(self.tissnm.keys()))
        for tissue in tissues:
            for e in graph.es:
                if self.tissnm[tissue] not in e['giant']:
                    e['giant'][self.tissnm[tissue]] = None
            self.iterate(graph, tissue, size=size)
        self.save()

    def from_cache(self, one, two, tissue):
        if tissue in self.cached:
            this_pair = [one, two]
            this_pair = tuple(sorted(this_pair))
            if this_pair in self.cached[tissue]:
                return self.cached[tissue][this_pair]
        return None

    def to_cache(self, data, tissue):
        if tissue not in self.cached:
            self.cached[tissue] = {}
        self.cached[tissue][tuple(sorted([data[0], data[1]]))] = data[2]

    def iterate(self, graph, tissue, size=100):
        if self.debug:
            self.trace['edges_from_cache'] = []
        eready = 0
        prev_ready = -1
        prg = progress.Progress(graph.ecount(),
                                'Loading tissue specificity data for %s' %
                                self.tissnm[tissue], 1)
        inum = 0
        while eready > prev_ready:
            inum += 1
            prg.step(eready - prev_ready, status='iteration %u' % inum)
            entlst = []
            if self.debug:
                self.trace['input_edges'] = []
            for e in graph.es:
                value = None
                if e['giant'][self.tissnm[tissue]] is None:
                    value = self.from_cache(graph.vs[e.source]['name'],
                                            graph.vs[e.target]['name'], tissue)
                if value is None:
                    if self.ventrez[e.source] is not None and self.ventrez[
                            e.source] not in entlst:
                        entlst.append(self.ventrez[e.source])
                    if self.ventrez[e.target] is not None and self.ventrez[
                            e.target] not in entlst:
                        entlst.append(self.ventrez[e.target])
                    #print 'length of query list: %u' % len(entlst)
                    if self.debug:
                        self.trace['input_edges'].append(e.index)
                    if len(entlst) > size or (e.index == graph.ecount() - 1 and
                                              len(entlst) > 0):
                        #print '\t:: starting query'
                        if self.debug:
                            self.trace['entrez_list'] = entlst
                        prg.step(
                            eready - prev_ready,
                            status='downloading (waiting for remote), iter: %u'
                            % inum)
                        data = self.query(tissue, entlst)
                        genes = self.genes_dict(data)
                        prg.step(
                            eready - prev_ready,
                            status='processing data, iter: %u' % inum)
                        if self.debug:
                            self.trace['edges_set'] = []
                        for gedge in data['edges']:
                            # here we get the Gene IDs from the list indexes:
                            epair = (str(genes[gedge['source']]),
                                     str(genes[gedge['target']]))
                            repair = tuple(reversed(list(epair)))
                            eid = self.eentrez[epair] if epair in self.eentrez \
                                else self.eentrez[repair] if repair in self.eentrez else None
                            if eid is not None:
                                # print '\tedge %s---%s [%s] in the network with id = %u; setting value to %s' % (
                                # graph.vs[self.ventrez.index(epair[0])]['name'],
                                # graph.vs[self.ventrez.index(epair[1])]['name'],
                                # str(epair), eid, str(gedge['weight']))
                                graph.es[eid]['giant'][self.tissnm[
                                    tissue]] = gedge['weight']
                                self.to_cache([
                                    graph.vs[graph.es[eid].source]['name'],
                                    graph.vs[graph.es[eid].target]['name'],
                                    gedge['weight']
                                ], tissue)
                                if self.debug:
                                    self.trace['edges_set'].append(eid)
                        prev_ready = eready
                        entlst = []
                else:
                    e['giant'][self.tissnm[tissue]] = value
                    if self.debug:
                        self.trace['edges_from_cache'].append(e.index)
                eready = len([
                    e for e in graph.es
                    if e['giant'][self.tissnm[tissue]] is not None
                ])
        prg.terminate()
