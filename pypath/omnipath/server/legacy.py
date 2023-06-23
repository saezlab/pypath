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

"""
Old server class working from a ``pypath.legacy.main.PyPath`` object.
"""

import json

import pypath.omnipath.server.run as server


class PypathServer(server.BaseServer):
    
    
    def __init__(self, pypath):
        
        self.p = pypath
        self.g = pypath.graph
        self.isLeaf = True
        
        server.BaseServer.__init__(self)
    
    
    def network(self, req):
        
        hdr = ['nodes', 'edges', 'is_directed', 'sources']
        val = [
            self.g.vcount(), self.g.ecount(), int(self.g.is_directed()),
            self.p.sources
        ]
        if b'format' in req.args and req.args[b'format'] == b'json':
            return json.dumps(dict(zip(hdr, val)))
        else:
            return '%s\n%s' % ('\t'.join(hdr), '\t'.join(
                [str(v) if type(v) is not list else ';'.join(v) for v in val]))


    def interactions(self, req):
        
        fields = [b'sources', b'references']
        result = []
        elist = self._get_eids(req)
        res = []
        hdr = [
            'source', 'target', 'is_directed', 'is_stimulation',
            'is_inhibition'
        ]
        
        if b'fields' in req.args:
            
            hdr += [
                f.decode('utf-8')
                for f in fields
                if f in b','.join(req.args[b'fields']).split(b',')
            ]
        
        if (
            b'genesymbols' in req.args and
            self._parse_arg(req.args[b'genesymbols'])
        ):
            genesymbols = True
            hdr.insert(2, 'source_genesymbol')
            hdr.insert(3, 'target_genesymbol')
        else:
            genesymbols = False
        
        all_sources = set([])
        for eid in elist:
            e = self.g.es[eid]
            all_sources = all_sources | e['sources']
            for d in ['straight', 'reverse']:
                uniprots = getattr(e['dirs'], d)
                
                if e['dirs'].dirs[uniprots]:
                    
                    thisEdge = [
                        uniprots[0], uniprots[1]
                    ]
                    
                    if genesymbols:
                        
                        thisEdge.extend([
                            self.g.vs[self.p.nodDct[uniprots[0]]]['label'],
                            self.g.vs[self.p.nodDct[uniprots[1]]]['label']
                        ])
                    
                    thisEdge.extend([
                        1,
                        int(e['dirs'].is_stimulation(uniprots)),
                        int(e['dirs'].is_inhibition(uniprots))
                    ])
                    dsources = e['dirs'].get_dir(uniprots, sources=True)
                    dsources = dsources | e['dirs'].get_dir(
                        'undirected', sources=True)
                    if 'sources' in hdr:
                        thisEdge.append(list(dsources))
                    if 'references' in hdr:
                        thisEdge.append([
                            r.pmid
                            for r in flat_list([
                                rs for s, rs in iteritems(e['refs_by_source'])
                                if s in dsources
                            ])
                        ])
                    thisEdge.append(self._dip_urls(e))
                    res.append(thisEdge)
            
            if not e['dirs'].is_directed():
                
                thisEdge = [e['dirs'].nodes[0], e['dirs'].nodes[1]]
                if genesymbols:
                    
                    thisEdge.extend([
                        self.g.vs[self.p.nodDct[e['dirs'].nodes[0]]]['label'],
                        self.g.vs[self.p.nodDct[e['dirs'].nodes[1]]]['label']
                    ])
                
                thisEdge.extend([0, 0, 0])
                if 'sources' in hdr:
                    thisEdge.append(list(e['sources']))
                if 'references' in hdr:
                    thisEdge.append([r.pmid for r in e['references']])
                thisEdge.append(self._dip_urls(e))
                res.append(thisEdge)
        
        if 'DIP' in all_sources:
            hdr.append('dip_url')
        else:
            res = map(lambda r: r[:-1], res)
        if b'format' in req.args and req.args[b'format'] == b'json':
            return json.dumps([dict(zip(hdr, r)) for r in res])
        else:
            return self._table_output(res, hdr, req)
    
    
    
    def _table_output(self, res, hdr, req):
        
        return '%s%s' % ('' if not bool(req.args[b'header']) else
                         '%s\n' % '\t'.join(hdr), '\n'.join([
                             '\t'.join([
                                 ';'.join(f) if type(f) is list or
                                 type(f) is set else str(f) for f in r
                             ]) for r in res
                         ]))
    
    
    
    def _get_eids(self, req):
        
        names = None if len(req.postpath) <= 1 else req.postpath[1].split(',')
        ids = range(
            0, self.g.vcount()) if names is None else self.p.names2vids(names)
        ids = set(ids)
        elist = set(range(0, self.g.ecount())) if names is None else set([])
        if names is not None:
            alist = dict(zip(ids, [self.g.neighbors(i) for i in ids]))
            for one, others in iteritems(alist):
                for two in others:
                    e = self.g.get_eid(one, two, directed=True, error=False)
                    if e != -1:
                        elist.add(e)
                    if self.g.is_directed():
                        e = self.g.get_eid(
                            two, one, directed=True, error=False)
                        if e != -1:
                            elist.add(e)
        return elist
    
    
    def ptms(self, req):
        fields = [
            b'is_stimulation', b'is_inhibition', b'sources', b'references'
        ]
        result = []
        elist = self._get_eids(req)
        res = []
        hdr = [
            'enzyme', 'substrate', 'residue_type', 'residue_offset',
            'modification'
        ]
        
        if (
            b'genesymbols' in req.args and
            self._parse_arg(req.args[b'genesymbols'])
        ):
            genesymbols = True
            hdr.insert(2, 'enzyme_genesymbol')
            hdr.insert(3, 'substrate_genesymbol')
        else:
            genesymbols = False
        
        if b'fields' in req.args:
            
            hdr += [
                f.decode('utf-8')
                for f in fields
                if f in b','.join(req.args[b'fields']).split(b',')
            ]
        
        if 'ptm' in self.g.es.attributes():
            for eid in elist:
                e = self.g.es[eid]
                for ptm in e['ptm']:
                    if 'ptmtype' not in req.args or ptm.ptm.typ in req.args[
                            'ptmtype']:
                        thisPtm = [
                            ptm.domain.protein, ptm.ptm.protein
                        ]
                        
                        if genesymbols:
                            
                            thisPtm.extend([
                                self.g.vs[self.p.nodDct[
                                    ptm.domain.protein]]['label'],
                                self.g.vs[self.p.nodDct[
                                    ptm.ptm.protein]]['label']
                            ])
                        
                        thisPtm.extend([
                            ptm.ptm.residue.name, ptm.ptm.residue.number,
                            ptm.ptm.typ
                        ])
                        
                        if 'is_stimulation' in hdr:
                            thisPtm.append(
                                int(e['dirs'].is_stimulation((
                                    ptm.domain.protein, ptm.ptm.protein))))
                        
                        if 'is_inhibition' in hdr:
                            thisPtm.append(
                                int(e['dirs'].is_inhibition((
                                    ptm.domain.protein, ptm.ptm.protein))))
                        
                        if 'sources' in hdr:
                            thisPtm.append(list(ptm.ptm.sources))
                        
                        if 'references' in hdr:
                            thisPtm.append(list(ptm.refs))
                        res.append(thisPtm)
        
        if b'format' in req.args and req.args[b'format'] == b'json':
            return json.dumps([dict(zip(hdr, r)) for r in res])
        else:
            return self._table_output(res, hdr, req)


    def resources(self, req):
        
        hdr = [
            'database', 'proteins', 'interactions', 'directions',
            'stimulations', 'inhibitions', 'signs'
        ]
        res = [[
            # database name
            s,
            # number of proteins
            len([1 for v in self.g.vs if s in v['sources']]),
            # number of interacting pairs
            len([1 for e in self.g.es if s in e['sources']]),
            # number of directions
            sum([s in e['dirs'].sources[d] \
                 for e in self.g.es for d in e['dirs'].which_dirs()]),
            # number of stimulations
            sum([s in e['dirs'].positive_sources[d] \
                 for e in self.g.es for d in e['dirs'].which_dirs()]),
            # number of inhibitions
            sum([s in e['dirs'].negative_sources[d] \
                 for e in self.g.es for d in e['dirs'].which_dirs()]),
            # number of signs
            sum([s in sg for e in self.g.es for d in e['dirs'].which_dirs() \
                 for sg in [e['dirs'].positive_sources[d],
                            e['dirs'].negative_sources[d]]]),
        ]
            for s in self.p.sources]
        if b'format' in req.args and req.args[b'format'] == b'json':
            return json.dumps([dict(zip(hdr, r)) for r in res])
        else:
            return self._table_output(res, hdr, req)


    def _dip_urls(self, e):
        result = []
        if 'dip_id' in e.attributes():
            for dip_id in e['dip_id']:
                try:
                    result.append(urls.urls['dip']['ik'] %
                                int(dip_id.split('-')[1][:-1]))
                except:
                    
                    print(dip_id)
                
        return ';'.join(result)
