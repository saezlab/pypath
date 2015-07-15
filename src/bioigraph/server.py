#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `bioigraph` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from twisted.web import server, resource
from twisted.internet import reactor
import urllib
import json
from common import *

class Rest(object):
    
    def __init__(self, biograph, port):
        self.port = port
        self.site = server.Site(RestResource(biograph))
        reactor.listenTCP(self.port, self.site)
        reactor.run()

class RestResource(resource.Resource):
    
    def __init__(self, biograph):
        self.b = biograph
        self.g = biograph.graph
        self.isLeaf = True
    
    def render_GET(self, request):
        self.set_defaults(request)
        if len(request.postpath) > 0 and hasattr(self, request.postpath[0]):
            toCall = getattr(self, request.postpath[0])
            if hasattr(toCall, '__call__'):
                return str(toCall(request))
        elif len(request.postpath) == 0:
            return self.root(request)
        return str(request.__dict__)
        return "Not found: %s%s" % ('/'.join(request.postpath), 
            '' if len(request.args) == 0 else \
            '?%s' % '&'.join(['%s=%s'%(k, v) for k, v in request.args.iteritems()]))
    
    def set_defaults(self, request):
        request.setHeader('Cache-Control', 'Public')
        if '' in request.postpath:
            request.postpath.remove('')
        if 'format' in request.args and request.args['format'][0] == 'json':
            request.args['format'] = 'json'
            request.setHeader('Content-Type', 'text/json')
        else:
            request.args['format'] = 'text'
            request.setHeader('Content-Type', 'text/plain')
        request.args['header'] = 1 if 'header' not in request.args \
            else int(request.args['header'][0])
        request.args['fields'] = [] if 'fields' not in request.args \
            else request.args['fields']
    
    def about(self, req):
        return 'Hello, this is the REST service of bioigraph %s. Welcome!' % __version__
    
    def root(self, req):
        hdr = ['nodes', 'edges', 'is_directed', 'sources']
        val = [self.g.vcount(), self.g.ecount(), 
            int(self.g.is_directed()), self.b.sources]
        if req.args['format'] == 'json':
            return json.dumps(dict(zip(hdr, val)))
        else:
            return '%s\n%s' % ('\t'.join(hdr), '\t'.join([str(v) \
                if type(v) is not list else ';'.join(v) for v in val]))
    
    def interactions(self, req):
        fields = ['sources', 'references']
        result = []
        elist = self._get_eids(req)
        res = []
        hdr = ['source', 'target', 'is_directed', 'is_stimulation', 'is_inhibition']
        hdr += [f for f in fields if f in req.args['fields']]
        for eid in elist:
            e = self.g.es[eid]
            for d in ['straight', 'reverse']:
                uniprots = getattr(e['dirs'], d)
                if e['dirs'].dirs[uniprots]:
                    thisEdge = [
                        uniprots[0],
                        uniprots[1],
                        1,
                        int(e['dirs'].is_stimulation(uniprots)),
                        int(e['dirs'].is_inhibition(uniprots))
                    ]
                    dsources = e['dirs'].get_dir(uniprots, sources = True)
                    dsources = dsources | e['dirs'].get_dir('undirected', sources = True)
                    if 'sources' in hdr:
                        thisEdge.append(list(dsources))
                    if 'references' in hdr:
                        thisEdge.append([r.pmid for r in flatList([rs for s, rs in \
                            e['refs_by_source'].iteritems() \
                            if s in dsources])])
                    res.append(thisEdge)
            if not e['dirs'].is_directed():
                thisEdge = [e['dirs'].nodes[0], e['dirs'].nodes[1], 0, 0, 0]
                if 'sources' in hdr:
                    thisEdge.append(e['sources'])
                if 'references' in hdr:
                    thisEdge.append([r.pmid for r in e['references']])
                res.append(thisEdge)
        if req.args['format'] == 'json':
            return json.dumps([dict(zip(hdr, r)) for r in res])
        else:
            return self._table_output(res, hdr, req)
    
    def _table_output(self, res, hdr, req):
        return '%s%s' % ('' if not bool(req.args['header']) else '%s\n'%'\t'.join(hdr), 
                '\n'.join(['\t'.join([';'.join(f) if type(f) is list else str(f) \
                for f in r]) for r in res]))
    
    def _get_eids(self, req):
        names = None if len(req.postpath) <= 1 else req.postpath[1].split(',')
        ids = range(0, self.g.vcount()) if names is None else self.b.names2vids(names)
        ids = set(ids) 
        elist = set(range(0, self.g.ecount())) if names is None else set([])
        if names is not None:
            alist = dict(zip(ids, [self.g.neighbors(i) for i in ids]))
            for one, others in alist.iteritems():
                for two in others:
                    e = self.g.get_eid(one, two, directed = True, error = False)
                    if e != -1:
                        elist.add(e)
                    if self.g.is_directed():
                        e = self.g.get_eid(two, one, directed = True, error = False)
                        if e != -1:
                            elist.add(e)
        return elist
    
    def ptms(self, req):
        fields = ['is_stimulation', 'is_inhibition', 'sources', 'references']
        result = []
        elist = self._get_eids(req)
        res = []
        hdr = ['enzyme', 'substrate', 'residue_type', 'residue_offset', 'modification']
        hdr += [f for f in fields if f in req.args['fields']]
        if 'ptm' in self.g.es.attributes():
            for eid in elist:
                e = self.g.es[eid]
                for ptm in e['ptm']:
                    thisPtm = [ptm.domain.protein, ptm.ptm.protein, ptm.ptm.residue.name, 
                        ptm.ptm.residue.number, ptm.ptm.typ]
                    if 'sources' in hdr:
                        thisPtm.append(ptm.ptm.sources)
                    if 'references' in hdr:
                        thisPtm.append(ptm.refs)
                    if 'is_stimulation' in hdr:
                        thisPtm.append(int(edge['dirs'].is_stimulation(\
                            (ptm.domain.protein, ptm.ptm.protein))))
                    if 'is_inhibition' in hdr:
                        thisPtm.append(int(edge['dirs'].is_stimulation(\
                            (ptm.domain.protein, ptm.ptm.protein))))
                    res.append(thisPtm)
        if req.args['format'] == 'json':
            return json.dumps([dict(zip(hdr, r)) for r in res])
        else:
            return self._table_output(res, hdr, req)