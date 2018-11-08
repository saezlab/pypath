#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2018
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from future.utils import iteritems

import sys
import os

try:
    from twisted.web import server, resource
    from twisted.internet import reactor
except:
    sys.stdout.write('\t:: No `twisted` available.\n')

import urllib
import json

import pandas as pd
import numpy as np

import pypath.descriptions as descriptions
import pypath._html as _html
import pypath.urls as urls
from pypath.common import flatList
from pypath._version import __version__

if 'unicode' not in __builtins__:
    unicode = str


def stop_server():
    
    reactor.removeAll()


class BaseServer(resource.Resource):
    
    def __init__(self):
        
        sys.stdout.write('BaseServer initialized\n')
        
        self.htmls = ['info', '']
        self.welcome_message = (
            'Hello, this is the REST service of pypath %s. Welcome!\n'\
            'For the descriptions of pathway resources go to `/info`.' % (
                __version__
            )
        )
        
        self.isLeaf = True
        resource.Resource.__init__(self)
    
    def render_GET(self, request):
        
        response = []
        
        request.postpath = [i.decode('utf-8') for i in request.postpath]
        
        html = len(request.postpath) == 0 or request.postpath[0] in self.htmls
        self._set_defaults(request, html=html)
        
        if (
            request.postpath and
            hasattr(self, request.postpath[0]) and
            request.postpath[0][0] != '_'
        ):
            
            self._process_postpath(request)
            
            toCall = getattr(self, request.postpath[0])
            
            if hasattr(toCall, '__call__'):
                
                response = toCall(request)
                response = (
                    response.encode('utf-8')
                    if type(response) is unicode else
                    response
                )
                response = [response]
            
        elif not request.postpath:
            
            response = [self.root(request)]
        
        if not response:
            
            response = [
                (
                    "Not found: %s%s" % (
                        '/'.join(request.postpath),
                        ''
                        if len(request.args) == 0 else
                        '?%s' %
                            '&'.join([
                                '%s=%s' % (
                                    k.decode('utf-8'),
                                    v[0].decode('utf-8')
                                )
                                for k, v in iteritems(request.args)
                                if v
                            ])
                    )
                ).encode('utf-8')
            ]
        
        request.write(response[0])
        
        request.finish()
        
        return server.NOT_DONE_YET
    
    def render_POST(self, request):
        
        if request.getHeader(b'content-type').startsWith(b'application/json'):
            
            args_raw = json.loads(request.content.getvalue())
            request.args = dict(
                (
                    k.encode('utf-8'),
                    [v.encode('utf-8')]
                    if type(v) is not list else
                    [','.join(v).encode('utf-8')]
                )
                for k, v in iteritems(args_raw)
            )
        
        return self.render_GET(request)
    
    def _set_defaults(self, request, html=False):
        
        for k, v in iteritems(request.args):
            
            request.args[k] = [b','.join(v)]
        
        request.setHeader('Cache-Control', 'Public')
        
        if '' in request.postpath:
            request.postpath.remove('')
        
        request.setHeader('Access-Control-Allow-Origin', '*')
        
        if html:
            request.setHeader('Content-Type', 'text/html; charset=utf-8')
        elif (
            b'format' in request.args and
            request.args[b'format'][0] == b'json'
        ):
            request.setHeader(
                'Content-Type',
                'application/json'
            )
        else:
            request.args[b'format'] = [b'text']
            request.setHeader('Content-Type', 'text/plain; charset=utf-8')
        
        request.args[b'header'] = [b'1'] if b'header' not in request.args \
            else request.args[b'header']
        
        request.args[b'fields'] = [] if b'fields' not in request.args \
            else request.args[b'fields']
    
    def _process_postpath(self, req):
        
        if len(req.postpath) > 1:
            
            ids_left = [req.postpath[1].encode('utf-8')]
            
            ids_right = (
                [req.postpath[2].encode('utf-8')]
                if (
                    len(req.postpath) > 2 and
                    req.postpath[2].lower() not in {'and', 'or'}
                ) else
                None
            )
            
            left_right = (
                [b'OR']
                if req.postpath[-1].lower() not in {'and', 'or'} else
                [req.postpath[-1].encode('utf-8')]
            )
            
            if ids_right:
                
                if req.postpath[0] == 'ptms':
                    
                    req.args[b'enzymes'] = ids_left
                    req.args[b'substrates'] = ids_right
                    
                else:
                    req.args[b'sources'] = ids_left
                    req.args[b'targets'] = ids_right
                
            else:
                req.args[b'partners'] = ids_left
            
            if req.postpath[0] == 'ptms':
                req.args[b'enzyme_substrate'] = left_right
            else:
                req.args[b'source_target'] = left_right
    
    def about(self, req):
        return self.welcome_message

    def info(self, req):
        return descriptions.gen_html()

    def root(self, req):
        return _html.main_page()
    
    def _parse_arg(self, arg):
        
        if type(arg) is list and len(arg):
            arg = arg[0]
        if hasattr(arg, 'decode'):
            arg = arg.decode('utf-8')
        if hasattr(arg, 'isdigit') and arg.isdigit():
            arg = int(arg)
        if arg == 'no':
            arg = False
        if arg == 'yes':
            arg = True
        
        return bool(arg)


class TableServer(BaseServer):
    
    list_fields = {
        'sources',
        'references',
        'isoforms'
    }
    
    int_list_fields = {
        'references',
        'isoforms'
    }
    
    args_reference = {
        'interactions': {
            'header': None,
            'format': {
                'json',
                'tab',
                'text',
                'tsv',
                'table'
            },
            'datasets': {
                'omnipath',
                'tfregulons',
                'kinaseextra',
                'mirnatarget'
            },
            'types': {
                'PPI',
                'TF',
                'MTI'
            },
            'sources':  None,
            'targets':  None,
            'partners': None,
            'genesymbols': {'1', '0', 'no', 'yes'},
            'fields': {
                'references',
                'sources',
                'tfregulons_level',
                'tfregulons_curated',
                'tfregulons_chipseq',
                'tfregulons_tfbs',
                'tfregulons_coexp',
                'type',
                'ncbi_tax_id',
                'databases',
                'organism'
            },
            'tfregulons_levels':  {'A', 'B', 'C', 'D', 'E'},
            'tfregulons_methods': {
                'curated',
                'chipseq',
                'coexp',
                'tfbs'
            },
            'organisms': {
                '9606',
                '10090',
                '10116'
            },
            'databases': None,
            'source_target': {
                'AND',
                'OR',
                'and',
                'or'
            },
            'directed': {'1', '0', 'no', 'yes'},
            'signed': {'1', '0', 'no', 'yes'},
        },
        'ptms': {
            'header':      None,
            'format': {
                'json',
                'tab',
                'text',
                'tsv',
                'table'
            },
            'enzymes':     None,
            'substrates':  None,
            'partners':    None,
            'genesymbols': {'1', '0', 'no', 'yes'},
            'organisms': {
                '9606',
                '10090',
                '10116'
            },
            'databases': None,
            'residues':  None,
            'modification': None,
            'types': None,
            'fields': {
                'sources',
                'references',
                'ncbi_tax_id',
                'organism',
                'databases',
                'isoforms'
            },
            'enzyme_substrate': {
                'AND',
                'OR',
                'and',
                'or'
            }
        }
    }
    
    datasets_ = {'omnipath', 'tfregulons', 'kinaseextra', 'mirnatarget'}
    tfregulons_methods = {'curated', 'coexp', 'chipseq', 'tfbs'}
    dataset2type = {
        'omnipath': 'PPI',
        'tfregulons': 'TF',
        'kinaseextra': 'PPI',
        'mirnatarget': 'MTI'
    }
    interaction_fields = {
        'references', 'sources', 'tfregulons_level',
        'tfregulons_curated', 'tfregulons_chipseq',
        'tfregulons_tfbs', 'tfregulons_coexp', 'type',
        'ncbi_tax_id', 'databases', 'organism'
    }
    ptms_fields = {
        'references', 'sources', 'databases',
        'isoforms', 'organism', 'ncbi_tax_id'
    }
    
    def __init__(self, a = 'hey', tbls = {
            'interactions': 'omnipath_webservice_interactions.tsv',
            'ptms': 'omnipath_webservice_ptms.tsv'
        }):
        
        sys.stdout.write('TableServer initialized\n')
        
        self.tbls = tbls
        self.data = {}
        self._read_tables()
        self._preprocess_interactions()
        self._preprocess_ptms()
        
        BaseServer.__init__(self)
    
    def _read_tables(self):
        
        sys.stdout.write('Loading data tables\n')
        
        for name, fname in iteritems(self.tbls):
            
            if not os.path.exists(fname):
                
                sys.stdout.write('\t:: Server: missing table: `%s`\n' % (
                    fname
                ))
                continue
            
            self.data[name] = pd.DataFrame.from_csv(
                fname,
                sep = '\t',
                index_col = None
            )
            
            sys.stdout.write(
                'Table `%s` loaded from file `%s`\n' % (name, fname)
            )
    
    def _network(self, req):
        
        hdr = ['nodes', 'edges', 'is_directed', 'sources']
        tbl = self.data['network'].field
        val = dict(zip(tbl.field, tbl.value))
        
        if b'format' in req.args and req.args[b'format'] == b'json':
            return json.dumps(val)
        else:
            return '%s\n%s' % ('\t'.join(hdr), '\t'.join(
                [str(val[h]) for h in hdr]))
    
    def _preprocess_interactions(self):
        
        sys.stdout.write('Preprocessing interactions\n')
        tbl = self.data['interactions']
        tbl['set_sources'] = pd.Series(
            [set(s.split(';')) for s in tbl.sources]
        )
        tbl['set_tfregulons_level'] = pd.Series(
            [
                set(s.split(';'))
                if not pd.isnull(s) else
                set([])
                for s in tbl.tfregulons_level
            ]
        )
    
    def _preprocess_ptms(self):
        
        sys.stdout.write('Preprocessing ptms\n')
        tbl = self.data['ptms']
        tbl['set_sources'] = pd.Series(
            [set(s.split(';')) for s in tbl.sources]
        )
    
    def _check_args(self, req):
        
        result = []
        ref = self.args_reference[req.postpath[0]]
        
        for arg, val in iteritems(req.args):
            
            arg = arg.decode('utf-8')
            
            if arg in ref:
                
                if not ref[arg] or not val:
                    
                    continue
                
                val = (
                    {val[0]}
                    if type(val[0]) is int else
                    set(val[0].decode('utf-8').split(','))
                )
                
                unknowns = val - ref[arg]
                
                if unknowns:
                    
                    result.append(
                        ' ==> Unknown values for argument `%s`: `%s`' % (
                            arg,
                            ', '.join(str(u) for u in unknowns)
                        )
                    )
                
            else:
                
                result.append(' ==> Unknown argument: `%s`' % arg)
        
        req.args[b'header'] = self._parse_arg(req.args[b'header'])
        
        if result:
            
            return (
                'Something is not entirely good:\n%s\n\n'
                'Please check the examples at\n'
                'https://github.com/saezlab/pypath\n'
                'and\n'
                'https://github.com/saezlab/DoRothEA\n'
                'If you still experiencing issues contact us at\n'
                'omnipath@googlegroups.com' % '\n'.join(result)
            )
    
    def queries(self, req):
        
        query_type = (
            req.postpath[1]
                if len(req.postpath) > 1 else
            'interactions'
        )
        
        query_param = (
            req.postpath[2]
                if len(req.postpath) > 2 else
            None
        )
        
        if query_type in self.args_reference:
            
            result = self.args_reference[query_type]
            
            if query_param is not None and query_param in result:
                
                result = {}
                result[query_param] = (
                    self.args_reference[query_type][query_param]
                )
            
        else:
            
            result = {}
            result[query_type] = (
                'No possible arguments defined for'
                'query `%s` or no such query available.' % query_type
            )
        
        if b'format' in req.args and req.args[b'format'][0] == b'json':
            
            return json.dumps(result)
            
        else:
            
            return 'argument\tvalues\n%s' % '\n'.join(
                '%s\t%s' % (
                    k,
                    ';'.join(v) if type(v) in {list, set} else str(v)
                )
                for k, v in iteritems(result)
            )
    
    def databases(self, req):
        
        query_type = (
            req.postpath[1]
                if len(req.postpath) > 1 else
            'interactions'
        )
        
        datasets = (
            set(req.postpath[2].split(','))
                if len(req.postpath) > 2 else
            None
        )
        
        tbl = (
            self.data[query_type]
                if query_type in self.data else
            self.data['interactions']
        )
        
        # filter for datasets
        if query_type == 'interactions':
            
            if datasets is not None:
                
                tbl = tbl[tbl.type.isin(datasets)]
                
            else:
                
                datasets = self._get_datasets()
            
            result = {}
            
            for dataset in datasets:
                
                result[dataset] = sorted(set.union(
                    *tbl[tbl.type == dataset].set_sources)
                )
            
        else:
            
            result = {}
            result['*'] = sorted(set.union(*tbl.set_sources))
        
        if b'format' in req.args and req.args[b'format'][0] == b'json':
            
            return json.dumps(result)
            
        else:
            
            return 'dataset\tdatabases\n%s' % '\n'.join(
                '%s\t%s' % (k, ';'.join(v)) for k, v in iteritems(result)
            )
    
    def _get_datasets(self):
        
        return list(self.data['interactions'].type.unique())
    
    def datasets(self, req):
        
        query_type = (
            req.postpath[1]
                if len(req.postpath) > 1 else
            'interactions'
        )
        
        if query_type == 'interactions':
            
            result = self._get_datasets()
        
        else:
            
            result = []
        
        if b'format' in req.args and req.args[b'format'][0] == b'json':
            
            return json.dumps(result)
            
        else:
            
            return ';'.join(result)
    
    def interactions(
            self,
            req,
            datasets  = {'omnipath'},
            databases = None,
            tfregulons_levels = {'A', 'B'},
            organisms = {9606},
            source_target = 'OR'
        ):
        
        bad_req = self._check_args(req)
        
        if bad_req:
            
            return bad_req
        
        hdr = [
            'source', 'target', 'is_directed', 'is_stimulation',
            'is_inhibition', 'dip_url'
        ]
        
        if b'source_target' in req.args:
            
            source_target = (
                req.args[b'source_target'][0].decode('utf-8').upper()
            )
        
        args = {}
        
        for arg in (
            'datasets', 'types', 'tfregulons_levels',
            'sources', 'targets', 'partners', 'databases',
            'tfregulons_methods', 'organisms'
        ):
            
            args[arg] = self._args_set(req, arg)
        
        # if user requested TF type interactions
        # they likely want the tfregulons dataset
        if 'TF' in args['types']:
            args['datasets'].add('tfregulons')
        if 'MTI' in args['types']:
            args['datasets'].add('mirnatarget')
        
        # here adjust on the defaults otherwise we serve empty
        # response by default
        args['datasets'] = args['datasets'] or datasets
        args['datasets'] = args['datasets'] & self.datasets_
        
        args['organisms'] = set(
            int(t) for t in args['organisms'] if t.isdigit()
        )
        args['organisms'] = args['organisms'] or organisms
        
        # do not allow impossible values
        # those would result KeyError later
        args['tfregulons_levels'] = (
            args['tfregulons_levels'] or tfregulons_levels
        )
        args['tfregulons_methods'] = (
            args['tfregulons_methods'] & self.tfregulons_methods
        )
        
        # provide genesymbols: yes or no
        if (
            b'genesymbols' in req.args and
            self._parse_arg(req.args[b'genesymbols'])
        ):
            genesymbols = True
            hdr.insert(2, 'source_genesymbol')
            hdr.insert(3, 'target_genesymbol')
        else:
            genesymbols = False
        
        # if user requested TF Regulons they likely want us
        # to serve TF-target interactions
        # but if they requested other types, then we
        # serve those as well
        if 'tfregulons' in args['datasets']:
            args['types'].add('TF')
        if 'mirnatarget' in args['datasets']:
            args['types'].add('MTI')
        
        # if no types provided we collect the types
        # for the datasets requested
        # or by default only the 'omnipath' dataset
        # which belongs to the 'PPI' type
        if not args['types'] or args['datasets']:
            args['types'].update(set(
                self.dataset2type[ds] for ds in args['datasets']
            ))
        
        # starting from the entire dataset
        tbl = self.data['interactions']
        
        # filter by type
        tbl = tbl[tbl.type.isin(args['types'])]
        
        # if partners provided those will overwrite
        # sources and targets
        args['sources'] = args['sources'] or args['partners']
        args['targets'] = args['targets'] or args['partners']
        
        # then we filter by source and target
        # which matched against both standard names
        # and gene symbols
        if args['sources'] and args['targets'] and source_target == 'OR':
            
            tbl = tbl[
                tbl.target.isin(args['targets']) |
                tbl.target_genesymbol.isin(args['targets']) |
                tbl.source.isin(args['sources']) |
                tbl.source_genesymbol.isin(args['sources'])
            ]
        
        else:
            
            if args['sources']:
                tbl = tbl[
                    tbl.source.isin(args['sources']) |
                    tbl.source_genesymbol.isin(args['sources'])
                ]
            
            if args['targets']:
                tbl = tbl[
                    tbl.target.isin(args['targets']) |
                    tbl.target_genesymbol.isin(args['targets'])
                ]
        
        # filter by datasets
        if args['datasets']:
            tbl = tbl.query(' or '.join(args['datasets']))
        
        # filter by organism
        tbl = tbl[
            tbl.ncbi_tax_id_source.isin(args['organisms']) |
            tbl.ncbi_tax_id_target.isin(args['organisms'])
        ]
        
        # filter by TG Regulons confidence levels
        if 'TF' in args['types'] and args['tfregulons_levels']:
            
            tbl = tbl[
                np.logical_not(tbl.tfregulons) |
                (tbl.set_tfregulons_level & args['tfregulons_levels'])
            ]
        
        # filter by databases
        if args['databases']:
            
            tbl = tbl[tbl.set_sources & args['databases']]
        
        # filtering by TF Regulons methods
        if 'TF' in args['types'] and args['tfregulons_methods']:
            
            q = ['tfregulons_%s' % m for m in args['tfregulons_methods']]
            
            tbl = tbl[
                tbl[q].any(1) | np.logical_not(tbl.tfregulons)
            ]
        
        # filter directed & signed
        if (
            b'directed' in req.args and
            self._parse_arg(req.args[b'directed'])
        ):
            
            tbl = tbl[tbl.is_directed == 1]
        
        if (
            b'signed' in req.args and
            self._parse_arg(req.args[b'signed'])
        ):
            
            tbl = tbl[np.logical_or(
                tbl.is_stimulation == 1,
                tbl.is_inhibition == 1
            )]
        
        if req.args[b'fields']:
            
            _fields = [
                f for f in
                req.args[b'fields'][0].decode('utf-8').split(',')
                if f in self.interaction_fields
            ]
            
            for f in _fields:
                
                if f == 'ncbi_tax_id' or f == 'organism':
                    
                    hdr.append('ncbi_tax_id_source')
                    hdr.append('ncbi_tax_id_target')
                    
                elif f == 'databases':
                    
                    hdr.append('sources')
                    
                else:
                    
                    hdr.append(f)
        
        tbl = tbl.loc[:,hdr]
        
        return self._serve_dataframe(tbl, req)
    
    def ptms(
            self,
            req,
            organisms = {9606},
            enzyme_substrate = 'OR'
        ):
        
        bad_req = self._check_args(req)
        
        if bad_req:
            
            return bad_req
        
        hdr = [
            'enzyme', 'substrate', 'residue_type',
            'residue_offset', 'modification'
        ]
        
        if b'enzyme_substrate' in req.args:
            
            enzyme_substrate = (
                req.args[b'enzyme_substrate'][0].decode('utf-8').upper()
            )
        
        args = {}
        
        for arg in (
            'enzymes', 'substrates', 'partners',
            'databases', 'organisms', 'types',
            'residues'
        ):
            
            args[arg] = self._args_set(req, arg)
        
        args['organisms'] = set(
            int(t) for t in args['organisms'] if t.isdigit()
        )
        args['organisms'] = args['organisms'] or organisms
        
        # provide genesymbols: yes or no
        if (
            b'genesymbols' in req.args and
            self._parse_arg(req.args[b'genesymbols'])
        ):
            genesymbols = True
            hdr.insert(2, 'enzyme_genesymbol')
            hdr.insert(3, 'substrate_genesymbol')
        else:
            genesymbols = False
        
        # starting from the entire dataset
        tbl = self.data['ptms']
        
        # filter by type
        if args['types']:
            tbl = tbl[tbl.modification.isin(args['types'])]
        
        # if partners provided those will overwrite
        # enzymes and substrates
        args['enzymes'] = args['enzymes'] or args['partners']
        args['substrates'] = args['substrates'] or args['partners']
        
        # then we filter by enzyme and substrate
        # which matched against both standard names
        # and gene symbols
        if args['enzymes'] and args['substrates'] and enzyme_substrate == 'OR':
            
            tbl = tbl[
                tbl.substrate.isin(args['substrates']) |
                tbl.substrate_genesymbol.isin(args['substrates']) |
                tbl.enzyme.isin(args['enzymes']) |
                tbl.enzyme_genesymbol.isin(args['enzymes'])
            ]
        
        else:
            
            if args['enzymes']:
                tbl = tbl[
                    tbl.enzyme.isin(args['enzymes']) |
                    tbl.enzyme_genesymbol.isin(args['enzymes'])
                ]
            
            if args['substrates']:
                tbl = tbl[
                    tbl.substrate.isin(args['substrates']) |
                    tbl.substrate_genesymbol.isin(args['substrates'])
                ]
        
        # filter by organism
        tbl = tbl[tbl.ncbi_tax_id.isin(args['organisms'])]
        
        # filter by databases
        if args['databases']:
            
            tbl = tbl[tbl.set_sources & args['databases']]
        
        if req.args[b'fields']:
            
            _fields = [
                f for f in
                req.args[b'fields'][0].decode('utf-8').split(',')
                if f in self.ptms_fields
            ]
            
            for f in _fields:
                
                if f == 'ncbi_tax_id' or f == 'organism':
                    
                    hdr.append('ncbi_tax_id')
                    
                elif f == 'databases':
                    
                    hdr.append('sources')
                    
                else:
                    
                    hdr.append(f)
        
        tbl = tbl.loc[:,hdr]
        
        return self._serve_dataframe(tbl, req)
    
    @classmethod
    def _serve_dataframe(cls, tbl, req):
        
        if b'format' in req.args and req.args[b'format'][0] == b'json':
            
            data_json = tbl.to_json(orient = 'records')
            # this is necessary because in the data frame we keep lists
            # as `;` separated strings but in json we is nicer to serve
            # them as lists
            data_json = json.loads(data_json)
            
            for i in data_json:
                
                for k, v in iteritems(i):
                    
                    if k in cls.list_fields:
                        
                        i[k] = [
                            int(f) if k in cls.int_list_fields else f
                            for f in v.split(';')
                        ]
            
            return json.dumps(data_json)
            
        else:
            
            return tbl.to_csv(
                sep = '\t',
                index = False,
                header = bool(req.args[b'header'])
            )
    
    @staticmethod
    def _args_set(req, arg):
        
        arg = arg.encode('utf-8')
        
        return (
            set(req.args[arg][0].decode('utf-8').split(','))
            if arg in req.args
            else set([])
        )


class PypathServer(BaseServer):
    
    def __init__(self, pypath):
        self.p = pypath
        self.g = pypath.graph
        self.isLeaf = True
        
        BaseServer.__init__(self)
    
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
                f.decode('utf-8') for f in fields if f in req.args[b'fields']
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
                            for r in flatList([
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
                f.decode('utf-8') for f in fields if f in req.args[b'fields']
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


class Rest(object):
    
    def __init__(self, port, serverclass = PypathServer, **kwargs):
        """
        Runs a webserver serving a `PyPath` instance listening
        to a custom port.
        
        Args:
        -----
        :param int port:
            The port to listen to.
        :param str serverclass'
            The class implementing the server.
        :param **kwargs:
            Arguments for initialization of the server class.
        """
        
        self.port = port
        self.site = server.Site(serverclass(kwargs))
        reactor.listenTCP(self.port, self.site)
        reactor.run()
