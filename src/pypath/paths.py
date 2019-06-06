#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import __main__

from future.utils import iteritems
from past.builtins import xrange, range, reduce

# external modules:
import sys
import imp
from functools import reduce
import itertools
import operator
import warnings

import igraph

class Path(object):
    
    """
    Represents a directed path in an igraph object.
    """
    
    def __init__(self, graph, sequence):
        
        self.graph   = graph
        self.sequence = sequence
        self.signs   = signs
        self._update_members()
    
    def reload(self):
        """
        Reloads the objects preserving the attributes of the instance.
        """
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def __iter__(self):
        """
        Iterates over members of the path.
        """
        for elem in self.sequence:
            yield elem
    
    def __iadd__(self, new):
        """
        Joins one element to the end of this path.
        """
        self.sequence.append(new)
        self._update_members()
    
    def __len__(self):
        """
        Returns the number of elements in the path.
        """
        return len(self.sequence)
    
    def _update_members(self):
        """
        Updates the members set.
        """
        self.members = set(self.sequence)
    
    def iterpairs(self):
        """
        Iterates over the links in the path as pairs of vertices.
        """
        for i in xrange(1, len(self)):
            yield (self.sequence[i - 1], self.sequence[i])
    
    def iterlinks(self):
        """
        Returns the links (edges) in the path.
        """
        
        for i in xrange(len(self) - 1):
            
            yield self.graph.get_eid(self.sequence[i], self.sequence[i+1])
    
    def links(self):
        """
        Returns the links (edges) in the path.
        """
        
        return list(self.iterlinks())
    
    def iternames(self):
        """
        Iterates over the names in the path.
        """
        
        for vi in self.sequence:
            
            yield self.graph[vi]['name']
    
    def names(self):
        """
        Returns the path as sequence of vertex names.
        """
        
        return list(self.iternames())
    
    def iterlabels(self):
        """
        Iterates over the labels in the path.
        """
        
        for vi in self.sequence:
            
            yield self.graph[vi]['label']
    
    def labels(self):
        """
        Returns the path as sequence of vertex labels.
        """
        
        return list(self.iterlabels())


class WeightedPath(Path):
    
    """
    Represents a directed path in an igraph object,
    with edge signs and weights and vertex weights optionally.
    """
    
    sign_notation = {
        'stimulation':  1,
         'activation':  1,
         'inhibition': -1,
                  '+':  1,
                  '-': -1,
                    1:  1,
                   -1: -1,
                    0:  0,
            'unknown':  0,
                  '*':  0
    }
    
    def __init__(self, graph, sequence,
                 signs = None,
                 vweights = None,
                 eweights = None):
        
        super(WeightedPath, self).__init__(graph, sequence)
        
        self.signs    = signs
        self.vweights = vweights
        self.eweights = eweights
        
        self._check()
    
    def has_signs(self):
        """
        Tells whether this path has signs (+/-) assigned to the edges.
        """
        return self.signs is not None
    
    def has_vweights(self):
        """
        Tells whether this path has weights assigned to its vertices.
        """
        return self.vweights is not None
    
    def has_eweights(self):
        """
        Tells whether this path has weights assigned to its vertices.
        """
        return self.vweights is not None
    
    def _check(self):
        """
        Checks if the length of sign and weight vectors matches the
        length of the path. Raises `ValueError` on mismatch.
        """
        if self.signs is not None and len(self.signs) != len(self) - 1:
            
            raise ValueError('Length of `signs` vector must match'\
                             'the number of edges')
        
        if self.eweights is not None and len(self.eweights) != len(self) - 1:
            
            raise ValueError('Length of `eweights` vector must match'\
                             'the number of edges')
        
        if self.vweights is not None and len(self.vweights) != len(self):
            
            raise ValueError('Length of `vweights` vector must match'\
                             'the number of vertices')
    
    def add_dict(self, new):
        """
        Adds a new element to the end of the path.
        New element defined by a dict,
        e.g. `{'id': 23, 'sign': -1, 'vweight': None, 'eweight': 4.6547}`.
        """
        sign    = new['sign'] if 'sign' in new else None
        vweight = new['vweight'] if 'vweight' in new else None
        eweight = new['eweight'] if 'eweight' in new else None
        self.__iadd__((new['id'], sign, vweight, eweight))
    
    def __iadd__(self, new):
        """
        Adds a new element to the end of the path.
        New element defined by a tuple,
        e.g. `(23, -1, None, 4.6547)`.
        Elements in the tuple correspond to id, sign,
        vweight and eweight, respectively.
        """
        
        super(WeightedPath, self).__iadd__(new[0])
        
        if self.signs is not None:
            
            sign = None if len(new) == 1 else new[1]
            self.signs.append(self._sign(sign))
        
        self._add_attr('signs', 1, new, 'Sign', self._sign)
        self._add_attr('vweights', 2, new, 'Vertex weight')
        self._add_attr('eweights', 3, new, 'Edge weight')
        
        self._check()
    
    def _add_attr(self, attr, iattr, new, name, transf = lambda x: x):
        """
        Adds new value to an attribute vector.
        """
        
        attr = getattr(self, attr)
        
        if attr is not None:
            
            val = None \
                if len(new) <= iattr or new[iattr] is None \
                else transf(new[iattr])
            
            attr.append(val)
        
        if len(new) > iattr and new[iattr] is not None and attr is None:
            
            warnings.warn('%s provided,'\
                          'but this path has no %ss' % (name, name.lower()))
    
    @classmethod
    def _sign(cls, sign):
        """
        Transforms signs to integer notation (-1 or +1).
        """
        return cls.sign_notation[sign]
    
    def vscore(self):
        """
        Returns the sum of vertex weights.
        """
        return 0.0 \
            if self.vweights is None \
            else sum(filter(lambda w: w is None, self.vweights))
    
    def escore(self):
        """
        Returns the sum of edge weights.
        """
        return 0.0 \
            if self.eweights is None \
            else sum(filter(lambda w: w is None, self.eweights))
    
    def sign(self):
        """
        Returns the overall sign of the path,
        which is the product of the sign vector.
        """
        return reduce(operator.mul, self.signs)
    
    


class PathwayPath(WeightedPath):
    
    """
    Represents a directed path in a PyPath object.
    """
    
    def __init__(self, pa, sequence):
        
        self.pa = pa
        
        super(PathwayPath, self).__init__(sequence)
    
    def uniprots(self):
        """
        Returns the path as a series of UniProt IDs.
        """
        
        return self.names()
    
    def genesymbols(self):
        """
        Returns the path as a series of GeneSymbols.
        """
        
        return self.labels()

class PathFinder(object):
    
    def __init__(self, graph, sources, targets = None,
                 maxlen = 2, mode = 'ALL'):
        """
        Finds all paths in an igraph.Graph object between arbitrary sets
        of vertices.
        
        :param igraph.Graph graph: The graph to look up paths in.
        :param list sources: The source vertices of the paths.
                             Either vertex indices or names.
        :param list targets: The target vertices.
                             If `None` sources will be used.
        :param int maxlen: The maximum length of the paths in number of steps.
        :param str mode: Direction of the paths: ALL, IN or OUT.
        """
        
        self.graph        = graph
        self.sources      = sources
        self.targets      = targets
        self.maxlen       = maxlen
        self.mode         = mode
        self.all_shortest = all_shortest
        
        self._name_dict()
        
        if self.sources and type(self.sources[0]) is tuple:
            
            self._by_pairs = True
            self._vid_pairs()
        
        else:
            
            self._by_pairs = False
            self._vid_sources()
            
            if self.targets is not None:
                
                self._vid_targets()
    
    def _name_dict(self):
        
        self.names = (
            dict(
                map(
                    lambda n:
                        (n[1], n[0]),
                    iteritems(self.graph['name'])
                )
            )
        )
    
    def _vid_pairs(self):
        
        self._pairs = (
            list(
                map(
                    lambda pair:
                        tuple(map(self._get_vid, pair)),
                    self.sources
                )
            )
        )
    
    def _vid_sources(self):
        
        self._sources = list(map(self._get_vid, self.sources))
    
    def _vid_targets(self):
        
        if self.targets is not None:
            
            self._targets = list(map(self._get_vid, self.targets))
    
    def _get_vid(self, v):
        
        return (
            v             if type(v) is int else
            self.names[v] if v in self.names else
            None
        )
    
    def iterpairs(self):
        """
        Iterates over pairs of vertex indices.
        Yields tuples where the first element is the source
        and the second is the target.
        """
        
        if self._by_pairs:
            
            return self._pairs.__iter__()
            
        elif self.targets is not None:
            
            return itertools.product(self._sources, self._targets)
            
        elif self.mode == 'ALL':
            
            return itertools.combinations(self._sources, 2)
            
        else:
            
            return (
                filter(
                    lambda p:
                        p[0] != p[1],
                    itertools.product(self._sources, repeat = 2)
                )
            )
    
    def __iter__(self):
        
        for pair in self.iterpairs():
            
            for path in self.paths_between_pair(pair):
                
                yield Path(self.graph, path)
    
    def paths_between_pair(self, pair):
    
        shortest_path_len = self.graph.shortest_paths(pair[0], pair[0])[0][0]
        
        if (
            self.maxlen is not None and
            self.all_shortest and
            self.maxlen < shortest_path_len
        ):
            
            # if the shortest is longer than maxlen,
            # we still find something with all_shortest
            _maxlen = shortest_path_len
        
        else:
            
            _maxlen = self.maxlen
        
        # this way if shortest == maxlen, we still benefit from the
        # speed of the igraph method
        if _maxlen and (not self.all_shortest or _maxlen > shortest_path_len):
            
            return self.bfs_all_paths(pair, maxlen = _maxlen, silent = True)
        
        # otherwise we include only all shortest
        if not maxlen or all_shortest:
            
            for path in self.graph.get_all_shortest_paths(pair[0], pair[1],
                                                          mode = self.mode):
                
                yield path
        
    def bfs_all_paths(self, pair, maxlen = 2):
        """
        Finds all paths up to length `maxlen` between a pair of
        vertices. This function is needed only becaues igraph`s
        get_all_shortest_paths() finds only the shortest, not any
        path up to a defined length.
        
        :param tuple pair: Tuple of vertex indices.
            Indices of the source and target nodes of the paths.
        :param int maxlen:
            Maximum length of paths in steps, i.e. if maxlen = 3, then
            the longest path may consist of 3 edges and 4 nodes.
        """

        def bfs(start, end, path0, maxlen = None):
            
            path0 = path0 + [start]
            
            if start == end:
                yield path0
            
            if len(path0) <= maxlen:
                
                for node in self.adjlist[start] - set(path0):
                    
                    for path1 in bfs(node, end, path0, maxlen):
                        
                        yield path1
        
        if not hasattr(self, 'adjlist'):
            
            self.update_adjlist()
        
        for path in bfs(pair[0], pair[1], [], maxlen):
            
            yield path
    
    def update_adjlist(self, mode = 'OUT'):
        
        self.adjlist = [
            set(self.graph.neighbors(vi, mode = self.mode))
            for vi in xrange(self.graph.vcount())
        ]
    

class PathwayFinder(PathFinder):
    
    def __init__(self, pa, sources, targets = None, name_type = 'uniprot', maxlen = 2):
        """
        
        """

#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# (c) Denes Turei EMBL 2017
# turei.denes@gmail.com
#

import pypath
import re
import imp
import itertools
import collections

fnDrugInfo = 'drug_info.csv'
fnOut = 'targets_paths_{0}_{1}{2}.tsv'
restars = re.compile(r'([A-Za-z0-9]{2,})([\*]?)([A-Za-z0-9]*)')
old_omnipath = True
maxlen = None
all_shortest = True

def read_drug_info(fn):
    
    result = []
    
    with open(fn, 'r') as fp:
        
        _ = fp.readline()
        
        for li in fp:
            
            li = li.strip().replace('"', '').split(',')
            
            result.append([li[0], process_names(li[1])])
    
    return result

def process_names(names):
    
    return (
        list(
            map(
                lambda rn:
                    (
                        ''.join(rn) if not rn[1]
                            else
                        re.compile(
                            '%s%s%s' % (
                                rn[0],
                                (
                                    '[A-Z]'
                                    if rn[0][-1].isdigit() else
                                    '[0-9]+'
                                ),
                                rn[2]
                            )
                        )
                    ),
                map(
                    lambda n:
                        restars.match(n.strip()).groups(),
                    names.split(';')
                )
            )
        )
    )

def map_targets(llDrugsRaw, mapper):
    
    return (
        list(
            map(
                lambda dt:
                    (
                        dt[0], # compound
                        set(
                            itertools.chain(
                                *map(
                                    lambda t: (
                                        mapper.map_name(t, 'genesymbol', 'uniprot')
                                        if type(t) is str else
                                        itertools.chain(
                                            *map(
                                                lambda gs:
                                                    mapper.map_name(gs, 'genesymbol', 'uniprot'),
                                                filter(
                                                    lambda gs:
                                                        t.match(gs),
                                                    mapper.tables[9606][
                                                        ('genesymbol', 'uniprot')
                                                    ].mapping['to'].keys()
                                                )
                                            )
                                        )
                                    ),
                                    dt[1]
                                )
                            )
                        ), # set of UniProts
                        dt[1] # target original name
                    ),
                llDrugsRaw
            )
        )
    )

def path_uniprots(path, pa):
    
    return list(map(lambda vid: pa.nodNam[vid], path))

def path_genesymbols(path, pa):
    
    return list(map(lambda vid: pa.nodLab[vid], path))

def path_to_edges(path, pa):
    
    return (
        list(
            map(
                lambda i:
                    pa.graph.get_eid(path[i], path[i+1]),
                range(len(path) - 1)
            )
        )
    )

def signed_paths(path, pa):
    
    return (
        list(
            itertools.product(
                *map(
                    lambda e:
                        pa.graph.es[e]['dirs'].consensus_edges(),
                    path_to_edges(path, pa)
                )
            )
        )
    )

def path_undirected(spath):
    
    return (
        any(
            map(
                lambda s:
                    s[2] != 'directed',
                spath
            )
        )
    )

def path_serial(spath):
    
    return (
        all(
            map(
                lambda i:
                    spath[i][1] == spath[i+1][0],
                range(len(spath) - 1)
            )
        ) or all(
            map(
                lambda i:
                    spath[i][0] == spath[i+1][1],
                range(len(spath) - 1)
            )
        )
    )

def path_downstream(spath, upath):
    
    return (
        path_bidirectional(spath, upath, '__gt__') and
        not path_bidirectional(spath, upath, '__lt__')
    )

def path_upstream(spath, upath):
    
    return (
        path_bidirectional(spath, upath, '__lt__') and
        not path_bidirectional(spath, upath, '__gt__')
    )

def pairwise(iterable):
    
    a, b = itertools.tee(iterable)
    next(b, None)
    
    return zip(a, b)

def path_bidirectional(spath, upath, op):
    
    return (
        any(
            map(
                lambda x:
                    getattr(x[0], op)(x[1]),
                pairwise(path_reldirs(spath, upath))
            )
        )
    )

def path_reldirs(spath, upath):
    
    return (
        list(
            map(
                lambda s:
                    ((s[1][0], s[1][1]) ==
                     (upath[s[0]], upath[s[0]+1])),
                enumerate(spath)
            )
        )
    )

def section_sign(reldirs, spath, i):
    
    return (
        sum_signs(
            map(
                lambda ir:
                    spath[ir[0]][3],
                filter(
                    lambda ir:
                        ir[1] == reldirs[i],
                    enumerate(reldirs)
                )
            )
        )
    )

def dirpath_sign(spath, upath):
    
    return sum_signs(map(lambda s: s[3], spath))

def path_bidirectional_signs(spath, upath):
    
    reldirs = path_reldirs(spath, upath)
    
    return (
        section_sign(reldirs, spath,  0), # A ---> o
        section_sign(reldirs, spath, -1) # o <--- B
    )

def path_signs(spath):
    
    return (
        list(
            map(
                lambda s:
                    (
                        1  if s[3] == 'positive' else
                        -1 if s[3] == 'negative' else
                        0
                    ),
                spath
            )
        )
    )

def sum_signs(signs):
    
    csigns = collections.Counter(signs)
    
    return (
        0 if (
            not csigns or
            'unknown' in csigns
        ) else  1 if (
            'negative' not in csigns or
            csigns['negative'] % 2 == 0
        ) else -1
    )

def path_str(spath, upath, gsymbol = False, pa = None):
    
    reldirs = path_reldirs(spath, upath)
    signs   = path_signs(spath)
    
    return (
        '%s%s' % (
            ''.join(
                map(
                    lambda s:
                        step_str(
                            upath[s[0]],
                            s[1][2] == 'directed',
                            not reldirs[s[0]],
                            signs[s[0]],
                            gsymbol = gsymbol,
                            pa = pa
                        ),
                    enumerate(spath)
                )
            ),
            (
                pa.nodLab[pa.nodDct[upath[-1]]]
                if gsymbol else
                upath[-1]
            )
        )
    )

def interactions_str(spath):
    
    return (
        ';'.join(
            map(
                lambda s:
                    ' '.join([
                        s[0],
                        '%s' % (
                            'x' if s[2] == 'undirected' else
                            '+' if s[3] == 'positive' else
                            '-' if s[3] == 'negative' else
                            '*'),
                        s[1]
                    ]),
                spath
            )
        )
    )

def step_data(name1, name2, d, r, s, pa):
    
    names = [name1, name2]
    
    if r:
        names = reversed(names)
    
    return [
        names[0],
        names[1],
        pa.nodLab[pa.nodDct[names[0]]],
        pa.nodLab[pa.nodDct[names[1]]],
        int(d),
        int(s)
    ]

def step_str(name, d, r, s, gsymbol = False, pa = None):
    
    return (
        '%s %s--%s--%s ' % (
            (
                pa.nodLab[pa.nodDct[name]]
                if gsymbol else
                name
            ),
            '<' if r else '-',
            (
                '[-]' if s == -1 else
                '[+]' if s ==  1 else
                '---'
            ),
            '>' if d and not r else '-'
        )
    )

def is_reverse(spath, upath):
    
    return len(spath) and (spath[0][0], spath[0][1]) == (upath[1], upath[0])

def path_data(spath, upath, pa):
    
    typ = (
            'serial' if path_serial(spath) else
            'downstream' if path_downstream(spath, upath) else
            'upstream' if path_upstream(spath, upath) else
            'mixed'
        )
    
    direction = (
        0  if typ != 'serial' else
        -1 if is_reverse(spath, upath) else
        1
    )
    
    directed = int(not path_undirected(spath))
    
    signs = (
        dirpath_sign(spath, upath),
        0
    ) if typ == 'serial' else (
        path_bidirectional_signs(spath, upath)
    )
    
    if not directed:
        typ = 'undirected'
    
    return [
        upath[0],
        upath[-1],
        pa.nodLab[pa.nodDct[upath[0]]],
        pa.nodLab[pa.nodDct[upath[-1]]],
        len(spath),
        typ,
        direction,
        directed,
        signs[0],
        signs[1]
    ]

def path_line(spath, upath, pa):
    
    return path_data(spath, upath, pa) + [path_str(spath, upath),
                                          path_str(spath, upath, True, pa),
                                          interactions_str(spath)]

def get_path_header():
    
    return [
        'uniprot_a',
        'uniprot_b',
        'genesymbol_a',
        'genesymbol_b',
        'length',
        'type',
        'direction',
        'is_directed',
        'sign1',
        'sign2',
        'path_uniprots',
        'path_genesymbols',
        'interactions'
    ]

def path_lines(spaths, upath, pa):
    
    return (
        list(
            map(
                lambda p: (
                    path_line(p, upath, pa)
                ),
                spaths
            )
        )
    )

def process_paths(p, pa):
    
    return signed_paths(p, pa), path_uniprots(p, pa)

def full_lines(pathd, pa, prg = None):
    
    if prg is not None:
        prg.step()
    
    result = []
    
    for path in pathd['paths']:
        
        spaths, upath = process_paths(path, pa)
        
        result.extend(
            list(
                map(
                    lambda pl:
                        [
                            pathd['drug1'],
                            pathd['drug2'],
                            pathd['target_name1'],
                            pathd['target_name2'],
                        ] + pl,
                    path_lines(spaths, upath, pa)
                )
            )
        )
    
    return result

def all_full_lines(ldPaths, pa):
    
    prg = pypath.progress.Progress(len(ldPaths),
                                   'Processing paths', 1,
                                   percent = False)
    
    result = (
        list(
            itertools.chain(
                map(
                    lambda pathd:
                        full_lines(pathd, pa, prg),
                    ldPaths
                )
            )
        )
    )
    
    prg.terminate()
    
    return result

def get_full_header():
    
    return [
        'drug_a',
        'drug_b',
        'target_name_b',
        'target_name_b'
    ] + get_path_header()

def list2str(l, sep = '\t'):
    
    return (
        sep.join(
            map(
                lambda f:
                    (
                        list2str(f, sep = ';')
                        if type(f) is list else
                        f.pattern
                        if hasattr(f, 'pattern') else
                        '{0}'.format(f)
                    ),
                l
            )
        )
    )

def write_table(llPaths, fnOut, sep = '\t'):
    
    hdr = get_full_header()
    
    with open(fnOut, 'w') as fp:
        
        fp.write('%s\n' % (list2str(hdr)))
        
        fp.write(
            '\n'.join(
                map(
                    lambda l:
                        list2str(l),
                    llPaths
                )
            )
        )

def one_drug_pair_paths(dt1, dt2, pa,
                        maxlen = 2,
                        all_shortest = True,
                        prg = None):
    
    if prg is not None:
        prg.step()
    
    return {
        'drug1': dt1[0],
        'drug2': dt2[0],
        'paths': list(
            itertools.chain(
                *map(
                    lambda tt:
                        paths_between_pair(tt[0], tt[1], pa,
                                           maxlen, all_shortest),
                    itertools.product(dt1[1], dt2[1])
                )
            )
        ),
        'target_name1': dt1[2],
        'target_name2': dt2[2]
    }

def paths_between_pair(u1, u2, pa, maxlen = 2, all_shortest = True):
    
    result = []
    
    v1 = pa.uniprot(u1)
    v2 = pa.uniprot(u2)
    
    if v1 and v2:
        
        shortest_path_len = pa.graph.shortest_paths(v1.index, v2.index)[0][0]
        
        if (
            maxlen is not None and
            all_shortest and
            maxlen < shortest_path_len
           ):
            
            # if the shortest is longer than maxlen,
            # we still find something with all_shortest
            maxlen = shortest_path_len
        
        # this way if shortest == maxlen, we still benefit from the
        # speed of the igraph method
        if maxlen and (not all_shortest or maxlen > shortest_path_len):
            
            result.extend(pa.find_all_paths(v1.index, v2.index,
                                            maxlen = maxlen,
                                            silent = True,
                                            mode = 'ALL'))
        
        # otherwise we include only all shortest
        if not maxlen or all_shortest:
            
            return pa.graph.get_all_shortest_paths(
                    v1.index, v2.index, mode = 'ALL')
    
    return result

def all_drug_pairs(llDrugs, pa, maxlen = 2, all_shortest = True):
    
    ndrugs = len(llDrugs)
    prg = pypath.progress.Progress(ndrugs * (ndrugs - 1) / 2,
            'Looking up paths between all drug target pairs',
            1, percent = False)
    
    result = (
        list(
            itertools.chain(
                *map(
                    lambda i:
                        list(
                            map(
                                lambda j:
                                    one_drug_pair_paths(llDrugs[i],
                                                        llDrugs[j],
                                                        pa, maxlen,
                                                        all_shortest,
                                                        prg),
                                range(i + 1, ndrugs)
                            )
                        ),
                    range(ndrugs)
                )
            )
        )
    )
    
    prg.terminate()
    
    return result

def reload_pypath(pa):
    pa.reload()
    imp.reload(pypath.common)
    imp.reload(pypath)

def preload(infile, old_omnipath = False):
    
    pa = pypath.PyPath()
    
    if old_omnipath:
        pa.load_old_omnipath(kinase_substrate_extra = True)
    else:
        pa.load_omnipath()
    
    llDrugsRaw = read_drug_info(fnDrugInfo)
    llDrugs = map_targets(llDrugsRaw, pa.mapper)
    
    return llDrugs, pa

def main(llDrugs, pa, fnOut, maxlen = 2, sep = '\t',
         return_result = True, old_omnipath = False,
         all_shortest = True):
    
    ldPaths = all_drug_pairs(llDrugs, pa, maxlen, all_shortest)
    lllPaths = all_full_lines(ldPaths, pa)
    outfile = fnOut.format(
        'omnipath0' if old_omnipath else 'omnipath1',
        '%u' % maxlen if maxlen else '0',
        '_all_shortest' if all_shortest else ''
    )
    write_table(itertools.chain(*lllPaths), outfile, sep = sep)
    
    if return_result:
        
        return llPaths

if __name__ == '__main__':
    
    llDrugs, pa = preload(fnDrugInfo,
                          old_omnipath = old_omnipath)
    main(llDrugs, pa, fnOut,
         return_result = False,
         old_omnipath = old_omnipath,
         maxlen = maxlen,
         all_shortest = all_shortest)
