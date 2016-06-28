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

import os
import sys
import math
import random
import textwrap
import hashlib

__all__ = ['ROOT', 'aacodes', 'aaletters', 'simpleTypes', 'numTypes', 'uniqList', 'addToList', 
           'gen_session_id', 'sorensen_index', 'console', 'wcl', 'flatList',
           'charTypes', 'delEmpty', '__version__', 'get_args',
           'something', 'rotate', 'cleanDict', 'igraph_graphics_attrs', 'md5', 'mod_keywords']

# get the location
ROOT = os.path.abspath(os.path.dirname(__file__))

def _get_version():
    with open(os.path.join(ROOT, '__version__'), 'r') as v:
        return tuple([int(i) for i in v.read().strip().split('.')])

_MAJOR, _MINOR, _MICRO = _get_version()
__version__ = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
__release__ = '%d.%d' % (_MAJOR, _MINOR)


default_name_type = {
    "protein": "uniprot", 
    "mirna": "mirbase", 
    "lncrna": "lncrnaname",
    "drug": "pubchem"
}

aacodes = {
    'A': 'ALA',
    'G': 'GLY',
    'M': 'MET',
    'W': 'TRP',
    'Y': 'TYR',
    'R': 'ARG',
    'C': 'CYS',
    'S': 'SER',
    'V': 'VAL',
    'I': 'ILE',
    'L': 'LEU',
    'P': 'PRO',
    'H': 'HIS',
    'T': 'THR',
    'E': 'GLU',
    'Q': 'GLN',
    'Z': 'GLX',
    'D': 'ASP',
    'N': 'ASN',
    'B': 'ASX',
    'K': 'LYS',
    'F': 'PHE',
    'J': 'XLE',
    'X': 'XAA'
}

aanames = {
    'alanine': 'A',
    'arginine': 'R',
    'asparagine': 'N',
    'aspartic acid': 'D',
    'cysteine': 'C',
    'glutamine': 'Q',
    'glutamic acid': 'E',
    'glycine': 'G',
    'histidine': 'H',
    'isoleucine': 'I',
    'leucine': 'L',
    'lysine': 'K',
    'methionine': 'M',
    'phenylalanine': 'F',
    'proline': 'P',
    'serine': 'S',
    'threonine': 'T',
    'tryptophan': 'W',
    'tyrosine': 'Y',
    'valine': 'V'
}

mod_keywords = {
    'Reactome': [
        ('phosphopantetheinylation', ['phosphopantet']),
        ('phosphorylation', ['phospho']),
        ('acetylneuraminylation', ['acetylneuraminyl']),
        ('acetylation', ['acetyl']),
        ('farnesylation', ['farnesyl']),
        ('palmitoylation', ['palmito']),
        ('methylation', ['methyl']),
        ('tetradecanoylation', ['tetradecanoyl']),
        ('decanoylation', ['decanoyl']),
        ('palmitoleylation', ['palmytoleil']),
        ('formylation', ['formyl']),
        ('ubiquitination', ['ubiquitin']),
        ('galactosylation', ['galactos']),
        ('glutamylation', ['glutamyl']),
        ('fucosylation', ['fucosyl']),
        ('myristoylation', ['myristoyl']),
        ('carboxylation', ['carboxyl']),
        ('biotinylation', ['biotinyl']),
        ('glycosylation', ['glycosyl']),
        ('octanoylation', ['octanoyl']),
        ('glycylation', ['glycyl']),
        ('hydroxylation', ['hydroxy']),
        ('sulfhydration', ['persulfid']),
        ('thiolation', ['thio']),
        ('amidation', ['amide']),
        ('selenation', ['seleno']),
        ('glucosylation', ['glucosyl']),
        ('neddylation', ['neddyl']),
        ('sumoylation', ['sumoyl'])
    ],
    'ACSN': [],
    'WikiPathways': []
}

if 'long' not in __builtins__:
    long = int

if 'unicode' not in __builtins__:
    unicode = str

aaletters = dict(zip(aacodes.values(),aacodes.keys()))

simpleTypes = set([int, long, float, str, unicode, bytes])

numTypes = set([int, long, float])

charTypes = set([str, unicode, bytes])

def uniqList(seq):
    # Not order preserving
    # from http://www.peterbe.com/plog/uniqifiers-benchmark
    keys = {}
    for e in seq:
        try:
            keys[e] = 1
        except:
            print('ERROR in pypath.common.uniqList(): unhashable type:')
            print(e)
    return list(keys.keys())

def flatList(lst):
    return [it for sl in lst for it in sl]

def delEmpty(lst):
    return [i for i in lst if len(i) > 0]

def uniqOrdList(seq, idfun = None): 
    # Order preserving
    # from http://www.peterbe.com/plog/uniqifiers-benchmark
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        try:
            seen[marker] = 1
        except:
            print('ERROR in pypath.common.uniqOrdList(): unhashable type:')
            print(marker)
        result.append(item)
    return result

def addToList(lst, toadd):
    if isinstance(toadd, list):
        lst += toadd
    else:
        lst.append(toadd)
    if None in lst:
        lst.remove(None)
    return uniqList(lst)

def something(anything):
    return not (anything is None or \
        (type(anything) in [list, set, dict, str, unicode] \
            and len(anything) == 0))

def gen_session_id(length=5):
    abc = '0123456789abcdefghijklmnopqrstuvwxyz'
    return ''.join(random.choice(abc) for i in range(length))

def simpson_index(a, b):
    a = set(a)
    b = set(b)
    ab = a & b
    return float(len(ab)) / float(min(len(a),len(b)))

def sorensen_index(a, b):
    a = set(a)
    b = set(b)
    ab = a & b
    return float(len(ab)) / float(len(a) + len(b))

def jaccard_index(a, b):
    a = set(a)
    b = set(b)
    ab = a & b
    return float(len(ab)) / float(len(a | b))

def console(message):
    message = '\n\t'.join(textwrap.wrap(message,80))
    sys.stdout.write(('\n\t'+message).ljust(80))
    sys.stdout.write('\n')
    sys.stdout.flush()

def wcl(f):
    toClose = type(f) is file
    f = f if type(f) is file else open(f, 'r')
    for i, l in enumerate(f):
        pass
    if toClose:
        f.seek(0)
    else:
        f.close()
    return i + 1

def get_args(loc_dict, remove = set([])):
    if type(remove) not in [set, list]: remove = set([remove])
    if type(remove) is list: remove = set(remove)
    remove.add('self')
    remove.add('kwargs')
    args = dict((k, v) for k, v in iteritems(loc_dict) if k not in remove)
    if 'kwargs' in loc_dict: args = dict(args.items() + loc_dict['kwargs'].items())
    return args

def rotate(point, angle, center = (0.0, 0.0)):
    """
    from http://stackoverflow.com/a/20024348/854988
    Rotates a point around center. Angle is in degrees.
    Rotation is counter-clockwise
    """
    angle = math.radians(angle)
    temp_point = point[0]-center[0] , point[1]-center[1]
    temp_point = ( temp_point[0]*math.cos(angle)-temp_point[1]*math.sin(angle), \
        temp_point[0]*math.sin(angle) + temp_point[1]*math.cos(angle))
    temp_point = temp_point[0]+center[0] , temp_point[1]+center[1]
    return temp_point

def cleanDict(dct):
    for k, v in dct.items():
        if v is None:
            del dct[k]
        else:
            dct[k] = str(v)
    return dct

def md5(value):
    try:
        string = str(value).encode('ascii')
    except:
        string = str(value).encode('ascii')
    return hashlib.md5(string).hexdigest()

igraph_graphics_attrs = {
    'vertex': ['size', ' color', 'frame_color', 'frame_width', 'shape', 'label', 'label_dist', 'label_color', 'label_size', 'label_angle'],
    'edge': ['curved', 'color', 'width', 'arrow_size', 'arrow_width']
}

def merge_dicts(d1, d2):
    """
    Merges dicts recursively
    """
    for k2, v2 in iteritems(d2):
        if k2 not in d1:
            d1[k2] = v2
        elif type(v2) is dict:
            d1[k2] = merge_dicts(d1[k2], v2)
        elif type(v2) is list:
            d1[k2].extend(v2)
        elif type(v2) is set:
            d1[k2] = d1[k2].union(v2)
    return d1

def dict_set_path(d, path):
    """
    In dict of dicts ``d`` looks up the keys following ``path``,
    creates new subdicts and keys if those do not exist yet,
    and sets/merges the leaf element according to simple heuristic.
    """
    val = path[-1]
    key = path[-2]
    subd = d
    for k in path[:-2]:
        if type(subd) is dict:
            if k in subd:
                subd = subd[k]
            else:
                subd[k] = {}
                subd = subd[k]
        else:
            return d
    if key not in subd:
        subd[key] = val
    elif type(val) is dict and type(subd[key]) is dict:
        subd[key].update(val)
    elif type(subd[key]) is list:
        if type(val) is list:
            subd[key].extend(val)
        else:
            subd[key].append(val)
    elif type(subd[key]) is set:
        if type(val) is set:
            subd[key] = subd.union(val)
        else:
            subd[key].add(val)
    return d
