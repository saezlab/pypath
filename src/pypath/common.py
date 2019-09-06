#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Contains helper functions shared by different modules.
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

#TODO requires cleaning, check what functions are not used and may be removed.
#Some parts can go to jsons.

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import os
import sys
import re
import math
import random
import textwrap
import hashlib

import numpy as np

__all__ = [
    'ROOT', 'aacodes', 'aaletters', 'simpleTypes', 'numTypes',
    'uniqList', 'addToList', 'addToSet', 'gen_session_id',
    'sorensen_index', 'simpson_index', 'simpson_index_counts',
    'jaccard_index', 'console', 'wcl', 'flatList', 'charTypes',
    'delEmpty', 'get_args', 'something', 'rotate', 'cleanDict',
    'igraph_graphics_attrs', 'md5', 'mod_keywords', 'Namespace', 'fun',
    'taxids', 'taxa', 'phosphoelm_taxids', 'dbptm_taxids',
    'uniqOrdList', 'dict_diff', 'to_set', 'to_list',
    'unique_list', 'basestring', 'amino_acids', 'aminoa_1_to_3_letter',
    'aminoa_3_to_1_letter', 'pmod_bel', 'pmod_other_to_bel',
    'pmod_bel_to_other', 'refloat', 'reint', 'is_float', 'is_int',
    'float_or_nan',
]

# get the location
ROOT = os.path.abspath(os.path.dirname(__file__))
DATA = os.path.join(ROOT, 'data')


try:
    basestring
except NameError:
    basestring = str


class _const:

    class ConstError(TypeError):

        pass

    def __setattr__(self, name, value):

        if name in self.__dict__:

            raise(self.ConstError, "Can't rebind const(%s)" % name)

        self.__dict__[name] = value


default_name_type = {"protein": "uniprot",
                     "mirna": "mirbase",
                     "lncrna": "lncrnaname",
                     "drug": "pubchem"}

aacodes = {'A': 'ALA', 'G': 'GLY', 'M': 'MET', 'W': 'TRP', 'Y': 'TYR',
           'R': 'ARG', 'C': 'CYS', 'S': 'SER', 'V': 'VAL', 'I': 'ILE',
           'L': 'LEU', 'P': 'PRO', 'H': 'HIS', 'T': 'THR', 'E': 'GLU',
           'Q': 'GLN', 'Z': 'GLX', 'D': 'ASP', 'N': 'ASN', 'B': 'ASX',
           'K': 'LYS', 'F': 'PHE', 'J': 'XLE', 'X': 'XAA'}

aanames = {'alanine': 'A',
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
           'valine': 'V'}

mod_keywords = {'Reactome': [('phosphopantetheinylation', ['phosphopantet']),
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
                             ('sumoylation', ['sumoyl']),
                             ('prenylation', ['prenyl'])],
                'ACSN': [('phosphorylation', ['phospho']),
                         ('glycylation', ['glycyl']),
                         ('ubiquitination', ['ubiquitin']),
                         ('acetylation', ['acetyl']),
                         ('myristoylation', ['myristoyl']),
                         ('prenylation', ['prenyl']),
                         ('hydroxylation', ['hydroxy'])],
                'WikiPathways': [],
                'NetPath': [('phosphorylation', ['phospho']),
                            ('glycylation', ['glycyl']),
                            ('ubiquitination', ['ubiquitin']),
                            ('acetylation', ['acetyl'])],
                'PANTHER': [('phosphorylation', ['phospho']),
                            ('acetylation', ['acetyl'])],
                'NCI-PID': [('phosphorylation', ['phospho']),
                            ('methylation', ['methyl']),
                            ('farnesylation', ['farnesyl']),
                            ('palmitoylation', ['palmito']),
                            ('myristoylation', ['myristoyl']),
                            ('glycylation', ['glycyl']),
                            ('ubiquitination', ['ubiquitin']),
                            ('acetylation', ['acetyl']),
                            ('glycosylation', ['glycosyl']),
                            ('geranylation', ['geranyl']),
                            ('hydroxylation', ['hydroxy'])],
                'KEGG': [('phosphorylation', ['hospho']),
                         ('methylation', ['methyl']),
                         ('ubiquitination', ['ubiquitin']),
                         ('acetylation', ['acetyl']),
                         ('hydroxylation', ['hydroxy']),
                         ('carboxyethylation', ['carboxyethyl']),
                         ('ribosylation', ['ribosyl']),
                         ('nitrosylation', ['nitrosyl']),
                         ('sulfoylation', ['ulfo']),
                         ('biotinylation', ['biotinyl']),
                         ('malonylation', ['malonyl']),
                         ('glutarylation', ['lutaryl'])]}


# XXX: This is also in main.py, could be directly added in __all__?'
# For compatibility with python 2, see https://docs.python.org/3/whatsnew/3.0.html
if 'long' not in dir(__builtins__):
    long = int

if 'unicode' not in dir(__builtins__):
    unicode = str

# Inverted `aacodes` dict
aaletters = dict(zip(aacodes.values(), aacodes.keys()))

# Type definitions
simpleTypes = (int, long, float, str, unicode, bytes, bool, type(None))
numTypes = (int, long, float)
charTypes = (str, unicode, bytes)


refloat = re.compile(r'\s*-?\s*[\s\.\d]+\s*')
reint   = re.compile(r'\s*-?\s*[\s\d]+\s*')


def is_float(num):
    """
    Tells if a string represents a floating point number,
    i.e. it can be converted by `float`.
    """
    
    return bool(refloat.match(num))


def is_int(num):
    """
    Tells if a string represents an integer,
    i.e. it can be converted by `int`.
    """
    
    return bool(reint.match(num))


def float_or_nan(num):
    """
    Returns `num` converted from string to float if `num` represents a
    float otherwise `numpy.nan`.
    """
    
    return float(num) if is_float(num) else np.nan


def to_set(var):
    """
    Makes sure the object `var` is a set, if it is a list converts it to set,
    otherwise it creates a single element set out of it.
    If `var` is None returns empty set.
    """
    
    if isinstance(var, set):
        
        return var
        
    elif var is None:
        
        return set()
        
    elif isinstance(var, list):
        
        return set(var)
        
    else:
        
        return set([var])


def to_list(var):
    """
    Makes sure `var` is a list otherwise creates a single element list
    out of it. If `var` is None returns empty list.
    """
    
    if isinstance(var, list):
        
        return var
        
    elif var is None:
        
        return []
        
    else:
        
        return [var]


# From http://www.peterbe.com/plog/uniqifiers-benchmark
def unique_list(seq):
    """Reduces a list to its unique elements.

    Takes any iterable and returns a list of unique elements on it. If
    the argument is a dictionary, returns a list of unique keys.
    **NOTE:** Does not preserve the order of the elements.

    :arg list seq:
        Sequence to be processed, can be any iterable type.

    :return:
        (*list*) -- List of unique elements in the sequence *seq*.

    **Examples:**
        >>> uniqList('aba')
        ['a', 'b']
        >>> uniqList([0, 1, 2, 1, 0])
        [0, 1, 2]
    """

    return list({}.fromkeys(seq).keys())


# old ugly name
uniqList = unique_list


def uniqList1(seq): # XXX: Not used
    """
    Not order preserving
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    """

    return list(set(seq))


def uniqList2(seq): # XXX: Not used
    """
    Not order preserving
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    """

    keys = {}

    for e in seq:
        keys[e] = 1

    return list(keys.keys())


def flatList(lst):
    """Coerces the elements of a list of iterables into a single list.

    :arg lsit lst:
        List to be flattened. Its elements can be also lists or any
        other iterable.

    :return:
        (*list*) -- Flattened list of *lst*.

    **Examples:**
        >>> flatList([(0, 1), (1, 1), (2, 1)])
        [0, 1, 1, 1, 2, 1]
        >> flatList(['abc', 'def'])
        ['a', 'b', 'c', 'd', 'e', 'f']
    """

    return [it for sl in lst for it in sl]


def delEmpty(lst): # XXX: Only used in main.py line: 1278
    """Removes empty entries of a list.

    It is assumed that elemenst of *lst* are iterables (e.g. [str] or
    [list]).

    :arg list lst:
        List from which empty elements will be removed.

    :return:
        (*list*) -- Copy of *lst* without elements whose length was
        zero.

    **Example:**
        >>> delEmpty(['a', '', 'b', 'c'])
        ['a', 'b', 'c']
    """

    return [i for i in lst if len(i) > 0]


# Order preserving
# From http://www.peterbe.com/plog/uniqifiers-benchmark
def uniqOrdList(seq, idfun=None): # XXX: Only used in plot.py line: 510
    """Reduces a list to its unique elements keeping their order.

    Returns a copy of *seq* without repeated elements. Preserves the
    order. If any element is repeated, the first instance is kept.

    :arg list seq:
        Or any other iterable type. The sequence from which repeated
        elements are to be removed.
    :arg function idfun:
        Optional, ``None`` by default. Identifier function, for each
        entry of *seq*, returns a identifier of that entry from which
        uniqueness is determined. Default behavior is f(x) = x. See
        examples below.

    :return:
        (*list*) -- Copy of *seq* without the repeated elements
        (according to *idfun*).

    **Examples:**
        >>> uniqOrdList([0, 1, 2, 1, 5])
        [0, 1, 2, 5]
        >>> uniqOrdList('abracadabra')
        ['a', 'b', 'r', 'c', 'd']
        >>> def f(x):
        ...     if x > 0:
        ...             return 0
        ...     else:
        ...             return 1
        >>> uniqOrdList([-32, -42, 1, 15, -12], idfun=f)
        [-32, 1]
        >>> def g(x): # Given a file name, return it without extension
        ...    return x.split('.')[0]
        >>> uniqOrdList(['a.png', 'a.txt', 'b.pdf'], idfun=g)
        ['a.png', 'b.pdf']
    """

    if idfun is None:

        def idfun(x):
            return x

    seen = {}
    result = []

    for item in seq:
        marker = idfun(item)

        if marker in seen:
            continue

        seen[marker] = 1
        result.append(item)

    return result


def addToList(lst, toadd):
    """ Adds elements to a list.

    Appends *toadd* to *lst*. Function differs from
    :py:func:`list.append` since is capable to handle different data
    types. This is, if *lst* is not a list, it will be converted to one.
    Similarly, if *toadd* is not a list, it will be converted and added.
    If *toadd* is or contains ``None``, these will be ommited. The
    returned list will only contain unique elements and does not
    necessarily preserve order.

    :arg list lst:
        List or any other type (will be converted into a list). Original
        sequence to which *toadd* will be appended.
    :arg any toadd:
        Element(s) to be added into *lst*.

    :return:
        (*list*) -- Contains the unique element(s) from the union of
        *lst* and *toadd*. **NOTE:** Makes use of
        :py:func:`common.uniqList`, does not preserve order of elements.

    **Examples:**
        >>> addToList('ab', 'cd')
        ['ab', 'cd']
        >>> addToList('ab', ['cd', None, 'ab', 'ef'])
        ['ab', 'ef', 'cd']
        >>> addToList((0, 1, 2), 4)
        [0, 1, 2, 4]
    """

    if type(lst) is not list:

        if type(lst) in simpleTypes:
            lst = [lst]

        else:
            lst = list(lst)

    if toadd is None:
        return lst

    if type(toadd) in simpleTypes:
        lst.append(toadd)

    else:

        if type(toadd) is not list:
            toadd = list(toadd)

        lst.extend(toadd)

    if None in lst:
        lst.remove(None)

    return uniqList(lst)


def addToSet(st, toadd):
    """Adds elements to a set.

    Appends *toadd* to *st*. Function is capable to handle different
    input data types. This is, if *toadd* is a list, it will be
    converted to a set and added.

    :arg set st:
        Original set to which *toadd* will be added.
    :arg any toadd:
        Element(s) to be added into *st*.

    :return:
        (*set*) -- Contains the element(s) from the union of *st* and
        *toadd*.

    **Examples:**
        >>> st = set([0, 1, 2])
        >>> addToSet(st, 3)
        set([0, 1, 2, 3])
        >>> addToSet(st, [4, 2, 5])
        set([0, 1, 2, 4, 5])
    """

    if type(toadd) in simpleTypes:
        st.add(toadd)

    if type(toadd) is list:
        toadd = set(toadd)

    if type(toadd) is set:
        st.update(toadd)

    return st


def something(anything):
    """Checks if argument is empty.

    Checks if *anything* is empty or ``None``.

    :arg any anything:
        Self-explanatory.

    :return:
        (*bool*) -- ``False`` if *anyhting* is ``None`` or any empty
        data type.

    **Examples:**
        >>> something(None)
        False
        >>> something(123)
        True
        >>> something('Hello world!')
        True
        >>> something('')
        False
        >>> something([])
        False
    """

    return not (anything is None
                or (type(anything) in [list, set, dict, str, unicode]
                    and len(anything) == 0))


def gen_session_id(length=5):
    """Generates a random alphanumeric string.

    :arg int length:
        Optional, ``5`` by default. Specifies the length of the random
        string.

    :return:
        (*str*) -- Random alphanumeric string of the specified length.
    """

    abc = '0123456789abcdefghijklmnopqrstuvwxyz'

    return ''.join(random.choice(abc) for i in xrange(length))


# XXX: Are you sure this is the way to compute Simpson's index?
def simpson_index(a, b):
    """Computes Simpson's index.

    Given two sets *a* and *b*, returns the Simpson index.

    :arg set a:
        Or any iterable type (will be converted to set).
    :arg set b:
        Or any iterable type (will be converted to set).

    :return:
        (*float*) -- The Simpson index between *a* and *b*.
    """

    a = set(a)
    b = set(b)
    ab = a & b

    return float(len(ab)) / float(min(len(a), len(b)))


# XXX: Related to comment above, what is this exactly?
def simpson_index_counts(a, b, c):
    """
    :arg a:
    :arg b:
    :arg c:

    :return:
        (*float*) --
    """

    return float(c) / float(min(a, b)) if min(a, b) > 0 else 0.0


def sorensen_index(a, b):
    """Computes the Sorensen index.

    Computes the Sorensen-Dice coefficient between two sets *a* and *b*.

    :arg set a:
        Or any iterable type (will be converted to set).
    :arg set b:
        Or any iterable type (will be converted to set).

    :return:
        (*float*) -- The Sorensen-Dice coefficient between *a* and *b*.
    """

    a = set(a)
    b = set(b)
    ab = a & b

    return float(len(ab)) / float(len(a) + len(b))


def jaccard_index(a, b):
    """Computes the Jaccard index.

    Computes the Jaccard index between two sets *a* and *b*.

    :arg set a:
        Or any iterable type (will be converted to set).
    :arg set b:
        Or any iterable type (will be converted to set).

    :return:
        (*float*) -- The Jaccard index between *a* and *b*.
    """

    a = set(a)
    b = set(b)
    ab = a & b

    return float(len(ab)) / float(len(a | b))


def console(message):
    """Prints a message on the terminal.

    Prints a *message* to the standard output (e.g. terminal) formatted
    to 80 characters per line plus first-level indentation.

    :arg str message:
        The message to be printed.
    """

    message = '\n\t'.join(textwrap.wrap(message, 80))
    sys.stdout.write(('\n\t' + message).ljust(80))
    sys.stdout.write('\n')
    sys.stdout.flush()


def wcl(f): # XXX: Not used (another function w/ same name defined in curl.py)
    """

    """

    toClose = type(f) is file
    f = f if type(f) is file else open(f, 'r')

    for i, l in enumerate(f):
        pass

    if toClose:
        f.seek(0)

    else:
        f.close()

    return i + 1


# XXX: Not very clear to me the purpose of this function.
def get_args(loc_dict, remove=set([])):
    """
    Given a dictionary of local variables, returns a copy of it without
    ``'self'``, ``'kwargs'`` (in the scope of a :py:obj:`class`) plus
    any other specified in the keyword argument *remove*.

    :arg dict loc_dict:
        Dictionary containing the local variables (e.g. a call to
        :py:func:`locals` in a given scope).
    :arg set remove:
        Optional, ``set([])`` by default. Can also be a list. Contains
        the keys of the elements in *loc_dict* that will be removed.

    :return:
        (*dict*) -- Copy of *loc_dict* without ``'self'``, ``'kwargs'``
        and any other element specified in *remove*.
    """

    if type(remove) not in [set, list]:
        remove = set([remove])

    if type(remove) is list:
        remove = set(remove)

    remove.add('self')
    remove.add('kwargs')
    args = dict((k, v) for k, v in iteritems(loc_dict) if k not in remove)

    if 'kwargs' in loc_dict:
        args = dict(args.items() + loc_dict['kwargs'].items())

    return args


# From http://stackoverflow.com/a/20024348/854988
def rotate(point, angle, center=(0.0, 0.0)): # XXX: Not used? Wrote the docs before checking XD
    """
    Rotates a point with respect to a center.

    Rotates a given *point* around a *center* according to the specified
    *angle* (in degrees) in a two-dimensional space. The rotation is
    counter-clockwise.

    :arg tuple point:
        Or list. Contains the two coordinates of the point to be
        rotated.
    :arg float angle:
        Angle (in degrees) from which the point will be rotated with
        respect to *center* (counter-clockwise).
    :arg tuple center:
        Optional, ``(0.0, 0.0)`` by default. Can also be a list.
        Determines the two coordinates of the center relative to which
        the point has to be rotated.

    :return:
        (*tuple*) -- Pair of coordinates of the rotated point.
    """

    angle = math.radians(angle)
    temp_point = point[0] - center[0], point[1] - center[1]
    temp_point = (temp_point[0] * math.cos(angle)
                  - temp_point[1] * math.sin(angle),
                  temp_point[0] * math.sin(angle)
                  + temp_point[1] * math.cos(angle))
    temp_point = temp_point[0] + center[0], temp_point[1] + center[1]

    return temp_point


def cleanDict(dct):
    """Cleans a dictionary of ``None`` values.

    Removes ``None`` values from  a dictionary *dct* and casts all other
    values to strings.

    :arg dict dct:
        Dictionary to be cleaned from ``None`` values.

    :return:
        (*dict*) -- Copy of *dct* without ``None`` value entries and all
        other values formatted to strings.
    """

    toDel = []

    for k, v in iteritems(dct):

        if v is None:
            toDel.append(k)

        else:
            dct[k] = str(v)

    for k in toDel:
        del dct[k]

    return dct


def md5(value):
    """
    Computes the sum of MD5 hash of a given string *value*.

    :arg str value:
        Or any other type (will be converted to string). Value for which
        the MD5 sum will be computed. Must follow ASCII encoding.

    :return:
        (*str*) -- Hash value resulting from the MD5 sum of the *value*
        string.
    """

    try:
        string = str(value).encode('ascii')

    except: # XXX: Bad practice to catch all exceptions
        string = str(value).encode('ascii') # XXX: Same as the try statement!?

    return hashlib.md5(string).hexdigest()


# XXX: Shouldn't we keep all functions and variables separated
#      (together among them)?
igraph_graphics_attrs = {'vertex': ['size', ' color', 'frame_color',
                                    'frame_width', 'shape', 'label',
                                    'label_dist', 'label_color', 'label_size',
                                    'label_angle'],
                         'edge': ['curved', 'color', 'width', 'arrow_size',
                                  'arrow_width']}


def merge_dicts(d1, d2):
    """Merges dictionaries recursively.

    If a key exists in both dictionaries, the values will be merged.

    :arg dict d1:
        Base dictionary where *d2* will be merged.
    :arg dict d1:
        Dictionary to be merged into *d1*.

    :return:
        (*dict*) -- Resulting dictionary from the merging.
    """

    for k2, v2 in iteritems(d2):
        t = type(v2)

        if k2 not in d1:
            d1[k2] = v2

        elif t is dict:
            d1[k2] = merge_dicts(d1[k2], v2)

        elif t is list:
            d1[k2].extend(v2)

        elif t is set:
            d1[k2].update(v2)

    return d1


# XXX: Not 100% clear, correct if I'm mistaken. What do you mean by
#      simple heuristic?
def dict_set_path(d, path):
    """
    Given a dictionary of dictionaries *d* looks up the keys according
    to *path*, creates new subdicts and keys if those do not exist yet,
    and sets/merges the leaf element according to simple heuristic.

    :arg dict d:
        Dictionary of dictionaries for which the path is to be set.
    :arg list path:
        Or tuple, contains the path of keys being the first element a
        key of *d* (if doesn't exist will be created), and the
        subsequent of the inner dictionaries. The last element is the
        value that will be set/merged on the specified path.

    :return:
        (*dict*) -- Copy of *d* including the specified *path*.

    **Example:**
        >>> dict_set_path(dict(), ['a', 'b', 1])
        {'a': {'b': 1}}
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
            subd[key].update(val)

        else:
            subd[key].add(val)

    return d


def dict_diff(d1, d2):
    """Compares two dictionaries.

    Compares two given dictionaries *d1* and *d2* whose values are sets
    or dictionaries (in such case the function is called recursively).
    **NOTE:** The comparison is only performed on the values of the
    keys that are common in *d1* and *d2* (see example below).

    :arg dict d1:
        First dictionary of the comparison.
    :arg dict d2:
        Second dictionary of the comparison.

    :return:
        * (*dict*) -- Unique elements of *d1* when compared to *d2*.
        * (*dict*) -- Unique elements of *d2* when compared to *d1*.

    **Examples:**
        >>> d1 = {'a': {1}, 'b': {2}, 'c': {3}} # 'c' is unique to d1
        >>> d2 = {'a': {1}, 'b': {3}}
        >>> dict_diff(d1, d2)
        ({'a': set([]), 'b': set([2])}, {'a': set([2]), 'b': set([3])})
    """

    ldiff = {}
    rdiff = {}
    keys = set(d1.keys()) & set(d2.keys())

    for k in keys:

        if type(d1[k]) is dict and type(d2[k]) is dict:
            ldiff[k], rdiff[k] = dict_diff(d1[k], d2[k])

        elif type(d1[k]) is set and type(d2[k]) is set:
            ldiff[k], rdiff[k] = (d1[k] - d2[k], d2[k] - d1[k])

    return ldiff, rdiff


def dict_sym_diff(d1, d2): # XXX: Not used
    """

    """

    diff = {}
    keys = set(d1.keys()) & set(d2.keys())

    for k in keys:

        if type(d1[k]) is dict and type(d2[k]) is dict:
            diff[k] = dict_sym_diff(d1[k], d2[k])

        elif type(d1[k]) is set and type(d2[k]) is set:
            diff[k] = d1[k] ^ d2[k]

    return diff


def swap_dict(d):
    """Swaps a dictionary.

    Interchanges the keys and values of a dictionary. If the values are
    lists (or any iterable type) and/or not unique, each unique element
    will be a key and values sets of the original keys of *d* (see
    example below).

    :arg dict d:
        Original dictionary to be swapped.

    :return:
        (*dict*) -- The swapped dictionary.

    **Examples:**
        >>> d = {'a': 1, 'b': 2}
        >>> swap_dict(d)
        {1: 'a', 2: 'b'}
        >>> d = {'a': 1, 'b': 1, 'c': 2}
        >>> swap_dict(d)
        {1: set(['a', 'b']), 2: set(['c'])}
        d = {'a': [1, 2, 3], 'b': [2, 3]}
        >>> swap_dict(d)
        {1: set(['a']), 2: set(['a', 'b']), 3: set(['a', 'b'])}
    """

    _d = {}

    for key, vals in iteritems(d):
        
        vals = [vals] if type(vals) in simpleTypes else vals
        
        for val in vals:
            
            if val not in _d:
                _d[val] = set([])
            
            _d[val].add(key)

    if all(len(v) <= 1 for v in _d.values()):
        
        _d = dict((k, list(v)[0]) for k, v in iteritems(_d) if len(v))

    return _d


def swap_dict_simple(d):
    """
    Swaps a dictionary.

    Interchanges the keys and values of a dictionary. Assumes the values
    are unique and hashable, otherwise overwrites duplicates or raises
    error.

    :arg dict d:
        Original dictionary to be swapped.

    :return:
        (*dict*) -- The swapped dictionary.
    """

    return dict((v, k) for k, v in iteritems(d))


# XXX: Not sure what this joins exactly... I tried:
#      >>> a = {'a': [1], 'b': [2]}
#      >>> b = {'a': [2], 'b': [4]}
#      >>> join_dicts(a, b)
#      and got an empty dictionary (?)

def join_dicts(d1, d2, _from='keys', to='values'): # TODO
    """
    Joins a pair of dictionaries.

    :arg dict d1:
        Dictionary to be merged with *d2*
    :arg dict d2:
        Dictionary to be merged with *d1*
    :arg str _from:
        Optional, ``'keys'`` by default.
    :arg str to:
        Optional, ``'values'`` by default.

    :return:
        (*dict*) -- .
    """

    result = {}

    if to == 'keys':
        d2 = swap_dict(d2)

    for key1, val1 in iteritems(d1):
        sources = ([key1] if _from == 'keys'
                          else [val1] if type(val1) in simpleTypes
                                      else val1)
        meds = ([key1] if _from == 'values'
                       else [val1] if type(val1) in simpleTypes
                                   else val1)
        targets = set([])

        for med in meds:

            if med in d2:

                if type(targets) is list:
                    targets.append(d2[med])

                elif type(d2[med]) in simpleTypes:
                    targets.add(d2[med])

                elif type(d2[med]) is list:
                    targets.update(set(d2[med]))

                elif type(d2[med]) is set:
                    targets.update(d2[med])

                elif d2[med].__hash__ is not None:
                    targets.add(d2[med])

                else:
                    targets = list(targets)
                    targets.append(d2[med])

        for source in sources:

            if type(targets) is list:

                if source not in result:
                    result[source] = []

                result[source].extend(targets)

            elif type(targets) is set:

                if source not in result:
                    result[source] = set([])

                result[source].update(targets)

    if all(len(x) <= 1 for x in result.values()):
        result = dict((k, list(v)[0]) for k, v in iteritems(result) if len(v))

    return result


class Namespace(object): # XXX: WHY? + Not used anywhere
    pass


def fun(): # XXX: Best name for a function. Not 100% sure if used anywhere...
    """

    """

    print(__name__)
    print(__name__ in globals())

    for n in __name__.split('.'):
        print(n, n in globals())

    return __name__


# XXX: Shouldn't we keep all functions and variables separated
#      (together among them)?
taxids = {
    9606: 'human',
    10090: 'mouse',
    10116: 'rat',
    9031: 'chicken',
    9913: 'cow',
    9986: 'rabbit',
    9940: 'sheep',
    10141: 'guinea pig',
    10036: 'hamster',
    7227: 'fruit fly',
    9615: 'dog',
    9823: 'pig',
    8355: 'frog',
    9091: 'quail',
    9796: 'horse',
    9925: 'goat',
    89462: 'water buffalo',
    9598: 'monkey',
    9103: 'turkey',
    9685: 'cat',
    7604: 'starfish',
    7609: 'spiny starfish',
    1213717: 'torpedo',
    9669: 'ferret',
    8839: 'duck',
    9593: 'gorilla',
    7460: 'honeybee',
    8407: 'european common frog',
    9544: 'rhesus macaque',
}

taxa = swap_dict(taxids)

taxa_synonyms = {
    'bovine': 'cow',
}

phosphoelm_taxids = {
    9606: 'Homo sapiens',
    10090: 'Mus musculus',
    9913: 'Bos taurus',
    9986: 'Oryctolagus cuniculus',
    9615: 'Canis familiaris',
    10029: 'Cricetulus griseus',
    9267: 'Didelphis virginiana',
    9031: 'Gallus gallus',
    10036: 'Mesocricetus auratus',
    9940: 'Ovis aries',
    10116: 'Rattus norvegicus',
    9823: 'Sus scrofa',
    8355: 'Xenopus laevis',
    10141: 'Cavia porcellus',
    9796: 'Equus caballus',
    7227: 'Drosophila melanogaster',
    487: 'Neisseria meningitidis',
    562: 'Escherichia coli',
    5207: 'Cryptococcus neoformans',
    470: 'Acinetobacter baumannii',
    1280: 'Staphylococcus aureus',
    4939: 'Saccharomyces cerevisiae',
    34: 'Myxococcus xanthus',
    1392: 'Bacillus anthracis',
    210: 'Helicobacter pylori',
    6239: 'Caenorhabditis elegans',
}

dbptm_taxids = {
    9606: 'HUMAN',
    10090: 'MOUSE',
    7227: 'DROME',
    10116: 'RAT',
    559292: 'YEAST',
    284812: 'SCHPO',
    4081: 'SOLLC',
    3702: 'ARATH',
    9940: 'SHEEP',
    9913: 'BOVIN',
    9925: 'CAPHI',
    44689: 'DICDI',
    4577: 'MAIZE',
    9823: 'PIG',
    9615: 'CANLF',
    6239: 'CAEEL',
    8455: 'XENLA',
    83333: 'ECOLI',
    1891767: 'SV40',
}


dbptm_to_ncbi_tax_id = swap_dict(dbptm_taxids)
latin_name_to_ncbi_tax_id = swap_dict(phosphoelm_taxids)

def taxid_from_common_name(taxon_name):
    
    taxon_name = taxon_name.lower().strip()
    
    if not taxon_name or taxon_name in {'none', 'unknown'}:
        
        return None
    
    if taxon_name in taxa_synonyms:
        
        taxon_name = taxa_synonyms[taxon_name]
    
    if taxon_name in taxa:
        
        return taxa[taxon_name]


def taxid_from_latin_name(taxon_name):
    
    if taxon_name in latin_name_to_ncbi_tax_id:
        
        return latin_name_to_ncbi_tax_id[taxon_name]


def taxid_from_dbptm_taxon_name(taxon_name):
    
    if taxon_name in dbptm_to_ncbi_tax_id:
        
        return dbptm_to_ncbi_tax_id[taxon_name]


def ensure_ncbi_tax_id(taxon_id):
    
    if isinstance(taxon_id, int):
        
        return taxon_id
        
    else:
        
        return (
            taxid_from_common_name(taxon_id) or
            taxid_from_latin_name(taxon_id) or
            taxid_from_dbptm_taxon_name(taxon_id)
        )




psite_mod_types = [('p', 'phosphorylation'),
                   ('ac', 'acetylation'),
                   ('ga', 'galactosylation'),
                   ('gl', 'glycosylation'),
                   ('sm', 'sumoylation'),
                   ('ub', 'ubiquitination'),
                   ('me', 'methylation')]

psite_mod_types2 = [('p', 'phosphorylation'),
                    ('ac', 'acetylation'),
                    ('ga', 'galactosylation'),
                    ('gl', 'glycosylation'),
                    ('ub', 'ubiquitination'),
                    ('me', 'methylation'),
                    ('sm', 'sumoylation'),
                    ('sc', 'succinylation'),
                    ('m1', 'mono-methylation'),
                    ('m2', 'di-methylation'),
                    ('m3', 'tri-methylation'),
                    ('ad', 'adenylation'),
                    ('pa', 'palmitoylation'),
                    ('ne', 'neddylation'),
                    ('sn', 'nitrosylation'),
                    ('ca', 'caspase-cleavage')]


pmod_bel = (
    ('Ac', ('acetylation',)),
    ('ADPRib', ('ADP-ribosylation', 'ADP-rybosylation', 'adenosine diphosphoribosyl',)),
    ('Farn', ('farnesylation',)),
    ('Gerger', ('geranylgeranylation',)),
    ('Glyco', ('glycosylation',)),
    ('Hy', ('hydroxylation',)),
    ('ISG', ('ISGylation', 'ISG15-protein conjugation',)),
    ('Me', ('methylation',)),
    ('Me1', ('monomethylation', 'mono-methylation',)),
    ('Me2', ('dimethylation', 'di-methylation',)),
    ('Me3', ('trimethylation', 'tri-methylation',)),
    ('Myr', ('myristoylation',)),
    ('Nedd', ('neddylation',)),
    ('NGlyco', ('N-linked glycosylation',)),
    ('OGlyco', ('O-linked glycosylation',)),
    ('Palm', ('palmitoylation',)),
    ('Ph', ('phosphorylation',)),
    ('Sumo', ('sumoylation',)),
    ('Ub', ('ubiquitination', 'ubiquitinylation', 'ubiquitylation',)),
    ('UbK48', ('Lysine 48-linked polyubiquitination',)),
    ('UbK63', ('Lysine 63-linked polyubiquitination',)),
    ('UbMono', ('monoubiquitination',)),
    ('UbPoly', ('polyubiquitination',)),
)

pmod_bel_to_other = dict(pmod_bel)
pmod_other_to_bel = dict(
    (other_name, bel_name)
    for bel_name, other_names in pmod_bel
    for other_name in other_names
)


mirbase_taxids = {9606: 'hsa',
                  10090: 'mmu',
                  10116: 'rno',
                  7227: 'dme'}

amino_acids = (
    ('alanine', 'Ala', 'A'),
    ('arginine', 'Arg', 'R'),
    ('asparagine', 'Asn', 'N'),
    ('aspartic acid', 'Asp', 'D'),
    ('asparagine or aspartic acid', 'Asx', 'B'),
    ('cysteine', 'Cys', 'C'),
    ('glutamic acid', 'Glu', 'E'),
    ('glutamine', 'Gln', 'Q'),
    ('glutamine or glutamic acid', 'Glx', 'Z'),
    ('glycine', 'Gly', 'G'),
    ('histidine', 'His', 'H'),
    ('isoleucine', 'Ile', 'I'),
    ('leucine', 'Leu', 'L'),
    ('lysine', 'Lys', 'K'),
    ('methionine', 'Met', 'M'),
    ('phenylalanine', 'Phe', 'F'),
    ('proline', 'Pro', 'P'),
    ('serine', 'Ser', 'S'),
    ('threonine', 'Thr', 'T'),
    ('tryptophan', 'Trp', 'W'),
    ('tyrosine', 'Tyr', 'Y'),
    ('valine', 'Val', 'V'),
)


aminoa_3_to_1_letter = dict(
    (code3, code1)
    for name, code3, code1 in amino_acids
)


aminoa_1_to_3_letter = dict(
    (code1, code3)
    for name, code3, code1 in amino_acids
)


class silent(object): # XXX: Never used
    """

    """

    def __init__(self):

        pass

    def __enter__(self):
        self.aux = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exception_type, exception_value, traceback):

        sys.stdout.close()
        sys.stdout = self.aux


def paginate(lst, size = 10):
    """
    Yields sections of length ``size`` from list ``lst``.
    The last section might be shorter than ``size``.
    Following https://stackoverflow.com/a/3744502/854988.
    """
    
    for i in xrange((len(lst) // size) + 1):
        
        yield lst[size * i:size * (i + 1)]
