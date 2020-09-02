#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Contains helper functions shared by different modules.
#
#  Copyright
#  2014-2020
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
import copy
import random
import operator
import collections
import itertools
import warnings
import textwrap
import hashlib
import inspect
import textwrap
import tabulate

import numpy as np


__all__ = [
    'ROOT',
    'aacodes', 'aaletters',
    'simple_types', 'numeric_types', 'list_like',
    'uniq_list', 'add_to_list', 'add_to_set', 'gen_session_id',
    'sorensen_index', 'simpson_index', 'simpson_index_counts',
    'jaccard_index', 'console', 'wcl', 'flat_list', 'char_types',
    'del_empty', 'get_args', 'something', 'rotate', 'clean_dict',
    'igraph_graphics_attrs', 'md5', 'mod_keywords',
    'uniq_ord_list', 'dict_diff', 'to_set', 'to_list',
    'unique_list', 'basestring', 'amino_acids', 'aminoa_1_to_3_letter',
    'aminoa_3_to_1_letter', 'pmod_bel', 'pmod_other_to_bel',
    'pmod_bel_to_other', 'refloat', 'reint', 'is_float', 'is_int',
    'float_or_nan', 'non_digit',
]

# get the location
ROOT = os.path.join(
    *os.path.split(
        os.path.abspath(os.path.dirname(__file__))
    )[:-1]
)
DATA = os.path.join(ROOT, 'data')


try:
    basestring
except NameError:
    basestring = str

non_digit = re.compile(r'[^\d.-]+')


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
simple_types = (int, long, float, str, unicode, bytes, bool, type(None))
numeric_types = (int, long, float)
char_types = (str, unicode, bytes)
list_like = (tuple, set, list)


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

    elif not isinstance(var, basestring) and hasattr(var, '__iter__'):

        return set(var)

    else:

        return {var}


def to_list(var):
    """
    Makes sure `var` is a list otherwise creates a single element list
    out of it. If `var` is None returns empty list.
    """

    if isinstance(var, list):

        return var

    elif var is None:

        return []

    elif not isinstance(var, basestring) and hasattr(var, '__iter__'):

        return list(var)

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
        >>> uniq_list('aba')
        ['a', 'b']
        >>> uniq_list([0, 1, 2, 1, 0])
        [0, 1, 2]
    """

    return list({}.fromkeys(seq).keys())


# old ugly name
uniq_list = unique_list


def uniq_list1(seq): # XXX: Not used
    """
    Not order preserving
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    """

    return list(set(seq))


def uniq_list2(seq): # XXX: Not used
    """
    Not order preserving
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    """

    keys = {}

    for e in seq:
        keys[e] = 1

    return list(keys.keys())


def flat_list(lst):
    """Coerces the elements of a list of iterables into a single list.

    :arg lsit lst:
        List to be flattened. Its elements can be also lists or any
        other iterable.

    :return:
        (*list*) -- Flattened list of *lst*.

    **Examples:**
        >>> flat_list([(0, 1), (1, 1), (2, 1)])
        [0, 1, 1, 1, 2, 1]
        >> flat_list(['abc', 'def'])
        ['a', 'b', 'c', 'd', 'e', 'f']
    """

    return [it for sl in lst for it in sl]


def del_empty(lst): # XXX: Only used in main.py line: 1278
    """Removes empty entries of a list.

    It is assumed that elemenst of *lst* are iterables (e.g. [str] or
    [list]).

    :arg list lst:
        List from which empty elements will be removed.

    :return:
        (*list*) -- Copy of *lst* without elements whose length was
        zero.

    **Example:**
        >>> del_empty(['a', '', 'b', 'c'])
        ['a', 'b', 'c']
    """

    return [i for i in lst if i or isinstance(i, (int, float))]


# Order preserving
# From http://www.peterbe.com/plog/uniqifiers-benchmark
def uniq_ord_list(seq, idfun=None): # XXX: Only used in plot.py line: 510
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
        >>> uniq_ord_list([0, 1, 2, 1, 5])
        [0, 1, 2, 5]
        >>> uniq_ord_list('abracadabra')
        ['a', 'b', 'r', 'c', 'd']
        >>> def f(x):
        ...     if x > 0:
        ...             return 0
        ...     else:
        ...             return 1
        >>> uniq_ord_list([-32, -42, 1, 15, -12], idfun=f)
        [-32, 1]
        >>> def g(x): # Given a file name, return it without extension
        ...    return x.split('.')[0]
        >>> uniq_ord_list(['a.png', 'a.txt', 'b.pdf'], idfun=g)
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


def add_to_list(lst, toadd):
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
        :py:func:`common.uniq_list`, does not preserve order of elements.

    **Examples:**
        >>> add_to_list('ab', 'cd')
        ['ab', 'cd']
        >>> add_to_list('ab', ['cd', None, 'ab', 'ef'])
        ['ab', 'ef', 'cd']
        >>> add_to_list((0, 1, 2), 4)
        [0, 1, 2, 4]
    """

    if type(lst) is not list:

        if type(lst) in simple_types:
            lst = [lst]

        else:
            lst = list(lst)

    if toadd is None:
        return lst

    if type(toadd) in simple_types:
        lst.append(toadd)

    else:

        if type(toadd) is not list:
            toadd = list(toadd)

        lst.extend(toadd)

    if None in lst:
        lst.remove(None)

    return uniq_list(lst)


def add_to_set(st, toadd):
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
        >>> add_to_set(st, 3)
        set([0, 1, 2, 3])
        >>> add_to_set(st, [4, 2, 5])
        set([0, 1, 2, 4, 5])
    """

    if type(toadd) in simple_types:
        st.add(toadd)

    if type(toadd) is list:
        toadd = set(toadd)

    if type(toadd) is set:
        st.update(toadd)

    return st


def upper0(string):
    """
    Ensures the first letter of a string is uppercase, except if the first
    word already contains uppercase letters, in order to avoid changing,
    for example, miRNA to MiRNA.
    """

    if not string:

        return string

    else:

        words = string.split(' ', maxsplit = 1)

        if words[0] and words[0].lower() == words[0]:

            words[0] = words[0][0].upper() + words[0][1:]

        return ' '.join(words)


def first(it, default = None):
    """
    Returns the first element of the iterable ``it`` or the value of
    ``default`` if the iterable is empty.
    """

    for i in it:

        return i

    return default


def swap_suffix(name, sep = '_', suffixes = None):
    """
    Changes the suffix of a string.

    suffixes : dict
        A mapping for the swap, by default is `{'a': 'b', 'b': 'a'}`.
    """

    suffixes = suffixes or {'a': 'b', 'b': 'a'}

    name_suffix = name.rsplit(sep, maxsplit = 1)

    if len(name_suffix) == 2 and name_suffix[1] in suffixes:

        name = '%s%s%s' % (name_suffix[0], sep, suffixes[name_suffix[1]])

    return name



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


def clean_dict(dct):
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

        vals = [vals] if type(vals) in simple_types else vals

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
                          else [val1] if type(val1) in simple_types
                                      else val1)
        meds = ([key1] if _from == 'values'
                       else [val1] if type(val1) in simple_types
                                   else val1)
        targets = set([])

        for med in meds:

            if med in d2:

                if type(targets) is list:
                    targets.append(d2[med])

                elif type(d2[med]) in simple_types:
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


def shared_unique(by_group, group, op = 'shared'):
    """
    For a *dict* of *set*s ``by_group`` and a particular key ``group``
    returns a *set* of all elements in the *set* belonging to the
    key ``group`` which either does or does not occure in any of the sets
    assigned to the other keys, depending on the operator ``op``.
    This method can be used among other things to find the shared and
    unique entities across resources.

    :arg str op:
        Either `shared` or `unique`.
    """

    if group not in by_group:

        warnings.warn(
            'Group `%s` missing from the dict of groups!' % group
        )

    _op = operator.sub if op == 'unique' else operator.and_

    return _op(
        by_group[group] if group in by_group else set(),
        set.union(*(
            elements
            for label, elements
            in iteritems(by_group)
            if label != group
        ), set())
    )


def shared_elements(by_group, group):
    """
    For a *dict* of *set*s ``by_group`` and a particular key ``group``
    returns a *set* of all elements in the *set* belonging to the
    key ``group`` which occure in any of the sets assigned to the other keys.
    This method can be used among other things to find the shared entities
    across resources.
    """

    return shared_unique(
        by_group = by_group,
        group = group,
        op = 'shared',
    )


def unique_elements(by_group, group):
    """
    For a *dict* of *set*s ``by_group`` and a particular key ``group``
    returns a *set* of all elements in the *set* belonging to the
    key ``group`` which don't occure in any of the sets assigned to the
    other keys. This method can be used among other things to find the
    unique entities across resources.
    """

    return shared_unique(
        by_group = by_group,
        group = group,
        op = 'unique',
    )


def n_shared_elements(by_group, group):
    """
    For a *dict* of *set*s ``by_group`` and a particular key ``group``
    returns the number of all elements in the *set* belonging to the
    key ``group`` which occure in any of the other sets.
    This method can be used among other things to count the shared entities
    across resources.
    """

    return len(shared_elements(by_group = by_group, group = group))


def n_unique_elements(by_group, group):
    """
    For a *dict* of *set*s ``by_group`` and a particular key ``group``
    returns the number of all elements in the *set* belonging to the
    key ``group`` which don't occure in any of the other sets.
    This method can be used among other things to count the unique entities
    across resources.
    """

    return len(unique_elements(by_group = by_group, group = group))


def shared_unique_foreach(by_group, op = 'shared', counts = False):
    """
    For a *dict* of *set*s ``by_group`` returns a *dict* of *set*s with
    either shared or unique elements across all *set*s, depending on
    the operation ``op``.
    """

    method = len if counts else lambda x: x

    return dict(
        (
            label,
            method(
                shared_unique(by_group = by_group, group = label, op = op)
            ),
        )
        for label in by_group.keys()
    )


def n_shared_unique_foreach(by_group, op = 'shared'):
    """
    For a *dict* of *set*s ``by_group`` returns a *dict* of numbers with
    the counts of either the shared or unique elements across all *set*s,
    depending on the operation ``op``.
    """

    return shared_unique_foreach(
        by_group = by_group,
        op = 'shared',
        counts = True,
    )


def shared_foreach(by_group):

    return shared_unique_foreach(by_group = by_group, op = 'shared')


def unique_foreach(by_group):

    return shared_unique_foreach(by_group = by_group, op = 'unique')


def n_shared_foreach(by_group):

    return n_shared_unique_foreach(by_group = by_group, op = 'shared')


def n_unique_foreach(by_group):

    return n_shared_unique_foreach(by_group = by_group, op = 'unique')


def dict_union(dict_of_sets):
    """
    For a *dict* of *set*s returns the union of the values.
    """

    return set.union(*dict_of_sets.values()) if dict_of_sets else set()


def dict_counts(dict_of_sets):
    """
    For a *dict* of *set*s or other values with ``__len__`` returns a
    *dict* of numbers with the length of each value in the original *dict*.

    This function is recursively works on dicts of dicts.
    """

    return dict(
        (
            key,
            (
                dict_counts(val)
                    if isinstance(val, dict) else
                len(val)
            )
        )
        for key, val in iteritems(dict_of_sets)
    )


def dict_expand_keys(dct, depth = 1, front = True):
    """
    From a *dict* with *tuple* keys builds a dict of dicts.

    :arg dict dct:
        A *dict* with tuple keys (if keys are not tuples ``dct`` will be
        returned unchanged).
    :arg int depth:
        Expand the keys up to this depth. If 0 *dct* will be returned
        unchanged, if 1 dict of dicts, if 2 dict of dict of dicts will be
        returned, and so on.
    :arg bool front:
        If ``True`` the tuple keys will be chopped from the front, otherwise
        from their ends.
    """

    if depth == 0:

        return dct

    new = {}

    for key, val in iteritems(dct):

        if not isinstance(key, tuple):

            new[key] = val

        elif len(key) == 1:

            new[key[0]] = val

        else:

            outer_key = key[0] if front else key[:-1]
            inner_key = key[1:] if front else key[-1]

            if len(inner_key) == 1:

                inner_key = inner_key[0]

            sub_dct = new.setdefault(outer_key, {})
            sub_dct[inner_key] = val

    if depth > 1:

        new = (
            dict(
                (
                    key,
                    dict_expand_keys(sub_dct, depth = depth - 1)
                )
                for key, sub_dct in iteritems(new)
            )
                if front else
            dict_expand_keys(new, depth = depth - 1, front = False)
        )

    return new


def dict_collapse_keys(
        dct,
        depth = 1,
        front = True,
        expand_tuple_keys = True,
    ):
    """
    From a dict of dicts builds a dict with tuple keys.

    :arg dict dct:
        A dict of dicts (if values are not dicts it will be returned
        unchanged).
    :arg int depth:
        Collapse the keys up to this depth. If 0 *dct* will be returned
        unchanged, if 1 tuple keys will have 2 elements, if 2 then
        2 elements, and so on.
    :arg bool front:
        If ``True`` the tuple keys will be collapsed first from the
        outermost dict going towards the innermost one until depth allows.
        Otherwise the method will start from the innermost ones.
    :arg bool expand_tuple_keys:
        If ``True`` the tuple keys of inner dicts will be concatenated with
        the outer key tuples. If ``False`` the inner tuple keys will be added
        as an element of the tuple key i.e. tuple in tuple.
    """

    if not front:

        # this is difficult to implement because we have no idea about
        # the depth; this version ensures an even key length for the
        # tuple keys; another alterntive would be to iterate recursively
        # over the dictionary tree
        dct = dict_collapse_keys(dct, depth = 9999999)
        maxdepth = max(
            len(k for k in dct.keys() if isinstance(k, tuple)),
            default = 0
        )
        return dict_expand_keys(dct, depth = maxdepth - depth, front = True)

    if not any(isinstance(val, dict) for val in dct.values()):

        return dct

    new = {}

    for key, val in iteritems(dct):

        key = key if isinstance(key, tuple) else (key,)

        if isinstance(val, dict):

            for key1, val1 in iteritems(val):

                _key = key + (
                    key1
                        if (
                            isinstance(key1, tuple) and
                            expand_tuple_keys
                        ) else
                    (key1,)
                )
                new[_key] = val1

        else:

            new[key] = val

    if depth > 1:

        new = dict_collapse_keys(new, depth = depth - 1)

    return new


def shared_unique_total(by_group, op = 'shared'):

    counts = collections.Counter(itertools.chain(*by_group.values()))
    _op = operator.eq if op == 'unique' else operator.gt

    return {key for key, val in iteritems(counts) if _op(val, 1)}


def shared_total(by_group):

    return shared_unique_total(by_group = by_group, op = 'shared')


def unique_total(by_group):

    return shared_unique_total(by_group = by_group, op = 'unique')


def n_shared_total(by_group):

    return len(shared_total(by_group))


def n_unique_total(by_group):

    return len(unique_total(by_group))


def dict_subtotals(dct):
    """
    For a dict of dicts of sets returns a dict with keys of the outer dict
    and values the union of the sets in each of the inner dicts.
    """

    return dict(
        (
            key,
            dict_union(sub_dct)
        )
        for key, sub_dct in iteritems(dct)
    )


def dict_percent(dict_of_counts, total):
    """
    For a *dict* of counts and a total count creates a *dict* of percentages.
    """

    return dict(
        (key, (val / total if total != 0 else 0) * 100)
        for key, val in iteritems(dict_of_counts)
    )


def dict_set_percent(dict_of_sets):

    total = len(dict_union(dict_of_sets))
    counts = dict_counts(dict_of_sets)

    return dict_percent(counts, total)


def df_memory_usage(df, deep = True):
    """
    Returns the memory usage of a ``pandas.DataFrame`` as a string.
    Modified from ``pandas.DataFrame.info``.
    """

    dtypes = {str(dt) for dt in df.dtypes}

    size_qualifier = (
        '+'
            if (
                'object' in dtypes or
                df.index._is_memory_usage_qualified()
            ) else
        ''
    )

    mem_usage = df.memory_usage(index = True, deep = deep).sum()

    for unit in ['bytes', 'KB', 'MB', 'GB', 'TB']:

        if mem_usage < 1024.0:

            return '%3.1f%s %s' % (mem_usage, size_qualifier, unit)

        mem_usage /= 1024.0

    return '%3.1f%s PB' % (mem_usage, size_qualifier)


def sum_dicts(*args):
    """
    For dicts of numbers returns a dict with the sum of the numbers from
    all dicts for all keys.
    """

    args = [
        collections.defaultdict(int, d)
        for d in args
    ]

    return dict(
        (
            key,
            sum(d[key] for d in args)
        )
        for key in set(itertools.chain(*(d.keys() for d in args)))
    )


def combine_attrs(attrs, num_method = max):
    """
    Combines multiple attributes into one. This method attempts
    to find out which is the best way to combine attributes.

        * If there is only one value or one of them is None, then
          returns the one available.
        * Lists: concatenates unique values of lists.
        * Numbers: returns the greater by default or calls
          *num_method* if given.
        * Sets: returns the union.
        * Dictionaries: calls :py:func:`pypath.common.merge_dicts`.
        * Direction: calls their special
          :py:meth:`pypath.main.Direction.merge` method.

    Works on more than 2 attributes recursively.

    :arg list attrs:
        List of one or more attribute values.
    :arg function num_method:
        Optional, ``max`` by default. Method to merge numeric attributes.
    """

    def list_or_set(one, two):

        if ((isinstance(one, list) and isinstance(two, set))
            or (isinstance(two, list) and isinstance(one, set))):

            try:
                return set(one), set(two)

            except TypeError:
                return list(one), list(two)

        else:
            return one, two


    # recursion:
    if len(attrs) > 2:
        attrs = [attrs[0], combine_attrs(attrs[1:], num_method = num_method)]

    # quick and simple cases:
    if len(attrs) == 0:
        return None

    if len(attrs) == 1:
        return attrs[0]

    if attrs[0] == attrs[1]:
        return attrs[0]

    if attrs[0] is None:
        return attrs[1]

    if attrs[1] is None:
        return attrs[0]

    # merge numeric values
    if type(attrs[0]) in numeric_types and type(attrs[1]) in numeric_types:
        return num_method(attrs)

    attrs = list(attrs)

    # in case one is list other is set
    attrs[0], attrs[1] = list_or_set(attrs[0], attrs[1])

    # merge lists:
    if isinstance(attrs[0], list) and isinstance(attrs[1], list):

        try:
            # lists of hashable elements only:
            return list(set(itertools.chain(attrs[0], attrs[1])))

        except TypeError:
            # if contain non-hashable elements:
            return list(itertools.chain(attrs[0], attrs[1]))

    # merge sets:
    if isinstance(attrs[0], set):
        return add_to_set(attrs[0], attrs[1])

    if isinstance(attrs[1], set):
        return add_to_set(attrs[1], attrs[0])

    # merge dicts:
    if isinstance(attrs[0], dict) and isinstance(attrs[1], dict):
        return merge_dicts(attrs[0], attrs[1])

    # 2 different strings: return a set with both of them
    if ((isinstance(attrs[0], str) or isinstance(attrs[0], unicode))
        and (isinstance(attrs[1], str) or isinstance(attrs[1], unicode))):

        if len(attrs[0]) == 0:
            return attrs[1]

        if len(attrs[1]) == 0:
            return attrs[0]

        return set([attrs[0], attrs[1]])

    # one attr is list, the other is simple value:
    if (isinstance(attrs[0], list) and type(attrs[1]) in simple_types):

        if attrs[1] in numeric_types or len(attrs[1]) > 0:
            return add_to_list(attrs[0], attrs[1])

        else:
            return attrs[0]

    if (isinstance(attrs[1], list) and type(attrs[0]) in simple_types):

        if attrs[0] in numeric_types or len(attrs[0]) > 0:
            return add_to_list(attrs[1], attrs[0])

        else:
            return attrs[1]

    # in case the objects have `__add__()` method:
    if hasattr(attrs[0], '__add__'):

        return attrs[0] + attrs[1]


def _add_method(cls, method_name, method, signature = None, doc = None):

    method.__name__ = method_name

    if signature and hasattr(inspect, 'Signature'): # Py2

        if not isinstance(signature, inspect.Signature):

            signature = inspect.Signature([
                inspect.Parameter(
                    name = param[0],
                    kind = inspect.Parameter.POSITIONAL_OR_KEYWORD,
                    default = (
                        param[1]
                            if len(param) > 1 else
                        inspect.Parameter.empty
                    )
                )
                for param in signature
            ])

        method.__signature__ = signature

    if doc:

        method.__doc__ = doc

    setattr(cls, method_name, method)


def at_least_in(n = 2):
    """
    Returns a method which is similar to the `intersection` operator on sets
    but requires the elements to present at least in ``n`` of the sets
    instead of in all of them. If ``n = 1`` it is equivalent with ``union``,
    if ``n`` is the same as the number of the sets it is equivalent with
    ``intersection``. The returned method accepts an arbitrary number of
    sets as non-keyword arguments.
    """

    def _at_least_in(*args):

        if len(args) < n:

            return set()

        counter = collections.Counter(itertools.chain(*args))

        return {
            key
            for key, count in iteritems(counter)
            if count >= n
        }

    return _at_least_in


def eqs(one, other):
    """
    Equality between ``one`` and ``other``. If any of them is type of `set`,
    returns True if it contains the other. If both of them are `set`,
    returns True if they share any element. Lists, tuples and similar objects
    will be converted to `set`.
    """

    one = one if isinstance(one, simple_types) else to_set(one)
    other = other if isinstance(other, simple_types) else to_set(other)

    if isinstance(one, set):

        if isinstance(other, set):

            return bool(one & other)

        else:

            return other in one

    elif isinstance(other, set):

        return one in other

    else:

        return one == other


def dict_str(dct):

    if not isinstance(dct, dict):

        return str(dct)

    return ', '.join(
        '%s=%s' % (str(key), str(val))
        for key, val in iteritems(dct)
    )


def none_or_len(value):

    return None if not hasattr(value, '__len__') else len(value)


def sets_to_sorted_lists(obj):

    if isinstance(obj, dict):

        return dict(
            (k, sets_to_sorted_list(v))
            for k, v in iteritems(obj)
        )

    elif isinstance(obj, (list, set, tuple)):

        return sorted(obj)

    else:

        return obj


def wrap_truncate(text, width = None, maxlen = None):

    if isinstance(text, list_like):

        text = ', '.join(text)

    if not isinstance(text, basestring):

        text = str(text)

    if maxlen:

        text = textwrap.shorten(text, width = maxlen)

    if width:

        text = textwrap.wrap(text, width = width)

    return os.linesep.join(text) if isinstance(text, list_like) else text


def table_add_row_numbers(tbl, **kwargs):

    nrows = len(next(tbl.values().__iter__())) if len(tbl) else 0

    return collections.OrderedDict(
        itertools.chain(
            (('No.', list(xrange(1, nrows + 1))),),
            tbl.items(),
        )
    )


def table_textwrap(tbl, width = None, maxlen = None):
    """
    Wraps and truncates the text content of cells in a table.
    The table is an ``OrderedDict`` with column titles as keys and column
    contents as lists.
    """

    def get_width(i):

        return (
            width[i]
                if (
                    isinstance(width, list_like) or
                    isinstance(width, dict) and
                    i in width
                ) else
            width['default']
                if (
                    isinstance(width, dict) and
                    'default' in width
                ) else
            width
        )


    return collections.OrderedDict(
        (
            wrap_truncate(title, width = get_width(i), maxlen = maxlen),
            [
                wrap_truncate(cell, width = width, maxlen = maxlen)
                for cell in column
            ]
        )
        for i, (title, column) in enumerate(iteritems(tbl))
    )


tabulate_defaults = {
    'numalign': 'right',
}


def table_format(
        tbl,
        width = None,
        maxlen = None,
        tablefmt = 'fancy_grid',
        wrap = True,
        lineno = True,
        **kwargs
    ):

    if wrap:

        tbl = table_textwrap(tbl, width = width, maxlen = maxlen)

    if lineno:

        tbl = table_add_row_numbers(tbl)

    tabulate_param = copy.deepcopy(tabulate_defaults)
    tabulate_param.update(kwargs)


    return tabulate.tabulate(
        zip(*tbl.values()),
        tbl.keys(),
        tablefmt = tablefmt,
        **tabulate_param
    )


def print_table(
        tbl,
        width = None,
        maxlen = None,
        tablefmt = 'fancy_grid',
        wrap = True,
        lineno = True,
        **kwargs
    ):

    sys.stdout.write(
        table_format(
            tbl,
            width = width,
            maxlen = wifth,
            tablefmt = tablefmt,
            wrape = wrap,
            lineno = lineno,
            **kwargs
        )
    )
    sys.stdout.write(os.linesep)
    sys.stdout.flush()


def tsv_table(
        tbl,
        path = None,
        maxlen = None,
        **kwargs,
    ):
    """
    From a table represented by an OrderedDict with column titles as keys
    and column contents as lists generates a tab separated string.
    If ``path`` provided writes out the tsv into a file, otherwise returns
    the string.
    """

    tbl = table_textwrap(tbl, width = None, maxlen = maxlen)
    tsv = []
    tsv.append('\t'.join(tbl.keys()))
    tsv.extend([
        '\t'.join(map(str, row))
        for row in zip(*tbl.values())
    ])
    tsv = os.linesep.join(tsv)

    if path:

        with open(path, 'w') as fp:

            fp.write(tsv)

    else:

        return tsv


def latex_table(
        tbl,
        colformat = None,
        maxlen = None,
        lineno = True,
        path = None,
        doc_template = True,
        booktabs = True,
        latex_compile = False,
        latex_executable = 'xelatex',
        latex_engine = 'xelatex',
        **kwargs
    ):
    """
    From a table represented by an OrderedDict with column titles as keys
    and column contents as lists generates LaTeX tabular.
    If ``path`` provided writes out the table into a file,
    if ``latex_compile`` is True compiles the document, otherwise returns
    it as a string.
    """

    maxlen = maxlen or 999999

    _xelatex_header = [
        r'\usepackage[no-math]{fontspec}',
        r'\usepackage{xunicode}',
        r'\usepackage{polyglossia}',
        r'\setdefaultlanguage{english}',
        r'\usepackage{xltxtra}',
    ]

    _pdflatex_header = [
        r'\usepackage[utf8]{inputenc}'
        r'\usepackage[T1]{fontenc}'
        r'\usepackage[english]{babel}'
    ]

    _doc_template_default = [
        r'\documentclass[9pt, a4paper, landscape]{article}',
    ]
    _doc_template_default.extend(
        _xelatex_header if latex_engine == 'xelatex' else _pdflatex_header
    )
    _doc_template_default.extend([
        r'\usepackage{array}',
        r'\usepackage{tabularx}'
        r'\usepackage{xltabular}',
        r'\usepackage{booktabs}',
        r'\usepackage[table]{xcolor}',
        (
            r'\usepackage[landscape,top=1cm,bottom=2cm,left=1cm,right=1cm]'
            r'{geometry}'
        ),
        r'\newcolumntype{L}{>{\raggedright\arraybackslash}X}',
        r'\newcolumntype{K}[1]{>{\raggedright\arraybackslash}p{#1}}'
        r'\renewcommand{\arraystretch}{1.5}',
        r'\begin{document}',
        r'\fontsize{4pt}{5pt}\selectfont'
        r'\rowcolors{2}{gray!25}{white}'
        r'',
        r'%s',
        r'',
        r'\end{document}',
    ])

    doc_template = (
        doc_template
            if isinstance(doc_template, basestring) else
        os.linesep.join(_doc_template_default)
            if doc_template or latex_compile else
        '%s'
    )

    _ = kwargs.pop('wrap', None)

    kwargs['tablefmt'] = 'latex_%s' % ('booktabs' if booktabs else 'raw')

    tbl = table_textwrap(tbl, width = None, maxlen = maxlen)
    tbl = collections.OrderedDict(
        (
            upper0(title.replace('_', ' ')),
            column,
        )
        for title, column in iteritems(tbl)
    )

    latex_table = table_format(
        tbl = tbl,
        maxlen = maxlen,
        lineno = lineno,
        wrap = False,
        **kwargs,
    )

    latex_table = latex_table.replace('tabular', 'xltabular')
    latex_table = latex_table.replace(
        r'\begin{xltabular',
        r'\begin{xltabular}{\linewidth'
    )
    recolformat = re.compile(r'(xltabular\}\{\\linewidth\}\{)(\w+)(\})')
    if not colformat:
        m = recolformat.search(latex_table)
        colformat = m.groups()[1].rsplit('r', maxsplit = 1)
        colformat = '%sr%s' % (colformat[0], colformat[1].replace('l', 'L'))

    latex_table = recolformat.sub(r'\g<1>%s\g<3>' % colformat, latex_table)
    latex_table_head, latex_table_body = (
        latex_table.split(r'\midrule', maxsplit = 1)
    )
    latex_table_head = os.linesep.join((
        latex_table_head,
        r'\midrule',
        r'\endhead',
        '',
    ))
    latex_full = doc_template % (latex_table_head + latex_table_body)
    latex_full = latex_full.replace(r'\ensuremath{<}', r'\textless ')
    latex_full = latex_full.replace(r'\ensuremath{>}', r'\textgreater ')
    latex_full = latex_full.replace(r'\_', '-')

    if not path and latex_compile:

        path = 'table-%s' % gen_session_id()

    if path and os.path.splitext(path)[1] != '.tex':

        path = '%s.tex' % path

    if path:

        with open(path, 'w') as fp:

            fp.write(latex_full)

    if latex_compile and doc_template:

        # doing twice to make sure it compiles all right
        os.system('%s %s' % (latex_executable, path))
        os.system('%s %s' % (latex_executable, path))

    if not path:

        return latex_full
