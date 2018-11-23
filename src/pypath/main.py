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
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import __main__ # XXX: Bad practice

from future.utils import iteritems
from past.builtins import xrange, range, reduce

# external modules:
import os
import sys
import re
import imp
import math
from functools import reduce

try:
    import cairo

except ImportError: # XXX: Catching any exception like this is bad practice
    sys.stdout.write(
        '\t:: Module `cairo` not available.\n'
        '\t   Some plotting functionalities won\'t be accessible.\n'
    ) # XXX: Single/double quotes are not consistent throughout the code

import igraph
import codecs
import random # XXX: This and other Python built-in modules shoud be up
import textwrap
import copy
import json
import operator
import locale
import heapq
import threading
import traceback
import itertools
from itertools import chain
import collections
from collections import Counter
from scipy import stats
import numpy as np

# XXX: Put together all import tries

try:
    import cPickle as pickle

except ImportError:
    import pickle

try:
    import pygraphviz as graphviz

except ImportError:#ModuleNotFoundError:
    sys.stdout.write(
        '\t:: Module `pygraphviz` not available.\n'
        '\t   You don\'t need it unless you want to export dot files.\n')
    sys.stdout.flush()

try:
    import pandas

except ImportError:#ModuleNotFoundError:
    sys.stdout.write('\t:: Module `pandas` not available.\n')
    sys.stdout.flush()

# from this module:
import pypath.logn as logn
import pypath.data_formats as data_formats
import pypath.mapping as mapping
import pypath.descriptions as descriptions
import pypath.chembl as chembl
import pypath.mysql as mysql
import pypath.dataio as dataio
import pypath.uniprot_input as uniprot_input
import pypath.curl as curl
import pypath.intera as intera
import pypath.seq as se
import pypath.go as go
import pypath.gsea as gsea
import pypath.drawing as bdrawing
import pypath.proteomicsdb as proteomicsdb
import pypath.reflists as reflists
import pypath.input_formats as input_formats
import pypath.refs as _refs
import pypath.plot as plot
import pypath.ptm
import pypath.export as export
import pypath.ig_drawing as ig_drawing
import pypath.common as common
import pypath._version as _version
from pypath.gr_plot import *
from pypath.progress import *
import pypath.settings as settings

omnipath = data_formats.omnipath

# XXX: The following aliases are already defined in common.py
if 'long' not in __builtins__:
    long = int

if 'unicode' not in __builtins__:
    unicode = str

# XXX: Referenced but not defined: __version__, a (what?), ReferenceList
__all__ = ['PyPath', 'Direction', '__version__', 'a',
           'AttrHelper', 'ReferenceList', 'omnipath']


class Direction(object):
    """
    Object storing directionality information of an edge. Also includes
    information about the reverse direction, mode of regulation and
    sources of that information.

    :arg str nameA:
        Name of the source node.
    :arg str nameB:
        Name of the target node.

    :var dict dirs:
        Dictionary containing the presence of directionality of the
        given edge. Keys are *straight*, *reverse* and ``'undirected'``
        and their values denote the presence/absence [bool].
    :var dict negative:
        Dictionary contianing the presence/absence [bool] of negative
        interactions for both :py:attr:`straight` and :py:attr:`reverse`
        directions.
    :var dict negative_sources:
        Contains the resource names [str] supporting a negative
        interaction on :py:attr:`straight` and :py:attr:`reverse`
        directions.
    :var list nodes:
        Contains the node names [str] sorted alphabetically (*nameA*,
        *nameB*).
    :var dict positive:
        Dictionary contianing the presence/absence [bool] of positive
        interactions for both :py:attr:`straight` and :py:attr:`reverse`
        directions.
    :var dict positive_sources:
        Contains the resource names [str] supporting a positive
        interaction on :py:attr:`straight` and :py:attr:`reverse`
        directions.
    :var tuple reverse:
        Contains the node names [str] in reverse order e.g. (*nameB*,
        *nameA*).
    :var dict sources:
        Contains the resource names [str] of a given edge for each
        directionality (:py:attr:`straight`, :py:attr:`reverse` and
        ``'undirected'``). Values are sets containing the names of those
        resources supporting such directionality.
    :var tuple straight:
        Contains the node names [str] in the original order e.g.
        (*nameA*, *nameB*).
    """

    __slots__ = ['nodes', 'straight', 'reverse', 'dirs', 'sources', 'positive',
                 'negative', 'positive_sources', 'negative_sources']

    def __init__(self, nameA, nameB):
        """Initializes the edge object between the given nodes."""

        self.nodes = [nameA, nameB]
        self.nodes.sort()

        self.straight = (self.nodes[0], self.nodes[1])
        self.reverse = (self.nodes[1], self.nodes[0])

        self.dirs = {self.straight: False,
                     self.reverse: False,
                     'undirected': False}
        self.sources = {self.straight: set([]),
                        self.reverse: set([]),
                        'undirected': set([])}

        self.positive = {self.straight: False, self.reverse: False}
        self.negative = {self.straight: False, self.reverse: False}

        self.positive_sources = {self.straight: set([]), self.reverse: set([])}
        self.negative_sources = {self.straight: set([]), self.reverse: set([])}

    def __str__(self):
        """Custom string/printing function for the object."""

    # XXX: Pretty sure this could be refactored and implemented more
    #      efficiently like using str.format()

        s = 'Directions and signs of interaction between %s and %s\n\n' % \
            (self.nodes[0], self.nodes[1])

        if self.dirs[self.straight]:
            s += '\t%s ===> %s :: %s\n' % \
                (self.nodes[0], self.nodes[1], ', '.join(
                    self.sources[self.straight]))

        if self.dirs[self.reverse]:
            s += '\t%s <=== %s :: %s\n' % \
                (self.nodes[0], self.nodes[1],
                 ', '.join(self.sources[self.reverse]))

        if self.dirs['undirected']:
            s += '\t%s ==== %s :: %s\n' % \
                (self.nodes[0], self.nodes[1],
                 ', '.join(self.sources['undirected']))

        if self.positive[self.straight]:
            s += '\t%s =+=> %s :: %s\n' % (
                self.nodes[0], self.nodes[1],
                ', '.join(self.positive_sources[self.straight]))

        if self.positive[self.reverse]:
            s += '\t%s <=+= %s :: %s\n' % (
                self.nodes[0], self.nodes[1],
                ', '.join(self.positive_sources[self.reverse]))

        if self.negative[self.straight]:
            s += '\t%s =-=> %s :: %s\n' % (
                self.nodes[0], self.nodes[1],
                ', '.join(self.negative_sources[self.straight]))

        if self.negative[self.reverse]:
            s += '\t%s <=-= %s :: %s\n' % (
                self.nodes[0], self.nodes[1],
                ', '.join(self.negative_sources[self.reverse]))

        return s

    def reload(self):
        """Reloads the object from the module level."""

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def check_nodes(self, nodes):
        """Checks if *nodes* is contained in the edge.

        :arg list nodes:
             Or [tuple], contains the names of the nodes to be checked.

        :return:
            (*bool*) -- ``True`` if all elements in *nodes* are
            contained in the object :py:attr:`nodes` list.
        """

        return not bool(len(set(self.nodes) - set(nodes)))

    def check_param(self, di):
        """
        Checks if *di* is ``'undirected'`` or contains the nodes of
        the current edge. Used internally to check that *di* is a valid
        key for the object attributes declared on dictionaries.

        :arg tuple di:
            Or [str], key to be tested for validity.

        :return:
            (*bool*) -- ``True`` if *di* is ``'undirected'`` or a tuple
              of node names contained in the edge, ``False`` otherwise.
        """

        return (di == 'undirected' or (isinstance(di, tuple) and
                                       self.check_nodes(di)))

    def set_dir(self, direction, source):
        """
        Adds directionality information with the corresponding data
        source named. Modifies self attributes :py:attr:`dirs` and
        :py:attr:`sources`.

        :arg tuple direction:
            Or [str], the directionality key for which the value on
            :py:attr:`dirs` has to be set ``True``.
        :arg set source:
            Contains the name(s) of the source(s) from which such
            information was obtained.
        """

        if self.check_param(direction) and len(source):
            self.dirs[direction] = True
            source = common.addToSet(set([]), source)
            self.sources[direction] = self.sources[direction] | source

    def get_dir(self, direction, sources=False):
        """
        Returns the state (or *sources* if specified) of the given
        *direction*.

        :arg tuple direction:
            Or [str] (if ``'undirected'``). Pair of nodes from which
            direction information is to be retrieved.
        :arg bool sources:
            Optional, ``'False'`` by default. Specifies if the
            :py:attr:`sources` information of the given direction is to
            be retrieved instead.

        :return:
            (*bool* or *set*) -- (if ``sources=True``). Presence/absence
            of the requested direction (or the list of sources if
            specified). Returns ``None`` if *direction* is not valid.
        """

        if self.check_param(direction):

            if sources:
                return self.sources[direction]

            else:
                return self.dirs[direction]

        else:
            return None

    def get_dirs(self, src, tgt, sources=False):
        """
        Returns all directions with boolean values or list of sources.

        :arg str src:
            Source node.
        :arg str tgt:
            Target node.
        :arg bool sources:
            Optional, ``False`` by default. Specifies whether to return
            the :py:attr:`sources` attribute instead of :py:attr:`dirs`.

        :return:
            Contains the :py:attr:`dirs` (or :py:attr:`sources` if
            specified) of the given edge.
        """

    # XXX: What's the point of using src and tgt if in the end straigth,
    #      reverse and undirected are returned? Also, in such case:
    #      `return self.sources/dirs.values()` is more straightforward.

        query = (src, tgt)

        if self.check_nodes(query):

            if sources:
                return [self.sources[query],
                        self.sources[(query[1], query[0])],
                        self.sources['undirected']]

            else:
                return [self.dirs[query],
                        self.dirs[(query[1], query[0])],
                        self.dirs['undirected']]

        else:
            return None

    def which_dirs(self):
        """
        Returns the pair(s) of nodes for which there is information
        about their directionality.

        :return:
            (*list*) -- List of tuples containing the nodes for which
            their attribute :py:attr:`dirs` is ``True``.
        """

        return [d for d, s in iteritems(self.dirs) if s and d != 'undirected']

    def unset_dir(self, direction, source=None):
        """
        Removes directionality and/or source information of the
        specified *direction*. Modifies attribute :py:attr:`dirs` and
        :py:attr:`sources`.

        :arg tuple direction:
            Or [str] (if ``'undirected'``) the pair of nodes specifying
            the directionality from which the information is to be
            removed.
        :arg set source:
            Optional, ``None`` by default. If specified, determines
            which specific source(s) is(are) to be removed from
            :py:attr:`sources` attribute in the specified *direction*.
        """

        if check_param(direction):

            if source is not None:

                try:
                    self.sources[direction].remove(source)

                except ValueError:
                    pass

            else:
                self.sources[direction] = []

            if len(self.sources[direction]) == 0:
                self.dirs[direction] = False

    def is_directed(self):
        """
        Checks if edge has any directionality information.

        :return:
            (*bool*) -- Returns ``True```if any of the :py:attr:`dirs`
            attribute values is ``True`` (except ``'undirected'``),
            ``False`` otherwise.
        """

        return bool(
            sum([v for k, v in iteritems(self.dirs) if k != 'undirected']))

    def is_stimulation(self, direction=None):
        """
        Checks if any (or for a specific *direction*) interaction is
        activation (positive interaction).

        :arg tuple direction:
            Optional, ``None`` by default. If specified, checks the
            :py:attr:`positive` attribute of that specific
            directionality. If not specified, checks both.

        :return:
            (*bool*) -- ``True`` if any interaction (or the specified
            *direction*) is activatory (positive).
        """

        if direction is None:
            return bool(sum(self.positive.values()))

        else:
            return self.positive[direction]

    def is_inhibition(self, direction=None):
        """
        Checks if any (or for a specific *direction*) interaction is
        inhibition (negative interaction).

        :arg tuple direction:
            Optional, ``None`` by default. If specified, checks the
            :py:attr:`negative` attribute of that specific
            directionality. If not specified, checks both.

        :return:
            (*bool*) -- ``True`` if any interaction (or the specified
            *direction*) is inhibitory (negative).
        """

        if direction is None:
            return bool(sum(self.negative.values()))

        else:
            return self.negative[direction]

    def has_sign(self, direction=None):
        """
        Checks whether the edge (or for a specific *direction*) has
        any signed information (about positive/negative interactions).

        :arg tuple direction:
            Optional, ``None`` by default. If specified, only the
            information of that direction is checked for sign.

        :return:
            (*bool*) -- ``True`` if there exist any information on the
              sign of the interaction, ``False`` otherwise.
        """

        if direction is None:
            return (bool(sum(self.positive.values()))
                    or bool(sum(self.negative.values())))

        else:
            return self.negative[direction] or self.positive[direction]

    def set_sign(self, direction, sign, source):
        """
        Sets sign and source information on a given direction of the
        edge. Modifies the attributes :py:attr:`positive` and
        :py:attr:`positive_sources` or :py:attr:`negative` and
        :py:attr:`negative_sources` depending on the sign. Direction is
        also updated accordingly, which also modifies the attributes
        :py:attr:`dirs` and :py:attr:`sources`.

        :arg tuple direction:
            Pair of edge nodes specifying the direction from which the
            information is to be set/updated.
        :arg str sign:
            Specifies the type of interaction. If ``'positive'``, is
            considered activation, otherwise, is assumed to be negative
            (inhibition).
        :arg set source:
            Contains the name(s) of the source(s) from which the
            information was obtained.
        """

        if self.check_nodes(direction) and len(source):
            self.set_dir(direction, source)
            source = common.addToSet(set([]), source)

            if sign == 'positive':
                self.positive[direction] = True
                self.positive_sources[direction] = \
                    self.positive_sources[direction] | source

            else:
                self.negative[direction] = True
                self.negative_sources[direction] = \
                    self.negative_sources[direction] | source

    def get_sign(self, direction, sign=None, sources=False):
        """
        Retrieves the sign information of the edge in the given
        diretion. If specified in *sign*, only that sign's information
        will be retrieved. If specified in *sources*, the sources of
        that information will be retrieved instead.

        :arg tuple direction:
            Contains the pair of nodes specifying the directionality of
            the edge from which th information is to be retrieved.
        :arg str sign:
            Optional, ``None`` by default. Denotes whether to retrieve
            the ``'positive'`` or ``'negative'`` specific information.
        :arg bool sources:
            Optional, ``False`` by default. Specifies whether to return
            the sources instead of sign.

        :return:
            (*list*) -- If ``sign=None`` containing [bool] values
            denoting the presence of positive and negative sign on that
            direction, if ``sources=True`` the [set] of sources for each
            of them will be returned instead. If *sign* is specified,
            returns [bool] or [set] (if ``sources=True``) of that
            specific direction and sign.
        """

        if self.check_nodes(direction):

            if sources:

                if sign == 'positive':
                    return self.positive_sources[direction]

                elif sign == 'negative':
                    return self.negative_sources[direction]

                else:
                    return [self.positive_sources[direction],
                            self.negative_sources[direction]]

            else:

                if sign == 'positive':
                    return self.positive[direction]

                elif sign == 'negative':
                    return self.negative[direction]

                else:
                    return [self.positive[direction], self.negative[direction]]

    def unset_sign(self, direction, sign, source=None):
        """
        Removes sign and/or source information of the specified
        *direction* and *sign*. Modifies attribute :py:attr:`positive`
        and :py:attr:`positive_sources` or :py:attr:`negative` and
        :py:attr:`negative_sources` (or
        :py:attr:`positive_attributes`/:py:attr:`negative_sources`
        only if ``source=True``).

        :arg tuple direction:
            The pair of nodes specifying the directionality from which
            the information is to be removed.
        :arg str sign:
            Sign from which the information is to be removed. Must be
            either ``'positive'`` or ``'negative'``.
        :arg set source:
            Optional, ``None`` by default. If specified, determines
            which source(s) is(are) to be removed from the sources in
            the specified *direction* and *sign*.
        """

        if self.check_nodes(direction):

            if source is not None:

                try:

                    if sign == 'positive':
                        self.positive_sources[direction].remove(source)

                    if sign == 'negative':
                        self.negative_sources[direction].remove(source)

                except:
                    pass

            else:

                if sign == 'positive':
                    self.positive_sources[direction] = []

                if sign == 'negative':
                    self.negative_sources[direction] = []

            if len(self.positive_sources[direction]) == 0:
                self.positive[direction] = False

            if len(self.negative_sources[direction]) == 0:
                self.negative[direction] = False

    # XXX: Not sure if intended or you noticed, but if undirected=True and
    #      self.dirs['undirected']=True, the returning list includes 'u'.
    def src(self, undirected=False):
        """
        Returns the name(s) of the source node(s) for each existing
        direction on the interaction.

        :arg bool undirected:
            Optional, ``False`` by default.

        :returns:
            (*list*) -- Contains the name(s) for the source node(s).
            This means if the interaction is bidirectional, the list
            will contain both identifiers on the edge. If the
            interaction is undirected, an empty list will be returned.
        """

        return [k[0] for k, v in iteritems(self.dirs)
                if (k != 'undirected' or undirected) and v]

    # XXX: Similarly, if undirected=True and self.dirs['undirected']=True,
    #      the returning list includes 'n'.
    def tgt(self, undirected=False):
        """
        Returns the name(s) of the target node(s) for each existing
        direction on the interaction.

        :arg bool undirected:
            Optional, ``False`` by default.

        :returns:
            (*list*) -- Contains the name(s) for the target node(s).
            This means if the interaction is bidirectional, the list
            will contain both identifiers on the edge. If the
            interaction is undirected, an empty list will be returned.
        """

        return [k[1] for k, v in iteritems(self.dirs)
                if (k != 'undirected' or undirected) and v]

    # FIXME: If the intended behavior of the functions above is to include
    #        'undirected' in the returned list if undirected=True and
    #        self.dirs['undirected']=True, then change return by:
    #
    #        `[k[0] if k != 'undirected' else k for k, v in
    #         iteritems(self.dirs) if (k != 'undirected' or undirected) and v]`
    #
    #        On the other hand, if the intended behavior is that returns an
    #        empty list if kwarg undirected=True, then add a simple if/else
    #        return

    def src_by_source(self, source):
        """
        Returns the name(s) of the source node(s) for each existing
        direction on the interaction for a specific *source*.

        :arg str source:
            Name of the source according to which the information is to
            be retrieved.

        :return:
            (*list*) -- Contains the name(s) for the source node(s)
            according to the specified *source*. This means if the
            interaction is bidirectional, the list will contain both
            identifiers on the edge. If the specified *source* is not
            found or invalid, an empty list will be returned.
        """

        return [k[0] for k, v in iteritems(self.sources)
                if k != 'undirected' and source in v]

    def tgt_by_source(self, source):
        """
        Returns the name(s) of the target node(s) for each existing
        direction on the interaction for a specific *source*.

        :arg str source:
            Name of the source according to which the information is to
            be retrieved.

        :return:
            (*list*) -- Contains the name(s) for the target node(s)
            according to the specified *source*. This means if the
            interaction is bidirectional, the list will contain both
            identifiers on the edge. If the specified *source* is not
            found or invalid, an empty list will be returned.
        """

        return [k[1] for k, v in iteritems(self.sources)
                if k != 'undirected' and source in v]

    def sources_straight(self):
        """
        Retrieves the list of sources for the :py:attr:`straight`
        direction.

        :return:
            (*set*) -- Contains the names of the sources supporting the
            :py:attr:`straight` directionality of the edge.
        """

        return self.sources[self.straight]

    def sources_reverse(self):
        """
        Retrieves the list of sources for the :py:attr:`reverse` direction.

        :return:
            (*set*) -- Contains the names of the sources supporting the
            :py:attr:`reverse` directionality of the edge.
        """

        return self.sources[self.reverse]

    def sources_undirected(self):
        """
        Retrieves the list of sources without directed information.

        :return:
            (*set*) -- Contains the names of the sources supporting the
            edge presence but without specific directionality
            information.
        """

        return self.sources['undirected']

    def positive_straight(self):
        """
        Checks if the :py:attr:`straight` directionality is a positive
        interaction.

        :return:
            (*bool*) -- ``True`` if there is supporting information on
            the :py:attr:`straight` direction of the edge as activation.
            ``False`` otherwise.
        """

        return self.positive[self.straight]

    def positive_reverse(self):
        """
        Checks if the :py:attr:`reverse` directionality is a positive
        interaction.

        :return:
            (*bool*) -- ``True`` if there is supporting information on
            the :py:attr:`reverse` direction of the edge as activation.
            ``False`` otherwise.
        """

        return self.positive[self.reverse]

    def negative_straight(self):
        """
        Checks if the :py:attr:`straight` directionality is a negative
        interaction.

        :return:
            (*bool*) -- ``True`` if there is supporting information on
            the :py:attr:`straight` direction of the edge as inhibition.
            ``False`` otherwise.
        """

        return self.negative[self.straight]

    def negative_reverse(self):
        """
        Checks if the :py:attr:`reverse` directionality is a negative
        interaction.

        :return:
            (*bool*) -- ``True`` if there is supporting information on
            the :py:attr:`reverse` direction of the edge as inhibition.
            ``False`` otherwise.
        """

        return self.negative[self.reverse]

    def negative_sources_straight(self):
        """
        Retrieves the list of sources for the :py:attr:`straight`
        direction and negative sign.

        :return:
            (*set*) -- Contains the names of the sources supporting the
            :py:attr:`straight` directionality of the edge with a
            negative sign.
        """

        return self.negative_sources[self.straight]

    def negative_sources_reverse(self):
        """
        Retrieves the list of sources for the :py:attr:`reverse`
        direction and negative sign.

        :return:
            (*set*) -- Contains the names of the sources supporting the
            :py:attr:`reverse` directionality of the edge with a
            negative sign.
        """

        return self.negative_sources[self.reverse]

    def positive_sources_straight(self):
        """
        Retrieves the list of sources for the :py:attr:`straight`
        direction and positive sign.

        :return:
            (*set*) -- Contains the names of the sources supporting the
            :py:attr:`straight` directionality of the edge with a
            positive sign.
        """

        return self.positive_sources[self.straight]

    def positive_sources_reverse(self):
        """
        Retrieves the list of sources for the :py:attr:`reverse`
        direction and positive sign.

        :return:
            (*set*) -- Contains the names of the sources supporting the
            :py:attr:`reverse` directionality of the edge with a
            positive sign.
        """

        return self.positive_sources[self.reverse]

    def majority_dir(self):
        """
        Infers which is the major directionality of the edge by number
        of supporting sources.

        :return:
            (*tuple*) -- Contains the pair of nodes denoting the
            consensus directionality. If the number of sources on both
            directions is equal, ``None`` is returned. If there is no
            directionality information, ``'undirected'``` will be
            returned.
        """

        if self.is_directed():

            if len(self.sources[self.straight]) == len(self.sources[
                    self.reverse]):
                return None

            elif len(self.sources[self.straight]) > len(self.sources[
                    self.reverse]):
                return self.straight

            else:
                return self.reverse

        else:
            return 'undirected'

    def majority_sign(self):
        """
        Infers which is the major sign (activation/inhibition) of the
        edge by number of supporting sources on both directions.

        :return:
            (*dict*) -- Keys are the node tuples on both directions
            (:py:attr:`straight`/:py:attr:`reverse`) and values can be
            either ``None`` if that direction has no sign information or
            a list of two [bool] elements corresponding to majority of
            positive and majority of negative support. In case both
            elements of the list are ``True``, this means the number of
            supporting sources for both signs in that direction is
            equal.
        """

        dirs = [self.straight, self.reverse]

        result = {(d, None) for d in dirs}

        for d in dirs:
            if self.has_sign(direction=d):
                p, n = map(len, [self.positive_sources[d],
                                 self.negative_sources[d]])

                result[d] = [p >= n, p <= n]

        return result

    def consensus_edges(self):
        """
        Infers the consensus edge(s) according to the number of
        supporting sources. This includes direction and sign.

        :return:
            (*list*) -- Contains the consensus edge(s) along with the
            consensus sign. If there is no major directionality, both
            are returned. The structure is as follows:
            ``['<source>', '<target>', '<(un)directed>', '<sign>']``
        """

        result = []
        d = self.majority_dir()
        s = self.majority_sign()

        # XXX: This case could actually return directly, would save some time
        if d == 'undirected':
            result.append([self.straight[0], self.straight[1],
                           'undirected', 'unknown'])

        else:
            dirs = [self.straight, self.reverse] if d is None else [d]

            for d in dirs:

                if s[d] is not None:

                    if s[d][0]:
                        result.append([d[0], d[1], 'directed', 'positive'])

    # XXX: Technically, if s[d] is not None, this could be an elif right?
                    if s[d][1]:
                        result.append([d[0], d[1], 'directed', 'negative'])

                else:
                    result.append([d[0], d[1], 'directed', 'unknown'])

        return result

    def merge(self, other):
        """
        Merges current edge with another (if and only if they are the
        same class and contain the same nodes). Updates the attributes
        :py:attr:`dirs`, :py:attr:`sources`, :py:attr:`positive`,
        :py:attr:`negative`, :py:attr:`positive_sources` and
        :py:attr:`negative_sources`.

        :arg pypath.main.Direction other:
            The new edge object to be merged with the current one.
        """

    # XXX: Not best way to check the class. Probably never happens, but there
    #      may be other class out there with the same name
    #      >  if (other.__class__.__name__ == 'Direction' and self.check_nodes(
    #                 other.nodes):
        if other.__class__ == self.__class__ and self.check_nodes(other.nodes):
            for k in [self.straight, self.reverse, 'undirected']:
                self.dirs[k] = self.dirs[k] or other.dirs[k]

                self.sources[k] = self.sources[k] | other.sources[k]

    # XXX: Is there a reason to only update positive with straight and negative
    #      only with reverse?
                if k == self.straight:
                    self.positive[k] = self.positive[k] or other.positive[k]
                    self.positive_sources[k] = (self.positive_sources[k]
                                                | other.positive_sources[k])

                elif k == self.reverse:
                    self.negative[k] = self.negative[k] or other.negative[k]
                    self.negative_sources[k] = (self.negative_sources[k]
                                                | other.negative_sources[k])

    def translate(self, ids):
        """
        Translates the node names/identifiers according to the
        dictionary *ids*.

        :arg dict ids:
            Dictionary containing (at least) the current names of the
            nodes as keys and their translation as values.

        :return:
            (*pypath.main.Direction*) -- The copy of current edge object
            with translated node names.
        """

        # new Direction object
        newd = Direction(ids[self.nodes[0]], ids[self.nodes[1]])

        # copying directions
        for k, v in iteritems(self.sources):
            di = (ids[k[0]], ids[k[1]]) if type(k) is tuple else k
            newd.set_dir(di, v)

        # copying signs
        for di in [self.straight, self.reverse]:
            pos, neg = self.get_sign(di, sources = True)

            newd.set_sign((ids[di[0]], ids[di[1]]), 'positive', pos)
            newd.set_sign((ids[di[0]], ids[di[1]]), 'negative', neg)

        return newd


# TODO: Ask Dénes about this and finish docstrings
class AttrHelper(object):
    """
    Attribute helper class.

    * Initialization arguments:
        - *value* [dict/str]?:
        - *name* [str]?: Optional, ``None`` by default.
        - *defaults* [dict]:

    * Attributes:
        - *value* [dict]?:
        - *name* [str]?:
        - *defaults* [dict]:
        - *id_type* [type]:

    * Call arguments:
        - *instance* []:
        - *thisDir* [tuple?]: Optional, ``None`` by default.
        - *thisSign* []: Optional, ``None`` by default.
        - *thisDirSources* []: Optional, ``None`` by default.
        - *thisSources* []: Optional, ``None`` by default.

    * Returns:
        -
    """

    def __init__(self, value, name=None, defaults={}):
        """
        """

        self.name = name
        self.value = value
        self.defaults = defaults

        if isinstance(self.value, dict):
            self.id_type = type(self.value.keys()[0])

    def __call__(self, instance, thisDir=None, thisSign=None,
                 thisDirSources=None, thisSources=None):
        """
        """

        _thisDir = 'directed' if isinstance(thisDir, tuple) else thisDir

        # user supplied callback function:
        if hasattr(self.value, '__call__'):
            return self.value(instance)

        # special cases #1: by direction/effect
        elif (self.value == 'DIRECTIONS' and self.defaults is not None
              and self.name is not None and self.name in self.defaults):

            if _thisDir in self.defaults[self.name]:

                if thisSign in self.defaults[self.name][_thisDir]:
                    return self.defaults[self.name][_thisDir][thisSign]

        # special cases #2: by source category
        elif self.value == 'RESOURCE_CATEGORIES':

            for resource_type in ['pathway', 'ptm', 'reaction', 'interaction']:

                if len(getattr(data_formats, '%s_resources' % resource_type)
                       &thisSources) > 0:

                    if (self.name in self.defaults
                        and resource_type in self.defaults[self.name]):
                        return self.defaults[self.name][resource_type]

            sys.stdout.wrtie('No category for %s\n' % thisSources)
            sys.stdout.flush()

        # if value is constant:
        elif type(self.value) in common.simpleTypes:
            return self.value

        # if a dictionary given to map some igraph attribute to values:
        elif hasattr(self.value, '__call__'):
            return self.value(instance)

        elif isinstance(self.value, dict) and self.attr_name is not None:

            if hasattr(instance, self.value['_name']):
                key_attr = getattr(instance, self.value['_name'])

            elif self.value['_name'] in instance.attributes():
                key_attr = instance[self.value['_name']]

            if key_attr in self.value:
                return self.value[key_attr]

        # if default value has been given for this attribute:
        elif (self.name is not None and self.defaults is not None
              and self.name in self.defaults):
            return self.defaults[self.name]

        # ultimately, return None
        else:
            return None


class _NamedVertexSeq(object):
    """
    Vertex sequence object. Combines the list of vertex objects, their
    UniProt IDs and corresponding GeneSymbols.

    :arg igraph.VertexSeq _vs:
        Collection of :py:class:`igraph.Vertex` objects.
    :arg list _nodNam:
        List of [str] containing the node names (UniProt IDs).
    :arg list _nodLab:
        List of [str] containing the node labels (GeneSymbols).

    :var igraph.VertexSeq _vs:
        Collection of :py:class:`igraph.Vertex` objects.
    :var lsit _nodNam:
        List of [str] containing the node names (UniProt IDs).
    :var list _nodLab:
        List of [str] containing the node labels (GeneSymbols).
    """

    __slots__ = ['_vs', '_nodNam', '_nodLab']

    def __init__(self, _vs, _nodNam, _nodLab):
        """
        Initializes the object and sets the attributes according to the
        passed arguments.
        """

        self._vs = _vs
        self._nodNam = _nodNam
        self._nodLab = _nodLab

    def __iter__(self):
        """
        Iterator function yielding the :py:class:`igraph.Vertex`
        instances contained in *_vs*. It's accessed through the alias
        *vs*.

        :yield:
            (*igraph.Vertex*) -- The :py:class:`igraph.Vertex` object.
        """

        for v in self._vs:
            yield v

    def genesymbol(self):
        """
        Iterator function yielding the GeneSymbols contained in
        *_nodLab*. It can be accessed through the alias *gs*.

        :yield:
            (*str*) -- The GeneSymbols from the node list.
        """

        for v in self._vs:
            yield self._nodLab[v.index]

    def uniprot(self):
        """
        Iterator function yielding the UniProt IDs contained in
        *_nodNam*. It can be accessed through the alias *up*.

        :yield:
            (*str*) -- The UniProt IDs from the node list.
        """

        for v in self._vs:
            yield self._nodNam[v.index]

    def ids(self):
        """
        Iterator function yielding the vertex indexes.

        :yield:
            (*int*) -- The node's index.
        """

        for v in self._vs:
            yield v.index

    # Aliases
    gs = genesymbol
    up = uniprot
    vs = __iter__

# TODO: Ask:
#       - Which/how many organisms are accepted/available?
#       - MySQL
class PyPath(object):
    """Main network object.

    :arg int ncbi_tax_id:
        Optional, ``9606`` (Homo sapiens) by default. NCBI Taxonomic
        identifier of the organism from which the data will be
        downloaded.
    :arg dict default_name_type:
        Optional, ``{'protein': 'uniprot', 'mirna': 'mirbase', 'drug':
        'chembl', 'lncrna': 'lncrna-genesymbol'}`` by default. Contains
        the default identifier types to which the downloaded data will
        be converted. If others are used, user may need to provide the
        format definitions for the conversion tables.
    :arg pypath.main.PyPath copy:
        Optional, ``None`` by default. Other
        :py:class:`pypath.main.PyPath` instance from which the data will
        be copied.
    :arg tuple mysql:
        Optional, ``(None, 'mapping')`` by default. Contains the MySQL
        parameters used by the :py:mod:`pypath.mapping` module to load
        the ID conversion tables.
    :arg tuple chembl_mysql:
        Optional, ``(None, 'chembl')`` by default. Contains the MySQL
        parameters used by the :py:mod:`pypath.mapping` module to load
        the ChEMBL ID conversion tables.
    :arg str name:
        Optional, ``'unnamed'`` by default. Session or project name
        (custom).
    :arg str outdir:
        Optional, ``'results'`` by default. Output directory where to
        store all output files.
    :arg str loglevel:
        Optional, ``'INFO'`` by default. Sets the level of the logger.
        Possible levels are: ``'DEBUG'``, ``'INFO'``, ``'WARNING'``,
        ``'ERROR'`` or ``'CRITICAL'``.
    :arg bool loops:
        Optional, ``False`` by default. Determines if self-loop edges
        are allowed in the graph.

    :var list adjlist:
        List of [set] containing the adjacency of each node. See
        :py:meth:`PyPath.update_adjlist` method for more information.
    :var pypath.chembl.Chembl chembl:
        Contains the ChEMBL data. See :py:mod:`pypath.chembl` module
        documentation for more information.
    :var tuple chembl_mysql:
        Contains the MySQL parameters used by the
        :py:mod:`pypath.mapping` module to load the ChEMBL ID conversion
        tables.
    :var dict data:
        Stores the loaded interaction and attribute table. See
        :py:meth:`PyPath.read_data_file` method for more information.
    :var dict db_dict:
        Dictionary of dictionaries. Outer-level keys are ``'nodes'`` and
        ``'edges'``, corresponding values are [dict] whose keys are the
        database sources with values of type [set] containing the
        edge/node indexes for which that database provided some
        information.
    :var igraph.Graph dgraph:
        Directed network graph object.
    :var str disclaimer:
        Disclaimer text.
    :var dict dlabDct:
        Maps the directed graph node labels [str] (keys) to their
        indices [int] (values).
    :var dict dnodDct:
        Maps the directed graph node names [str] (keys) to their indices
        [int] (values).
    :var set dnodInd:
        Stores the directed graph node names [str].
    :var dict dnodLab:
        Maps the directed graph node indices [int] (keys) to their
        labels [str] (values).
    :var dict dnodNam:
        Maps the directed graph node indices [int] (keys) to their names
        [str] (values).
    :var dict edgeAttrs:
        Stores the edge attribute names [str] as keys and their
        corresponding types (e.g.: ``set``, ``list``, ``str``, ...) as
        values.
    :var pandas.DataFrame exp:
        Stores the expression data for the nodes (if loaded).
    :var pandas.DataFrame exp_prod:
        Stores the edge expression data (as the product of the
        normalized expression between the pair of nodes by default). For
        more details see :py:meth:`PyPath.edges_expression`.
    :var set exp_samples:
        Contains a list of tissues as downloaded by ProteomicsDB. See
        :py:meth:`PyPath.get_proteomicsdb` for more information.
    :var list failed_edges:
        List of lists containing information about the failed edges.
        Each failed edge sublist contains (in this order): [tuple] with
        the node IDs, [str] names of nodes A and B, [int] IDs of nodes
        A and B and [int] IDs of the edges in both directions.
    :var dict go:
        Contains the organism(s) NCBI taxonomy ID as key [int] and
        :py:class:`pypath.go.GOAnnotation` object as value, which
        contains the GO annotations for the nodes in the graph. See
        :py:class:`pypath.go.GOAnnotation` for more information.
    :var igraph.Graph graph:
        Undirected network graph object.
    :var pypath.gsea.GSEA gsea:
        Contains the loaded gene-sets from MSigDB. See
        :py:class:`pypath.gsea.GSEA` for more information.
    :var set has_cats:
        Contains the categories (e.g.: resources) [str] loaded in the
        current network.
    :var dict htp:
        Contains information about high-throughput data of the network
        for different thresholds [int] (keys). Values are [dict]
        containing the number of references (``'rnum'``) [int], number
        of edges (``'enum'``) [int], number of sources (``'snum'``)
        [int] and list of PMIDs of the most common references above the
        given threshold (``'htrefs'``) [set].
    :var dict labDct:
        Maps the undirected graph node labels [str] (keys) to their
        indices [int] (values).
    :var dict lists:
        Contains specific lists of nodes (values) for different
        categories [str] (keys). These can to be loaded from a file or
        a resource. Some methods include :py:meth:`PyPath.receptor_list`
        (``'rec'``), :py:meth:`PyPath.druggability_list` (``'dgb'``),
        :py:meth:`PyPath.kinases_list` (``'kin'``),
        :py:meth:`PyPath.tfs_list` (``'tf'``),
        :py:meth:`PyPath.disease_genes_list` (``'dis'``),
        :py:meth:`PyPath.signaling_proteins_list` (``'sig'``),
        :py:meth:`PyPath.proteome_list` (``'proteome'``) and
        :py:meth:`PyPath.cancer_drivers_list` (``'cdv'``).
    :var str loglevel:
        The level of the logger.
    :var bool loops:
        Whether if self-loop edges are allowed in the graph.
    :var pypath.mapping.Mapper mapper:
        :py:class:`pypath.mapper.Mapper` object for ID conversion and
        other ID-related operations across resources.
    :var list mutation_samples:
        DEPRECATED
    :var tuple mysql_conf:
        Contains the MySQL parameters used by the
        :py:mod:`pypath.mapping` module to load the ID conversion
        tables.
    :var str name:
        Session or project name (custom).
    :var int ncbi_tax_id:
        NCBI Taxonomic identifier of the organism from which the data
        will be downloaded.
    :var dict negatives:
        Contains a list of negative interactions according to a given
        source (e.g.: Negatome database). See
        :py:meth:`PyPath.apply_negative` for more information.
    :var dict nodDct:
        Maps the undirected graph node names [str] (keys) to their
        indices [int] (values).
    :var set nodInd:
        Stores the undirected graph node names [str].
    :var dict nodLab:
        Maps the undirected graph node indices [int] (keys) to their
        labels [str] (values).
    :var dict nodNam:
        Maps the directed graph node indices [int] (keys) to their names
        [str] (values).
    :var str outdir:
        Output directory where to store all output files.
    :var pypath.logn.logw ownlog:
        Logger class instance, see :py:class:`pypath.logn.logw` for more
        information.
    :var list palette:
        Contains a list of hexadecimal [str] of colors. Used for
        plotting purposes.
    :var list pathway_types:
        Contains the names of all the loaded pathway resources [str].
    :var dict pathways:
        Contains the list of pathways (values) for each resource (keys)
        loaded in the network.
    :var dict plots:
        DEPRECATED (?)
    :var pypath.proteomicsdb.ProteomicsDB proteomicsdb:
        Contains a :py:class:`pypath.proteomicsdb.ProteomicsDB`
        instance, see the class documentation for more information.
    :var list raw_data:
        Contains a list of loaded edges [dict] from a data file. See
        :py:meth:`PyPath.read_data_file` for more information.
    :var dict reflists:
        Contains the reference list(s) loaded. Keys are [tuple]
        containing the node name type [str] (e.g.: ``'uniprot'``), type
        [str] (e.g.: ``'protein'``) and taxonomic ID [int] (e.g.:
        ``'9606'``). Values are the corresponding
        :py:class:`pypath.reflists.ReferenceList` instance (see class
        documentation for more information).
    :var dict seq:
        (?)
    :var str session:
        Session ID, a five random alphanumeric characters.
    :var str session_name:
        Session name and ID (e.g. ``'unnamed-abc12'``).
    :var igraph.Graph sourceNetEdges:
        (?)
    :var igraph.Graph sourceNetNodes:
        (?)
    :var list sources:
        List contianing the names of the loaded resources [str].
    :var dict u_pfam:
        Dictionary of dictionaries, contains the mapping of UniProt IDs
        to their respective protein families and other information.
    :var list uniprot_mapped:
        DEPRECATED (?)
    :var list unmapped:
        Contains the names of unmapped items [str]. See
        :py:meth:`PyPath.map_item` for more information.
    :var dict vertexAttrs:
        Stores the node attribute names [str] as keys and their
        corresponding types (e.g.: ``set``, ``list``, ``str``, ...) as
        values.
    """

    default_name_type = {'protein': 'uniprot',
                         'mirna': 'mirbase',
                         'drug': 'chembl',
                         'lncrna': 'lncrna-genesymbol'}

    def __init__(self, ncbi_tax_id=9606, default_name_type=default_name_type,
                 copy=None, mysql=(None, 'mapping'),
                 chembl_mysql=(None, 'chembl'), name='unnamed',
                 cache_dir = None,
                 outdir='results', loglevel='INFO', loops=False):
        """Initializes the network object.

        **NOTE:** Only the instance is created, no data is donwloaded
        until the corresponding function is called (e.g.:
        :py:meth:`PyPath.init_network`).
        """

        self.__version__ = _version.__version__

        # Setting up the working directory
        for d in ['results', 'log']:

            if not os.path.exists(d):
                os.makedirs(d)

        self.cache_dir = cache_dir or settings.get('cachedir')

        if not os.path.exists(self.cache_dir):

            os.makedirs(self.cache_dir)

        if copy is None:
            # Setting up graph object
            self.graph = igraph.Graph(0)
            g = self.graph
            g['entity_types'] = {}
            g['ncbi_tax_id'] = ncbi_tax_id
            g['name'] = name
            g['sources'] = {}
            g['references'] = {}
            g['directed'] = False
            g.vs['type'] = []
            g.vs['name'] = []
            g.vs['nameType'] = []
            g.vs['originalNames'] = [[] for _ in xrange(self.graph.vcount())]
            g.vs['ncbi_tax_id'] = []
            g.vs['exp'] = [{}]
            g.es['sources'] = [set([]) for _ in xrange(self.graph.ecount())]
            g.es['type'] = [[] for _ in xrange(self.graph.ecount())]
            g.es['references'] = [[] for _ in xrange(self.graph.ecount())]
            g.es['refs_by_source'] = [{} for _ in xrange(self.graph.ecount())]
            g.es['refs_by_dir'] = [{} for _ in xrange(self.graph.ecount())]
            g.es['refs_by_type'] = [{} for _ in xrange(self.graph.ecount())]
            g.es['sources_by_type'] = [{} for _ in xrange(self.graph.ecount())]
            g.es['negative_refs'] = [[] for _ in xrange(self.graph.ecount())]
            g.es['negative'] = [[] for _ in xrange(self.graph.ecount())]
            g.es['dirs'] = [None]
            g['layout_type'] = None
            g['layout_data'] = None
            g['only_directed'] = False

            self.loops = loops
            self.dgraph = None
            self._undirected = self.graph
            self._directed = None
            self.failed_edges = []

            self.uniprot_mapped = [] # XXX: Not used anywhere
            self.mysql_conf = mysql
            self.set_chembl_mysql(chembl_mysql[1], chembl_mysql[0])
            # self.mysql = mysql.MysqlRunner(self.mysql_conf)
            self.unmapped = []
            self.name = name
            self.outdir = outdir
            self.ncbi_tax_id = ncbi_tax_id
            self.data = {}
            self.reflists = {}
            self.negatives = {}
            self.raw_data = None
            self.lists = {}
            self.plots = {} # XXX: Not used anywhere
            self.proteomicsdb = None
            self.exp_samples = set([])
            self.sources = []
            self.has_cats = set([])
            self.db_dict = {}
            self.pathway_types = []
            self.pathways = {}
            self.vertexAttrs = {}
            self.edgeAttrs = {}
            self.u_pfam = None
            self.seq = None
            self.palette = ['#6EA945', '#007B7F', '#FCCC06', '#DA0025',
                            '#000000']

            # Session and log
            self.session = common.gen_session_id()
            self.session_name = ''.join([self.name, '-', self.session]) # XXX: What about '-'.join([self.name, self.session])?
            self.loglevel = loglevel
            self.ownlog = logn.logw(self.session, self.loglevel)
            self.mapper = mapping.Mapper(self.ncbi_tax_id,
                                         mysql_conf=self.mysql_conf,
                                         log=self.ownlog)
            self.disclaimer = '\n\t=== d i s c l a i m e r ===\n\n'\
                '\tAll data coming with this module\n'\
                '\teither as redistributed copy or downloaded using the\n'\
                '\tprogrammatic interfaces included in the present module\n'\
                '\tare available under public domain, are free to use at\n'\
                '\tleast for academic research or education purposes.\n'\
                '\tPlease be aware of the licences of all the datasets\n'\
                '\tyou use in your analysis, and please give appropriate\n'\
                '\tcredits for the original sources when you publish your\n'\
                '\tresults. To find out more about data sources please\n'\
                '\tlook at `pypath.descriptions` and\n'\
                '\t`pypath.data_formats.urls`.\n\n'
            self.licence()
            self.ownlog.msg(1, "PyPath has been initialized")
            self.ownlog.msg(1, "Beginning session '%s'" % self.session)
            sys.stdout.write(
                """\t> New session started,\n\tsession ID: '%s'\n\tlogfile: """
                """'./%s'\n\tpypath version: %s\n""" % (
                    self.session, self.ownlog.logfile, _version.__version__))

        else:
            self.copy(copy)

    def set_chembl_mysql(self, title, config_file=None):
        """
        Sets the ChEMBL MySQL configuration according to the *title*
        section in *config_file* ini file configuration.

        :arg str title:
            Section title of the ini file.
        :arg str config_file:
            Optional, ``None`` by default. Specifies the configuration
            file name if none is passed, ``mysql_config/defaults.mysql``
            will be used.
        """

        self.chembl_mysql = (config_file, title)

    def copy(self, other):
        """
        Copies another :py:class:`pypath.main.PyPath` instance into the
        current one.

        :arg pypath.main.PyPath other:
            The instance to be copied from.
        """

        self.__dict__ = other.__dict__
        self.ownlog.msg(1, "Reinitialized", 'INFO')

    def init_network(self, lst=None, exclude=[], cache_files={},
                     pfile=False, save=False, reread=False, redownload=False,
                     **kwargs): # XXX: kwargs is not used anywhere
        """Loads the network data.

        This is a lazy way to start the module, load data and build the
        high confidence, literature curated part of the signaling
        network.

        :arg dict lst:
            Optional, ``None`` by default. Specifies the data input
            formats for the different resources (keys) [str]. Values
            are :py:class:`pypath.input_formats.ReadSettings` instances
            containing the information. By default uses the set of
            resources of OmniPath.
        :arg list exclude:
            Optional, ``[]`` by default. List of resources [str] to
            exclude from the network.
        :arg dict cache_files:
            Optional, ``{}`` by default. Contains the resource name(s)
            [str] (keys) and the corresponding cached file name [str].
            If provided (and file exists) bypasses the download of the
            data for that resource and uses the cache file instead.
        :arg str pfile:
            Optional, ``False`` by default. If any, provides the file
            name or path to a previously saved network pickle file.
            If ``True`` is passed, takes the default path from
            :py:meth:`PyPath.save_network`
            (``'cache/default_network.pickle'``).
        :arg bool save:
            Optional, ``False`` by default. If set to ``True``, saves
            the loaded network to its default location
            (``'cache/default_network.pickle'``).
        :arg bool reread:
            Optional, ``False`` by default. Specifies whether to reread
            the data files from the cache or omit them (similar to
            *redownload*).
        :arg bool redownload:
            Optional, ``False`` by default. Specifies whether to
            re-download the data and ignore the cache.
        :arg \*\*kwargs:
            Not used.
        """

        if lst is None:
            lst = omnipath

        if pfile:
            pfile = pfile if not isinstance(pfile, bool) \
                else os.path.join(self.cache_dir, 'default_network.pickle')

            if os.path.exists(pfile):
                sys.stdout.write(
                    '\t:: Loading igraph object from file `%s`...' % pfile)
                sys.stdout.flush()
                graph = pickle.load(open(pfile, 'rb'))

                if isinstance(graph, igraph.Graph) and graph.vcount() > 0:
                    self.graph = graph
                    sys.stdout.write(
                        '\r%s\r\t:: Network loaded from `%s`. %u nodes, '
                        '%u edges.\n' % (' ' * 90, pfile, self.graph.vcount(),
                                         self.graph.ecount()))
                    sys.stdout.flush()
                    self.update_vname()
                    self.update_vindex()
                    self.update_sources()

                    return None

        self.load_reflists() # XXX: This is redundant (see line 4565 in load_resources)
        self.load_resources(
            lst=lst, exclude=exclude, reread=reread, redownload=redownload,
            cache_files=cache_files)
        if save:
            sys.stdout.write('\t:: Saving igraph object to file `%s`...' %
                             pfile)
            sys.stdout.flush()
            self.save_network()
            sys.stdout.write(
                '\r%s\r\t:: Network saved successfully to file `%s`.\n' %
                (' ' * 90, pfile))
            sys.stdout.flush()

    def save_network(self, pfile=None):
        """Saves the network object.

        Stores the instance into a pickle (binary) file which can be
        reloaded in the future.

        :arg str pfile:
            Optional, ``None`` by default. The path/file name where to
            store the pcikle file. If not specified, saves the network
            to its default location
            (``'cache/default_network.pickle'``).
        """

        pfile = pfile if pfile is not None \
            else os.path.join(self.cache_dir, 'default_network.pickle')
        pickle.dump(self.graph, open(pfile, 'wb'), -1)

    ###
    # functions to read networks from text files or mysql
    ###

    def get_max(self, attrList): # TODO
        """
        """

        maxC = 0

        for val in attrList.values():

            if val.__class__ is tuple:
                val = val[0]

            if val > maxC:
                maxC = val

        return maxC

    def get_attrs(self, line, spec, lnum): # TODO
        """
        """

        attrs = {}

        for col in spec.keys():
            # extraEdgeAttrs and extraNodeAttrs are dicts
            # of additional parameters assigned to edges and nodes respectively;
            # key is the name of the parameter, value is the col number,
            # or a tuple of col number and the separator,
            # if the column contains additional subfields e.g. (5, ";")

            try:

                if spec[col].__class__ is tuple:

                    if hasattr(spec[col][1], '__call__'):
                        fieldVal = spec[col][1](line[spec[col][0]])

                    else:
                        fieldVal = line[spec[col][0]].split(spec[col][1])

                else:
                    fieldVal = line[spec[col]]

            except:
                self.ownlog.msg(2, (
                    """Wrong column index (%s) in extra attributes?
                    Line #%u\n""" % (str(col), lnum)), 'ERROR')
                readError = 1
                break

            fieldName = col
            attrs[fieldName] = fieldVal

        return attrs

    def get_taxon(self, tax_dict, fields): # TODO
        """
        """

        if 'A' in tax_dict and 'B' in tax_dict:

            return (self.get_taxon(tax_dict['A'], fields),
                    self.get_taxon(tax_dict['B'], fields))

        else:

            if 'dict' not in tax_dict:
                return int(fields[tax_dict['col']])

            elif fields[tax_dict['col']] in tax_dict['dict']:
                return tax_dict['dict'][fields[tax_dict['col']]]

            else:
                return None

    def numof_references(self):
        """Counts the number of reference on the network.

        Counts the total number of unique references in the edges of the
        network.

        :return:
            (*int*) -- Number of unique references in the network.
        """

        return len(
            common.uniqList(
                common.flatList(
                    list(map(lambda e: e['references'], self.graph.es)))))

    def mean_reference_per_interaction(self):
        """
        Computes the mean number of references per interaction of the
        network.

        :return:
            (*float*) -- Mean number of interactions per edge.
        """

        return np.mean(
            list(map(lambda e: len(e['references']), self.graph.es)))

    def numof_reference_interaction_pairs(self): # XXX: Not really sure about this one
        """
        Returns the total of unique references per interaction.

        :return:
            (*int*) -- Total number of unique references per
            interaction.
        """

        return len(common.uniqList(common.flatList(
            list(map(lambda e:
                     list(map(lambda r:
                              (e.index, r), e['references'])),
                     self.graph.es)))))

    def curators_work(self):
        """
        Computes and prints an estimation of how many years of curation
        took to achieve the amount of information on the network.
        """

        curation_effort = self.numof_reference_interaction_pairs()
        sys.stdout.write(
            '\t:: Curators worked %.01f-%.01f years to accomplish '
            'what currently you have incorporated in this network!'
            '\n\n\tAmazing, isn\'t it?\n' %
            (curation_effort * 15 / 60.0 / 2087.0,
             curation_effort * 60 / 60.0 / 2087.0))
        sys.stdout.flush()

    def reference_edge_ratio(self):
        """
        Computes the average number of references per edge (as in the
        undirected graph).

        :return:
            (*float*) -- Average number of references per edge.
        """

        return self.numof_references() / float(self.graph.ecount())

    def get_giant(self, replace=False, graph=None):
        """
        Returns the giant component of the *graph*, or replaces the
        :py:class:`igraph.Graph` instance with only the giant component
        if specified.

        :arg bool replace:
            Optional, ``False`` by default. Specifies whether to replace
            the :py:class:`igraph.Graph` instance. This can be either
            the undirected network of the current
            :py:class:`pypath.main.PyPath` instance (default) or the one
            passed under the keyword argument *graph*.
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The graph object from which
            the giant component is to be computed. If none is specified,
            takes the undirected network of the current
            :py:class:`pypath.main.PyPath` instance.

        :return:
            (*igraph.Graph*) -- If ``replace=False``, returns a copy of
            the giant component graph.
        """

        g = graph if graph is not None else self.graph
        gg = g if replace else copy.deepcopy(g)

        cl = gg.components(mode='WEAK')
        cl_sizes = cl.sizes()
        giant_component_index = cl_sizes.index(max(cl_sizes))
        in_giant = [x == giant_component_index for x in cl.membership]

        common.console(':: Nodes in giant component: %u' %
                       in_giant.count(True))

        toDel = [i for i in xrange(0, gg.vcount()) if not in_giant[i]]
        gg.delete_vertices(toDel)

        common.console(':: Giant component size: %u edges, %u nodes' %
                       (gg.ecount(), gg.vcount()))

        if not replace:
            return gg

    def update_vname(self):
        """
        Fast lookup of node names and indexes, these are hold in a
        [list] and a [dict] as well. However, every time new nodes are
        added, these should be updated. This function is automatically
        called after all operations affecting node indices.
        """

        self.genesymbol_labels()
        graph = self._get_undirected()
        self._already_has_directed()
        dgraph = self._directed

        if graph is not None:
            self.nodInd = set(graph.vs['name'])
            self.nodDct = dict(zip(graph.vs['name'], xrange(graph.vcount())))
            self.labDct = dict(zip(graph.vs['label'], xrange(graph.vcount())))
            self.nodNam = dict(zip(xrange(graph.vcount()), graph.vs['name']))
            self.nodLab = dict(zip(xrange(graph.vcount()), graph.vs['label']))

        if dgraph is not None:
            self.dnodInd = set(dgraph.vs['name'])
            self.dnodDct = dict(
                zip(dgraph.vs['name'], xrange(dgraph.vcount())))
            self.dlabDct = dict(
                zip(dgraph.vs['label'], xrange(dgraph.vcount())))
            self.dnodNam = dict(
                zip(xrange(dgraph.vcount()), dgraph.vs['name']))
            self.dnodLab = dict(
                zip(xrange(dgraph.vcount()), dgraph.vs['label']))

    def vsgs(self):
        """
        Returns a generator sequence of the node names as GeneSymbols
        [str] (from the undirected graph).

        :return:
            (*generator*) -- Sequence containing the node names as
            GeneSymbols [str].
        """

        return _NamedVertexSeq(self.graph.vs, self.nodNam, self.nodLab).gs()

    def vsup(self):
        """
        Returns a generator sequence of the node names as UniProt IDs
        [str] (from the undirected graph).

        :return:
            (*generator*) -- Sequence containing the node names as
            UniProt IDs [str].
        """

        return _NamedVertexSeq(self.graph.vs, self.nodNam, self.nodLab).up()

    def update_vindex(self): # XXX: If so, shouldn't it be removed?
        """
        This is deprecated.
        """

        self.nodNam = dict(
            zip(range(0, self.graph.vcount()), self.graph.vs['name']))

    def vertex_pathways(self):
        """
        Some resources assignes interactions some others proteins to
        pathways. This function copies pathway annotations from edge
        attributes to vertex attributes.
        """

        for eattr in self.graph.es.attributes():

            if eattr.endswith('pathways'):

                if eattr not in self.graph.vs.attributes():
                    self.graph.vs[eattr] = [[] for _ in self.graph.vs]

                for e in self.graph.es:
                    self.graph.vs[e.source][eattr] = e[eattr]
                    self.graph.vs[e.target][eattr] = e[eattr]

    # XXX: Not very clear for me what this function is actually doing...
    #      I mean, returns True/False and True as soon as the first
    #      len(thisVal & filtrVal) > 0
    def filters(self, line, positiveFilters=[], negativeFilters=[]): # TODO
        """
        """

        for filtr in negativeFilters:

            if len(filtr) > 2:
                sep = filtr[2]
                thisVal = set(line[filtr[0]].split(sep))

            else:
                thisVal = set([line[filtr[0]]])

            filtrVal = set(filtr[1]
                           if isinstance(filtr[1], list) else [filtr[1]])

            if len(thisVal & filtrVal) > 0:
                return True

        for filtr in positiveFilters:

            if len(filtr) > 2:
                sep = filtr[2]
                thisVal = set(line[filtr[0]].split(sep))

            else:
                thisVal = set([line[filtr[0]]])

            filtrVal = set(filtr[1]
                           if isinstance(filtr[1], list) else [filtr[1]])

            if len(thisVal & filtrVal) == 0:
                return True

        return False

    def lookup_cache(self, name, cache_files, int_cache, edges_cache):
        """
        Checks up the cache folder for the files of a given resource.
        First checks if *name* is on the *cache_files* dictionary.
        If so, loads either the interactions or edges otherwise. If
        not, checks *edges_cache* or *int_cache* otherwise.

        :arg str name:
            Name of the resource (lower-case).
        :arg dict cache_files:
            Contains the resource name(s) [str] (keys) and the
            corresponding cached file name [str] (values).
        :arg str int_cache:
            Path to the interactions cache file of the resource.
        :arg str edges_cache:
            Path to the edges cache file of the resource.

        :return:
            * (*file*) -- The loaded pickle file from the cache if the
              file is contains the interactions. ``None`` otherwise.
            * (*list*) -- List of mapped edges if the file contains the
              information from the edges. ``[]`` otherwise.
        """

        infile = None
        edgeListMapped = []
        cache_file = cache_files[name] if name in cache_files else None

        if cache_file is not None and os.path.exists(cache_file):
            cache_type = cache_file.split('.')[-2]

            if cache_type == 'interactions':
                infile = self.read_from_cache(int_cache)

            elif cache_type == 'edges':
                edgeListMapped = self.read_from_cache(edges_cache)

        elif os.path.exists(edges_cache):
            edgeListMapped = self.read_from_cache(edges_cache)

        else: # XXX: You could use another elif statement here

            if os.path.exists(int_cache):
                infile = self.read_from_cache(int_cache)

        return infile, edgeListMapped

    def read_from_cache(self, cache_file):
        """
        Reads a pickle file from the cache and returns it. It is assumed
        that the subfolder ``cache/`` is on the supplied path.

        :arg str cache_file:
            Path to the cache file that is to be loaded.

        :return:
            (*file*) -- The loaded pickle file from the cache. Type will
            depend on the file itself (e.g.: if the pickle was saved
            from a dictionary, the type will be [dict]).
        """

        sys.stdout.write('\t:: Reading from cache: %s\n' % cache_file)
        sys.stdout.flush()
        self.ownlog.msg(2, 'Data have been read from cache: %s' % cache_file)

        return pickle.load(open(cache_file, 'rb'))

    def process_sign(self, signData, signDef):
        """
        Processes the sign of an interaction, used when processing an
        input file.

        :arg str signData:
            Data regarding the sign to be processed.
        :arg tuple signDef:
            Contains information about how to process *signData*. This
            is defined in :py:mod:`pypath.data_formats`. First element
            determines the position on the direction information of each
            line on the data file [int], second element is either [str]
            or [list] and defines the terms for which an interaction is
            defined as stimulation, third element is similar but for the
            inhibition and third (optional) element determines the
            separator for *signData* if contains more than one element.

        :return:
            * (*bool*) -- Determines whether the processed interaction
              is considered stimulation or not.
            * (*bool*) -- Determines whether the processed interaction
              is considered inhibition or not.
        """

        stim = False
        inh = False
        signSep = signDef[3] if len(signDef) > 3 else None
        signData = set(str(signData).split(signSep))
        pos = set(signDef[1] if isinstance(signDef[1], list) else [signDef[1]])
        neg = set(signDef[2] if isinstance(signDef[2], list) else [signDef[2]])

        # XXX: Isn't using elif here bias the choice to stimulatory interactions
        #      even though there can also be negative sources?

        if len(signData & pos) > 0:
            stim = True

        elif len(signData & neg) > 0:
            inh = True

        return stim, inh

    def process_direction(self, line, dirCol, dirVal, dirSep):
        """
        Processes the direction information of an interaction according
        to a data file from a source.

        :arg list line:
            The stripped and separated line from the resource data file
            containing the information of an interaction.
        :arg int dirCol:
            The column/position number where the information about the
            direction is to be found (on *line*).
        :arg list dirVal:
            Contains the terms [str] for which that interaction is to be
            considered directed.
        :arg str dirSep:
            Separator for the field in *line* containing the direction
            information (if any).

        :return:
            (*bool*) -- Determines whether the given interaction is
            directed or not.
        """

        if dirCol is None or dirVal is None:
            return False

        else:
            thisDir = set(line[dirCol].split(dirSep))
            return len(thisDir & dirVal) > 0

    def read_data_file(self, settings, keep_raw=False, cache_files={},
                       reread=False, redownload=False):
        """
        Reads interaction data file containing node and edge attributes
        that can be read from simple text based files and adds it to the
        networkdata. This function works not only with files, but with
        lists as well. Any other function can be written to download and
        preprocess data, and then give it to this function to finally
        attach to the network.

        :arg pypath.input_formats.ReadSettings settings:
            :py:class:`pypath.input_formats.ReadSettings` instance
            containing the detailed definition of the input format of
            the file. Instead of the file name (on the
            :py:attr:`pypath.input_formats.ReadSettings.inFile`
            attribute) you can give a custom function name, which will
            be executed, and the returned data will be used instead.
        :arg bool keep_raw:
            Optional, ``False`` by default. Whether to keep the raw data
            read by this function, in order for debugging purposes, or
            further use.
        :arg dict cache_files:
            Optional, ``{}`` by default. Contains the resource name(s)
            [str] (keys) and the corresponding cached file name [str].
            If provided (and file exists) bypasses the download of the
            data for that resource and uses the cache file instead.
        :arg bool reread:
            Optional, ``False`` by default. Specifies whether to reread
            the data files from the cache or omit them (similar to
            *redownload*).
        :arg bool redownload:
            Optional, ``False`` by default. Specifies whether to
            re-download the data and ignore the cache.
        """

        listLike = set([list, tuple])
        edgeList = []
        nodeList = []
        edgeListMapped = []
        infile = None
        _name = settings.name.lower()
        int_cache = os.path.join(
            self.cache_dir,
            '%s.interactions.pickle' % _name
        )
        edges_cache = os.path.join(
            self.cache_dir,
            '%s.edges.pickle' % _name
        )
        if not reread and not redownload:
            infile, edgeListMapped = self.lookup_cache(_name, cache_files,
                                                       int_cache, edges_cache)

        if not len(edgeListMapped):

            if infile is None:

                if settings.__class__.__name__ != "ReadSettings":
                    self.ownlog.msg(2, (
                        """No proper input file definition!\n\'settings\'
                        should be a \'ReadSettings\' instance\n"""), 'ERROR')
                    return None

                if settings.huge:
                    sys.stdout.write(
                        '\n\tProcessing %s requires huge memory.\n'
                        '\tPlease hit `y` if you have at least 2G free memory,\n'
                        '\tor `n` to omit %s.\n'
                        '\tAfter processing once, it will be saved in \n'
                        '\t%s, so next time can be loaded quickly.\n\n'
                        '\tProcess %s now? [y/n]\n' %
                        (settings.name, settings.name, edges_cache,
                         settings.name))
                    sys.stdout.flush()

                    while True:
                        answer = raw_input().lower()

                        if answer == 'n':
                            return None

                        elif answer == 'y':
                            break

                        else:
                            sys.stdout.write(
                                '\n\tPlease answer `y` or `n`:\n\t')
                            sys.stdout.flush()

                inputFunc = self.get_function(settings.inFile)

                if inputFunc is None and hasattr(dataio, settings.inFile):
                    inputFunc = getattr(dataio, settings.inFile)

                # reading from remote or local file, or executing import
                # function:
                if settings.inFile.startswith('http') or \
                        settings.inFile.startswith('ftp'):
                    curl_use_cache = not redownload
                    c = curl.Curl(
                        settings.inFile,
                        silent=False,
                        large=True,
                        cache=curl_use_cache)
                    infile = c.result.read()

                    if type(infile) is bytes:

                        try:
                            infile = infile.decode('utf-8')

                        except:

                            try:
                                infile = infile.decode('iso-8859-1')

                            except:
                                pass

                    infile = [
                        x for x in infile.replace('\r', '').split('\n')
                        if len(x) > 0
                    ]
                    self.ownlog.msg(2, "Retrieving data from%s ..." %
                                    settings.inFile)

                # elif hasattr(dataio, settings.inFile):
                elif inputFunc is not None:
                    self.ownlog.msg(2, "Retrieving data by dataio.%s() ..." %
                                    inputFunc.__name__)
                    _store_cache = curl.CACHE
                    curl.CACHE = not redownload

                    # this try-except needs to be removed
                    # once correct exception handling will
                    # be implemented in every input function
                    try:
                        infile = inputFunc(**settings.inputArgs)

                    except Exception as e:
                        sys.stdout.write(
                            '\n\t:: Error in `pypath.dataio.%s()`. '
                            'Skipping to next resource.\n' %
                            inputFunc.__name__)
                        sys.stdout.write('\t:: %s\n' % str(e.args))
                        sys.stdout.flush()

                        try:
                            traceback.print_tb(
                                e.__traceback__, file=sys.stdout)

                        except Exception as e:
                            sys.stdout.write(
                                '\t:: Failed handling exception.\n')
                            sys.stdout.write('\t%s\n' % str(e.args))
                            sys.stdout.flush()

                    curl.CACHE = _store_cache

                elif os.path.isfile(settings.inFile):
                    infile = codecs.open(
                        settings.inFile, encoding='utf-8', mode='r')
                    self.ownlog.msg(2, "%s opened..." % settings.inFile)

                if infile is None:
                    self.ownlog.msg(2, "%s: No such file or "
                                    "dataio function! :(\n" %
                                    (settings.inFile), 'ERROR')
                    return None

            # finding the largest referred column number,
            # to avoid references out of range
            isDir = settings.isDirected
            sign = settings.sign
            refCol = settings.refs[0] if isinstance(settings.refs, tuple) \
                else settings.refs if isinstance(settings.refs, int) else None
            refSep = settings.refs[1] if isinstance(settings.refs,
                                                    tuple) else ';'
            sigCol = None if not isinstance(sign, tuple) else sign[0]
            dirCol = None
            dirVal = None
            dirSep = None

            if isinstance(isDir, tuple):
                dirCol = isDir[0]
                dirVal = isDir[1]
                dirSep = isDir[2] if len(isDir) > 2 else None

            elif isinstance(sign, tuple):
                dirCol = sign[0]
                dirVal = sign[1:3]
                dirVal = dirVal if type(dirVal[
                    0]) in common.simpleTypes else common.flatList(dirVal)
                dirSep = sign[3] if len(sign) > 3 else None

            dirVal = set(dirVal if isinstance(dirVal, list) else [dirVal])
            maxCol = max(
                filter(
                    lambda i: i is not None, [
                        settings.nameColA, settings.nameColB, self.get_max(
                            settings.extraEdgeAttrs),
                        self.get_max(settings.extraNodeAttrsA), self.get_max(
                            settings.extraNodeAttrsB), refCol, dirCol, sigCol,
                        max(itertools.chain(
                            map(lambda x: x[0],
                                settings.positiveFilters),
                            [0])),
                        max(itertools.chain(
                            map(lambda x: x[0],
                                settings.negativeFilters),
                            [0]))
                    ]))
            # iterating lines from input file
            lnum = 0
            lFiltered = 0
            rFiltered = 0
            tFiltered = 0
            readError = 0

            for line in infile: # XXX: here could be used enumerate for lnum
                lnum += 1

                if len(line) <= 1 or (lnum == 1 and settings.header):
                    # empty lines
                    # or header row
                    continue

                if type(line) not in listLike:

                    if hasattr(line, 'decode'):
                        line = line.decode('utf-8')

                    # XXX: Maybe str.strip() instead of two str.replace()?
                    line = line.replace('\n', '').replace('\r', '').\
                        split(settings.separator)

                else:
                    line = [
                        x.replace('\n', '').replace('\r', '')
                        if hasattr(x, 'replace') else x for x in line
                    ]

                # in case line has less fields than needed
                if len(line) < maxCol:
                    self.ownlog.msg(2, ('Line #%u has less than %u fields,'
                                        ' skipping! :(\n' % (lnum, maxCol)),
                                    'ERROR')
                    readError = 1
                    continue

                else:

                    # applying filters:
                    if self.filters(line, settings.positiveFilters,
                                    settings.negativeFilters):
                        lFiltered += 1
                        continue

                    # reading names and attributes:
                    if isDir and not isinstance(isDir, tuple):
                        thisEdgeDir = True

                    else:
                        thisEdgeDir = self.process_direction(line, dirCol,
                                                             dirVal, dirSep)

                    refs = []

                    if refCol is not None:
                        refs = common.delEmpty(
                            list(set(line[refCol].split(refSep))))

                    refs = dataio.only_pmids([r.strip() for r in refs])

                    if len(refs) == 0 and settings.must_have_references:
                        rFiltered += 1
                        continue

                    # to give an easy way:
                    if isinstance(settings.ncbiTaxId, int):
                        taxA = settings.ncbiTaxId
                        taxB = settings.ncbiTaxId

                    # to enable more sophisticated inputs:
                    elif isinstance(settings.ncbiTaxId, dict):

                        taxx = self.get_taxon(settings.ncbiTaxId, line)

                        if isinstance(taxx, tuple):
                            taxA = taxx[0]
                            taxB = taxx[1]

                        else:
                            taxA = taxB = taxx

                        taxdA = (
                            settings.ncbiTaxId['A']
                            if 'A' in settings.ncbiTaxId else
                            settings.ncbiTaxId)
                        taxdB = (
                            settings.ncbiTaxId['B']
                            if 'B' in settings.ncbiTaxId else
                            settings.ncbiTaxId)

                        if (('include' in taxdA and
                            taxA not in taxdA['include']) or
                            ('include' in taxdB and
                            taxB not in taxdB['include']) or
                            ('exclude' in taxdA and
                            taxA in taxdA['exclude']) or
                            ('exclude' in taxdB and
                            taxB in taxdB['exclude'])):

                            tFiltered += 1
                            continue

                    else:
                        taxA = taxB = self.ncbi_tax_id

                    if taxA is None or taxB is None:
                        tFiltered += 1
                        continue

                    stim = False
                    inh = False

                    if isinstance(sign, tuple):
                        stim, inh = self.process_sign(line[sign[0]], sign)

                    resource = (
                        [line[settings.resource]]
                        if type(settings.resource) is int else
                        line[settings.resource[0]].split(settings.resource[1])
                        if type(settings.resource) is tuple else
                        [settings.resource]
                    )
                    newEdge = {
                        "nameA": line[settings.nameColA].strip(),
                        "nameB": line[settings.nameColB].strip(),
                        "nameTypeA": settings.nameTypeA,
                        "nameTypeB": settings.nameTypeB,
                        "typeA": settings.typeA,
                        "typeB": settings.typeB,
                        "source": resource,
                        "isDirected": thisEdgeDir,
                        "references": refs,
                        "stim": stim,
                        "inh": inh,
                        "taxA": taxA,
                        "taxB": taxB,
                        "type": settings.intType
                    }
                    # except:
                    # self.ownlog.msg(2,("""Wrong name column indexes (%u and %u),
                    # or wrong separator (%s)? Line #%u\n"""
                    #% (
                    #settings.nameColA, settings.nameColB,
                    # settings.separator, lnum)), 'ERROR')
                    #readError = 1
                    # break
                    # getting additional edge and node attributes
                    attrsEdge = self.get_attrs(line, settings.extraEdgeAttrs,
                                               lnum)
                    attrsNodeA = self.get_attrs(line, settings.extraNodeAttrsA,
                                                lnum)
                    attrsNodeB = self.get_attrs(line, settings.extraNodeAttrsB,
                                                lnum)
                    # merging dictionaries
                    nodeAttrs = {
                        "attrsNodeA": attrsNodeA,
                        "attrsNodeB": attrsNodeB,
                        "attrsEdge": attrsEdge
                    }
                    newEdge.update(nodeAttrs)

                if readError != 0:
                    self.ownlog.msg(2, (
                        'Errors occured, certain lines skipped.'
                        'Trying to read the remaining.\n'), 'ERROR')
                    readError = 1

                edgeList.append(newEdge)

            if hasattr(infile, 'close'):
                infile.close()

            ### !!!! ##
            edgeListMapped = self.map_list(edgeList)
            self.ownlog.msg(
                2, "%u lines have been read from %s,"
                "%u links after mapping; \n\t\t"
                "%u lines filtered by filters;\n\t\t"
                "%u lines filtered because lack of references;\n\t\t"
                "%u lines filtered by taxon filters." %
                (lnum - 1, settings.inFile, len(edgeListMapped), lFiltered,
                 rFiltered, tFiltered))

            if reread or redownload:
                pickle.dump(edgeListMapped, open(edges_cache, 'wb'), -1)
                self.ownlog.msg(2,
                                'Mapped edge list saved to %s' % edges_cache)
        if keep_raw:
            self.data[settings.name] = edgeListMapped

        self.raw_data = edgeListMapped

    def load_list(self, lst, name): # XXX: Not used anywhere
        """
        Loads a custom list to the object's node data lists. See
        :py:attr:`pypath.main.PyPath.lists` attribute for more
        information.

        :arg list lst:
            The list containing the node names [str] from the given
            category (*name*).
        :arg str name:
            The category or identifier for the list of nodes provided.
        """

        self.lists[name] = lst

    def receptors_list(self):
        """
        Loads the Human Plasma Membrane Receptome as a list. This
        resource is human only.
        """

        self.lists['rec'] = common.uniqList(common.flatList([
            self.mapper.map_name(rec, 'genesymbol', 'uniprot',
                                 ncbi_tax_id = 9606)
            for rec in dataio.get_hpmr()]))

    def druggability_list(self):
        """
        Loads the list of druggable proteins from DgiDB. This resource
        is human only.
        """

        self.lists['dgb'] = common.uniqList(common.flatList([
            self.mapper.map_name(dgb, 'genesymbol', 'uniprot', 9606)
            for dgb in dataio.get_dgidb()]))

    def kinases_list(self):
        """
        Loads the list of all known kinases in the proteome from
        kinase.com. This resource is human only.
        """

        self.lists['kin'] = common.uniqList(common.flatList([
            self.mapper.map_name(kin, 'genesymbol', 'uniprot', 9606)
            for kin in dataio.get_kinases()]))

    def tfs_list(self):
        """
        Loads the list of all known transcription factors from TF census
        (Vaquerizas 2009). This resource is human only.
        """

        tfs = dataio.get_tfcensus()
        utfs = [self.mapper.map_name(tf, 'ensg', 'uniprot', 9606)
                for tf in tfs['ensg']]
        utfs += [self.mapper.map_name(h, 'hgnc', 'uniprot', 9606)
                 for h in tfs['hgnc']]

        self.lists['tf'] = common.uniqList(common.flatList(utfs))

    def disease_genes_list(self, dataset='curated'):
        """
        Loads the list of all disease related genes from DisGeNet. This
        resource is human only.
        """

        diss = dataio.get_disgenet(dataset=dataset)
        dis = []

        for di in diss:
            dis.extend(self.mapper.map_name(di['entrez'], 'entrez', 'uniprot',
                                            9606))

        self.lists['dis'] = common.uniqList(dis)

    def signaling_proteins_list(self):
        """
        Compiles a list of signaling proteins (as opposed to other
        proteins like metabolic enzymes, matrix proteins, etc), by
        looking up a few simple keywords in short description of GO
        terms.
        """

        goq = dataio.get_go_quick()

        gosig = set([])

        for term, name in iteritems(goq['names']):

            if 'signal' in name or 'regulat' in name:
                gosig.add(term)

        upsig = set([])

        if 'proteome' not in self.lists:
            self.proteome_list()

        for up, term in iteritems(goq['terms']['P']):

            if len(term & gosig):
                upsig.add(up)

        spsig = set([])

        for u in upsig:
            spsig.update(set(self.mapper.map_name(
                u, 'uniprot', 'uniprot', ncbi_tax_id = self.ncbi_tax_id)))

        upsig = spsig & set(self.lists['proteome'])

        self.lists['sig'] = list(upsig)

    def proteome_list(self, swissprot=True):
        """
        Loads the whole proteome as a list.

        :arg bool swissprot:
            Optional, ``True`` by default. Determines whether to use
            also the information from SwissProt.
        """

        swissprot = 'yes' if swissprot else None
        self.lists['proteome'] = dataio.all_uniprots(self.ncbi_tax_id,
                                                     swissprot=swissprot)

    def cancer_gene_census_list(self):
        """
        Loads the list of cancer driver proteins from the COSMIC Cancer
        Gene Census.
        """

        self.read_list_file(data_formats.cgc)

    def intogen_cancer_drivers_list(self, intogen_file):
        """
        Loads the list of cancer driver proteins from IntOGen data.

        :arg str intogen_file:
            Path to the data file. Can also be [function] that provides
            the data. In general, anything accepted by
            :py:attr:`pypath.input_formats.ReadSettings.inFile`.
        """

        data_formats.intogen_cancer.inFile = intogen_file
        self.read_list_file(data_formats.intogen_cancer)

    def cancer_drivers_list(self, intogen_file=None):
        """
        Loads the list of cancer drivers. Contains information from
        COSMIC (needs user log in credentials) and IntOGen (if provided)
        and adds the attribute to the undirected network nodes.

        :arg str intogen_file:
            Optional, ``None`` by default. Path to the data file. Can
            also be [function] that provides the data. In general,
            anything accepted by
            :py:attr:`pypath.input_formats.ReadSettings.inFile`.
        """

        self.cancer_gene_census_list()

        if intogen_file is not None:
            self.intogen_cancer_drivers_list(intogen_file=intogen_file)
            self.lists['cdv'] = list(
                set(self.lists['cgc']) | set(self.lists['IntOGen']))

        else:
            self.lists['cdv'] = self.lists['cgc']

        self.graph.vs['cdv'] = list(
            map(lambda v: True if v['name'] in self.lists['cdv'] else False,
                self.graph.vs))

    def coverage(self, lst):
        """
        Computes the coverage (range [0, 1]) of a list of nodes against
        the current (undirected) network.

        :arg set lst:
            Can also be [list] (will be converted to [set]) or [str]. In
            the latter case it will retrieve the list with that name (if
            such list exists in :py:attr:`pypath.main.PyPath.lists`).
        """

        lst = lst if isinstance(lst, set) \
            else set(lst) if isinstance(lst, list) \
            else set(self.lists[lst]) \
            if isinstance(lst, str) and lst in self.lists \
            else set([])

        return len(set(self.graph.vs['name']) & lst) / float(len(lst))

    def fisher_enrichment(self, lst, attr, ref='proteome'):
        """
        Computes an enrichment analysis using Fisher's exact test. The
        contingency table is built as follows:
        First row contains the number of nodes in the *ref* list (such
        list is considered to be loaded in
        :py:attr:`pypath.main.PyPath.lists`) and the number of nodes in
        the current (undirected) network. Second row contains the number
        of nodes in *lst* list (also considered to be already loaded)
        and the number of nodes in the network with a non-empty
        attribute *attr*. Uses :py:func:`scipy.stats.fisher_exact`, see
        the documentation of the corresponding package for more
        information.

        :arg str lst:
            Name of the list in :py:attr:`pypath.main.PyPath.lists`
            whose number of elements will be the first element in the
            second row of the contingency table.
        :arg str attr:
            The node attribute name for which the number of nodes in the
            network with such attribute will be the second element of
            the second row of the contingency table.
        :arg str ref:
            Optional, ``'proteome'`` by default. The name of the list in
            :py:attr:`pypath.main.PyPath.lists` whose number of elements
            will be the first element of the first row of the
            contingency table.

        :return:
            * (*float*) -- Prior odds ratio.
            * (*float*) -- P-value or probability of obtaining a
              distribution as extreme as the observed, assuming that the
              null hypothesis is true.
        """

        cont = np.array([[len(self.lists[ref]), self.graph.vcount()],
                         [len(self.lists[lst]),
                          len([1 for v in self.graph.vs if len(v[attr]) > 0])]])

        return stats.fisher_exact(cont)

    def read_list_file(self, settings, **kwargs):
        """
        Reads a list from a file and adds it to
        :py:attr:`pypath.main.PyPath.lists`.

        :arg pypath.input_formats.ReadList settings:
            :py:class:`python.data_formats.ReadList` instance specifying
            the settings of the file to be read. See the class
            documentation for more details.
        :arg \*\*kwargs:
            Extra arguments passed to the file reading function. Such
            function name is outlined in the
            :py:attr:`python.data_formats.ReadList.inFile` attribute and
            defined in :py:mod:`pypath.dataio`.
        """

        _input = None

        if settings.__class__.__name__ != "ReadList":
            self.ownlog.msg(2,
                            ("""No proper input file definition!\n\'settings\'
                             should be a \'readList\' instance\n"""), 'ERROR')
            return None

        if hasattr(dataio, settings.inFile):
            toCall = getattr(dataio, settings.inFile)
            _input = toCall(**kwargs)

        elif not os.path.isfile(settings.inFile):
            self.ownlog.msg(2, "%s: No such file! :(\n" % (settings.inFile),
                            'ERROR')
            return None

        else:
            _input = settings.inFile

        originalNameType = settings.nameType
        defaultNameType = self.default_name_type[settings.typ]
        mapTbl = ''.join([originalNameType, "_", defaultNameType])

        if type(_input) in common.charTypes and os.path.isfile(_input):
            _input = codecs.open(_input, encoding='utf-8', mode='r')

        if _input is None:
            self.ownlog.msg(2, ("""Could not find '\
                                'file or dataio function.\n"""), 'ERROR')
            return None

        self.ownlog.msg(2, "%s opened..." % settings.inFile)
        # finding the largest referred column number,
        # to avoid references out of index
        maxCol = max([settings.nameCol, self.get_max(settings.extraAttrs)])
        # iterating lines from input file
        lnum = 1
        readError = 0
        itemList = []

        for line in _input: # XXX: Could use enumerate(_input) instead of lnum

            if len(line) == 0 or (lnum == 1 and settings.header):
                # empty lines
                # or header row
                lnum += 1
                continue

            if type(line) in common.charTypes:
                line = line.rstrip().split(settings.separator)

            # in case line has less fields than needed
            if len(line) < maxCol:
                self.ownlog.msg(2, ("Line #%u has less than %u fields! :(\n" %
                                    (lnum, maxCol)), 'ERROR')
                readError = 1
                break

            else:

                # reading names and attributes
                try:
                    newItem = {"name": line[settings.nameCol],
                               "nameType": settings.nameType,
                               "type": settings.typ,
                               "source": settings.name}

                except:
                    self.ownlog.msg(2,
                                    ("""Wrong name column indexes (%u and %u),
                                     or wrong separator (%s)? Line #%u\n""" %
                                     (settings.nameCol, settings.separator,
                                      lnum)), 'ERROR')
                    readError = 1
                    break

                # getting additional attributes
                attrsItem = self.get_attrs(line, settings.extraAttrs, lnum)
                # merging dictionaries
                newItem.update(attrsItem)

            if readError != 0:
                break

            itemList.append(newItem)
            lnum += 1

        if hasattr(_input, 'close'):
            _input.close()

        itemListMapped = self.map_list(itemList, singleList=True)
        itemListMapped = list(set(itemListMapped))
        self.ownlog.msg(2, "%u lines have been read from %s, %u '\
                        items after mapping" %
                        (lnum, settings.inFile, len(itemListMapped)))
        self.lists[settings.name] = itemListMapped

    def map_list(self, lst, singleList=False):
        """
        Maps the names from a list of edges or items (molecules).

        :arg list lst:
            List of items or edge dictionaries whose names have to be
            mapped.
        :arg bool singleList:
            Optional, ``False`` by default. Determines whether the
            provided elements are items or edges. This is, either calls
            :py:meth:`pypath.main.PyPath.map_edge` or
            :py:meth:`pypath.main.PyPath.map_item` to map the item
            names.

        :return:
            (*list*) -- Copy of *lst* with their elements' names mapped.
        """

        listMapped = []

        if singleList:

            for item in lst:
                listMapped += self.map_item(item)

        else:

            for edge in lst:
                listMapped += self.map_edge(edge)

        return listMapped

    def map_item(self, item):
        """
        Translates the name in *item* representing a molecule. Default
        name types are defined in
        :py:attr:`pypath.main.PyPath.default_name_type` If the mapping
        is unsuccessful, the item will be added to
        :py:attr:`pypath.main.PyPath.unmapped` list.

        :arg dict item:
            Item whose name is to be mapped to a default name type.

        :return:
            (*list*) -- The default mapped name(s) [str] of *item*.
        """

        # TODO: include
        defaultNames = self.mapper.map_name(
            item['name'], item['nameType'],
            self.default_name_type[item['type']])

        if len(defaultNames) == 0:
            self.unmapped.append(item['name'])

        return defaultNames

    def map_edge(self, edge):
        """
        Translates the name in *edge* representing an edge. Default
        name types are defined in
        :py:attr:`pypath.main.PyPath.default_name_type` If the mapping
        is unsuccessful, the item will be added to
        :py:attr:`pypath.main.PyPath.unmapped` list.

        :arg dict edge:
            Item whose name is to be mapped to a default name type.

        :return:
            (*list*) -- Contains the edge(s) [dict] with default mapped
            names.
        """

        edgeStack = []

        defNameA = self.mapper.map_name(edge['nameA'], edge['nameTypeA'],
                                        self.default_name_type[edge['typeA']],
                                        ncbi_tax_id = edge['taxA'])
        # print 'mapped %s to %s' % (str(edge['nameA']), str(defaultNameA))

        defNameB = self.mapper.map_name(edge['nameB'], edge['nameTypeB'],
                                        self.default_name_type[edge['typeB']],
                                        ncbi_tax_id = edge['taxB'])
        # print 'mapped %s to %s' % (str(edge['nameB']), str(defaultNameB))

        # this is needed because the possibility ambigous mapping
        # one name can be mapped to multiple ones
        # this multiplies the nodes and edges
        # in case of proteins this does not happen too often

        # XXX: I refactored this into a single for loop (using itertools.product):

        #for dnA in defNameA:
            #for dnB in defNameB:
        for dnA, dnB in itertools.product(defNameA, defNameB):
            edge['defaultNameA'] = dnA
            edge['defaultNameTypeA'] = self.default_name_type[edge['typeA']]

            edge['defaultNameB'] = dnB
            edge['defaultNameTypeB'] = self.default_name_type[edge['typeB']]
            edgeStack.append(edge)
            # print 'new edge: %s' % str(edge)

        return edgeStack

    def combine_attr(self, lst, num_method=max):
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

        :arg list lst:
            List of one or two attribute values.
        :arg function num_method:
            Optional, ``max`` by default. Method to merge numeric
            attributes.
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
        if len(lst) > 2:
            lst = [lst[0], self.combine_attr(lst[1:], num_method=num_method)]

        # quick and simple cases:
        if len(lst) == 0:
            return None

        if len(lst) == 1:
            return lst[0]

        if lst[0] == lst[1]:
            return lst[0]

        if lst[0] is None:
            return lst[1]

        if lst[1] is None:
            return lst[0]

        # merge numeric values
        if type(lst[0]) in common.numTypes and type(lst[1]) in common.numTypes:
            return num_method(lst)

        # in case one is list other is set
        lst[0], lst[1] = list_or_set(lst[0], lst[1])

        # merge lists:
        if isinstance(lst[0], list) and isinstance(lst[1], list):

            try:
                # lists of hashable elements only:
                return list(set(itertools.chain(lst[0], lst[1])))

            except TypeError:
                # if contain non-hashable elements:
                return list(itertools.chain(lst[0], lst[1]))

        # merge sets:
        if isinstance(lst[0], set):
            return common.addToSet(lst[0], lst[1])

        if isinstance(lst[1], set):
            return common.addToSet(lst[1], lst[0])

        # merge dicts:
        if isinstance(lst[0], dict) and isinstance(lst[1], dict):
            return common.merge_dicts(lst[0], lst[1])

        # 2 different strings: return a set with both of them
        if ((isinstance(lst[0], str) or isinstance(lst[0], unicode))
            and (isinstance(lst[1], str) or isinstance(lst[1], unicode))):

            if len(lst[0]) == 0:
                return lst[1]

            if len(lst[1]) == 0:
                return lst[0]

            return set([lst[0], lst[1]])

        # one attr is list, the other is simple value:
        if (isinstance(lst[0], list) and type(lst[1]) in common.simpleTypes):

            if lst[1] in common.numTypes or len(lst[1]) > 0:
                return common.addToList(lst[0], lst[1])

            else:
                return lst[0]

        if (isinstance(lst[1], list) and type(lst[0]) in common.simpleTypes):

            if lst[0] in common.numTypes or len(lst[0]) > 0:
                return common.addToList(lst[1], lst[0])

            else:
                return lst[1]

        # special: merging directions
        if (lst[0].__class__.__name__ == 'Direction'
            and lst[1].__class__.__name__ == 'Direction'):

            lst[0].merge(lst[1])
            return lst[0]

        # in case the objects have `__add__()` method:
        if hasattr(lst[0], '__add__'):
            return lst[0] + lst[1]

    def uniq_node_list(self, lst):
        """
        Returns a given list of nodes containing only the unique
        elements.

        :arg list lst:
            List of nodes.

        :return:
            (*list*) -- Copy of *lst* containing only unique nodes.
        """

        uniqLst = {}

        for n in lst:

            if n[0] not in uniqLst:
                uniqLst[n[0]] = n[1]

            else:
                uniqLst[n[0]] = self.merge_attrs(uniqLst[n[0]], n[1])

        return uniqLst

    def collapse_by_name(self, graph=None):
        """
        Collapses nodes with the same name by copying and merging
        all edges and attributes. Operates directly on the provided
        network object.

        :arg igraph.Graph graph:
            Optional, ``None`` by default. The network for which the
            nodes are to be collapsed. If none is provided, takes
            :py:attr:`pypath.main.PyPath.graph` (undirected network) by
            default.
        """

        graph = self.graph if graph is None else graph

        dupli = Counter(graph.vs['name'])

        for name, count in iteritems(dupli):

            if count > 1:
                nodes = graph.vs.select(name = name)

                # the number of nodes might have changed
                if len(nodes) > 1:
                    self.merge_nodes(nodes)

    def merge_nodes(self, nodes, primary=None, graph=None):
        """
        Merges all attributes and edges of selected nodes and assigns
        them to the primary node (by default the one with lowest index).

        :arg list nodes:
            List of node indexes [int] that are to be collapsed.
        :arg int primary:
            Optional, ``None`` by default. ID of the primary edge, if
            none is passed, the node with lowest index on *nodes* is
            selected.
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The network graph object from
            which the nodes are to be merged. If none is passed, takes
            the undirected network graph.
        """

        graph = self.graph if graph is None else graph
        nodes = sorted(list(map(lambda n:
                            n.index if type(n) is not int else n, nodes)))
        nodes = sorted(nodes) # XXX: Isn't it already sorted just above?
        primary = nodes[0] if primary is None else primary
        primary = primary.index if type(primary) is not int else primary
        nonprimary = list(filter(lambda n: n != primary, nodes))
        graph.vs['id_merge'] = list(range(graph.vcount()))

        # combining vertex attributes:
        vprim = graph.vs[primary]

        for attr in vprim.attributes():

            if attr != 'name':
                vprim[attr] = self.combine_attr(list(map(
                                lambda vid: graph.vs[vid][attr],
                                # combining from all nodes
                                nodes)))

        # moving edges of non primary vertices to the primary one
        self.copy_edges(nonprimary, primary, move = True, graph = graph)

        # deleting non primary vertices:
        toDel = list(map(lambda i: graph.vs.select(id_merge=i)[0].index,
                         nonprimary))

        graph.delete_vertices(toDel)
        del graph.vs['id_merge']

    def copy_edges(self, sources, target, move=False, graph=None):
        """
        Copies edges from *sources* node(s) to another one (*target*),
        keeping attributes and directions.

        :arg list sources:
            Contains the vertex index(es) [int] of the node(s) to be
            copied or moved.
        :arg int target:
            Vertex index where edges and attributes are to be copied to.
        :arg bool move:
            Optional, ``False`` by default. Whether to perform copy or
            move (remove or keep the source edges).
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The network graph object from
            which the nodes are to be merged. If none is passed, takes
            the undirected network graph.
        """

        toDel = set([])
        graph = self.graph if graph is None else graph
        graph.vs['id_old'] = list(range(graph.vcount()))
        graph.es['id_old'] = list(range(graph.ecount()))

        # preserve a permanent marker of the target vertex
        ovidt = graph.vs[target]['id_old']

        # collecting the edges of all source vertices into dict
        ses = dict(map(lambda s: (
                        # id_old of source vertices:
                        s,
                        # edges of current source node:
                        set(map(lambda e: e.index,
                                itertools.chain(graph.es.select(_source=s),
                                                graph.es.select(_target=s))))),
                       sources))

        # collecting edges to be newly created
        toAdd = set([])

        for s, es in iteritems(ses):

            for eid in es:
                # the source edge:
                e = graph.es[eid]
                # looking up if target edge already exists:
                vid1 = target if e.source == s else e.source
                vid2 = target if e.target == s else e.target
                te = graph.get_eid(vid1, vid2, error = False)

                if te == -1:
                    # target edge not found, needs to be added:
                    toAdd.add((vid1, vid2))

        # creating new edges
        graph.add_edges(toAdd)

        # copying attributes:
        for ovids, es in iteritems(ses):

            for oeid in es:
                # this is the index of the current source node:
                s = graph.vs.select(id_old = ovids)[0].index
                # this is the index of the current target node:
                t = graph.vs.select(id_old = ovidt)[0].index
                # this is the current source edge:
                e = graph.es.select(id_old = oeid)[0]
                # looking up target edge and peer vertex:
                vid1 = t if e.source == s else e.source
                vid2 = t if e.target == s else e.target
                vid_peer = e.source if e.target == s else e.target
                te = graph.es[graph.get_eid(vid1, vid2)]

                # old direction:
                d = e['dirs']
                # dict from old names to new ones
                # the peer does no change, only s->t
                ids = {graph.vs[s]['name']: graph.vs[t]['name'],
                       graph.vs[vid_peer]['name']: graph.vs[vid_peer]['name']}

                # copying directions and signs:
                te['dirs'] = (d.translate(ids).merge(te['dirs'])
                              if isinstance(te['dirs'], Direction)
                              else d.translate(ids))

                # copying `refs_by_dir`
                te['refs_by_dir'] = self.translate_refsdir(e['refs_by_dir'],
                                                           ids)

                # copying further attributes:
                for eattr in e.attributes():

                    if eattr != 'dirs' and eattr != 'refs_by_dir':
                        te[eattr] = self.combine_attr([te[eattr], e[eattr]])

                # in case we want to delete old edges:
                toDel.add(e.index)

        if move:
            graph.delete_edges(list(toDel))

        # removing temporary attributes
        del graph.es['id_old']
        del graph.vs['id_old']

    def delete_by_taxon(self, tax):
        """
        Removes the proteins of all organisms which are not given in
        *tax*.

        :arg list tax:
            List of NCBI Taxonomy IDs [int] of the organism(s) that are
            to be kept.
        """
        g = self.graph
        toDel = []

        for v in g.vs:

            if v['ncbi_tax_id'] not in tax:
                toDel.append(v.index)

        g.delete_vertices(toDel)
        self.update_vname()
        self.update_db_dict()

    def delete_unknown(self, tax, typ='protein', defaultNameType=None):
        """
        Removes those items which are not in the list of all default
        IDs of the organisms. By default, it means to remove all protein
        nodes not having a human UniProt ID.

        :arg list tax:
            List of NCBI Taxonomy IDs [int] of the organism(s) of
            interest.
        :arg str typ:
            Optional, ``'protein'`` by default. Determines the molecule
            type. These can be ``'protein'``, ``'drug'``, ``'lncrna'``,
            ``'mirna'`` or any other type defined in
            :py:attr:`pypath.main.PyPath.default_name_type`.
        :arg str defaultNameType:
            Optional, ``None`` by default. The default name type for the
            given molecular species. If none is specified takes it from
            :py:attr:`pypath.main.PyPath.default_name_type` (e.g.: for
            ``'protein'``, default is ``'uniprot'``).
        """

        g = self.graph

        if not defaultNameType:
            defaultNameType = self.default_name_type[typ]

        toDel = []
        reflists = {}
        self.update_vname()

        for t in tax:
            idx = (defaultNameType, typ, t)

            if idx in self.reflists:
                reflists[t] = self.reflists[idx].lst

            else:
                msg = (
                    'Missing reference list for %s (default name type: %s), in taxon %u'
                ) % (idx[1], idx[0], t)
                self.ownlog.msg(2, msg, 'ERROR')
                sys.stdout.write(''.join(['\t', msg, '\n']))

                return False

        sys.stdout.write(' :: Comparing with reference lists...')

        for t in tax:
            nt = g.vs['nameType']
            nt = [i for i, j in enumerate(nt) if j == defaultNameType]
            ty = g.vs['type']
            ty = [i for i, j in enumerate(ty) if j == typ]
            tx = g.vs['ncbi_tax_id']
            tx = [i for i, j in enumerate(tx) if j == t]
            vs = list((set(nt) & set(ty)) & set(tx))
            vn = [g.vs[i]['name'] for i in vs]
            toDelNames = list(set(vn) - set(reflists[t]))
            toDel += [self.nodDct[n] for n in toDelNames]

        g.delete_vertices(toDel)
        sys.stdout.write(' done.\n')

    def clean_graph(self):
        """
        Removes multiple edges, unknown molecules and those from wrong
        taxon. Multiple edges will be combined by
        :py:meth:`pypath.main.PyPath.combine_attr` method.
        Loops will be deleted unless the attribute
        :py:attr:`pypath.main.PyPath.loops` is set to ``True``.
        """

        self.ownlog.msg(1, "Removing duplicate edges...", 'INFO')
        g = self.graph

        if not g.is_simple():
            g.simplify(loops=not self.loops, multiple=True,
                       combine_edges=self.combine_attr)

        self.delete_unmapped()

        ## TODO: multiple taxons ##
        if len(self.reflists) != 0:
            self.delete_by_taxon([self.ncbi_tax_id])
            self.delete_unknown([self.ncbi_tax_id])

        x = g.vs.degree()
        zeroDeg = [i for i, j in enumerate(x) if j == 0]
        g.delete_vertices(zeroDeg)
        self.update_vname()

    ###
    # functions to integrate new data into the main igraph network object
    ###

    def count_sol(self): # XXX: Not used anywhere
        """
        Counts the number of nodes with zero degree.

        :return:
            (*int*) -- The number of nodes with zero degree.
        """

        # XXX: Refactor option 1:
        #          return len([1 for i in self.graph.vs.degree() if i == 0])
        #      Refactor option 2:
        #          return Counter(self.graph.vs.degree()[0])

        s = 0

        for i in self.graph.vs.degree():

            if i == 0:
                s += 1

        return s

    def add_update_vertex(self, defAttrs, originalName, originalNameType,
                          extraAttrs={}, add=False):
        """
        Updates the attributes of one node in the (undirected) network.
        Optionally it creates a new node and sets the attributes, but it
        is not efficient as :py:mod:`igraph` needs to reindex vertices
        after this operation, so better to create new nodes in batches.

        :arg dict defAttrs:
            The attribute dictionary of the node to be updated/created.
        :arg str originalName:
            Original node name (e.g.: UniProt ID).
        :arg str originalNameType:
            The original node name type (e.g.: for the previous example,
            this would be ``'uniprot'``).
        :arg dict extraAttrs:
            Optional, ``{}`` by default. Contains any extra attributes
            for the node to be updated.
        :arg bool add:
            Optional, ``False`` by default. If set to ``True`` and the
            node is not in the network, it will be created. Otherwise,
            in such case it will raise an error message.
        """

        g = self.graph

        if not defAttrs["name"] in g.vs["name"]:

            if not add:
                self.ownlog.msg(2, 'Failed to add some vertices', 'ERROR')
                return False

            n = g.vcount()
            g.add_vertices(1)
            g.vs[n]['originalNames'] = {originalName: originalNameType}
            thisNode = g.vs.find(name = defAttrs["name"])

        else:
            thisNode = g.vs.find(name = defAttrs["name"])

            if thisNode["originalNames"] is None:
                thisNode["originalNames"] = {}

            thisNode["originalNames"][originalName] = originalNameType

        for key, value in iteritems(defAttrs):
            thisNode[key] = value

        for key, value in iteritems(extraAttrs):

            if key not in g.vs.attributes():
                g.vs[key] = ([[] for _ in xrange(self.graph.vcount())]
                             if isinstance(value, list) else [None])

            thisNode[key] = self.combine_attr([thisNode[key], value])

    def add_update_edge(self, nameA, nameB, source, isDir, refs, stim, inh,
                        taxA, taxB, typ, extraAttrs={}, add=False):
        """
        Updates the attributes of one edge in the (undirected) network.
        Optionally it creates a new edge and sets the attributes, but it
        is not efficient as :py:mod:`igraph` needs to reindex edges
        after this operation, so better to create new edges in batches.

        :arg str nameA:
            Name of the source node of the edge to be added/updated.
        :arg str nameB:
            Name of the source node of the edge to be added/updated.
        :arg set source:
            Or [list], contains the names [str] of the resources
            supporting that edge.
        :arg bool isDir:
            Whether if the edge is directed or not.
        :arg set refs:
            Or [list], contains the instances of the references
            :py:class:`pypath.refs.Reference` for that edge.
        :arg bool stim:
            Whether the edge is stimulatory or not.
        :arg bool inh:
            Whether the edge is inhibitory or note
        :arg int taxA:
            NCBI Taxonomic identifier of the source molecule.
        :arg int taxB:
            NCBI Taxonomic identifier of the target molecule.
        :arg str typ:
            The type of interaction (e.g.: ``'PPI'``)
        :arg dict extraAttrs:
            Optional, ``{}`` by default. Contains any extra attributes
            for the edge to be updated.
        :arg bool add:
            Optional, ``False`` by default. If set to ``True`` and the
            edge is not in the network, it will be created. Otherwise,
            in such case it will raise an error message.
        """

        g = self.graph

        if not hasattr(self, 'nodDct') or len(self.nodInd) != g.vcount():
            self.update_vname()

        edge = self.edge_exists(nameA, nameB)

        if isinstance(edge, list):

            if not add:
                sys.stdout.write('\tERROR: Failed to add some edges\n')
                self.ownlog.msg(2, 'Failed to add some edges', 'ERROR')
                aid = self.nodDct[nameA]
                bid = self.nodDct[nameB]
                a = g.get_eid(aid, bid, error=False)
                b = g.get_eid(aid, bid, error=False)
                self.failed_edges.append([edge, nameA, nameB, aid, bid, a, b])

                return False

            g.add_edge(edge[0], edge[1])
            edge = self.edge_exists(nameA, nameB)

        # assigning source:
        self.add_set_eattr(edge, 'sources', source)
        # adding references:
        # if len(refs) > 0:
        refs = [_refs.Reference(pmid) for pmid in refs]
        self.add_list_eattr(edge, 'references', refs)
        # updating references-by-source dict:
        sources = source if type(source) in {tuple, set, list} else (source,)

        for src in sources:
            self.add_grouped_set_eattr(edge, 'refs_by_source', src, refs)

        # updating refrences-by-type dict:
        self.add_grouped_set_eattr(edge, 'refs_by_type', typ, refs)

        # setting directions:
        if not g.es[edge]['dirs']:
            g.es[edge]['dirs'] = Direction(nameA, nameB)

        if isDir:
            g.es[edge]['dirs'].set_dir((nameA, nameB), source)
            # updating references-by-direction dict:
            self.add_grouped_set_eattr(edge, 'refs_by_dir', (nameA, nameB),
                                       refs)
        else:
            g.es[edge]['dirs'].set_dir('undirected', source)
            self.add_grouped_set_eattr(edge, 'refs_by_dir', 'undirected', refs)

        # setting signs:
        if stim:
            g.es[edge]['dirs'].set_sign((nameA, nameB), 'positive', source)

        if inh:
            g.es[edge]['dirs'].set_sign((nameA, nameB), 'negative', source)

        # updating sources-by-type dict:
        self.add_grouped_set_eattr(edge, 'sources_by_type', typ, source)
        # adding type:
        self.add_list_eattr(edge, 'type', typ)

        # adding extra attributes:
        for key, value in iteritems(extraAttrs):

            if key not in g.es.attributes():
                g.es[key] = ([[] for _ in xrange(self.graph.ecount())]
                             if isinstance(value, list) else [None])

            g.es[edge][key] = self.combine_attr([g.es[edge][key], value])
###############################################################################
    def add_list_eattr(self, edge, attr, value):
        """

        """

        value = value if isinstance(value, list) else [value]
        e = self.graph.es[edge]

        if attr not in self.graph.es.attributes():
            self.graph.es[attr] = [[] for _ in xrange(0, self.graph.ecount())]

        if e[attr] is None:
            e[attr] = []

        elif not isinstance(e[attr], list):
            e[attr] = [e[attr]]

        e[attr] = common.uniqList(e[attr] + value)

    def add_set_eattr(self, edge, attr, value):
        """

        """

        value = (value if isinstance(value, set) else set(value)
                 if isinstance(value, list) else set([value]))
        e = self.graph.es[edge]

        if attr not in self.graph.es.attributes():
            self.graph.es[attr] = [set([]) for _ in
                                   xrange(0, self.graph.ecount())]
        if e[attr] is None:
            e[attr] = set([])

        elif not isinstance(e[attr], set):
            e[attr] = set(e[attr]) if isinstance(e[attr],
                                                 list) else set([e[attr]])

        e[attr].update(value)

    def add_grouped_eattr(self, edge, attr, group, value):
        """

        """

        value = value if isinstance(value, list) else [value]
        e = self.graph.es[edge]

        if attr not in self.graph.es.attributes():
            self.graph.es[attr] = [{} for _ in xrange(0, self.graph.ecount())]

        if not isinstance(e[attr], dict):
            e[attr] = {}

        if group not in e[attr] or isinstance(e[attr][group], type(None)):
            e[attr][group] = []

        elif not isinstance(e[attr][group], list):
            e[attr][group] = [e[attr][group]]

        e[attr][group] = common.uniqList(e[attr][group] + value)

    def add_grouped_set_eattr(self, edge, attr, group, value):
        """

        """

        value = (value if isinstance(value, set) else set(value)
                 if isinstance(value, list) else set([value]))
        e = self.graph.es[edge]

        if attr not in self.graph.es.attributes():
            self.graph.es[attr] = [{} for _ in xrange(0, self.graph.ecount())]

        if not isinstance(e[attr], dict):
            e[attr] = {}

        if group not in e[attr] or isinstance(e[attr][group], type(None)):
            e[attr][group] = set([])

        elif not isinstance(e[attr][group], set):
            e[attr][group] = (set(e[attr][group])
                              if isinstance(e[attr][group], list)
                              else set([e[attr][group]]))

        e[attr][group].update(value)

    def get_directed(self, graph=False, conv_edges=False, mutual=False,
                     ret=False):
        """
        Converts ``graph`` undirected ``igraph.Graph`` object to a directed one.
        By default it converts the graph in ``PyPath.graph`` and places the directed
        instance in ``PyPath.dgraph``.

        @graph : igraph.Graph
            Undirected graph object.

        @conv_edges : bool
            Whether to convert undirected edges (those without explicit
            direction information) to an arbitrary direction edge or
            a pair of opposite edges.
            Otherwise those will be deleted. Default is ``False``.

        @mutual : bool
            If ``conv_edges`` is ``True``, whether to convert the
            undirected edges to a single, arbitrary directed edge,
            or a pair of opposite directed edges. Default is ``False``.

        @ret : bool
            Return the directed graph instance, or return ``None``.
            Default is ``False`` (returns ``None``).
        """

        toDel = []
        g = self.graph if not graph else graph
        d = g.as_directed(mutual=True)
        self.update_vname()
        d.es['directed_sources'] = [[] for _ in xrange(g.ecount())]
        d.es['undirected_sources'] = [[] for _ in xrange(g.ecount())]
        d.es['directed'] = [False for _ in xrange(g.ecount())]
        prg = Progress(
            total=g.ecount(), name="Setting directions", interval=17)

        for e in g.es:
            """
            This works because in directed graphs get_eid() defaults to
            directed = True, so the source -> target edge is returned.
            """
            dir_one = (g.vs['name'][e.source], g.vs['name'][e.target])
            dir_two = (g.vs['name'][e.target], g.vs['name'][e.source])
            dir_edge_one = d.get_eid(
                d.vs['name'].index(g.vs['name'][e.source]),
                d.vs['name'].index(g.vs['name'][e.target]))
            dir_edge_two = d.get_eid(
                d.vs['name'].index(g.vs['name'][e.target]),
                d.vs['name'].index(g.vs['name'][e.source]))

            if not e['dirs'].get_dir(dir_one):

                if not conv_edges or e['dirs'].get_dir(dir_two):
                    toDel.append(dir_edge_one)

            else:
                d.es[dir_edge_one]['directed'] = True
                d.es[dir_edge_one]['directed_sources'] += \
                    e['dirs'].get_dir(dir_one, sources=True)
                d.es[dir_edge_one]['undirected_sources'] += \
                    e['dirs'].get_dir('undirected', sources=True)

            if not e['dirs'].get_dir(dir_two):

                if not conv_edges or e['dirs'].get_dir(dir_one):
                    toDel.append(dir_edge_two)

            else:
                d.es[dir_edge_two]['directed'] = True
                d.es[dir_edge_two]['directed_sources'] += \
                    e['dirs'].get_dir(dir_two, sources=True)
                d.es[dir_edge_two]['undirected_sources'] += \
                    e['dirs'].get_dir('undirected', sources=True)

            if e['dirs'].get_dir('undirected') and \
                not e['dirs'].get_dir(dir_one) and \
                    not e['dirs'].get_dir(dir_two):

                if conv_edges:
                    d.es[dir_edge_one]['undirected_sources'] += \
                        e['dirs'].get_dir('undirected', sources=True)

                    if mutual:
                        d.es[dir_edge_two]['undirected_sources'] += \
                            e['dirs'].get_dir('undirected', sources=True)

                    else:
                        toDel.append(dir_edge_two)

                else:
                    toDel += [dir_edge_one, dir_edge_two]

            prg.step()

        d.delete_edges(list(set(toDel)))
        prg.terminate()
        deg = d.vs.degree()
        toDel = []

        for v in d.vs:

            if deg[v.index] == 0:
                toDel.append(v.index)

        del self.nodInd

        if len(toDel) > 0:
            d.delete_vertices(list(set(toDel)))

        if not graph:
            self.dgraph = d
            self._directed = self.dgraph

        self._get_directed()
        self._get_undirected()
        self.update_vname()

        if graph or ret:
            return d

    def new_edges(self, edges):
        """
        """

        self.graph.add_edges(list(edges))

    def new_nodes(self, nodes):
        """
        """

        self.graph.add_vertices(list(nodes))

    def edge_exists(self, nameA, nameB):
        """
        Returns a tuple of vertice indices if edge doesn't exists,
        otherwise edge id. Not sensitive to direction.
        """

        if not hasattr(self, 'nodDct'):
            self.update_vname()

        nodes = [self.nodDct[nameA], self.nodDct[nameB]]
        edge = self.graph.get_eid(nodes[0], nodes[1], error=False)

        if edge != -1:
            return edge

        else:
            nodes.sort()
            return nodes

    def edge_names(self, e):
        """
        """

        if isinstance(e, int):
            e = self.graph.es[e]

        return (self.graph.vs[e.source]['name'],
                self.graph.vs[e.target]['name'])

    def node_exists(self, name):
        """
        """

        if not hasattr(self, 'nodInd'):
            self.update_vname()

        return name in self.nodInd

    def names2vids(self, names):
        """
        """

        vids = []

        if not hasattr(self, 'nodInd'):
            self.update_vname()

        for n in names:

            if n in self.nodInd:
                vids.append(self.nodDct[n])

        return vids

    def _get_edge(self, nodes):
        """
        Returns the edge id only if there is an edge from nodes[0] to nodes[1],
        returns False if edge exists in opposite direction, or no edge exists
        between the two vertices, or any of the vertice ids doesn't exist.
        To find edges without regarding their direction, see edge_exists().
        """

        g = self.graph

        try:
            e = g.get_eid(nodes[0], nodes[1])
            return e

        except:
            return False

    def straight_between(self, nameA, nameB):
        """
        This does actually the same as get_edge(), but by names
        instead of vertex ids.
        """

        nodNm = sorted([nameA, nameB])
        nodes = [
            self.graph.vs['name'].index(nodNm[0]),
            self.graph.vs['name'].index(nodNm[1])
        ]
        edge = self._get_edge(nodes)

        if isinstance(edge, int):
            return edge

        else:
            return nodes

    def all_between(self, nameA, nameB):
        """
        Returns all edges between two given vertex names. Similar to
        straight_between(), but checks both directions, and returns
        list of edge ids in [undirected, straight, reversed] format,
        for both nameA -> nameB and nameB -> nameA edges.
        """

        g = self.graph
        edges = {'ab': [None, None, None], 'ba': [None, None, None]}
        eid = self.edge_exists(nameA, nameB)

        if isinstance(eid, int):

            if g.es[eid]['dirs'].get_dir('undirected'):
                edges['ab'][0] = eid
                edges['ba'][0] = eid

            if g.es[eid]['dirs'].get_dir((nameA, nameB)):
                edges['ab'][1] = eid
                edges['ba'][2] = eid

            if g.es[eid]['dirs'].get_dir((nameB, nameA)):
                edges['ab'][2] = eid
                edges['ba'][1] = eid

        return edges

    def get_node_pair(self, nameA, nameB, directed = False):
        """
        """

        if not hasattr(self, 'nodDct'):
            self.update_vname()

        g = self._directed if directed else self._undirected
        nodDct = self.dnodDct if directed else self.nodDct
        nodes = [nameA, nameB] if not directed else sorted([nameA, nameB])

        try:
            nodeA = nodDct[nodes[0]]
            nodeB = nodDct[nodes[1]]
            return (nodeA, nodeB)

        except:
            return False

    def update_attrs(self):
        """
        """

        for attr in self.graph.vs.attributes():
            types = list(
                set([type(x) for x in self.graph.vs[attr] if x is not None]))

            if len(types) > 1:
                self.ownlog.msg(2,
                                'Vertex attribute `%s` has multiple types of'
                                ' values: %s' % (attr, ', '.join(
                                    [x.__name__ for x in types])), 'WARNING')

            elif len(types) == 0:
                self.ownlog.msg(2, 'Vertex attribute `%s` has only None values'
                                % (attr), 'WARNING')

            if len(types) > 0:

                if list in types:
                    self.vertexAttrs[attr] = list

                else:
                    self.vertexAttrs[attr] = types[0]

                self.init_vertex_attr(attr)

        for attr in list(set(self.graph.es.attributes()) - set(['dirs'])):
            types = list(
                set([type(x) for x in self.graph.es[attr] if x is not None]))

            if len(types) > 1:
                self.ownlog.msg(2, 'Edge attribute `%s` has multiple types of'
                                ' values: %s' % (attr, ', '.join(
                                    [x.__name__ for x in types])), 'WARNING')

            elif len(types) == 0:
                self.ownlog.msg(2, 'Edge attribute `%s` has only None values' %
                                (attr), 'WARNING')

            if len(types) > 0:

                if set in types:
                    self.edgeAttrs[attr] = set

                elif list in types:
                    self.edgeAttrs[attr] = list

                else:
                    self.edgeAttrs[attr] = types[0]

                self.init_edge_attr(attr)

    def init_vertex_attr(self, attr):
        """
        Fills vertex attribute with its default values, creates
        lists if in `vertexAttrs` the attribute is registered as list.
        """

        for v in self.graph.vs:

            if v[attr] is None:
                v[attr] = self.vertexAttrs[attr]()

            if self.vertexAttrs[attr] is list and type(v[
                    attr]) in common.simpleTypes:
                v[attr] = [v[attr]] if len(v[attr]) > 0 else []

    def init_edge_attr(self, attr):
        """
        Fills edge attribute with its default values, creates
        lists if in `edgeAttrs` the attribute is registered as list.
        """

        for e in self.graph.es:

            if e[attr] is None:
                e[attr] = self.edgeAttrs[attr]()

            if (self.edgeAttrs[attr] is list or
                self.edgeAttrs[attr] is set) and type(e[
                    attr]) in common.simpleTypes:

                e[attr] = [e[attr]] if (
                    type(e[attr]) not in common.charTypes or
                    len(e[attr]) > 0) else []

            if self.edgeAttrs[attr] is set and type(e[attr]) is list:

                e[attr] = set(e[attr])

    def attach_network(self, edgeList=False, regulator=False):
        """
        Adds edges to the network from edgeList obtained from file or
        other input method.
        """

        g = self.graph

        if not edgeList:

            if self.raw_data is not None:
                edgeList = self.raw_data

            else:
                self.ownlog.msg(2, "attach_network(): No data, nothing to do.",
                                'INFO')
                return True

        if isinstance(edgeList, str):

            if edgeList in self.data:
                edgeList = self.data[edgeList]

            else:
                self.ownlog.msg(2, "`%s' looks like a source name, but no data"
                                "available under this name." % (edgeList),
                                'ERROR')
                return False

        nodes = []
        edges = []
        # adding nodes and edges first in bunch,
        # to avoid multiple reindexing by igraph
        self.update_vname()
        prg = Progress(
            total=len(edgeList), name="Processing nodes", interval=50)

        for e in edgeList:
            aexists = self.node_exists(e["defaultNameA"])
            bexists = self.node_exists(e["defaultNameB"])

            if not aexists and (not regulator or bexists):
                nodes.append(e["defaultNameA"])

            if not bexists and not regulator:
                nodes.append(e["defaultNameB"])

            prg.step()

        prg.terminate()
        self.new_nodes(set(nodes))
        self.ownlog.msg(2, 'New nodes have been created', 'INFO')
        self.update_vname()
        prg = Progress(
            total=len(edgeList), name='Processing edges', interval=50)

        for e in edgeList:
            aexists = self.node_exists(e["defaultNameA"])
            bexists = self.node_exists(e["defaultNameB"])

            if aexists and bexists:
                edge = self.edge_exists(e["defaultNameA"], e["defaultNameB"])

                if isinstance(edge, list):
                    edges.append(tuple(edge))

                prg.step()

        prg.terminate()
        self.new_edges(set(edges))
        self.ownlog.msg(2, "New edges have been created", 'INFO')
        self.ownlog.msg(2, ("""Introducing new node and edge attributes..."""),
                        'INFO')
        prg = Progress(
            total=len(edgeList), name="Processing attributes", interval=30)
        nodes_updated = []
        self.update_vname()

        for e in edgeList:
            # adding new node attributes

            if e["defaultNameA"] not in nodes_updated:
                defAttrs = {
                    "name": e["defaultNameA"],
                    "label": e["defaultNameA"],
                    "nameType": e["defaultNameTypeA"],
                    "type": e["typeA"],
                    "ncbi_tax_id": e["taxA"]
                }
                self.add_update_vertex(defAttrs, e["nameA"], e["nameTypeA"],
                                       e["attrsNodeA"])
                nodes_updated.append(e["defaultNameA"])

            if e["defaultNameB"] not in nodes_updated:
                defAttrs = {
                    "name": e["defaultNameB"],
                    "label": e["defaultNameB"],
                    "nameType": e["defaultNameTypeB"],
                    "type": e["typeB"],
                    "ncbi_tax_id": e["taxB"]
                }
                self.add_update_vertex(defAttrs, e["nameB"], e["nameTypeB"],
                                       e["attrsNodeB"])
                nodes_updated.append(e["defaultNameB"])

            # adding new edge attributes
            self.add_update_edge(e["defaultNameA"], e["defaultNameB"],
                                 e["source"], e["isDirected"], e["references"],
                                 e["stim"], e["inh"], e["taxA"], e["taxB"],
                                 e["type"], e["attrsEdge"])
            prg.step()

        prg.terminate()
        self.raw_data = None
        self.update_attrs()

    def apply_list(self, name, node_or_edge="node"):
        """
        Creates vertex or edge attribute based on a list.
        """

        if name not in self.lists:
            self.ownlog.msg(1, ("""No such list: %s""" % name), 'ERROR')
            return None

        g = self.graph

        if node_or_edge == "edge":
            g.es[name] = [None]

        else:
            g.vs[name] = [None]

        if isinstance(self.lists[name], dict):

            if node_or_edge == "edge":

                for e in g.es:

                    if (v[e.source]["name"], v[e.target]["name"]
                        ) in self.lists[name]:
                        e[name] = self.lists[name][(v[e.source]["name"],
                                                    v[e.target]["name"])]

                    if (v[e.target]["name"], v[e.source]["name"]
                        ) in self.lists[name]:
                        e[name] = self.lists[name][(v[e.target]["name"],
                                                    v[e.source]["name"])]

            else:

                for v in g.vs:

                    if v["name"] in self.lists[name]:
                        v[name] = self.lists[name][v["name"]]

        if isinstance(self.lists[name], list):

            if node_or_edge == "edge":

                for e in g.es:

                    if (v[e.source]["name"], v[e.target]["name"]
                        ) in self.lists[name]:
                        e[name] = True

                    else:
                        e[name] = False

            else:

                for v in g.vs:

                    if v["name"] in self.lists[name]:
                        v[name] = True

                    else:
                        v[name] = False

    def merge_lists(self, nameA, nameB, name=None, and_or="and", delete=False,
                    func="max"):
        """
        Merges two lists in `lists`.
        """

        if nameA not in self.lists:
            self.ownlog.msg(1, ("""No such list: %s""" % nameA), 'ERROR')
            return None

        if nameB not in self.lists:
            self.ownlog.msg(1, ("""No such list: %s""" % nameB), 'ERROR')
            return None

        name = '_'.join([nameA, nameB]) if name is None else name

        if isinstance(self.lists[nameA], list) and isinstance(
                self.lists[nameB], list):

            if and_or == "and":
                self.lists[name] = list(
                    set(self.lists[nameA]) | set(self.lists[nameB]))

            if and_or == "or":
                self.lists[name] = list(
                    set(self.lists[nameA]) & set(self.lists[nameB]))

        if isinstance(self.lists[nameA], dict) and isinstance(
                self.lists[nameB], dict):
            self.lists[name] = {}

            if and_or == "and":
                keys = list(
                    set(self.lists[nameA].keys) | set(self.lists[nameB].keys(
                    )))

                for k in keys:

                    if k in self.lists[nameA]:
                        self.lists[name][k] = self.lists[nameA][k]

                    if k in self.lists[nameB]:
                        self.lists[name][k] = self.lists[nameB][k]

                    if k in self.lists[nameA] and k in self.lists[nameB]:
                        self.lists[name][k] = self.combine_attr(
                            [self.lists[nameA][k], self.lists[nameB][k]])

            if and_or == "or":
                keys = list(
                    set(self.lists[nameA].keys) & set(self.lists[nameB].keys(
                    )))

                for k in keys:
                    self.lists[name][k] = self.combine_attr(
                        [self.lists[nameA][k], self.lists[nameB][k]])

        if delete:
            del self.lists[nameA]
            del self.lists[nameB]

    def save_session(self):
        """
        Save current state into pickle dump.
        """

        pickleFile = "pypath-" + self.session + ".pickle"
        self.ownlog.msg(1, ("""Saving session to %s... """ % pickleFile),
                        'INFO')

        with open(pickleFile, "wb") as f:
            pickle.dump(self, f, -1)

    ###
    # functions for plotting // with custom typeface ;)
    ###

    #
    # functions to compare networks and pathways
    #

    def databases_similarity(self, index='simpson'):
        """
        """

        g = self.graph
        edges = {}
        nodes = {}
        self.update_sources()
        nodes = dict([(s, [v.index for v in g.vs if s in v['sources']])
                      for s in self.sources])
        edges = dict([(s, [e.index for e in g.es if s in e['sources']])
                      for s in self.sources])
        sNodes = self.similarity_groups(nodes, index=index)
        sEdges = self.similarity_groups(edges, index=index)
        return {'nodes': sNodes, 'edges': sEdges}

    def similarity_groups(self, groups, index='simpson'):
        """
        """

        index_func = '%s_index' % index
        # these are imported into globals() from common:

        if hasattr(sys.modules['%s.common' % self.__module__.split('.')[0]],
                   index_func):
            to_call = getattr(sys.modules['%s.common' %
                                          self.__module__.split('.')[0]],
                              index_func)
            grs = sorted(groups.keys())
            sor = dict([(g, {}) for g in grs])

            for g1 in xrange(0, len(grs)):

                for g2 in xrange(g1, len(grs)):
                    sor[grs[g1]][grs[g2]] = to_call(groups[grs[g1]],
                                                    groups[grs[g2]])
                    sor[grs[g2]][grs[g1]] = sor[grs[g1]][grs[g2]]

            return sor

        else:
            self.ownlog.msg(2, 'No such function: %s()' % index_func, 'ERROR')

    def sorensen_pathways(self, pwlist=None):
        """
        """

        g = self.graph

        if pwlist is None:
            self.update_pathway_types()
            pwlist = self.pathway_types

        for p in pwlist:

            if p not in g.vs.attributes():
                self.ownlog.msg(2, ("No such vertex attribute: %s" % p),
                                'ERROR')

        edges = {}
        nodes = {}

        for e in g.es:
            indA = e.source
            indB = e.target
            pwsA = []
            pwsB = []

            for p in pwlist:

                if g.vs[indA][p] is not None:

                    for pw in g.vs[indA][p]:
                        thisPw = p.replace("_pathways", "__") + pw

                        if thisPw not in nodes:
                            nodes[thisPw] = []

                        nodes[thisPw].append(indA)
                        pwsA.append(thisPw)

                if g.vs[indB][p] is not None:

                    for pw in g.vs[indB][p]:
                        thisPw = p.replace("_pathways", "__") + pw

                        if thisPw not in nodes:
                            nodes[thisPw] = []

                        nodes[thisPw].append(indB)
                        pwsB.append(thisPw)

            pwsE = set(pwsA).intersection(set(pwsB))

            for pw in pwsE:

                if pw not in edges:
                    edges[pw] = []

                edges[pw].append(e.index)

        sNodes = self.similarity_groups(nodes, index='sorensen')
        sEdges = self.similarity_groups(edges, index='sorensen')

        return {"nodes": sNodes, "edges": sEdges}

    def write_table(self, tbl, outfile, sep="\t", cut=None, colnames=True,
                    rownames=True):
        """
        """

        out = ''
        rn = list(tbl.keys())

        if "header" in rn:
            cn = tbl["header"]
            del tbl["header"]
            rn.remove("header")

        else:
            cn = [str(i) for i in xrange(0, len(tbl[rn[0]]))]

        if colnames:

            if rownames:
                out += sep

            out += sep.join(cn) + "\n"

        for r in rn:

            if rownames:
                out += str(r)[0:cut] + sep

            thisRow = [str(i) for i in tbl[r]]
            out += sep.join(thisRow) + "\n"

        f = open(os.path.join(self.outdir, outfile), 'w')
        f.write(out)
        f.close()

    def search_attr_or(self, obj, lst):
        """
        """

        if len(lst) == 0:
            return True

        for a, v in iteritems(lst):

            if (isinstance(v, list) and
                    len(set(obj[a]).intersection(v)) > 0) or (
                        not isinstance(v, list) and obj[a] == v):
                return True

        return False

    def search_attr_and(self, obj, lst):
        """
        """

        for a, v in iteritems(lst):

            if (isinstance(v, list) and
                    len(set(obj[a]).intersection(v)) == 0) or (
                        not isinstance(v, list) and obj[a] != v):
                return False

        return True

    def get_sub(self, crit, andor="or", graph=None):
        """
        """

        g = self.graph if graph is None else graph
        keepV = []
        delE = []

        if andor == "and":

            for e in g.es:
                keepThis = self.search_attr_and(e, crit["edge"])

                if keepThis:
                    keepA = self.search_attr_and(g.vs[e.source], crit["node"])

                    if keepA:
                        keepV += [e.source]

                    keepB = self.search_attr_and(g.vs[e.target], crit["node"])

                    if keepB:
                        keepV += [e.target]

                    if not keepA or not keepB:
                        delE += [e.index]

                else:
                    delE += [e.index]

        else:

            for e in g.es:
                keepThis = self.search_attr_or(e, crit["edge"])

                if keepThis:
                    keepV += [e.source, e.target]
                    continue

                else:
                    delE += [e.index]

                if len(crit["node"]) > 0:
                    keepA = self.search_attr_or(g.vs[e.source], crit["node"])

                    if keepA:
                        keepV += [e.source]

                    keepB = self.search_attr_or(g.vs[e.target], crit["node"])

                    if keepB:
                        keepV += [e.target]

        return {"nodes": list(set(keepV)), "edges": list(set(delE))}

    def edgeseq_inverse(self, edges):
        """
        """

        g = self.graph
        inv = []

        for e in g.es:

            if e.index not in set(edges):
                inv.append(e.index)

        return inv

    def get_network(self, crit, andor="or", graph=None):
        """
        """

        g = self.graph if graph is None else graph
        sub = self.get_sub(crit, andor=andor, graph=g)
        new = g.copy()
        new.delete_edges(sub["edges"])

        return new.induced_subgraph(sub["nodes"])

    def separate(self):
        """
        Separates networks from different sources.
        Returns dict of igraph objects.
        """

        return dict([(s, self.get_network({
            'edge': {
                'sources': [s]
            },
            'node': {}
        })) for s in self.sources])

    def separate_by_category(self):
        """
        Separate networks based on categories.
        Returns dict of igraph objects.
        """

        cats = \
            dict(
                list(
                    map(
                        lambda c:
                            (
                                c,
                                list(
                                    filter(
                                        lambda s:
                                            s in self.sources,
                                        map(
                                            lambda cs:
                                                cs[0],
                                            filter(
                                                lambda cs:
                                                    cs[1] == c,
                                                iteritems(
                                                    data_formats.categories)
                                                    )
                                        )
                                    )
                                )
                            ),
                        self.has_cats
                    )
                )
            )
        return dict([(c, self.get_network({
            'edge': {
                'sources': s
            },
            'node': {}
        })) for c, s in iteritems(cats)])

    def update_pathway_types(self):
        """
        """

        g = self.graph
        pwTyp = []

        for i in g.vs.attributes():

            if i.find("_pathways") > -1:
                pwTyp.append(i)

        self.pathway_types = pwTyp

    def source_similarity(self, outfile=None):
        """
        """

        if outfile is None:
            outfile = ''.join(["pwnet-", self.session, "-sim-src"])

        res = self.database_similarity(index='sorensen')
        self.write_table(res["nodes"], outfile + "-nodes")
        self.write_table(res["edges"], outfile + "-edges")

    def pathway_similarity(self, outfile=None):
        """
        """

        if outfile is None:
            outfile = ''.join(["pwnet-", self.session, "-sim-pw"])

        res = self.sorensen_pathways()
        self.write_table(res["nodes"], outfile + "-nodes", cut=20)
        self.write_table(res["edges"], outfile + "-edges", cut=20)

    def update_sources(self):
        """
        Makes sure that the `sources` attribute is an up to date
        list of all sources in the current network.
        """

        g = self.graph
        src = []

        for e in g.es:
            src += e["sources"]

        self.sources = list(set(src))
        self.update_cats()

    def update_cats(self):
        """
        Makes sure that the `has_cats` attribute is an up to date
        set of all categories in the current network.
        """

        self.has_cats = set(
            list(
                map(lambda s: data_formats.categories[s],
                    filter(lambda s: s in data_formats.categories,
                           self.sources))))

    def update_pathways(self):
        """
        """

        g = self.graph
        pws = {}

        for v in g.vs:

            for k in g.vs.attributes():

                if k.find('_pathways') > -1:

                    if k not in pws:
                        pws[k] = []

                    if v[k] is not None:

                        for p in v[k]:

                            if p not in set(pws[k]):
                                pws[k].append(p)

        self.pathways = pws

    def delete_unmapped(self):
        """
        """

        if "unmapped" in self.graph.vs["name"]:
            self.graph.delete_vertices(
                self.graph.vs.find(name="unmapped").index)
            self.update_db_dict()
            self.update_vname()

    def genesymbol_labels(self, graph=None, remap_all=False):
        """
        Creats vertex attribute ``label`` and fills up with Gene Symbols
        of all proteins where the Gene Symbol can be looked up based on
        the default name of the protein vertex.
        If the attribute ``label`` had been already initialized,
        updates this attribute or recreates if ``remap_all``
        is ``True``.
        """

        self._already_has_directed()

        if graph is None and self.dgraph is not None:
            self.genesymbol_labels(graph=self.dgraph, remap_all=remap_all)

        g = self.graph if graph is None else graph
        defaultNameType = self.default_name_type["protein"]
        labelNameTypes = {'protein': 'genesymbol',
                          'mirna': 'mir-mat-name'}

        if 'label' not in g.vs.attributes():
            remap_all = True

        labels = [
            None if remap_all or v['label'] == v['name'] else v['label']
            for v in g.vs
        ]

        for v, l, i in zip(g.vs, labels, xrange(g.vcount())):

            if l is None:
                label = []

                if v['type'] in labelNameTypes:
                    label = self.mapper.map_name(v['name'],
                                                self.default_name_type[v['type']],
                                                labelNameTypes[v['type']],
                                                ncbi_tax_id = v['ncbi_tax_id'])

                if not len(label):
                    labels[i] = v['name']

                else:
                    labels[i] = label[0]

        g.vs['label'] = labels

    def network_stats(self, outfile=None):
        """
        Calculates basic statistics for the whole network
        and each of sources. Writes the results in a tab file.
        """

        if outfile is None:
            outfile = '-'.join(["pwnet", self.session, "stats"])

        stats = {}
        stats['header'] = [
            "vnum", "enum", "deg_avg", "diam", "trans", "adh", "coh"
        ]

        for k in xrange(0, len(self.sources) + 1):
            s = "All" if k == len(self.sources) else self.sources[k]
            g = self.graph if k == len(self.sources) else self.get_network({
                "edge": {
                    "sources": [s]
                },
                "node": {}
            })

            if g.vcount() > 0:
                stats[s] = [
                    g.vcount(), g.ecount(),
                    sum(g.vs.degree()) / float(len(g.vs)), g.diameter(),
                    g.transitivity_undirected(), g.adhesion(), g.cohesion()
                ]
        self.write_table(stats, outfile)

    def degree_dists(self):
        """
        """

        dds = {}

        for s in self.sources:
            g = self.get_network({"edge": {"sources": [s]}, "node": {}})

            if g.vcount() > 0:
                dds[s] = g.degree_distribution()

        for k, v in iteritems(dds):
            filename = os.path.join(self.outdir, ''.join(["pwnet-", self.session, "-degdist-", k]))
            bins = []
            vals = []

            for i in v.bins():
                bins.append(int(i[0]))
                vals.append(int(i[2]))

            out = ''.join([
                ';'.join(str(x)
                         for x in bins), "\n", ';'.join(str(x)
                                                        for x in vals), "\n"
            ])
            f = codecs.open(filename, encoding='utf-8', mode='w')
            f.write(out)
            f.close()

    def intergroup_shortest_paths(self, groupA, groupB, random=False):
        """
        """

        self.update_sources()

        if groupA not in self.graph.vs.attributes():
            self.ownlog.msg(2, ("No such attribute: %s" % groupA), 'ERROR')
            return False

        if groupB not in self.graph.vs.attributes():
            self.ownlog.msg(2, ("No such attribute: %s" % groupB), 'ERROR')
            return False

        deg_pathlen = {}
        rat_pathlen = {}
        rand_pathlen = {}
        diam_pathlen = {}

        for k in xrange(0, len(self.sources) + 1):
            s = "All" if k == len(self.sources) else self.sources[k]
            outfile = '-'.join([s, groupA, groupB, "paths"])
            f = self.graph if k == len(self.sources) else self.get_network({
                "edge": {
                    "sources": [s]
                },
                "node": {}
            })
            paths = []
            grA = []
            grB = []

            for v in f.vs:

                if v[groupA]:
                    grA.append(v.index)

                if v[groupB]:
                    grB.append(v.index)

            for v in f.vs:

                if v[groupB]:
                    pt = f.get_shortest_paths(v.index, grA, output="epath")

                    for p in pt:
                        l = len(p)

                        if l > 0:
                            paths.append(l)

                    if (v.index in grA):
                        paths.append(0)

            self.write_table(
                {
                    "paths": paths
                },
                outfile,
                sep=";",
                colnames=False,
                rownames=False)
            deg = f.vs.degree()
            mean_pathlen = sum(paths) / float(len(paths))
            deg_pathlen[s] = [mean_pathlen, sum(deg) / float(len(deg))]
            rat_pathlen[s] = [
                mean_pathlen, f.vcount() / float(len(list(set(grA + grB))))
            ]
            diam_pathlen[s] = [mean_pathlen, f.diameter()]

            if random:
                groupA_random = groupA + "_random"
                groupB_random = groupB + "_random"
                random_pathlen = []

                for i in xrange(0, 100):
                    f.vs[groupA_random] = copy.copy(f.vs[groupA])
                    f.vs[groupB_random] = copy.copy(f.vs[groupB])
                    random.shuffle(f.vs[groupA_random])
                    random.shuffle(f.vs[groupB_random])
                    paths = []
                    grA = []
                    grB = []

                    for v in f.vs:

                        if v[groupA_random]:
                            grA.append(v.index)

                        if v[groupB_random]:
                            grB.append(v.index)

                    for v in f.vs:

                        if v[groupB_random]:
                            pt = f.get_shortest_paths(
                                v.index, grA, output="epath")

                            for p in pt:
                                l = len(p)

                                if l > 0:
                                    paths.append(l)

                            if (v.index in grA):
                                paths.append(0)

                    if len(paths) > 0:
                        random_pathlen.append(sum(paths) / float(len(paths)))

                if len(random_pathlen) > 0:
                    rand_pathlen[s] = [
                        mean_pathlen,
                        sum(random_pathlen) / float(len(random_pathlen))
                    ]

                else:
                    rand_pathlen[s] = [mean_pathlen, 0.0]

        deg_pathlen["header"] = ["path_len", "degree"]
        self.write_table(deg_pathlen, "deg_pathlen", sep=";")
        diam_pathlen["header"] = ["path_len", "diam"]
        self.write_table(diam_pathlen, "diam_pathlen", sep=";")
        rat_pathlen["header"] = ["path_len", "ratio"]
        self.write_table(rat_pathlen, "rat_pathlen", sep=";")

        if random:
            rand_pathlen["header"] = ["path_len", "random"]
            self.write_table(rand_pathlen, "rand_pathlen", sep=";")

    def update_vertex_sources(self):
        """
        """

        g = self.graph

        for attr in ['sources', 'references']:
            g.vs[attr] = [set([]) for _ in g.vs]

            for e in g.es:
                g.vs[e.source][attr].update(e[attr])
                g.vs[e.target][attr].update(e[attr])

    def set_categories(self):
        """
        """

        self.graph.vs['cat'] = [set([]) for _ in self.graph.vs]
        self.graph.es['cat'] = [set([]) for _ in self.graph.es]
        self.graph.es['refs_by_cat'] = [{} for _ in self.graph.es]

        for v in self.graph.vs:

            for s in v['sources']:

                if s in data_formats.categories:
                    v['cat'].add(data_formats.categories[s])

        for e in self.graph.es:

            for s in e['sources']:

                if s in data_formats.categories:
                    cat = data_formats.categories[s]
                    e['cat'].add(cat)

                    if cat not in e['refs_by_cat']:
                        e['refs_by_cat'][cat] = set([])

                    if s in e['refs_by_source']:
                        e['refs_by_cat'][cat].update(e['refs_by_source'][s])

    def basic_stats_intergroup(self, groupA, groupB, header=None):
        """
        """

        result = {}
        g = self.graph

        for k in xrange(0, len(self.sources) + 1):
            s = "All" if k == len(self.sources) else self.sources[k]
            f = self.graph if k == len(self.sources) else self.get_network({
                "edge": {
                    "sources": set([s])
                },
                "node": {}
            })
            deg = f.vs.degree()
            bw = f.vs.betweenness()
            vnum = f.vcount()
            enum = f.ecount()
            cancerg = 0
            drugt = 0
            cdeg = []
            tdeg = []
            ddeg = []
            cbw = []
            tbw = []
            dbw = []
            cg = []
            dt = []

            for v in f.vs:
                tdeg.append(deg[v.index])
                tbw.append(bw[v.index])

                if v['name'] in self.lists[groupA]:
                    cg.append(v['name'])
                    cancerg += 1
                    cdeg.append(deg[v.index])
                    cbw.append(bw[v.index])

                if v['name'] in self.lists[groupB]:
                    dt.append(v['name'])
                    drugt += 1
                    ddeg.append(deg[v.index])
                    dbw.append(bw[v.index])

            cpct = cancerg * 100 / float(len(self.lists[groupA]))
            dpct = drugt * 100 / float(len(self.lists[groupB]))
            tdgr = sum(tdeg) / float(len(tdeg))
            cdgr = sum(cdeg) / float(len(cdeg))
            ddgr = sum(ddeg) / float(len(ddeg))
            tbwn = sum(tbw) / float(len(tbw))
            cbwn = sum(cbw) / float(len(cbw))
            dbwn = sum(dbw) / float(len(dbw))
            src = []
            csrc = []
            dsrc = []

            for e in f.es:
                src.append(len(e["sources"]))

                if f.vs[e.source][groupA] or f.vs[e.target][groupA]:
                    csrc.append(len(e["sources"]))

                if f.vs[e.source][groupB] or f.vs[e.target][groupB]:
                    dsrc.append(len(e["sources"]))

            snum = sum(src) / float(len(src))
            csnum = sum(csrc) / float(len(csrc))
            dsnum = sum(dsrc) / float(len(dsrc))
            result[s] = [
                s, str(vnum), str(enum), str(cancerg), str(drugt), str(cpct),
                str(dpct), str(tdgr), str(cdgr), str(ddgr), str(tbwn),
                str(cbwn), str(dbwn), str(snum), str(csnum), str(dsnum)
            ]
        outfile = '-'.join([groupA, groupB, "stats"])

        if header is None:
            self.write_table(result, outfile, colnames=False)

        else:
            result["header"] = header
            self.write_table(result, outfile, colnames=True)

    def sources_venn_data(self, fname=None, return_data=False):
        """
        """

        result = {}
        self.update_sources()
        g = self.graph

        for i, j in itertools.product(self.sources, self.sources):

            ini = [e.index for e in g.es if i in e['sources']]
            inj = [e.index for e in g.es if j in e['sources']]

            onlyi = str(len(list(set(ini) - set(inj))))
            onlyj = str(len(list(set(inj) - set(ini))))
            inter = str(len(list(set(ini) & set(inj))))
            result[i + "-" + j] = [i, j, onlyi, onlyj, inter]

        if fname:
            self.write_table(result, fname)

        if return_data:
            return result

    def sources_hist(self):
        """
        """

        srcnum = [len(e['sources']) for e in self.graph.es]

        self.write_table({"srcnum": srcnum}, "source_num", sep=";",
                         rownames=False, colnames=False)

    def degree_dist(self, prefix, g, group=None):
        """
        """

        deg = g.vs.degree()
        self.write_table({"deg": deg}, prefix + "-whole-degdist", sep=";",
                         rownames=False, colnames=False)

        if group is not None:

            if len(set(group) - set(self.graph.vs.attributes())) > 0:
                self.ownlog.msg(2, ("Missing vertex attribute!"), 'ERROR')
                return False

            if not isinstance(group, list):
                group = [group]

            for gr in group:
                dgr = [deg[i] for i, v in enumerate(g.vs) if v[gr]]

                self.write_table({"deg": dgr}, prefix + "-" + gr + "-degdist",
                                 sep=";", rownames=False, colnames=False)

    def delete_by_source(self, source, vertexAttrsToDel=None,
                         edgeAttrsToDel=None):
        """
        """

        self.update_vertex_sources()
        g = self.graph
        verticesToDel = []

        for v in g.vs:

            if len(v['sources'] - set([source])) == 0:
                verticesToDel.append(v.index)

        g.delete_vertices(verticesToDel)
        edgesToDel = []

        for e in g.es:

            if len(e['sources'] - set([source])) == 0:
                edgesToDel.append(e.index)

            else:
                e['sources'] = set(e['sources']) - set([source])

        g.delete_edges(edgesToDel)

        if vertexAttrsToDel is not None:

            for vAttr in vertexAttrsToDel:

                if vAttr in g.vs.attributes():
                    del g.vs[vAttr]

        if edgeAttrsToDel is not None:

            for eAttr in edgeAttrsToDel:

                if eAttr in g.vs.attributes():
                    del g.vs[eAttr]

        self.update_vertex_sources()

    def reference_hist(self, filename=None):
        """
        """

        g = self.graph
        tbl = []

        for e in g.es:
            tbl.append((g.vs[e.source]['name'], g.vs[e.target]['name'],
                        str(len(e['references'])), str(len(e['sources']))))

        if filename is None:
            filename = self.outdir + "/" + self.session + "-refs-hist"

        out = ''

        for i in tbl:
            out += "\t".join(list(i)) + "\n"

        outf = codecs.open(filename, encoding='utf-8', mode='w')
        outf.write(out[:-1])
        outf.close()

    def load_resources(self, lst=None, exclude=[], cache_files={},
                       reread=False, redownload=False):
        """
        Loads multiple resources, and cleans up after.
        Looks up ID types, and loads all ID conversion
        tables from UniProt if necessary. This is much
        faster than loading the ID conversion and the
        resources one by one.
        """

        if lst is None:
            lst = omnipath

        self.load_reflists()
        huge = dict(
            (k, v) for k, v in iteritems(lst)
            if v.huge and k not in exclude and v.name not in cache_files)
        nothuge = dict(
            (k, v) for k, v in iteritems(lst)
            if (not v.huge or v.name in cache_files) and k not in exclude)

        for lst in [huge, nothuge]:

            for k, v in iteritems(lst):
                self.load_resource(
                    v,
                    clean=False,
                    cache_files=cache_files,
                    reread=reread,
                    redownload=redownload)

                # try:
                #    self.load_resource(v, clean = False,
                #        cache_files = cache_files,
                #        reread = reread,
                #        redownload = redownload)
                # except:
                #    sys.stdout.write('\t:: Could not load %s, unexpected error '\
                #        'occurred, see %s for error.\n'%(k, self.ownlog.logfile))
                #    self.ownlog.msg(1, 'Error at loading %s: \n%s\n, \t%s, %s\n' % \
                #            (k, sys.exc_info()[1],
                #            sys.exc_info()[2],
                #            sys.exc_info()[0]),
                #        'ERROR')
                #    sys.stdout.flush()

        sys.stdout.write('\n')

        self.clean_graph()
        self.update_sources()
        self.update_vertex_sources()
        self.update_pathways()
        self.update_pathway_types()

        sys.stdout.write(
            """\n > %u interactions between %u nodes\n from %u"""
            """ resources have been loaded,\n for details see the log: ./%s\n"""
            % (self.graph.ecount(), self.graph.vcount(), len(self.sources),
               self.ownlog.logfile))

    def load_mappings(self):
        """
        """

        self.mapper.load_mappings(maps=data_formats.mapList)

    def load_resource(self, settings, clean=True, cache_files={}, reread=False,
                      redownload=False):
        """
        """

        sys.stdout.write(' > ' + settings.name + '\n')
        self.read_data_file(
            settings,
            cache_files=cache_files,
            reread=reread,
            redownload=redownload)
        self.attach_network()

        if clean:
            self.clean_graph()

        self.update_sources()
        self.update_vertex_sources()

    def load_reflists(self, reflst=None):
        """
        """

        if reflst is None:
            reflst = reflists.get_reflists()

        for rl in reflst:
            self.load_reflist(rl)

    def load_reflist(self, reflist):
        """
        """

        reflist.load()
        idx = (reflist.nameType, reflist.typ, reflist.tax)
        self.reflists[idx] = reflist

    def load_negatives(self):
        """
        """

        for k, v in iteritems(negative):
            sys.stdout.write(' > ' + v.name + '\n')
            self.apply_negative(v)

    def load_tfregulons(self, levels = {'A', 'B'}, only_curated = False):
        """
        Adds TF-target interactions from TF regulons to the network.

        :param set levels:
            Confidence levels to be used.
        :param bool only_curated:
            Retrieve only literature curated interactions.

        Details
        -------
        TF regulons is a comprehensive resource of TF-target interactions
        combining multiple lines of evidences: literature curated databases,
        ChIP-Seq data, PWM based prediction using HOCOMOCO and JASPAR matrices
        and prediction from GTEx expression data by ARACNe.

        For details see https://github.com/saezlab/DoRothEA.

        Example
        -------
        >>> import pypath
        >>> pa = pypath.PyPath()
        >>> pa.load_tfregulons(levels = {'A'})
        """

        settings = copy.deepcopy(data_formats.transcription['tfregulons'])
        settings.inputArgs = {
            'levels': levels,
            'only_curated': only_curated
        }

        self.load_resources({'tfregulons': settings})

    def list_resources(self):
        """
        """

        sys.stdout.write(' > omnipath\n')

        for k, v in iteritems(omnipath):
            sys.stdout.write('\t:: %s (%s)\n' % (v.name, k))

        sys.stdout.write(' > good\n')

        for k, v in iteritems(good):
            sys.stdout.write('\t:: %s (%s)\n' % (v.name, k))

    def info(self, name):
        """
        """

        d = descriptions.descriptions

        if name not in d:
            sys.stdout.write(' :: Sorry, no description available about %s\n' %
                             name)
            return None

        dd = d[name]
        out = '\n\t::: %s :::\n\n' % name

        if 'urls' in dd:

            if 'webpages' in dd['urls']:
                out += '\t:: Webpages:\n'

                for w in dd['urls']['webpages']:
                    out += '\t    > %s\n' % w

            if 'articles' in dd['urls']:
                out += '\t:: Articles:\n'

                for w in dd['urls']['articles']:
                    out += '\t    > %s\n' % w

        if 'taxons' in dd:
            out += '\t:: Taxons: %s\n' % ', '.join(dd['taxons'])

        if 'descriptions' in dd and len(dd['descriptions']) > 0:
            out += '\n\t:: From the authors:\n'
            txt = dd['descriptions'][0].split('\n')
            txt = '\n\t'.join(['\n\t'.join(textwrap.wrap(t, 50)) for t in txt])
            out += '\t\t %s\n' % txt.replace('\t\t', '\t    ')

        if 'notes' in dd and len(dd['notes']) > 0:
            out += '\n\t:: Notes:\n'
            txt = dd['notes'][0].split('\n')
            txt = '\n\t'.join(['\n\t'.join(textwrap.wrap(t, 50)) for t in txt])
            out += '\t\t %s\n\n' % txt.replace('\t\t', '\t    ')

        sys.stdout.write(out)

    #
    # functions to make life easier
    #

    def having_attr(self, attr, graph=None, index=True, edges=True):
        """
        """

        graph = graph or self.graph
        es_or_vs = getattr(graph, 'es' if edges else 'vs')

        if attr in es_or_vs.attributes():

            for i in es_or_vs:

                if something(i[attr]):
                    yield i.index if index else i

    def having_eattr(self, attr, graph=None, index=True):
        """
        """

        return self.having_attr(attr, graph, index)

    def having_vattr(self, attr, graph=None, index=True):
        """
        """

        return self.having_attr(attr, graph, index, False)

    def having_ptm(self, index=True, graph=None):
        """
        """

        return self.having_eattr('ptm', graph, index)

    def loop_edges(self, index=True, graph=None):
        """
        """

        graph = graph or self.graph

        for e in graph.es:

            if e.source == e.target:
                yield e.index if index else e

    #
    # functions to make topological analysis on the graph
    #

    def first_neighbours(self, node, indices=False):
        """
        """

        g = self.graph
        lst = []

        if not isinstance(node, int):
            node = g.vs.select(name=node)

            if len(node) > 0:
                node = node[0].index

            else:
                return lst

        lst = list(set(g.neighborhood(node)) - set([node]))

        if indices:
            return lst

        else:
            nlst = []

            for v in lst:
                nlst.append(g.vs[v]['name'])

            return nlst

    def second_neighbours(self, node, indices=False, with_first=False):
        """
        """

        g = self.graph
        lst = []

        if isinstance(node, int):
            node_i = node
            node_n = g.vs[node_i]['name']

        else:
            node_i = g.vs.select(name=node)

            if len(node_i) > 0:
                node_i = node_i[0].index
                node_n = node

            else:
                return lst

        first = self.first_neighbours(node_i, indices=indices)

        for n in first:
            lst += self.first_neighbours(n, indices=indices)

        if with_first:
            lst += first

        else:
            lst = list(set(lst) - set(first))

        if indices:
            return list(set(lst) - set([node_i]))

        else:
            return list(set(lst) - set([node_n]))

    def all_neighbours(self, indices=False):
        """
        """

        g = self.graph
        g.vs['neighbours'] = [[] for _ in xrange(g.vcount())]
        prg = Progress(
            total=g.vcount(), name="Searching neighbours", interval=30)

        for v in g.vs:
            v['neighbours'] = self.first_neighbours(v.index, indices=indices)
            prg.step()

        prg.terminate()

    def jaccard_edges(self):
        """
        """

        g = self.graph
        self.all_neighbours(indices=True)
        metaEdges = []
        prg = Progress(
            total=g.vcount(), name="Calculating Jaccard-indices", interval=11)

        for v in xrange(0, g.vcount() - 1):

            for w in xrange(v + 1, g.vcount()):
                vv = g.vs[v]
                vw = g.vs[w]
                ja = (len(set(vv['neighbours']) & set(vw['neighbours'])) /
                      float(len(vv['neighbours']) + len(vw['neighbours'])))
                metaEdges.append((vv['name'], vw['name'], ja))

            prg.step()

        prg.terminate()

        return metaEdges

    def jaccard_meta(self, jedges, critical):
        """
        """

        edges = []

        for e in jedges:

            if e[2] > critical:
                edges.append((e[0], e[1]))

        return igraph.Graph.TupleList(edges)

    def apply_negative(self, settings):
        """
        """

        g = self.graph

        if settings.name not in self.negatives:
            self.raw_data = None
            self.read_data_file(settings)
            self.negatives[settings.name] = self.raw_data

        neg = self.negatives[settings.name]
        prg = Progress(
            total=len(neg), name="Matching interactions", interval=11)
        matches = 0

        g.es['negative'] = [
            set([]) if e['negative'] is None else e['negative'] for e in g.es
        ]
        g.es['negative_refs'] = [
            set([]) if e['negative_refs'] is None else e['negative_refs']
            for e in g.es
        ]

        for n in neg:
            aexists = n["defaultNameA"] in g.vs['name']
            bexists = n["defaultNameB"] in g.vs['name']

            if aexists and bexists:
                edge = self.edge_exists(n["defaultNameA"], n["defaultNameB"])

                if isinstance(edge, int):
                    g.es[edge]['negative'].add(settings.name)
                    refs = set(
                        list(
                            map(lambda r: _refs.Reference(int(r)), n[
                                'attrsEdge']['references'])))
                    g.es[edge]['negative_refs'].add(refs)
                    matches += 1

            prg.step()

        prg.terminate()
        sys.stdout.write('\t%u matches found with negative set\n' % matches)

    def negative_report(self, lst=True, outFile=None):
        """
        """

        if outFile is None:
            outFile = self.outdir + self.session + '-negatives'

        self.genesymbol_labels()
        out = ''
        neg = []
        g = self.graph

        for e in g.es:

            if len(e['negative']) > 0:

                if outFile:
                    out += '\t'.join([
                        g.vs[e.source]['name'], g.vs[e.target]['name'],
                        g.vs[e.source]['label'], g.vs[e.target]['label'],
                        ';'.join(list(e['sources'])),
                        ';'.join(map(lambda r: r.pmid, e['references'])),
                        ';'.join(e['negative']), ';'.join(e['negative_refs'])
                    ]) + '\n'

                if lst:
                    neg.append(e)

        if outFile:
            outf = codecs.open(outFile, encoding='utf-8', mode='w')
            outf.write(out)
            outf.close()

        if lst:
            return neg

    def export_ptms_tab(self, outfile=None):
        """
        """

        if outfile is None:
            outfile = os.path.join(self.outdir,
                                   'network-' + self.session + '.tab')

        self.genesymbol_labels()
        self.update_vname()
        g = self.graph

        if 'ddi' not in g.es.attributes():
            g.es['ddi'] = [[] for _ in g.es]

        if 'ptm' not in g.es.attributes():
            g.es['ptm'] = [[] for _ in g.es]

        header = [
            'UniProt_A', 'UniProt_B', 'GeneSymbol_A', 'GeneSymbol_B',
            'Databases', 'PubMed_IDs', 'Stimulation', 'Inhibition',
            'Substrate-isoform', 'Residue_number', 'Residue_letter', 'PTM_type'
        ]
        stripJson = re.compile(r'[\[\]{}\"]')
        # first row is header
        outl = [header]

        with codecs.open(outfile, encoding='utf-8', mode='w') as f:
            f.write('\t'.join(header) + '\n')
            prg = Progress(
                total=self.graph.ecount(), name='Writing table', interval=31)
            uniqedges = []

            for e in g.es:
                prg.step()
                # only directed

                for di in e['dirs'].which_dirs():
                    src = self.nodDct[di[0]]
                    tgt = self.nodDct[di[1]]
                    # uniprot names
                    row = list(di)
                    # genesymbols
                    row += [g.vs[src]['label'], g.vs[tgt]['label']]
                    # sources
                    dbs = e['dirs'].get_dir(di, sources=True)
                    row.append(';'.join(dbs))
                    # references
                    row.append(';'.join([
                        r
                        for rs in [
                            refs for db, refs in iteritems(e['refs_by_source'])
                            if db in dbs
                        ] for r in rs
                    ]))
                    # signs
                    row += [str(int(x)) for x in e['dirs'].get_sign(di)]
                    # domain-motif
                    # row.append('#'.join([x.print_residues() for x in e['ptm'] \
                    #    if x.__class__.__name__ == 'DomainMotif']))

                    for dmi in e['ptm']:

                        if dmi.__class__.__name__ == 'DomainMotif':

                            if dmi.ptm.residue is not None:

                                if dmi.ptm.residue.protein == di[1]:
                                    uniqedges.append(e.index)
                                    r = row + [
                                        '%s-%u' %
                                        (dmi.ptm.protein, dmi.ptm.isoform
                                         ), str(dmi.ptm.residue.number),
                                        dmi.ptm.residue.name, dmi.ptm.typ
                                    ]
                                    # here each ptm in separate row:
                                    outl.append(r)
                    # row complete, appending to main list
                    # outl.append(row)
                    # tabular text from list of lists; writing to file

            out = '\n'.join(['\t'.join(row) for row in outl])

            with codecs.open(outfile, encoding='utf-8', mode='w') as f:
                f.write(out)

            prg.terminate()
            console(':: Data has been written to %s' % outfile)

            return list(set(uniqedges))

    def export_struct_tab(self, outfile=None):
        """
        """

        if outfile is None:
            outfile = os.path.join(self.outdir,
                                   'network-' + self.session + '.tab')

        self.genesymbol_labels()
        self.update_vname()
        g = self.graph

        if 'ddi' not in g.es.attributes():
            g.es['ddi'] = [[] for _ in g.es]

        if 'ptm' not in g.es.attributes():
            g.es['ptm'] = [[] for _ in g.es]

        header = [
            'UniProt_A', 'UniProt_B', 'GeneSymbol_B', 'GeneSymbol_A',
            'Databases', 'PubMed_IDs', 'Stimulation', 'Inhibition',
            'Domain-domain', 'Domain-motif-PTM', 'PTM-regulation'
        ]
        stripJson = re.compile(r'[\[\]{}\"]')
        # first row is header
        outl = [header]

        with codecs.open(outfile, encoding='utf-8', mode='w') as f:
            f.write('\t'.join(header) + '\n')
            prg = Progress(
                total=self.graph.ecount(), name='Writing table', interval=31)

            for e in g.es:
                prg.step()
                # only directed

                for di in e['dirs'].which_dirs():
                    src = self.nodDct[di[0]]
                    tgt = self.nodDct[di[1]]
                    # uniprot names
                    row = list(di)
                    # genesymbols
                    row += [g.vs[src]['label'], g.vs[tgt]['label']]
                    # sources
                    dbs = e['dirs'].get_dir(di, sources=True)
                    row.append(';'.join(dbs))
                    # references
                    row.append(';'.join([
                        r
                        for rs in [
                            refs for db, refs in iteritems(e['refs_by_source'])
                            if db in dbs
                        ] for r in rs
                    ]))
                    # signs
                    row += [str(int(x)) for x in e['dirs'].get_sign(di)]
                    # domain-domain
                    row.append('#'.join([x.serialize() for x in e['ddi']]))
                    # domain-motif
                    row.append('#'.join([
                        x.serialize() for x in e['ptm']
                        if x.__class__.__name__ == 'Ptm'
                    ]))
                    # domain-motif
                    row.append('#'.join([
                        x.serialize() for x in e['ptm']
                        if x.__class__.__name__ == 'Regulation'
                    ]))
                    # row complete, appending to main list
                    outl.append(row)

            # tabular text from list of lists; writing to file
            out = '\n'.join(['\t'.join(row) for row in outl])

            with codecs.open(outfile, encoding='utf-8', mode='w') as f:
                f.write(out)

            prg.terminate()
            console(':: Data has been written to %s' % outfile)

    def export_tab(self, outfile=None, extra_node_attrs={},
                   extra_edge_attrs={}, unique_pairs=True, **kwargs):
        """
        Exports the network in a tabular format.

        By default UniProt IDs, Gene Symbols, source databases, literature
        references, directionality and sign information and interaction type
        are included.

        Args:
        -----
        :param str outfile:
            Name of the output file. If `None` a file name
            "netrowk-<session id>.tab" is used.
        :param dict extra_node_attrs:
            Additional node attributes to be included in the exported table.
            Keys are column ames used in the header while values are names
            of vertex attributes. In the header `_A` and `_B` suffixes will
            be appended to the column names so the values can be assigned to
            A and B side interaction partners.
        :param dict extra_edge_attrs:
            Additional edge attributes to be included in the exported table.
            Keys are column ames used in the header while values are names
            of edge attributes.
        """

        e = export.Export(
            pa = self,
            extra_node_attrs = extra_node_attrs,
            extra_edge_attrs = extra_edge_attrs,
            **kwargs
        )
        e.write_tab(unique_pairs = unique_pairs, outfile = outfile)

    def export_sif(self, outfile=None):
        """
        """

        outfile = outfile if outfile is not None \
            else 'network-%s.sif' % self.session

        with open(outfile, 'w') as f:

            for e in self.graph.es:

                for d in [d for d, b in iteritems(e['dirs'].dirs) if b]:

                    if e['dirs'].is_directed() and d == 'undirected':
                        continue

                    sign = '' if d == 'undirected' \
                        else ''.join([['+', '-'][i]
                                      for i, v in enumerate(e['dirs'].get_sign(d)) if v])
                    dirn = '=' if d == 'undirected' else '>'
                    source = self.graph.vs[e.source]['name'] \
                        if d == 'undirected' else d[0]
                    target = self.graph.vs[e.target]['name'] \
                        if d == 'undirected' else d[1]
                    f.write('\t'.join([source, sign + dirn, target]) + '\n')

    def export_graphml(self, outfile=None, graph=None, name='main'):
        """
        """

        self.genesymbol_labels()
        g = self.graph if graph is None else graph

        if outfile is None:
            outfile = os.path.join(self.outdir,
                                   'network-' + self.session + '.graphml')

        isDir = 'directed' if g.is_directed() else 'undirected'
        isDirB = 'true' if g.is_directed() else 'false'
        nodeAttrs = [('UniProt', 'string'), ('GeneSymbol', 'string'),
                     ('Type', 'string')]
        edgeAttrs = [('Databases', 'string'), ('PubMedIDs', 'string'),
                     ('Undirected', 'string'), ('DirectionAB', 'string'),
                     ('DirectionBA', 'string'), ('StimulatoryAB', 'string'),
                     ('InhibitoryAB', 'string'), ('StimulatoryBA', 'string'),
                     ('InhibitoryBA', 'string'), ('Category', 'string')]
        header = """<?xml version="1.0" encoding="UTF-8"?>
                    <graphml xmlns="http://graphml.graphdrawing.org/xmlns"
                    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
                    http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">\n\n
                """

        with codecs.open(outfile, encoding='utf-8', mode='w') as f:
            f.write(header)

            for attr in nodeAttrs:
                f.write(
                    '\t<key id="%s" for="node" attr.name="%s" attr.type="%s" />\n'
                    % (attr[0], attr[0], attr[1]))

            for attr in edgeAttrs:
                f.write(
                    '\t<key id="%s" for="edge" attr.name="%s" attr.type="%s" />\n'
                    % (attr[0], attr[0], attr[1]))

            f.write("""\n<graph id="%s" edgedefault="%s"
                        parse.nodeids="free" parse.edgeids="canonical"
                        parse.order="nodesfirst">\n\n""" % (name, isDir))
            prg = Progress(total=g.vcount(), name='Writing nodes', interval=31)
            f.write('\n\t<!-- graph properties -->\n\n')
            f.write('\n\t<!-- vertices -->\n\n')

            for v in g.vs:
                f.write('<node id="%s">\n' % (v['name']))
                f.write('\t<data key="UniProt">%s</data>\n' % (v['name']))
                f.write('\t<data key="GeneSymbol">%s</data>\n' %
                        (v['label'].replace(' ', '')))
                f.write('\t<data key="Type">%s</data>\n' % (v['type']))
                f.write('</node>\n')

            prg.terminate()
            prg = Progress(total=g.ecount(), name='Writing edges', interval=31)
            f.write('\n\t<!-- edges -->\n\n')

            for e in g.es:
                f.write(
                    '<edge id="%s_%s" source="%s" target="%s" directed="%s">\n'
                    % (g.vs[e.source]['name'], g.vs[e.target]['name'],
                       g.vs[e.source]['name'], g.vs[e.target]['name'], isDirB))
                f.write('\t<data key="Databases">%s</data>\n' %
                        (';'.join(list(e['sources']))))
                f.write(
                    '\t<data key="PubMedIDs">%s</data>\n' %
                    (';'.join(list(map(lambda r: r.pmid, e['references'])))))
                f.write('\t<data key="Undirected">%s</data>\n' %
                        (';'.join(common.uniqList(e['dirs_by_source'][0]))))
                f.write('\t<data key="DirectionAB">%s</data>\n' %
                        (';'.join(common.uniqList(e['dirs_by_source'][1]))))
                f.write('\t<data key="DirectionBA">%s</data>\n' %
                        (';'.join(common.uniqList(e['dirs_by_source'][2]))))
                f.write('\t<data key="StimulatoryAB">%s</data>\n' %
                        (';'.join(common.uniqList(e['signs'][0][0]))))
                f.write('\t<data key="InhibitoryAB">%s</data>\n' %
                        (';'.join(common.uniqList(e['signs'][0][1]))))
                f.write('\t<data key="StimulatoryBA">%s</data>\n' %
                        (';'.join(common.uniqList(e['signs'][1][0]))))
                f.write('\t<data key="InhibitoryBA">%s</data>\n' %
                        (';'.join(common.uniqList(e['signs'][1][1]))))
                f.write('\t<data key="InhibitoryBA">%s</data>\n' % (e['type']))
                f.write('</edge>\n')
                prg.step()

            f.write('\n\t</graph>\n</graphml>')

        prg.terminate()

    def compounds_from_chembl(self, chembl_mysql=None, nodes=None, crit=None,
                              andor="or", assay_types=['B', 'F'],
                              relationship_types=['D', 'H'], multi_query=False,
                              **kwargs):
        """
        """

        chembl_mysql = chembl_mysql or self.chembl_mysql
        self.chembl = chembl.Chembl(
            chembl_mysql, self.ncbi_tax_id, mapper=self.mapper)

        if nodes is None:

            if crit is None:
                nodes = xrange(0, self.graph.vcount())

            else:
                sub = self.get_sub(crit=crit, andor=andor)
                nodes = sub['nodes']

        uniprots = []

        for v in self.graph.vs:

            if v.index in nodes and v['nameType'] == 'uniprot':
                uniprots.append(v['name'])

        self.chembl.compounds_targets(
            uniprots,
            assay_types=assay_types,
            relationship_types=relationship_types,
            **kwargs)
        self.chembl.compounds_by_target()
        self.update_vname()
        self.graph.vs['compounds_chembl'] = [
            [] for _ in xrange(self.graph.vcount())
        ]
        self.graph.vs['compounds_names'] = [
            [] for _ in xrange(self.graph.vcount())
        ]
        self.graph.vs['compounds_data'] = [[]
                                           for _ in xrange(self.graph.vcount())
                                           ]
        num = len(self.chembl.compounds)
        prg = Progress(total=num, name='Adding compounds', interval=11)
        hascomp = 0

        for target, data in iteritems(self.chembl.compounds):
            prg.step()

            if self.node_exists(target):
                hascomp += 1
                node = self.graph.vs[self.graph.vs['name'].index(target)]
                node['compounds_data'] = data

                for comp in data:
                    node['compounds_chembl'].append(comp['chembl'])
                    node['compounds_names'] += comp['compound_names']

                node['compounds_chembl'] = common.uniqList(node[
                    'compounds_chembl'])
                node['compounds_names'] = common.uniqList(node[
                    'compounds_names'])

        prg.terminate()
        percent = hascomp / float(self.graph.vcount())
        sys.stdout.write(
            '\n\tCompounds found for %u targets, (%.2f%% of all proteins).\n\n'
            % (hascomp, percent * 100.0))

    def network_filter(self, p=2.0):
        """
        This function aims to cut the number of edges in the network,
        without loosing nodes, to make the network less connected,
        less hairball-like, more usable for analysis.
        """

        ref_freq = {}

        for s in self.sources:
            ref_freq[s] = {}

            for e in self.graph.es:

                if s in e['sources']:

                    for r in e['refs_by_source'][s]:

                        if r not in ref_freq[s]:
                            ref_freq[s][r] = 1

                        else:
                            ref_freq[s][r] += 1

        self.graph.es['score'] = [0.0]
        deg = self.graph.vs.degree()
        avdeg = sum(deg) / len(deg)
        prg = Progress(self.graph.ecount(), 'Calculating scores', 11)

        for e in self.graph.es:
            score = 0.0

            for s, rs in iteritems(e['refs_by_source']):

                for r in rs:
                    score += 1.0 / ref_freq[s][r]

            mindeg = min(self.graph.vs[e.source].degree(),
                         self.graph.vs[e.target].degree())

            if mindeg < avdeg:
                score *= pow((mindeg - avdeg), p)

            e['score'] = score
            prg.step()

        prg.terminate()

    def shortest_path_dist(self, graph=None, subset=None, outfile=None,
                           **kwargs):
        """
        subset is a tuple of two lists if you wish to look for
        paths between elements of two groups, or a list if you
        wish to look for shortest paths within this group
        """

        graph = graph if graph is not None else self.graph
        shortest_paths = []
        subset = subset if isinstance(subset, tuple) or subset is None else (
            subset, subset)
        prg = Progress(graph.vcount(), 'Calculating paths', 1)

        for i in xrange(0, graph.vcount() - 1):

            if subset is None or i in subset[0] or i in subset[1]:
                paths = graph.get_shortest_paths(i,
                                                 xrange(i + 1, graph.vcount()),
                                                 **kwargs)

                for j in xrange(0, len(paths)):

                    if subset is None or (
                            i in subset[0] and i + j + 1 in subset[1]) or (
                                i in subset[1] and i + j + 1 in subset[0]):
                        shortest_paths.append(len(paths[j]))

            prg.step()

        prg.terminate()

        if outfile is not None:
            out = '\n'.join([str(i) for i in shortest_paths])

            with codecs.open(outfile, encoding='utf-8', mode='w') as f:
                f.write(out)

        return shortest_paths

    def load_pdb(self, graph=None):
        """
        """

        graph = graph if graph is not None else self.graph
        u_pdb, pdb_u = dataio.get_pdb()

        if u_pdb is None:
            self.ownlog.msg(2, 'Failed to download UniProt-PDB dictionary',
                            'ERROR')

        else:
            graph.vs['pdb'] = [None]

            for v in graph.vs:
                v['pdb'] = {}

                if v['name'] in u_pdb:

                    for pdb in u_pdb[v['name']]:
                        v['pdb'][pdb[0]] = (pdb[1], pdb[2])

            self.ownlog.msg(2, 'PDB IDs for proteins has been retrieved.',
                            'INFO')

    def load_pfam(self, graph=None):
        """
        """

        graph = graph if graph is not None else self.graph
        u_pfam, pfam_u = dataio.get_pfam(graph.vs['name'])

        if u_pfam is None:
            self.ownlog.msg(2, 'Failed to download Pfam data from UniProt',
                            'ERROR')

        else:
            graph.vs['pfam'] = [None]

            for v in graph.vs:
                v['pfam'] = []

                if v['name'] in u_pfam:
                    v['pfam'] += u_pfam[v['name']]

            self.ownlog.msg(2, 'Pfam domains has been retrieved.', 'INFO')

    def load_pfam2(self):
        """
        """

        self.pfam_regions()

        if self.u_pfam is None:
            self.ownlog.msg(2, 'Failed to download data from Pfam', 'ERROR')

        else:
            self.graph.vs['pfam'] = [{} for _ in self.graph.vs]

            for v in self.graph.vs:

                if v['name'] in self.u_pfam:
                    v['pfam'] = self.u_pfam[v['name']]

            self.ownlog.msg(2, 'Pfam domains has been retrieved.', 'INFO')

    def load_pfam3(self):
        """
        """

        self.pfam_regions()

        if self.u_pfam is None:
            self.ownlog.msg(2, 'Failed to download data from Pfam', 'ERROR')

        else:
            self.graph.vs['doms'] = [[] for _ in self.graph.vs]

            for v in self.graph.vs:

                if v['name'] in self.u_pfam:

                    for pfam, regions in iteritems(self.u_pfam[v['name']]):

                        for region in regions:
                            v['doms'].append(
                                intera.Domain(
                                    protein=v['name'],
                                    domain=pfam,
                                    start=region['start'],
                                    end=region['end'],
                                    isoform=region['isoform']))

            self.ownlog.msg(2, 'Pfam domains has been retrieved.', 'INFO')

    def load_corum(self, graph=None):
        """
        Loads complexes from CORUM database. Loads data into vertex attribute
        `graph.vs['complexes']['corum']`.
        This resource is human only.
        """

        graph = graph if graph is not None else self.graph
        complexes, members = dataio.get_corum()

        if complexes is None:
            self.ownlog.msg(2, 'Failed to download data from CORUM', 'ERROR')

        else:
            self.init_complex_attr(graph, 'corum')

            for u, cs in iteritems(members):
                sw = self.mapper.map_name(u, 'uniprot', 'uniprot', 9606)

                for s in sw:

                    if s in graph.vs['name']:

                        for c in cs:
                            others = []

                            for memb in complexes[c[0]][0]:
                                others += self.mapper.map_name(memb,
                                                               'uniprot',
                                                               'uniprot',
                                                               9606)

                            graph.vs.select(
                                name=s)[0]['complexes']['corum'][c[1]] = {
                                    'full_name': c[0],
                                    'all_members': others,
                                    'all_members_original': complexes[c[0]][0],
                                    'references': c[2],
                                    'functions': c[4],
                                    'diseases': c[5],
                                }

            self.ownlog.msg(2, 'Complexes from CORUM have been retrieved.',
                            'INFO')

    def init_complex_attr(self, graph, name):
        """
        """

        if 'complexes' not in graph.vs.attributes():
            graph.vs['complexes'] = [None]

        for v in graph.vs:

            if v['complexes'] is None:
                v['complexes'] = {}

            if name not in v['complexes']:
                v['complexes'][name] = {}

    def load_havugimana(self, graph=None):
        """
        Loads complexes from Havugimana 2012. Loads data into vertex attribute
        `graph.vs['complexes']['havugimana']`.
        This resource is human only.
        """

        graph = graph if graph is not None else self.graph
        complexes = dataio.read_complexes_havugimana()

        if complexes is None:
            self.ownlog.msg(2, 'Failed to read data from Havugimana', 'ERROR')

        else:
            self.init_complex_attr(graph, 'havugimana')

            for c in complexes:
                membs = []
                names = []

                for memb in c:
                    membs += self.mapper.map_name(memb,
                                                  'uniprot',
                                                  'uniprot',
                                                  9606)

                for u in membs:
                    names += self.mapper.map_name(u,
                                                  'uniprot',
                                                  'genesymbol',
                                                  9606)

                names = sorted(set(names))
                name = ':'.join(names)

                for u in membs:

                    if u in graph.vs['name']:
                        graph.vs.select(
                            name=u)[0]['complexes']['havugimana'][name] = {
                                'full_name': '-'.join(names) + ' complex',
                                'all_members': membs,
                                'all_members_original': c
                            }

                self.ownlog.msg(
                    2, 'Complexes from Havugimana have been retrieved.',
                    'INFO')

    def load_compleat(self, graph=None):
        """
        Loads complexes from Compleat. Loads data into vertex attribute
        `graph.vs['complexes']['compleat']`.
        This resource is human only.
        """

        graph = graph if graph is not None else self.graph
        complexes = dataio.get_compleat()

        if complexes is None:
            self.ownlog.msg(2, 'Failed to load data from COMPLEAT', 'ERROR')

        else:
            self.init_complex_attr(graph, 'compleat')

            for c in complexes:
                c['uniprots'] = []
                c['gsymbols'] = []

                for e in c['entrez']:
                    c['uniprots'] += self.mapper.map_name(e,
                                                          'entrez',
                                                          'uniprot',
                                                          9606)

                for u in c['uniprots']:
                    c['gsymbols'] += self.mapper.map_name(u,
                                                          'uniprot',
                                                          'genesymbol',
                                                          9606)

                c['gsymbols'] = list(
                    set([gs.replace('; ', '') for gs in c['gsymbols']]))

                if len(c['uniprots']) > 0:

                    for u in c['uniprots']:

                        if u in graph.vs['name']:
                            name = '-'.join(c['gsymbols'])
                            graph.vs.select(
                                name=u)[0]['complexes']['compleat'][name] = {
                                    'full_name': name + ' complex',
                                    'all_members': c['uniprots'],
                                    'all_members_original': c['entrez'],
                                    'source': c['source'],
                                    'references': c['pubmeds']
                                }

            self.ownlog.msg(2, 'Complexes from COMPLEAT have been retrieved.',
                            'INFO')

    def load_complexportal(self, graph=None):
        """
        Loads complexes from ComplexPortal. Loads data into vertex attribute
        `graph.vs['complexes']['complexportal']`.
        This resource is human only.
        """

        graph = graph if graph is not None else self.graph
        # TODO: handling species
        complexes = dataio.get_complexportal()

        if complexes is None:
            self.ownlog.msg(2, 'Failed to read data from Havugimana', 'ERROR')

        else:

            if 'complexes' not in graph.vs.attributes():
                graph.vs['complexes'] = [None]

            for v in graph.vs:

                if v['complexes'] is None:
                    v['complexes'] = {}

                if 'complexportal' not in v['complexes']:
                    v['complexes']['complexportal'] = {}

            for c in complexes:
                swprots = []

                if 'complex recommended name' in c['names']:
                    name = c['names']['complex recommended name']

                else:
                    name = c['fullname']

                for u in c['uniprots']:
                    swprots += self.mapper.map_name(u,
                                                    'uniprot',
                                                    'uniprot',
                                                    9606)

                swprots = list(set(swprots))

                for sw in swprots:

                    if sw in graph.vs['name']:
                        v = graph.vs.select(name=sw)[0]
                        v['complexes']['complexportal'][name] = {
                            'full_name': c['fullname'],
                            'all_members': swprots,
                            'all_members_original': c['uniprots'],
                            'pdbs': c['pdbs'],
                            'references': c['pubmeds'],
                            'synonyms': c['names'],
                            'description': c['description']
                        }

            self.ownlog.msg(
                2, 'Complexes from Complex Portal have been retrieved.',
                'INFO')

    def load_3dcomplexes(self, graph=None):
        """
        """

        graph = graph if graph is not None else self.graph
        c3d = dataio.get_3dcomplexes()

        if c3d is None:
            self.ownlog.msg(
                2, 'Failed to download data from 3DComplexes and PDB', 'ERROR')

        else:

            if 'complexes' not in graph.vs.attributes():
                graph.vs['complexes'] = [None]

            for v in graph.vs:

                if v['complexes'] is None:
                    v['complexes'] = {}

                if '3dcomplexes' not in v['complexes']:
                    v['complexes']['3dcomplexes'] = {}

            if '3dstructure' not in graph.es.attributes():
                graph.es['3dstructure'] = [None]

            for e in graph.es:

                if e['3dstructure'] is None:
                    e['3dstructure'] = {}

                if '3dcomplexes' not in e['3dstructure']:
                    e['3dstructure']['3dcomplexes'] = []

            compl_names = {}
            compl_membs = {}
            compl_links = []
            prg = Progress(len(c3d), 'Processing complexes', 7)

            for cname, v in iteritems(c3d):
                swprots = []
                uprots = []

                for ups, contact in iteritems(v):
                    uprots += list(ups)
                    up1s = self.mapper.map_name(ups[0], 'uniprot', 'uniprot')

                    for up1 in up1s:
                        inRefLists = False

                        for tax, lst in iteritems(self.reflists):

                            if up1 in lst.lst:
                                up2s = self.mapper.map_name(ups[1], 'uniprot',
                                                            'uniprot')

                                for up2 in up2s:

                                    for tax, lst in iteritems(self.reflists):

                                        if up2 in lst.lst:
                                            inRefLists = True
                                            thisPair = sorted([up1, up2])
                                            swprots += thisPair
                                            compl_links.append(
                                                (thisPair, contact, cname))
                                            break

                swprots = list(set(swprots))
                uprots = list(set(uprots))

                if len(swprots) > 0:
                    name = []

                    for sp in swprots:
                        name += self.mapper.map_name(sp,
                                                     'uniprot',
                                                     'genesymbol',
                                                     9606)

                    compl_names[cname] = '-'.join(name) + ' complex'
                    compl_membs[cname] = (swprots, uprots)

                prg.step()

            prg.terminate()
            prg = Progress(len(compl_membs), 'Processing nodes', 3)

            for cname, uniprots in iteritems(compl_membs):

                for sp in uniprots[0]:

                    if sp in graph.vs['name']:
                        graph.vs.select(
                            name=sp)[0]['complexes']['3dcomplexes'][cname] = {
                                'full_name': compl_names[cname],
                                'all_members': uniprots[0],
                                'all_members_original': uniprots[1]
                            }

                prg.step()

            prg.terminate()
            prg = Progress(len(compl_links), 'Processing edges', 3)

            for links in compl_links:
                one = links[0][0]
                two = links[0][1]

                if one in graph.vs['name'] and two in graph.vs['name']:
                    nodeOne = graph.vs.select(name=one)[0].index
                    nodeTwo = graph.vs.select(name=two)[0].index

                    try:
                        edge = graph.es.select(_between=((nodeOne, ),
                                                         (nodeTwo, )))[0]
                        edge['3dstructure']['3dcomplexes'].append({
                            'pdb': links[2].split('_')[0],
                            'numof_residues': links[1]
                        })

                    except:
                        pass

                prg.step()
            prg.terminate()

    def load_pisa(self, graph=None):
        """
        """

        graph = graph if graph is not None else self.graph

        if 'complexes' not in graph.vs.attributes() \
                or '3dcomplexes' not in graph.vs[0]['complexes']:
            self.load_3dcomplexes(graph=graph)

        if 'complexportal' not in graph.vs[0]['complexes']:
            self.load_complexportal(graph=graph)

        pdblist = []

        for v in graph.vs:

            for n, d in iteritems(v['complexes']['3dcomplexes']):
                pdblist.append(n.split('_')[0])

            for n, d in iteritems(v['complexes']['complexportal']):
                pdblist += d['pdbs']

        pdblist = list(set(pdblist))
        pisa, unmapped = dataio.get_pisa(pdblist)

        if 'interfaces' not in graph.es.attributes():
            graph.es['interfaces'] = [None]

        for e in graph.es:

            if e['interfaces'] is None:
                e['interfaces'] = {}

            if 'pisa' not in e['interfaces']:
                e['interfaces']['pisa'] = {}

        for pdb, iflist in iteritems(pisa):

            for uniprots, intf in iteritems(iflist):

                if uniprots[0] in graph.vs['name'] and uniprots[1] in graph.vs[
                        'name']:
                    e = self.edge_exists(uniprots[0], uniprots[1])

                    if isinstance(e, int):

                        if pdb not in graph.es[e]['interfaces']['pisa']:
                            graph.es[e]['interfaces']['pisa'][pdb] = []

                        graph.es[e]['interfaces']['pisa'][pdb].append(intf)

        return unmapped

    def find_complex(self, search):
        """
        Finds complexes by their non standard names.
        E.g. to find DNA polymerases you can use the search
        term `DNA pol` which will be tested against complex names
        in CORUM.
        """

        result = []
        comp, memb = dataio.get_corum()

        for cname, cdata in comp.items():

            if search.lower() in cname.lower():
                result.append((cname, cdata[0]))

        return result

    def genesymbol(self, genesymbol):
        """
        Returns ``igraph.Vertex()`` object if the GeneSymbol
        can be found in the default undirected network,
        otherwise ``None``.

        @genesymbol : str
            GeneSymbol.
        """

        graph = self._get_undirected()
        return graph.vs[self.labDct[genesymbol]] \
            if genesymbol in self.labDct else None

    gs = genesymbol

    def dgenesymbol(self, genesymbol):
        """
        Returns ``igraph.Vertex()`` object if the GeneSymbol
        can be found in the default directed network,
        otherwise ``None``.

        @genesymbol : str
            GeneSymbol.
        """

        dgraph = self._get_directed()
        return dgraph.vs[self.dlabDct[genesymbol]] \
            if genesymbol in self.dlabDct else None

    dgs = dgenesymbol

    def genesymbols(self, genesymbols):
        """
        """

        return filter(lambda v: v is not None,
                      map(self.genesymbol, genesymbols))

    gss = genesymbols

    def dgenesymbols(self, genesymbols):
        """
        """

        return filter(lambda v: v is not None,
                      map(self.dgenesymbol, genesymbols))

    dgss = dgenesymbols

    def uniprot(self, uniprot):
        """
        Returns ``igraph.Vertex()`` object if the UniProt
        can be found in the default undirected network,
        otherwise ``None``.

        @uniprot : str
            UniProt ID.
        """

        graph = self._get_undirected()
        return graph.vs[self.nodDct[uniprot]] \
            if uniprot in self.nodDct else None

    up = uniprot

    def duniprot(self, uniprot):
        """
        Same as ``PyPath.uniprot(), just for directed graph.
        Returns ``igraph.Vertex()`` object if the UniProt
        can be found in the default directed network,
        otherwise ``None``.

        @uniprot : str
            UniProt ID.
        """

        dgraph = self._get_directed()
        return dgraph.vs[self.dnodDct[uniprot]] \
            if uniprot in self.dnodDct else None

    dup = duniprot

    def uniprots(self, uniprots):
        """
        Returns list of ``igraph.Vertex()`` object
        for a list of UniProt IDs omitting those
        could not be found in the default
        undirected graph.
        """

        return filter(lambda v: v is not None, map(self.uniprot, uniprots))

    ups = uniprots

    def duniprots(self, uniprots):
        """
        Returns list of ``igraph.Vertex()`` object
        for a list of UniProt IDs omitting those
        could not be found in the default
        directed graph.
        """

        return filter(lambda v: v is not None, map(self.duniprot, uniprots))

    dups = duniprots

    def get_node(self, identifier):
        """
        Returns ``igraph.Vertex()`` object if the identifier
        is a valid vertex index in the default undirected graph,
        or a UniProt ID or GeneSymbol which can be found in the
        default undirected network, otherwise ``None``.

        @identifier : int, str
            Vertex index (int) or GeneSymbol (str) or UniProt ID (str) or
            ``igraph.Vertex`` object.
        """

        graph = self._get_undirected()

        if isinstance(identifier, igraph.Vertex):
            identifier = identifier['name']

        return graph.vs[identifier] \
            if isinstance(identifier, int) and identifier < graph.vcount() \
            else graph.vs[self.nodDct[identifier]] \
            if identifier in self.nodDct else \
            graph.vs[self.labDct[identifier]] \
            if identifier in self.labDct else None

    # synonyms
    v = get_node
    protein = get_node
    p = get_node

    def get_node_d(self, identifier):
        """
        Same as ``PyPath.get_node``, just for the directed graph.
        Returns ``igraph.Vertex()`` object if the identifier
        is a valid vertex index in the default directed graph,
        or a UniProt ID or GeneSymbol which can be found in the
        default directed network, otherwise ``None``.

        @identifier : int, str
            Vertex index (int) or GeneSymbol (str) or UniProt ID (str) or
            ``igraph.Vertex`` object.
        """

        dgraph = self._get_directed()

        if isinstance(identifier, igraph.Vertex):
            identifier = identifier['name']

        return dgraph.vs[identifier] \
            if isinstance(identifier, int) and identifier < dgraph.vcount() \
            else dgraph.vs[self.dnodDct[identifier]] \
            if identifier in self.dnodDct else \
            dgraph.vs[self.dlabDct[identifier]] \
            if identifier in self.dlabDct else None

    # synonyms
    dv = get_node_d
    dp = get_node_d
    protein = get_node_d

    def get_nodes(self, identifiers):
        """
        """

        return filter(lambda v: v is not None, map(self.get_node, identifiers))

    vs = get_nodes
    proteins = get_nodes
    ps = get_nodes

    def get_nodes_d(self, identifiers):
        """
        """

        return filter(lambda v: v is not None, map(self.get_node_d, identifiers))

    # these are just synonyms
    dvs = get_nodes_d
    dps = get_nodes_d
    dproteins = get_nodes_d

    def up_edge(self, source, target, directed=True):
        """
        Returns ``igraph.Edge`` object if an edge exist between
        the 2 proteins, otherwise ``None``.

        @source : str
            UniProt ID
        @target : str
            UniProt ID
        @directed : bool
            To be passed to igraph.Graph.get_eid()
        """

        v_source = self.uniprot(source)
        v_target = self.uniprot(target)

        if v_source is not None and v_target is not None:
            eid = self.graph.get_eid(
                v_source.index, v_target.index, directed=directed, error=False)

            if eid != -1:
                return self.graph.es[eid]

        return None

    def gs_edge(self, source, target, directed=True):
        """
        Returns ``igraph.Edge`` object if an edge exist between
        the 2 proteins, otherwise ``None``.

        @source : str
            GeneSymbol
        @target : str
            GeneSymbol
        @directed : bool
            To be passed to igraph.Graph.get_eid()
        """

        v_source = self.genesymbol(source)
        v_target = self.genesymbol(target)

        if v_source is not None and v_target is not None:
            eid = self.graph.get_eid(
                v_source.index, v_target.index, directed=directed, error=False)

            if eid != -1:
                return self.graph.es[eid]

        return None

    def get_edge(self, source, target, directed=True):
        """
        Returns ``igraph.Edge`` object if an edge exist between
        the 2 proteins, otherwise ``None``.

        :param int,str source:
            Vertex index or UniProt ID or GeneSymbol or ``igraph.Vertex``
            object.
        :param int,str target:
            Vertex index or UniProt ID or GeneSymbol or ``igraph.Vertex``
            object.
        :param bool directed:
            To be passed to igraph.Graph.get_eid()
        """

        v_source = self.get_node(source) \
            if not self.graph.is_directed() else self.get_node_d(source)
        v_target = self.get_node(target) \
            if not self.graph.is_directed() else self.get_node_d(target)

        if source is not None and target is not None:
            eid = self.graph.get_eid(
                source.index,
                target.index,
                directed=directed,
                error=False
            )

            if eid != -1:
                return self.graph.es[eid]

        return None

    # synonyms
    protein_edge = get_edge

    def get_edges(self, sources, targets, directed=True):
        """
        Returns a generator with all edges between source and target vertices.

        :param iterable sources:
            Source vertex IDs, names, labels, or any iterable yielding
            ``igraph.Vertex`` objects.
        :param iterable targets:
            Target vertec IDs, names, labels, or any iterable yielding
            ``igraph.Vertex`` objects.
        :param bool directed:
            Passed to ``igraph.get_eid()``.
        """

        return (e for e in (
            self.get_edge(s, t, directed)
            for s in sources
            for t in targets)
            if e is not None)

    def _has_directed(self):
        if self._directed is None:
            if self.graph.is_directed():
                self._directed = self.graph
            else:
                self.get_directed()

    def _already_has_directed(self):
        """
        """

        if self._directed is None:

            if self.graph is not None and self.graph.is_directed():
                self._directed = self.graph

            elif self.dgraph is not None and self.dgraph.is_directed():
                self._directed = self.dgraph

    def _get_directed(self):
        """
        Returns the directed instance of the graph.
        If not available, creates one at ``PyPath.dgraph``.
        """

        self._has_directed()
        return self._directed

    # conversion between directed and undirected vertices

    def _get_undirected(self):
        """
        """

        if self._undirected != self.graph and not self.graph.is_directed():
            self._undirected = self.graph

        if self.graph.is_directed():
            self._undirected = None

        return self._undirected

    def up_in_directed(self, uniprot):
        """
        """

        self._has_directed()
        return self.dnodDct[uniprot] if uniprot in self.dnodDct else None

    def up_in_undirected(self, uniprot):
        """
        """

        self._has_undirected()
        return self.nodDct[uniprot] if uniprot in self.nodDct else None

    def gs_in_directed(self, genesymbol):
        """
        """

        self._has_directed()
        return self.dlabDct[genesymbol] if genesymbol in self.dlabDct else None

    def gs_in_undirected(self, genesymbol):
        """
        """

        self._has_undirected()
        return self.labDct[genesymbol] if genesymbol in self.labDct else None

    def in_directed(self, vertex):
        """
        """

        return self.up_in_directed(vertex['name'])

    def in_undirected(self, vertex):
        """
        """

        return self.up_in_undirected(vertex['name'])

    # affects and affected_by

    def _affected_by(self, vertex):
        """
        """

        dgraph = self._get_directed()

        if dgraph != vertex.graph:
            vertex = self.in_directed(vertex['name'])

        vs = vertex.neighbors(mode='IN') \
            if vertex is not None else []

        return _NamedVertexSeq(vs, self.dnodNam, self.dnodLab)

    def _affects(self, vertex):
        """
        """

        dgraph = self._get_directed()

        if dgraph != vertex.graph:
            vertex = self.in_directed(vertex['name'])

        vs = vertex.neighbors(mode='OUT') \
            if vertex is not None else []

        return _NamedVertexSeq(vs, self.dnodNam, self.dnodLab)

    def affected_by(self, identifier):
        vrtx = self.get_node_d(identifier)
        if vrtx is not None:
            return self._affected_by(vrtx)
        return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    def affects(self, identifier):
        """
        """

        vrtx = self.get_node_d(identifier)

        if vrtx is not None:
            return self._affects(vrtx)

        return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    def up_affects(self, uniprot):
        """
        """

        vrtx = self.duniprot(uniprot)

        if vrtx is not None:
            return self._affects(vrtx)

        return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    def up_affected_by(self, uniprot):
        """
        """

        vrtx = self.duniprot(uniprot)

        if vrtx is not None:
            return self._affected_by(vrtx)

        return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    def gs_affects(self, genesymbol):
        """
        """

        vrtx = self.dgenesymbol(genesymbol)

        if vrtx is not None:
            return self._affects(vrtx)

        return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    def gs_affected_by(self, genesymbol):
        """
        """

        vrtx = self.dgenesymbol(genesymbol)

        if vrtx is not None:
            return self._affected_by(vrtx)

        return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    def up_stimulated_by(self, uniprot):
        """
        """

        dgraph = self._get_directed()

        if uniprot in self.dnodDct:
            vid = self.dnodDct[uniprot]
            vs = self.up_affected_by(uniprot)
            return _NamedVertexSeq(
                filter(
                    lambda v: dgraph.es[dgraph.get_eid(v.index, vid)]['dirs'].positive[v['name'], uniprot],
                    vs), self.dnodNam, self.dnodLab)

        else:
            return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    def up_inhibited_by(self, uniprot):
        """
        """

        dgraph = self._get_directed()

        if uniprot in self.dnodDct:
            vid = self.dnodDct[uniprot]
            vs = self.up_affected_by(uniprot)
            return _NamedVertexSeq(
                filter(
                    lambda v: dgraph.es[dgraph.get_eid(v.index, vid)]['dirs'].negative[v['name'], uniprot],
                    vs), self.dnodNam, self.dnodLab)

        else:
            return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    def up_stimulates(self, uniprot):
        """
        """

        dgraph = self._get_directed()

        if uniprot in self.dnodDct:
            vid = self.dnodDct[uniprot]
            vs  = self.up_affects(uniprot)
            return _NamedVertexSeq(
                filter(
                    lambda v: dgraph.es[dgraph.get_eid(vid, v.index)]['dirs'].positive[uniprot, v['name']],
                    vs), self.dnodNam, self.dnodLab)

        else:
            return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    def up_inhibits(self, uniprot):
        """
        """

        dgraph = self._get_directed()

        if uniprot in self.dnodDct:
            vid = self.dnodDct[uniprot]
            vs = self.up_affects(uniprot)
            return _NamedVertexSeq(
                filter(
                    lambda v: dgraph.es[dgraph.get_eid(vid, v.index)]['dirs'].negative[uniprot, v['name']],
                    vs), self.dnodNam, self.dnodLab)

        else:
            return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    def gs_stimulated_by(self, genesymbol):
        """
        """

        dgraph = self._get_directed()
        uniprot = self.dnodNam[self.dlabDct[genesymbol]] \
            if genesymbol in self.dlabDct else None

        return self.up_stimulated_by(uniprot)

    def gs_stimulates(self, genesymbol):
        """
        """

        dgraph = self._get_directed()
        uniprot = self.dnodNam[self.dlabDct[genesymbol]] \
            if genesymbol in self.dlabDct else None

        return self.up_stimulates(uniprot)

    def gs_inhibited_by(self, genesymbol):
        """
        """

        dgraph = self._get_directed()
        uniprot = self.dnodNam[self.dlabDct[genesymbol]] \
            if genesymbol in self.dlabDct else None

        return self.up_inhibited_by(uniprot)

    def gs_inhibits(self, genesymbol):
        """
        """

        dgraph = self._get_directed()
        uniprot = self.dnodNam[self.dlabDct[genesymbol]] \
            if genesymbol in self.dlabDct else None

        return self.up_inhibits(uniprot)

    # neighbors variations

    def up_neighbors(self, uniprot, mode='ALL'):
        """
        """

        vrtx = self.uniprot(uniprot)

        if vrtx is not None:
            return _NamedVertexSeq(
                vrtx.neighbors(mode=mode), self.nodNam, self.nodLab)

        return _NamedVertexSeq([], self.nodNam, self.nodLab)

    def gs_neighbors(self, genesymbol, mode='ALL'):
        """
        """

        vrtx = self.genesymbol(genesymbol)

        if vrtx is not None:
            return _NamedVertexSeq(
                vrtx.neighbors(mode=mode), self.nodNam, self.nodLab)

        return _NamedVertexSeq([], self.nodNam, self.nodLab)

    def neighbors(self, identifier, mode='ALL'):
        """
        """

        vrtx = self.get_node(identifier)

        if vrtx is not None:
            return _NamedVertexSeq(
                vrtx.neighbors(mode=mode), self.nodNam, self.nodLab)

        return _NamedVertexSeq([], self.nodNam, self.nodLab)

    def dneighbors(self, identifier, mode='ALL'):
        """
        """

        vrtx = self.get_node_d(identifier)

        if vrtx is not None:
            return _NamedVertexSeq(
                vrtx.neighbors(mode=mode), self.dnodNam, self.dnodLab)

        return _NamedVertexSeq([], self.dnodNam, self.dnodLab)

    # neighborhood variations:

    def _neighborhood(self, vs, order=1, mode='ALL'):
        """
        """

        return _NamedVertexSeq(
            map(
                lambda vi: self.graph.vs[vi],
                reduce(
                    lambda a, b: a.extend(b),
                    self.graph.neighborhood(
                        vs, order=order, mode=mode), [])),
            self.nodNam,
            self.nodLab)

    def up_neighborhood(self, uniprot, order=1, mode='ALL'):
        """
        """

        if type(uniprots) in common.simpleTypes:
            uniprots = [uniprots]

        vs = self.uniprots(uniprots)
        return self._neighborhood(vs, order=order, mode=mode)

    def gs_neighborhood(self, genesymbols, order=1, mode='ALL'):
        """
        """

        if type(genesymbols) in common.simpleTypes:
            genesymbols = [genesymbols]

        vs = self.genesymbols(genesymbols)
        return self._neighborhood(vs, order=order, mode=mode)

    def neighborhood(self, identifiers, order=1, mode='ALL'):
        """
        """

        if type(identifiers) in common.simpleTypes:
            identifiers = [identifiers]

        vs = self.get_nodes(identifiers)
        return self._neighborhood(vs, order=order, mode=mode)

    # complexes

    def complexes_in_network(self, csource='corum', graph=None):
        """
        """

        graph = self.graph if graph is None else graph
        cdict = {}
        allv = set(graph.vs['name'])

        for v in graph.vs:

            for c, cdata in iteritems(v['complexes'][csource]):

                if c not in cdict:
                    cdict[c] = set(cdata['all_members'])

        return [c for c, memb in iteritems(cdict) if len(memb - allv) == 0]

    def complex_comembership_network(self, graph=None, resources=None):
        """
        """

        graph = graph if graph is not None else self.graph

        if resources is None:
            resources = graph.vs[0]['complexes'].keys()

        cnet = igraph.Graph()
        cdict = {}

        for v in graph.vs:

            for src in resources:

                if src not in cdict:
                    cdict[src] = {}

                if src in c['complexes']:

                    for cname, cannot in iteritems(c['complexes'][src]):

                        if cname not in cdict[src]:
                            cdict[src][cname] = []

                        cdict[src][cname].append(v['name'])

        pass

    def edges_in_comlexes(self, csources=['corum'], graph=None):
        """
        Creates edge attributes ``complexes`` and ``in_complex``.
        These are both dicts where the keys are complex resources.
        The values in ``complexes`` are the list of complex names
        both the source and the target vertices belong to.
        The values ``in_complex`` are boolean values whether there
        is at least one complex in the given resources both the
        source and the target vertex of the edge belong to.

        @csources : list
            List of complex resources. Should be already loaded.
        @graph : igraph.Graph()
            The graph object to do the calculations on.
        """

        graph = graph if graph is not None else self.graph

        if 'complexes' not in graph.es.attributes():
            graph.es['complexes'] = [{} for _ in graph.es]

        if 'in_complex' not in graph.es.attributes():
            graph.es['in_complex'] = [{} for _ in graph.es]

        def _common_complexes(e, cs):
            return set(graph.vs[e.source]['complexes'][cs].keys()) & \
                set(graph.vs[e.source]['complexes'][cs].keys())

        nul = map(lambda cs:
                  map(lambda e:
                      e['complexes'].__setitem__(cs, _common_complexes(e, cs)),
                      graph.es),
                  csources)

        nul = map(lambda cs:
                  map(lambda e:
                      e['in_complex'].__setitem__(
                          cs, bool(len(e['complexes'][cs]))),
                      graph.es),
                  csources)

    def sum_in_complex(self, csources=['corum'], graph=None):
        """
        Returns the total number of edges in the network falling
        between two members of the same complex.
        Returns as a dict by complex resources.
        Calls :py:func:pypath.pypath.Pypath.edges_in_comlexes()
        to do the calculations.

        @csources : list
            List of complex resources. Should be already loaded.
        @graph : igraph.Graph()
            The graph object to do the calculations on.
        """

        graph = graph if graph is not None else self.graph
        self.edges_in_comlexes(csources=csources, graph=graph)

        return dict(map(lambda cs:
                        (cs, sum(map(lambda e:
                                     e['in_complex'][cs], graph.es)
                                 )),
                        csources)
                    )

    #
    #
    #

    def get_function(self, fun):
        """
        """

        if hasattr(fun, '__call__'):
            return fun

        fun = fun.split('.')
        fun0 = fun.pop(0)
        toCall = None

        if fun0 in globals():
            toCall = globals()[fun0]

        elif fun0 in locals():
            toCall = locals()[fun0]

        elif fun0 == __name__.split('.')[0]:
            toCall = __main__

        elif fun0 == 'self' or fun0 == __name__.split('.')[-1]:
            toCall = self

        elif fun0 in dir(self):
            toCall = getattr(self, fun0)

        else:

            for subm in globals().keys():

                if hasattr(globals()[subm], '__name__') and globals(
                )[subm].__name__.split('.')[0] == __name__.split('.')[-1]:

                    if hasattr(globals()[subm], fun0):
                        toCall = getattr(globals()[subm], fun0)

        if toCall is None:
            return None

        for fun0 in fun:

            if hasattr(toCall, fun0):
                toCall = getattr(toCall, fun0)

        if hasattr(toCall, '__call__'):
            return toCall

        else:
            return None

    def edges_3d(self, methods=['dataio.get_instruct', 'dataio.get_i3d']):
        """
        """

        all_3d = []
        self.update_vname()

        for m in methods:
            fun = self.get_function(m)

            if fun is not None:
                all_3d += fun()

        # initializing edge attributes:
        if 'ddi' not in self.graph.es.attributes():
            self.graph.es['ddi'] = [None]

        for e in self.graph.es:

            if e['ddi'] is None:
                e['ddi'] = []

        # assigning annotations:
        for i in all_3d:
            u1 = i['uniprots'][0]
            u2 = i['uniprots'][1]

            if u1 != u2 and self.node_exists(u1) and self.node_exists(u2):
                edge = self.edge_exists(u1, u2)

                if isinstance(edge, int):

                    for seq1 in i[u1]['seq']:

                        for seq2 in i[u2]['seq']:
                            dom1 = intera.Domain(
                                protein=u1,
                                domain=i[u1]['pfam'],
                                start=seq1[0],
                                end=seq1[1],
                                chains={i['pdb'][0]: i[u1]['chain']})
                            dom2 = intera.Domain(
                                protein=u2,
                                domain=i[u2]['pfam'],
                                start=seq2[0],
                                end=seq2[1],
                                chains={i['pdb'][0]: i[u2]['chain']})
                            domdom = intera.DomainDomain(
                                dom1,
                                dom2,
                                refs=i['references'],
                                sources=i['source'],
                                pdbs=i['pdb'])
                            self.graph.es[edge]['ddi'].append(domdom)

    def edges_ptms(self):
        """
        """

        self.pfam_regions()
        self.load_pepcyber()
        self.load_psite_reg()
        self.load_ielm()
        self.load_phosphoelm()
        self.load_elm()
        self.load_3did_dmi()

    def pfam_regions(self):
        """
        """

        if self.u_pfam is None:
            self.u_pfam = dataio.get_pfam_regions(
                uniprots=self.graph.vs['name'], dicts='uniprot', keepfile=True)

    def complexes(self, methods=['3dcomplexes', 'havugimana', 'corum',
                                 'complexportal', 'compleat']):
        """
        """

        for m in methods:
            m = 'load_' + m

            if hasattr(self, m):
                toCall = getattr(self, m)

                if hasattr(toCall, '__call__'):
                    toCall()

    def load_domino_dmi(self, organism=None):
        """
        """

        organism = organism if organism is not None else self.ncbi_tax_id
        #all_unip = dataio.uniprot_input.all_uniprots(organism)
        domi = dataio.get_domino_ptms()

        if domi is None:
            self.ownlog.msg(
                2, 'Failed to load domain-motif interaction data from DOMINO',
                'ERROR')
            return None

        self.update_vname()

        if 'ptm' not in self.graph.es.attributes():
            self.graph.es['ptm'] = [[] for _ in self.graph.es]

        prg = Progress(
            len(domi['dmi']), 'Loading domain-motif interactions', 11)

        for dm in domi['dmi']:
            prg.step()
            uniprot1 = dm.domain.protein
            uniprot2 = dm.ptm.protein

            if self.node_exists(uniprot1) and self.node_exists(uniprot2):
                e = self.edge_exists(uniprot1, uniprot2)

                if isinstance(e, int):

                    if isinstance(e, int):

                        if not isinstance(self.graph.es[e]['ptm'], list):
                            self.graph.es[e]['ptm'] = []

                        self.graph.es[e]['ptm'].append(dm)

        prg.terminate()

    def load_3did_dmi(self):
        """
        """

        dmi = dataio.process_3did_dmi()

        if dmi is None:
            self.ownlog.msg(
                2, 'Failed to load domain-motif interaction data from 3DID',
                'ERROR')
            return None

        self.update_vname()

        if 'ptm' not in self.graph.es.attributes():
            self.graph.es['ptm'] = [[] for _ in self.graph.es]

        prg = Progress(len(dmi), 'Loading domain-motif interactions', 11)

        for uniprots, dmi_list in iteritems(dmi):
            prg.step()

            if self.node_exists(uniprots[0]) and self.node_exists(uniprots[1]):
                e = self.edge_exists(uniprots[0], uniprots[1])

                if isinstance(e, int):

                    if not isinstance(self.graph.es[e]['ptm'], list):
                        self.graph.es[e]['ptm'] = []

                    self.graph.es[e]['ptm'] += dmi_list

        prg.terminate()

    def load_ddi(self, ddi):
        """
        ddi is either a list of intera.DomainDomain objects,
        or a function resulting this list
        """

        data = ddi if not hasattr(ddi, '__call__') else ddi()

        if data is None:

            if ddi.__module__.split('.')[1] == 'dataio':
                self.ownlog.msg(2, 'Function %s() failed' % ddi, 'ERROR')

            return None

        if 'ddi' not in self.graph.es.attributes():
            self.graph.es['ddi'] = [[] for _ in self.graph.es]

        prg = Progress(len(data), 'Loading domain-domain interactions', 99)
        in_network = 0

        for dd in data:
            prg.step()
            uniprot1 = dd.domains[0].protein
            uniprot2 = dd.domains[1].protein

            if self.node_exists(uniprot1) and self.node_exists(uniprot2):
                e = self.edge_exists(uniprot1, uniprot2)

                if isinstance(e, int):

                    if not isinstance(self.graph.es[e]['ddi'], list):
                        self.graph.es[e]['ddi'] = []

                    in_network += 1
                    self.graph.es[e]['ddi'].append(dd)

        prg.terminate()

    def load_dmi(self, dmi):
        """
        dmi is either a list of intera.DomainMotif objects,
        or a function resulting this list
        """

        data = dmi if not hasattr(dmi, '__call__') else dmi()

        if data is None:

            if dmi.__module__.split('.')[1] == 'dataio':
                self.ownlog.msg(2, 'Function %s() failed' % dmi, 'ERROR')

            return None

        if 'ptm' not in self.graph.es.attributes():
            self.graph.es['ptm'] = [[] for _ in self.graph.es]

    def run_batch(self, methods, toCall=None):
        """
        """

        if toCall is not None:
            toCall = self.get_function(toCall)

        for m in methods:
            fun = self.get_function(m)

            if fun is not None:

                if hasattr(toCall, '__call__'):
                    toCall(fun)

                else:
                    fun()

    def load_ddis(self, methods=['dataio.get_3dc_ddi', 'dataio.get_domino_ddi',
                                 'self.load_3did_ddi2']):
        """
        """

        self.run_batch(methods, toCall=self.load_ddi)

    def load_dmis(self, methods=['self.pfam_regions', 'self.load_depod_dmi',
                                 'self.load_dbptm', 'self.load_mimp_dmi',
                                 'self.load_pnetworks_dmi',
                                 'self.load_domino_dmi', 'self.load_pepcyber',
                                 'self.load_psite_reg', 'self.load_psite_phos',
                                 'self.load_ielm', 'self.load_phosphoelm',
                                 'self.load_elm', 'self.load_3did_dmi']):
        """
        """

        self.run_batch(methods)
        self.uniq_ptms()
        self.phosphorylation_directions()

    def load_interfaces(self):
        """
        """

        self.load_3did_ddi2(ddi=False, interfaces=True)
        unm = self.load_pisa()

    def load_3did_ddi(self):
        """
        """

        g = self.graph
        ddi, interfaces = dataio.get_3did_ddi(residues=True)

        if ddi is None or interfaces is None:
            self.ownlog.msg(
                2, 'Failed to load domain-domain interaction data from 3DID',
                'ERROR')
            return None

        if 'ddi' not in g.es.attributes():
            g.es['ddi'] = [None]

        for e in g.es:

            if e['ddi'] is None:
                e['ddi'] = {}

            if '3did' not in e['ddi']:
                e['ddi']['3did'] = []

        prg = Progress(len(ddi), 'Loading domain-domain interactions', 99)

        for k, v in iteritems(ddi):
            uniprot1 = k[0]
            uniprot2 = k[1]
            swprots = self.mapper.swissprots([uniprot1, uniprot2])

            for swprot1 in swprots[uniprot1]:

                for swprot2 in swprots[uniprot2]:

                    if swprot1 in g.vs['name'] and swprot2 in g.vs['name']:
                        e = self.edge_exists(swprot1, swprot2)

                        if isinstance(e, int):

                            for domains, structures in iteritems(v):

                                for pdb, pdb_uniprot_pairs in \
                                        iteritems(structures['pdbs']):

                                    for pdbuniprots in pdb_uniprot_pairs:
                                        pdbswprots = self.mapper.swissprots(
                                            pdbuniprots)

                                        for pdbswprot1 in pdbswprots[
                                                pdbuniprots[0]]:

                                            for pdbswprot2 in pdbswprots[
                                                    pdbuniprots[1]]:
                                                this_ddi = {}
                                                this_ddi[swprot1] = {
                                                    'pfam': domains[0]
                                                }
                                                this_ddi[swprot2] = {
                                                    'pfam': domains[1]
                                                }
                                                this_ddi['uniprots'] = [
                                                    swprot1, swprot2
                                                ]
                                                this_ddi[
                                                    'uniprots_original'] = k
                                                this_ddi['pdb'] = {
                                                    swprot1: pdbswprot1,
                                                    swprot2: pdbswprot2,
                                                    'pdb': pdb
                                                }
                                                g.es[e]['ddi']['3did'].append(
                                                    this_ddi)

            prg.step()

        prg.terminate()

    def load_3did_interfaces(self):
        """
        """

        self.load_3did_ddi2(ddi=False, interfaces=True)

    def load_3did_ddi2(self, ddi=True, interfaces=False):
        """
        """

        self.update_vname()
        ddi, intfs = dataio.get_3did()

        if ddi is None or interfaces is None:
            self.ownlog.msg(
                2, 'Failed to load domain-domain interaction data from 3DID',
                'ERROR')
            return None

        # ERROR
        if interfaces:

            if 'interfaces' not in self.graph.es.attributes():
                self.graph.es['interfaces'] = [[] for _ in self.graph.es]

            for intf in intfs:
                uniprot1 = intf.domains[0].protein
                uniprot2 = intf.domains[1].protein

                if self.node_exists(uniprot1) and self.node_exists(uniprot2):
                    e = self.edge_exists(uniprot1, uniprot2)

                    if isinstance(e, int):

                        if not isinstance(self.graph.es[e]['interfaces'],
                                          list):
                            self.graph.es[e]['interfaces'] = []

                        self.graph.es[e]['interfaces'].append(intf)

        if ddi:
            return ddi

    def load_ielm(self):
        """
        """

        self.update_vname()

        if 'ptm' not in self.graph.es.attributes():
            self.graph.es['ptm'] = [[] for _ in self.graph.es]

        ppi = []

        for e in self.graph.es:
            ppi.append((self.graph.vs[e.source]['name'],
                        self.graph.vs[e.target]['name']))

        ielm = dataio.get_ielm(ppi)
        elm = dataio.get_elm_classes()
        prg = Progress(len(ielm), 'Processing domain-motif interactions', 13)
        noelm = []

        for l in ielm:
            prg.step()
            nodes = self.get_node_pair(l[1], l[9])

            if nodes:
                e = self.graph.get_eid(nodes[0], nodes[1], error=False)

                if e != -1:

                    if self.graph.es[e]['ptm'] is None:
                        self.graph.es[e]['ptm'] = []

                    if l[2] not in elm:
                        noelm.append(l[2])

                    motif = [None] * 5 if l[2] not in elm else elm[l[2]]
                    mot = intera.Motif(
                        l[1],
                        int(l[3]),
                        int(l[4]),
                        regex=motif[3],
                        prob=None if motif[4] is None else float(motif[4]),
                        description=motif[2],
                        elm=motif[0],
                        motif_name=l[2])
                    dom = intera.Domain(
                        protein=l[9],
                        start=int(l[11]),
                        end=int(l[12]),
                        domain=l[10],
                        domain_id_type='ielm_hmm')

                    if self.u_pfam is not None and l[9] in self.u_pfam:

                        for pfam, regions in iteritems(self.u_pfam[l[9]]):

                            for region in regions:

                                if int(l[11]) == region['start'] and \
                                        int(l[12]) == region['end']:
                                    dom.pfam = pfam
                                    dom.isoform = region['isoform']

                    ptm = intera.Ptm(protein=l[1],
                                     motif=mot,
                                     typ=l[2].split('_')[0])
                    domMot = intera.DomainMotif(dom, ptm, 'iELM')
                    self.graph.es[e]['ptm'].append(domMot)

        prg.terminate()
        # find deprecated elm classes:
        # list(set(noelm))

    def load_elm(self):
        """
        """

        self.update_vname()

        if 'ptm' not in self.graph.es.attributes():
            self.graph.es['ptm'] = [[] for _ in self.graph.es]

        elm = dataio.get_elm_interactions()
        prg = Progress(len(elm), 'Processing domain-motif interactions', 11)

        self.sequences()

        for l in elm:
            prg.step()

            if len(l) > 7:
                uniprot_elm = l[2]
                uniprot_dom = l[3]

                if self.node_exists(uniprot_elm) and self.node_exists(
                        uniprot_dom):
                    nodes = self.get_node_pair(uniprot_dom, uniprot_elm)

                    if nodes:
                        e = self.graph.get_eid(nodes[0], nodes[1], error=False)

                        if e != -1:

                            if self.graph.es[e]['ptm'] is None:
                                self.graph.es[e]['ptm'] = []

                            start = int(l[4])
                            end   = int(l[5])

                            if uniprot_elm in self.seq:
                                inst = self.seq[uniprot_elm].get_region(start = start,
                                                                        end = end)[2]

                            else:
                                inst = None

                            mot = intera.Motif(
                                protein=uniprot_elm,
                                motif_name=l[0],
                                start=start,
                                end=end,
                                instance = inst)
                            ptm = intera.Ptm(protein=uniprot_elm, motif=mot)
                            dstart = start = None if l[6] == 'None' else int(l[
                                6])
                            dend = None if l[6] == 'None' else int(l[7])

                            doms = []

                            if self.u_pfam is not None and uniprot_dom in self.u_pfam:

                                if l[1] in self.u_pfam[uniprot_dom]:

                                    for region in self.u_pfam[uniprot_dom][l[
                                            1]]:

                                        if dstart is None or dend is None or \
                                            (dstart == region['start'] and
                                             dend == region['end']):
                                            doms.append(
                                                intera.Domain(
                                                    protein=uniprot_dom,
                                                    domain=l[1],
                                                    start=region['start'],
                                                    end=region['end'],
                                                    isoform=region['isoform']))

                            else:
                                doms = [
                                    intera.Domain(
                                        protein=uniprot_dom,
                                        domain=l[1]
                                    )
                                ]

                            for dom in doms:
                                self.graph.es[e]['ptm'].append(
                                    intera.DomainMotif(dom, ptm, 'ELM'))

        prg.terminate()

    def load_lmpid(self, method):
        """
        """

        pass

    def process_dmi(self, source, **kwargs):
        """
        This is an universal function
        for loading domain-motif objects
        like load_phospho_dmi() for PTMs.
        TODO this will replace load_elm, load_ielm, etc
        """

        functions = {'LMPID': 'lmpid_dmi'}
        motif_plus = {'LMPID': []}
        self.update_vname()
        toCall = getattr(dataio, functions[source])
        data = toCall(**kwargs)

        if self.seq is None:
            self.sequences(self.ncbi_tax_id)

        if self.seq is None:
            self.sequences(isoforms=True)

        if 'ptm' not in self.graph.es.attributes():
            self.graph.es['ptm'] = [[] for _ in self.graph.es]

        prg = Progress(
            len(data), 'Processing domain-motif interactions from %s' % source,
            7)

        for d in data:
            prg.step()
            domain_ups = []
            motif_upd = []

            if source == 'LMPID':
                domain_ups = self.mapper.map_name(d['domain_protein'],
                                                  'uniprot', 'uniprot')
                motif_ups = self.mapper.map_name(d['motif_protein'], 'uniprot',
                                                 'uniprot')
            for du in domain_ups:
                dom = intera.Domain(
                    du,
                    domain=None
                    if 'domain_name' not in d else d['domain_name'],
                    domain_id_type=None
                    if 'domain_name_type' not in d else d['domain_name_type'])

                for mu in motif_ups:

                    if mu in self.seq and mu in self.nodInd and du in self.nodInd:
                        edge = self._get_edge(
                            (self.nodDct[du], self.nodDct[mu]))

                        if edge:
                            mse = self.seq[mu]
                            isos = []

                            if d['instance'] is None:
                                start, end, inst = mse.get_region(
                                    start=d['motif_start'],
                                    end=d['motif_end'],
                                    isoform=1)
                                isos.append((start, end, 1, inst))

                            for iso in mse.isoforms():
                                start, end, inst = mse.get_region(
                                    start=d['motif_start'],
                                    end=d['motif_end'],
                                    isoform=iso)

                                if inst == d['instance']:
                                    isos.append((start, end, iso, inst))

                            for iso in isos:
                                motargs = dict([(k, d[k])
                                                for k in motif_plus[source]])
                                mot = intera.Motif(
                                    mu,
                                    iso[0],
                                    iso[1],
                                    isoform=iso[2],
                                    instance=iso[3],
                                    **motargs)
                                ptm = intera.Ptm(mu,
                                                 motif=mot,
                                                 isoform=iso[2],
                                                 source=[source])
                                dommot = intera.DomainMotif(
                                    domain=dom,
                                    ptm=ptm,
                                    sources=[source],
                                    refs=d['refs'])
                                self.graph.es[edge]['ptm'].append(dommot)

        prg.terminate()

    def load_pepcyber(self):
        """
        """

        self.update_vname()

        if 'ptm' not in self.graph.es.attributes():
            self.graph.es['ptm'] = [[] for _ in self.graph.es]

        pepc = dataio.get_pepcyber()
        prg = Progress(len(pepc), 'Processing domain-motif interactions', 13)

        for l in pepc:
            prg.step()
            uniprot1 = [l[9]] if len(l[9]) > 0 else []

            if len(l[9]) == 0 and len(l[10]) > 0:
                uniprot1 = self.mapper.map_name(l[10],
                                                'refseq',
                                                'uniprot',
                                                9606)

            uniprot2 = [l[11]] if len(l[11]) > 0 else []

            # ptm on u2,
            # u1 interacts with u2 depending on ptm
            for u1 in uniprot1:

                for u2 in uniprot2:
                    nodes = self.get_node_pair(u1, u2)

                    if nodes:
                        e = self.graph.get_eid(nodes[0], nodes[1], error=False)

                        if e != -1:
                            res = intera.Residue(int(l[4]), l[8], u2)
                            ptm = intera.Ptm(protein=u2,
                                             residue=res,
                                             typ='phosphorylation')
                            reg = intera.Regulation(
                                source=u2,
                                target=u1,
                                sources=['PepCyber'],
                                effect='induces',
                                ptm=ptm)
                            self.graph.es[e]['ptm'].append(reg)

        prg.terminate()

    def load_psite_reg(self):
        """
        """

        modtyp = {
            'p': 'phosphorylation',
            'ac': 'acetylation',
            'ub': 'ubiquitination',
            'me': 'methylation',
            'sm': 'sumoylation',
            'ga': 'o-galnac',
            'gl': 'o-glcnac',
            'sc': 'succinylation',
            'm3': 'tri-methylation',
            'm1': 'mono-methylation',
            'm2': 'di-methylation',
            'pa': 'palmitoylation'
        }
        self.update_vname()

        if 'ptm' not in self.graph.es.attributes():
            self.graph.es['ptm'] = [[] for _ in self.graph.es]

        preg = dataio.get_psite_reg()
        prg = Progress(len(preg), 'Processing regulatory effects', 11)

        for src, tgts in iteritems(preg):
            prg.step()

            # ptm on src
            # tgt: interactor of src, depending on ptm
            if self.node_exists(src):

                for tgt in tgts:

                    for effect in ['induces', 'disrupts']:

                        for ind in tgt[effect]:
                            uniprots = self.mapper.map_name(ind, 'genesymbol',
                                                            'uniprot')

                            for u in uniprots:

                                if self.node_exists(u):
                                    nodes = self.get_node_pair(src, u)

                                    if nodes:
                                        e = self.graph.get_eid(
                                            nodes[0], nodes[1], error=False)

                                        if e != -1:

                                            if self.graph.es[e]['ptm'] is None:
                                                self.graph.es[e]['ptm'] = []

                                            res = intera.Residue(
                                                int(tgt['res']), tgt['aa'],
                                                src)
                                            ptm = intera.Ptm(
                                                protein=src,
                                                residue=res,
                                                typ=modtyp[tgt['modt']])
                                            reg = intera.Regulation(
                                                source=src,
                                                target=u,
                                                sources=['PhosphoSite'],
                                                refs=tgt['pmids'],
                                                effect=effect,
                                                ptm=ptm)
                                            self.graph.es[e]['ptm'].append(reg)

        prg.terminate()

    def load_comppi(self, graph=None):
        """
        """

        graph = graph if graph is not None else self.graph
        self.update_vname()

        if 'comppi' not in graph.es.attributes():
            graph.es['comppi'] = [None for _ in graph.es]

        if 'comppi' not in graph.vs.attributes():
            graph.vs['comppi'] = [{} for _ in graph.vs]

        comppi = dataio.get_comppi()
        prg = Progress(len(comppi), 'Processing localizations', 33)

        for c in comppi:
            prg.step()
            uniprots1 = self.mapper.map_name(c['uniprot1'], 'uniprot',
                                             'uniprot', 9606)
            uniprots2 = self.mapper.map_name(c['uniprot2'], 'uniprot',
                                             'uniprot', 9606)

            for u1 in uniprots1:

                if self.node_exists(u1):

                    for loc in [x.split(':') for x in c['loc1'].split('|')]:
                        graph.vs[self.nodDct[u1]]\
                            ['comppi'][loc[0]] = float(loc[1])

            for u2 in uniprots1:

                if self.node_exists(u2):

                    for loc in [x.split(':') for x in c['loc2'].split('|')]:
                        graph.vs[self.nodDct[u2]]\
                            ['comppi'][loc[0]] = float(loc[1])

            for u1 in uniprots1:

                for u2 in uniprots2:

                    if self.node_exists(u1) and self.node_exists(u2):
                        nodes = self.get_node_pair(u1, u2)

                        if nodes:
                            e = graph.get_eid(nodes[0], nodes[1], error=False)

                            if e != -1:
                                graph.es[e]['comppi'] = float(c['loc_score'])

        prg.terminate()

    def edge_loc(self, graph=None, topn=2):
        """
        """

        graph = graph if graph is not None else self.graph

        if 'comppi' not in graph.vs.attributes():
            self.load_comppi()

        graph.es['loc'] = [[] for _ in graph.es]

        for e in graph.es:
            e['loc'] = list(
                set([
                    i[0]
                    for i in heapq.nlargest(2,
                                            iteritems(graph.vs[e.source][
                                                'comppi']),
                                            operator.itemgetter(1))
                ]) & set([
                    i[0]
                    for i in heapq.nlargest(2,
                                            iteritems(graph.vs[e.target][
                                                'comppi']),
                                            operator.itemgetter(1))
                ]))

    def sequences(self, isoforms=True, update=False):
        """
        """

        if self.seq is None or update:
            self.seq = se.swissprot_seq(self.ncbi_tax_id, isoforms)

    def load_ptms2(self, input_methods=None, map_by_homology_from=[9606],
                   homology_only_swissprot=True, ptm_homology_strict=False,
                   nonhuman_direct_lookup=True, inputargs={}):
        """
        This is a new method which will replace `load_ptms`.
        It uses `pypath.ptm.PtmAggregator`, a newly introduced
        module for combining enzyme-substrate data from multiple
        resources using homology translation on users demand.

        :param list input_methods: Resources to collect enzyme-substrate
            interactions from. E.g. `['Signor', 'phosphoELM']`. By default
            it contains Signor, PhosphoSitePlus, HPRD, phosphoELM, dbPTM,
            PhosphoNetworks, Li2012 and MIMP.
        :param list map_by_homology_from: List of NCBI Taxonomy IDs of
            source taxons used for homology translation of enzyme-substrate
            interactions. If you have a human network and you add here
            `[10090, 10116]` then mouse and rat interactions from the source
            databases will be translated to human.
        :param bool homology_only_swissprot: `True` by default which means
            only SwissProt IDs are accepted at homology translateion, Trembl
            IDs will be dropped.
        :param bool ptm_homology_strict: For homology translation use
            PhosphoSite's PTM homology table. This guarantees that only
            truely homologous sites will be included. Otherwise we only
            check if at the same numeric offset in the homologous sequence
            the appropriate residue can be find.
        :param bool nonhuman_direct_lookup: Fetch also directly nonhuman
            data from the resources whereever it's available. PhosphoSite
            contains mouse enzyme-substrate interactions and it is possible
            to extract these directly beside translating the human ones
            to mouse.
        :param dict inputargs: Additional arguments passed to `PtmProcessor`.
            A `dict` can be supplied for each resource, e.g.
            `{'Signor': {...}, 'PhosphoSite': {...}, ...}`.
            Those not used by `PtmProcessor` are forwarded to the
            `pypath.dataio` methods.
        """

        ptma = pypath.ptm.PtmAggregator(
            input_methods = input_methods,
            ncbi_tax_id = self.ncbi_tax_id,
            map_by_homology_from = map_by_homology_from,
            # here we don't share the mapper as later many
            # tables which we don't need any more would
            # just occupy memory
            mapper = self.mapper,
            homology_only_swissprot = homology_only_swissprot,
            ptm_homology_strict = ptm_homology_strict,
            nonhuman_direct_lookup = nonhuman_direct_lookup,
            inputargs = inputargs
        )
        ptma.assign_to_network(self)

        if self.ncbi_tax_id == 9606 and (
                input_methods is None or
                'DEPOD' in input_methods
            ):
            self.load_depod_dmi()

        self.uniq_ptms()
        self.phosphorylation_directions()

    def load_ptms(self):
        """
        """

        self.load_depod_dmi()
        self.load_signor_ptms()
        self.load_li2012_ptms()
        self.load_hprd_ptms()
        self.load_mimp_dmi()
        self.load_pnetworks_dmi()
        self.load_phosphoelm()
        self.load_dbptm()
        self.load_psite_phos()
        self.uniq_ptms()
        self.phosphorylation_directions()

    def load_mimp_dmi(self, non_matching=False, trace=False, **kwargs):
        """
        """

        trace = self.load_phospho_dmi(source='MIMP', trace=trace, **kwargs)

        if trace:
            return trace

    def load_pnetworks_dmi(self, trace=False, **kwargs):
        """
        """

        trace = self.load_phospho_dmi(
            source='PhosphoNetworks', trace=trace, **kwargs)

        if trace:
            return trace

    def load_phosphoelm(self, trace=False, **kwargs):
        """
        """

        trace = self.load_phospho_dmi(
            source='phosphoELM', trace=trace, **kwargs)

        if trace:
            return trace

    def load_dbptm(self, non_matching=False, trace=False, **kwargs):
        """
        """

        trace = self.load_phospho_dmi(source='dbPTM', trace=trace, **kwargs)

        if trace:
            return trace

    def load_hprd_ptms(self, non_matching=False, trace=False, **kwargs):
        """
        """

        trace = self.load_phospho_dmi(source='HPRD', trace=trace, **kwargs)

        if trace:
            return trace

    def load_li2012_ptms(self, non_matching=False, trace=False, **kwargs):
        """
        """

        trace = self.load_phospho_dmi(source='Li2012', trace=trace, **kwargs)

        if trace:
            return trace

    def load_signor_ptms(self, non_matching=False, trace=False, **kwargs):
        """
        """

        trace = self.load_phospho_dmi(source='Signor', trace=trace, **kwargs)

        if trace:
            return trace

    def load_psite_phos(self, trace=False, **kwargs):
        """
        """

        trace = self.load_phospho_dmi(
            source='PhosphoSite', trace=trace, **kwargs)

        if trace:
            return trace

    def load_phospho_dmi(self, source, trace=False, return_raw=False,
                         **kwargs):
        """
        """

        functions = {
            'Signor': 'load_signor_ptms',
            'MIMP': 'get_mimp',
            'PhosphoNetworks': 'get_phosphonetworks',
            'phosphoELM': 'get_phosphoelm',
            'dbPTM': 'get_dbptm',
            'PhosphoSite': 'get_psite_phos',
            'HPRD': 'get_hprd_ptms',
            'Li2012': 'li2012_phospho'
        }

        if source == 'PhosphoSite' and self.ncbi_tax_id in common.taxids:
            kwargs['organism'] = common.taxids[self.ncbi_tax_id]

            if 'strict' not in kwargs:
                kwargs['strict'] = False
                kwargs['mapper'] = self.mapper

        if source == 'Signor':
            kwargs['organism'] = self.ncbi_tax_id

        if source == 'phosphoELM':
            kwargs['organism'] = self.ncbi_tax_id

            if self.ncbi_tax_id != 9606 and 'ltp_only' not in kwargs:
                kwargs['ltp_only'] = False

        if source == 'dbPTM':
            kwargs['organism'] = self.ncbi_tax_id

        if self.ncbi_tax_id != 9606 and source not in \
            ['Signor', 'PhosphoSite', 'phosphoELM', 'dbPTM']:
            return None

        self.update_vname()
        toCall = getattr(dataio, functions[source])
        data = toCall(**kwargs)

        if self.seq is None:
            self.sequences(isoforms=True)

        if 'ptm' not in self.graph.es.attributes():
            self.graph.es['ptm'] = [[] for _ in self.graph.es]

        nomatch = []
        kin_ambig = {}
        sub_ambig = {}
        prg = Progress(len(data), 'Processing PTMs from %s' % source, 23)
        raw = []

        for p in data:
            prg.step()

            if p['kinase'] is not None and len(p['kinase']) > 0:

                # database specific id conversions
                if source in ['PhosphoSite', 'phosphoELM', 'Signor']:
                    kinase_ups = self.mapper.map_name(p['kinase'], 'uniprot',
                                                      'uniprot')

                else:

                    if not isinstance(p['kinase'], list):
                        p['kinase'] = [p['kinase']]
                    kinase_ups = [
                        i
                        for ii in [
                            self.mapper.map_name(k, 'genesymbol', 'uniprot')
                            for k in p['kinase']
                        ] for i in ii
                    ]

                    kinase_ups = list(set(kinase_ups))

                if not isinstance(p['kinase'], list):
                    p['kinase'] = [p['kinase']]

                if p['substrate'].startswith('HLA'):
                    continue

                if source in ['MIMP', 'PhosphoNetworks', 'Li2012']:
                    substrate_ups_all = self.mapper.map_name(
                        p['substrate'], 'genesymbol', 'uniprot')

                if source == 'MIMP':
                    substrate_ups_all += self.mapper.map_name(
                        p['substrate_refseq'], 'refseq', 'uniprot')
                    substrate_ups_all = list(set(substrate_ups_all))

                if source in ['phosphoELM', 'dbPTM', 'PhosphoSite', 'Signor']:
                    substrate_ups_all = self.mapper.map_name(
                        p['substrate'], 'uniprot', 'uniprot')

                if source == 'HPRD':
                    substrate_ups_all = self.mapper.map_name(
                        p['substrate_refseqp'], 'refseqp', 'uniprot')
                    substrate_ups_all = common.uniqList(substrate_ups_all)

                for k in p['kinase']:

                    if k not in kin_ambig:
                        kin_ambig[k] = kinase_ups

                # looking up sequences in all isoforms:
                substrate_ups = []

                for s in substrate_ups_all:

                    if s in self.seq:

                        for isof in self.seq[s].isoforms():

                            if p['instance'] is not None:

                                if self.seq[s].match(
                                        p['instance'],
                                        p['start'],
                                        p['end'],
                                        isoform=isof):
                                    substrate_ups.append((s, isof))

                            else:

                                if self.seq[s].match(
                                        p['resaa'], p['resnum'], isoform=isof):
                                    substrate_ups.append((s, isof))

                if p['substrate'] not in sub_ambig:
                    sub_ambig[p['substrate']] = substrate_ups

                # generating report on non matching substrates
                if len(substrate_ups) == 0:

                    for s in substrate_ups_all:

                        if s[0] in self.seq:
                            nomatch.append((s[0], s[1], p['substrate_refseq'],
                                            s, p['instance'], self.seq[s].get(
                                                p['start'], p['end'])))

                # adding kinase-substrate interactions
                for k in kinase_ups:

                    for s in substrate_ups:
                        nodes = self.get_node_pair(k, s[0],
                            directed = self.graph.is_directed())

                        if nodes or return_raw:
                            e = None

                            if nodes:
                                e = self.graph.get_eid(
                                    nodes[0], nodes[1], error=False)

                            if (isinstance(e, int) and e > 0) or return_raw:
                                res = intera.Residue(
                                    p['resnum'],
                                    p['resaa'],
                                    s[0],
                                    isoform=s[1])

                                if p['instance'] is None:
                                    reg = self.seq[s[0]].get_region(
                                        p['resnum'],
                                        p['start'],
                                        p['end'],
                                        isoform=s[1])

                                    if reg is not None:
                                        p['instance'] = reg[2]
                                        p['start'] = reg[0]
                                        p['end'] = reg[1]

                                if 'typ' not in p:
                                    p['typ'] = 'phosphorylation'

                                mot = intera.Motif(
                                    s[0],
                                    p['start'],
                                    p['end'],
                                    instance=p['instance'],
                                    isoform=s[1])
                                ptm = intera.Ptm(s[0],
                                                 motif=mot,
                                                 residue=res,
                                                 typ=p['typ'],
                                                 source=[source],
                                                 isoform=s[1])
                                dom = intera.Domain(protein=k)

                                if 'references' not in p:
                                    p['references'] = []
                                dommot = intera.DomainMotif(
                                    domain=dom,
                                    ptm=ptm,
                                    sources=[source],
                                    refs=p['references'])

                                if source == 'MIMP':
                                    dommot.mimp_sources = ';'.split(p[
                                        'databases'])
                                    dommot.npmid = p['npmid']

                                elif source == 'PhosphoNetworks':
                                    dommot.pnetw_score = p['score']

                                elif source == 'dbPTM':
                                    dommot.dbptm_sources = [p['source']]

                                if isinstance(e, int) and e > 0:

                                    if self.graph.es[e]['ptm'] is None:
                                        self.graph.es[e]['ptm'] = []

                                    self.graph.es[e]['ptm'].append(dommot)

                                if return_raw:
                                    raw.append(dommot)

        prg.terminate()

        if trace:
            return {
                'non_matching': nomatch,
                'kinase_ambiguousity': kin_ambig,
                'substrate_ambiguousity': sub_ambig
            }

        if return_raw:
            return raw

    def load_depod_dmi(self):
        """
        """

        reres = re.compile(r'([A-Z][a-z]+)-([0-9]+)')
        non_digit = re.compile(r'[^\d.-]+')
        data = dataio.get_depod()
        aadict = dict(
            zip([a.lower().capitalize() for a in common.aaletters.keys()],
                common.aaletters.values()))

        if self.seq is None:
            self.sequences()

        if 'ptm' not in self.graph.es.attributes():
            self.graph.es['ptm'] = [[] for _ in self.graph.es]

        prg = Progress(
            len(data), 'Loading dephosphorylation data from DEPOD', 11)
        mapp = []

        for l in data:
            prg.step()
            reslist = [(aadict[r[0]], int(non_digit.sub('', r[1])))
                       for r in reres.findall(l[2]) if r[0] in aadict]

            if len(l[4]) != 0 and len(l[5]) != 0:
                enzymes = self.mapper.map_name(l[4][0], 'uniprot', 'uniprot')
                substrates = self.mapper.map_name(l[5][0], 'uniprot',
                                                  'uniprot')
                substrates = [s for s in substrates if s in self.seq]
                mapp.append([l[0], l[1], enzymes, substrates, reslist])

                for r in reslist:
                    substrates = [
                        s for s in substrates if self.seq[s].match(r[0], r[1])
                    ]

                mapp[-1].append(substrates)

                for d in enzymes:

                    for s in substrates:
                        nodes = self.get_node_pair(d, s)

                        if nodes:
                            e = self.graph.get_eid(
                                nodes[0], nodes[1], error=False)

                            if isinstance(e, int) and e > 0:

                                for r in reslist:
                                    res = intera.Residue(r[1], r[0], s)
                                    ptm = intera.Ptm(s,
                                                     residue=res,
                                                     typ='dephosphorylation',
                                                     source=['DEPOD'])
                                    dom = intera.Domain(protein=d)
                                    dommot = intera.DomainMotif(
                                        domain=dom,
                                        ptm=ptm,
                                        sources=['DEPOD'],
                                        refs=[
                                            x.strip() for x in l[3].split(',')
                                        ])

                                    if self.graph.es[e]['ptm'] is None:
                                        self.graph.es[e]['ptm'] = []

                                    self.graph.es[e]['ptm'].append(dommot)

        prg.terminate()
        # return mapp

    def uniq_ptms(self):
        """
        """

        if 'ptm' in self.graph.es.attributes():
            self.graph.es['ptm'] = [
                self.uniq_ptm(e['ptm']) for e in self.graph.es
            ]

    def uniq_ptm(self, ptms):
        """
        """

        ptms_uniq = []

        for ptm in ptms:
            merged = False

            for i, ptmu in enumerate(ptms_uniq):

                if ptm == ptmu:
                    ptms_uniq[i].merge(ptm)
                    merged = True

            if not merged:
                ptms_uniq.append(ptm)

        return ptms_uniq

    def phosphorylation_directions(self):
        """
        """

        self.uniq_ptms()
        isdir = 0

        for e in self.graph.es:

            if e['dirs'].is_directed():
                isdir += 1

        for e in self.graph.es:

            for ptm in e['ptm']:

                if ptm.__class__.__name__ == 'DomainMotif' and \
                        ptm.ptm.typ in ['phosphorylation', 'dephosphorylation']:
                    e['dirs'].set_dir((ptm.domain.protein, ptm.ptm.protein),
                                      ptm.ptm.sources)

        isdir2 = 0

        for e in self.graph.es:

            if e['dirs'].is_directed():
                isdir2 += 1

        sys.stdout.write('\t:: Directionality set for %u interactions\n'
                         '\t   based on known (de)phosphorylation events.\n' %
                         (isdir2 - isdir))
        sys.stdout.flush()

    def phosphorylation_signs(self):
        """
        """

        self.update_vname()
        new_signs = 0
        dephos = {}

        # looking up effects of dephosphorylations:
        for e in self.graph.es:

            for ptm in e['ptm']:
                sign = None

                if ptm.ptm.typ == 'dephosphorylation':
                    direction = (ptm.domain.protein, ptm.ptm.protein)
                    signs = e['dirs'].get_sign(direction)

                    if signs[0] and not signs[1]:
                        sign = 'positive'

                    elif signs[1] and not signs[0]:
                        sign = 'negative'

                    if sign:

                        if ptm.ptm.protein not in dephos:
                            dephos[ptm.ptm.protein] = []

                        dephos[ptm.ptm.protein].append({
                            'protein': ptm.ptm.protein,
                            'resaa': ptm.ptm.residue.name,
                            'resnum': ptm.ptm.residue.number,
                            'sign': sign
                        })

        # set the opposite sign for posphorylations as it is
        # at corresponding dephosphorylations:
        for e in self.graph.es:

            for ptm in e['ptm']:

                if ptm.ptm.typ == 'phosphorylation':

                    if ptm.ptm.protein in dephos:

                        for de in dephos[protein]:

                            if ptm.ptm.residue.number == de['resnum'] and \
                                    ptm.ptm.residue.name == de['resaa']:

                                direction = (ptm.domain.protein,
                                             ptm.ptm.protein)
                                signs = e['dirs'].get_sign(direction)
                                sign = 'positive' if de['sign'] == 'negative' \
                                    else 'negative'

                                if not bool(sum(signs)):
                                    e['dirs'].get_sign(direction, sign,
                                                       ptm.sources)
                                    new_signs += 1

        sys.stdout.write('\t:: Signes set based on phosphorylation-'
                         'dephosphorylation pairs: %u\n' % new_signs)
        sys.stdout.flush()

    def kinase_stats(self):
        """
        """

        counts = {}
        pcounts = {}
        ks_pairs = 0
        psite_num = 0

        for e in self.graph.es:
            ks_srcs = []

            for p in e['ptm']:

                if p.__class__.__name__ == 'DomainMotif' and p.ptm.typ == 'phosphorylation':
                    ks_srcs += p.ptm.sources + p.sources
                    p_srcs = p.ptm.sources + p.sources
                    p_srcs = sorted(set(p_srcs) - set(['Swiss-Prot']))
                    p_srcs = tuple(p_srcs)

                    if p_srcs not in pcounts:
                        pcounts[p_srcs] = 0

                    pcounts[p_srcs] += 1
                    psite_num += 1

            if len(ks_srcs) > 0:
                ks_srcs = sorted(set(ks_srcs) - set(['Swiss-Prot']))
                ks_srcs = tuple(ks_srcs)

                if ks_srcs not in counts:
                    counts[ks_srcs] = 0

                counts[ks_srcs] += 1
                ks_pairs += 1

        return {'phosphorylations': pcounts, 'kinase_substrate': counts}

    def update_db_dict(self):
        """
        """

        self.db_dict = {'nodes': {}, 'edges': {}}
        self.update_vertex_sources()
        self.update_sources()

        for e in self.graph.es:

            for s in e['sources']:

                if s not in self.db_dict['edges']:
                    self.db_dict['edges'][s] = set()

                self.db_dict['edges'][s].add(e.index)

        for v in self.graph.vs:

            for s in v['sources']:

                if s not in self.db_dict['nodes']:
                    self.db_dict['nodes'][s] = set()

                self.db_dict['nodes'][s].add(v.index)

    def sources_overlap(self, diagonal=False):
        """
        """

        self.update_db_dict()
        result = {
            'single': {
                'nodes': {},
                'edges': {}
            },
            'overlap': {
                'nodes': {},
                'edges': {}
            }
        }

        for s in self.sources:
            result['single']['nodes'][s] = len(self.db_dict['nodes'][s])
            result['single']['edges'][s] = len(self.db_dict['edges'][s])

        for s1 in self.sources:

            for s2 in self.sources:

                if diagonal or s1 != s2:
                    result['overlap']['nodes'][(s1, s2)] = \
                        len(self.db_dict['nodes'][s1] &
                            self.db_dict['nodes'][s2])
                    result['overlap']['edges'][(s1, s2)] = \
                        len(self.db_dict['edges'][s1] &
                            self.db_dict['edges'][s2])

        return result

    def source_stats(self):
        """
        """

        stats = self.sources_overlap()
        degrees = self.graph.vs.degree()
        bwness = self.graph.vs.betweenness()
        ebwness = self.graph.es.edge_betweenness()

        for s in stats['single']['nodes'].keys():
            pass

    def source_diagram(self, outf=None, **kwargs):
        """
        """

        outf = outf if outf is not None else os.path.join('pdf',
                                                          'databases.pdf')
        stats = self.sources_overlap(diagonal=True)
        sources = []
        overlaps = {}

        for s, n in iteritems(stats['single']['nodes']):
            sources.append((s, (n, None), (stats['single']['edges'][s], None)))

        for s, n in iteritems(stats['overlap']['nodes']):
            overlaps[s] = {'size': (n, stats['overlap']['edges'][s])}

        diagram = bdrawing.InterSet(sources, overlaps, outf=outf, **kwargs)
        diagram.draw()

    def get_dirs_signs(self):
        """
        """

        result = {}

        for dbs in [
                data_formats.ptm.values(),
                data_formats.interaction_htp.values(),
                data_formats.pathway.values(),
                data_formats.transcription.values()
        ]:

            for db in dbs:
                result[db.name] = [bool(db.isDirected), bool(db.sign)]

        return result

    def basic_stats(self, latex=False, caption='', latex_hdr=True, fontsize=8,
                    font='HelveticaNeueLTStd-LtCn', fname=None,
                    header_format='%s', row_order=None, by_category=True,
                    use_cats=['p', 'm', 'i', 'r'], urls=True, annots=False):
        """
        Returns basic numbers about the network resources, e.g. edge and
        node counts.

        latex
            Return table in a LaTeX document. This can be compiled by
            PDFLaTeX:
            latex stats.tex
        """

        url_strip = {
            'SPIKE': 4,
            'PDZbase': 4,
            'DEPOD': 4,
            'LMPID': 5,
            'DIP': 4,
            'MPPI': 5
        }
        outf = fname if fname is not None else os.path.join(
            'results',
            'databases-%s.tex' % self.session) if latex else os.path.join(
                'results', 'databases-%s.tsv' % self.session)
        stats = self.sources_overlap()
        dirs = self.get_dirs_signs()
        header = ['Database', 'Nodes', 'Edges', 'Directions', 'Signs']

        if urls:
            header.append('URL')

        if annots:
            header.append('Annotations')

        header.append('Notes')
        out = [header]

        desc = descriptions.descriptions

        for s in sorted(stats['single']['nodes'].keys()):
            d = desc[s] if s in desc else desc[
                'HPRD'] if s == 'HPRD-phos' else desc[
                'PhosphoSite'] if s == 'PhosphoSite_noref' else {}
            annot = ', '.join(d['annot']) if 'annot' in d else ''
            recom = d['recommend'] if 'recommend' in d else ''
            url = d['urls']['webpages'][0] \
                if 'urls' in d \
                and 'webpages' in d['urls'] \
                and len(d['urls']['webpages']) \
                else ''

            if len(url):
                _url = url.split('/')[:(3 if s not in url_strip else url_strip[
                    s])]
                url = '/'.join(_url[:3])

                if len(_url) > 3:
                    url += r'/\-'
                    url += '/'.join(_url[3:])

            url = url.replace('~', r'\textasciitilde ')
            row = [
                s, stats['single']['nodes'][s], stats['single']['edges'][s],
                int(dirs[s][0]) if s in dirs else '', int(dirs[s][1])
                if s in dirs else ''
            ]

            if urls:
                row.append(url)

            if annots:
                row.append(annot)

            row.append(recom)
            out.append(row)

        if not latex:
            out = '\n'.join(['\t'.join(x) for x in out])

            with open(outf, 'w') as f:
                f.write(out)

            sys.stdout.write('\t:: Output written to file `%s`\n' % outf)
            sys.stdout.flush()

        else:
            cats = list(
                reduce(lambda e1, e2: e1 | e2['cat'], self.graph.es, set(
                    []))) if by_category else []
            cats = set(cats) & set(use_cats)
            use_cats = list(filter(lambda c: c in cats, use_cats))
            out = dict(map(lambda l: (l[0], l[1:]), out[1:]))

            if not by_category:
                row_order = sorted(self.sources, key=lambda s: s.lower())

            else:
                row_order = []

                for cat in use_cats:
                    row_order.append((data_formats.catnames[cat], 'subtitle'))
                    row_order.extend(
                        sorted(
                            filter(lambda s: data_formats.categories[s] == cat,
                                   self.sources),
                            key=lambda s: s.lower()))

            texnewline = r'\\' + '\n'
            _latex_tab = r"""%s
                    \begin{tabularx}{0.95\textwidth}{%s}
                    \toprule
                        %s
                    \midrule
                        %s
                    \bottomrule
                    \end{tabularx}%s"""
            _latex_hdr = r"""\documentclass[a4paper,%upt]{extarticle}
                    \usepackage{fontspec}
                    \usepackage{xunicode}
                    \usepackage{polyglossia}
                    \setdefaultlanguage{english}
                    \usepackage{xltxtra}
                    \usepackage{microtype}
                    \usepackage[margin=5pt,portrait,paperwidth=23cm,paperheight=30cm]{geometry}
                    \usepackage{amsmath}
                    \usepackage{amssymb}
                    \usepackage{textcomp}
                    \usepackage[usenames,dvipsnames,svgnames,table]{xcolor}
                    \usepackage{color}
                    \usepackage{booktabs}
                    \usepackage{tabularx}
                    \setmainfont{%s}
                    \definecolor{grey875}{gray}{0.125}
                    \begin{document}
                    \color{grey875}
                    \thispagestyle{empty}
                    \vfill
                """ % (fontsize, font) if latex_hdr else ''
            _latex_end = r"""
                    \end{document}
                """ if latex_hdr else ''

            _hdr_row = ' & '.join(
                [header_format % h.replace(r'%', r'\%')
                 for h in header]) + '\\\\'
            formatter = lambda x: locale.format('%.2f', x, grouping=True) \
                if isinstance(x, float) \
                else locale.format('%d', x, grouping=True) \
                if isinstance(x, int) \
                else x
            intfloat = lambda x: float(x) if '.' in x else int(x)
            _rows = ''

            for i, k in enumerate(row_order):

                if isinstance(k, tuple) and k[1] == 'subtitle':
                    row = '\n'.join([
                        r'\midrule'
                        if i > 0 else '', r'\multicolumn{%u}{l}{%s}\\' %
                        (6 + int(annots) + int(urls), k[0]), r'\midrule', ''
                    ])

                else:
                    out[k][2] = r'$\bigstar$' if bool(out[k][2]) else ''
                    out[k][3] = r'$\bigstar$' if bool(out[k][3]) else ''
                    k0 = (
                        r'PhosphoSite\textsubscript{noref}'
                        if k == 'PhosphoSite_noref'
                        else k
                    )
                    row = ' & '.join([k0] + [
                        xx.replace('_', r'\_') if not isinstance(xx, int) else
                        '{:,}'.format(xx) for xx in out[k]
                    ]) + '\\\\\n'

                _rows += row

            _latex_tab = _latex_tab % \
                (
                    _latex_hdr,
                    'r' * (len(header) - 2) +
                    (r'p{3.0cm}<{\raggedright}' if urls else '') +
                    (r'p{2.7cm}<{\raggedright}' if annots else '') +
                    r'X<{\raggedright}',
                    _hdr_row,
                    _rows,
                    _latex_end
                )

            with open(outf, 'w') as f:
                f.write(_latex_tab)

            sys.stdout.write('\t:: Output written to file `%s`\n' % outf)
            sys.stdout.flush()

    # XXX: is this used? `PlotParam` is not defined in the whole package.
    #      Does not work.
    def source_network(self, font='HelveticaNeueLTStd'):
        """
        For EMBL branding, use Helvetica Neue Linotype Standard light
        """

        stats = self.sources_overlap()
        tupleListNodes = []
        tupleListEdges = []

        for dbs, overlap in iteritems(stats['overlap']['nodes']):

            if overlap > 0:
                tupleListNodes.append(tuple(sorted(list(dbs)) + [overlap]))

        for dbs, overlap in iteritems(stats['overlap']['edges']):

            if overlap > 0:
                tupleListEdges.append(tuple(sorted(list(dbs)) + [overlap]))

        tupleListNodes = list(set(tupleListNodes))
        tupleListEdges = list(set(tupleListEdges))
        self.sourceNetNodes = igraph.Graph.TupleList(
            tupleListNodes, weights=True)
        self.sourceNetEdges = igraph.Graph.TupleList(
            tupleListEdges, weights=True)
        # edge widths proportional to overlaps in unique nodes/edges
        self.sourceNetNodes.es['edge_width'] = [
            math.log(stats['overlap']['nodes'][(
                self.sourceNetNodes.vs[x.source]['name'],
                self.sourceNetNodes.vs[x.target]['name'])], 10)
            for x in self.sourceNetNodes.es
        ]
        self.sourceNetEdges.es['edge_width'] = [
            math.log(stats['overlap']['edges'][(
                self.sourceNetEdges.vs[x.source]['name'],
                self.sourceNetEdges.vs[x.target]['name'])], 10)
            for x in self.sourceNetEdges.es
        ]
        # node sizes proportional to number of nodes/edges
        self.sourceNetNodes.vs['vertex_size'] = [
            math.log(len(self.db_dict['nodes'][x['name']]), 10) * 5
            for x in self.sourceNetNodes.vs
        ]
        self.sourceNetEdges.vs['vertex_size'] = [
            math.log(len(self.db_dict['edges'][x['name']]), 10) * 5
            for x in self.sourceNetEdges.vs
        ]
        # numbers in vertex labels:
        v['label'] = [
            v['name'] + ' (%u)' % len(self.db_dict['nodes'][v['name']])
            for v in self.sourceNetNodes.vs
        ]
        v['label'] = [
            v['name'] + ' (%u)' % len(self.db_dict['edges'][v['name']])
            for v in self.sourceNetEdges.vs
        ]
        plotParamNodes = PlotParam( # XXX: global name PlotParam is not defined anywhere
            graph=self.sourceNetNodes,
            filename='db_net_nodes.pdf',
            layout='circular',
            vertex_label_color='#4d4d4d',
            vertex_color='#73b360',
            edge_color='#006666',
            vertex_fill_alpha='ff',
            vertex_label_font=font,
            vertex_size=self.sourceNetNodes.vs['vertex_size'],
            bbox=igraph.drawing.utils.BoundingBox(20, 20, 994, 994),
            dimensions=(1024, 1024))
        plotParamEdges = PlotParam(
            graph=self.sourceNetEdges,
            filename='db_net_nodes.pdf',
            layout='circular',
            vertex_label_color='#4d4d4d',
            vertex_color='#73b360',
            edge_color='#006666',
            vertex_fill_alpha='ff',
            vertex_label_font=font,
            vertex_size=self.sourceNetEdges.vs['vertex_size'],
            bbox=igraph.drawing.utils.BoundingBox(20, 20, 994, 994),
            dimensions=(1024, 1024))

        return plotParamNodes, plotParamEdges

    #
    # Load data from GDSC
    #
    # XXX: global name 'gdsc' is not defined (nor here or the whole package)
    def load_mutations(self, attributes=None, gdsc_datadir=None,
                       mutation_file=None):
        """
        Mutations are listed in vertex attributes. Mutation() objects
        offers methods to identify residues and look up in Ptm(), Motif()
        and Domain() objects, to check if those residues are
        modified, or are in some short motif or domain.
        """

        self.mutation_samples = []
        attributes = attributes if attributes is not None else [
            'Consequence', 'COSMIC_ID', 'SAMPLE_NAME', 'Tissue_TCGA',
            'ZYGOSITY'
        ]
        self.update_vname()
        g = gdsc.GDSC(datadir=gdsc_datadir)
        data = g.read_mutations(attributes=attributes, infile=mutation_file)

        if 'mut' not in self.graph.vs.attributes():
            self.graph.vs['mut'] = [{} for _ in self.graph.vs]

        prg = Progress(len(data), 'Processing mutations', 33)

        for uniprot, muts in iteritems(data):
            prg.step()

            if uniprot in self.nodDct:

                for mu in muts:

                    if self.graph.vs[self.nodDct[uniprot]]['mut'] is None:
                        self.graph.vs[self.nodDct[uniprot]]['mut'] = {}

                    if mu.sample not in \
                            self.graph.vs[self.nodDct[uniprot]]['mut']:
                        self.graph.vs[self.nodDct[uniprot]]['mut']\
                            [mu.sample] = []

                    self.graph.vs[self.nodDct[uniprot]]['mut']\
                        [int(mu.sample)].append(mu)

                    if mu.sample not in self.mutation_samples:
                        self.mutation_samples.append(int(mu.sample))

        prg.terminate()
        data = None

    def load_expression(self, array=False):
        """
        Expression data can be loaded into vertex attributes,
        or into a pandas DataFrame – the latter offers faster
        ways to process and use these huge matrices.
        """

        self.update_vname()
        data = gdsc.read_transcriptomics()

        if 'exp' not in self.graph.vs.attributes():
            self.graph.vs['exp'] = [{} for _ in self.graph.vs]

        prg = Progress(len(data), 'Processing transcriptomics', 33)

        if array:

            for gsymb in data.keys():
                prg.step()

                for up in self.mapper.map_name(gsymb, 'genesymbol', 'uniprot'):

                    if up in self.nodDct:
                        data[up] = data[gsymb]

                del data[gsymb]

            self.exp = pandas.DataFrame(data)

        else:

            for gsymb, expr in iteritems(data):
                prg.step()
                uniprots = self.mapper.map_name(gsymb, 'genesymbol', 'uniprot')

                for up in uniprots:

                    if up in self.nodDct:
                        self.graph.vs[self.nodDct[up]]['exp'] = \
                            dict(expr.items() +
                                 self.graph.vs[self.nodDct[up]]['exp'].items())

        prg.terminate()

    def edges_expression(self, func=lambda x, y: x * y):
        """
        Executes function `func` for each pairs of connected proteins in the
        network, for every expression dataset. By default, `func` simply
        gives the product the (normalized) expression values.

        func : callable
            Function to handle 2 vectors (pandas.Series() objects), should
            return one vector of the same length.
        """

        self.update_vindex()
        self.exp_prod = pandas.DataFrame(index=self.exp.index, columns=[])
        prg = Progress(self.graph.ecount(),
                       'Weigh edges based on protein expression', 37)

        for e in self.graph.es:
            prg.step()
            nodes = self.edge_names(e)

            if nodes[0] in self.exp.columns and nodes[1] in self.exp.columns:
                self.exp_prod[e.index] = func(self.exp[nodes[0]],
                                              self.exp[nodes[1]])

            else:
                self.exp_prod[e.index] = pandas.Series([None] *
                                                       self.exp_prod.shape[0])

        prg.terminate()

    #
    # Find and remove edges where mutations disrupt PTMs
    #

    def mutated_edges(self, sample):
        """
        Compares the mutated residues and the modified residues in PTMs.
        Interactions are marked as mutated if the target residue in the
        underlying PTM is mutated.
        """

        toDel = []
        disrupted = {}

        for e in self.graph.es:

            for v in [e.source, e.target]:

                for smpl, mut in iteritems(self.graph.vs[v]['mut']):

                    if smpl == sample:

                        for ptm in e['ptm']:

                            if mut.original == ptm.ptm.residue:
                                toDel.append(e.index)

                                if e.index not in disrupted:
                                    disrupted[e.index] = \
                                        {'ptms': len(e['ptm']), 'disr': 0}

                                disrupted[e.index]['disr'] += 1

        return toDel, disrupted

    def select_by_go(self, go_terms, go_desc=None, aspects=('C', 'F', 'P'),
                     method='ANY'):
        """
        Selects the nodes annotated by certain GO terms.

        Returns set of vertex IDs.

        :param str method:
            If `ANY` nodes annotated with any of the terms returned.
            If `ALL` nodes annotated with all the terms returned.
        """

        _method = (
            lambda s1, s2: not s2.difference(s1)
            if method == 'ALL' else
            lambda s1, s2: s1.intersection(s2)
        )

        def _method(s1, s2):
            return s1.intersection(s2)

        if go_desc is None:
            go_desc = dataio.go_descendants_goose(aspects = aspects)

        go_terms = (
            set(go_terms)
            if type(go_terms) in {set, list, tuple} else
            {go_terms}
        )
        all_desc = set.union(
            *(go_desc[term] for term in go_terms if term in go_desc)
        )

        if 'go' not in self.graph.vs.attributes():
            self.go_annotate(aspects = aspects)

        vids = set(
            i for i, v in enumerate(self.graph.vs)
            if any(
                _method(v['go'][a], all_desc)
                for a in aspects
            )
        )

        return vids

    def label_by_go(self, label, go_terms, **kwargs):
        """
        Assigns a boolean vertex attribute to nodes which tells whether
        the node is annotated by all or any (see ``method`` parameter of
        ``select_by_go``) the GO terms.
        """

        vids = self.select_by_go(go_terms, **kwargs)
        self.graph.vs[label] = [False for _ in xrange(self.graph.vcount())]

        for vid in vids:
            self.graph.vs[vid][label] = True

    def load_ligand_receptor_network(self, lig_rec_resources=True,
                                     inference_from_go=True,
                                     sources=None,
                                     keep_undirected=False,
                                     keep_rec_rec=False,
                                     keep_lig_lig=False):
        """
        Initializes a ligand-receptor network.
        """

        if sources is None:
            sources = data_formats.pathway

        CC_EXTRACELL   = 'GO:0005576'
        CC_PLASMAMEM   = 'GO:0005887'
        MF_RECBINDING  = 'GO:0005102'
        MF_RECACTIVITY = 'GO:0038023'

        if inference_from_go:
            go_desc = dataio.go_descendants_goose(aspects = ('C', 'F'))
            self.work(sources)

            if 'go' not in self.graph.vs.attributes():
                self.go_annotate()

            vids_extracell   = self.select_by_go(CC_EXTRACELL,   go_desc)
            vids_plasmamem   = self.select_by_go(CC_PLASMAMEM,   go_desc)
            vids_recbinding  = self.select_by_go(MF_RECBINDING,  go_desc)
            vids_recactivity = self.select_by_go(MF_RECACTIVITY, go_desc)

            receptors = vids_plasmamem & vids_recactivity
            ligands   = vids_extracell & vids_recbinding

            lig_with_interactions = set()
            rec_with_interactions = set()
            lig_rec_edges = set()
            lig_lig_edges = set()
            rec_rec_edges = set()

            for e in self.graph.es:
                srcs = set(
                    self.up(u).index
                    for u in e['dirs'].src(keep_undirected)
                )
                tgts = set(
                    self.up(u).index
                    for u in e['dirs'].tgt(keep_undirected)
                )

                if srcs & ligands and tgts & receptors:
                    lig_rec_edges.add(e.index)
                    lig_with_interactions.update(srcs)
                    rec_with_interactions.update(tgts)

                elif keep_lig_lig and srcs & ligands and tgts & ligands:
                    lig_lig_edges.add(e.index)

                elif keep_rec_rec and srcs & receptors and tgts & receptors:
                    rec_rec_edges.add(e.index)

            self.set_boolean_vattr('ligand_go',   lig_with_interactions)
            self.set_boolean_vattr('receptor_go', rec_with_interactions)
            edges_to_delete = (
                set(xrange(self.graph.ecount())) -
                set.union(
                    lig_rec_edges,
                    rec_rec_edges,
                    lig_lig_edges,
                )
            )

            self.graph.delete_edges(edges_to_delete)
            self.graph.delete_vertices(
                [i for i, d in enumerate(self.graph.degree()) if not d]
            )

            for e in self.graph.es:
                e['sources'].add('GO_lig_rec')

        self.update_vname()
        self.update_sources()

        if lig_rec_resources:
            datasets = copy.deepcopy(data_formats.ligand_receptor)
            datasets['cellphonedb'].inputArgs = {
                'ligand_ligand':     keep_lig_lig,
                'receptor_receptor': keep_rec_rec,
            }
            self.load_resources(datasets)

    def set_boolean_vattr(self, attr, vids, negate = False):
        """
        """

        vids = set(vids)

        if negate:
            vids = set(xrange(self.graph.vcount())) - vids

        self.graph.vs[attr] = [i in vids for i in xrange(self.graph.vcount())]

    def go_annotate(self, aspects=('C', 'F', 'P')):
        """
        Annotates protein nodes with GO terms. In the ``go`` vertex
        attribute each node is annotated by a dict of sets where keys are
        one letter codes of GO aspects and values are sets of GO accessions.
        """

        go.annotate(self.graph, aspects = aspects)

    # old name as synonym
    load_go = go_annotate

    def go_dict(self, organism=9606):
        """
        Creates a ``pypath.go.GOAnnotation`` object for one organism in the
        dict under ``go`` attribute.

        :param int organism:
            NCBI Taxonomy ID of the organism.
        """

        if not hasattr(self, 'go'):
            self.go = {}

        self.go[organism] = go.GOAnnotation(organism)

    def go_enrichment(self, proteins=None, aspect='P', alpha=0.05,
                      correction_method='hommel', all_proteins=None):
        """
        """

        if not hasattr(self, 'go') or self.ncbi_tax_id not in self.go:
            self.go_dict()

        all_proteins = (
            set(all_proteins)
                if isinstance(all_proteins, list) else
            all_proteins
                if isinstance(all_proteins, set) else
            set(self.graph.vs['name'])
        )

        annotation = dict(
            (up, g)
            for up, g in iteritems(
                getattr(
                    self.go[self.ncbi_tax_id],
                    aspect.lower()
                )
            )
            if up in all_proteins
        )

        enr = go.GOEnrichmentSet(
            aspect=aspect,
            organism=self.ncbi_tax_id,
            basic_set=annotation,
            alpha=alpha,
            correction_method=correction_method)

        if proteins is not None:
            enr.new_set(set_names=proteins)

        return enr

    def init_gsea(self, user):
        """
        """

        self.gsea = gsea.GSEA(user=user, mapper=self.mapper)
        sys.stdout.write('\n :: GSEA object initialized, use '
                         'load_genesets() to load some of the collections.\n')
        sys.stdout.write('      e.g. load_genesets([\'H\'])\n\n')
        sys.stdout.flush()
        self.gsea.show_collections()

    def add_genesets(self, genesets):
        """
        """

        for gsetid in genesets:

            if gsetid in self.gsea.collections:
                self.gsea.load_collection(gsetid)

    def geneset_enrichment(self, proteins, all_proteins=None, geneset_ids=None,
                           alpha=0.05, correction_method='hommel'):
        """
        """

        all_proteins = self.graph.vs['name'] \
            if all_proteins is None else all_proteins
        enr = gsea.GSEABinaryEnrichmentSet(
            basic_set=all_proteins,
            gsea=self.gsea,
            geneset_ids=geneset_ids,
            alpha=alpha,
            correction_method=correction_method)
        enr.new_set(proteins)

        return enr

    def update_adjlist(self, graph=None, mode='ALL'):
        """
        Creates an adjacency list in a list of sets format.
        """

        graph = graph or self.graph
        self.adjlist = [
            set(graph.neighbors(
                node, mode=mode)) for node in xrange(graph.vcount())
        ]

    def find_all_paths(self, start, end, mode='OUT', maxlen=2,
                       graph=None, silent=False):
        """
        Finds all paths up to length `maxlen` between groups of
        vertices. This function is needed only becaues igraph`s
        get_all_shortest_paths() finds only the shortest, not any
        path up to a defined length.

        @start : int or list
            Indices of the starting node(s) of the paths.
        @end : int or list
            Indices of the target node(s) of the paths.
        @mode : 'IN', 'OUT', 'ALL'
            Passed to igraph.Graph.neighbors()
        @maxlen : int
            Maximum length of paths in steps, i.e. if maxlen = 3, then
            the longest path may consist of 3 edges and 4 nodes.
        @graph : igraph.Graph object
            The graph you want to find paths in. self.graph by default.
        """

        def find_all_paths_aux(start, end, path, maxlen=None):
            path = path + [start]

            if start == end:
                return [path]

            paths = []

            if len(path) < maxlen + 1:

                for node in self.adjlist[start] - set(path):
                    paths.extend(
                        find_all_paths_aux(node, end, path, maxlen))

            return paths

        graph = graph or self.graph

        if not hasattr(self, 'adjlist'):
            self.update_adjlist(graph, mode = mode)

        all_paths = []
        start = start if isinstance(start, list) else [start]
        end = end if isinstance(end, list) else [end]

        if not silent:
            prg = Progress(
                len(start) * len(end),
                'Looking up all paths up to length %u' % maxlen, 1)

        for s in start:

            for e in end:

                if not silent:
                    prg.step()

                all_paths.extend(find_all_paths_aux(s, e, [], maxlen))

        if not silent:
            prg.terminate()

        return all_paths

    def find_all_paths2(self, graph, start, end, mode='OUT', maxlen=2,
                        psize=100):
        """
        """

        def one_step(paths, adjlist):
            # extends all paths by one step using all neighbors in adjacency
            # list
            return [
                i for ii in [[s + [a] for a in adjlist[s[-1]] if a not in s]
                             for s in paths] for i in ii
            ]

        def parts(paths, targets, adjlist, maxlen=2, psize=100, depth=1):
            complete_paths = [p for p in paths if p[-1] in targets]

            if len(paths) > 0 and len(paths[0]) <= maxlen:

                for i in xrange(0, len(paths), psize):
                    new_paths = one_step(paths[i:i + psize], adjlist)
                    complete_paths += parts(
                        new_paths,
                        targets,
                        adjlist,
                        maxlen,
                        psize,
                        depth=depth + 1)
                    sys.stdout.write("\r" + " " * 90)
                    sys.stdout.write('\r\tDepth: %u :: Paths found: %u' %
                                     (depth, len(complete_paths)))
                    sys.stdout.flush()

            return complete_paths

        graph = self.graph if graph is None else graph
        all_paths = []
        start = start if isinstance(start, list) else [start]
        end = end if isinstance(end, list) else [end]
        send = set(end)
        sys.stdout.write('\n')
        adjlist = [
            set(graph.neighbors(
                node, mode=mode)) for node in xrange(graph.vcount())
        ]
        paths = [[s] for s in start]
        all_paths = parts(paths, end, adjlist, maxlen, psize)
        sys.stdout.write('\n')

        return all_paths

    def transcription_factors(self):
        """
        """

        return common.uniqList([
            j
            for ssl in [
                i for sl in [[
                    e['dirs'].src_by_source(s)
                    for s in e['sources_by_type']['TF']
                ] for e in self.graph.es if 'TF' in e['type']] for i in sl
            ] for j in ssl
        ])

    def neighbourhood_network(self, center, second=False):
        """
        """

        center = center if isinstance(
            center, int) else self.graph.vs['name'].index(center)

        if second:
            nodes = self.second_neighbours(
                center, indices=True, with_first=True)

        else:
            nodes = self.first_neighbours(center, indices=True)

        nodes.append(center)

        return self.graph.induced_subgraph(
            nodes, implementation='create_from_scratch')

    def get_proteomicsdb(self, user, passwd, tissues=None, pickle=None):
        """
        """

        self.proteomicsdb = proteomicsdb.ProteomicsDB(user, passwd)
        self.proteomicsdb.load(pfile=pickle)
        self.proteomicsdb.get_tissues()
        self.proteomicsdb.tissues_x_proteins(tissues=tissues)
        self.exp_samples = self.proteomicsdb.tissues_loaded

    def prdb_tissue_expr(self, tissue, prdb=None, graph=None, occurrence=1,
                         group_function=lambda x: sum(x) / float(len(x)),
                         na_value=0.0):
        """
        """

        graph = self.graph if graph is None else graph
        prdb = self.proteomicsdb if prdb is None else prdb
        samples = (
            prdb.samples[tissue]
            if tissue in prdb.samples
            else [tissue] if tissue in prdb.expression
            else []
        )

        nsamples = len(samples)
        occurrence = min(nsamples, occurrence) \
            if isinstance(occurrence, int) \
            else nsamples * occurrence

        proteins_present = set([
            uniprot
            for uniprot, cnt in iteritems(Counter(itertools.chain(*[
                        prdb.expression[sample].keys()
                        for sample in samples
                    ])))
            if cnt >= occurrence
        ]) & set(graph.vs['name'])

        expressions = dict([(uniprot, group_function([
            prdb.expression[sample][uniprot] for sample in samples
            if uniprot in prdb.expression[sample]
        ])) for uniprot in proteins_present])

        graph.vs[tissue] = [
            na_value if v['name'] not in expressions else expressions[v['name']]
            for v in graph.vs
        ]

    def load_hpa(self, normal=True, pathology=True, cancer=True,
                 summarize_pathology=True, tissues=None,
                 quality=set(['Supported', 'Approved']),
                 levels={'High': 3, 'Medium': 2, 'Low': 1, 'Not detected': 0},
                 graph=None, na_value=0):
        """
        Loads Human Protein Atlas data into vertex attributes.
        """

        graph = graph or self.graph
        hpa = dataio.get_proteinatlas(
            normal = normal, pathology = pathology or cancer)

        graph.vs['hpa'] = [{} for _ in xrange(graph.vcount())]

        if normal:

            for tissue, data in iteritems(hpa['normal']):

                if tissues is not None and tissue not in tissues:
                    continue

                for v in graph.vs:
                    v['hpa'][tissue] = (
                        levels[data[v['name']][0]]
                        if v['name'] in data and data[v['name']][1] in quality
                        else na_value)

        if cancer or pathology:
            counts = set(['Not detected', 'Low', 'Medium', 'High'])

            for tissue, data in iteritems(hpa['pathology']):

                if tissues is not None and tissue not in tissues:
                    continue

                graph.vs[tissue] = [{} if summarize_pathology else na_value
                                    for v in xrange(graph.vcount())]

                for v in graph.vs:

                    if v['name'] in data:

                        if summarize_pathology:
                            patho_sum = collections.Counter(dict(
                                i for i in data[v['name']].items()
                                if i[0] in counts))

                            if len(patho_sum):
                                v[tissue] = levels[patho_sum.most_common()[0][0]]

                        else:
                            v[tissue] = data[v['name']]

    def tissue_network(self, tissue, graph=None):
        """
        Returns a network which includes the proteins expressed in
        certain tissue according to ProteomicsDB.

        Args:
        -----
        :param str tissue:
            Tissue name as used in ProteomicsDB.
        :param igraph.Graph graph:
            A graph object, by default the `graph` attribute of
            the current instance.
        """

        graph = self.graph if graph is None else graph

        if tissue not in graph.vs.attributes():
            self.prdb_tissue_expr(tissue, graph=graph)

        return graph.induced_subgraph(
            [v.index for v in graph.vs if v[tissue] > 0.0])

    def small_plot(self, graph, **kwargs):
        """
        This method is deprecated, do not use it.
        """

        arrow_size = []
        arrow_width = []
        edge_color = []
        dgraph = graph if graph.is_directed() else self.get_directed(
            graph, True, False, True)

        toDel = []

        for e in dgraph.es:

            if not e['dirs'].is_directed and e.index not in toDel:
                opp = dgraph.get_eid(e.target, e.soure, error=False)

                if opp != -1 and not dgraph.es[opp]['dirs'].is_directed():
                    toDel.append(opp)

        dgraph.delete_edges(list(set(toDel)))

        for e in dgraph.es:
            src = dgraph.vs[e.source]['name']
            tgt = dgraph.vs[e.target]['name']

            if not e['dirs'].is_directed():
                edge_color.append('#646567')
                arrow_size.append(0.0001)
                arrow_width.append(0.0001)

            else:

                if e['dirs'].is_stimulation((src, tgt)):
                    edge_color.append('#6EA945')
                    arrow_size.append(0.5)
                    arrow_width.append(0.7)

                elif e['dirs'].is_inhibition((src, tgt)):
                    edge_color.append('#DA0025')
                    arrow_size.append(0.5)
                    arrow_width.append(0.7)

                else:
                    edge_color.append('#007B7F')
                    arrow_size.append(0.5)
                    arrow_width.append(0.7)

        p = bdrawing.Plot(
            graph=dgraph,
            edge_color=edge_color,
            edge_arrow_size=arrow_size,
            edge_arrow_width=arrow_width,
            vertex_label=dgraph.vs['label'],
            vertex_color='ltp',
            vertex_alpha='FF',
            edge_alpha='FF',
            edge_label=[', '.join(l) for l in dgraph.es['loc']],
            vertex_label_family='HelveticaNeueLT Std Lt',
            edge_label_family='HelveticaNeueLT Std Lt',
            **kwargs)

        return p.draw()

    def communities(self, method, **kwargs):
        """
        """

        graph = self.graph

        if method == 'spinglass':
            giant = self.get_giant()
            graph = giant

            if 'spins' not in kwargs:
                kwargs['spins'] = 255

        elif method == 'fastgreedy':
            undgr = self.as_undirected(combine_edges='ignore')
            graph = undgr

        if hasattr(graph, '%s_community' % method):
            to_call = getattr(graph, '%s_community' % method)

            if hasattr(to_call, '__call__'):
                comm = to_call(**kwargs)

            if comm.__class__.__name__ == 'VertexDendrogram':
                comm = comm.as_clustering()

        return comm

    def set_receptors(self):
        """
        Creates a vertex attribute `rec` with value *True* if
        the protein is a receptor, otherwise *False*.
        """

        self.update_vname()
        self.graph.vs['rec'] = [False for _ in self.graph.vs]

        if 'rec' not in self.lists:
            self.receptors_list()

        for rec in self.lists['rec']:

            if rec in self.nodDct:
                self.graph.vs[self.nodDct[rec]]['rec'] = True

    def set_kinases(self):
        """
        Creates a vertex attribute `kin` with value *True* if
        the protein is a kinase, otherwise *False*.
        """

        self.update_vname()
        self.graph.vs['kin'] = [False for _ in self.graph.vs]

        if 'kin' not in self.lists:
            self.kinases_list()

        for kin in self.lists['kin']:

            if kin in self.nodDct:
                self.graph.vs[self.nodDct[kin]]['kin'] = True

    def set_signaling_proteins(self):
        """
        Creates a vertex attribute `kin` with value *True* if
        the protein is a kinase, otherwise *False*.
        """

        self.update_vname()
        self.graph.vs['sig'] = [False for _ in self.graph.vs]

        if 'kin' not in self.lists:
            self.signaling_proteins_list()

        for sig in self.lists['sig']:

            if sig in self.nodDct:
                self.graph.vs[self.nodDct[sig]]['sig'] = True

    def set_druggability(self):
        """
        Creates a vertex attribute `dgb` with value *True* if
        the protein is druggable, otherwise *False*.
        """

        self.update_vname()
        self.graph.vs['dgb'] = [False for _ in self.graph.vs]

        if 'dgb' not in self.lists:
            self.druggability_list()

        for dgb in self.lists['dgb']:

            if dgb in self.nodDct:
                self.graph.vs[self.nodDct[dgb]]['dgb'] = True

    def set_drugtargets(self, pchembl=5.0):
        """
        Creates a vertex attribute `dtg` with value *True* if
        the protein has at least one compound binding with
        affinity higher than `pchembl`, otherwise *False*.

        Args:
        -----
        :param float pchembl:
            Pchembl threshold.
        """

        self.update_vname()
        self.graph.vs['dtg'] = [False for _ in self.graph.vs]

        if 'compounds_data' not in self.graph.vs.attributes():
            self.compounds_from_chembl(assay_types=['B'], pchembl=True)

        for v in self.graph.vs:

            for d in v['compounds_data']:
                pc = [float(p) for p in d['pchembl'] if len(p) > 0]

                if len(pc) > 0:

                    if max(pc) > pchembl:
                        v['dtg'] = True

    def set_transcription_factors(self, classes=['a', 'b', 'other']):
        """
        Creates a vertex attribute `tf` with value *True* if
        the protein is a transcription factor, otherwise *False*.

        Args:
        -----
        :param list classes:
            Classes to use from TF Census. Default is `['a', 'b', 'other']`.
        """

        self.update_vname()
        self.graph.vs['tf'] = [False for _ in self.graph.vs]

        if 'tf' not in self.lists:
            self.tfs_list()

        for tf in self.lists['tf']:

            if tf in self.nodDct:
                self.graph.vs[self.nodDct[tf]]['tf'] = True

    def set_disease_genes(self, dataset='curated'):
        """
        Creates a vertex attribute named `dis` with boolean values *True*
        if the protein encoded by a disease related gene according to
        DisGeNet.

        Args:
        -----
        :param str dataset:
            Which dataset to use from DisGeNet. Default is `curated`.
        """

        self.update_vname()
        self.graph.vs['dis'] = [False for _ in self.graph.vs]
        self.disease_genes_list(dataset=dataset)

        for tf in self.lists['dis']:

            if tf in self.nodDct:
                self.graph.vs[self.nodDct[tf]]['dis'] = True

    def get_pathways(self, source):
        """
        """

        attrname = '%s_pathways' % source
        proteins_pws = None
        interactions_pws = None

        if hasattr(dataio, attrname):
            fun = getattr(dataio, attrname)
            proteins_pws, interactions_pws = fun(mapper=self.mapper)

        return proteins_pws, interactions_pws

    def pathway_members(self, pathway, source):
        """
        Returns an iterator with the members of a single pathway.
        Apart from the pathway name you need to supply its source
        database too.
        """

        attr = '%s_pathways' % source
        if attr in self.graph.vs.attribute_names():
            return _NamedVertexSeq(
                filter(lambda v: pathway in v[attr], self.graph.vs),
                self.nodNam, self.nodLab)

        else:
            return _NamedVertexSeq([], self.nodNam, self.nodLab)

    def pathway_names(self, source, graph = None):
        """
        Returns the names of all pathways having at least one member
        in the current graph.
        """

        graph = graph or self.graph

        attr = '%s_pathways' % source

        if attr in graph.vertex_attributes():
            return set.union(*graph.vs[attr])

    def load_all_pathways(self, graph=None):
        """
        """

        self.kegg_pathways(graph=graph)
        self.signor_pathways(graph=graph)
        self.pathway_attributes(graph=graph)

    def load_pathways(self, source, graph=None):
        """
        Generic method to load pathway annotations from a resource.
        We don't recommend calling this method but either specific
        methods for a single source e.g. `kegg_pathways()`
        or `sirnor_pathways()` or call `load_all_pathways()` to
        load all resources.

        Args:
        -----

        :param str source:
            Name of the source, this need to match a method in the dict
            in `get_pathways()` method and the edge and vertex attributes
            with pathway annotations will be called "<source>_pathways".
        :param igraph.Graph graph:
            A graph, by default the default the `graph` attribute of the
            current instance.
        """

        attrname = '%s_pathways' % source
        g = self.graph if graph is None else graph
        nodDct = dict(zip(g.vs['name'], xrange(g.vcount())))
        proteins_pws, interactions_pws = self.get_pathways(source)
        g.vs[attrname] = [set([]) for _ in g.vs]
        g.es[attrname] = [set([]) for _ in g.es]

        if isinstance(proteins_pws, dict):

            for pw, proteins in iteritems(proteins_pws):

                for protein in proteins:

                    if protein in proteins:
                        uniprots = self.mapper.map_name(protein, 'uniprot',
                                                        'uniprot')

                        for u in uniprots:

                            if u in nodDct:
                                g.vs[nodDct[u]][attrname].add(pw)

        if isinstance(interactions_pws, dict):

            for pw, ia in iteritems(interactions_pws):

                for pair in ia:
                    usrcs = self.mapper.map_name(pair[0], 'uniprot', 'uniprot')
                    utgts = self.mapper.map_name(pair[1], 'uniprot', 'uniprot')

                    for usrc in usrcs:

                        for utgt in utgts:

                            if usrc in nodDct and utgt in nodDct:
                                eid = g.get_eid(
                                    nodDct[usrc], nodDct[utgt], error=False)

                                if eid != -1:
                                    g.es[eid][attrname].add(pw)
        self.update_pathway_types()
        self.update_pathways()

    def signor_pathways(self, graph=None):
        """
        """

        self.load_pathways('signor', graph=graph)

    def kegg_pathways(self, graph=None):
        """
        """

        self.load_pathways('kegg', graph=graph)

    def pathway_attributes(self, graph=None):
        """
        """

        g = self.graph if graph is None else graph

        if 'netpath_pathways' in g.es.attributes():
            g.vs['netpath_pathways'] = [set([]) for _ in xrange(g.vcount())]

            for e in g.es:
                g.vs[e.source]['netpath_pathways'] = \
                    g.vs[e.source]['netpath_pathways'] | set(
                        e['netpath_pathways'])

                g.vs[e.target]['netpath_pathways'] = \
                    g.vs[e.target]['netpath_pathways'] | set(
                        e['netpath_pathways'])

        if 'slk_pathways' in g.vs.attributes():
            g.vs['signalink_pathways'] = [set(v['slk_pathways']) for v in g.vs]

            for v in g.vs:

                if v['atg']:
                    v['signalink_pathways'].add('autophagy')

            g.es['signalink_pathways'] = [
                g.vs[e.source]['signalink_pathways'] |
                g.vs[e.target]['signalink_pathways'] for e in g.es]

    def pathways_table(self, filename='genes_pathways.list',
                       pw_sources=['signalink', 'signor', 'netpath', 'kegg'],
                       graph=None):
        """
        """

        result = []
        hdr = ['UniProt', 'GeneSymbol', 'Database', 'Pathway']
        g = self.graph if graph is None else graph
        self.genesymbol_labels(graph=g)

        for v in g.vs:

            for src in pw_sources:
                pw_attr = '%s_pathways' % src

                for pw in v[pw_attr]:
                    result.append([v['name'], v['label'], src, pw])

        with open(filename, 'w') as f:
            f.write('\t'.join(hdr))
            f.write('\n'.join(map(lambda l: '\t'.join(l), result)))

    def guide2pharma(self):
        """
        """

        result = []
        data = dataio.get_guide2pharma()

        for d in data:
            ulig = []
            urec = self.mapper.map_name(d['receptor_uniprot'], 'uniprot',
                                        'uniprot')

            if len(d['ligand_uniprot']) > 0:
                ulig = self.mapper.map_name(d['ligand_uniprot'], 'uniprot',
                                            'uniprot')

            if len(d['ligand_genesymbol']) > 0:
                ulig += self.mapper.map_name(d['ligand_genesymbol'],
                                             'genesymbol', 'uniprot')

            if len(ulig) > 0 and len(urec) > 0:

                for ur in urec:

                    for ul in ulig:
                        result.append([ul, ur, d['effect'], d['pubmed'], True])

        return result

    def set_tfs(self, classes=['a', 'b', 'other']):
        """
        """

        self.set_transcription_factors(classes)

    def load_disgenet(self, dataset='curated', score=0.0, umls = False,
                      full_data = False):
        """
        Assigns DisGeNet disease-gene associations to the proteins
        in the network. Disease annotations will be added to the `dis`
        vertex attribute.

        :param float score: Confidence score from DisGeNet. Only associations
            above the score provided will be considered.
        :param bool ulms: By default we assign a list of disease names to
            each protein. To use Unified Medical Language System IDs instead
            set this to `True`.
        :param bool full_data: By default we load only disease names. Set this
            to `True` if you wish to load additional annotations like number
            of PubMed IDs, number of SNPs and original sources.
        """

        self.update_vname()
        data = dataio.get_disgenet(dataset=dataset)
        self.graph.vs['dis'] = [[] for _ in self.graph.vs]

        for d in data:

            if d['score'] >= score:
                uniprots = self.mapper.map_name(d['entrez'], 'entrez',
                                                'uniprot')

                for up in uniprots:

                    if up in self.nodInd:

                        if full_data:
                            _ = d.pop('entrez', None)
                            _ = d.pop('genesymbol', None)
                            self.graph.vs[self.nodDct[up]]['dis'].append(d)

                        elif umls:
                            self.graph.vs[self.nodDct[up]]['dis'].append(d[
                                'umls'])

                        else:
                            self.graph.vs[self.nodDct[up]]['dis'].append(d[
                                'disease'])

    def curation_stats(self, by_category=True):
        """
        """

        result = {}
        all_refs = len(
            set(
                common.flatList([[r.pmid for r in e['references']]
                                 for e in self.graph.es])))

        cats = list(
            reduce(lambda e1, e2: e1 | e2['cat'], self.graph.es, set(
                []))) if by_category else []

        for s in list(self.sources) + cats:
            sattr = 'cat' if s in data_formats.catnames else 'sources'
            rattr = 'refs_by_cat' if s in data_formats.catnames else 'refs_by_source'

            cat = None if s in data_formats.catnames \
                or s not in data_formats.categories \
                else data_formats.categories[s]

            catmembers = set(data_formats.catnames.keys()) \
                if s in data_formats.catnames \
                else set(self.sources) if not hasattr(data_formats, cat) \
                else getattr(data_formats, cat)

            src_nodes = len([v for v in self.graph.vs if s in v[sattr]])
            cat_nodes = self.graph.vcount() if cat is None else len(
                [v for v in self.graph.vs if len(v[sattr] & catmembers)])

            src_nodes_pct = src_nodes / float(cat_nodes) * 100.0
            only_src_nodes = len([
                v for v in self.graph.vs
                if s in v[sattr] and len(v[sattr] & catmembers) == 1
            ])

            only_src_nodes_pct = only_src_nodes / float(cat_nodes) * 100.0
            shared_nodes = len([
                v for v in self.graph.vs
                if s in v[sattr] and len(v[sattr] & catmembers) > 1
            ])

            src_edges = len([e for e in self.graph.es if s in e[sattr]])
            cat_edges = self.graph.ecount() if cat is None else len(
                [e for e in self.graph.es if len(e[sattr] & catmembers)])

            src_edges_pct = src_edges / \
                float(cat_edges) * 100.0

            only_src_edges = len([
                e for e in self.graph.es
                if s in e[sattr] and len(e[sattr] & catmembers) == 1
            ])

            only_src_edges_pct = only_src_edges / float(cat_edges) * 100.0
            shared_edges = len([
                e.index for e in self.graph.es
                if s in e[sattr] and len(e[sattr] & catmembers) > 1
            ])

            src_refs = set(
                common.uniqList([
                    r.pmid
                    for r in common.flatList(
                        [e[rattr][s] for e in self.graph.es if s in e[rattr]])
                ]))

            other_refs = set(
                common.flatList([[
                    r.pmid
                    for r in common.flatList([
                        rr for sr, rr in iteritems(e[rattr])
                        if sr != s and sr in catmembers
                    ])
                ] for e in self.graph.es]))

            only_src_refs = len(src_refs - other_refs)
            only_src_refs_pct = len(src_refs - other_refs) / float(
                all_refs) * 100.0

            src_refs_pct = len(src_refs) / float(all_refs) * 100.0
            shared_refs = len(src_refs & other_refs)

            shared_curation_effort = sum([
                len(x)
                for x in [
                    set(
                        common.flatList([[(r.pmid, e.index) for r in rr]
                                         for sr, rr in iteritems(e[rattr])
                                         if sr != s and sr in catmembers])) &
                    set([(rl.pmid, e.index) for rl in e[rattr][s]])
                    for e in self.graph.es if s in e[rattr]
                ]
            ])

            src_only_curation_effort = sum([
                len(x)
                for x in [
                    set([(rl.pmid, e.index) for rl in e[rattr][s]]) - set(
                        common.flatList([[(r.pmid, e.index) for r in rr]
                                         for sr, rr in iteritems(e[rattr])
                                         if sr != s if sr in catmembers]))
                    for e in self.graph.es if s in e[rattr]
                ]
            ])

            src_curation_effort = sum([
                len(x)
                for x in [
                    set([(rl.pmid, e.index) for rl in e[rattr][s]])
                    for e in self.graph.es if s in e[rattr]
                ]
            ])
            ratio = len(src_refs) / float(src_edges)

            if s in data_formats.catnames:
                s = data_formats.catnames[s]

            result[s] = {
                'source_nodes': src_nodes,
                'source_nodes_percentage': src_nodes_pct,
                'specific_nodes': only_src_nodes,
                'specific_nodes_percentage': only_src_nodes_pct,
                'shared_nodes': shared_nodes,
                'source_edges': src_edges,
                'source_edges_percentage': src_edges_pct,
                'specific_edges': only_src_edges,
                'specific_edges_percentage': only_src_edges_pct,
                'shared_edges': shared_edges,
                'specific_refs': only_src_refs,
                'specific_refs_percentage': only_src_refs_pct,
                'source_refs': len(src_refs),
                'source_refs_percentage': src_refs_pct,
                'shared_refs': shared_refs,
                'source_curation_effort': src_curation_effort,
                'source_specific_curation_effort': src_only_curation_effort,
                'shared_curation_effort': shared_curation_effort,
                'refs_edges_ratio': ratio,
                'corrected_curation_effort': src_curation_effort * ratio}

        return result

    def table_latex(self, fname, header, data, sum_row=True, row_order=None,
                    latex_hdr=True, caption='', font='HelveticaNeueLTStd-LtCn',
                    fontsize=8, sum_label='Total', sum_cols=None,
                    header_format='%s', by_category=True):
        """
        """

        non_digit = re.compile(r'[^\d.-]+')
        row_order = sorted(data.keys(), key=lambda x: x.upper()) \
            if row_order is None else row_order

        _latex_tab = r"""%s
                \begin{tabularx}{\textwidth}{%s}
                \toprule
                    %s
                \midrule
                    %s
                %s\bottomrule
                \end{tabularx}
                %s
            """

        _latex_hdr = r"""\documentclass[a4wide,%upt]{extarticle}
                \usepackage{fontspec}
                \usepackage{xunicode}
                \usepackage{polyglossia}
                \setdefaultlanguage{english}
                \usepackage{xltxtra}
                \usepackage{microtype}
                \usepackage[margin=5pt,landscape]{geometry}
                \usepackage{amsmath}
                \usepackage{amssymb}
                \usepackage[usenames,dvipsnames,svgnames,table]{xcolor}
                \usepackage{color}
                \usepackage{booktabs}
                \usepackage{tabularx}
                \setmainfont{%s}
                \definecolor{grey875}{gray}{0.125}
                \begin{document}
                \color{grey875}
                \thispagestyle{empty}
                \vfill
                \begin{table}[h]
            """ % (fontsize, font) if latex_hdr else ''

        _latex_end = r"""
                \caption{%s}
                \end{table}
                \end{document}
            """ % caption if latex_hdr else ''

        _hdr_row = ' & '.join([''] + [
            header_format % h[1].replace(r'%', r'\%') for h in header
        ]) + '\\\\'

        formatter = lambda x: locale.format('%.2f', x, grouping=True) \
            if isinstance(x, float) \
            else locale.format('%d', x, grouping=True) \
            if isinstance(x, int) \
            else x

        intfloat = lambda x: float(x) if '.' in x else int(x)
        _rows = ''

        for i, k in enumerate(row_order):

            row = ' & '.join([k[0] if isinstance(k, tuple) else k] + [
                formatter(data[k][h[0]]) for h in header
            ]) + '\\\\'

            if isinstance(k, tuple) and k[1] == 'subtitle':
                row = '\n'.join(
                    [r'\midrule' if i > 0 else '', row, r'\midrule', ''])

            else:
                row = row + '\n'

            _rows += row

        sum_cols = xrange(len(header)) if sum_cols is None \
            else sum_cols
        _sum_row = ' & '.join([sum_label] + [
            '' if i not in sum_cols else formatter(
                sum([
                    intfloat(
                        list(
                            filter(lambda x: len(x) > 0, [
                                non_digit.sub('', str(v[h[0]])), '0'
                            ]))[0]) for v in data.values()
                ])) for i, h in enumerate(header)
        ]) + '\\\\\n'
        _latex_tab = _latex_tab % (_latex_hdr, 'X' + 'r' * len(header),
                                   _hdr_row, _rows,
                                   r'\midrule' + '\n%s' % _sum_row
                                   if sum_row else '', _latex_end)

        with open(fname, 'w') as f:
            f.write(_latex_tab)

    def curation_tab(self, fname='curation_stats.tex', by_category=True,
                     use_cats=['p', 'm', 'i', 'r'], header_size='normalsize',
                     **kwargs):
        """
        """

        if by_category and 'cat' not in self.graph.es.attributes():
            self.set_categories()

        cats = list(
            reduce(lambda e1, e2: e1 | e2['cat'], self.graph.es, set(
                []))) if by_category else []

        cats = set(cats) & set(use_cats)
        use_cats = list(filter(lambda c: c in cats, use_cats))

        header = [('source_nodes', 'Proteins'),
                  ('source_nodes_percentage', r'Proteins [%]'),
                  ('shared_nodes', 'Shared proteins'),
                  ('specific_nodes', 'Specific proteins'),
                  ('specific_nodes_percentage', r'Specific proteins [%]'),
                  ('source_edges', 'Edges'),
                  ('source_edges_percentage', r'Edges [%]'),
                  ('shared_edges', 'Shared edges'),
                  ('specific_edges', 'Specific edges'),
                  ('specific_edges_percentage', r'Specific edges [%]'),
                  ('source_refs', 'References'),
                  ('source_refs_percentage', r'References [%]'),
                  ('shared_refs', 'Shared references'),
                  ('specific_refs', 'Specific references'),
                  ('specific_refs_percentage', r'Specific references [%]'),
                  ('source_curation_effort', 'Curation effort'),
                  ('shared_curation_effort', 'Shared curation effort'), (
                      'source_specific_curation_effort',
                      'Specific curation effort'), ('refs_edges_ratio',
                                                    'References-edges ratio'),
                  ('corrected_curation_effort', 'Corrected curation effort')]
        header_format = (r'\rotatebox{90}{\%s ' % header_size) + '%s}'

        cs = self.curation_stats(by_category=by_category)

        for name, data in iteritems(cs):

            if data['specific_refs'] == 0 and data['shared_refs'] == 0:
                data['source_refs'] = 'N/A'
                data['specific_refs'] = 'N/A'
                data['shared_refs'] = 'N/A'
                data['source_refs_percentage'] = 'N/A'
                data['specific_refs_percentage'] = 'N/A'
                data['source_curation_effort'] = 'N/A'
                data['shared_curation_effort'] = 'N/A'
                data['source_specific_curation_effort'] = 'N/A'
                data['refs_edges_ratio'] = 'N/A'
                data['corrected_curation_effort'] = 'N/A'

        for key in use_cats:

            if key in data_formats.catnames:
                name = data_formats.catnames[key]

                if by_category:
                    cs[(name, 'subtitle')] = cs[name]

                else:

                    if name in cs:
                        del cs[name]

        row_order = []

        if by_category:

            for cat in use_cats:
                row_order.append((data_formats.catnames[cat], 'subtitle'))
                row_order.extend(sorted(filter(lambda name:
                            (name in data_formats.categories
                             and data_formats.categories[name] == cat),
                             cs.keys())))

        self.table_latex(fname, header, cs, header_format=header_format,
                         row_order=row_order if by_category else None,
                         by_category=by_category, sum_row=False, **kwargs)

    def load_old_omnipath(self, kinase_substrate_extra = False,
                          remove_htp = False, htp_threshold = 1,
                          keep_directed = False, min_refs_undirected = 2):
        """
        Loads the OmniPath network as it was before August 2016.
        Furthermore it gives some more options.
        """

        self.load_omnipath(**locals())

    def load_omnipath(self, kinase_substrate_extra = False, remove_htp = True,
                      htp_threshold = 1, keep_directed = True,
                      min_refs_undirected = 2, old_omnipath_resources=False):
        """
        Loads the OmniPath network.
        """

        # XXX: According to the alias above omnipath = data_formats.omnipath already

        if old_omnipath_resources:
            omnipath = copy.deepcopy(data_formats.omnipath)
            omnipath['biogrid'] = data_formats.interaction['biogrid']
            omnipath['alz'] = data_formats.interaction['alz']
            omnipath['netpath'] = data_formats.interaction['netpath']
            exclude = ['intact', 'hprd']

        else:
            omnipath = data_formats.omnipath
            exclude = []

        self.load_resources(omnipath, exclude=exclude)

        if kinase_substrate_extra:
            self.load_resources(data_formats.ptm_misc)

        self.third_source_directions()

        if remove_htp:
            self.remove_htp(threshold=htp_threshold, keep_directed=keep_directed)

        if not keep_directed:
            self.remove_undirected(min_refs=min_refs_undirected)

    def remove_htp(self, threshold=50, keep_directed=False):
        """
        """

        self.htp_stats()
        vcount_before = self.graph.vcount()
        ecount_before = self.graph.ecount()
        htedgs = [
            e.index for e in self.graph.es
            if len(
                set([r.pmid for r in e['references']]) - self.htp[threshold][
                    'htrefs']) == 0 and (not keep_directed or not e['dirs']
                                         .is_directed())
        ]

        self.graph.delete_edges(htedgs)
        zerodeg = [v.index for v in self.graph.vs if v.degree() == 0]
        self.graph.delete_vertices(zerodeg)
        self.update_vname()
        sys.stdout.write(
            '\t:: Interactions with only high-throughput references '
            'have been removed.\n\t   %u interactions removed.\n\t   Number of edges '
            'decreased from %u to %u, number of vertices from %u to %u.\n' %
            (len(htedgs), ecount_before, self.graph.ecount(), vcount_before,
             self.graph.vcount()))

        sys.stdout.flush()

    def remove_undirected(self, min_refs=None):
        """
        """

        vcount_before = self.graph.vcount()
        ecount_before = self.graph.ecount()
        udedgs = [
            e.index for e in self.graph.es
            if not e['dirs'].is_directed() and (min_refs is None or len(
                set([r.pmid for r in e['references']])) < min_refs)
        ]

        self.graph.delete_edges(udedgs)
        zerodeg = [v.index for v in self.graph.vs if v.degree() == 0]
        self.graph.delete_vertices(zerodeg)
        self.update_vname()
        sys.stdout.write(
            '\t:: Undirected interactions %s '
            'have been removed.\n\t   %u interactions removed.\n\t   Number of edges '
            'decreased from %u to %u, number of vertices from %u to %u.\n' %
            ('' if min_refs is None else 'with less than %u references' %
             min_refs, len(udedgs), ecount_before, self.graph.ecount(),
             vcount_before, self.graph.vcount()))

        sys.stdout.flush()

    def numof_directed_edges(self):
        """
        """

        return len(list(filter(lambda e: e['dirs'].is_directed(),
                               self.graph.es)))

    def numof_undirected_edges(self):
        """
        """

        return len(list(filter(lambda e: not e['dirs'].is_directed(),
                               self.graph.es)))

    def htp_stats(self):
        """
        """

        htdata = {}
        refc = Counter(
            common.flatList((r.pmid for r in e['references'])
                            for e in self.graph.es))

        for htlim in reversed(xrange(1, 201)):
            htrefs = set([i[0] for i in refc.most_common() if i[1] > htlim])
            htedgs = [
                e.index for e in self.graph.es
                if len(set([r.pmid for r in e['references']]) - htrefs) == 0
            ]
            htsrcs = common.uniqList(
                common.flatList([self.graph.es[e]['sources'] for e in htedgs]))
            htdata[htlim] = {'rnum': len(htrefs), 'enum': len(htedgs),
                             'snum': len(htsrcs), 'htrefs': htrefs}

        self.htp = htdata

    def third_source_directions(self, graph=None, use_string_effects=False,
                                use_laudanna_data=False):
        """
        This method calls a series of methods to get
        additional direction & effect information
        from sources having no literature curated references,
        but giving sufficient evidence about the directionality
        for interactions already supported by literature
        evidences from other sources.
        """

        if use_string_effects:
            self.string_effects(graph = graph)

        self.kegg_directions(graph=graph)

        if use_laudanna_data:
            self.laudanna_effects(graph=graph)
            self.laudanna_directions(graph=graph)

        self.wang_effects(graph=graph)
        self.acsn_effects(graph=graph)
        self.phosphosite_directions(graph=graph)
        self.phosphopoint_directions(graph=graph)
        self.phosphonetworks_directions(graph=graph)
        self.mimp_directions(graph=graph)

    def kegg_directions(self, graph=None):
        """
        """

        keggd = dataio.get_kegg()
        self.process_directions(keggd, 'KEGG', stimulation='activation',
                                inhibition='inhibition', graph=graph)

    def phosphosite_directions(self, graph=None):
        """
        """

        psite = dataio.phosphosite_directions()
        self.process_directions(psite, 'PhosphoSite_dir', dirs_only=True,
                                id_type='uniprot', graph=graph)

    def phosphopoint_directions(self, graph=None):
        """
        """

        ppoint = dataio.phosphopoint_directions()
        self.process_directions(ppoint, 'PhosphoPoint', dirs_only=True,
                                id_type='genesymbol', graph=graph)

    def phosphonetworks_directions(self, graph=None):
        """
        """

        pnet = dataio.pnetworks_interactions()
        self.process_directions(pnet, 'PhosphoNetworks', dirs_only=True,
                                id_type='genesymbol', graph=graph)

    def mimp_directions(self, graph=None):
        """
        """

        mimp = dataio.mimp_interactions()
        self.process_directions(mimp, 'MIMP', dirs_only=True,
                                id_type='genesymbol', graph=graph)

    def laudanna_directions(self, graph=None):
        """
        """

        laud = dataio.get_laudanna_directions()
        self.process_directions(laud, 'Laudanna_sigflow', dirs_only=True,
                                id_type='genesymbol', graph=graph)

    def laudanna_effects(self, graph=None):
        """
        """

        laud = dataio.get_laudanna_effects()
        self.process_directions(laud, 'Laudanna_effects',
                                stimulation='activation',
                                inhibition='inhibition', directed='docking',
                                id_type='genesymbol', graph=graph)

    def string_effects(self, graph=None):
        """
        """

        string = dataio.get_string_effects()
        self.process_directions(string, 'STRING', stimulation='+',
                                inhibition='-', directed='*', id_type='ensp',
                                graph=graph)

    def acsn_effects(self, graph=None):
        """
        """

        acsnd = dataio.get_acsn_effects()
        self.process_directions(acsnd, 'ACSN', stimulation='+', inhibition='-',
                                directed='*', id_type='genesymbol',
                                graph=graph)

    def wang_effects(self, graph=None):
        """
        """

        wangd = dataio.get_wang_effects()
        self.process_directions(wangd, 'Wang', stimulation='+', inhibition='-',
                                directed='0', id_type='genesymbol',
                                graph=graph)

    def process_directions(self, dirs, name, directed=None, stimulation=None,
                           inhibition=None, graph=None, id_type=None,
                           dirs_only=False):
        """
        """

        g = graph if graph is not None else self.graph
        nodes = set(g.vs['name'])
        sourcedirs = 0
        newdirs = 0
        newsigns = 0

        for k in dirs:

            if dirs_only or k[2] == stimulation or k[2] == inhibition or k[
                    2] == directed:
                src = [k[0]] if id_type is None \
                    else self.mapper.map_name(k[0], id_type, 'uniprot')
                tgt = [k[1]] if id_type is None \
                    else self.mapper.map_name(k[1], id_type, 'uniprot')

                for s in src:

                    for t in tgt:

                        if s in nodes and t in nodes:
                            v1 = g.vs.find(name=s)
                            v2 = g.vs.find(name=t)
                            e = g.get_eid(v1, v2, error=False)

                            if e != -1:
                                sourcedirs += 1

                                if not g.es[e]['dirs'].get_dir((s, t)):
                                    newdirs += 1

                                if dirs_only or k[2] == directed:
                                    g.es[e]['dirs'].set_dir((s, t), name)

                                elif k[2] == stimulation:

                                    if not g.es[e]['dirs'].get_sign(
                                        (s, t), 'positive'):
                                        newsigns += 1

                                    g.es[e]['dirs'].set_sign((s, t),
                                                             'positive', name)

                                elif k[2] == inhibition:

                                    if not g.es[e]['dirs'].get_sign(
                                        (s, t), 'negative'):
                                        newsigns += 1

                                    g.es[e]['dirs'].set_sign((s, t),
                                                             'negative', name)
        sys.stdout.write(
            '\t:: Directions and signs set for %u edges based on %s,'
            ' %u new directions, %u new signs.\n' %
            (sourcedirs, name, newdirs, newsigns))

    def curation_effort(self, sum_by_source=False):
        """
        Returns the total number of reference-interactions pairs.

        @sum_by_source : bool
            If True, counts the refrence-interaction pairs by
            sources, and returns the sum of these values.
        """

        if sum_by_source:
            return sum(
                map(sum,
                    map(lambda rs: map(len, rs.values()), self.graph.es[
                        'refs_by_source'])))

        else:
            return sum(map(len, self.graph.es['references']))

    def export_dot(self, nodes=None, edges=None, directed=True,
                   labels='genesymbol', edges_filter=lambda e: True,
                   nodes_filter=lambda v: True, edge_sources=None,
                   dir_sources=None, graph=None, return_object=False,
                   save_dot=None, save_graphics=None, prog='neato',
                   format=None, hide=False, font=None, auto_edges=False,
                   hide_nodes=[], defaults={}, **kwargs):
        """
        Builds a pygraphviz.AGraph() object with filtering the edges
        and vertices along arbitrary criteria.
        Returns the Agraph object if requesred, or exports the dot
        file, or saves the graphics.

        @nodes : list
        List of vertex ids to be included.
        @edges : list
        List of edge ids to be included.
        @directed : bool
        Create a directed or undirected graph.
        @labels : str
        Name type to be used as id/label in the dot format.
        @edges_filter : function
        Function to filter edges, accepting igraph.Edge as argument.
        @nodes_filter : function
        Function to filter vertices, accepting igraph.Vertex as argument.
        @edge_sources : list
        Sources to be included.
        @dir_sources : list
        Direction and effect sources to be included.
        @graph : igraph.Graph
        The graph object to export.
        @return_object : bool
        Whether to return the pygraphviz.AGraph object.
        @save_dot : str
        Filename to export the dot file to.
        @save_graphics : str
        Filename to export the graphics, the extension defines the format.
        @prog : str
        The graphviz layout algorithm to use.
        @format : str
        The graphics format passed to pygraphviz.AGrapg().draw().
        @hide : bool
        Hide filtered edges instead of omit them.
        @hide nodes : list
        Nodes to hide. List of vertex ids.
        @auto_edges : str
        Automatic, built-in style for edges.
        'DIRECTIONS' or 'RESOURCE_CATEGORIES' are supported.
        @font : str
        Font to use for labels.
        For using more than one fonts refer to graphviz attributes with constant values
        or define callbacks or mapping dictionaries.
        @defaults : dict
        Default values for graphviz attributes, labeled with the entity, e.g.
        `{'edge_penwidth': 0.2}`.
        @**kwargs : constant, callable or dict
        Graphviz attributes, labeled by the target entity. E.g. `edge_penwidth`,
        'vertex_shape` or `graph_label`.
        If the value is constant, this value will be used.
        If the value is dict, and has `_name` as key, for every instance of the
        given entity, the value of the attribute defined by `_name` will be looked
        up in the dict, and the corresponding value will be given to the graphviz
        attribute. If the key `_name` is missing from the dict, igraph vertex and
        edge indices will be looked up among the keys.
        If the value is callable, it will be called with the current instance of
        the entity and the returned value will be used for the graphviz attribute.
        E.g. `edge_arrowhead(edge)` or `vertex_fillcolor(vertex)`
        Example:
            import pypath
            from pypath import data_formats
            net = pypath.PyPath()
            net.init_network(pfile = 'cache/default.pickle')
            #net.init_network({'arn': data_formats.omnipath['arn']})
            tgf = [v.index for v in net.graph.vs if 'TGF' in v['slk_pathways']]
            dot = net.export_dot(nodes = tgf, save_graphics = 'tgf_slk.pdf', prog = 'dot',
                main_title = 'TGF-beta pathway', return_object = True,
                label_font = 'HelveticaNeueLTStd Med Cn',
                edge_sources = ['SignaLink3'],
                dir_sources = ['SignaLink3'], hide = True)
        """

        _attrs = {}
        _custom_attrs = kwargs
        graph_attrs, vertex_attrs, edge_attrs = dataio.get_graphviz_attrs()
        _defaults = {'edge_color': {'undirected': {'unknown': '#CCCCCC'},
                                    'directed': {'stimulation': '#00AA00',
                                                 'inhibition': '#CC0000',
                                                 'unknown': '#0000CC'},
                                    'pathway': '#7AA0A177', 'ptm': '#C6909C77',
                                    'reaction': '#C5B26E77',
                                    'interaction': '#9D8BB777'},
                     'edge_arrowhead': {'undirected': {'unknown': 'none'},
                                        'directed': {'stimulation': 'normal',
                                                     'inhibition': 'tee',
                                                     'unknown': 'diamond'}},
                     'vertex_fillcolor': '#AAAAAA',
                     'vertex_fontcolor': '#000000'}

        #
        for k, v in iteritems(defaults):
            _defaults[k] = v
        #

        g = self.graph if graph is None else graph
        labels = 'name' if labels == 'uniprot' else 'label'
        edge_sources = edge_sources if edge_sources is None else set(
            edge_sources)
        dir_sources = dir_sources if dir_sources is None else set(dir_sources)
        labels = 'name' if labels == 'uniprot' else labels
        labels = 'label' if labels == 'genesymbol' else labels

        if labels == 'label' and g.vs['label'] == g.vs['name']:
            self.genesymbol_labels()

        if nodes is None and edges is not None:
            nodes = [
                n for e in g.es for n in (e.source, e.target)
                if e.index in edges
            ]

        elif nodes is not None and edges is None:
            edges = [
                e.index for e in g.es
                if e.source in nodes and e.target in nodes
            ]

        else:
            edges = xrange(g.ecount())
            nodes = xrange(g.vcount())

        dNodes = dict((v.index, v[labels]) for v in g.vs
                      if nodes is None or v.index in nodes)
        hide_nodes = set(hide_nodes)

        if font is not None:
            _custom_attrs['graph_fontname'] = font
            _custom_attrs['vertex_fontname'] = font
            _custom_attrs['edge_fontname'] = font

        # attribute callbacks
        for entity in ['graph', 'vertex', 'edge']:
            callbacks_dict = '%s_callbacks' % entity
            _attrs[callbacks_dict] = {}
            callbacks = _attrs[callbacks_dict]

            if entity == 'edge':

                if (auto_edges == 'RESOURCE_CATEGORIES' or
                        auto_edges == 'DIRECTIONS') \
                        and 'edge_color' not in _custom_attrs \
                        and 'edge_arrowhead' not in _custom_attrs:
                    callbacks['color'] = \
                        AttrHelper(auto_edges, 'edge_color', _defaults)

                    if auto_edges == 'DIRECTIONS':
                        callbacks['arrowhead'] = \
                            AttrHelper(auto_edges, 'edge_arrowhead', _defaults)

                    else:
                        callbacks['arrowhead'] = \
                            AttrHelper('none', 'edge_arrowhead', _defaults)

            for attr in locals()['%s_attrs' % entity].keys():
                callback_name = '%s_%s' % (entity, attr)

                if callback_name in _custom_attrs:
                    callback_value = _custom_attrs[callback_name]

                    if isinstance(callback_value, dict):

                        if '_name' not in callback_value:
                            callback_value['_name'] = 'index'
                    callbacks[attr] = AttrHelper(
                        value=callback_value, name=attr, defaults=_defaults)

        # graph
        dot = graphviz.AGraph(directed=directed)
        attrs = {}

        for gattr, fun in iteritems(_attrs['graph_callbacks']):
            attrs[gattr] = fun(g)

        attrs = common.cleanDict(attrs)

        for gattr, value in iteritems(attrs):
            dot.graph_attr[gattr] = value

        # vertices
        for vid, node in iteritems(dNodes):
            attrs = {}

            for vattr, fun in iteritems(_attrs['vertex_callbacks']):
                attrs[vattr] = fun(g.vs[vid])

            if vid in hide_nodes:
                attrs['style'] = 'invis'

            attrs = common.cleanDict(attrs)
            dot.add_node(node, **attrs)

        # edges
        edge_callbacks = _attrs['edge_callbacks']

        for eid in edges:
            s = g.es[eid].source
            t = g.es[eid].target
            sn = g.vs[s]['name']
            tn = g.vs[t]['name']
            sl = g.vs[s][labels]
            tl = g.vs[t][labels]
            d = g.es[eid]['dirs']
            # this edge should be visible at all?
            evis = edges_filter(g.es[eid]) and \
                nodes_filter(g.vs[s]) and nodes_filter(g.vs[t]) and \
                (edge_sources is None or
                    len(g.es[eid]['sources'] & edge_sources) > 0) and \
                g.es[eid].source not in hide_nodes and g.es[
                    eid].target not in hide_nodes

            if evis or hide:
                drawn_directed = False
                thisSign = 'unknown'
                thisDir = 'undirected'
                thisSources = set([])

                if directed:

                    if d.get_dir((sn, tn)):
                        sdir = d.get_dir((sn, tn), sources=True)
                        thisDir = (sn, tn)
                        vis = (dir_sources is None or
                               len(sdir & dir_sources) > 0) and evis

                        if vis or hide:
                            attrs = {}
                            sign = d.get_sign((sn, tn))
                            ssign = d.get_sign((sn, tn), sources=True)
                            drawn_directed = True

                            if vis and (sign[0] and dir_sources is None or
                                        dir_sources is not None and
                                        len(ssign[0] & dir_sources) > 0):
                                thisSign = 'stimulation'
                                thisSources = ssign[0] if dir_sources is None else \
                                    ssign[0] & dir_sources

                                for eattr, fun in iteritems(edge_callbacks):
                                    attrs[eattr] = fun(g.es[eid], thisDir,
                                                       thisSign, thisSources,
                                                       g.es[eid]['sources'])

                            elif vis and (sign[1] and dir_sources is None or
                                          dir_sources is not None and
                                          len(ssign[1] & dir_sources) > 0):
                                thisSign = 'inhibition'
                                thisSoures = ssign[1] if dir_sources is None else \
                                    ssign[1] & dir_sources

                                for eattr, fun in iteritems(edge_callbacks):
                                    attrs[eattr] = fun(g.es[eid], thisDir,
                                                       thisSign, thisSources,
                                                       g.es[eid]['sources'])

                            elif vis:
                                thisSign = 'unknown'
                                thisSources = sdir

                                for eattr, fun in iteritems(edge_callbacks):
                                    attrs[eattr] = fun(g.es[eid], thisDir,
                                                       thisSign, thisSources,
                                                       g.es[eid]['sources'])

                            else:
                                attrs['style'] = 'invis'
                                drawn_directed = False
                            attrs = common.cleanDict(attrs)
                            dot.add_edge(sl, tl, **attrs)

                    if d.get_dir((tn, sn)):
                        sdir = d.get_dir((tn, sn), sources=True)
                        thisDir = (tn, sn)
                        vis = (dir_sources is None or
                               len(sdir & dir_sources) > 0) and evis

                        if vis or hide:
                            attrs = {}
                            sign = d.get_sign((tn, sn))
                            ssign = d.get_sign((tn, sn), sources=True)
                            drawn_directed = True

                            if vis and (sign[0] and dir_sources is None or
                                        dir_sources is not None and
                                        len(ssign[0] & dir_sources) > 0):
                                thisSign = 'stimulation'
                                thisSources = ssign[0] if dir_sources is None else \
                                    ssign[0] & dir_sources

                                for eattr, fun in iteritems(edge_callbacks):
                                    attrs[eattr] = fun(g.es[eid], thisDir,
                                                       thisSign, thisSources,
                                                       (g.es[eid]['sources']))

                            elif vis and (sign[1] and dir_sources is None or
                                          dir_sources is not None and
                                          len(ssign[1] & dir_sources) > 0):
                                thisSign = 'inhibition'
                                thisSources = ssign[1] if dir_sources is None else \
                                    ssign[1] & dir_sources

                                for eattr, fun in iteritems(edge_callbacks):
                                    attrs[eattr] = fun(g.es[eid], thisDir,
                                                       thisSign, thisSources,
                                                       g.es[eid]['sources'])

                            elif vis:
                                thisSign = 'unknown'
                                thisSources = sdir

                                for eattr, fun in iteritems(edge_callbacks):
                                    attrs[eattr] = fun(g.es[eid], thisDir,
                                                       thisSign, thisSources,
                                                       g.es[eid]['sources'])

                            else:
                                attrs['style'] = 'invis'
                                drawn_directed = False
                            attrs = common.cleanDict(attrs)
                            dot.add_edge(tl, sl, **attrs)

                if not directed or d.get_dir('undirected'):
                    attrs = {}
                    thisDir = 'undirected'
                    thisSign = 'unknown'
                    thisSources = d.get_dir('undirected', sources=True)

                    for eattr, fun in iteritems(edge_callbacks):
                        attrs[eattr] = fun(g.es[eid], thisDir, thisSign,
                                           thisSources, g.es[eid]['sources'])

                    if (not evis and hide) or drawn_directed:
                        attrs['style'] = 'invis'

                    if dot.has_neighbor(sl, tl):

                        if dot.get_edge(sl, tl).attr['style'] is not None and \
                                'invis' in dot.get_edge(sl, tl).attr['style']:
                            dot.delete_edge(sl, tl)

                    if not dot.has_neighbor(sl, tl):
                        attrs = common.cleanDict(attrs)
                        dot.add_edge((sl, tl), **attrs)

        if type(save_dot) in set([str, unicode]):

            with open(save_dot, 'w') as f:
                f.write(dot.to_string())

        if type(save_graphics) in set([str, unicode]):
            dot.draw(save_graphics, format=format, prog=prog)

        if return_object:
            return dot

    def consistency(self):
        """
        """

        con = dict(map(lambda c: (c, dict(
            map(lambda t: (t, dict(
                ((s1, s2), dict(
                    map(lambda a: (a, set([]) if t.endswith('edges') else 0),
                        ['total', 'minor', 'major'])))
                for s1 in self.sources for s2 in self.sources)),
                ['directions', 'directions_edges', 'signs', 'signs_edges']))),
            ['consistency', 'inconsistency']))

        # inconsistency #
        prg = Progress(len(self.sources)**2, 'Counting inconsistency', 1)

        for s1 in self.sources:

            for s2 in self.sources:
                prg.step()
                s12 = set([s1, s2])

                for e in self.graph.es:
                    d = e['dirs']

                    if s1 in d.sources_straight() and \
                            s1 not in d.sources_reverse() and \
                            s2 in d.sources_reverse() and \
                            s2 not in d.sources_straight() or \
                            s2 in d.sources_straight() and \
                            s2 not in d.sources_reverse() and \
                            s1 in d.sources_reverse() and \
                            s1 not in d.sources_straight():
                        con['inconsistency']['directions'][(s1,
                                                            s2)]['total'] += 1
                        con['inconsistency']['directions_edges'][(
                            s1, s2)]['total'].add(e.index)

                        if s1 in d.sources_straight() and s1 != s2:

                            if len(d.sources_straight()) > len(
                                    d.sources_reverse()):
                                con['inconsistency']['directions'][(
                                    s1, s2)]['major'] += 1
                                con['inconsistency']['directions'][(
                                    s2, s1)]['minor'] += 1
                                con['inconsistency']['directions_edges'][(
                                    s1, s2)]['major'].add(e.index)
                                con['inconsistency']['directions_edges'][(
                                    s2, s1)]['minor'].add(e.index)

                            elif len(d.sources_straight()) < \
                                    len(d.sources_reverse()):
                                con['inconsistency']['directions'][(
                                    s1, s2)]['minor'] += 1
                                con['inconsistency']['directions'][(
                                    s2, s1)]['major'] += 1
                                con['inconsistency']['directions_edges'][(
                                    s1, s2)]['minor'].add(e.index)
                                con['inconsistency']['directions_edges'][(
                                    s2, s1)]['major'].add(e.index)

                            elif len(d.sources_straight()) == 1 and \
                                    len(d.sources_reverse()) == 1:
                                con['inconsistency']['directions'][(
                                    s1, s2)]['minor'] += 1
                                con['inconsistency']['directions'][(
                                    s2, s1)]['minor'] += 1
                                con['inconsistency']['directions_edges'][(
                                    s1, s2)]['minor'].add(e.index)
                                con['inconsistency']['directions_edges'][(
                                    s2, s1)]['minor'].add(e.index)

                    if s1 in d.positive_sources_straight() and \
                            s2 in d.negative_sources_straight() or \
                            s1 in d.negative_sources_straight() and \
                            s2 in d.positive_sources_straight():
                        con['inconsistency']['signs'][(s1, s2)]['total'] += 1
                        con['inconsistency']['signs_edges'][(
                            s1, s2)]['total'].add(e.index)

                        if s1 in d.positive_sources_straight():

                            if len(d.positive_sources_straight()) > \
                                    len(d.negative_sources_straight()):
                                con['inconsistency']['signs'][(
                                    s1, s2)]['major'] += 1
                                con['inconsistency']['signs'][(
                                    s2, s1)]['minor'] += 1
                                con['inconsistency']['signs_edges'][(
                                    s1, s2)]['major'].add(e.index)
                                con['inconsistency']['signs_edges'][(
                                    s2, s1)]['minor'].add(e.index)

                            elif len(d.positive_sources_straight()) < \
                                    len(d.negative_sources_straight()):
                                con['inconsistency']['signs'][(
                                    s1, s2)]['minor'] += 1
                                con['inconsistency']['signs'][(
                                    s2, s1)]['major'] += 1
                                con['inconsistency']['signs_edges'][(
                                    s1, s2)]['minor'].add(e.index)
                                con['inconsistency']['signs_edges'][(
                                    s2, s1)]['major'].add(e.index)

                            elif len(d.positive_sources_straight()) == 1 and \
                                    len(d.negative_sources_straight()) == 1:
                                con['inconsistency']['signs'][(
                                    s1, s2)]['minor'] += 1
                                con['inconsistency']['signs'][(
                                    s2, s1)]['minor'] += 1
                                con['inconsistency']['signs_edges'][(
                                    s1, s2)]['minor'].add(e.index)
                                con['inconsistency']['signs_edges'][(
                                    s2, s1)]['minor'].add(e.index)

                    if s1 in d.positive_sources_reverse() and \
                            s2 in d.negative_sources_reverse() or \
                            s1 in d.negative_sources_reverse() and \
                            s2 in d.positive_sources_reverse():
                        con['inconsistency']['signs'][(s1, s2)]['total'] += 1
                        con['inconsistency']['signs_edges'][(
                            s1, s2)]['total'].add(e.index)

                        if s1 in d.positive_sources_reverse():

                            if len(d.positive_sources_reverse()) > \
                                    len(d.negative_sources_reverse()):
                                con['inconsistency']['signs'][(
                                    s1, s2)]['major'] += 1
                                con['inconsistency']['signs'][(
                                    s2, s1)]['minor'] += 1
                                con['inconsistency']['signs_edges'][(
                                    s1, s2)]['major'].add(e.index)
                                con['inconsistency']['signs_edges'][(
                                    s2, s1)]['minor'].add(e.index)

                            elif len(d.positive_sources_reverse()) < \
                                    len(d.negative_sources_reverse()):
                                con['inconsistency']['signs'][(
                                    s1, s2)]['minor'] += 1
                                con['inconsistency']['signs'][(
                                    s2, s1)]['major'] += 1
                                con['inconsistency']['signs_edges'][(
                                    s1, s2)]['minor'].add(e.index)
                                con['inconsistency']['signs_edges'][(
                                    s2, s1)]['major'].add(e.index)

                            elif len(d.positive_sources_reverse()) == 1 and \
                                    len(d.negative_sources_reverse()) == 1:
                                con['inconsistency']['signs'][(
                                    s1, s2)]['minor'] += 1
                                con['inconsistency']['signs'][(
                                    s2, s1)]['minor'] += 1
                                con['inconsistency']['signs_edges'][(
                                    s1, s2)]['minor'].add(e.index)
                                con['inconsistency']['signs_edges'][(
                                    s2, s1)]['minor'].add(e.index)
        prg.terminate()

        # consistency #
        prg = Progress(len(self.sources)**2, 'Counting consistency', 1)

        for s1 in self.sources:

            for s2 in self.sources:
                prg.step()
                s12 = set([s1, s2])

                for e in self.graph.es:
                    d = e['dirs']

                    if s12 <= d.sources_straight():
                        con['consistency']['directions'][(s1,
                                                          s2)]['total'] += 1
                        con['consistency']['directions_edges'][(
                            s1, s2)]['total'].add(e.index)

                        if len(d.sources_straight()) > len(d.sources_reverse(
                        )) and len(s12 & d.sources_reverse()) == 0:
                            con['consistency']['directions'][(
                                s1, s2)]['major'] += 1
                            con['consistency']['directions_edges'][(
                                s1, s2)]['major'].add(e.index)

                        elif len(d.sources_straight()) < len(d.sources_reverse()) and \
                                len(s12 & d.sources_reverse()) == 0:
                            con['consistency']['directions'][(
                                s1, s2)]['minor'] += 1
                            con['consistency']['directions_edges'][(
                                s1, s2)]['minor'].add(e.index)

                    if s12 <= d.sources_reverse():
                        con['consistency']['directions'][(s1,
                                                          s2)]['total'] += 1
                        con['consistency']['directions_edges'][(
                            s1, s2)]['total'].add(e.index)

                        if len(d.sources_reverse()) > len(d.sources_straight(
                        )) and len(s12 & d.sources_straight()) == 0:
                            con['consistency']['directions'][(
                                s1, s2)]['major'] += 1
                            con['consistency']['directions_edges'][(
                                s1, s2)]['major'].add(e.index)

                        elif len(d.sources_reverse()) < len(d.sources_straight()) and \
                                len(s12 & d.sources_straight()) == 0:
                            con['consistency']['directions'][(
                                s1, s2)]['minor'] += 1
                            con['consistency']['directions_edges'][(
                                s1, s2)]['minor'].add(e.index)

                    if s12 <= d.positive_sources_straight():
                        con['consistency']['signs'][(s1, s2)]['total'] += 1
                        con['consistency']['signs_edges'][(
                            s1, s2)]['total'].add(e.index)

                        if len(d.positive_sources_straight()) > len(
                                d.positive_sources_reverse()) and len(
                                    s12 & d.positive_sources_reverse()) == 0:
                            con['consistency']['signs'][(s1, s2)]['major'] += 1
                            con['consistency']['signs_edges'][(
                                s1, s2)]['major'].add(e.index)

                        elif len(d.positive_sources_straight()) < len(d.positive_sources_reverse()) and \
                                len(s12 & d.positive_sources_reverse()) == 0:
                            con['consistency']['signs'][(s1, s2)]['minor'] += 1
                            con['consistency']['signs_edges'][(
                                s1, s2)]['minor'].add(e.index)

                    if s12 <= d.negative_sources_straight():
                        con['consistency']['signs'][(s1, s2)]['total'] += 1
                        con['consistency']['signs_edges'][(
                            s1, s2)]['total'].add(e.index)

                        if len(d.negative_sources_straight()) > len(
                                d.negative_sources_reverse()) and len(
                                    s12 & d.negative_sources_reverse()) == 0:
                            con['consistency']['signs'][(s1, s2)]['major'] += 1
                            con['consistency']['signs_edges'][(
                                s1, s2)]['major'].add(e.index)

                        elif len(d.negative_sources_straight()) < len(d.negative_sources_reverse()) and \
                                len(s12 & d.negative_sources_reverse()) == 0:
                            con['consistency']['signs'][(s1, s2)]['minor'] += 1
                            con['consistency']['signs_edges'][(
                                s1, s2)]['minor'].add(e.index)

                    if s12 <= d.positive_sources_reverse():
                        con['consistency']['signs'][(s1, s2)]['total'] += 1
                        con['consistency']['signs_edges'][(
                            s1, s2)]['total'].add(e.index)

                        if len(d.positive_sources_reverse()) > len(
                                d.positive_sources_straight()) and len(
                                    s12 & d.positive_sources_straight()) == 0:
                            con['consistency']['signs'][(s1, s2)]['major'] += 1
                            con['consistency']['signs_edges'][(
                                s1, s2)]['major'].add(e.index)

                        elif len(d.positive_sources_reverse()) < len(d.positive_sources_straight()) and \
                                len(s12 & d.positive_sources_straight()) == 0:
                            con['consistency']['signs'][(s1, s2)]['minor'] += 1
                            con['consistency']['signs_edges'][(
                                s1, s2)]['minor'].add(e.index)

                    if s12 <= d.negative_sources_reverse():
                        con['consistency']['signs'][(s1, s2)]['total'] += 1
                        con['consistency']['signs_edges'][(
                            s1, s2)]['total'].add(e.index)

                        if len(d.negative_sources_reverse()) > len(
                                d.negative_sources_straight()) and len(
                                    s12 & d.negative_sources_straight()) == 0:
                            con['consistency']['signs'][(s1, s2)]['major'] += 1
                            con['consistency']['signs_edges'][(
                                s1, s2)]['major'].add(e.index)

                        elif len(d.negative_sources_reverse()) < len(d.negative_sources_straight()) and \
                                len(s12 & d.negative_sources_straight()) == 0:
                            con['consistency']['signs'][(s1, s2)]['minor'] += 1
                            con['consistency']['signs_edges'][(
                                s1, s2)]['minor'].add(e.index)

        prg.terminate()
        return con

    def export_edgelist(self, fname, graph=None, names=['name'],
                        edge_attributes=[], sep='\t'):
        """
        Write edge list to text file with attributes

        :param fname:
            the name of the file or a stream to read from.
        :param graph:
            the igraph object containing the network
        :param names:
            list with the vertex attribute names to be printed
            for source and target vertices
        :param edge_attributes:
            list with the edge attribute names
            to be printed
        :param sep:
            string used to separate columns
        """

        # from Luis Tobalina
        graph = self.graph if graph is None else graph
        # check that input 'names' and 'edge_attributes' exist
        names = filter(lambda name: name in graph.vs.attribute_names(), names)
        edge_attributes = filter(lambda attr:
                        attr in graph.es.attribute_names(), edge_attributes)

        # write file
        with open(fname, 'wt') as fid:

            # write header
            for iname in names:
                fid.write('%s%s' % (sep.join([
                    '{}_{}'.format(st, iname) for st in ('source', 'target')
                ]), sep))

            fid.write('%s\n' % sep.join(eattr for eattr in edge_attributes))

            # write data
            for edge in graph.es:

                for iname in names:
                    fid.write('%s%s' % (sep.join([graph.vs[v][iname] for v
                                                  in edge.tuple]), sep))

                fid.write('%s\n' % sep.join(['{}'.format(edge[eattr]) for eattr
                                             in edge_attributes]))

    def in_complex(self, csources=['corum']):
        """
        """

        self.graph.es['in_complex'] = \
            [sum([len(set(self.graph.vs[e.source]['complexes'][cs].keys()) &
                      set(self.graph.vs[e.target]['complexes'][cs].keys()))
                  for cs in csources]) > 0 for e in self.graph.es]

    #
    # Methods for translating network to other organism
    #

    def translate_refsdir(self, rd, ids):
        """
        """

        new_refsdir = {}

        for k, v in iteritems(rd):
            di = (ids[k[0]], ids[k[1]]) if type(k) is tuple else k
            new_refsdir[di] = v

            return new_refsdir

    def orthology_translation(self, target, source=None, only_swissprot=True,
                              graph=None):
        """
        Translates the current object to another organism by orthology.
        Proteins without known ortholog will be deleted.

        :param int target: NCBI Taxonomy ID of the target organism.
            E.g. 10090 for mouse.
        """

        return_graph = graph is not None
        graph = self.graph if graph is None else graph
        source = self.ncbi_tax_id if source is None else source
        orto = dataio.homologene_uniprot_dict(source, target,
                                              only_swissprot=only_swissprot,
                                              mapper=self.mapper)

        vcount_before = graph.vcount()
        ecount_before = graph.ecount()

        # nodes could not be mapped are to be deleted
        vids = dict(map(lambda v: (v[1], v[0]), enumerate(graph.vs['name'])))

        orto = dict(filter(lambda i: i[0] in vids and len(i[1]),
                           iteritems(orto)))

        # print(list(iteritems(vdict))[:10])

        toDel = list(map(lambda v: v[1],
                         filter(lambda v: ((v[0] not in orto
                                    or not len(orto[v[0]])) and
                                # nodes of other species or compounds ignored
                                    graph.vs[v[1]]['ncbi_tax_id'] == source),
                                iteritems(vids))))

        ndel = len(toDel)
        graph.delete_vertices(toDel)

        # this for permanent identification of nodes:
        graph.vs['id_old'] = list(range(graph.vcount()))
        # a dict of these permanent ids and the orthologs:
        ovid_orto = dict(map(lambda v: (v['id_old'], orto[v['name']]),
                             graph.vs))

        # renaming vertices
        newnames = list(map(lambda v:(orto[v['name']][0]
                                      if v['ncbi_tax_id'] == source
                                # nodes of other species or compounds ignored
                                      else v['name']), graph.vs))

        graph.vs['name'] = newnames

        # the new nodes to be added because of ambiguous mapping
        toAdd = list(set(itertools.chain(*map(lambda v: v[1:], orto.values())))
                    # except those already exist:
                     - set(graph.vs['name']))

        graph += toAdd

        # this for permanent identification of nodes:
        graph.vs['id_new'] = list(range(graph.vcount()))

        # this is a dict of vertices to be multiplied:
                                # key is id_new
        vmul = dict(map(lambda v: (graph.vs.select(id_old = v[0])[0]['id_new'],
                                    # id_new of all new orthologs
                                   list(map(lambda vv:
                                        graph.vs.select(name=vv)[0]['id_new'],
                                        v[1]))),
                        iteritems(ovid_orto)))

        # compiling a dict of new edges to be added due to ambigous mapping

        # this is for unambiguously identify edges both at directed and
        # undirected graphs after reindexing at adding new edges:
        graph.es['id_old'] = list(range(graph.ecount()))
        graph.es['s_t_old'] = list(map(lambda e: (graph.vs[e.source]['id_new'],
                                                  graph.vs[e.target]['id_new']),
                                       graph.es))
        graph.es['u_old'] = list(map(lambda e: (graph.vs[e.source]['name'],
                                                graph.vs[e.target]['name']),
                                     graph.es))

                                # the parent edge original id as key
        edgesToAdd = dict(map(lambda epar: (epar[0], list(filter(lambda enew:
                                # removing the parent edge itself
                                (enew[0] != epar[1][0] or
                                 enew[1] != epar[1][1]),
                                itertools.product(vmul[epar[1][0]],
                                                  vmul[epar[1][1]])))),
                              map(lambda e: (e['id_old'], e['s_t_old']),
                                  graph.es)))

        # translating the dict values to vertex indices
        edgesToAddVids = list(set(map(lambda e: (
                            graph.vs.select(id_new = e[0])[0].index,
                            graph.vs.select(id_new = e[1])[0].index),
                        itertools.chain(*edgesToAdd.values()))))

        # creating new edges
        graph += edgesToAddVids

        # id_new > current index
        vids = dict(map(lambda v: (v['id_new'], v.index), graph.vs))
        # id_new > uniprot
        vnms = dict(map(lambda v: (v['id_new'], v['name']), graph.vs))
        prg = Progress(graph.ecount(), 'Translating network by homology', 21)

        # setting attributes on old and new edges:
        for e in graph.es:
            prg.step()
            d = e['dirs']

            # this lookup is appropriate as old node names are certainly
            # unique; for newly added edges `dirs` will be None
            if d is not None and d.nodes[0] in orto and d.nodes[1] in orto:
                # translation of direction object attached to original edges
                ids = {d.nodes[0]: orto[d.nodes[0]][0],
                       d.nodes[1]: orto[d.nodes[1]][0]}
                e['dirs'] = d.translate(ids)
                e['refs_by_dir'] = self.translate_refsdir(e['refs_by_dir'],
                                                          ids)

                # if new edges have been introduced
                # based on this specific edge
                if e['id_old'] in edgesToAdd:

                    # iterating new edges between orthologs
                    for enew in edgesToAdd[e['id_old']]:
                        vid1 = vids[enew[0]]
                        vid2 = vids[enew[1]]
                        # in case of directed graphs this will be correct:
                        es = graph.es.select(_source=vid1, _target=vid2)

                        if not len(es):
                            # at undirected graphs
                            # source/target might be opposite:
                            es = graph.es.select(_source=vid2, _target=vid1)

                        if not len(es):
                            sys.stdout.write('\t:: Could not find edge '\
                                             'between %s and %s!\n' % (
                                                graph.vs[vid1]['name'],
                                                graph.vs[vid2]['name']
                                                ))
                            continue

                        # this is a new edge between orthologs
                        eenew = es[0]

                        ids = {e['u_old'][0]: graph.vs[vid1]['name'],
                               e['u_old'][1]: graph.vs[vid2]['name']}

                        eenew['dirs'] = e['dirs'].translate(ids)
                        eenew['refs_by_dir'] = \
                            self.translate_refsdir(e['refs_by_dir'], ids)

                        # copying the remaining attributes
                        for eattr in e.attributes():

                            if eattr != 'dirs' and eattr != 'refs_by_dir':
                                eenew[eattr] = copy.deepcopy(e[eattr])

        prg.terminate()
        # id_new > current index
        vids = dict(map(lambda v: (v['id_new'], v.index), graph.vs))

        # setting attributes of vertices
        for vn0, vns in iteritems(vmul):
            # the first ortholog:
            v0 = graph.vs[vids[vn0]]
            # now setting its taxon to the target:
            v0['ncbi_tax_id'] = target

            for vn in vns:
                # iterating further orthologs:
                v = graph.vs[vids[vn]]

                # copying attributes:
                for vattr in v0.attributes():

                    if vattr != 'name':
                        v[vattr] = copy.deepcopy(v0[vattr])

        # removing temporary edge attributes
        del self.graph.es['id_old']
        del self.graph.vs['id_old']
        del self.graph.vs['id_new']
        del self.graph.es['s_t_old']
        del self.graph.es['u_old']

        self.collapse_by_name(graph=graph)

        if not return_graph:
            self.ncbi_tax_id = target

            self.update_vname()
            self.update_vindex()
            self.genesymbol_labels(remap_all=True)
            refl = ('uniprot', 'protein', target)

            if refl not in self.reflists:
                self.load_reflists([
                    reflists.ReferenceList(*refl, inFile='all_uniprots')])

            sys.stdout.write('\n')
            self.clean_graph()

        sys.stdout.write(' > Network successfully translated from `%u` to'\
            ' `%u`.\n   Nodes before: %u, after: %u\n   Edges before: %u,'\
            ' after %u\n' % (source, target, vcount_before, graph.vcount(),
                            ecount_before, graph.ecount()))

        if return_graph:
            return graph

    homology_translation = orthology_translation

    def random_walk_with_return(self, q, graph=None, c=.5, niter=1000):
        """
        Random walk with return (RWR) starting from one or more query nodes.
        Returns affinity (probability) vector of all nodes in the graph.

        Args:
        -----
            :param int,list q:
                Vertex IDs of query nodes.
            :param igraph.Graph graph:
                An `igraph.Graph` object.
            :param float c:
                Probability of restart.
            :param int niter:
                Number of iterations.

        Example:
        --------
            >>> import igraph
            >>> import pypath
            >>> pa = pypath.PyPath()
            >>> pa.init_network({
                    'signor': pypath.data_formats.pathway['signor']
                })
            >>> q = [
                    pa.gs('EGFR').index,
                    pa.gs('ATG4B').index
                ]
            >>> rwr = pa.random_walk_with_return(q = q)
            >>> palette = igraph.RainbowPalette(n = 100)
            >>> colors  = [palette.get(int(round(i))) for i in rwr / max(rwr) * 99]
            >>> igraph.plot(pa.graph, vertex_color = colors)
        """

        graph = graph or self._get_directed()

        if not graph.is_directed():
            sys.stdout.write('\t:: Warning: undirected graph provided\n')

        # making q a set of vertex IDs
        q = (q if type(q) is set else set(q) if type(q) is list else {q}
             if type(q) is int else None)

        if not q:
            sys.stdout.write('\t:: Warning: no starting node(s)\n')
            return np.array([0.] * graph.vcount())

        # probability per start node
        cp = c / len(q)
        # vector with restarting at starting nodes and 0 at all other nodes
        _q = np.array([cp if v in q else .0 for v in xrange(graph.vcount())])

        # vector of probabilities;
        # this will be subject of iteration
        # and its final state will be the result
        _p = copy.copy(_q)

        # transition matrix
        __A = np.array(list(graph.get_adjacency()), dtype=np.float64).T
        __A = np.nan_to_num(__A / __A.sum(0), copy=False)

        #return __A, _p, _q

        # iteration converges to affinity vector
        for _ in xrange(niter):
            _p = (1. - c) * __A.dot(_p) + _q

        return _p

    def random_walk_with_return2(self, q, c=.5, niter=1000):
        """
        Literally does random walks.
        Only for testing of the other method, to be deleted later.
        """

        # making q a set of vertex IDs
        q = (q if type(q) is set else set(q) if type(q) is list else {q}
             if type(q) is int else None)

        graph = self._get_directed()

        p = np.array([0] * graph.vcount())

        for qq in q:
            current = qq

            for _ in xrange(niter):
                p[current] += 1

                if random.random() <= c:
                    current = qq

                else:
                    successors = graph.successors(current)

                    if not successors:
                        current = qq

                    else:
                        current = random.choice(successors)

        return p

    def reload(self):
        """
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def _disclaimer(self):
        """
        """

        sys.stdout.write(self.disclaimer)

    def licence(self):
        """
        """

        self._disclaimer()
