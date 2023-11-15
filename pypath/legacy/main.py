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

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import pypath.share.session as session_mod

_logger = session_mod.Logger(name = 'main')

# external modules:
import os
import sys
import re
import importlib as imp
import math
import codecs
import random
import textwrap
import copy as copy_mod
import json
import operator
import locale
import heapq
import threading
import traceback
import itertools
import collections

from scipy import stats
import numpy as np

# Py 2/3
try:
    basestring

except NameError:
    basestring = str

# Py 2/3
try:
    input = raw_input
except NameError:
    pass

try:
    import cPickle as pickle

except ModuleNotFoundError:
    import pickle

try:
    import pygraphviz as graphviz

except ModuleNotFoundError:
    _logger._log(
        'Module `pygraphviz` not available. '
        'You don\'t need it unless you want to export dot files.'
    )

try:
    import pandas

except ModuleNotFoundError:
    _logger._log('Module `pandas` not available.')

try:
    import igraph

except ModuleNotFoundError:
    _logger._log(
        'Module `igraph` not available. '
        'The legacy `pypath.main.PyPath` class requires `igraph`.'
    )

# importing cairo here only to test if the module is available
try:
    import cairo

except ImportError: # XXX: Catching any exception like this is bad practice
                    # YYY: Why? We only want to send a log message so one can
                    # find out if this causes any trouble later.
    _logger._log(
        'Module `cairo` not available. '
        'Some plotting functionalities won\'t work.'
    )

# from this module:
import pypath
import pypath.share.cache as cache_mod
import pypath.resources.data_formats as data_formats
import pypath.utils.mapping as mapping
import pypath.resources.descriptions as descriptions
import pypath.inputs.pubmed as pubmed_input
import pypath.inputs.kegg as kegg_input
import pypath.inputs.phosphosite as phosphosite_input
import pypath.inputs.phosphonetworks as phosphonetworks_input
import pypath.inputs.mimp as mimp_input
import pypath.inputs.laudanna as laudanna_input
import pypath.inputs.string as string_input
import pypath.inputs.go as go_input
import pypath.inputs.pdb as pdb_input
import pypath.inputs.corum as corum
import pypath.inputs.havugimana as havugimana
import pypath.inputs.compleat as compleat
import pypath.inputs.complexportal as complexportal
import pypath.inputs.proteinatlas as proteinatlas
import pypath.inputs.surfaceome as surfaceome_input
import pypath.inputs.matrisome as matrisome_input
import pypath.inputs.phosphopoint as phosphopoint_input
import pypath.inputs.graphviz as graphviz_input
import pypath.inputs.disgenet as disgenet_input
import pypath.inputs.membranome as membranome_input
import pypath.inputs.wang as wang_input
import pypath.inputs.acsn as acsn_input
import pypath.inputs.exocarta as exocarta_input
import pypath.inputs.guide2pharma as g2p_input
import pypath.inputs.comppi as comppi_input
import pypath.inputs.pepcyber as pepcyber_input
import pypath.inputs.domino as domino_input
import pypath.inputs.threedid as threedid
import pypath.inputs.threedcomplex as threedcomplex
import pypath.inputs.elm as elm_input
import pypath.inputs.ielm as ielm_input
import pypath.inputs.pisa as pisa
import pypath.inputs as inputs
import pypath.core.network as network
import pypath.utils.orthology as orthology
import pypath.inputs.uniprot as uniprot_input
import pypath.inputs.pfam as pfam_input
import pypath.share.curl as curl
import pypath.internals.intera as intera
import pypath.utils.seq as se
import pypath.utils.go as go
import pypath.visual.drawing as bdrawing
import pypath.utils.proteomicsdb as proteomicsdb
import pypath.utils.reflists as reflists
import pypath.internals.input_formats as input_formats
import pypath.internals.refs as _refs
import pypath.visual.plot as plot
import pypath.core.enz_sub
import pypath.omnipath.export as export
import pypath.share.common as common
import pypath_common._constants as _const
from pypath.share.progress import *
import pypath.share.settings as settings
import pypath.core.entity as entity_mod
import pypath.utils.taxonomy as taxonomy
import pypath.legacy.db_categories as db_categories
import pypath.resources.network as network_resources
import pypath.core.evidence as evidence
import pypath.core.interaction as interaction

# to make it accessible directly from the module
omnipath = network_resources.omnipath

# XXX: The following aliases are already defined in common.py
# For compatibility with python 2, see https://docs.python.org/3/whatsnew/3.0.html
if 'long' not in dir(__builtins__):
    long = int

if 'unicode' not in dir(__builtins__):
    unicode = str

__all__ = ['PyPath', 'Direction', '_AttrHelper', 'omnipath']


NetworkEntityCollection = collections.namedtuple(
    'NetworkEntityCollection',
    [
        'total',
        'by_resource',
        'by_category',
        'shared',
        'unique',
        'shared_res_cat',
        'unique_res_cat',
        'shared_cat',
        'unique_cat',
        'resource_cat',
        'cat_resource',
        'method',
        'label',
    ],
)
NetworkEntityCollection.__new__.__defaults__ = (None,) * 8


NetworkStatsRecord = collections.namedtuple(
    'NetworkStatsRecord',
    [
        'total',
        'by_resource',
        'by_category',
        'shared',
        'unique',
        'percent',
        'shared_res_cat',
        'unique_res_cat',
        'percent_res_cat',
        'shared_cat',
        'unique_cat',
        'percent_cat',
        'resource_cat',
        'cat_resource',
        'method',
        'label',
    ],
)
NetworkStatsRecord.__new__.__defaults__ = (None,) * 11


class Direction(object):
    """
    This is a *legacy* object for handling directionality information
    associated with unique pairs of interacting molecular entities.
    The py:class:`pypath.interaction.Interaction` available in the ``attrs``
    edge attribute of the legacy :py:class:`pypath.main.PyPath`` object
    provides a clearer and much more versatile interface. This object will
    be removed at some point, we don't recommend to build applications by
    using it.

    Object storing directionality information of an edge. Also includes
    information about the reverse direction, mode of regulation and
    sources of that information.

    :arg str id_a:
        Name of the source node.
    :arg str id_b:
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
        Contains the node names [str] sorted alphabetically (*id_a*,
        *id_b*).
    :var dict positive:
        Dictionary contianing the presence/absence [bool] of positive
        interactions for both :py:attr:`straight` and :py:attr:`reverse`
        directions.
    :var dict positive_sources:
        Contains the resource names [str] supporting a positive
        interaction on :py:attr:`straight` and :py:attr:`reverse`
        directions.
    :var tuple reverse:
        Contains the node names [str] in reverse order e.g. (*id_b*,
        *id_a*).
    :var dict sources:
        Contains the resource names [str] of a given edge for each
        directionality (:py:attr:`straight`, :py:attr:`reverse` and
        ``'undirected'``). Values are sets containing the names of those
        resources supporting such directionality.
    :var tuple straight:
        Contains the node names [str] in the original order e.g.
        (*id_a*, *id_b*).
    """

    __slots__ = [
        'nodes',
        'straight',
        'reverse',
        'dirs',
        'sources',
        'positive',
        'negative',
        'positive_sources',
        'negative_sources',
    ]

    def __init__(self, id_a, id_b):
        """Initializes the edge object between the given nodes."""

        self.nodes = [id_a, id_b]
        self.nodes.sort()

        self.straight = (self.nodes[0], self.nodes[1])
        self.reverse  = (self.nodes[1], self.nodes[0])

        self.dirs = {
            self.straight: False,
            self.reverse:  False,
            'undirected':  False,
        }
        self.sources = {
            self.straight: set(),
            self.reverse:  set(),
            'undirected':  set(),
        }

        self.positive = {
            self.straight: False,
            self.reverse:  False,
        }
        self.negative = {
            self.straight: False,
            self.reverse:  False,
        }

        self.positive_sources = {
            self.straight: set(),
            self.reverse:  set(),
        }
        self.negative_sources = {
            self.straight: set(),
            self.reverse:  set(),
        }


    def reload(self):
        """Reloads the object from the module level."""

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


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


    def set_direction(self, direction, source):
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
            source = common.add_to_set(set([]), source)
            self.sources[direction] = self.sources[direction] | source


    # synonym: old name
    set_dir = set_direction


    def get_direction(self, direction, sources = False):
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


    # synonym: old name
    get_dir = get_direction


    def get_directions(self, src, tgt, sources=False):
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
                return [
                    self.sources[query],
                    self.sources[(query[1], query[0])],
                    self.sources['undirected'],
                ]

            else:
                return [
                    self.dirs[query],
                    self.dirs[(query[1], query[0])],
                    self.dirs['undirected'],
                ]

        else:
            return None


    # synonym: old name
    get_dirs = get_directions


    def which_directions(self, resources = None, effect = None):
        """
        Returns the pair(s) of nodes for which there is information
        about their directionality.

        :param str effect:
            Either *positive* or *negative*.
        :param str,set resources:
            Limits the query to one or more resources. Optional.

        :return:
            (*tuple*) -- Tuple of tuples with pairs of nodes where the
            first element is the source and the second is the target
            entity, according to the given resources and limited to the
            effect.
        """

        resources = self._resources_set(resources)
        effect = self._effect_synonyms(effect)

        return tuple(
            _dir
            for _dir, _resources in iteritems(self.sources)
            if _dir != 'undirected' and
            _resources and (
                not resources or
                resources & _resources
            ) and (
                not effect
                or (
                    not resources and
                    getattr(self, '%s_sources' % effect)
                ) or
                getattr(self, '%s_sources' % effect) & resources
            )
        )


    # synonym: old name
    which_dirs = which_directions


    def which_signs(self, resources = None, effect = None):
        """
        Returns the pair(s) of nodes for which there is information
        about their effect signs.

        :param str,set resources:
            Limits the query to one or more resources. Optional.
        :param str effect:
            Either *positive* or *negative*, limiting the query to positive
            or negative effects; for any other values effects of both
            signs will be returned.

        :return:
            (*tuple*) -- Tuple of tuples with pairs of nodes where the
            first element is a tuple of the source and the target entity,
            while the second element is the effect sign, according to
            the given resources. E.g. ((('A', 'B'), 'positive'),)
        """

        resources = self._resources_set(resources)
        effect = self._effect_synonyms(effect)
        effects = (effect,) if effect else ('positive', 'negative')

        return tuple(
            (_dir, _effect)
            for _effect in effects
            for _dir, _resources
            in iteritems(getattr(self, '%s_sources' % _effect))
            if _resources and (
                not resources or
                resources & _resources
            )
        )



    @staticmethod
    def _effect_synonyms(effect):

        if not effect:

            return

        if effect in {'positive', 'stimulation', 'stimulatory'}:

            return 'positive'

        if effect in {'negative', 'inhibition', 'inhibitory'}:

            return 'negative'


    def unset_direction(self, direction, source = None):
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

        if self.check_param(direction):

            if source is not None:

                try:
                    self.sources[direction].remove(source)

                except ValueError:
                    pass

            else:
                self.sources[direction] = []

            if len(self.sources[direction]) == 0:
                self.dirs[direction] = False


    # synonym: old name
    unset_dir = unset_direction


    def _resources_set(self, resources = None):

        return common.to_set(resources)


    def is_directed(self):
        """
        Checks if edge has any directionality information.

        :return:
            (*bool*) -- Returns ``True`` if any of the :py:attr:`dirs`
            attribute values is ``True`` (except ``'undirected'``),
            ``False`` otherwise.
        """

        return self.dirs[self.straight] or self.dirs[self.reverse]


    def is_directed_by_resources(self, resources = None):
        """
        Checks if edge has any directionality information from some
        resource(s).

        :return:
            (*bool*) -- Returns ``True`` if any of the :py:attr:`dirs`
            attribute values is ``True`` (except ``'undirected'``),
            ``False`` otherwise.
        """

        return self._by_resource(resources, op = operator.or_)


    def is_mutual(self, resources = None):
        """
        Checks if the edge has mutual directions (both A-->B and B-->A).
        """

        return (
            self.dirs[self.straight] and self.dirs[self.reverse]
                if not resources else
            self.is_mutual_by_resources(resources = resources)
        )


    def is_mutual_by_resources(self, resources = None):
        """
        Checks if the edge has mutual directions (both A-->B and B-->A)
        according to some resource(s).
        """

        return self._by_resource(resources, op = operator.and_)


    def _by_resource(self, resources = None, op = operator.or_):

        resources = self._resources_set(resources)

        return op(
            self.sources_straight() & resources,
            self.sources_reverse() & resources
        )


    def is_stimulation(self, direction = None, resources = None):
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

        return self._is_effect(
            sign = 'positive',
            direction = direction,
            resources = resources,
        )


    def is_inhibition(self, direction = None, resources = None):
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

        return self._is_effect(
            sign = 'negative',
            direction = direction,
            resources = resources,
        )


    def _is_effect(self, sign, direction = None, resources = None):

        _sign = (
            self.negative_sources
                if sign == 'negative' else
            self.positive_sources
        )
        _resources = self._resources_set(resources)

        return (
            any(
                bool(
                    val
                        if not _resources else
                    val & _resources
                )
                for _dir, val in iteritems(_sign)
                if not direction or direction == _dir
            )
        )


    def has_sign(self, direction = None, resources = None):
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

        return (
            self.is_stimulation(direction = direction, resources = resources)
                or
            self.is_inhibition(direction = direction, resources = resources)
        )


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
            source = common.add_to_set(set([]), source)

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


    def source(self, undirected = False, resources = None):
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

        return self._partner(
            source_target = 'source',
            undirected = undirected,
            resources = resources,
        )


    # synonym: old name
    src = source


    def target(self, undirected = False, resources = None):
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

        return self._partner(
            source_target = 'target',
            undirected = undirected,
            resources = resources,
        )


    # synonym: old name
    tgt = target


    def _partner(self, source_target, undirected = False, resources = None):

        resources = self._resources_set(resources)
        _slice = slice(0, 1) if source_target == 'source' else slice(1, 2)

        return tuple(itertools.chain(
            (
                _dir[_slice]
                    if _dir != 'undirected' else
                self.nodes
                    if undirected else
                ()
            )
            for _dir, _resources in iteritems(self.sources)
            if (
                (not resources and _resources) or
                (resources & _resources)
            )
        ))


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

        result = dict((d, None) for d in dirs)

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


class _AttrHelper(object):
    """
    *Legacy* object for internal use. Will be removed.
    Assists in assigning visual attributes for plotting and export methods.

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
        - *this_directed* [tuple?]: Optional, ``None`` by default.
        - *thisSign* []: Optional, ``None`` by default.
        - *this_directedSources* []: Optional, ``None`` by default.
        - *thisSources* []: Optional, ``None`` by default.
    """

    def __init__(self, value, name=None, defaults={}):
        """
        """

        self.name = name
        self.value = value
        self.defaults = defaults

        if isinstance(self.value, dict):
            self.id_type = type(self.value.keys()[0])

    def __call__(self, instance, this_directed=None, thisSign=None,
                 this_directedSources=None, thisSources=None):
        """
        """

        _this_directed = 'directed' if isinstance(this_directed, tuple) else this_directed

        # user supplied callback function:
        if hasattr(self.value, '__call__'):
            return self.value(instance)

        # special cases #1: by direction/effect
        elif (self.value == 'DIRECTIONS' and self.defaults is not None
              and self.name is not None and self.name in self.defaults):

            if _this_directed in self.defaults[self.name]:

                if thisSign in self.defaults[self.name][_this_directed]:
                    return self.defaults[self.name][_this_directed][thisSign]

        # special cases #2: by source category
        elif self.value == 'RESOURCE_CATEGORIES':

            for resource_type in ['pathway', 'ptm', 'reaction', 'interaction']:

                if len(getattr(db_categories, '%s_resources' % resource_type)
                       &thisSources) > 0:

                    if (self.name in self.defaults
                        and resource_type in self.defaults[self.name]):
                        return self.defaults[self.name][resource_type]

            sys.stdout.wrtie('No category for %s\n' % thisSources)
            sys.stdout.flush()

        # if value is constant:
        elif type(self.value) in _const.SIMPLE_TYPES:
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
    *Legacy* object, will be removed.
    A more versatile replacement in the new API is
    :py:class:`pypath.entity.EntityList`.
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


class PyPath(session_mod.Logger):
    """
    This is the a *legacy* object representing a molecular interaction
    network. At some point it will be removed, we don't recommend to rely
    on it when you build your applications.
    The :py:class:`pypath.network.Network` object offers a much clearer and
    more versatile API. As of end of 2019 not all functionalities have been
    migrated to the new API. For this reason we offer an intermediate
    solution: in this `igraph` based object the `attrs` edge attribute
    holds instances of :py:class:`pypath.interaction.Interaction` objects,
    the same type of object we use to represent interactions in the new
    :py:class:`pypath.network.Network`.
    At the same time we will keep supporting `igraph` with a method for
    converting :py:class:`pypath.network.Network` to a
    :py:class:`igraph.Graph` object, however this won't provide all the
    methods available here but will serve only the purpose to make it
    possible to use the graph theory methods from the `igraph` library
    on networks built with `pypath`.

    An object representing a molecular interaction network.

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
    :arg str name:
        Optional, ``'unnamed'`` by default. Session or project name
        (custom).
    :arg str outdir:
        Optional, ``'results'`` by default. Output directory where to
        store all output files.
    :arg int loglevel:
        Optional, 0 by default. Sets the level of the logger.
        The higher the level the more messages will be written to the log.
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
        DEPRECATED Contains the MySQL parameters used by the
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
        more details see :py:meth:`pypath.main.PyPath.edges_expression`.
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
        Contains the categories (e.g.: resource types) [str] loaded in
        the current network. Possible categories are: ``'m'`` for
        PTM/enzyme-substrate resources, ``'p'`` for pathway/activity
        flow resources, ``'i'`` for undirected/PPI resources, ``'r'``
        for process description/reaction resources and ``'t'`` for
        transcription resources.
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
        a resource. Some methods include
        :py:meth:`pypath.main.PyPath.receptor_list` (``'rec'``),
        :py:meth:`pypath.main.PyPath.druggability_list` (``'dgb'``),
        :py:meth:`pypath.main.PyPath.kinases_list` (``'kin'``),
        :py:meth:`pypath.main.PyPath.tfs_list` (``'tf'``),
        :py:meth:`pypath.main.PyPath.disease_genes_list` (``'dis'``),
        :py:meth:`pypath.main.PyPath.signaling_proteins_list`
        (``'sig'``), :py:meth:`pypath.main.PyPath.proteome_list`
        (``'proteome'``) and
        :py:meth:`pypath.main.PyPath.cancer_drivers_list` (``'cdv'``).
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
        DEPRECATED Contains the MySQL parameters used by the
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
        :py:meth:`pypath.main.PyPath.apply_negative` for more
        information.
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
        :py:meth:`pypath.main.PyPath.map_item` for more information.
    :var dict vertexAttrs:
        Stores the node attribute names [str] as keys and their
        corresponding types (e.g.: ``set``, ``list``, ``str``, ...) as
        values.
    """


    def __init__(
            self,
            ncbi_tax_id = None,
            copy = None,
            name = 'unnamed',
            cache_dir  =  None,
            outdir = 'results',
            loglevel = 0,
            loops = False,
        ):
        """Initializes the network object.

        **NOTE:** Only the instance is created, no data is donwloaded
        until the corresponding function is called (e.g.:
        :py:meth:`PyPath.init_network`).
        """

        session_mod.Logger.__init__(self, name = 'network')

        self.__version__ = pypath.__version__

        self.cache_dir = cache_mod.get_cachedir(cachedir = cache_dir)

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
            g.vs['id_type'] = []
            g.vs['original_names'] = [[] for _ in xrange(self.graph.vcount())]
            g.vs['ncbi_tax_id'] = []
            g.vs['exp'] = [{}]
            g.es['sources'] = [set([]) for _ in xrange(self.graph.ecount())]
            g.es['attrs'] = [None]
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
            self.unmapped = []
            self.name = name
            self.outdir = outdir
            self.ncbi_tax_id = ncbi_tax_id or settings.get('default_organism')
            self.default_name_type = settings.get('default_name_types')
            self.data = {}
            self.negatives = {}
            self.raw_data = None
            self.lists = {}
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

            self._log('PyPath has been initialized')

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

        self.__dict__ = copy_mod.deepcopy(other.__dict__)
        self.update_vname()
        self._log('PyPath object has been copied and reinitialized.')

    def __copy__(self):

        new = PyPath(copy = self)

        return new

    def init_network(
            self,
            lst = None,
            exclude = [],
            cache_files = {},
            pickle_file = None,
            pfile = False,
            save = False,
            reread = None,
            redownload = None,
            keep_raw = False,
            **kwargs
        ): # XXX: kwargs is not used anywhere
        """Loads the network data.

        This is a lazy way to start the module, load data and build the
        high confidence, literature curated part of the signaling
        network.

        :arg dict lst:
            Optional, ``None`` by default. Specifies the data input
            formats for the different resources (keys) [str]. Values
            are :py:class:`pypath.input_formats.NetworkInput` instances
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

        pfile = pickle_file or pfile

        if pfile and os.path.exists(pfile):

            self._log(
                'Loading igraph object from pickle `%s`...' % pfile
            )
            graph = pickle.load(open(pfile, 'rb'))

            if isinstance(graph, igraph.Graph) and graph.vcount() > 0:
                self.graph = graph
                self._log(
                    'Network loaded from `%s`. %u nodes, %u edges.' % (
                        pfile,
                        self.graph.vcount(),
                        self.graph.ecount(),
                    )
                )
                self.update_vname()
                self.update_vindex()
                self.update_sources()

                return None

        self.load_resources(
            lst = lst,
            exclude = exclude,
            reread = reread,
            redownload = redownload,
            cache_files = cache_files,
            keep_raw = keep_raw,
        )

        if save:

            self._log('Saving igraph object to file `%s`...' % pfile)
            self.save_network(pfile = pfile)
            self._log('Network saved successfully to file `%s`.' % pfile)


    def load_from_pickle(self, pickle_file):
        """
        Shortcut for loading a network from a pickle dump.
        """

        self.init_network(pfile = pickle_file)


    def save_to_pickle(self, pickle_file = None, pfile = None):
        """Saves the network object.

        Stores the instance into a pickle (binary) file which can be
        reloaded in the future.

        :arg str pickle_file:
            Optional, ``None`` by default. The path/file name where to
            store the pcikle file. If not specified, saves the network
            to its default location
            (``'cache/default_network.pickle'``).
        """

        self._log('Saving to pickle `%s`.' % pickle_file)

        pfile = (
            pickle_file or
            pfile or
            os.path.join(self.cache_dir, 'default_network.pickle')
        )
        pickle.dump(self.graph, open(pfile, 'wb'), -1)

        self._log('Saved to pickle `%s`.' % pickle_file)


    # synonym for old name
    save_network = save_to_pickle

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
            # extra_edge_attrs and extraNodeAttrs are dicts
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
                self._log(
                    'Wrong column index (%s) in extra attributes? '
                    'Line #%u' % (str(col), lnum),
                    -5,
                )
                readError = 1
                break

            fieldName = col
            attrs[fieldName] = fieldVal

        return attrs

    def get_taxon(self, tax_dict, fields): # TODO
        """
        """

        if 'A' in tax_dict and 'B' in tax_dict:

            return (
                self.get_taxon(tax_dict['A'], fields),
                self.get_taxon(tax_dict['B'], fields),
            )

        else:

            if 'dict' not in tax_dict:
                return int(fields[tax_dict['col']])

            elif fields[tax_dict['col']] in tax_dict['dict']:
                return tax_dict['dict'][fields[tax_dict['col']]]

            else:
                return None


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
        gg = g if replace else copy_mod.deepcopy(g)

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

        self._log('Updating network component lookup dictionaries.')

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


    def _filters(
            self,
            line,
            positive_filters = None,
            negative_filters = None,
        ):
        """
        Applies negative and positive filters on a line (record from an
        interaction database). If returns ``True`` the interaction will be
        discarded, if ``False`` the interaction will be further processed
        and if all other criteria fit then will be added to the network
        after identifier translation.
        """

        negative_filters = negative_filters or ()

        for filtr in negative_filters:

            if len(filtr) > 2:
                sep = filtr[2]
                thisVal = set(line[filtr[0]].split(sep))

            else:
                thisVal = set([line[filtr[0]]])

            filtrVal = common.to_set(filtr[1])

            if thisVal & filtrVal:
                return True

        positive_filters = positive_filters or ()

        for filtr in positive_filters:

            if len(filtr) > 2:
                sep = filtr[2]
                thisVal = set(line[filtr[0]].split(sep))

            else:
                thisVal = {line[filtr[0]]}

            filtrVal = common.to_set(filtr[1])

            if not thisVal & filtrVal:
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
        edge_list_mapped = []
        cache_file = cache_files[name] if name in cache_files else None

        if cache_file is not None and os.path.exists(cache_file):
            cache_type = cache_file.split('.')[-2]

            if cache_type == 'interactions':
                infile = self.read_from_cache(int_cache)

            elif cache_type == 'edges':
                edge_list_mapped = self.read_from_cache(edges_cache)

        elif os.path.exists(edges_cache):
            edge_list_mapped = self.read_from_cache(edges_cache)

        else: # XXX: You could use another elif statement here

            if os.path.exists(int_cache):
                infile = self.read_from_cache(int_cache)

        return infile, edge_list_mapped

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

        self._log(
            'Reading edge list pickle dump from cache: %s\n' % cache_file
        )

        data = pickle.load(open(cache_file, 'rb'))

        self._log('Data have been read from cache: `%s`' % cache_file)

        return data


    def _process_sign(self, signData, signDef):
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
        pos = common.to_set(signDef[1])
        neg = common.to_set(signDef[2])

        # XXX: Isn't using elif here bias the choice to stimulatory interactions
        #      even though there can also be negative sources?

        if len(signData & pos) > 0:
            stim = True

        elif len(signData & neg) > 0:
            inh = True

        return stim, inh


    def _process_direction(self, line, dir_col, dir_val, dir_sep):
        """
        Processes the direction information of an interaction according
        to a data file from a source.

        :arg list line:
            The stripped and separated line from the resource data file
            containing the information of an interaction.
        :arg int dir_col:
            The column/position number where the information about the
            direction is to be found (on *line*).
        :arg list dir_val:
            Contains the terms [str] for which that interaction is to be
            considered directed.
        :arg str dir_sep:
            Separator for the field in *line* containing the direction
            information (if any).

        :return:
            (*bool*) -- Determines whether the given interaction is
            directed or not.
        """

        if dir_col is None or dir_val is None:
            return False

        elif isinstance(line[dir_col], bool):

            return line[dir_col]

        else:

            this_directed = set(line[dir_col].split(dir_sep))
            return len(this_directed & dir_val) > 0


    def _read_network_data(
            self,
            param,
            keep_raw = False,
            cache_files = {},
            reread = None,
            redownload = None,
        ):
        """
        Reads interaction data file containing node and edge attributes
        that can be read from simple text based files and adds it to the
        networkdata. This function works not only with files, but with
        lists as well. Any other function can be written to download and
        preprocess data, and then give it to this function to finally
        attach to the network.

        :arg pypath.input_formats.NetworkInput param:
            :py:class:`pypath.input_formats.NetworkInput` instance
            containing the detailed definition of the input format of
            the file. Instead of the file name (on the
            :py:attr:`pypath.input_formats.NetworkInput.input`
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

        self._log('Reading network data from `%s`.' % param.name)

        # workaround in order to make it work with both NetworkInput
        # and NetworkResource type param
        _resource = (
            param
                if isinstance(
                    param,
                    network_resources.resource.NetworkResource
                ) else
            network_resources.resource.NetworkResource(
                name = param.name,
                interaction_type = param.interaction_type,
                networkinput = param,
                data_model = param.data_model or 'unknown',
            )
        )

        param = _resource.networkinput

        _resources_secondary = ()

        expand_complexes = (
            param.expand_complexes
                if isinstance(param.expand_complexes, bool) else
            settings.get('network_expand_complexes')
        )
        reread = (
            reread
                if isinstance(reread, bool) else
            not settings.get('network_pickle_cache')
        )

        self._log('Expanding complexes for `%s`: %s' % (
            param.name, str(expand_complexes),
        ))

        listLike = {list, tuple}
        edge_list = []
        nodeList = []
        edge_list_mapped = []
        infile = None
        _name = param.name.lower()

        int_cache = os.path.join(
            self.cache_dir,
            '%s.interactions.pickle' % _name
        )
        edges_cache = os.path.join(
            self.cache_dir,
            '%s.edges.pickle' % _name
        )

        if not reread and not redownload:

            infile, edge_list_mapped = self.lookup_cache(
                _name,
                cache_files,
                int_cache,
                edges_cache,
            )

        if not len(edge_list_mapped):

            if infile is None:

                if not isinstance(
                    param,
                    (
                        data_formats.input_formats.NetworkInput,
                        network_resources.resource.NetworkResource,
                    )
                ):

                    self._log(
                        '_read_network_data: No proper input file '
                        'definition. `param` should be either '
                        'a `pypath.input_formats.NetworkInput` or a '
                        '`pypath.resource.NetworkResource` instance.',
                        -5,
                    )

                    return None

                if param.huge:

                    sys.stdout.write(
                        '\n\tProcessing %s requires huge memory.\n'
                        '\tPlease hit `y` if you have at '
                        'least 2G free memory,\n'
                        '\tor `n` to omit %s.\n'
                        '\tAfter processing once, it will be saved in \n'
                        '\t%s, so next time can be loaded quickly.\n\n'
                        '\tProcess %s now? [y/n]\n' %
                        (param.name, param.name, edges_cache,
                         param.name))
                    sys.stdout.flush()

                    while True:
                        answer = input().lower()

                        if answer == 'n':
                            return None

                        elif answer == 'y':
                            break

                        else:
                            sys.stdout.write(
                                '\n\tPlease answer `y` or `n`:\n\t')
                            sys.stdout.flush()

                input_func = inputs.get_method(param.input)

                # reading from remote or local file, or executing import
                # function:
                if (
                    isinstance(param.input, basestring) and (
                        param.input.startswith('http') or
                        param.input.startswith('ftp')
                    )
                ):

                    curl_use_cache = not redownload
                    c = curl.Curl(
                        param.input,
                        silent=False,
                        large=True,
                        cache=curl_use_cache
                    )
                    infile = c.fileobj.read()

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
                    self._log(
                        "Retrieving data from%s ..." % param.input
                    )

                elif input_func is not None:

                    self._log("Retrieving data by `%s` ..." %
                                    input_func.__name__)

                    _store_cache = curl.CACHE

                    if isinstance(redownload, bool):

                        curl.CACHE = not redownload

                    # this try-except needs to be removed
                    # once correct exception handling will
                    # be implemented in every input function
                    try:
                        infile = input_func(**param.input_args)

                    except Exception as e:
                        self._log(
                            'Error in `%s`. '
                            'Skipping to next resource. '
                            'See below the traceback.' % input_func.__name__
                        )
                        self._log_traceback()

                        try:
                            traceback.print_tb(
                                e.__traceback__, file = sys.stdout)

                        except Exception as e:
                            self._log('Failed handling exception.')
                            self._log_traceback()

                    curl.CACHE = _store_cache

                elif os.path.isfile(param.input):

                    infile = curl.Curl(
                        param.input,
                        large = True,
                        silent = False
                    ).result

                    self._log('%s opened...' % param.input)

                if infile is None:

                    self._log(
                        '`%s`: Could not find file or input function '
                        'or failed preprocessing.' %
                        param.input,
                        -5,
                    )
                    return None

            # at this point we can be sure we have something to iterate
            # rows (either lists or strings which can be split into lists)

            # finding the largest referred column number,
            # to avoid list indices out of range
            # TODO: this should not be necessary, if an input gives an
            # IndexError that should be fixed elsewhere
            is_directed = param.is_directed
            sign = param.sign
            refCol = param.refs[0] if isinstance(param.refs, tuple) \
                else param.refs if isinstance(param.refs, int) else None
            refSep = param.refs[1] if isinstance(param.refs,
                                                    tuple) else ';'
            sigCol = None if not isinstance(sign, tuple) else sign[0]
            dir_col = None
            dir_val = None
            dir_sep = None

            if isinstance(is_directed, tuple):

                dir_col = is_directed[0]
                dir_val = is_directed[1]
                dir_sep = is_directed[2] if len(is_directed) > 2 else None

            elif isinstance(sign, tuple):

                dir_col = sign[0]
                dir_val = sign[1:3]
                dir_val = dir_val if type(dir_val[
                    0]) in _const.SIMPLE_TYPES else common.flat_list(dir_val)
                dir_sep = sign[3] if len(sign) > 3 else None

            dir_val = common.to_set(dir_val)
            max_col = max(
                filter(
                    lambda i: i is not None, [
                        param.id_col_a,
                        param.id_col_b,
                        self.get_max(param.extra_edge_attrs),
                        self.get_max(param.extra_node_attrs_a),
                        self.get_max(param.extra_node_attrs_b),
                        refCol,
                        dir_col,
                        sigCol,
                        max(itertools.chain(
                            map(lambda x: x[0],
                                param.positive_filters),
                            [0])),
                        max(itertools.chain(
                            map(lambda x: x[0],
                                param.negative_filters),
                            [0]))
                    ]
                )
            )

            must_have_references = (
                settings.get('keep_noref') or
                param.must_have_references
            )
            self._log(
                'Resource `%s` %s have literature references '
                'for all interactions. Interactions without references '
                'will be dropped. You can alter this condition globally by '
                '`pypath.settings.keep_noref` or for individual resources '
                'by the `must_have_references` attribute of their '
                '`NetworkInput` object.' % (
                    param.name,
                    'must' if must_have_references else 'does not need to'
                ),
                1,
            )
            self._log(
                '`%s` must have references: %s' % (
                    param.name,
                    str(must_have_references)
                )
            )

            # iterating lines from input file
            lFiltered = 0
            rFiltered = 0
            tFiltered = 0
            readError = 0
            lnum = 0 # we need to define it here to avoid errors if the
                     # loop below runs zero cycles

            for lnum, line in enumerate(infile):

                if len(line) <= 1 or (lnum == 1 and param.header):
                    # empty lines
                    # or header row
                    continue

                if not isinstance(line, (list, tuple)):

                    if hasattr(line, 'decode'):
                        line = line.decode('utf-8')

                    line = line.strip('\n\r').split(param.separator)

                else:
                    line = [
                        x.replace('\n', '').replace('\r', '')
                        if hasattr(x, 'replace') else x for x in line
                    ]

                # in case line has less fields than needed
                if len(line) < max_col:

                    self._log(
                        'Line #%u has less than %u fields,'
                        ' skipping! :(\n' % (lnum, max_col),
                        5,
                    )

                    readError = 1
                    continue

                else:
                    # from here we are committed to process this row

                    # applying filters:
                    if self._filters(
                        line,
                        param.positive_filters,
                        param.negative_filters
                    ):

                        lFiltered += 1
                        continue

                    # reading names and attributes:
                    if is_directed and not isinstance(is_directed, tuple):
                        this_edge_dir = True

                    else:
                        this_edge_dir = self._process_direction(
                            line,
                            dir_col,
                            dir_val,
                            dir_sep,
                        )

                    refs = []
                    if refCol is not None:

                        if isinstance(line[refCol], (list, set, tuple)):

                            refs = line[refCol]

                        elif isinstance(line[refCol], int):

                            refs = (line[refCol],)

                        else:

                            refs = line[refCol].split(refSep)

                        refs = common.del_empty(list(set(refs)))

                    refs = pubmed_input.only_pmids(
                        [str(r).strip() for r in refs]
                    )

                    if len(refs) == 0 and must_have_references:
                        rFiltered += 1
                        continue

                    # to give an easy way for input definition:
                    if isinstance(param.ncbi_tax_id, int):
                        taxon_a = param.ncbi_tax_id
                        taxon_b = param.ncbi_tax_id

                    # to enable more sophisticated inputs:
                    elif isinstance(param.ncbi_tax_id, dict):

                        taxx = self.get_taxon(param.ncbi_tax_id, line)

                        if isinstance(taxx, tuple):
                            taxon_a = taxx[0]
                            taxon_b = taxx[1]

                        else:
                            taxon_a = taxon_b = taxx

                        taxdA = (
                            param.ncbi_tax_id['A']
                            if 'A' in param.ncbi_tax_id else
                            param.ncbi_tax_id
                        )
                        taxdB = (
                            param.ncbi_tax_id['B']
                            if 'B' in param.ncbi_tax_id else
                            param.ncbi_tax_id
                        )

                        if (('include' in taxdA and
                            taxon_a not in taxdA['include']) or
                            ('include' in taxdB and
                            taxon_b not in taxdB['include']) or
                            ('exclude' in taxdA and
                            taxon_a in taxdA['exclude']) or
                            ('exclude' in taxdB and
                            taxon_b in taxdB['exclude'])):

                            tFiltered += 1
                            continue

                    else:
                        taxon_a = taxon_b = self.ncbi_tax_id

                    if taxon_a is None or taxon_b is None:
                        tFiltered += 1
                        continue

                    stim = False
                    inh = False

                    if isinstance(sign, tuple):
                        stim, inh = self._process_sign(line[sign[0]], sign)

                    resource = (
                        line[param.resource]
                        if isinstance(param.resource, int) else
                        line[param.resource[0]].split(param.resource[1])
                        if isinstance(param.resource, tuple) else
                        param.resource
                    )

                    resource = common.to_set(resource)

                    _resources_secondary = tuple(
                        network_resources.resource.NetworkResource(
                            name = sec_res,
                            interaction_type = _resource.interaction_type,
                            data_model = _resource.data_model,
                            via = _resource.name,
                        )
                        for sec_res in resource
                        if sec_res != _resource.name
                    )

                    resource.add(param.name)

                    id_a = line[param.id_col_a]
                    id_b = line[param.id_col_b]
                    id_a = id_a.strip() if hasattr(id_a, 'strip') else id_a
                    id_b = id_b.strip() if hasattr(id_b, 'strip') else id_b

                    evidences = evidence.Evidences(
                        evidences = (
                            evidence.Evidence(
                                resource = _res,
                                references = refs,
                            )
                            for _res in
                            _resources_secondary + (_resource,)
                        )
                    )


                    new_edge = {
                        'id_a': id_a,
                        'id_b': id_b,
                        'id_type_a': param.id_type_a,
                        'id_type_b': param.id_type_b,
                        'entity_type_a': param.entity_type_a,
                        'entity_type_b': param.entity_type_b,
                        'source': resource,
                        'is_directed': this_edge_dir,
                        'references': refs,
                        'stim': stim,
                        'inh': inh,
                        'taxon_a': taxon_a,
                        'taxon_b': taxon_b,
                        'type': param.interaction_type,
                        'evidences': evidences,
                    }

                    # getting additional edge and node attributes
                    attrs_edge = self.get_attrs(
                        line,
                        param.extra_edge_attrs,
                        lnum,
                    )
                    attrs_node_a = self.get_attrs(
                        line,
                        param.extra_node_attrs_a,
                        lnum,
                    )
                    attrs_node_b = self.get_attrs(
                        line,
                        param.extra_node_attrs_b,
                        lnum,
                    )

                    if param.mark_source:

                        attrs_node_a[param.mark_source] = this_edge_dir

                    if param.mark_target:

                        attrs_node_b[param.mark_target] = this_edge_dir

                    # merging dictionaries
                    node_attrs = {
                        'attrs_node_a': attrs_node_a,
                        'attrs_node_b': attrs_node_b,
                        'attrs_edge': attrs_edge,
                    }
                    new_edge.update(node_attrs)

                if readError != 0:

                    self._log(
                        'Errors occured, certain lines skipped.'
                        'Trying to read the remaining.\n',
                        5,
                    )
                    readError = 1

                edge_list.append(new_edge)

            if hasattr(infile, 'close'):

                infile.close()

            # ID translation of edges
            edge_list_mapped = self._map_list(
                edge_list,
                expand_complexes = expand_complexes,
            )

            self._log(
                '%u lines have been read from %s, '
                '%u links after mapping; '
                '%u lines filtered by filters; '
                '%u lines filtered because lack of references; '
                '%u lines filtered by taxon filters.' %
                (
                    lnum - 1,
                    param.input,
                    len(edge_list_mapped),
                    lFiltered,
                    rFiltered,
                    tFiltered,
                )
            )

            if reread or redownload:
                pickle.dump(edge_list_mapped, open(edges_cache, 'wb'), -1)
                self._log('ID translated edge list saved to %s' % edges_cache)

        else:

            self._log(
                'Previously ID translated edge list '
                'has been loaded from `%s`.' % edges_cache
            )

        if keep_raw:

            self.data[param.name] = edge_list_mapped

        self.raw_data = edge_list_mapped


    def signaling_proteins_list(self):
        """
        Compiles a list of signaling proteins (as opposed to other
        proteins like metabolic enzymes, matrix proteins, etc), by
        looking up a few simple keywords in short description of GO
        terms.
        """

        terms   = go_input.go_terms_quickgo()
        goannot = go_input.go_annotations_goa()

        gosig = set([])

        for term, name in iteritems(terms['P']):

            if 'signal' in name or 'regulat' in name:
                gosig.add(term)

        upsig = set([])

        if 'proteome' not in self.lists:
            self.proteome_list()

        for up, term in iteritems(goannot['P']):

            if len(term & gosig):
                upsig.add(up)

        spsig = set([])

        for u in upsig:
            spsig.update(set(mapping.map_name(
                u, 'uniprot', 'uniprot', ncbi_tax_id = self.ncbi_tax_id)))

        upsig = spsig & set(self.lists['proteome'])

        self.lists['sig'] = list(upsig)

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
            :py:attr:`pypath.input_formats.NetworkInput.input`.
        """

        data_formats.intogen_cancer.input = intogen_file
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
            :py:attr:`pypath.input_formats.NetworkInput.input`.
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


    def entities_by_resources(self):
        """
        Returns a dict of sets with resources as keys and sets of entity IDs
        as values.
        """

        results = collections.defaultdict(set)

        for v in self.graph.vs:

            for resource in v['sources']:

                result[resource].add(v['name'])

        return result


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
            :py:attr:`python.data_formats.ReadList.input` attribute and
            defined in :py:mod:`pypath.inputs`.
        """

        _input = None

        if settings.__class__.__name__ != "ReadList":
            self._log(
                'No proper input file definition! `settings`'
                'should be a `readList` instance',
                -5,
            )
            return None

        if (
            isinstance(settings.input, str) and
            inputs.get_method(settings.input)
        ):

            toCall = inputs.get_method(settings.input)
            _input = toCall(**kwargs)

        elif not os.path.isfile(settings.input):

            self._log('%s: No such file! :(\n' % settings.input, -5)
            return None

        else:
            _input = settings.input

        original_name_type = settings.id_type
        default_name_type = self.default_name_type[settings.entity_type]
        mapTbl = ''.join([original_name_type, "_", default_name_type])

        if type(_input) in _const.CHAR_TYPES and os.path.isfile(_input):

            _input = curl.Curl(_input, large = True).result

            #codecs.open(_input, encoding='utf-8', mode='r')

        if _input is None:
            self._log('Could not find file or input function.', -5)
            return None

        self._log("%s opened..." % settings.input)
        # finding the largest referred column number,
        # to avoid references out of index
        max_col = max([settings.id_col, self.get_max(settings.extra_attrs)])
        # iterating lines from input file
        lnum = 1
        readError = 0
        item_list = []

        for line in _input: # XXX: Could use enumerate(_input) instead of lnum

            if len(line) == 0 or (lnum == 1 and settings.header):
                # empty lines
                # or header row
                lnum += 1
                continue

            if type(line) in _const.CHAR_TYPES:
                line = line.rstrip().split(settings.separator)

            # in case line has less fields than needed
            if len(line) < max_col:
                self._log(
                    'Line #%u has less than %u fields! :(\n' % (lnum, max_col)
                )
                readError = 1
                break

            else:

                # reading names and attributes
                try:
                    newItem = {
                        "name": line[settings.id_col],
                       "id_type": settings.id_type,
                       "type": settings.entity_type,
                       "source": settings.name
                    }

                except:

                    print(line)
                    self._log(
                        'Wrong name column indexes (%u and %u), '
                        'or wrong separator (%s)? Line #%u' % (
                            settings.id_col,
                            settings.separator,
                            lnum,
                        )
                    )
                    readError = 1
                    break

                # getting additional attributes
                attrsItem = self.get_attrs(line, settings.extra_attrs, lnum)
                # merging dictionaries
                newItem.update(attrsItem)

            if readError != 0:
                break

            item_list.append(newItem)
            lnum += 1

        if hasattr(_input, 'close'):
            _input.close()

        item_list_mapped = self.map_list(item_list, single_list=True)
        item_list_mapped = list(set(item_list_mapped))
        self._log(
            '%u lines have been read from %s, %u  items after mapping' % (
                lnum,
                settings.input,
                len(item_list_mapped)
            )
        )
        self.lists[settings.name] = item_list_mapped


    def _map_list(
            self,
            lst,
            single_list = False,
            expand_complexes = True,
        ):
        """
        Maps the names from a list of edges or items (molecules).

        :arg list lst:
            List of items or edge dictionaries whose names have to be
            mapped.
        :arg bool single_list:
            Optional, ``False`` by default. Determines whether the
            provided elements are items or edges. This is, either calls
            :py:meth:`pypath.main.PyPath.map_edge` or
            :py:meth:`pypath.main.PyPath.map_item` to map the item
            names.
        :arg bool expand_complexes:
            Expand complexes, i.e. create links between each member of
            the complex and the interacting partner.

        :return:
            (*list*) -- Copy of *lst* with their elements' names mapped.
        """

        list_mapped = []

        if single_list:

            for item in lst:
                list_mapped += self._map_item(
                    item,
                    expand_complexes = expand_complexes,
                )

        else:

            for edge in lst:
                list_mapped += self._map_edge(
                    edge,
                    expand_complexes = expand_complexes,
                )

        return list_mapped


    def _map_item(self, item, expand_complexes = True):
        """
        Translates the name in *item* representing a molecule. Default
        name types are defined in
        :py:attr:`pypath.main.PyPath.default_name_type` If the mapping
        is unsuccessful, the item will be added to
        :py:attr:`pypath.main.PyPath.unmapped` list.

        :arg dict item:
            Item whose name is to be mapped to a default name type.
        :arg bool expand_complexes:
            Expand complexes, i.e. create links between each member of
            the complex and the interacting partner.

        :return:
            (*list*) -- The default mapped name(s) [str] of *item*.
        """

        # TODO: include
        default_id = mapping.map_name(
            item['name'], item['id_type'],
            self.default_name_type[item['type']],
            expand_complexes = expand_complexes,
        )

        if len(default_id) == 0:

            self.unmapped.append(item['name'])

        return default_id


    def _map_edge(self, edge, expand_complexes = True):
        """
        Translates the identifiers in *edge* representing an edge. Default
        name types are defined in
        :py:attr:`pypath.main.PyPath.default_name_type` If the mapping
        is unsuccessful, the item will be added to
        :py:attr:`pypath.main.PyPath.unmapped` list.

        :arg dict edge:
            Item whose name is to be mapped to a default name type.
        :arg bool expand_complexes:
            Expand complexes, i.e. create links between each member of
            the complex and the interacting partner.

        :return:
            (*list*) -- Contains the edge(s) [dict] with default mapped
            names.
        """

        edge_stack = []

        default_id_a = mapping.map_name(
            edge['id_a'],
            edge['id_type_a'],
            self.default_name_type[edge['entity_type_a']],
            ncbi_tax_id = edge['taxon_a'],
            expand_complexes = expand_complexes,
        )

        default_id_b = mapping.map_name(
            edge['id_b'],
            edge['id_type_b'],
            self.default_name_type[edge['entity_type_b']],
            ncbi_tax_id = edge['taxon_b'],
            expand_complexes = expand_complexes,
        )

        # this is needed because the possibility ambigous mapping
        # and expansion of complexes
        # one name can be mapped to multiple ones
        # this multiplies the nodes and edges
        # in case of proteins this does not happen too often
        for id_a, id_b in itertools.product(default_id_a, default_id_b):

            this_edge = copy_mod.copy(edge)
            this_edge['default_name_a'] = id_a
            this_edge['default_name_type_a'] = (
                self.default_name_type[edge['entity_type_a']]
            )

            this_edge['default_name_b'] = id_b
            this_edge['default_name_type_b'] = (
                self.default_name_type[edge['entity_type_b']]
            )
            edge_stack.append(this_edge)

        return edge_stack


    def combine_attr(self, lst, num_method = max):
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
        if type(lst[0]) in common.numeric_types and type(lst[1]) in common.numeric_types:
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
            return common.add_to_set(lst[0], lst[1])

        if isinstance(lst[1], set):
            return common.add_to_set(lst[1], lst[0])

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
        if (isinstance(lst[0], list) and type(lst[1]) in _const.SIMPLE_TYPES):

            if lst[1] in common.numeric_types or len(lst[1]) > 0:
                return common.add_to_list(lst[0], lst[1])

            else:
                return lst[0]

        if (isinstance(lst[1], list) and type(lst[0]) in _const.SIMPLE_TYPES):

            if lst[0] in common.numeric_types or len(lst[0]) > 0:
                return common.add_to_list(lst[1], lst[0])

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

    def collapse_by_name(self, graph = None):
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

        multiple = collections.Counter(graph.vs['name'])

        for name, count in iteritems(multiple):

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
        graph.vs['id_merge'] = list(xrange(graph.vcount()))

        # combining vertex attributes:
        vprim = graph.vs[primary]

        for attr in vprim.attributes():

            if attr not in {'name', 'id_merge', 'label'}:

                vprim[attr] = self.combine_attr(list(map(
                                lambda vid: graph.vs[vid][attr],
                                # combining from all nodes
                                nodes)))

        # moving edges of non primary vertices to the primary one
        self.copy_edges(nonprimary, primary, move = True, graph = graph)

        # deleting non primary vertices:
        toDel = [graph.vs.select(id_merge = i)[0].index for i in nonprimary]

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
        ses = dict(
            (
                # id_old of source vertices:
                s,
                # edges of current source node:
                set(
                    e.index for e in
                    itertools.chain(
                        graph.es.select(_source = s),
                        graph.es.select(_target = s),
                    )
                )
            )
            for s in sources
        )

        # collecting edges to be newly created
        toAdd = set([])

        for s, es in iteritems(ses):

            for eid in es:
                # the source edge:
                e = graph.es[eid]

                new_source = target if e.source == s else e.source
                new_target = target if e.target == s else e.target

                es0 = graph.es.select(
                    _source = new_source,
                    _target = new_target,
                )
                es1 = ()
                if not graph.is_directed():
                    es1 = graph.es.select(
                        _source = new_target,
                        _target = new_source,
                    )
                # looking up if target edge already exists:

                if not len(list(itertools.chain(es0, es1))):

                    toAdd.add((new_source, new_target))

        # creating new edges
        graph.add_edges(toAdd)
        nvidt = graph.vs.select(id_old = target)[0].index

        # copying attributes:

        for ovids, es in iteritems(ses):

            for oeid in es:

                # this is the current source edge:
                e = graph.es.select(id_old = oeid)[0]
                # this is the index of the other (peer) node:
                new_source = (
                    nvidt
                        if graph.vs[e.source]['id_old'] == ovids else
                    e.source
                )
                new_target = (
                    nvidt
                        if graph.vs[e.target]['id_old'] == ovids else
                    e.target
                )
                nvids = graph.vs.select(id_old = ovids)[0].index
                nvid_other = e.target if e.source == nvids else e.source

                # looking up new edge:
                es0 = graph.es.select(
                    _source = new_source,
                    _target = new_target,
                )
                es1 = ()

                if not graph.is_directed():

                    es1 = graph.es.select(
                        _source = new_target,
                        _target = new_source,
                    )

                es_all = list(itertools.chain(es0, es1))

                for new_edge in es_all:

                    # old direction:
                    d = e['dirs']
                    # dict from old names to new ones
                    # the peer does no change, only s->t
                    ids = {
                        graph.vs[nvids]['name']:
                            graph.vs[nvidt]['name'],
                        graph.vs[nvid_other]['name']:
                            graph.vs[nvid_other]['name'],
                    }

                    # copying directions and signs:
                    new_dirs = d.translate(ids)

                    if new_edge['dirs']:
                        new_dirs.merge(new_edge['dirs'])

                    new_edge['dirs'] = new_dirs

                    # copying `refs_by_dir`
                    new_edge['refs_by_dir'] = self._translate_refsdir(
                        e['refs_by_dir'], ids,
                    )

                    # copying further attributes:
                    for eattr in e.attributes():

                        if eattr not in {'dirs', 'refs_by_dir', 'id_old'}:

                            new_edge[eattr] = self.combine_attr([
                                new_edge[eattr],
                                e[eattr]
                            ])

                # in case we want to delete old edges:
                toDel.add(e.index)

        if move:

            graph.delete_edges(list(toDel))

        # removing temporary attributes
        del graph.es['id_old']
        del graph.vs['id_old']


    def delete_by_organism(self, organisms_allowed = None):
        """
        Removes the proteins of all organisms which are not given in
        *tax*.

        :arg list,set organisms_allowed:
            List of NCBI Taxonomy IDs [int] of the organism(s) that are
            to be kept.
        """

        g = self.graph

        organisms_allowed = organisms_allowed or {self.ncbi_tax_id}

        to_delete = [
            v.index for v in g.vs
            if v['ncbi_tax_id'] not in organisms_allowed
        ]

        g.delete_vertices(to_delete)
        self.update_vname()
        self.update_db_dict()


    def delete_unknown(
            self,
            organisms_allowed = None,
            entity_type = 'protein',
            default_name_type = None,
        ):
        """
        Removes those items which are not in the list of all default
        IDs of the organisms. By default, it means to remove all protein
        nodes not having a human UniProt ID.

        :arg str typ:
            Optional, ``'protein'`` by default. Determines the molecule
            type. These can be ``'protein'``, ``'drug'``, ``'lncrna'``,
            ``'mirna'`` or any other type defined in
            :py:attr:`pypath.main.PyPath.default_name_type`.
        :arg str default_name_type:
            Optional, ``None`` by default. The default name type for the
            given molecular species. If none is specified takes it from
            :py:attr:`pypath.main.PyPath.default_name_type` (e.g.: for
            ``'protein'``, default is ``'uniprot'``).
        :arg set organisms_allowed:
            NCBI Taxonomy identifiers [int] of the organisms allowed in
            the network.
        """

        g = self.graph

        if not default_name_type:

            default_name_type = self.default_name_type[entity_type]

        organisms_allowed = organisms_allowed or {self.ncbi_tax_id}

        self.update_vname()

        self._log('Checking network components against reference lists.')

        names = g.vs['id_type']
        names = [i for i, j in enumerate(names) if j == default_name_type]
        entity_types = g.vs['type']
        entity_types = [
            i for i, j in enumerate(entity_types)
            if j == entity_type
        ]
        organisms = g.vs['ncbi_tax_id']
        organisms = [
            i for i, j in enumerate(organisms)
            if j in organisms_allowed
        ]
        vertices = list(set(names) & set(entity_types) & set(organisms))
        names_selected = [g.vs[i]['name'] for i in vertices]

        names_to_delete = set.intersection(
            *(
                reflists.is_not(
                    names = names_selected,
                    id_type = default_name_type,
                    ncbi_tax_id = ncbi_tax_id,
                )
                for ncbi_tax_id in organisms_allowed
            )
        )

        vertices_to_delete = [self.nodDct[n] for n in names_to_delete]

        g.delete_vertices(vertices_to_delete)

        self.update_vname()


    def clean_graph(self, organisms_allowed = None):
        """
        Removes multiple edges, unknown molecules and those from wrong
        taxon. Multiple edges will be combined by
        :py:meth:`pypath.main.PyPath.combine_attr` method.
        Loops will be deleted unless the attribute
        :py:attr:`pypath.main.PyPath.loops` is set to ``True``.

        :arg set organisms_allowed:
            NCBI Taxonomy identifiers [int] of the organisms allowed
            in the network.
        """

        self._log('Removing duplicate edges...')
        g = self.graph

        if not g.is_simple():
            g.simplify(loops=not self.loops, multiple=True,
                       combine_edges=self.combine_attr)

        self._log(
            'After duplicate edge removal: '
            'number of nodes: %u, edges: %u' % (
                self.graph.vcount(),
                self.graph.ecount(),
            )
        )

        self.delete_unmapped()

        self._log(
            'After removing unmapped nodes: '
            'number of nodes: %u, edges: %u' % (
                self.graph.vcount(),
                self.graph.ecount(),
            )
        )

        self.delete_by_organism(organisms_allowed = organisms_allowed)

        self._log(
            'After removing unknown organism nodes: '
            'number of nodes: %u, edges: %u' % (
                self.graph.vcount(),
                self.graph.ecount(),
            )
        )

        self.delete_unknown(organisms_allowed = organisms_allowed)

        x = g.vs.degree()
        zeroDeg = [i for i, j in enumerate(x) if j == 0]
        g.delete_vertices(zeroDeg)

        self._log(
            'After removing zero degree nodes: '
            'number of node: %u, edges: %u' % (
                self.graph.vcount(),
                self.graph.ecount(),
            )
        )

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


    def _add_update_vertex(
            self,
            default_attrs,
            original_name,
            original_name_type,
            extra_attrs = {},
            add = False,
        ):
        """
        Updates the attributes of one node in the (undirected) network.
        Optionally it creates a new node and sets the attributes, but it
        is not efficient as :py:mod:`igraph` needs to reindex vertices
        after this operation, so better to create new nodes in batches.

        :arg dict default_attrs:
            The attribute dictionary of the node to be updated/created.
        :arg str original_name:
            Original node name (e.g.: UniProt ID).
        :arg str original_name_type:
            The original node name type (e.g.: for the previous example,
            this would be ``'uniprot'``).
        :arg dict extra_attrs:
            Optional, ``{}`` by default. Contains any extra attributes
            for the node to be updated.
        :arg bool add:
            Optional, ``False`` by default. If set to ``True`` and the
            node is not in the network, it will be created. Otherwise,
            in such case it will raise an error message.
        """

        keep_original_names = settings.get('network_keep_original_names')

        g = self.graph
        g.vs._reindex_names()

        if not default_attrs['name'] in g.vs['name']:

            if not add:
                self._log('Failed to add some vertices', -5)
                return False

            n = g.vcount()
            g.add_vertices(1)

            # only keep track of original names if they are strings
            # not, for example, complexes
            if (
                keep_original_names and
                isinstance(original_name, str)
            ):
                g.vs[n]['original_names'] = {
                    original_name: original_name_type,
                }

            this_node = g.vs[g.vs._name_index[default_attrs['name']]]

        else:

            this_node = g.vs[g.vs._name_index[default_attrs['name']]]

            if this_node['original_names'] is None:

                this_node['original_names'] = {}

            # only keep track of original names if they are strings
            # not, for example, complexes
            if (
                keep_original_names and
                isinstance(original_name, str)
            ):

                this_node['original_names'][original_name] = (
                    original_name_type
                )

        if isinstance(default_attrs['name'], intera.Complex):

            default_attrs['type'] = 'complex'

        for key, value in iteritems(default_attrs):

            this_node[key] = value

        for key, value in iteritems(extra_attrs):

            if key not in g.vs.attributes():

                g.vs[key] = (
                    [[] for _ in xrange(self.graph.vcount())]
                        if isinstance(value, list) else
                    [None]
                )

            this_node[key] = self.combine_attr([this_node[key], value])


    def _add_update_edge(
            self,
            id_a,
            id_b,
            id_type_a,
            id_type_b,
            entity_type_a,
            entity_type_b,
            source,
            evidences,
            is_directed,
            refs,
            stim,
            inh,
            taxon_a,
            taxon_b,
            typ,
            extra_attrs = {},
            add = False,
        ):
        """
        Updates the attributes of one edge in the (undirected) network.
        Optionally it creates a new edge and sets the attributes, but it
        is not efficient as :py:mod:`igraph` needs to reindex edges
        after this operation, so better to create new edges in batches.

        :arg str id_a:
            Name of the source node of the edge to be added/updated.
        :arg str id_b:
            Name of the source node of the edge to be added/updated.
        :arg set source:
            Or [list], contains the names [str] of the resources
            supporting that edge.
        :arg pypath.evidence.Evidence evidence:
            A ``pypath.evidence.Evidence`` object.
        :arg bool is_directed:
            Whether if the edge is directed or not.
        :arg set refs:
            Or [list], contains the instances of the references
            :py:class:`pypath.refs.Reference` for that edge.
        :arg bool stim:
            Whether the edge is stimulatory or not.
        :arg bool inh:
            Whether the edge is inhibitory or note
        :arg int taxon_a:
            NCBI Taxonomic identifier of the source molecule.
        :arg int taxon_b:
            NCBI Taxonomic identifier of the target molecule.
        :arg str typ:
            The type of interaction (e.g.: ``'PPI'``)
        :arg dict extra_attrs:
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

        edge = self.edge_exists(id_a, id_b)

        if isinstance(edge, list):

            if not add:
                # normally we create all new edges before calling this method
                # for efficiency reasons; hence if ``add`` is False and the
                # edge does not exist something must have gone wrong earlier
                self._log('Failed to add some edges', -5)
                aid = self.nodDct[id_a]
                bid = self.nodDct[id_b]
                a = g.get_eid(aid, bid, error=False)
                b = g.get_eid(aid, bid, error=False)
                self.failed_edges.append([edge, id_a, id_b, aid, bid, a, b])

                return False

            g.add_edge(edge[0], edge[1])
            edge = self.edge_exists(id_a, id_b)

        # assigning source:
        self.add_set_eattr(edge, 'sources', source)

        # adding references:
        # if len(refs) > 0:
        refs = [_refs.Reference(pmid) for pmid in refs]
        self.add_list_eattr(edge, 'references', refs)

        entity_a = entity_mod.Entity(
            identifier = id_a,
            id_type = id_type_a,
            entity_type = entity_type_a,
            taxon = taxon_a,
        )
        entity_b = entity_mod.Entity(
            identifier = id_b,
            id_type = id_type_b,
            entity_type = entity_type_b,
            taxon = taxon_b,
        )

        attrs = interaction.Interaction(
            a = entity_a,
            b = entity_b,
        )

        # updating references-by-source dict:
        sources = (
            source
                if isinstance(source, (tuple, set, list)) else
            (source,)
        )

        for src in sources:
            self.add_grouped_set_eattr(edge, 'refs_by_source', src, refs)

        # updating refrences-by-type dict:
        self.add_grouped_set_eattr(edge, 'refs_by_type', typ, refs)

        # setting directions:
        if not g.es[edge]['dirs']:
            g.es[edge]['dirs'] = Direction(id_a, id_b)

        if is_directed:
            g.es[edge]['dirs'].set_dir((id_a, id_b), source)
            # updating references-by-direction dict:
            self.add_grouped_set_eattr(
                edge,
                'refs_by_dir',
                (id_a, id_b),
                refs,
            )
            attrs.add_evidence(
                evidence = evidences,
                direction = (entity_a, entity_b),
            )
        else:
            g.es[edge]['dirs'].set_dir('undirected', source)
            self.add_grouped_set_eattr(
                edge,
                'refs_by_dir',
                'undirected',
                refs,
            )
            attrs.add_evidence(
                evidence = evidences,
                direction = 'undirected',
            )

        # setting signs:
        if stim:
            g.es[edge]['dirs'].set_sign((id_a, id_b), 'positive', source)
            attrs.add_evidence(
                evidence = evidences,
                direction = (entity_a, entity_b),
                effect = 1,
            )

        if inh:
            g.es[edge]['dirs'].set_sign((id_a, id_b), 'negative', source)
            attrs.add_evidence(
                evidence = evidences,
                direction = (entity_a, entity_b),
                effect = -1,
            )

        # adding interaction attributes (this new kind of object either will
        # replace the igraph based network representation or is a temporary
        # solution and something else will replace them):
        if not isinstance(g.es[edge]['attrs'], interaction.Interaction):

            g.es[edge]['attrs'] = attrs

        else:

            g.es[edge]['attrs'] += attrs

        # updating sources-by-type dict:
        self.add_grouped_set_eattr(edge, 'sources_by_type', typ, source)
        # adding type:
        self.add_list_eattr(edge, 'type', typ)

        # adding extra attributes:
        for key, value in iteritems(extra_attrs):

            if key not in g.es.attributes():
                g.es[key] = ([[] for _ in xrange(self.graph.ecount())]
                             if isinstance(value, list) else [None])

            g.es[edge][key] = self.combine_attr([g.es[edge][key], value])


    def add_list_eattr(self, edge, attr, value):
        """
        Merges (or creates) a given edge attribute as [list].

        :arg int edge:
            Edge index where the given attribute value is to be merged
            or created.
        :arg str attr:
            The name of the attribute. If such attribute does not exist
            in the network edges, it will be created on all edges (as an
            empty [list], *value* will only be assigned to the given
            *edge*).
        :arg list value:
            The value of the attribute to be assigned/merged.
        """

        value = value if isinstance(value, list) else [value]
        e = self.graph.es[edge]

        if attr not in self.graph.es.attributes():
            self.graph.es[attr] = [[] for _ in xrange(0, self.graph.ecount())]

        if e[attr] is None:
            e[attr] = []

        elif not isinstance(e[attr], list):
            e[attr] = [e[attr]]

        e[attr] = common.unique_list(e[attr] + value)

    def add_set_eattr(self, edge, attr, value):
        """
        Merges (or creates) a given edge attribute as [set].

        :arg int edge:
            Edge index where the given attribute value is to be merged
            or created.
        :arg str attr:
            The name of the attribute. If such attribute does not exist
            in the network edges, it will be created on all edges (as an
            empty [set], *value* will only be assigned to the given
            *edge*).
        :arg set value:
            The value of the attribute to be assigned/merged.
        """

        value = common.to_set(value)

        e = self.graph.es[edge]

        if attr not in self.graph.es.attributes():

            self.graph.es[attr] = [
                set()
                for _ in xrange(self.graph.ecount())
            ]
        if e[attr] is None:
            e[attr] = set()

        elif not isinstance(e[attr], set):

            e[attr] = common.to_set(e[attr])

        e[attr].update(value)

    def add_grouped_eattr(self, edge, attr, group, value):
        """
        Merges (or creates) a given edge attribute as [dict] of [list]
        values.

        :arg int edge:
            Edge index where the given attribute value is to be merged
            or created.
        :arg str attr:
            The name of the attribute. If such attribute does not exist
            in the network edges, it will be created on all edges (as an
            empty [dict], *value* will only be assigned to the given
            *edge* and *group*).
        :arg str group:
            The key of the attribute dictionary where *value* is to be
            assigned.
        :arg list value:
            The value of the attribute to be assigned/merged.
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

        e[attr][group] = common.unique_list(e[attr][group] + value)

    def add_grouped_set_eattr(self, edge, attr, group, value):
        """
        Merges (or creates) a given edge attribute as [dict] of [set]
        values.

        :arg int edge:
            Edge index where the given attribute value is to be merged
            or created.
        :arg str attr:
            The name of the attribute. If such attribute does not exist
            in the network edges, it will be created on all edges (as an
            empty [dict], *value* will only be assigned to the given
            *edge* and *group*).
        :arg str group:
            The key of the attribute dictionary where *value* is to be
            assigned.
        :arg set value:
            The value of the attribute to be assigned/merged.
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
        Converts a copy of *graph* undirected *igraph.Graph* object to a
        directed one. By default it converts the current network
        instance in :py:attr:`pypath.main.PyPath.graph` and places the
        copy of the directed instance in
        :py:attr:`pypath.main.PyPath.dgraph`.

        :arg igraph.Graph graph:
            Optional, ``None`` by default. Undirected graph object. If
            none is passed, takes the current undirected network
            instance and saves the directed network under the attribute
            :py:attr:`pypath.main.PyPath.dgraph`. Otherwise, the
            directed graph will be returned instead.
        :arg bool conv_edges:
            Optional, ``False`` by default. Whether to convert
            undirected edges (those without explicit direction
            information) to an arbitrary direction edge or
            a pair of opposite edges. Otherwise those will be deleted.
        :arg bool mutual:
            Optional, ``False`` by default. If *conv_edges* is ``True``,
            whether to convert the undirected edges to a single,
            arbitrary directed edge, or a pair of opposite directed
            edges.
        :arg bool ret:
            Optional, ``False`` by default. Whether to return the
            directed graph instance, or not. If a *graph* is provided,
            its directed version will be returned anyway.

        :return:
            (*igraph.Graph*) -- If *graph* is passed or *ret* is
            ``True``, returns the copy of the directed graph. otherwise
            returns ``None``.
        """

        self._log('Creating directed network object.')

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

        self._log('Directed igraph object created.')

        if graph or ret:
            return d


    def new_edges(self, edges):
        """
        Adds new edges from any iterable of edges to the undirected
        graph. Basically, calls :py:meth:`igraph.Graph.add_edges`.

        :arg list edges:
            Contains the edges that are to be added to the network.
        """

        self.graph.add_edges(list(edges))


    def new_nodes(self, nodes):
        """
        Adds new nodes from any iterable of nodes to the undirected
        graph. Basically, calls :py:meth:`igraph.Graph.add_vertices`.

        :arg list nodes:
            Contains the nodes that are to be added to the network.
        """

        self.graph.add_vertices(list(nodes))

    def edge_exists(self, id_a, id_b):
        """
        Returns a tuple of vertex indices if edge doesn't exist,
        otherwise, the edge ID. Not sensitive to direction.

        :arg str id_a:
            Name of the source node.
        :arg str id_b:
            Name of the target node.

        :return:
            (*int*) -- The edge index, if exists such edge. Otherwise,
            [tuple] of [int] corresponding to the node IDs.
        """

        if not hasattr(self, 'nodDct'):
            self.update_vname()

        nodes = [self.nodDct[id_a], self.nodDct[id_b]]
        edge = self.graph.get_eid(nodes[0], nodes[1], error = False)

        if edge != -1:
            return edge

        else:
            nodes.sort()
            return nodes


    def edge_names(self, e):
        """
        Returns the node names of a given edge.

        :arg int e:
            The edge index.

        :return:
            (*tuple*) -- Contains the source and target node names of
            the edge [str].
        """

        if isinstance(e, int):
            e = self.graph.es[e]

        return (self.graph.vs[e.source]['name'],
                self.graph.vs[e.target]['name'])


    def node_exists(self, name):
        """
        Checks if a node exists in the (undirected) network.

        :arg str name:
            The name of the node to be searched.

        :return:
            (*bool*) -- Whether the node exists in the network or not.
        """

        if not hasattr(self, 'nodInd'):
            self.update_vname()

        return name in self.nodInd

    def names2vids(self, names):
        """
        From a list of node names, returns their corresponding indices.

        :arg list names:
            Contains the node names [str] for which the IDs are to be
            searched.

        :return:
            (*list*) -- The queried node IDs [int].
        """

        vids = []

        if not hasattr(self, 'nodInd'):
            self.update_vname()

        for n in names:

            if n in self.nodInd:
                vids.append(self.nodDct[n])

        return vids

    # XXX: this function name may lead to some confusion with get_edge() method
    #      also, leading underscore makes the method private (is this intended?)
    def _get_edge(self, nodes):
        """
        Returns the edge index only if there is such an edge from
        *nodes*[0] to *nodes*[1], returns ``False```otherwise (e.g.: if
        edge exists in the opposite direction, no such edge exists or
        any of the vertex ids doesn't exist). To find edges regardless
        of their direction, see
        :py:meth:`pypath.main.PyPath.edge_exists`.

        :arg tuple nodes:
            Or [list], contains the node IDs [int] where the first
            element is the source and the second one the target.

        :return:
            (*int*) -- The edge ID if it exists, ``False`` otherwise.
        """

        g = self.graph

        try:
            e = g.get_eid(nodes[0], nodes[1])
            return e

        except:
            return False

    def straight_between(self, id_a, id_b):
        """
        Finds an edge between the provided node names.

        :arg str id_a:
            The name of the source node.
        :arg str id_b:
            The name of the target node.

        :return:
            (*int*) -- The edge ID. If the edge doesn't exist, returns
            [list] with the node indices [int].
        """

        nodNm = sorted([id_a, id_b])
        nodes = [self.graph.vs['name'].index(nodNm[0]),
                 self.graph.vs['name'].index(nodNm[1])]
        edge = self._get_edge(nodes)

        if isinstance(edge, int):
            return edge

        else:
            return nodes

    # XXX: Not sure if the intended behavior, according to old description:
    # """
    # Returns all edges between two given vertex names. Similar to
    # straight_between(), but checks both directions, and returns
    # list of edge ids in [undirected, straight, reversed] format,
    # for both id_a -> id_b and id_b -> id_a edges.
    # """
    # Just returns A SINGLE edge ID assigned according to the 'dirs' attribute
    # on a position of the dict and list

    def all_between(self, id_a, id_b):
        """
        Checks for any edges (in any direction) between the provided
        nodes.

        :arg str id_a:
            The name of the source node.
        :arg str id_b:
            The name of the target node.

        :return:
            (*dict*) -- Contains information on the directionality of
            the requested edge. Keys are ``'ab'`` and ``'ba'``, denoting
            the straight/reverse directionalities respectively. Values
            are [list] whose elements are the edge ID or ``None``
            according to the existance of that edge in the following
            categories: undirected, straight and reverse (in that
            order).
        """

        g = self.graph
        edges = {'ab': [None, None, None], 'ba': [None, None, None]}
        eid = self.edge_exists(id_a, id_b)

        if isinstance(eid, int):

            if g.es[eid]['dirs'].get_dir('undirected'):
                edges['ab'][0] = eid
                edges['ba'][0] = eid

            if g.es[eid]['dirs'].get_dir((id_a, id_b)):
                edges['ab'][1] = eid
                edges['ba'][2] = eid

            if g.es[eid]['dirs'].get_dir((id_b, id_a)):
                edges['ab'][2] = eid
                edges['ba'][1] = eid

        return edges

    def get_node_pair(self, id_a, id_b, directed=False):
        """
        Retrieves the node IDs from a pair of node names.

        :arg str id_a:
            Name of the source node.
        :arg str id_b:
            Name of the target node.
        :arg bool directed:
            Optional, ``False`` by default. Whether to return the node
            indices from the directed or undirected graph.

        :return:
            (*tuple*) -- The pair of node IDs of the selected graph.
            If not found, returns ``False``.
        """

        if not hasattr(self, 'nodDct'):
            self.update_vname()

        g = self._directed if directed else self._undirected
        nodDct = self.dnodDct if directed else self.nodDct
        nodes = [id_a, id_b] if not directed else sorted([id_a, id_b])

        try:
            nodeA = nodDct[nodes[0]]
            nodeB = nodDct[nodes[1]]
            return (nodeA, nodeB)

        except:
            return False

    def update_attrs(self):
        """
        Updates the node and edge attributes. Note that no data is
        donwloaded, mainly updates the dictionaries of attributes
        :py:attr:`pypath.main.PyPath.edgeAttrs` and
        :py:attr:`pypath.main.PyPath.vertexAttrs` containing the
        attributes names and their correspoding types and initializes
        such attributes in the network nodes/edges if they weren't.
        """

        for attr in self.graph.vs.attributes():
            types = list(
                set([type(x) for x in self.graph.vs[attr] if x is not None]))

            if len(types) > 1:
                self._log(
                    'Vertex attribute `%s` has multiple types of'
                    ' values: %s' % (
                        attr,
                        ', '.join(x.__name__ for x in types)
                    )
                )

            elif len(types) == 0:
                self._log('Vertex attribute `%s` has only None values' % attr)

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
                self._log(
                    'Edge attribute `%s` has multiple types of'
                    ' values: %s' % (
                        attr,
                        ', '.join(x.__name__ for x in types)
                    )
                )

            elif len(types) == 0:

                self._log('Edge attribute `%s` has only None values' % attr)

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
        Fills all vertices attribute *attr* with its default type (if
        such attribute value is ``None``), creates [list] if in
        :py:attr:`pypath.main.PyPath.vertexAttrs` such attribute is
        registered as [list].

        :arg str attr:
            The attribute name to be initialized on the network
            vertices.
        """

        # XXX: Doesn't handle potential KeyError
        for v in self.graph.vs:

            if v[attr] is None:
                v[attr] = self.vertexAttrs[attr]()

            if self.vertexAttrs[attr] is list and type(v[
                    attr]) in _const.SIMPLE_TYPES:
                v[attr] = [v[attr]] if len(v[attr]) > 0 else []

    def init_edge_attr(self, attr):
        """
        Fills all edges attribute *attr* with its default type (if
        such attribute value is ``None``), creates [list] if in
        :py:attr:`pypath.main.PyPath.edgeAttrs` such attribute is
        registered as [list].

        :arg str attr:
            The attribute name to be initialized on the network edges.
        """

        # XXX: Doesn't handle potential KeyError
        for e in self.graph.es:

            if e[attr] is None:
                e[attr] = self.edgeAttrs[attr]()

            if (self.edgeAttrs[attr] is list or
                self.edgeAttrs[attr] is set) and type(e[
                    attr]) in _const.SIMPLE_TYPES:

                e[attr] = [e[attr]] if (
                    type(e[attr]) not in _const.CHAR_TYPES or
                    len(e[attr]) > 0) else []

            if self.edgeAttrs[attr] is set and type(e[attr]) is list:

                e[attr] = set(e[attr])

    def _add_network(self, edge_list = False, regulator = False):
        """
        Adds edges to the network from *edge_list* obtained from file or
        other input method. If none is passed, checks for such data in
        :py:attr:`pypath.main.PyPath.raw_data`.

        :arg str edge_list:
            Optional, ``False`` by default. The source name of the list
            of edges to be added. This must have been loaded previously
            (e.g.: with :py:meth:`pypath.main.PyPath.read_data_file`).
            If none is passed, loads the data directly from
            :py:attr:`pypath.main.PyPath.raw_data`.
        :arg bool regulator:
            Optional, ``False`` by default. If set to ``True``, non
            previously existing nodes, will not be added (and hence, the
            edges involved).
        """

        self._log('Adding preprocessed edge list to existing network.')

        g = self.graph

        if not edge_list:

            if self.raw_data is not None:
                edge_list = self.raw_data

            else:
                self._log('_add_network(): No data, nothing to do.')
                return True

        if isinstance(edge_list, str):

            if edge_list in self.data:
                edge_list = self.data[edge_list]

            else:
                self._log(
                    '`%s` looks like a source name, but no data '
                    'available under this name.' % edge_list
                )

                return False

        nodes = []
        edges = []
        # adding nodes and edges first in bunch,
        # to avoid multiple reindexing by igraph
        self.update_vname()
        prg = Progress(
            total=len(edge_list), name="Processing nodes", interval=50)

        for e in edge_list:
            aexists = self.node_exists(e["default_name_a"])
            bexists = self.node_exists(e["default_name_b"])

            if not aexists and (not regulator or bexists):
                nodes.append(e["default_name_a"])

            if not bexists and not regulator:
                nodes.append(e["default_name_b"])

            prg.step()

        prg.terminate()
        self.new_nodes(set(nodes))
        self._log('New nodes have been created (%u)' % len(nodes))
        self.update_vname()
        prg = Progress(
            total=len(edge_list), name='Processing edges', interval=50)

        for e in edge_list:
            aexists = self.node_exists(e["default_name_a"])
            bexists = self.node_exists(e["default_name_b"])

            if aexists and bexists:

                edge = self.edge_exists(
                    e["default_name_a"],
                    e["default_name_b"],
                )

                if isinstance(edge, list):
                    edges.append(tuple(edge))

                prg.step()

        prg.terminate()
        self.new_edges(set(edges))
        self._log('New edges have been created')
        self._log('Introducing new node and edge attributes...')
        prg = Progress(
            total = len(edge_list),
            name = 'Processing attributes',
            interval = 30,
        )
        nodes_updated = []
        self.update_vname()

        for e in edge_list:
            # adding new node attributes

            if e['default_name_a'] not in nodes_updated:
                default_attrs = {
                    'name': e['default_name_a'],
                    'label': e['default_name_a'],
                    'id_type': e['default_name_type_a'],
                    'type': e['entity_type_a'],
                    'ncbi_tax_id': e['taxon_a'],
                }
                self._add_update_vertex(
                    default_attrs,
                    e['id_a'],
                    e['id_type_a'],
                    e['attrs_node_a'],
                )
                nodes_updated.append(e['default_name_a'])

            if e['default_name_b'] not in nodes_updated:
                default_attrs = {
                    'name': e['default_name_b'],
                    'label': e['default_name_b'],
                    'id_type': e['default_name_type_b'],
                    'type': e['entity_type_b'],
                    'ncbi_tax_id': e['taxon_b']
                }
                self._add_update_vertex(
                    default_attrs,
                    e['id_b'],
                    e['id_type_b'],
                    e['attrs_node_b'],
                )
                nodes_updated.append(e['default_name_b'])

            # adding new edge attributes
            self._add_update_edge(
                e['default_name_a'],
                e['default_name_b'],
                e['default_name_type_a'],
                e['default_name_type_b'],
                e['entity_type_a'],
                e['entity_type_b'],
                e['source'],
                e['evidences'],
                e['is_directed'],
                e['references'],
                e['stim'],
                e['inh'],
                e['taxon_a'],
                e['taxon_b'],
                e['type'],
                e['attrs_edge'],
            )

            prg.step()

        self._log(
            'New network resource added, current number '
            'of nodes: %u, edges: %u.' % (
                self.graph.vcount(),
                self.graph.ecount()
            )
        )

        prg.terminate()
        self.raw_data = None
        self.update_attrs()


    def apply_list(self, name, node_or_edge = 'node'):
        """
        Creates vertex or edge attribute based on a list.

        :arg str name:
            The name of the list to be added as attribute. Must have
            been previously loaded with
            :py:meth:`pypath.main.PyPath.load_list` or other methods.
            See description of :py:attr:`pypath.main.PyPath.lists`
            attribute for more information.
        :arg str node_or_edge:
            Optional, ``'node'`` by default. Whether the attribute list
            is to be added to the nodes or to the edges.
        """

        if name not in self.lists:
            self._log('No such list: %s' % name, -5)
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

    def merge_lists(self, id_a, id_b, name=None, and_or='and', delete=False,
                    func="max"): # XXX: kwarg func not used
        """
        Merges two lists from :py:attr:`pypat.main.PyPath.lists`.

        :arg str id_a:
            Name of the first list to be merged.
        :arg str id_b:
            Name of the second list to be merged.
        :arg str name:
            Optional, ``None`` by default. Specifies a new name for the
            merged list. If none is passed, name will be set to
            *id_a*_*id_b*.
        :arg str and_or:
            Optional, ``'and'`` by default. The logic operation perfomed
            in the merging: ``'and'`` performs an union, ``'or'`` for
            the intersection.
        :arg bool delete:
            Optional, ``False`` by default. Whether to delete the
            former lists or not.
        :arg str func:
            Optional, ``'max'`` by default. Not used.
        """

        if id_a not in self.lists:
            self._log('No such list: %s' % id_a, -5)
            return None

        if id_b not in self.lists:
            self._log('No such list: %s' % id_b, -5)
            return None

        name = '_'.join([id_a, id_b]) if name is None else name

        if isinstance(self.lists[id_a], list) and isinstance(
                self.lists[id_b], list):

            if and_or == "and":
                self.lists[name] = list(
                    set(self.lists[id_a]) | set(self.lists[id_b]))

            if and_or == "or":
                self.lists[name] = list(
                    set(self.lists[id_a]) & set(self.lists[id_b]))

        if isinstance(self.lists[id_a], dict) and isinstance(
                self.lists[id_b], dict):
            self.lists[name] = {}

            if and_or == "and":
                keys = list(
                    set(self.lists[id_a].keys) | set(self.lists[id_b].keys(
                    )))

                for k in keys:

                    if k in self.lists[id_a]:
                        self.lists[name][k] = self.lists[id_a][k]

                    if k in self.lists[id_b]:
                        self.lists[name][k] = self.lists[id_b][k]

                    if k in self.lists[id_a] and k in self.lists[id_b]:
                        self.lists[name][k] = self.combine_attr(
                            [self.lists[id_a][k], self.lists[id_b][k]])

            if and_or == "or":
                keys = list(
                    set(self.lists[id_a].keys) & set(self.lists[id_b].keys(
                    )))

                for k in keys:
                    self.lists[name][k] = self.combine_attr(
                        [self.lists[id_a][k], self.lists[id_b][k]])

        if delete:
            del self.lists[id_a]
            del self.lists[id_b]

    def save_session(self):
        """
        Save the current session state into pickle dump. The file will
        be saved in the current working directory as
        ``'pypath-<session_id>.pickle'``.
        """

        pickle_file = (
            'pypath-%s.pickle' % self.session_mod.get_session().label
        )
        self._log("Saving session to %s... " % pickle_file)

        with open(pickle_file, "wb") as f:
            pickle.dump(self, f, -1)

    ###
    # functions for plotting // with custom typeface ;)
    ###

    #
    # functions to compare networks and pathways
    #

    def databases_similarity(self, index='simpson'):
        """
        Computes the similarity across databases according to a given
        index metric. Computes the similarity across the loaded
        resources (listed in :py:attr:`pypath.main.PyPath.sources` in
        terms of nodes and edges separately.

        :arg str index:
            Optional, ``'simpson'`` by default. The type of index metric
            to use to compute the similarity. Options are ``'simpson'``,
            ``'sorensen'`` and ``'jaccard'``.

        :return:
            (*dict*) -- Nested dictionaries (three levels). First-level
            keys are ``'nodes'`` and ``'edges'``, then second and third
            levels correspond to sources names which map to the
            similarity index between those sources [float].
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
        Computes the similarity index across the given *groups*.

        :arg dict groups:
            Contains the different group names [str] as keys and their
            corresponding elements [set].
        :arg str index:
            Optional, ``'simpson'`` by default. The type of index metric
            to use to compute the similarity. Options are ``'simpson'``,
            ``'sorensen'`` and ``'jaccard'``.

        :return:
            (*dict*) -- Dictionary of dictionaries containing the groups
            names [str] as keys (for both inner and outer dictionaries)
            and the index metric as inner value [float] between those
            groups.
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
            self._log('No such function: %s()' % index_func, -5)

    def sorensen_pathways(self, pwlist=None):
        """
        Computes the Sorensen's similarity index across nodes and edges
        for the given list of pathway sources (all loaded pathway
        sources by default).

        :arg list pwlist:
            Optional, ``None`` by default. The list of pathway sources
            to be compared.

        :return:
            (*dict*) -- Nested dictionaries (three levels). First-level
            keys are ``'nodes'`` and ``'edges'``, then second and third
            levels correspond to ``<source>__<patwhay>`` names which map
            to the similarity index between those pathways [float].
        """

        g = self.graph

        if pwlist is None:
            self.update_pathway_types()
            pwlist = self.pathway_types

        for p in pwlist:

            if p not in g.vs.attributes():
                self._log('No such vertex attribute: %s' % p, -5)

        edges = {} # Keys = <source>__<pathway>, values = lsit of edge IDs
        nodes = {} # Keys = <source>__<pathway>, values = lsit of node IDs

        for e in g.es:
            indA = e.source
            indB = e.target
            pwsA = []
            pwsB = []

            for p in pwlist: # p is '<source>_pathways'

                if g.vs[indA][p] is not None:
                    # pw is '<patwhay>' from in node[p]
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
        Writes a given table to a file.

        :arg dict tbl:
            Contains the data of the table. It is assumed that keys are
            the row names [str] and the values, well, values. Column
            names (if any) are defined with the key ``'header'``.
        :arg str outfile:
            File name where to save the table. The file will be saved
            under the object's :py:attr:`pypath.main.PyPath.outdir`
            (``'results'`` by default).
        :arg str sep:
            Optional, ``'\t'`` (tab) by default. Specifies the separator
            for the file.
        :arg int cut:
            Optional, ``None`` by default. Specifies the maximum number
            of characters for the row names.
        :arg bool colnames:
            Optional, ``True`` by default. Specifies whether to write
            the column names in the file or not.
        :arg bool rownames:
            Optional, ``True`` by default. Specifies whether to write
            the row names in the file or not.
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
        Searches a given collection of attributes in a given object. As
        soon as one item is found, returns ``True``, if none could be
        found then returns ``False``.

        :arg object obj:
            Object (dictionary-like) where to search for elements of
            *lst*.
        :arg dict lst:
            Keys are the attribute names [str] and values the collection
            of elements to be searched in such attribute [set].

        :return:
            (*bool*) -- ``True`` if *lst* is empty or any of its
            elements is found in *obj*. Returns only ``False`` if cannot
            find anything.
        """

        if len(lst) == 0:
            return True

        for a, v in iteritems(lst): # XXX: Why call it lst if it's dict?

            if ((isinstance(v, list) and len(set(obj[a]).intersection(v)) > 0)
                    or (not isinstance(v, list) and obj[a] == v)):
                return True

        return False

    def search_attr_and(self, obj, lst):
        """
        Searches a given collection of attributes in a given object.
        Only returns ``True``, if all elements of *lst* can be found in
        *obj*.

        :arg object obj:
            Object (dictionary-like) where to search for elements of
            *lst*.
        :arg dict lst:
            Keys are the attribute names [str] and values the collection
            of elements to be searched in such attribute [set].

        :return:
            (*bool*) -- ``True`` only if *lst* is empty or all of its
            elements are found in *obj*. Returns ``False`` otherwise (as
            soon as one element of *lst* is not found).
        """

        for a, v in iteritems(lst): # XXX: Why call it lst if it's dict?

            if ((isinstance(v, list) and len(set(obj[a]).intersection(v)) == 0)
                    or (not isinstance(v, list) and obj[a] != v)):
                return False

        return True

    def get_sub(self, crit, andor="or", graph=None):
        """
        Selects the nodes from *graph* (and edges to be removed)
        according to a set of user-defined attributes.

        :arg dict crit:
            Defines the critical attributes to generate the subnetwork.
            Keys are ``'edge'`` and ``'node'`` and values are [dict]
            containing the critical attribute names [str] and values
            are [set] containing those attributes of the nodes/edges
            that are to be kept.
        :arg str andor:
            Optional, ``'or'`` by default. Determines the search mode.
            See :py:meth:`pypath.main.PyPath.search_attr_or` and
            :py:meth:`pypath.main.PyPath.search_attr_and` for more
            details.
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The graph object where to
            extract the subnetwork. If none is passed, takes the current
            network (undirected) graph
            (:py:attr:`pypath.main.PyPath.graph`).

        :return:
            (*dict*) -- Keys are ``'nodes'`` and ``'edges'`` whose
            values are [lst] of elements (as indexes [int]). Nodes are
            those to be kept and edges to be removed on the extracted
            subnetwork.
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
        Returns the sequence of all edge indexes that are not in
        the argument *edges*.

        :arg set edges:
            Sequence of edge indices [int] that will not be returned.

        :return:
            (*list*) -- Contains all edge indices [int] of the
            undirected network except the ones on *edges* argument.
        """

        g = self.graph
        inv = []

        # XXX: This could be refactored with a list comprehension:
        # return [e.index for e in g.es if e.index not in set(edges)]
        # also, assuming edges index is always == range(g.ecount()):
        # return set(range(g.ecount())) - set(edges)

        for e in g.es:

            if e.index not in set(edges):
                inv.append(e.index)

        return inv

    def get_network(self, crit, andor="or", graph=None):
        """
        Retrieves a subnetwork according to a set of user-defined
        attributes. Basically applies
        :py:meth:`pypath.main.PyPath.get_sub` on a given *graph*.

        :arg dict crit:
            Defines the critical attributes to generate the subnetwork.
            Keys are ``'edge'`` and ``'node'`` and values are [dict]
            containing the critical attribute names [str] and values
            are [set] containing those attributes of the nodes/edges
            that are to be kept.
        :arg str andor:
            Optional, ``'or'`` by default. Determines the search mode.
            See :py:meth:`pypath.main.PyPath.search_attr_or` and
            :py:meth:`pypath.main.PyPath.search_attr_and` for more
            details.
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The graph object where to
            extract the subnetwork. If none is passed, takes the current
            network (undirected) graph
            (:py:attr:`pypath.main.PyPath.graph`).

        :return:
            (*igraph.Graph*) -- The subgraph obtained from filtering
            according to the attributes defined in *crit*.
        """

        g = self.graph if graph is None else graph
        sub = self.get_sub(crit, andor=andor, graph=g)
        new = g.copy()
        new.delete_edges(sub["edges"])

        return new.induced_subgraph(sub["nodes"])

    def separate(self):
        """
        Separates the undirected network according to the different
        sources. Basically applies
        :py:meth:`pypath.main.PyPath.get_network` for each resource.

        :return:
            (*dict*) -- Keys are resource names [str] whose values are
            the subnetwork [igraph.Graph] containing the elements of
            that source.
        """

        return dict([(s, self.get_network({'edge': {'sources': [s]},
                                           'node': {}}))
                     for s in self.sources])

    def separate_by_category(self):
        """
        Separates the undirected network according to resource
        categories. Possible categories are:

            * ``'m'``: PTM/enzyme-substrate resources.
            * ``'p'``: Pathway/activity flow resources.
            * ``'i'``: Undirected/PPI resources.
            * ``'r'``: Process description/reaction resources.
            * ``'t'``: Transcription resources.

        Works in the same way as :py:meth:`pypath.main.PyPath.separate`.

        :return:
            (*dict*) -- Keys are category names [str] whose values are
            the subnetwork [igraph.Graph] containing the elements of
            those resources corresponding to that category.
        """

        cats = dict(list(map(lambda c:
            (c, list(filter(lambda s:
                s in self.sources, map(lambda cs:
                    cs[0], filter(lambda cs:
                        cs[1] == c, iteritems(db_categories.categories)))))),
                             self.has_cats)))

        return dict([(c, self.get_network({'edge': {'sources': s},
                                           'node': {}}))
                     for c, s in iteritems(cats)])

    def update_pathway_types(self):
        """
        Updates the pathway types attribute
        (:py:attr:`pypath.main.PyPath.pathway_types`) according to the
        loaded resources of the undirected network.
        """

        g = self.graph
        pwTyp = []

        for i in g.vs.attributes():

            if i.find("_pathways") > -1:
                pwTyp.append(i)

        self.pathway_types = pwTyp

    def source_similarity(self, outfile=None):
        """
        Computes the Sorensen's similarity index across nodes and edges
        for all the sources available (already loaded in the network)
        and saves them into table files. Files are stored in
        :py:attr:`pypath.main.PyPath.outdir` (``'results'`` by default).
        See :py:meth:`pypath.main.PyPath.databases_similarity` for more
        information.

        :arg str outfile:
            Optional, ``None`` by default. Specifies the file name
            prefix (suffixes will be ``'-nodes'`` and ``'-edges'``). If
            none is specified, this will be
            ``'pwnet-<session_id>-sim-src'``.
        """

        if outfile is None:
            outfile = ''.join(["pwnet-", self.session, "-sim-src"])

        res = self.database_similarity(index='sorensen')
        self.write_table(res["nodes"], outfile + "-nodes")
        self.write_table(res["edges"], outfile + "-edges")

    def pathway_similarity(self, outfile=None):
        """
        Computes the Sorensen's similarity index across nodes and edges
        for all the available pathway sources (already loaded in the
        network) and saves them into table files. Files are stored in
        :py:attr:`pypath.main.PyPath.outdir` (``'results'`` by default).
        See :py:meth:`pypath.main.PyPath.sorensen_pathways` for more
        information..

        :arg str outfile:
            Optional, ``None`` by default. Specifies the file name
            prefix (suffixes will be ``'-nodes'`` and ``'-edges'``). If
            none is specified, this will be
            ``'pwnet-<session_id>-sim-pw'``.
        """

        if outfile is None:
            outfile = ''.join(["pwnet-", self.session, "-sim-pw"])

        res = self.sorensen_pathways()
        self.write_table(res["nodes"], outfile + "-nodes", cut=20)
        self.write_table(res["edges"], outfile + "-edges", cut=20)


    def update_sources(self):
        """
        Makes sure that the :py:attr:`pypath.main.PyPath.sources`
        attribute is an up to date [list] of all sources in the current
        network.
        """

        g = self.graph
        src = []

        for e in g.es:
            src += e["sources"]

        self.sources = list(set(src))
        self.sources2 = set(
            ev.resource
            for attr in self.graph.es['attrs']
            for ev in attr.evidences
        )
        self.update_cats()


    def update_cats(self):
        """
        Makes sure that the :py:attr:`pypath.main.PyPath.has_cats`
        attribute is an up to date [set] of all categories in the
        current network.
        """

        self.has_cats = {
            db_categories.catnames[catletter]
            for src in self.sources
            for catletter in db_categories.get_categories(src)
        }
        self.has_cats2 = {
            ev.resource.data_model_label
            for attr in self.graph.es['attrs']
            for ev in attr.evidences
        }


    def update_pathways(self):
        """
        Makes sure that the :py:attr:`pypath.main.PyPath.pathways`
        attribute is an up to date [dict] of all pathways and their
        sources in the current network.
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
        Checks the network for any existing unmapped node and removes
        it.
        """

        g = self.graph

        if "unmapped" in g.vs["name"]:
            g.delete_vertices(g.vs.find(name="unmapped").index)
            self.update_db_dict()
            self.update_vname()

    def genesymbol_labels(self, graph=None, remap_all=False):
        """
        Creats vertex attribute ``'label'`` and fills up with the
        corresponding GeneSymbols of all proteins where the GeneSymbol
        can be looked up based on the default name of the protein
        vertex (UniProt ID by default). If the attribute ``'label'`` has
        been already initialized, updates this attribute or recreates if
        *remap_all* is set to ``True``.

        :arg igraph.Graph graph:
            Optional, ``None`` by default. The network graph object
            where the GeneSymbol labels are to be set/updated. If none
            is passed, takes the current network undirected graph by
            default (:py:attr:`pypath.main.PyPath.graph`).
        :arg bool remap_all:
            Optional, ``False`` by default. Whether to map anew the
            GeneSymbol labels if those were already initialized.
        """

        self._log('Updating vertex labels.')

        # XXX: What's the purpose of this? I mean attribute _directed is not
        #      accessed in this function (?)
        self._already_has_directed()

        if graph is None and self.dgraph is not None:
            self.genesymbol_labels(graph=self.dgraph, remap_all=remap_all)

        g = self.graph if graph is None else graph
        default_name_types = settings.get('default_name_types')
        label_name_types = {
            'protein': 'genesymbol',
            'mirna': 'mir-mat-name',
        }

        if 'label' not in g.vs.attributes():
            remap_all = True

        labels = [
            (
                None
                    if remap_all or v['label'] == v['name'] else
                v['label']
            )
            for v in g.vs
        ]

        for v, l, i in zip(g.vs, labels, xrange(g.vcount())):

            if l is None:

                label = None

                if isinstance(v['name'], intera.Complex):

                    label = v['name'].genesymbol_str

                elif (
                    v['type'] in label_name_types and
                    v['type'] in default_name_types
                ):

                    label = mapping.map_name0(
                        v['name'],
                        default_name_types[v['type']],
                        label_name_types[v['type']],
                        ncbi_tax_id=v['ncbi_tax_id'],
                    )

                if label:

                    labels[i] = label

                else:
                    labels[i] = v['name']

        g.vs['label'] = labels

    def network_stats(self, outfile=None):
        """
        Calculates basic statistics for the whole network and each of
        sources (node and edge counts, average node degree, graph
        diameter, transitivity, adhesion and cohesion). Writes the
        results in a tab file. File is stored in
        :py:attr:`pypath.main.PyPath.outdir` (``'results'`` by default).

        :arg str outfile:
            Optional, ``None`` by default. Specifies the file name. If
            none is specified, this will be
            ``'pwnet-<session_id>-stats'``.
        """

        if outfile is None:
            outfile = '-'.join(["pwnet", self.session, "stats"])

        stats = {}
        stats['header'] = ["vnum", "enum", "deg_avg", "diam", "trans", "adh",
                           "coh"]

        for k in xrange(0, len(self.sources) + 1):
            s = "All" if k == len(self.sources) else self.sources[k]
            g = (self.graph if k == len(self.sources)
                 else self.get_network({"edge": {"sources": [s]}, "node": {}}))

            if g.vcount() > 0:
                stats[s] = [g.vcount(), g.ecount(),
                            sum(g.vs.degree()) / float(len(g.vs)), # XXX: g.vcount()?
                            g.diameter(), g.transitivity_undirected(),
                            g.adhesion(), g.cohesion()]

        self.write_table(stats, outfile)

    def degree_dists(self):
        """
        Computes the degree distribution for all the different network
        sources. This is, for each source, the subnetwork comprising all
        interactions coming from it is extracted and the degree
        distribution information is computed and saved into a file.
        A file is created for each resource under the name
        ``'pwnet-<session_id>-degdist-<resource>''``. Files are stored
        in :py:attr:`pypath.main.PyPath.outdir` (``'results'`` by
        default).
        """

        dds = {}

        for s in self.sources:
            g = self.get_network({"edge": {"sources": [s]}, "node": {}})

            if g.vcount() > 0:
                dds[s] = g.degree_distribution()

        for k, v in iteritems(dds):
            filename = os.path.join(self.outdir,
                                    ''.join(["pwnet-", self.session,
                                             "-degdist-", k]))
            bins = []
            vals = []

            for i in v.bins():
                bins.append(int(i[0]))
                vals.append(int(i[2]))

            out = ''.join([";".join(str(x) for x in bins),
                           "\n", ";".join(str(x) for x in vals), "\n"])
            f = codecs.open(filename, encoding='utf-8', mode='w')
            f.write(out)
            f.close()

    def intergroup_shortest_paths(self, groupA, groupB, random=False): # TODO
        """

        """

        self.update_sources()

        if groupA not in self.graph.vs.attributes():
            self._log('No such attribute: %s' % groupA, -5)
            return False

        if groupB not in self.graph.vs.attributes():
            self._log('No such attribute: %s' % groupB, -5)
            return False

        deg_pathlen = {}
        rat_pathlen = {}
        rand_pathlen = {}
        diam_pathlen = {}

        for k in xrange(0, len(self.sources) + 1):
            s = "All" if k == len(self.sources) else self.sources[k]
            outfile = '-'.join([s, groupA, groupB, "paths"])
            f = (self.graph if k == len(self.sources)
                 else self.get_network({"edge": {"sources": [s]}, "node": {}}))
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

            self.write_table({"paths": paths}, outfile, sep=";",
                             colnames=False, rownames=False)
            deg = f.vs.degree()
            mean_pathlen = sum(paths) / float(len(paths))
            deg_pathlen[s] = [mean_pathlen, sum(deg) / float(len(deg))]
            rat_pathlen[s] = [mean_pathlen,
                              f.vcount() / float(len(list(set(grA + grB))))]
            diam_pathlen[s] = [mean_pathlen, f.diameter()]

            if random:
                groupA_random = groupA + "_random"
                groupB_random = groupB + "_random"
                random_pathlen = []

                for i in xrange(0, 100):
                    f.vs[groupA_random] = copy_mod.copy(f.vs[groupA])
                    f.vs[groupB_random] = copy_mod.copy(f.vs[groupB])
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
                            pt = f.get_shortest_paths(v.index, grA,
                                                      output="epath")

                            for p in pt:
                                l = len(p)

                                if l > 0:
                                    paths.append(l)

                            if (v.index in grA):
                                paths.append(0)

                    if len(paths) > 0:
                        random_pathlen.append(sum(paths) / float(len(paths)))

                if len(random_pathlen) > 0:
                    rand_pathlen[s] = [mean_pathlen,
                                       sum(random_pathlen)
                                       / float(len(random_pathlen))]

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
        Updates the all the vertex attributes ``'sources'`` and
        ``'references'`` according to their related edges (on the
        undirected graph).
        """

        g = self.graph

        for attr in ['sources', 'references']:
            g.vs[attr] = [set([]) for _ in g.vs]

            for e in g.es:
                g.vs[e.source][attr].update(e[attr])
                g.vs[e.target][attr].update(e[attr])


    def set_categories(self):
        """
        Sets the category attribute on the network nodes and edges
        (``'cat'``) as well the edge attribute coercing the references
        by category (``'refs_by_cat'``). The possible categories are
        as follows:

            * ``'m'``: PTM/enzyme-substrate resources.
            * ``'p'``: Pathway/activity flow resources.
            * ``'i'``: Undirected/PPI resources.
            * ``'r'``: Process description/reaction resources.
            * ``'t'``: Transcription resources.
        """

        self.graph.vs['cat'] = [set([]) for _ in self.graph.vs]
        self.graph.es['cat'] = [set([]) for _ in self.graph.es]
        self.graph.es['refs_by_cat'] = [{} for _ in self.graph.es]

        for v in self.graph.vs:

            for s in v['sources']:

                vcats = db_categories.get_categories(s)

                for cat in vcats:

                    v['cat'].add(cat)

        for e in self.graph.es:

            for s in e['sources']:

                ecats = db_categories.get_categories(s)

                for cat in ecats:

                    e['cat'].add(cat)

                    if cat not in e['refs_by_cat']:
                        e['refs_by_cat'][cat] = set()

                    if s in e['refs_by_source']:

                        e['refs_by_cat'][cat].update(e['refs_by_source'][s])


    def basic_stats_intergroup(self, groupA, groupB, header=None): # TODO
        """

        """

        result = {}
        g = self.graph

        for k in xrange(0, len(self.sources) + 1):
            s = "All" if k == len(self.sources) else self.sources[k]
            f = (self.graph if k == len(self.sources)
                 else self.get_network({"edge": {"sources": set([s])},
                                        "node": {}}))
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

        # FIXME: Here raises AttributeError. From above, I understood that
        #        arguments groupA and groupB were names for self.lists but
        #        here seems they are vertex attributes??

            for e in f.es:
                src.append(len(e["sources"]))

                if f.vs[e.source][groupA] or f.vs[e.target][groupA]:
                    csrc.append(len(e["sources"]))

                if f.vs[e.source][groupB] or f.vs[e.target][groupB]:
                    dsrc.append(len(e["sources"]))

            snum = sum(src) / float(len(src))
            csnum = sum(csrc) / float(len(csrc))
            dsnum = sum(dsrc) / float(len(dsrc))
            result[s] = [s, str(vnum), str(enum), str(cancerg), str(drugt),
                         str(cpct), str(dpct), str(tdgr), str(cdgr), str(ddgr),
                         str(tbwn), str(cbwn), str(dbwn), str(snum),
                         str(csnum), str(dsnum)]

        outfile = '-'.join([groupA, groupB, "stats"])

        if header is None:
            self.write_table(result, outfile, colnames=False)

        else:
            result["header"] = header
            self.write_table(result, outfile, colnames=True)

    def sources_venn_data(self, fname=None, return_data=False):
        """
        Computes the overlap in number of interactions for all pairs of
        sources.

        :arg str fname:
            Optional, ``None`` by default. If provided, saves the
            results into a table file. File is stored in
            :py:attr:`pypath.main.PyPath.outdir` (``'results'`` by
            default).
        :arg bool return_data:
            Optional, ``False`` by default. Whether to return the
            results as a [list].

        :return:
            (*list*) -- Only if *return_data* is set to ``True``. List
            of lists containing the counts for each pair of resources.
            This is, for instance, number of interactions only in
            resource A, number of interactions only in resource B and
            number of common interactions between A and B.
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
        Counts the number of sources per interaction in the graph and
        saves them into a file named ``source_num``. File is stored in
        :py:attr:`pypath.main.PyPath.outdir` (``'results'`` by
        default).
        """

        srcnum = [len(e['sources']) for e in self.graph.es]

        self.write_table({"srcnum": srcnum}, "source_num", sep=";",
                         rownames=False, colnames=False)

    def degree_dist(self, prefix, g=None, group=None):
        """
        Computes the degree distribution over all nodes of the network.
        If *group* is provided, also across nodes of that group(s).

        :arg str prefix:
            Prefix for the file name(s).
        :arg igraph.Graph g:
            Optional, ``None`` by default. The network over which to
            compute the degree distribution. If none is passed, takes
            the undirected network of the current instance.
        :arg list group:
            Optional, ``None`` by default. Additional group(s) name(s)
            [str] of node attributes to subset the network and compute
            its degree distribution.
        """

        if g is None:
            g = self.graph

        deg = g.vs.degree()
        self.write_table({"deg": deg}, prefix + "-whole-degdist", sep=";",
                         rownames=False, colnames=False)

        if group is not None:

            if not isinstance(group, list):
                group = [group]

            if len(set(group) - set(self.graph.vs.attributes())) > 0:
                self._log('Missing vertex attribute!', -5)
                return False

            for gr in group:
                dgr = [deg[i] for i, v in enumerate(g.vs) if v[gr]]

                self.write_table({"deg": dgr}, prefix + "-" + gr + "-degdist",
                                 sep=";", rownames=False, colnames=False)

    def delete_by_source(self, source, vertexAttrsToDel=None,
                         edgeAttrsToDel=None):
        """
        Deletes nodes and edges from the network according to a provided
        source name. Optionally can also remove the given list of
        attributes from nodes and/or edges.

        :arg str source:
            Name of the source from which the nodes and edges have to be
            removed.
        :arg list vertexAttrsToDel:
            Optional, ``None`` by default. Contains the names [str] of
            the attributes to be removed from the nodes.
        :arg list edgeAttrsToDel:
            Optional, ``None`` by default. Contains the names [str] of
            the attributes to be removed from the edges.
        """

        self.update_vertex_sources()
        g = self.graph
        verticesToDel = []

        # XXX: Refactor to?:
        #      verticesToDel = [v.index for v in g.vs
        #                if len(v['sources']-set([source])) == 0]

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
        Generates a file containing a table with information about the
        network's edges. First column contains the source node ID,
        followed by the target's ID, third column contains the number of
        references for that interaction and finally the number of
        sources. Writes the results in a tab file.

        :arg str filename:
            Optional, ``None`` by default. Specifies the file name and
            path to save the table. If none is passed, file will be
            saved in :py:attr:`pypath.main.PyPath.outdir` (``'results'``
            by default) with the name ``'<session_id>-refs-hist'``.
        """

        g = self.graph
        tbl = []

        for e in g.es:
            tbl.append((g.vs[e.source]['name'], g.vs[e.target]['name'],
                        str(len(e['references'])), str(len(e['sources']))))

        if filename is None:
            filename = os.path.join(self.outdir, self.session + "-refs-hist")

        out = ''

        for i in tbl:
            out += "\t".join(list(i)) + "\n"

        outf = codecs.open(filename, encoding='utf-8', mode='w')
        outf.write(out[:-1])
        outf.close()

    def load_resources(self, lst=None, exclude=[], cache_files={},
                       reread=False, redownload=None, keep_raw = False):
        """
        Loads multiple resources, and cleans up after. Looks up ID
        types, and loads all ID conversion tables from UniProt if
        necessary. This is much faster than loading the ID conversion
        and the resources one by one.

        :arg dict lst:
            Optional, ``None`` by default. Specifies the data input
            formats for the different resources (keys) [str]. Values
            are :py:class:`pypath.input_formats.NetworkInput` instances
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
        :arg bool reread:
            Optional, ``False`` by default. Specifies whether to reread
            the data files from the cache or omit them (similar to
            *redownload*).
        :arg bool redownload:
            Optional, ``False`` by default. Specifies whether to
            re-download the data and ignore the cache.
        """

        if lst is None:
            lst = omnipath

        huge = dict(
            (k, v)
            for k, v in iteritems(lst)
            if (
                (
                    (hasattr(v, 'huge') and v.huge) or (
                        hasattr(v, 'networkinput') and
                        v.networkinput.huge
                    ) or (
                        hasattr(v, 'huge') and
                        v.huge
                    )
                ) and
                k not in exclude and
                v.name not in exclude and
                v.name not in cache_files
            )
        )
        nothuge = dict(
            (k, v)
            for k, v in iteritems(lst)
            if (
                (
                    (hasattr(v, 'huge') and not v.huge) or (
                        hasattr(v, 'networkinput') and
                        not v.networkinput.huge
                    ) or (
                        hasattr(v, 'huge') and
                        not v.huge
                    ) and
                    v.name in cache_files
                ) and
                k not in exclude and
                v.name not in exclude
            )
        )

        for _inputs in [huge, nothuge]:

            for this_input in _inputs.values():

                self.load_resource(
                    this_input,
                    clean = False,
                    cache_files = cache_files,
                    reread = reread,
                    redownload = redownload,
                    keep_raw = keep_raw,
                )

        self._log(
            'load_resources(): all resources have been loaded, '
            'current number of nodes: %u, edges: %u' % (
                self.graph.vcount(),
                self.graph.ecount(),
            )
        )

        self.clean_graph()
        self.update_sources()
        self.update_vertex_sources()
        self.update_pathways()
        self.update_pathway_types()

        self._log(
            'Network has been built with %u interactions '
            'between %u nodes from %u resources' % (
                self.graph.ecount(),
                self.graph.vcount(),
                len(self.sources),
            )
        )


    def load_resource(
            self,
            settings,
            clean = True,
            cache_files = {},
            reread = None,
            redownload = None,
            keep_raw = False,
        ):
        """
        Loads the data from a single resource and attaches it to the
        network

        :arg pypath.input_formats.NetworkInput settings:
            :py:class:`pypath.input_formats.NetworkInput` instance
            containing the detailed definition of the input format to
            the downloaded file.
        :arg bool clean:
            Optional, ``True`` by default. Whether to clean the graph
            after importing the data or not. See
            :py:meth:`pypath.main.PyPath.clean_graph` for more
            information.
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

        self._log('Loading network data from resource `%s`.' % settings.name)

        self._read_network_data(
            settings,
            cache_files = cache_files,
            reread = reread,
            redownload = redownload,
            keep_raw = keep_raw,
        )
        self._add_network()

        if clean:
            self.clean_graph()

        self.update_sources()
        self.update_vertex_sources()


    def load_negatives(self):
        """
        """

        for k, v in iteritems(data_formats.negative):

            self._log(
                'Loading resource of negative interactions: `%s`.' % v.name
            )
            self.apply_negative(v)


    def load_dorothea(self, levels = {'A', 'B'}, only_curated = False):
        """
        Adds TF-target interactions from TF regulons to the network.
        DoRothEA is a comprehensive resource of TF-target
        interactions combining multiple lines of evidences: literature
        curated databases, ChIP-Seq data, PWM based prediction using
        HOCOMOCO and JASPAR matrices and prediction from GTEx expression
        data by ARACNe.

        For details see https://github.com/saezlab/DoRothEA.

        :arg set levels:
            Optional, ``{'A', 'B'}`` by default. Confidence levels to be
            loaded (from A to E) [str].
        :arg bool only_curated:
            Optional, ``False`` by default. Whether to retrieve only the
            literature curated interactions or not.
        """

        settings = copy_mod.deepcopy(
            network_resources.transcription['dorothea']
        )
        settings.networkinput.input_args = {
            'levels': levels,
            'only_curated': only_curated,
        }

        self.load_resources({'dorothea': settings})


    load_tfregulons = load_dorothea

    # XXX: Wouldn't it be better if only printed the resources loaded in
    #      **the current network** rather than omnipath? (e.g. if the user
    #      just removed one or several of them because he doesn't like/want them
    #      and then this function still prints them, as if they were still
    #      part of the network can be confusing).
    @staticmethod
    def list_resources():
        """
        Prints the list of resources through the standard output.
        """

        sys.stdout.write(' > omnipath\n')

        for k, v in iteritems(omnipath):
            sys.stdout.write('\t:: %s (%s)\n' % (v.name, k))

        sys.stdout.write(' > good\n')

        for k, v in iteritems(good): # FIXME: global name 'good' is not defined
            sys.stdout.write('\t:: %s (%s)\n' % (v.name, k))

    def info(self, name):
        """
        Given the name of a resource, prints out the information about
        that source/database. You can check the list of available
        resource descriptions in
        :py:func:`ypath.descriptions.descriptions.keys`.

        :arg str name:
            The name of the resource from which to print the
            information.
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

    # XXX: all these having_* functions are actually not used anywhere

    def having_attr(self, attr, graph=None, index=True, edges=True):
        """
        Checks if edges or nodes of the network have a specific
        attribute and returns an iterator of the indices (or the
        edge/node instances) of edges/nodes having such attribute.

        :arg str attr:
            The name of the attribute to look for.
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The graph object where the
            edge/node attribute is to be searched. If none is passed,
            takes the undirected network of the current instance.
        :arg bool index:
            Optional, ``True`` by default. Whether to return the
            iterator of the indices or the node/edge instances.
        :arg bool edges:
            Optional, ``True`` by default. Whether to look for the
            attribute in the networks edges or nodes instead.

        :return:
            (*generator*) -- Generator object containing the edge/node
            indices (or instances) having the specified attribute.
        """

        graph = graph or self.graph
        es_or_vs = getattr(graph, 'es' if edges else 'vs')

        if attr in es_or_vs.attributes():

            for i in es_or_vs:

                if common.something(i[attr]):
                    yield i.index if index else i

    def having_eattr(self, attr, graph=None, index=True):
        """
        Checks if edges of the network have a specific attribute and
        returns an iterator of the indices (or the edge instances) of
        edges having such attribute.

        :arg str attr:
            The name of the attribute to look for.
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The graph object where the
            edge/node attribute is to be searched. If none is passed,
            takes the undirected network of the current instance.
        :arg bool index:
            Optional, ``True`` by default. Whether to return the
            iterator of the indices or the node/edge instances.

        :return:
            (*generator*) -- Generator object containing the edge
            indices (or instances) having the specified attribute.
        """

        return self.having_attr(attr, graph, index)

    def having_vattr(self, attr, graph=None, index=True):
        """
        Checks if nodes of the network have a specific attribute and
        returns an iterator of the indices (or the node instances) of
        nodes having such attribute.

        :arg str attr:
            The name of the attribute to look for.
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The graph object where the
            edge/node attribute is to be searched. If none is passed,
            takes the undirected network of the current instance.
        :arg bool index:
            Optional, ``True`` by default. Whether to return the
            iterator of the indices or the node/edge instances.

        :return:
            (*generator*) -- Generator object containing the node
            indices (or instances) having the specified attribute.
        """

        return self.having_attr(attr, graph, index, False)

    def having_ptm(self, index=True, graph=None):
        """
        Checks if edges of the network have the ``'ptm'`` attribute and
        returns an iterator of the indices (or the edge instances) of
        edges having such attribute.

        :arg bool index:
            Optional, ``True`` by default. Whether to return the
            iterator of the indices or the node/edge instances.
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The graph object where the
            edge/node attribute is to be searched. If none is passed,
            takes the undirected network of the current instance.

        :return:
            (*generator*) -- Generator object containing the edge
            indices (or instances) having the ``ptm''`` attribute.
        """

        return self.having_eattr('ptm', graph, index)

    def loop_edges(self, index=True, graph=None):
        """
        Returns an iterator of the indices (or the edge instances) of
        the edges which represent a loop (whose source and target node
        are the same).

        :arg bool index:
            Optional, ``True`` by default. Whether to return the
            iterator of the indices or the edge instances.
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The graph object where the
            edge loops are to be searched. If none is passed, takes the
            undirected network of the current instance.

        :return:
            (*generator*) -- Generator object containing the edge
            indices (or instances) containing loops.
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
        Looks for the first neighbours of a given node and returns a
        list of their UniProt IDs.

        :arg str node:
            The UniProt ID of the node of interest. Can also be the
            index of such node [int].
        :arg bool indices:
            Optional, ``False`` by default. Whether to return the
            neighbour nodes indices or their UniProt IDs.

        :return:
            (*list*) -- The list containing the first neighbours of the
            queried node.
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
        Looks for the (first and) second neighbours of a given node and
        returns a list of their UniProt IDs.

        :arg str node:
            The UniProt ID of the node of interest. Can also be the
            index of such node [int].
        :arg bool indices:
            Optional, ``False`` by default. Whether to return the
            neighbour nodes indices or their UniProt IDs.
        :arg bool wit_first:
            Optional, ``False`` by default. Whether to return also the
            first neighbours or not.

        :return:
            (*list*) -- The list containing the second neighbours of the
            queried node (including the first ones if specified).
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
        Looks for the first neighbours of all the nodes and creates an
        attribute (``'neighbours'``) on each one of them containing a
        list of their UniProt IDs.

        :arg bool indices:
            Optional, ``False`` by default. Whether to list the
            neighbour nodes indices or their UniProt IDs.
        """

        g = self.graph
        g.vs['neighbours'] = [[] for _ in xrange(g.vcount())]
        prg = Progress(total=g.vcount(), name="Searching neighbours",
                       interval=30)

        for v in g.vs:
            v['neighbours'] = self.first_neighbours(v.index, indices=indices)
            prg.step()

        prg.terminate()

    def jaccard_edges(self):
        """
        Computes the Jaccard similarity index between the sets of first
        neighbours of all node pairs. **NOTE:** this method can take a
        while to compute, e.g.: if the network has 10K nodes, the total
        number of possible pairs to compute is:

        .. math::
          \\binom{10^4}{2} = 49995000

        :return:
            (*list*) -- Large list of [tuple] elements containing the
            node pair names [str] and their corresponding first
            neighbours Jaccard index [float].
        """

        g = self.graph
        self.all_neighbours(indices=True)
        metaEdges = []
        prg = Progress(total=g.vcount(), name="Calculating Jaccard-indices",
                       interval=11)

        # XXX: could reduce to a single for loop using itertools.combinations()
        for v in xrange(0, g.vcount() - 1):

            for w in xrange(v + 1, g.vcount()):
                vv = g.vs[v]
                vw = g.vs[w]

        # XXX: Why not using the function from common.jaccard_index()?

                ja = (len(set(vv['neighbours']) & set(vw['neighbours'])) /
                      float(len(vv['neighbours']) + len(vw['neighbours'])))
                metaEdges.append((vv['name'], vw['name'], ja))

            prg.step()

        prg.terminate()

        return metaEdges

    # XXX: Should consider returning a generator instead of list in the
    #      function above... that list and pa.graph are overflowing 16GB
    #      of RAM + 8GB of swap (maybe Chrome and Atom also helping, but
    #      you get the point)

    def jaccard_meta(self, jedges, critical):
        """
        Creates a (undirected) graph from a list of edges filtering by
        their Jaccard index.

        :arg list jedges:
            List of [tuple] containing the edges node names [str] and
            their Jaccard index. Basically, the output of
            :py:meth:`pypath.main.PyPath.jaccard_edges`.
        :arg float critical:
            Specifies the threshold of the Jaccard index from above
            which an edge will be included in the graph.

        :return:
            (*igraph.Graph*) -- The Undirected graph instance containing
            only the edges whose Jaccard similarity index is above the
            threshold specified by *critical*.
        """

        edges = []

        for e in jedges:

            if e[2] > critical:
                edges.append((e[0], e[1]))

        return igraph.Graph.TupleList(edges)

    def apply_negative(self, settings):
        """
        Loads a negative interaction source (e.g.: Negatome) into the
        current network.

        :arg pypath.input_formats.NetworkInput settings:
            :py:class:`pypath.input_formats.NetworkInput` instance
            containing the detailed definition of the input format to
            the downloaded file. For instance
            :py:data:`pypath.data_formats.negative['negatome']`
        """

        g = self.graph

        if settings.name not in self.negatives:
            self.raw_data = None
            self._read_network_data(settings)
            self.negatives[settings.name] = self.raw_data

        neg = self.negatives[settings.name]
        prg = Progress(total=len(neg), name="Matching interactions",
                       interval=11)
        matches = 0

        g.es['negative'] = [set([]) if e['negative'] is None else e['negative']
                            for e in g.es]
        g.es['negative_refs'] = [set([]) if e['negative_refs'] is None
                                 else e['negative_refs'] for e in g.es]

        for n in neg:
            aexists = n["default_name_a"] in g.vs['name']
            bexists = n["default_name_b"] in g.vs['name']

            if aexists and bexists:
                edge = self.edge_exists(n["default_name_a"], n["default_name_b"])

                if isinstance(edge, int):
                    g.es[edge]['negative'].add(settings.name)
                    refs = set(
                        list(
                            map(lambda r: _refs.Reference(int(r)), n[
                                'attrs_edge']['references'])))
                    g.es[edge]['negative_refs'].update(refs)
                    matches += 1

            prg.step()

        prg.terminate()
        sys.stdout.write('\t%u matches found with negative set\n' % matches)

    def negative_report(self, lst=True, outFile=None):
        """
        Generates a report file with the negative interactions (assumed
        to be already loaded).

        :arg bool lst:
            Optional, ``True`` by default. Whether to retun a list of
            edges containing the edge instances which have negative
            references.
        :arg str outFile:
            Optional, ``None`` by default. The output file name/path. If
            none is passed, the default is
            ``'results/<session_id>-negatives'``

        :return:
            (*list*) -- If *lst* is set to ``True``, returns a [list]
            is returned with the :py:class:`igraph.Edge` instances that
            contain at least a negative reference.
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
                    out += '\t'.join([g.vs[e.source]['name'],
                                      g.vs[e.target]['name'],
                                      g.vs[e.source]['label'],
                                      g.vs[e.target]['label'],
                                      ';'.join(list(e['sources'])),
                                      ';'.join(map(lambda r: r.pmid,
                                                   e['references'])),
                                      ';'.join(e['negative']),
                                      ';'.join(e['negative_refs'])]) + '\n'

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
        Exports a tab file containing the PTM interaction information
        loaded in the network.

        :arg str outfile:
            Optional, ``None`` by default. The output file nama/path to
            store the PTM information. If none is provided, the default
            is ``'results/network-<session_id>.tab'``.

        :return:
            (*list*) -- Contains the edge indices [int] of all PTM
            interactions.
        """

        # XXX: Shouldn't we use another default name for the file? too generic
        #      plus the same one is used in the functions below
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

        header = ['UniProt_A', 'UniProt_B', 'GeneSymbol_A', 'GeneSymbol_B',
                  'Databases', 'PubMed_IDs', 'Stimulation', 'Inhibition',
                  'Substrate-isoform', 'Residue_number', 'Residue_letter',
                  'PTM_type']

        stripJson = re.compile(r'[\[\]{}\"]')
        # first row is header
        outl = [header]

        with codecs.open(outfile, encoding='utf-8', mode='w') as f:
            f.write('\t'.join(header) + '\n')
            prg = Progress(total=self.graph.ecount(), name='Writing table',
                           interval=31)
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
                    row.append(';'.join([r for rs in
                                         [refs for db, refs
                                          in iteritems(e['refs_by_source'])
                                          if db in dbs] for r in rs]))
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
                                    r = row + ['%s-%u' % (dmi.ptm.protein,
                                                          dmi.ptm.isoform),
                                               str(dmi.ptm.residue.number),
                                               dmi.ptm.residue.name,
                                               dmi.ptm.typ]
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
        Exports a tab file containing the domain interaction information
        and PTM regulation loaded in the network.

        :arg str outfile:
            Optional, ``None`` by default. The output file nama/path to
            store the PTM information. If none is provided, the default
            is ``'results/network-<session_id>.tab'``.

        :return:
            (*list*) -- Contains the edge indices [int] of all PTM
            interactions.
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

        header = ['UniProt_A', 'UniProt_B', 'GeneSymbol_B', 'GeneSymbol_A',
                  'Databases', 'PubMed_IDs', 'Stimulation', 'Inhibition',
                  'Domain-domain', 'Domain-motif-PTM', 'PTM-regulation']

        stripJson = re.compile(r'[\[\]{}\"]')
        # first row is header
        outl = [header]

        with codecs.open(outfile, encoding='utf-8', mode='w') as f:
            f.write('\t'.join(header) + '\n')
            prg = Progress(total=self.graph.ecount(), name='Writing table',
                           interval=31)

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
                    row.append(';'.join([r for rs in
                                         [refs for db, refs in
                                          iteritems(e['refs_by_source'])
                                          if db in dbs] for r in rs]))
                    # signs
                    row += [str(int(x)) for x in e['dirs'].get_sign(di)]
                    # domain-domain
                    row.append('#'.join([x.serialize() for x in e['ddi']]))
                    # domain-motif
                    row.append('#'.join([x.serialize() for x in e['ptm']
                                         if x.__class__.__name__ == 'Ptm']))
                    # domain-motif
                    row.append('#'.join([x.serialize() for x in e['ptm'] if
                                         x.__class__.__name__ == 'Regulation']))
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
        Exports the network in a tabular format. By default UniProt IDs,
        Gene Symbols, source databases, literature references,
        directionality and sign information and interaction type are
        included.

        :arg str outfile:
            Optional, ``None`` by default. Name/path of the output file.
            If none is passed, ``'results/netrowk-<session_id>.tab'``
            is used.
        :arg dict extra_node_attrs:
            Optional, ``{}`` by default. Additional node attributes to
            be included in the exported table. Keys are column names
            used in the header while values are names of vertex
            attributes. In the header ``'_A'`` and ``'_B'`` suffixes
            will be appended to the column names so the values can be
            assigned to A and B side interaction partners.
        :arg dict extra_edge_attrs:
            Optional, ``{}`` by default. Additional edge attributes to
            be included in the exported table. Keys are column names
            used in the header while values are names of edge
            attributes.
        :arg bool unique_pairs:
            Optional, ``True`` by default. If set to ``True`` each line
            corresponds to a unique pair of molecules, all
            directionality and sign information are covered in other
            columns. If ``False``, order of ``'A'`` and ``'B'`` IDs
            corresponds to the direction while sign covered in further
            columns.
        :arg \*\* kwargs:
            Additional keyword arguments passed to
            :py:class:`pypath.export.Export`.
        """

        e = export.Export(pa=self, extra_node_attrs=extra_node_attrs,
                          extra_edge_attrs=extra_edge_attrs, **kwargs)
        e.write_tab(unique_pairs=unique_pairs, outfile=outfile)

    def export_sif(self, outfile=None):
        """
        Exports the network interactions in ``.sif`` format (Simple
        Interaction Format).

        :arg str outfile:
            Optional, ``None`` by default. Name/path of the output file.
            If none is passed, ``'results/netrowk-<session_id>.sif'``
            is used.
        """

        outfile = (outfile if outfile is not None
                   else 'network-%s.sif' % self.session)

        with open(outfile, 'w') as f:

            for e in self.graph.es:

                for d in [d for d, b in iteritems(e['dirs'].dirs) if b]:

                    if e['dirs'].is_directed() and d == 'undirected':
                        continue

                    sign = ('' if d == 'undirected' else
                            ''.join([['+', '-'][i] for i, v in
                                     enumerate(e['dirs'].get_sign(d)) if v]))
                    dirn = '=' if d == 'undirected' else '>'
                    source = (self.graph.vs[e.source]['name']
                              if d == 'undirected' else d[0])
                    target = (self.graph.vs[e.target]['name']
                              if d == 'undirected' else d[1])
                    f.write('\t'.join([source, sign + dirn, target]) + '\n')

    # TODO: Further description when function works/know what to do.
    def export_graphml(self, outfile=None, graph=None, name='main'):
        """
        Saves the network in a ``.graphml`` file.

        :arg str outfile:
            Optional, ``None`` by default. Name/path of the output file.
            If none is passed,
            ``'results/netrowk-<session_id>.graphml'`` is used.
        :arg igraph.Graph graph:
            Optional, ``None`` by default. The network object to be
            saved. If none is passed, takes the undirected network of
            the current instance.
        :arg str name:
            Optional, ``'main'`` by default. The graph name for the
            output file.
        """

        self.genesymbol_labels()
        g = self.graph if graph is None else graph

        if outfile is None:
            outfile = os.path.join(self.outdir,
                                   'network-' + self.session + '.graphml')

        is_directed = 'directed' if g.is_directed() else 'undirected'
        is_directedB = 'true' if g.is_directed() else 'false'
        node_attrs = [('UniProt', 'string'), ('GeneSymbol', 'string'),
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

            for attr in node_attrs:
                f.write(
                    '\t<key id="%s" for="node" attr.name="%s" attr.type="%s" />\n'
                    % (attr[0], attr[0], attr[1]))

            for attr in edgeAttrs:
                f.write(
                    '\t<key id="%s" for="edge" attr.name="%s" attr.type="%s" />\n'
                    % (attr[0], attr[0], attr[1]))

            f.write("""\n<graph id="%s" edgedefault="%s"
                        parse.nodeids="free" parse.edgeids="canonical"
                        parse.order="nodesfirst">\n\n""" % (name, is_directed))
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
                       g.vs[e.source]['name'], g.vs[e.target]['name'], is_directedB))
                f.write('\t<data key="Databases">%s</data>\n' %
                        (';'.join(list(e['sources']))))
                f.write(
                    '\t<data key="PubMedIDs">%s</data>\n' %
                    (';'.join(list(map(lambda r: r.pmid, e['references'])))))
    # XXX: attribute 'dirs_by_source' does not exist (nor created anywhere)
                f.write(
                    '\t<data key="Undirected">%s</data>\n' % (
                        ';'.join(sorted(
                            e['dirs'].sources['undirected']
                        ))
                    )
                )

                f.write(
                    '\t<data key="DirectionAB">%s</data>\n' % (
                        ';'.join(sorted(
                            e['dirs'].sources[e['dirs'].straight]
                        ))
                    )
                )
                f.write('\t<data key="DirectionAB">%s</data>\n' %
                        (';'.join(common.unique_list(e['dirs_by_source'][1]))))
                f.write('\t<data key="DirectionBA">%s</data>\n' %
                        (';'.join(common.unique_list(e['dirs_by_source'][2]))))
                f.write('\t<data key="StimulatoryAB">%s</data>\n' %
                        (';'.join(common.unique_list(e['signs'][0][0]))))
                f.write('\t<data key="InhibitoryAB">%s</data>\n' %
                        (';'.join(common.unique_list(e['signs'][0][1]))))
                f.write('\t<data key="StimulatoryBA">%s</data>\n' %
                        (';'.join(common.unique_list(e['signs'][1][0]))))
                f.write('\t<data key="InhibitoryBA">%s</data>\n' %
                        (';'.join(common.unique_list(e['signs'][1][1]))))
                f.write('\t<data key="InhibitoryBA">%s</data>\n' % (e['type']))
                f.write('</edge>\n')
                prg.step()

            f.write('\n\t</graph>\n</graphml>')

        prg.terminate()

    # XXX: Consider to remove kwarg `multi_query` since it's not used.
    def compounds_from_chembl(self, chembl_mysql=None, nodes=None, crit=None,
                              andor="or", assay_types=['B', 'F'],
                              relationship_types=['D', 'H'], multi_query=False,
                              **kwargs):
        """
        Loads compound data from ChEMBL to the network.

        :arg tuple chebl_mysql:
            Optional, ``None`` by default. Contains the MySQL parameters
            used by the :py:mod:`pypath.mapping` module to loadthe
            ChEMBL ID conversion tables. If none is passed, takes the
            current instance :py:attr:`pypath.main.PyPath.chembl_mysql`
            attribute.
        :arg list nodes:
            Optional, ``None`` by default. List of node indices for
            which the information is to be loaded. If none is provided
            calls the method :py:meth:`pypath.main.PyPath.get_sub` with
            the provided *crit* parameter.
        :arg dict crit:
            Optional, ``None`` by default. Defines the critical
            attributes to generate a subnetwork to extract the nodes in
            case *nodes* is not provided. Keys are ``'edge'`` and
            ``'node'`` and values are [dict] containing the critical
            attribute names [str] and values are [set] containing those
            attributes of the nodes/edges that are to be kept in the
            subnetwork. If none is provided, takes the whole network.
        :arg str andor:
            Optional, ``'or'`` by default. Determines the search mode
            for the subnetwork generation (if ``nodes=None``). See
            :py:meth:`pypath.main.PyPath.search_attr_or` and
            :py:meth:`pypath.main.PyPath.search_attr_and` for more
            details.
        :arg list assay_types:
            Optional, ``['B', 'F']`` by default. Types of assay to query
            Options are: ``'A'`` (ADME), ``'B'`` (Binding),``'F'``
            (Functional), ``'P'`` (Physicochemical), ``'T'`` (Toxicity)
            and/or ``'U'`` (Unassigned).
        :arg list relationship_types:
            Optional, ``['D', 'H']`` by default. Assay relationship
            types to query. Possible values are: ``'D'`` (Direct protein
            target assigned), ``'H'`` (Homologous protein target
            assigned), ``'M'`` (Molecular target other than protein
            assigned), ``'N'`` (Non-molecular target assigned), ``'S'``
            (Subcellular target assigned) and/or ``'U'`` (Default value,
            target has yet to be curated).
        :arg bool multi_query:
            Optional, ``False`` by default. Not used.
        :arg \*\*kwargs:
            Additional keyword arguments for
            :py:meth:`pypath.chembl.Chembl.compounds_targets`.
        """

        chembl_mysql = chembl_mysql or self.chembl_mysql
        self.chembl = chembl.Chembl(
            chembl_mysql, self.ncbi_tax_id)

        if nodes is None:

            if crit is None:
                nodes = xrange(0, self.graph.vcount())

            else:
                sub = self.get_sub(crit=crit, andor=andor)
                nodes = sub['nodes']

        uniprots = []

        for v in self.graph.vs:

            if v.index in nodes and v['id_type'] == 'uniprot':
                uniprots.append(v['name'])

        self.chembl.compounds_targets(uniprots, assay_types=assay_types,
                                      relationship_types=relationship_types,
                                      **kwargs)
        self.chembl.compounds_by_target()
        self.update_vname()
        self.graph.vs['compounds_chembl'] = [[] for _ in
                                             xrange(self.graph.vcount())]
        self.graph.vs['compounds_names'] = [[] for _ in
                                            xrange(self.graph.vcount())]
        self.graph.vs['compounds_data'] = [[] for _ in
                                           xrange(self.graph.vcount())]
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

                node['compounds_chembl'] = common.unique_list(node['compounds_chembl'])
                node['compounds_names'] = common.unique_list(node['compounds_names'])

        prg.terminate()
        percent = hascomp / float(self.graph.vcount())
        sys.stdout.write(
            '\n\tCompounds found for %u targets, (%.2f%% of all proteins).\n\n'
            % (hascomp, percent * 100.0))

    # XXX: Is it me or this function does not actually filter but only
    #      computes the score for the edges and that's all?
    def network_filter(self, p=2.0): # TODO
        """
        This function aims to cut the number of edges in the network,
        without losing nodes, to make the network less connected,
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
        Computes the distribution of shortest paths for each pair of
        nodes in the network (or between group(s) of nodes if *subset*
        is provided). **NOTE:** this method can take a while to compute,
        e.g.: if the network has 10K nodes, the total number of possible
        pairs to compute is:

        .. math::
          \\binom{10^4}{2} = 49995000

        :arg igraph.Graph graph:
            Optional, ``None`` by default. The network object for which
            the shortest path distribution is to be computed. If none is
            passed, takes the undirected network of the current
            instance.
        :arg tuple susbet:
            Optional, ``None`` by default. Contains two lists of node
            indices defining two groups between which the distribution
            is to be computed. Can also be [list] if the shortest paths
            are to be searched whithin the group. If none is passed, the
            whole network is taken by default.
        :arg str outfile:
            Optional, ``None`` by default. File name/path to save the
            shortest path distribution. If none is passed, no file is
            generated.
        :arg \*\*kwargs:
            Additional keyword arguments passed to
            :py:meth:`igraph.Graph.get_shortest_paths`.

        :return:
            (*list*) -- The length of the shortest paths for each pair
            of nodes of the network (or whithin/between group/s if
            *subset* is provided).
        """

        graph = graph if graph is not None else self.graph
        shortest_paths = []
        subset = (subset if isinstance(subset, tuple) or subset is None
                  else (subset, subset))
        prg = Progress(graph.vcount(), 'Calculating paths', 1)

        for i in xrange(0, graph.vcount() - 1): # i = node index

            if subset is None or i in subset[0] or i in subset[1]:
                paths = graph.get_shortest_paths(i,
                                                 xrange(i + 1, graph.vcount()),
                                                 **kwargs)

                for j in xrange(0, len(paths)): # j = `paths` list index
    # XXX: Or I'm missing something or something's wrong here... you're
    #      adding a node index to a index of the `paths` list and checking
    #      if such number is in a subset?
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
        Loads the 3D structure information from PDB into the network.
        Creates the node attribute ``'pdb'`` containing a [dict] whose
        keys are the PDB identifier [str] and values are [tuple] of two
        elements denoting the experimental method [str] (e.g.:
        ``'X-ray'``, ``'NMR'``, ...) and the resolution [float] (if
        applicable).

        :arg igraph.Graph graph:
            Optional, ``None`` by default. The network object for which
            the information is to be loaded. If none is passed, takes
            the undirected network of the current instance.
        """

        graph = graph if graph is not None else self.graph
        u_pdb, pdb_u = pdb_input.pdb_uniprot()

        if u_pdb is None:
            self._log('Failed to download UniProt-PDB dictionary')

        else:
            graph.vs['pdb'] = [None]

            for v in graph.vs:
                v['pdb'] = {}

                if v['name'] in u_pdb:

                    for pdb in u_pdb[v['name']]:
                        v['pdb'][pdb[0]] = (pdb[1], pdb[2])

            self._log('PDB IDs for proteins has been retrieved.')

    def load_pfam(self, graph=None):
        """
        Loads the protein family information from UniProt into the
        network. Creates the node attribute ``'pfam'`` containing a
        [list] of protein family identifier(s) [str].

        :arg igraph.Graph graph:
            Optional, ``None`` by default. The network object for which
            the information is to be loaded. If none is passed, takes
            the undirected network of the current instance.
        """

        graph = graph if graph is not None else self.graph
        u_pfam, pfam_u = pfam_input.pfam_uniprot(graph.vs['name'])

        if u_pfam is None:
            self._log('Failed to download Pfam data from UniProt')

        else:
        # FIXME: should assign `[None]` if the attribute is empty not otherwise
            graph.vs['pfam'] = [None]

            for v in graph.vs:
                v['pfam'] = []

                if v['name'] in u_pfam:
                    v['pfam'] += u_pfam[v['name']]

            self._log('Pfam domains has been retrieved.')

    def load_pfam2(self):
        """
        Loads the protein family information from Pfam into the network.
        Creates the node attribute ``'pfam'`` containing a [list] of
        [dict] whose keys are protein family identifier(s) [str] and
        corresponding values are [list] of [dict] containing detailed
        information about the protein family(ies) for regions and
        isoforms of the protein.

        :arg igraph.Graph graph:
            Optional, ``None`` by default. The network object for which
            the information is to be loaded. If none is passed, takes
            the undirected network of the current instance.
        """

        self.pfam_regions()

        if self.u_pfam is None:
            self._log('Failed to download data from Pfam', -5)

        else:
            self.graph.vs['pfam'] = [{} for _ in self.graph.vs]

            for v in self.graph.vs:

                if v['name'] in self.u_pfam:
                    v['pfam'] = self.u_pfam[v['name']]

            self._log('Pfam domains has been retrieved.')

    def load_pfam3(self):
        """
        Loads the protein domain information from Pfam into the network.
        Creates the node attribute ``'doms'`` containing a [list] of
        :py:class:`pypath.intera.Domain` instances with information
        about each domain of the protein (see the corresponding class
        documentation for more information).
        """

        self.pfam_regions()

        if self.u_pfam is None:
            self._log('Failed to download data from Pfam', -5)

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

            self._log('Pfam domains has been retrieved.')

    def load_corum(self, graph=None):
        """
        Loads complexes from CORUM database. Loads data into vertex attribute
        `graph.vs['complexes']['corum']`.
        This resource is human only.
        """

        graph = graph if graph is not None else self.graph
        complexes, members = corum.corum_complexes()

        if complexes is None:
            self._log('Failed to download data from CORUM', -5)

        else:
            self.init_complex_attr(graph, 'corum')

            for u, cs in iteritems(members):
                sw = mapping.map_name(u, 'uniprot', 'uniprot', 9606)

                for s in sw:

                    if s in graph.vs['name']:

                        for c in cs:
                            others = []

                            for memb in complexes[c[0]][0]:
                                others += mapping.map_name(memb,
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

            self._log('Complexes from CORUM have been retrieved.')

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
        complexes = havugimana.havugimana_complexes()

        if complexes is None:

            self._log('Failed to read data from Havugimana')

        else:
            self.init_complex_attr(graph, 'havugimana')

            for c in complexes:
                membs = []
                names = []

                for memb in c:
                    membs += mapping.map_name(memb,
                                                  'uniprot',
                                                  'uniprot',
                                                  9606)

                for u in membs:
                    names += mapping.map_name(u,
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

                self._log('Complexes from Havugimana have been retrieved.')


    def load_compleat(self, graph=None):
        """
        Loads complexes from Compleat. Loads data into vertex attribute
        `graph.vs['complexes']['compleat']`.
        This resource is human only.
        """

        graph = graph if graph is not None else self.graph
        complexes = compleat.compleat_complexes()

        if complexes is None:
            self._log('Failed to load data from COMPLEAT', -5)

        else:
            self.init_complex_attr(graph, 'compleat')

            for c in complexes:
                c['uniprots'] = []
                c['gsymbols'] = []

                for e in c['entrez']:
                    c['uniprots'] += mapping.map_name(e,
                                                          'entrez',
                                                          'uniprot',
                                                          9606)

                for u in c['uniprots']:
                    c['gsymbols'] += mapping.map_name(u,
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

            self._log('Complexes from COMPLEAT have been retrieved.')

    def load_complexportal(self, graph=None):
        """
        Loads complexes from ComplexPortal. Loads data into vertex attribute
        `graph.vs['complexes']['complexportal']`.
        This resource is human only.
        """

        graph = graph if graph is not None else self.graph
        # TODO: handling species
        complexes = complexportal.complexportal_complexes()

        if complexes is None:
            self._log('Failed to read data from ComplexPortal', -5)

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
                    swprots += mapping.map_name(u,
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

            self._log('Complexes from Complex Portal have been retrieved.')

    def load_3dcomplexes(self, graph=None):
        """
        """

        graph = graph if graph is not None else self.graph
        c3d = threedcomplex.threedcomplex_ddi()

        if c3d is None:
            self._log('Failed to download data from 3DComplex and PDB', -5)

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
                    up1s = mapping.map_name(ups[0], 'uniprot', 'uniprot')

                    for up1 in up1s:
                        inRefLists = False

                        for tax, lst in iteritems(self.reflists):

                            if up1 in lst.lst:
                                up2s = mapping.map_name(ups[1], 'uniprot',
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
                        name += mapping.map_name(sp,
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
        pisa, unmapped = pisa.pisa_interfaces(pdblist)

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
        comp, memb = corum.corum_complexes()

        for cname, cdata in comp.items():

            if search.lower() in cname.lower():
                result.append((cname, cdata[0]))

        return result


    @staticmethod
    def vertex_name(v):

        return v['name'] if isinstance(v, igraph.Vertex) else v


    def get_entity_type(self, entity):

        entity = self.vertex_name(entity)

        return entity_mod.Entity._get_entity_type(entity)


    def is_protein(self, entity):

        entity = self.vertex_name(entity)

        return entity_mod.Entity._is_protein(entity)


    def is_complex(self, entity):

        entity = self.vertex_name(entity)

        return entity_mod.Entity._is_complex(entity)


    def is_mirna(self, entity):

        entity = self.vertex_name(entity)

        return entity_mod.Entity._is_mirna(entity)


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

        if v_source is not None and v_target is not None:
            eid = self.graph.get_eid(
                v_source.index,
                v_target.index,
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
            (
                self.graph.vs[vi] for vi in
                itertools.chain(
                    *self.graph.neighborhood(
                        vs, order = order, mode = mode
                    )
                )
            ),
            self.nodNam,
            self.nodLab,
        )

    def up_neighborhood(self, uniprots, order = 1, mode = 'ALL'):
        """
        """

        if type(uniprots) in _const.SIMPLE_TYPES:
            uniprots = [uniprots]

        vs = self.uniprots(uniprots)
        return self._neighborhood(vs, order=order, mode=mode)

    def gs_neighborhood(self, genesymbols, order=1, mode='ALL'):
        """
        """

        if type(genesymbols) in _const.SIMPLE_TYPES:
            genesymbols = [genesymbols]

        vs = self.genesymbols(genesymbols)
        return self._neighborhood(vs, order=order, mode=mode)

    def neighborhood(self, identifiers, order=1, mode='ALL'):
        """
        """

        if type(identifiers) in _const.SIMPLE_TYPES:
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

    def edges_in_complexes(self, csources=['corum'], graph=None):
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

        #elif fun0 == __name__.split('.')[0]:
        #    toCall = __main__

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

    def edges_3d(
        self,
        methods = [
            'inputs.instruct.get_instruct',
            'inputs.i3d.get_i3d',
        ]
    ):
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

            self.u_pfam = pfam_input.pfam_regions(
                uniprots = self.graph.vs['name'],
                value = 'uniprot',
                keepfile = True,
            )

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
        domi = domino.domino_enzsub()

        if domi is None:
            self._log(
                'Failed to load domain-motif interaction data from DOMINO',
                -5
            )
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

        dmi = inputs.threedid.process_3did_dmi()

        if dmi is None:
            self._log(
                'Failed to load domain-motif interaction data from 3DID',
                -5,
            )
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

            if ddi.__module__.split('.')[1] == 'inputs':

                self._log('Function %s() failed' % ddi, -5)

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

            if dmi.__module__.split('.')[1] == 'inputs':

                self._log('Function %s() failed' % dmi, -5)

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


    def load_ddis(
        self,
        methods = [
            'inputs.threedcomplex.threedcomplex_ddi',
            'inputs.domino.domino_ddi',
            'self.load_3did_ddi2',
        ]
    ):
        """
        """

        self.run_batch(methods, toCall=self.load_ddi)


    def load_dmis(
        self,
        methods = [
            'self.pfam_regions',
            'self.load_depod_dmi',
            'self.load_dbptm',
            'self.load_mimp_dmi',
            'self.load_pnetworks_dmi',
            'self.load_domino_dmi',
            'self.load_pepcyber',
            'self.load_psite_reg',
            'self.load_psite_phos',
            'self.load_ielm',
            'self.load_phosphoelm',
            'self.load_elm',
            'self.load_3did_dmi',
        ]
    ):
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
        ddi, interfaces = threedid.threedid_ddi(residues=True)

        if ddi is None or interfaces is None:
            self._log(
                'Failed to load domain-domain interaction data from 3DID',
                -5,
            )
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
            swprots = mapping.swissprots([uniprot1, uniprot2])

            for swprot1 in swprots[uniprot1]:

                for swprot2 in swprots[uniprot2]:

                    if swprot1 in g.vs['name'] and swprot2 in g.vs['name']:
                        e = self.edge_exists(swprot1, swprot2)

                        if isinstance(e, int):

                            for domains, structures in iteritems(v):

                                for pdb, pdb_uniprot_pairs in \
                                        iteritems(structures['pdbs']):

                                    for pdbuniprots in pdb_uniprot_pairs:
                                        pdbswprots = mapping.swissprots(
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
        ddi, intfs = inputs.threedid.get_3did()

        if ddi is None or interfaces is None:
            self._log(
                'Failed to load domain-domain interaction data from 3DID',
                -5,
            )
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

        ielm = ielm_input.get_ielm(ppi)
        elm = elm_input.elm_classes()
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

        elm = elm_input.elm_interactions()
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
        toCall = inputs.get_method(functions[source])
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
                domain_ups = mapping.map_name(d['domain_protein'],
                                                  'uniprot', 'uniprot')
                motif_ups = mapping.map_name(d['motif_protein'], 'uniprot',
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

        pepc = pepcyber.pepcyber_interactions()
        prg = Progress(len(pepc), 'Processing domain-motif interactions', 13)

        for l in pepc:
            prg.step()
            uniprot1 = [l[9]] if len(l[9]) > 0 else []

            if len(l[9]) == 0 and len(l[10]) > 0:
                uniprot1 = mapping.map_name(l[10],
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

        preg = phosphosite.phosphosite_regsites()
        prg = Progress(len(preg), 'Processing regulatory effects', 11)

        for src, tgts in iteritems(preg):
            prg.step()

            # ptm on src
            # tgt: interactor of src, depending on ptm
            if self.node_exists(src):

                for tgt in tgts:

                    for effect in ['induces', 'disrupts']:

                        for ind in tgt[effect]:
                            uniprots = mapping.map_name(ind, 'genesymbol',
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

        comppi = comppi_input.comppi_interaction_locations()
        prg = Progress(len(comppi), 'Processing localizations', 33)

        for c in comppi:
            prg.step()
            uniprots1 = mapping.map_name(c['uniprot1'], 'uniprot',
                                             'uniprot', 9606)
            uniprots2 = mapping.map_name(c['uniprot2'], 'uniprot',
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

    def load_ptms2(
            self,
            input_methods = None,
            map_by_orthology_from = [9606],
            orthology_only_swissprot = True,
            ptm_orthology_strict = False,
            nonhuman_direct_lookup = True,
            inputargs = {},
            database = None,
            force_load = False,
        ):
        """
        This is a new method which will replace `load_ptms`.
        It uses `pypath.enz_sub.EnzymeSubstrateAggregator`, a newly
        introduced module for combining enzyme-substrate data from multiple
        resources using orthology translation on users demand.

        :param list input_methods: Resources to collect enzyme-substrate
            interactions from. E.g. `['Signor', 'phosphoELM']`. By default
            it contains Signor, PhosphoSitePlus, HPRD, phosphoELM, dbPTM,
            PhosphoNetworks, Li2012 and MIMP.
        :param list map_by_orthology_from: List of NCBI Taxonomy IDs of
            source taxons used for orthology translation of enzyme-substrate
            interactions. If you have a human network and you add here
            `[10090, 10116]` then mouse and rat interactions from the source
            databases will be translated to human.
        :param bool orthology_only_swissprot: `True` by default which means
            only SwissProt IDs are accepted at orthology translateion, Trembl
            IDs will be dropped.
        :param bool ptm_orthology_strict: For orthology translation use
            PhosphoSite's PTM orthology table. This guarantees that only
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
            `pypath.inputs` methods.
        :param database:
            A ``PtmAggregator`` object. If provided no new database will be
            created.
        :param bool force_load:
            If ``True`` the database will be loaded with the parameters
            provided here; otherwise if the ``ptm`` module already has a
            database no new database will be created. This means the
            parameters specified in other arguments might have no effect.
        """

        if database:

            ptma = database

        else:

            method = 'init_db' if force_load else 'get_db'

            _ = getattr(pypath.core.enz_sub, 'method')(
                input_methods = input_methods,
                ncbi_tax_id = self.ncbi_tax_id,
                map_by_orthology_from = map_by_orthology_from,
                orthology_only_swissprot = orthology_only_swissprot,
                ptm_orthology_strict = ptm_orthology_strict,
                nonhuman_direct_lookup = nonhuman_direct_lookup,
                inputargs = inputargs
            )

            ptma = pypath.core.enz_sub.get_db()

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

        if source == 'PhosphoSite' and self.ncbi_tax_id in taxonomy.taxids:
            kwargs['organism'] = taxonomy.taxids[self.ncbi_tax_id]

            if 'strict' not in kwargs:
                kwargs['strict'] = False

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
        toCall = inputs.get_method(functions[source])
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

            if (
                p['kinase'] is not None and
                not isinstance(p['kinase'], intera.Complex) and
                not isinstance(p['substrate'], intera.Complex) and
                len(p['kinase']) > 0
            ):

                # database specific id conversions
                if source in ['PhosphoSite', 'phosphoELM', 'Signor']:
                    kinase_ups = mapping.map_name(p['kinase'], 'uniprot',
                                                      'uniprot')

                else:

                    if not isinstance(p['kinase'], list):
                        p['kinase'] = [p['kinase']]
                    kinase_ups = [
                        i
                        for ii in [
                            mapping.map_name(k, 'genesymbol', 'uniprot')
                            for k in p['kinase']
                        ] for i in ii
                    ]

                    kinase_ups = list(set(kinase_ups))

                if not isinstance(p['kinase'], list):
                    p['kinase'] = [p['kinase']]

                if p['substrate'].startswith('HLA'):
                    continue

                if source in ['MIMP', 'PhosphoNetworks', 'Li2012']:
                    substrate_ups_all = mapping.map_name(
                        p['substrate'], 'genesymbol', 'uniprot')

                if source == 'MIMP':
                    substrate_ups_all.update(
                        mapping.map_name(
                            p['substrate_refseq'],
                            'refseq',
                            'uniprot',
                        )
                    )
                    substrate_ups_all = list(set(substrate_ups_all))

                if source in ['phosphoELM', 'dbPTM', 'PhosphoSite', 'Signor']:
                    substrate_ups_all = mapping.map_name(
                        p['substrate'], 'uniprot', 'uniprot')

                if source == 'HPRD':
                    substrate_ups_all = mapping.map_name(
                        p['substrate_refseqp'], 'refseqp', 'uniprot')
                    substrate_ups_all = common.unique_list(substrate_ups_all)

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
        data = depod_input.depod_enzyme_substrate()
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
                enzymes = mapping.map_name(l[4][0], 'uniprot', 'uniprot')
                substrates = mapping.map_name(l[5][0], 'uniprot',
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

        self._log(
            'Directionality set for %u interactions '
            'based on known (de)phosphorylation events.' % (isdir2 - isdir)
        )


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

        self._log(
            'Signes set based on phosphorylation-'
            'dephosphorylation pairs: %u' % new_signs
        )


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

                if (
                    p.__class__.__name__ == 'DomainMotif' and
                    p.ptm.typ == 'phosphorylation'
                ):

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
                network_resources.ptm.values(),
                network_resources.interaction_htp.values(),
                network_resources.pathway.values(),
                network_resources.transcription.values()
        ]:

            for db in dbs:
                result[db.name] = [bool(db.is_directed), bool(db.sign)]

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
                    row_order.append((db_categories.catnames[cat], 'subtitle'))
                    row_order.extend(
                        sorted(
                            filter(
                                lambda s:
                                    list(
                                        db_categories.get_categories(s)
                                    )[0] == cat,
                                self.sources
                            ),
                            key = lambda s: s.lower()
                        )
                    )

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

                for up in mapping.map_name(gsymb, 'genesymbol', 'uniprot'):

                    if up in self.nodDct:
                        data[up] = data[gsymb]

                del data[gsymb]

            self.exp = pandas.DataFrame(data)

        else:

            for gsymb, expr in iteritems(data):
                prg.step()
                uniprots = mapping.map_name(gsymb, 'genesymbol', 'uniprot')

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

    def edges_between(self, group1, group2, directed = True, strict = False):
        """
        Selects edges between two groups of vertex IDs.
        Returns set of edge IDs.

        :param set group1,group2:
            List, set or tuple of vertex IDs.
        :param bool directed:
            Only edges with direction `group1 -> group2` selected.
        :param bool strict:
            Edges with no direction information still selected even if
            ``directed`` is `False`.
        """

        edges = set()

        for e in self.graph.es:

            if (
                (e.source in group1 and e.target in group2) or
                (e.target in group1 and e.target in group2)
            ):

                if not directed or (not e['dirs'].is_directed and not strict):

                    edges.add(e.index)
                    continue

                up1 = self.up(e.source)
                up2 = self.up(e.target)

                if (
                    (e.source in group1 and e['dirs'].get_dir((up1, up2))) or
                    (e.target in group1 and e['dirs'].get_dir((up2, up1)))
                ):

                    edges.add(e.index)

        return edges

    def label(self, label, idx, what = 'vertices'):
        """
        Creates a boolean attribute ``label`` True for the
        vertex or edge IDs in the set ``idx``.
        """

        seq = self.graph.es if what == 'edges' else self.graph.vs
        cnt = self.graph.ecount() if what == 'edges' else self.graph.vcount()

        seq[label] = [False for _ in xrange(cnt)]

        for i in idx:

            seq[i][label] = True

    def label_vertices(self, label, vertices):
        """
        Creates a boolean vertex attribute ``label`` True for the
        vertex IDs in the set ``vertices``.
        """

        self.label(label, vertices)

    def label_edges(self, label, edges):
        """
        Creates a boolean edge attribute ``label`` True for the
        edge IDs in the set ``edges``.
        """

        self.label(label, edges, 'edges')

    def select_by_go_all(self, go_terms):
        """
        Selects the nodes annotated by all GO terms in ``go_terms``.

        Returns set of vertex IDs.

        :param list go_terms:
            List, set or tuple of GO terms.
        """

        annot = self.get_go()

        return annot.select_by_all(
            terms = go_terms,
            uniprots = self.graph.vs['name'],
        )

    def _select_by_go(self, go_terms):
        """
        Retrieves the vertex IDs of all vertices annotated with any of the
        Gene Ontology terms or their descendants.
        """

        annot = self.get_go()

        return annot.select_by_term(
            terms = go_terms,
            uniprots = self.graph.vs['name'],
        )

    def select_by_go(self, go_terms):
        """
        Retrieves the vertex IDs of all vertices annotated with any
        Gene Ontology terms or their descendants, or evaluates string
        expression (see ``select_by_go_expr``).

        :param str,set go_terms:
            A single GO term, a set of GO terms or an expression with
            GO terms.
        """

        annot = self.get_go()

        return annot.select(
            terms = go_terms,
            uniprots = self.graph.vs['name'],
        )

    def select_by_go_expr(self, go_expr):
        """
        Selects vertices based on an expression of Gene Ontology terms.
        Operator precedence not considered, please use parentheses.

        :param str go_expr:
            An expression of Gene Ontology terms. E.g.
            ``'(GO:0005576 and not GO:0070062) or GO:0005887'``. Parentheses
            and operators ``and``, ``or`` and ``not`` can be used.
        """

        annot = self.get_go()

        return annot.select_by_expr(
            expr = go_expr,
            uniprots = self.graph.vs['name'],
        )

    def label_by_go(self, label, go_terms, method = 'ANY'):
        """
        Assigns a boolean vertex attribute to nodes which tells whether
        the node is annotated by all or any of the GO terms.
        """

        _method = (
            self.select_by_go_all
                if method.upper() == 'ALL'
            else self.select_by_go
        )

        vids = _method(go_terms)

        self.label_vertices(label, vids)

    def network_by_go(
            self,
            node_categories,
            network_sources = None,
            include = None,
            exclude = None,
            directed = False,
            keep_undirected = False,
            prefix  = 'GO',
            delete  = True,
            copy = False,
            vertex_attrs = True,
            edge_attrs = True,
        ):
        """
        Creates or filters a network based on Gene Ontology annotations.

        :param dict node_categories:
            A dict with custom category labels as keys and expressions of
            GO terms as values. E.g.
            ``{'extracell': 'GO:0005576 and not GO:0070062',
               'plasmamem': 'GO:0005887'}``.
        :param dict network_sources:
            A dict with anything as keys and network input format definintions
            (``input_formats.NetworkInput`` instances) as values.
        :param list include:
            A list of tuples of category label pairs. By default we keep all
            edges connecting proteins annotated with any of the defined
            categories. If ``include`` is defined then only edges between
            category pairs defined here will be kept and all others deleted.
        :param list exclude:
            Similarly to include, all edges will be kept but the ones listed
            in ``exclude`` will be deleted.
        :param bool directed:
            If True ``include`` and ``exclude`` relations will be processed
            with directed (source, target) else direction won't be considered.
        :param bool keep_undirected:
            If True the interactions without direction information will be
            kept even if ``directed`` is True. Passed to ``edges_between``
            as ``strict`` argument.
        :param str prefix:
            Prefix for all vertex and edge attributes created in this
            operation. E.g. if you have a category label 'bar' and prefix
            is 'foo' then you will have a new vertex attribute 'foo__bar'.
        :param bool delete:
            Delete the vertices and edges which don't belong to any of the
            categories.
        :param bool copy:
            Return a copy of the entire ``PyPath`` object with the graph
            filtered by GO terms. By default the object is modified in place
            and ``None`` is returned.
        :param bool vertex_attrs:
            Create vertex attributes.
        :param bool edge_attrs:
            Create edge attributes.
        """

        if network_sources:

            self.init_network(network_sources)

        if self.graph.vcount() == 0:

            self.load_omnipath()

        if copy:

            graph_original = copy_mod.deepcopy(self.graph)

        vselections = dict(
            (
                label,
                self.select_by_go(definition)
            )
            for label, definition in iteritems(node_categories)
        )

        categories = tuple(vselections.keys())

        eselections = {}

        for label in itertools.product(categories, categories):

            if include and label not in include:

                continue

            elif exclude and label in exclude:

                continue

            eselections[label] = self.edges_between(
                    vselections[label[1]],
                    vselections[label[2]],
                    directed = directed,
                    strict = keep_undirected,
            )

        if vertex_attrs:

            for label, vertices in iteritems(vselections):

                self.label_vertices('%s__%s' % (prefix, label), vertices)

        if edge_attrs:

            for (label1, label2), edges in iteritems(eselections):

                self.label_edges(
                    '%s__%s__%s' % (prefix, label1, label2),
                    edges,
                )

        if delete:

            edges_to_delete = (
                set(xrange(self.graph.ecount())) -
                set.union(*eselections.values())
            )

            self.graph.delete_edges(edges_to_delete)

            self.graph.delete_vertices(
                np.where(np.array(self.graph.degree()) == 0)
            )

            self.update_vname()

        if copy:

            pypath_new = PyPath(copy = self)
            self.graph = graph_original
            self.update_vname()
            return pypath_new

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

            sources = network_resources.pathway

        CC_EXTRACELL    = 'GO:0005576' # select all extracellular
        CC_EXOSOME      = 'GO:0070062' # remove exosome localized proteins
        CC_PLASMAMEM    = 'GO:0005887' # select plasma membrane proteins
        MF_RECBINDING   = 'GO:0005102' # select ligands
        MF_RECACTIVITY  = 'GO:0038023' # select receptors
        MF_ECM_STRUCT   = 'GO:0005201' # select matrix structure proteins
                                       # e.g. collagene
        MF_CATALYTIC    = 'GO:0003824' # select enzymes, e.g. MMPs

        BP_SURF_REC_SIG = 'GO:0007166' # cell surface receptor signaling pw.
        CC_ECM_COMP     = 'GO:0044420' # ECM component


        if inference_from_go:

            go_desc = go_input.go_descendants_quickgo(aspects = ('C', 'F'))
            self.init_network(sources)

            if 'go' not in self.graph.vs.attributes():
                self.go_annotate()

            vids_extracell   = self.select_by_go(CC_EXTRACELL,   go_desc)
            vids_exosome     = self.select_by_go(CC_EXOSOME,     go_desc)
            #vids_extracell   = vids_extracell - vids_exosome
            vids_plasmamem   = self.select_by_go(CC_PLASMAMEM,   go_desc)
            vids_recbinding  = self.select_by_go(MF_RECBINDING,  go_desc)
            vids_recactivity = self.select_by_go(MF_RECACTIVITY, go_desc)

            receptors = vids_plasmamem & vids_recactivity
            ligands   = vids_extracell & vids_recbinding

            ureceptors = set(self.nodNam[i] for i in receptors)
            uligands = set(self.nodNam[i] for i in ligands)

            lig_with_interactions = set()
            rec_with_interactions = set()
            lig_rec_edges = set()
            lig_lig_edges = set()
            rec_rec_edges = set()

            # vertex attributes to mark ligands and receptors
            self.graph.vs['GO_ligand'] = [
                False for _ in xrange(self.graph.vcount())
            ]
            self.graph.vs['GO_receptor'] = [
                False for _ in xrange(self.graph.vcount())
            ]

            for e in self.graph.es:

                # set boolean vertex attributes
                if e.source in ligands:

                    self.graph.vs[e.source]['GO_ligand'] = True

                if e.target in ligands:

                    self.graph.vs[e.target]['GO_ligand'] = True

                if e.source in receptors:

                    self.graph.vs[e.source]['GO_receptors'] = True

                if e.target in receptors:

                    self.graph.vs[e.target]['GO_receptors'] = True

                # set edge attributes and collect edges to keep
                di = e['dirs']

                srcs = set(
                    self.up(u).index
                    for u in di.src(keep_undirected)
                )
                tgts = set(
                    self.up(u).index
                    for u in di.tgt(keep_undirected)
                )

                # ligand-receptor interaction
                if srcs & ligands and tgts & receptors:
                    lig_rec_edges.add(e.index)
                    lig_with_interactions.update(srcs)
                    rec_with_interactions.update(tgts)
                    e['sources'].add('GO_lig_rec')

                    if (
                        di.straight[0] in uligands and
                        di.straight[1] in ureceptors
                    ):

                        di.sources[di.straight].add('GO_lig_rec')

                    if (
                        di.reverse[0] in uligands and
                        di.reverse[1] in ureceptors
                    ):

                        di.sources[di.reverse].add('GO_lig_rec')

                # ligand-ligand interaction
                elif keep_lig_lig and srcs & ligands and tgts & ligands:
                    lig_lig_edges.add(e.index)
                    e['sources'].add('GO_lig_lig')

                    for d in di.sources.values():

                        if d:

                            d.add('GO_lig_lig')

                # receptor-receptor interaction
                elif keep_rec_rec and srcs & receptors and tgts & receptors:
                    rec_rec_edges.add(e.index)
                    e['sources'].add('GO_rec_rec')

                    for d in di.sources.values():

                        if d:

                            d.add('GO_rec_rec')

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

        self.update_vname()
        self.update_sources()

        if lig_rec_resources:

            datasets = copy_mod.deepcopy(network_resources.ligand_receptor)
            datasets['cellphonedb'].networkinput.input_args = {
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


    def go_annotate_graph(self, aspects = ('C', 'F', 'P')):
        """
        Annotates protein nodes with GO terms. In the ``go`` vertex
        attribute each node is annotated by a dict of sets where keys are
        one letter codes of GO aspects and values are sets of GO accessions.
        """

        go.annotate(self.graph, aspects = aspects)


    def load_go(self, organism = None):
        """
        Creates a ``pypath.go.GOAnnotation`` object for one organism in the
        dict under ``go`` attribute.

        :param int organism:
            NCBI Taxonomy ID of the organism.
        """

        organism = organism or self.ncbi_tax_id

        if not hasattr(self, 'go'):
            self.go = {}

        self.go[organism] = go.GOAnnotation(organism)


    def get_go(self, organism = None):
        """
        Returns the ``GOAnnotation`` object for the organism requested
        (or the default one).
        """

        organism = organism or self.ncbi_tax_id

        if not hasattr(self, 'go') or organism not in self.go:

            self.load_go(organism = organism)

        return self.go[organism]

    def go_enrichment(self, proteins=None, aspect='P', alpha=0.05,
                      correction_method='hommel', all_proteins=None):
        """
        Does not work at the moment because cfisher module should be
        replaced with scipy.
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
        Initializes a ``pypath.gsea.GSEA`` object and shows the list of the
        collections in MSigDB.
        """

        self.gsea = gsea.GSEA(user=user)
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
        Does not work at the moment because cfisher module should be
        replaced with scipy.
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
            set(graph.neighbors(node, mode = mode))
            for node in xrange(graph.vcount())
        ]

    def find_all_paths(
            self,
            start,
            end,
            attr = None,
            mode = 'OUT',
            maxlen = 2,
            graph = None,
            silent = False,
            update_adjlist = True,
        ):
        """
        Finds all paths up to length `maxlen` between groups of
        vertices. This function is needed only becaues igraph`s
        get_all_shortest_paths() finds only the shortest, not any
        path up to a defined length.

        start : int or list
            Indices of the starting node(s) of the paths.
        end : int or list
            Indices of the target node(s) of the paths.
        attr : str
            Name of the vertex attribute to identify the vertices by.
            Necessary if ``start`` and ``end`` are not igraph vertex ids
            but for example vertex names or labels.
        mode : 'IN', 'OUT', 'ALL'
            Passed to igraph.Graph.neighbors()
        maxlen : int
            Maximum length of paths in steps, i.e. if maxlen = 3, then
            the longest path may consist of 3 edges and 4 nodes.
        graph : igraph.Graph object
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

        if update_adjlist or not hasattr(self, 'adjlist'):

            self.update_adjlist(graph, mode = mode)

        all_paths = []

        start = start if isinstance(start, list) else [start]
        end = end if isinstance(end, list) else [end]

        if attr:

            attr_to_id = dict(reversed(i) for i in enumerate(graph.vs[attr]))
            start = [attr_to_id[a] for a in start]
            end = [attr_to_id[a] for a in end]

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

        if attr:

            all_paths = [
                [graph.vs[i][attr] for i in path]
                for path in all_paths
            ]

        return all_paths


    def find_all_paths2(
            self,
            graph,
            start,
            end,
            mode = 'OUT',
            maxlen = 2,
            psize = 100,
            update_adjlist = True,
        ):
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
                node, mode = mode
            ))
            for node in xrange(graph.vcount())
        ]
        paths = [[s] for s in start]
        all_paths = parts(paths, end, adjlist, maxlen, psize)
        sys.stdout.write('\n')

        return all_paths


    def transcription_factors(self):
        """
        """

        return common.unique_list([
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
            for uniprot, cnt in iteritems(
                collections.Counter(
                    itertools.chain(*(
                        prdb.expression[sample].keys()
                        for sample in samples
                    ))
                )
            )
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
        hpa = proteinatlas.get_proteinatlas(
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

        Args
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

    def set_plasma_membrane_proteins_cspa(self):
        """
        Creates a vertex attribute `cspa` with value *True* if
        the protein is a plasma membrane protein according to CPSA,
        otherwise *False*.
        """

        self.update_vname()
        self.graph.vs['cspa'] = [False for _ in self.graph.vs]

        if 'cspa' not in self.lists:
            self.cspa_list()

        for sp in self.lists['cspa']:

            if sp in self.nodDct:
                self.graph.vs[self.nodDct[sp]]['cspa'] = True

    def set_plasma_membrane_proteins_surfaceome(self, score_threshold = .0):
        """
        Creates a vertex attribute `ishs` with value *True* if
        the protein is a plasma membrane protein according to the In Silico
        Human Surfaceome, otherwise *False*.
        """

        self.update_vname()
        self.graph.vs['ishs'] = [False for _ in self.graph.vs]

        # always load this to apply score threshold
        self.surfaceome_list(score_threshold = score_threshold)

        for sp in self.lists['ishs']:

            if sp in self.nodDct:

                self.graph.vs[self.nodDct[sp]]['ishs'] = True

    def set_plasma_membrane_proteins_cspa_surfaceome(
            self,
            score_threshold = .0,
        ):
        """
        Creates a vertex attribute ``surf`` with value *True* if
        the protein is a plasma membrane protein according either to the
        Cell Surface Protein Atlas or the In Silico Human Surfaceome.
        """

        self.set_plasma_membrane_proteins_cspa()
        self.set_plasma_membrane_proteins_surfaceome(
            score_threshold = score_threshold
        )

        self.graph.vs['surf'] = [
            cspa or ishs
            for cspa, ishs in
            zip(self.graph.vs['cspa'], self.graph.vs['ishs'])
        ]

    def load_surfaceome_attrs(self):
        """
        Loads vertex attributes from the In Silico Human Surfaceome.
        Attributes are ``surfaceome_score``, ``surfaceome_class`` and
        ``surfaceome_subclass``.
        """

        attrs = (
            'surfaceome_score',
            'surfaceome_class',
            'surfaceome_subclass',
        )

        self.update_vname()

        sf = surfaceome_input.get_surfaceome()

        for i, attr in enumerate(attrs):

            self.graph.vs[attr] = [
                sf[v['name']][i]
                for v in self.graph.vs
                if v['name'] in sf
            ]

    def load_matrisome_attrs(self, organism = None):
        """
        Loads vertex attributes from MatrisomeDB 2.0. Attributes are
        ``matrisome_class``, ``matrisome_subclass`` and ``matrisome_notes``.
        """

        attrs = (
            'matrisome_class',
            'matrisome_subclass',
            'matrisome_notes',
        )

        organism = organism or self.ncbi_tax_id

        matrisome = matrisome_input.get_matrisome(organism = organism)

        for i, attr in enumerate(attrs):

            self.graph.vs[attr] = [
                matrisome[v['name']][i]
                for v in self.graph.vs
                if v['name'] in matrisome
            ]

    def load_membranome_attrs(self):
        """
        Loads attributes from Membranome, a database of single-helix
        transmembrane proteins.
        """

        self.update_vname()

        self.graph.vs['membranome_location'] = [
            None for _ in xrange(self.graph.vcount())
        ]

        m = membranome_input.get_membranome()

        for uniprot, mem, side in m:

            if uniprot in self.nodDct:

                self.graph.vs[
                        self.nodDct[uniprot]
                    ][
                        'membranome_location'
                    ] = (mem, side)

    def load_exocarta_attrs(
            self,
            load_samples = False,
            load_refs = False,
        ):
        """
        Creates vertex attributes from ExoCarta data. Creates a boolean
        attribute ``exocarts_exosomal`` which tells whether a protein is
        in ExoCarta i.e. has been found in exosomes. Optionally creates
        attributes ``exocarta_samples`` and ``exocarta_refs`` listing the
        sample tissue and the PubMed references, respectively.
        """

        self._load_exocarta_vesiclepedia_attrs(
            database = 'exocarta',
            load_samples = load_samples,
            load_refs = load_refs,
        )

    def load_vesiclepedia_attrs(
            self,
            load_samples = False,
            load_refs = False,
            load_vesicle_type = False,
        ):
        """
        Creates vertex attributes from Vesiclepedia data. Creates a boolean
        attribute ``vesiclepedia_in_vesicle`` which tells whether a protein is
        in ExoCarta i.e. has been found in exosomes. Optionally creates
        attributes ``vesiclepedia_samples``, ``vesiclepedia_refs`` and
        ``vesiclepedia_vesicles`` listing the sample tissue, the PubMed
        references and the vesicle types, respectively.
        """

        self._load_exocarta_vesiclepedia_attrs(
            database = 'vesiclepedia',
            load_samples = load_samples,
            load_refs = load_refs,
            load_vesicle_type = load_vesicle_type,
        )

    def _load_exocarta_vesiclepedia_attrs(
            self,
            database = 'exocarta',
            load_samples = False,
            load_refs = False,
            load_vesicle_type = False,
        ):

        database = database.lower()

        exo = exocarta_input._get_exocarta_vesiclepedia(
            database = database,
            organism = self.ncbi_tax_id
        )

        bool_attr = (
            'exocarta_exosomal'
                if database == 'exocarta' else
            'vesiclepedia_in_vesicle'
        )

        self.graph.vs[bool_attr] = [
            False for _ in xrange(self.graph.vcount())
        ]
        if load_samples:

            self.graph.vs['%s_samples' % database] = [
                set() for _ in xrange(self.graph.vcount())
            ]
        if load_refs:

            self.graph.vs['%s_refs' % database] = [
                set() for _ in xrange(self.graph.vcount())
            ]
        if database == 'vesiclepedia' and load_vesicle_type:

            self.graph.vs['%s_vesicles' % database] = [
                set() for _ in xrange(self.graph.vcount())
            ]

        for e in exo:

            uniprots = mapping.map_name(
                e[1],
                'genesymbol',
                'uniprot',
                ncbi_tax_id = self.ncbi_tax_id,
            )

            for u in uniprots:

                if u not in self.nodInd:

                    continue

                v = self.graph.vs[self.nodDct[u]]

                v[bool_attr] = True

                if load_samples:

                    v['%s_samples' % database].add(e[3][2])

                if load_refs and e[3][0] is not None:

                    v['%s_refs' % database].add(e[3][0])

                if database == 'vesiclepedia' and load_vesicle_type:

                    v['%s_vesicles' % database].update(set(e[3][3]))

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

        if 'sig' not in self.lists:
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

        Args
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

        Args
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

        Args
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

        if (
            isinstance(attrname, str) and
            inputs.get_method(attrname)
        ):

            fun = inputs.get_method(attrname)
            proteins_pws, interactions_pws = fun()

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

        Args
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
                        uniprots = mapping.map_name(protein, 'uniprot',
                                                        'uniprot')

                        for u in uniprots:

                            if u in nodDct:
                                g.vs[nodDct[u]][attrname].add(pw)

        if isinstance(interactions_pws, dict):

            for pw, ia in iteritems(interactions_pws):

                for pair in ia:
                    usrcs = mapping.map_name(pair[0], 'uniprot', 'uniprot')
                    utgts = mapping.map_name(pair[1], 'uniprot', 'uniprot')

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
        data = g2p_input.guide2pharma_download()

        for d in data:
            ulig = []
            urec = mapping.map_name(d['receptor_uniprot'], 'uniprot',
                                        'uniprot')

            if len(d['ligand_uniprot']) > 0:
                ulig = mapping.map_name(d['ligand_uniprot'], 'uniprot',
                                            'uniprot')

            if len(d['ligand_genesymbol']) > 0:
                ulig += mapping.map_name(d['ligand_genesymbol'],
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
        data = disgenet_input.disgenet_annotations(dataset=dataset)
        self.graph.vs['dis'] = [[] for _ in self.graph.vs]

        for d in data:

            if d['score'] >= score:
                uniprots = mapping.map_name(d['entrez'], 'entrez',
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
                common.flat_list([[r.pmid for r in e['references']]
                                 for e in self.graph.es])))

        cats = list(
            reduce(lambda e1, e2: e1 | e2['cat'], self.graph.es, set(
                []))) if by_category else []

        for s in list(self.sources) + cats:
            sattr = 'cat' if s in db_categories.catnames else 'sources'
            rattr = 'refs_by_cat' if s in db_categories.catnames else 'refs_by_source'

            cat = db_categories.get_category(s)

            catmembers = set(db_categories.catnames.keys()) \
                if s in db_categories.catnames \
                else set(self.sources) if not hasattr(db_categories, cat) \
                else getattr(db_categories, cat)

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
                common.unique_list([
                    r.pmid
                    for r in common.flat_list(
                        [e[rattr][s] for e in self.graph.es if s in e[rattr]])
                ]))

            other_refs = set(
                common.flat_list([[
                    r.pmid
                    for r in common.flat_list([
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
                        common.flat_list([[(r.pmid, e.index) for r in rr]
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
                        common.flat_list([[(r.pmid, e.index) for r in rr]
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

            if s in db_categories.catnames:
                s = db_categories.catnames[s]

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

            if key in db_categories.catnames:
                name = db_categories.catnames[key]

                if by_category:
                    cs[(name, 'subtitle')] = cs[name]

                else:

                    if name in cs:
                        del cs[name]

        row_order = []

        if by_category:

            for cat in use_cats:
                row_order.append((db_categories.catnames[cat], 'subtitle'))
                row_order.extend(
                    sorted(
                        filter(
                            lambda name:
                                cat in db_categories.get_categories(name),
                             cs.keys()
                        )
                    )
                )

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

    def load_omnipath(
            self,
            omnipath = None,
            kinase_substrate_extra = False,
            ligand_receptor_extra = False,
            pathway_extra = False,
            remove_htp = True,
            htp_threshold = 1,
            keep_directed = True,
            min_refs_undirected = 2,
            old_omnipath_resources = False,
            exclude = None,
            pickle_file = None,
        ):
        """
        Loads the OmniPath network.
        Note, if ``pickle_file`` provided the network will be loaded directly
        from there regardless of its content.
        """

        # XXX: According to the alias above omnipath = data_formats.omnipath
        # already
        # YYY: Ok, but here the user has a chance to override it, is it bad?

        if pickle_file and os.path.exists(pickle_file):

            self.init_network(pickle_file = pickle_file)
            return


        def reference_constraint(formats, extra, cat):
            """
            If we anyways load extra interactions without references it does
            not make sense to throw away the records without references from
            the default OmniPath sources.
            """

            if not extra:

                return formats

            formats_noref = {}

            for name, fmt in iteritems(formats):

                fmt_noref = copy_mod.deepcopy(fmt)

                if fmt.name in getattr(db_categories, cat):

                    fmt_noref.networkinput.must_have_references = False

                formats_noref[name] = fmt_noref

            return formats_noref


        exclude = exclude or []

        if omnipath is None:

            if old_omnipath_resources:
                omnipath = copy_mod.deepcopy(network_resources.omnipath)
                omnipath['biogrid'] = network_resources.interaction['biogrid']
                omnipath['alz'] = network_resources.interaction['alz']
                omnipath['netpath'] = network_resources.interaction['netpath']
                exclude.extend(['intact', 'hprd'])

            else:
                omnipath = network_resources.omnipath

        omnipath = reference_constraint(omnipath, pathway_extra, 'p')
        omnipath = reference_constraint(omnipath, kinase_substrate_extra, 'm')
        omnipath = reference_constraint(omnipath, ligand_receptor_extra, 'l')

        self.load_resources(omnipath, exclude = exclude)

        if kinase_substrate_extra:
            self.load_resources(network_resources.ptm_misc)

        if ligand_receptor_extra:
            self.load_resources(network_resources.ligand_receptor)

        if pathway_extra:
            self.load_resources(network_resources.pathway_noref)

        self.third_source_directions()

        if remove_htp:
            self.remove_htp(
                threshold = htp_threshold,
                keep_directed = keep_directed,
            )

        if not keep_directed:
            self.remove_undirected(min_refs = min_refs_undirected)


    def remove_htp(self, threshold = 50, keep_directed = False):
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
        self._log(
            'Interactions with only high-throughput references '
            'have been removed. %u interactions removed. '
            'Number of edges decreased from %u to %u, '
            'number of vertices from %u to %u.' % (
                len(htedgs),
                ecount_before,
                self.graph.ecount(),
                vcount_before,
                self.graph.vcount()
            )
        )


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
        self._log(
            'Undirected interactions %s have been removed. '
            '%u interactions removed. Number of edges '
            'decreased from %u to %u, number of vertices '
            'from %u to %u.' % (
                ''
                    if min_refs is None else
                'with less than %u references' % min_refs,
                len(udedgs),
                ecount_before,
                self.graph.ecount(),
                vcount_before,
                self.graph.vcount()
            )
        )


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
        refc = collections.Counter(
            common.flat_list((r.pmid for r in e['references'])
                            for e in self.graph.es))

        for htlim in reversed(xrange(1, 201)):
            htrefs = set([i[0] for i in refc.most_common() if i[1] > htlim])
            htedgs = [
                e.index for e in self.graph.es
                if len(set([r.pmid for r in e['references']]) - htrefs) == 0
            ]
            htsrcs = common.unique_list(
                common.flat_list([self.graph.es[e]['sources'] for e in htedgs]))
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

        keggd = kegg_input.kegg_interactions()
        self.process_directions(keggd, 'KEGG', stimulation='activation',
                                inhibition='inhibition', graph=graph)

    def phosphosite_directions(self, graph=None):
        """
        """

        psite = phosphosite_input.phosphosite_directions()
        self.process_directions(psite, 'PhosphoSite_dir', dirs_only=True,
                                id_type='uniprot', graph=graph)

    def phosphopoint_directions(self, graph=None):
        """
        """

        ppoint = phosphopoint_input.phosphopoint_directions()
        self.process_directions(ppoint, 'PhosphoPoint', dirs_only=True,
                                id_type='genesymbol', graph=graph)

    def phosphonetworks_directions(self, graph=None):
        """
        """

        pnet = phosphonetworks_input.phosphonetworks_interactions()
        self.process_directions(pnet, 'PhosphoNetworks', dirs_only=True,
                                id_type='genesymbol', graph=graph)

    def mimp_directions(self, graph=None):
        """
        """

        mimp = mimp_input.mimp_interactions()
        self.process_directions(mimp, 'MIMP', dirs_only=True,
                                id_type='genesymbol', graph=graph)

    def laudanna_directions(self, graph=None):
        """
        """

        laud = laudanna_input.laudanna_directions()
        self.process_directions(laud, 'Laudanna_sigflow', dirs_only=True,
                                id_type='genesymbol', graph=graph)

    def laudanna_effects(self, graph=None):
        """
        """

        laud = laudanna_input.laudanna_effects()
        self.process_directions(laud, 'Laudanna_effects',
                                stimulation='activation',
                                inhibition='inhibition', directed='docking',
                                id_type='genesymbol', graph=graph)

    def string_effects(self, graph=None):
        """
        """

        string = string_input.string_effects()
        self.process_directions(string, 'STRING', stimulation='+',
                                inhibition='-', directed='*', id_type='ensp',
                                graph=graph)

    def acsn_effects(self, graph=None):
        """
        """

        acsnd = acsn_input.get_acsn_effects()
        self.process_directions(acsnd, 'ACSN', stimulation='+', inhibition='-',
                                directed='*', id_type='genesymbol',
                                graph=graph)

    def wang_effects(self, graph=None):
        """
        """

        wangd = wang_input.wang_interactions()
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
                    else mapping.map_name(k[0], id_type, 'uniprot')
                tgt = [k[1]] if id_type is None \
                    else mapping.map_name(k[1], id_type, 'uniprot')

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

        self._log(
            'Directions and signs set for %u edges based on %s,'
            ' %u new directions, %u new signs.' % (
                sourcedirs,
                name,
                newdirs,
                newsigns
            )
        )


    def update_summaries(self):
        """
        Creates a dict with many summarizing and comparative statistics
        about the resources in the current network.
        The result will be assigned to the attribute ``summaries``.
        """

        def one_entity(stats, resource):

            key = stats.method
            label = stats.label or stats.method.replace('_', ' ').capitalize()

            result = []

            for ext, attr in (
                ('', 'by_resource'),
                ('unique', 'unique'),
                ('shared', 'shared'),
                (('percent', '[%]'), 'percent'),
            ):

                key_sep, lab_sep = ('_', ' ') if ext else ('', '')
                key_ext, lab_ext = (
                    ext if isinstance(ext, tuple) else (ext, ext)
                )

                result.append(
                    (
                        'n_%s%s%s' % (key, key_sep, key_ext),
                        '%s%s%s' % (label, lab_sep, lab_ext),
                        getattr(stats.counts, attr)[resource],
                    )
                )

            return result


        summaries = {}
        summaries_labels = {}

        resources = self.resources
        references = self.references_stats()
        entities = self.entities_stats()
        interactions_all = self.interactions_all_stats()
        interactions_undirected = self.interactions_undirected_stats()
        interactions_directed = self.interactions_directed_stats()
        interactions_signed = self.interactions_signed_stats()
        interactions_stimulatory = self.interactions_stimulatory_stats()
        interactions_inhibitory = self.interactions_inhibitory_stats()
        interactions_mutual = self.interactions_mutual_stats()
        curation_effort = self.curation_effort_stats()

        for resource in itertools.chain(resources, ('Total',)):

            summaries[resource] = [
                ('name', 'Resource', resource),
            ]

            for stats in (
                entities,
                interactions_all,
                interactions_undirected,
                interactions_directed,
                interactions_signed,
                interactions_stimulatory,
                interactions_inhibitory,
                interactions_mutual,
                references,
                curation_effort,
            ):

                summaries[resource].extend(
                    one_entity(
                        stats = stats,
                        resource = resource,
                    )
                )

            summaries_labels = (
                summaries_labels or
                collections.OrderedDict(
                    (key, label)
                    for key, label, value in summaries[resource]
                )
            )

            summaries[resource] = collections.OrderedDict(
                (key, value)
                for key, label, value in summaries[resource]
            )

        self.summaries = summaries
        self.summaries_labels = summaries_labels


    def mean_reference_per_interaction(self, resources = None):
        """
        Computes the mean number of references per interaction of the
        network.

        :return:
            (*float*) -- Mean number of interactions per edge.
        """

        return (
            self.numof_references(resources = resources) /
            self.numof_edges(resources = resources)
        )


    def mean_reference_per_interaction_by_resource(self, resources = None):
        """
        Computes the mean number of references per interaction of the
        network.

        :return:
            (*float*) -- Mean number of interactions per edge.
        """

        return self._by_resource(
            method = self.mean_reference_per_interaction,
            resources = resources,
        )


    def numof_edges(self, resources = None):
        """
        Number of edges optionally limited to certain resources.
        """

        return len(list(self.iter_edges(resources = resources)))


    def iter_edges(self, resources = None):
        """
        Iterates the edges in the graph optionally limited to certain
        resources. Yields ``igraph.Edge`` objects.
        """

        resources = common.to_set(resources)

        for e in self.graph.es:

            if not resources or resources & e['sources']:

                yield e


    def numof_reference_interaction_pairs(self): # XXX: Not really sure about this one
        """
        Returns the total of unique references per interaction.

        :return:
            (*int*) -- Total number of unique references per
            interaction.
        """

        return len(common.unique_list(common.flat_list(
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


    @classmethod
    def _remove_cp(cls, resources):
        """
        Removes the *CellPhoneDB* derived resources from a set in order
        to avoid them to be compared as they were primary resources.
        """
        return cls._remove_by_label(resources, labels = {'CP'})


    @staticmethod
    def _remove_by_label(elements, labels):

        return {
            elem
            for elem in elements
            if next(reversed(elem.split('_'))) not in labels
        }


    @property
    def resources(self):
        """
        All network resources. Returns *set* of strings.
        """

        return self._remove_cp(self._collect_across_edges('sources'))


    def _collect_across_edges(self, attr):
        """
        For a *set* type edge attribute ``attr`` collects the elements
        across all edges.
        """

        return set(itertools.chain(*self.graph.es[attr]))


    def _by_resource(self, method, resources = None, **kwargs):
        """
        Calls a method for each resource (by default for all resources).
        Returns dict with resources as keys and the output of the method
        as values.
        """

        resources = common.to_set(resources) or self.resources
        method = method if callable(method) else getattr(self, method)

        return dict(
            (
                resource,
                method(resources = resource, **kwargs)
            )
            for resource in resources
        )


    @staticmethod
    def resources_by_category():

        return dict(
            (
                label,
                {
                    resource
                    for resource, cat in iteritems(db_categories.categories)
                    if cat == letter
                }
            )
            for label, letter in iteritems(db_categories.catletters)
        )


    def _by_category(self, method, **kwargs):

        cat_res = self.resources_by_category()
        method = method if callable(method) else getattr(self, method)

        result = {}

        for cat in db_categories.catletters.keys():

            if not cat_res[cat] & self.resources:

                continue

            entities = method(resources = cat_res[cat], **kwargs)

            if entities:

                result[cat] = entities

        return result


    @property
    def all_references(self):

        return self._collect_across_edges('references')


    def numof_references(self, resources = None, **kwargs):
        """
        Counts the number of reference on the network.

        Counts the total number of unique references in the edges of the
        network.

        resources : None,str,set
            Limits the query to one or more resources.

        :return:
            (*int*) -- Number of unique references in the network.
        """

        return len(self.references(resources = resources))


    def references(self, resources = None, **kwargs):
        """
        Returns a set of references for all edges.

        resources : None,str,set
            Limits the query to one or more resources.
        """

        resources = common.to_set(resources)

        return set.union(*(
            refs
            for e in self.graph.es
            for res, refs in iteritems(e['refs_by_source'])
            if not resources or res in resources
        ))


    def numof_references(self, resources = None, **kwargs):

        return len(self.references(resources = resources))


    def references_by_resource(self, resources = None, **kwargs):
        """
        Creates a dict with resources as keys and sets of references
        as values.
        """

        return self._by_resource(self.references, resources = resources)


    def references_by_category(self, **kwargs):

        return self._by_category(self.references)


    def numof_references_by_resource(self, resources = None, **kwargs):
        """
        Counts the references for each resource, optionally limited
        to certain resources.
        """

        return self._by_resource(self.numof_references, resources = resources)


    def numof_references_by_category(self, **kwargs):

        return self._by_category(selsf.numof_references)


    def curation_effort(self, resources = None, **kwargs):
        """
        Returns a *set* of reference-interactions pairs.
        """

        resources = common.to_set(resources) or self.resources

        return {
            (ref, self.nodNam[e.source], self.nodNam[e.target])
            for e in self.graph.es
            for resource, refs in iteritems(e['refs_by_source'])
            for ref in refs
            if not resources or resource in resources
        }


    def curation_effort_by_resource(self, resources = None, **kwargs):
        """
        A *dict* with resources as keys and *set*s of curation items
        (interaction-reference pairs) as values.
        """

        return self._by_resource(self.curation_effort, resources = resources)


    def curation_effort_by_category(self, **kwargs):

        return self._by_category(self.curation_effort)


    def iter_vertices(
            self,
            resources = None,
            categories = None,
            entity_type = None,
        ):

            for _id in self.iter_entities(
                resources = resources,
                categories = categories,
                entity_type = entity_type,
            ):

                yield self.graph.vs(self.nodDct[_id])


    def iter_entities(
            self,
            resources = None,
            categories = None,
            entity_type = None,
        ):
        """
        Iterates nodes optionally only for certain resources. Yields
        ``igraph.Vertex`` objects.
        """

        resources = common.to_set(resources)
        categories = common.to_set(categories)

        for iattr in self.graph.es['attrs']:

            if (
                (
                    not resources or
                    iattr.evidences & resources
                ) and (
                    not categories or
                    iattr.data_models % categories
                )
            ):

                if (
                    not entity_type or
                    iattr.entity_type_a == entity_type
                ):

                    yield iattr_id_a

                if (
                    not entity_type or
                    iattr.entity_type_b == entity_type
                ):

                    yield iattr_id_b


    def iter_entities(self, resources = None, entity_type = None):

        for v in self.iter_vertices(
            resources = resources,
            entity_type = entity_type,
        ):

            yield v['name']


    def entities(self, resources = None, entity_type = None, **kwargs):

        return set(
            self.iter_entities(
                resources = resources,
                entity_type = entity_type,
            )
        )


    def entities_by_resource(
            self,
            resources = None,
            entity_type = None,
            **kwargs
        ):
        """
        Returns a *dict* of *set*s with resources as keys and sets of
        entities as values.
        """

        return self._by_resource(
            method = self.entities,
            resources = resources,
            entity_type = entity_type,
        )


    def entities_by_category(self, entity_type = None, **kwargs):

        return self._by_category(self.entities, entity_type = entity_type)


    def protein_entities(self, resources = None, **kwargs):

        return self.entities(resources = resources, entity_type = 'protein')


    def protein_entities_by_resource(self, resources = None, **kwargs):

        return self.entities_by_resource(
            resources = resources,
            entity_type = 'protein',
        )


    def protein_entities_by_category(self, resources = None, **kwargs):

        return self.entities_by_category(
            resources = resources,
            entity_type = 'protein',
        )


    def mirna_entities(self, resources = None, **kwargs):

        return self.entities(resources = resources, entity_type = 'mirna')


    def mirna_entities_by_resource(self, resources = None, **kwargs):

        return self.entities_by_resource(
            resources = resources,
            entity_type = 'mirna',
        )


    def mirna_entities_by_category(self, resources = None, **kwargs):

        return self.entities_by_category(
            resources = resources,
            entity_type = 'mirna',
        )


    def complex_entities(self, resources = None):

        return self.entities(resources = resources, entity_type = 'complex')


    def complex_entities_by_resource(self, resources = None, **kwargs):

        return self.entities_by_resource(
            resources = resources,
            entity_type = 'complex',
        )


    def complex_entities_by_category(self, resources = None, **kwargs):

        return self.entities_by_category(
            resources = resources,
            entity_type = 'complex',
        )

    #
    # interactions undirected
    #

    def interactions_undirected(self, resources = None, **kwargs):
        """
        Returns a *set* of tuples of node name pairs without being aware
        of the directions.
        Pairs of node names will be sorted alphabetically.
        """

        resources = common.to_set(resources)

        return {
            tuple(sorted(
                (
                    self.nodNam[e.source],
                    self.nodNam[e.target],
                )
            ))
            for e in self.graph.es
            if (
                not resources and
                not e['dirs'].is_directed()
            ) or (
                e['dirs'].sources['undirected'] & resources
            )
        }


    def interactions_undirected_by_resource(self, resources = None, **kwargs):
        """
        Returns a *dict* of *set*s of tuples of node name pairs without
        being aware of the directions.
        Pairs of node names will be sorted alphabetically.
        """

        return self._by_resource(
            method = self.interactions_undirected,
            resources = resources,
        )


    def interactions_undirected_by_category(self, **kwargs):

        return self._by_category(self.interactions_undirected)

    #
    # interactions directed
    #

    def interactions_directed(self, resources = None, **kwargs):
        """
        Returns a *set* of tuples of node name pairs with being aware
        of the directions.
        Undirected interactions will be discarded.
        Pairs of node names represent the directions: first is the source,
        second is the target.
        """

        return self._interactions_directed(resources = resources, **kwargs)


    def interactions_directed_by_resource(
            self,
            resources = None,
            effect = None,
            **kwargs
        ):
        """
        Returns a *dict* of *set*s of tuples of node name pairs with being
        aware of the directions.
        Undirected interactions will be discarded.
        Pairs of node names represent the directions: first is the source,
        second is the target.
        """

        return self._by_resource(
            method = self._interactions_directed,
            resources = resources,
            effect = effect,
            **kwargs
        )


    def interactions_directed_by_category(self, effect = None, **kwargs):

        return self._by_category(
            method = self._interactions_directed,
            effect = effect,
        )


    def _interactions_directed(
            self,
            resources = None,
            effect = None,
            **kwargs
        ):

        resources = common.to_set(resources)
        method = 'which_directions' if not effect else 'which_signs'
        args = {} if not effect else {'effect': effect}

        return set(
            itertools.chain(*(
                getattr(e['dirs'], method)(resources = resources, **args)
                for e in self.graph.es
            )
        ))


    #
    # interactions signed
    #

    def interactions_signed(self, resources = None, **kwargs):
        """
        Returns a *set* of tuples of node name pairs only for signed
        interactions.
        Pairs of node names represent the directions: first is the source,
        second is the target.
        """

        return self._interactions_directed(
            resources = resources,
            effect = True,
        )


    def interactions_signed_by_resource(self, resources = None, **kwargs):
        """
        Returns a *dict* of *set*s of tuples of node name pairs with being
        aware of the directions.
        Undirected interactions will be discarded.
        Pairs of node names represent the directions: first is the source,
        second is the target.
        """

        return self._by_resource(
            method = self._interactions_directed,
            resources = resources,
            effect = True,
        )


    def interactions_signed_by_category(self, **kwargs):

        return self.interactions_directed_by_category(effect = True)


    #
    # interactions stimulatory
    #

    def interactions_stimulatory(self, resources = None, **kwargs):
        """
        Returns a *set* of tuples of node name pairs only for stimulatory
        interactions.
        Pairs of node names represent the directions: first is the source,
        second is the target.
        """

        return self._interactions_directed(
            resources = resources,
            effect = 'stimulation',
        )


    def interactions_stimulatory_by_resource(
            self,
            resources = None,
            **kwargs
        ):
        """
        Returns a *dict* of *set*s of tuples of node name pairs with being
        aware of the directions.
        Undirected interactions will be discarded.
        Pairs of node names represent the directions: first is the source,
        second is the target.
        """

        return self._by_resource(
            method = self._interactions_directed,
            resources = resources,
            effect = 'stimulation',
        )


    def interactions_stimulatory_by_category(self, **kwargs):

        return self.interactions_directed_by_category(effect = 'stimulation')


    #
    # interactions inhibitory
    #

    def interactions_inhibitory(self, resources = None, **kwargs):
        """
        Returns a *set* of tuples of node name pairs only for inhibitory
        interactions.
        Pairs of node names represent the directions: first is the source,
        second is the target.
        """

        return self._interactions_directed(
            resources = resources,
            effect = 'inhibition',
        )


    def interactions_inhibitory_by_resource(self, resources = None, **kwargs):
        """
        Returns a *dict* of *set*s of tuples of node name pairs with being
        aware of the directions.
        Undirected interactions will be discarded.
        Pairs of node names represent the directions: first is the source,
        second is the target.
        """

        return self._by_resource(
            method = self._interactions_directed,
            resources = resources,
            effect = 'inhibition',
        )


    def interactions_inhibitory_by_category(self, **kwargs):

        return self.interactions_directed_by_category(effect = 'inhibition')

    #
    # interactions mutual
    #

    def interactions_mutual(self, resources = None, **kwargs):
        """
        Returns a *set* of tuples of node name pairs representing
        mutual interactions (i.e. A-->B and B-->A).
        Pairs of node names will be sorted alphabetically.
        """

        return {
            tuple(e['dirs'].nodes) for e in self.graph.es
            if e['dirs'].is_mutual(resources = resources)
        }


    def interactions_mutual_by_resource(self, resources = None, **kwargs):
        """
        Returns a *dict* of *set*s of tuples of node name pairs representing
        mutual interactions (i.e. A-->B and B-->A).
        Pairs of node names will be sorted alphabetically.
        """

        return self._by_resource(
            method = self.interactions_mutual,
            resources = resources,
        )


    def interactions_mutual_by_category(self, **kwargs):

        return self._by_category(self.interactions_mutual)


    def interactions_all(self, resources = None, **kwargs):
        """
        Returns a *set* of tuples of node name pairs representing
        interactions, both directed and undirected. Directed interactions
        will be present according to their direction, mutual directed
        interactions are represented by two tuples. The directed and
        undirected interactions are not distinguished here.
        """

        return (
            self.interactions_undirected(resources = resources, **kwargs) |
            self.interactions_directed(resources = resources, **kwargs)
        )


    def interactions_all_by_resource(self, resources = None, **kwargs):

        return self._by_resource(
            method = self.interactions_all,
            resources = resources,
        )


    def interactions_all_by_category(self, **kwargs):

        return self._by_category(method = self.interactions_all)


    #
    # methods for collecting and counting entities
    #

    def collect(self, method, **kwargs):
        """
        Collects various entities over the network according to ``method``.
        """

        def cat_shared_unique(method):

            return (
                dict(
                    itertools.chain(
                        *(
                            iteritems(
                                getattr(
                                    common,
                                    '%s_foreach' % method
                                )(
                                    dict(
                                        (
                                            resource,
                                            entities
                                        )
                                        for resource, entities
                                        in iteritems(by_resource)
                                        if resource in resources
                                    )
                                )
                            )
                            for cat, resources in iteritems(cat_resource)
                            if resources
                        )
                    )
                )
            )


        self._log('Collecting `%s`.' % method)

        total = getattr(self, method)(**kwargs)
        by_resource = getattr(self, '%s_by_resource' % method)(**kwargs)
        by_category = getattr(self, '%s_by_category' % method)(**kwargs)
        shared = common.shared_foreach(by_resource)
        unique = common.unique_foreach(by_resource)

        cat_resource = dict(
            (
                cat,
                {
                    resource
                    for resource in by_resource.keys()
                    if cat in {
                        db_categories.catnames[c]
                        for c in db_categories.get_categories(resource)
                    }
                }
            )
            for cat, resources in iteritems(by_category)
            if resources
        )
        resource_cat = common.swap_dict(cat_resource)

        shared_res_cat = cat_shared_unique(method = 'shared')
        unique_res_cat = cat_shared_unique(method = 'unique')

        shared_cat = common.shared_foreach(by_category)
        unique_cat = common.shared_foreach(by_category)

        return NetworkEntityCollection(
            total = total,
            by_resource = by_resource,
            by_category  = by_category,
            shared = shared,
            unique = unique,
            shared_res_cat = shared_res_cat,
            unique_res_cat = unique_res_cat,
            shared_cat = shared_cat,
            unique_cat = unique_cat,
            resource_cat = resource_cat,
            cat_resource = cat_resource,
            method = method,
        )


    def counts(
            self,
            collection_method,
            add_total = True,
            add_percent = True,
            add_cat_total = True,
            **kwargs
        ):
        """
        Collects various entities over the network according to ``method``
        and counts them in total and by resources.
        """

        coll = (
            collection_method
                if isinstance(
                    collection_method,
                    NetworkEntityCollection
                ) else
            self.collect(collection_method, **kwargs)
        )

        self._log('Counting `%s`.' % coll.method)

        n_total = len(coll.total)
        n_by_resource = common.dict_counts(coll.by_resource)
        n_by_category = common.dict_counts(coll.by_category)
        n_shared = common.dict_counts(coll.shared)
        n_unique = common.dict_counts(coll.unique)
        _percent = (
            common.dict_percent(n_by_resource, n_total)
                if add_percent else
            None
        )
        n_shared_res_cat = common.dict_counts(coll.shared_res_cat)
        n_unique_res_cat = common.dict_counts(coll.unique_res_cat)

        for resource in coll.by_resource.keys():

            if not db_categories.get_category(resource):

                self._log('Category not known for resource `%s`.' % resource)

        percent_res_cat = dict(
            (
                resource,
                (
                    n_by_resource[resource] /
                    n_by_category[
                        db_categories.catnames[
                            db_categories.get_category(resource)
                        ]
                    ] * 100
                        if n_by_resource[resource] else
                    .0
                )
            )
            for resource in coll.by_resource.keys()
        )
        n_shared_cat = common.dict_counts(coll.shared_cat)
        n_unique_cat = common.dict_counts(coll.shared_cat)
        percent_cat = (
            common.dict_percent(n_by_category, n_total)
                if add_percent else
            None
        )

        if add_total:

            if _percent:

                _percent['Total'] = 100.

            n_by_resource['Total'] = n_total
            n_shared['Total'] = common.n_shared_total(coll.by_resource)
            n_unique['Total'] = common.n_unique_total(coll.by_resource)
            n_shared_res_cat['Total'] = (
                common.n_shared_total(coll.by_category)
            )
            n_unique_res_cat['Total'] = (
                common.n_unique_total(coll.by_category)
            )
            percent_res_cat['Total'] = 100.

        if add_cat_total:

            for cat in n_by_category.keys():

                n_by_resource[cat] = n_by_category[cat]
                n_shared[cat] = n_shared_cat[cat]
                n_unique[cat] = n_unique_cat[cat]
                _percent[cat] = percent_cat[cat]

                this_cat_by_resource = dict(
                    it
                    for it in iteritems(coll.by_resource)
                    if it[0] in coll.cat_resource[cat]
                )

                n_shared_res_cat[cat] = common.n_shared_total(
                    this_cat_by_resource
                )
                n_unique_res_cat[cat] = common.n_unique_total(
                    this_cat_by_resource
                )

                coll.resource_cat[cat] = cat
                coll.cat_resource[cat].add(cat)

        return NetworkStatsRecord(
            total = n_total,
            by_resource = n_by_resource,
            by_category = n_by_category,
            shared = n_shared,
            unique = n_unique,
            percent = _percent,
            shared_res_cat = n_shared_res_cat,
            unique_res_cat = n_unique_res_cat,
            percent_res_cat = percent_res_cat,
            shared_cat = n_shared_cat,
            unique_cat = n_unique_cat,
            percent_cat = percent_cat,
            resource_cat = coll.resource_cat,
            cat_resource = coll.cat_resource,
            method = coll.method,
        )


    def stats(self, method, keep_collection = False, **kwargs):
        """
        Creates a collection of entities over the network according to
        ``method`` and counts them. By default the collection won't be
        returned but only the counts.
        """

        NetworkEntities = collections.namedtuple(
            'NetworkEntities',
            [
                'counts',
                'entities',
                'method',
                'label',
            ],
        )
        NetworkEntities.__new__.__defaults__ = (None,)

        collection = self.collect(method = method, **kwargs)
        counts = self.counts(collection_method = collection, **kwargs)

        return NetworkEntities(
            counts = counts,
            entities = collection if keep_collection else None,
            method = method,
        )


    def references_stats(self, **kwargs):

        return self.stats('references', **kwargs)


    def interactions_undirected_stats(self, **kwargs):

        return self.stats('interactions_undirected', **kwargs)


    def interactions_all_stats(self, **kwargs):

        return self.stats('interactions_all', **kwargs)


    def interactions_directed_stats(self, **kwargs):

        return self.stats('interactions_directed', **kwargs)


    def interactions_mutual_stats(self, **kwargs):

        return self.stats('interactions_mutual', **kwargs)


    def interactions_signed_stats(self, **kwargs):

        return self.stats('interactions_signed', **kwargs)


    def interactions_stimulatory_stats(self, **kwargs):

        return self.stats('interactions_stimulatory', **kwargs)


    def interactions_inhibitory_stats(self, **kwargs):

        return self.stats('interactions_inhibitory', **kwargs)


    def entities_stats(self, **kwargs):

        return self.stats('entities', **kwargs)


    def curation_effort_stats(self, **kwargs):

        return self.stats('curation_effort', **kwargs)


    #
    # exporting resource vs. entity counts
    #

    def summaries_tab(self, outfile = None, return_table = False):
        """
        Creates a table from resource vs. entity counts and optionally
        writes it to ``outfile`` and returns it.
        """

        tab = []
        tab.append(self.summaries_labels.values())

        tab.extend([
            [
                str(value)
                for value in self.summaries[resource].values()
            ]
            for resource in sorted(
                self.summaries.keys(),
                key = lambda s: (1 if s == 'Total' else 0, s.lower())
            )
        ])

        if outfile:

            with open(outfile, 'w') as fp:

                fp.write('\n'.join('\t'.join(row) for row in tab))

        if return_table:

            return tab


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
        graph_attrs, vertex_attrs, edge_attrs = (
            graphviz_input.get_graphviz_attrs()
        )
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
        for _entity in ['graph', 'vertex', 'edge']:
            callbacks_dict = '%s_callbacks' % _entity
            _attrs[callbacks_dict] = {}
            callbacks = _attrs[callbacks_dict]

            if _entity == 'edge':

                if (auto_edges == 'RESOURCE_CATEGORIES' or
                        auto_edges == 'DIRECTIONS') \
                        and 'edge_color' not in _custom_attrs \
                        and 'edge_arrowhead' not in _custom_attrs:
                    callbacks['color'] = \
                        _AttrHelper(auto_edges, 'edge_color', _defaults)

                    if auto_edges == 'DIRECTIONS':
                        callbacks['arrowhead'] = \
                            _AttrHelper(auto_edges, 'edge_arrowhead', _defaults)

                    else:
                        callbacks['arrowhead'] = \
                            _AttrHelper('none', 'edge_arrowhead', _defaults)

            for attr in locals()['%s_attrs' % _entity].keys():
                callback_name = '%s_%s' % (_entity, attr)

                if callback_name in _custom_attrs:
                    callback_value = _custom_attrs[callback_name]

                    if isinstance(callback_value, dict):

                        if '_name' not in callback_value:
                            callback_value['_name'] = 'index'
                    callbacks[attr] = _AttrHelper(
                        value=callback_value, name=attr, defaults=_defaults)

        # graph
        dot = graphviz.AGraph(directed=directed)
        attrs = {}

        for gattr, fun in iteritems(_attrs['graph_callbacks']):
            attrs[gattr] = fun(g)

        attrs = common.clean_dict(attrs)

        for gattr, value in iteritems(attrs):
            dot.graph_attr[gattr] = value

        # vertices
        for vid, node in iteritems(dNodes):
            attrs = {}

            for vattr, fun in iteritems(_attrs['vertex_callbacks']):
                attrs[vattr] = fun(g.vs[vid])

            if vid in hide_nodes:
                attrs['style'] = 'invis'

            attrs = common.clean_dict(attrs)
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
                this_directed = 'undirected'
                thisSources = set([])

                if directed:

                    if d.get_dir((sn, tn)):
                        sdir = d.get_dir((sn, tn), sources=True)
                        this_directed = (sn, tn)
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
                                    attrs[eattr] = fun(g.es[eid], this_directed,
                                                       thisSign, thisSources,
                                                       g.es[eid]['sources'])

                            elif vis and (sign[1] and dir_sources is None or
                                          dir_sources is not None and
                                          len(ssign[1] & dir_sources) > 0):
                                thisSign = 'inhibition'
                                thisSoures = ssign[1] if dir_sources is None else \
                                    ssign[1] & dir_sources

                                for eattr, fun in iteritems(edge_callbacks):
                                    attrs[eattr] = fun(g.es[eid], this_directed,
                                                       thisSign, thisSources,
                                                       g.es[eid]['sources'])

                            elif vis:
                                thisSign = 'unknown'
                                thisSources = sdir

                                for eattr, fun in iteritems(edge_callbacks):
                                    attrs[eattr] = fun(g.es[eid], this_directed,
                                                       thisSign, thisSources,
                                                       g.es[eid]['sources'])

                            else:
                                attrs['style'] = 'invis'
                                drawn_directed = False
                            attrs = common.clean_dict(attrs)
                            dot.add_edge(sl, tl, **attrs)

                    if d.get_dir((tn, sn)):
                        sdir = d.get_dir((tn, sn), sources=True)
                        this_directed = (tn, sn)
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
                                    attrs[eattr] = fun(g.es[eid], this_directed,
                                                       thisSign, thisSources,
                                                       (g.es[eid]['sources']))

                            elif vis and (sign[1] and dir_sources is None or
                                          dir_sources is not None and
                                          len(ssign[1] & dir_sources) > 0):
                                thisSign = 'inhibition'
                                thisSources = ssign[1] if dir_sources is None else \
                                    ssign[1] & dir_sources

                                for eattr, fun in iteritems(edge_callbacks):
                                    attrs[eattr] = fun(g.es[eid], this_directed,
                                                       thisSign, thisSources,
                                                       g.es[eid]['sources'])

                            elif vis:
                                thisSign = 'unknown'
                                thisSources = sdir

                                for eattr, fun in iteritems(edge_callbacks):
                                    attrs[eattr] = fun(g.es[eid], this_directed,
                                                       thisSign, thisSources,
                                                       g.es[eid]['sources'])

                            else:
                                attrs['style'] = 'invis'
                                drawn_directed = False
                            attrs = common.clean_dict(attrs)
                            dot.add_edge(tl, sl, **attrs)

                if not directed or d.get_dir('undirected'):
                    attrs = {}
                    this_directed = 'undirected'
                    thisSign = 'unknown'
                    thisSources = d.get_dir('undirected', sources=True)

                    for eattr, fun in iteritems(edge_callbacks):
                        attrs[eattr] = fun(g.es[eid], this_directed, thisSign,
                                           thisSources, g.es[eid]['sources'])

                    if (not evis and hide) or drawn_directed:
                        attrs['style'] = 'invis'

                    if dot.has_neighbor(sl, tl):

                        if dot.get_edge(sl, tl).attr['style'] is not None and \
                                'invis' in dot.get_edge(sl, tl).attr['style']:
                            dot.delete_edge(sl, tl)

                    if not dot.has_neighbor(sl, tl):
                        attrs = common.clean_dict(attrs)
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
        Deprecated, will be removed.
        """

        self.graph.es['in_complex'] = \
            [sum([len(set(self.graph.vs[e.source]['complexes'][cs].keys()) &
                      set(self.graph.vs[e.target]['complexes'][cs].keys()))
                  for cs in csources]) > 0 for e in self.graph.es]

    #
    # Methods for translating network to other organism
    #

    def _translate_refsdir(self, rd, ids):
        """
        Orthology translation of the `references by direction` dictionaries.
        """

        new_refsdir = {}

        for k, v in iteritems(rd):
            di = (ids[k[0]], ids[k[1]]) if type(k) is tuple else k
            new_refsdir[di] = v

            return new_refsdir


    def orthology_translation(
            self,
            target,
            source = None,
            only_swissprot = True,
            graph = None,
        ):
        """
        Translates the current object to another organism by orthology.
        Proteins without known ortholog will be deleted.

        :param int target:
            NCBI Taxonomy ID of the target organism. E.g. 10090 for mouse.
        """

        return_graph = graph is not None
        graph = self.graph if graph is None else graph
        source = self.ncbi_tax_id if source is None else source

        self._log(
            'Translating network by orthology '
            'from taxon `%u` to taxon `%u`.' % (source, target)
        )

        name_old__name_new = orthology.homologene_uniprot_dict(
            source = source,
            target = target,
            only_swissprot = only_swissprot,
        )

        self._log(
            'UniProt to UniProt orthology dictionary obtained from '
            'NCBI Homologene. Contains %u UniProt IDs for taxon `%u`.' % (
                len(name_old__name_new),
                source,
            )
        )

        vcount_before = graph.vcount()
        ecount_before = graph.ecount()

        self._log(
            'Starting translation. Original network consists of %u '
            'nodes and %u edges.' % (vcount_before, ecount_before)
        )

        # nodes could not be mapped are to be deleted
        name_old__vid_old = dict(
            (name, vid)
            for vid, name in enumerate(graph.vs['name'])
        )

        name_old__name_new = dict(
            (name_old, name_new)
            for name_old, name_new in iteritems(name_old__name_new)
            if name_old in name_old__vid_old and name_new
        )

        delete_vids = [
            vid_old
            for name_old, vid_old in iteritems(name_old__vid_old)
            if (
                (
                    name_old not in name_old__name_new or
                    not len(name_old__name_new[name_old])
                ) and
                # nodes of other species or compounds ignored
                graph.vs[vid_old]['ncbi_tax_id'] == source
            )
        ]

        ndel = len(delete_vids)

        self._log(
            'Found orthologues for %u node IDs, %u nodes will be '
            'deleted because no ortholog is known.' % (
                len(name_old__name_new),
                ndel,
            )
        )

        graph.delete_vertices(delete_vids)

        self._log(
            'Number of nodes reduced from %u to %u after '
            'deletion of nodes with no known orthologs.' % (
                vcount_before,
                graph.vcount(),
            )
        )

        graph.vs['old_name'] = copy_mod.deepcopy(graph.vs['name'])

        del delete_vids
        del name_old__vid_old

        # this for permanent identification of nodes:
        graph.vs['id_old'] = list(xrange(graph.vcount()))
        # a dict of these permanent ids and the orthologs:
        vid_old__name_new = dict(
            (
                v['id_old'],
                name_old__name_new[v['name']]
            )
            for v in graph.vs
        )

        # renaming vertices using the first ortholog
        new_names = [
            (
                name_old__name_new[v['name']][0]
                    if v['ncbi_tax_id'] == source else
                # nodes of other species or compounds ignored
                v['name']
            )
            for v in graph.vs
        ]

        graph.vs['name'] = new_names

        self._log(
            '%u nodes renamed to the name of the first ortholog. '
            'New nodes will be created for genes mapping to multiple '
            'orthologues.' % len(new_names)
        )

        del new_names

        # the new nodes to be added because of ambiguous mapping
        new_nodes = list(
            set(
                itertools.chain(*(
                    names_new[1:]
                    for names_new in name_old__name_new.values()
                ))
            ) -
            # except those already exist:
            set(graph.vs['name'])
        )

        graph += new_nodes

        self._log(
            '%u new nodes have been added due to genes mapped to more '
            'than one orthologues. Total node count is %u.' % (
                len(new_nodes),
                graph.vcount(),
            )
        )

        del new_nodes

        # this for permanent identification of nodes:
        graph.vs['id_new'] = list(xrange(graph.vcount()))

        # this is a dict of vertices to be multiplied:
        vid_new_orig__vid_new_all = dict(
            (
                # keys are id_new of the original node
                graph.vs.select(id_old = vid_old)[0]['id_new'],
                # id_new of all orthologs
                [
                    graph.vs.select(name = name_new)[0]['id_new']
                    for name_new in names_new
                ]
            )
            for vid_old, names_new in iteritems(vid_old__name_new)
        )

        self._log(
            'Dictionary of nodes created: '
            '%u nodes are affected by ambiguous mapping.' % (
                len([
                    1
                    for v in vid_new_orig__vid_new_all.values()
                    if len(v) > 1
                ]),
            )
        )

        # compiling a dict of new edges to be added due to ambigous mapping

        # this is for unambiguously identify edges both at directed and
        # undirected graphs after reindexing at adding new edges:
        graph.es['id_orig'] = list(range(graph.ecount()))
        # pairs of new vertex IDs and names for the original edges
        graph.es['vids_new__e_orig'] = [
            (
                graph.vs[e.source]['id_new'],
                graph.vs[e.target]['id_new'],
            )
            for e in graph.es
        ]
        graph.es['names_new__e_orig'] = [
            (
                graph.vs[e.source]['name'],
                graph.vs[e.target]['name'],
            )
            for e in graph.es
        ]

        # the parent edge original id as key
        eid_orig__vids_new_not_orig = dict(
            (
                eid_orig,
                [
                    vids_new__e_new
                    for vids_new__e_new in itertools.product(
                        vid_new_orig__vid_new_all[vids_new__e_orig[0]],
                        vid_new_orig__vid_new_all[vids_new__e_orig[1]],
                    )
                    # not the parent edge itself
                    if not (
                        vids_new__e_new[0] == vids_new__e_orig[0] and
                        vids_new__e_new[1] == vids_new__e_orig[1]
                    )
                ]
            )
            for eid_orig, vids_new__e_orig in (
                (e['id_orig'], e['vids_new__e_orig']) for e in graph.es
            )
        )

        self._log(
            '%u edges have to be mapped to %u other edges due to '
            'ambiguous mapping of their endpoints. The network '
            'consists of %u edges in total.' % (
                len([
                    1
                    for v in eid_orig__vids_new_not_orig.values()
                    if len(v)
                ]),
                sum(len(v) for v in eid_orig__vids_new_not_orig.values()),
                graph.ecount(),
            )
        )

        # translating the dict values to new vertex indices
        vids_missing = list(
            set(
                (
                    graph.vs.select(id_new = s_vid_new)[0].index,
                    graph.vs.select(id_new = t_vid_new)[0].index,
                )
                for s_vid_new, t_vid_new in (
                    itertools.chain(*eid_orig__vids_new_not_orig.values())
                )
            ) -
            # but not any existing edge
            set((e.source, e.target) for e in graph.es) -
            (
                set()
                    if graph.is_directed() else
                set((e.target, e.source) for e in graph.es)
            )
        )

        # creating new edges
        graph += vids_missing

        self._log(
            '%u new edges have been added, number of edges '
            'increased to %u.' % (len(vids_missing), graph.ecount())
        )

        # id_new > current index
        vid_new__vid = dict((v['id_new'], v.index) for v in graph.vs)
        # id_new > name
        vid_new__name = dict((v['id_new'], v['name']) for v in graph.vs)

        prg = Progress(graph.ecount(), 'Translating network by orthology', 21)

        self._log(
            'Copying edge attributes from original edges '
            'to corresponding new edges between orthologues. '
            'Number of edges without attributes: %u.' % (
                len([1 for e in graph.es if e['dirs'] is None]),
            )
        )

        # setting attributes on old and new edges:
        for e in graph.es:

            prg.step()

            d = e['dirs']

            # this lookup is appropriate as old node names are certainly
            # unique; for newly added edges `dirs` will be None
            if (
                d is not None and
                d.nodes[0] in name_old__name_new and
                d.nodes[1] in name_old__name_new
            ):

                # translation of direction object attached to original edges

                ids = {
                    d.nodes[0]: name_old__name_new[d.nodes[0]][0],
                    d.nodes[1]: name_old__name_new[d.nodes[1]][0],
                }

                e['dirs'] = d.translate(ids)

                e['refs_by_dir'] = (
                    self._translate_refsdir(e['refs_by_dir'], ids)
                )

                # if no new edges have been introduced
                # based on this specific edge
                if e['id_orig'] not in eid_orig__vids_new_not_orig:

                    continue

                # iterating new edges between orthologs
                for vid_new_s, vid_new_t in (
                    eid_orig__vids_new_not_orig[e['id_orig']]
                ):

                    # the current vertex indices
                    vid_s = vid_new__vid[vid_new_s]
                    vid_t = vid_new__vid[vid_new_t]
                    # in case of directed graphs this will be correct:
                    new_edges_0 = graph.es.select(
                        _source = vid_s,
                        _target = vid_t,
                    )
                    new_edges_1 = ()

                    if not graph.is_directed():
                        # at undirected graphs
                        # source/target might be opposite:
                        new_edges_1 = graph.es.select(
                        _source = vid_t,
                        _target = vid_s,
                    )

                    if not len(new_edges_0) and not len(new_edges_1):

                        self._log(
                            'Orthology translation: could not find edge '
                            'between %s and %s!' % (
                                graph.vs[vid_s]['name'],
                                graph.vs[vid_t]['name']
                            ),
                            -5,
                        )
                        continue

                    # this is a new edge between orthologs
                    for new_edge in itertools.chain(new_edges_0, new_edges_1):

                        # this mapping supposed to be correct because
                        # the vid_s and vid_t pairs built at the same time
                        # as the `names_new__e_orig` edge attr

                        ids = {
                            e['names_new__e_orig'][0]:
                                graph.vs[vid_s]['name'],
                            e['names_new__e_orig'][1]:
                                graph.vs[vid_t]['name'],
                        }

                        new_edge['dirs'] = e['dirs'].translate(ids)

                        new_edge['refs_by_dir'] = (
                            self._translate_refsdir(e['refs_by_dir'], ids)
                        )
                        new_edge['attrs'] = e['attrs'].translate(ids)

                        # copying the remaining attributes
                        for eattr in e.attributes():

                            if (
                                eattr not in
                                {'dirs', 'refs_by_dir', 'attrs'}
                            ):

                                new_edge[eattr] = copy_mod.deepcopy(e[eattr])

        prg.terminate()

        self._log(
            'Copying edge attributes finished. '
            'Number of edges without attributes (should be zero): %u.' % (
                len([1 for e in graph.es if e['dirs'] is None]),
            )
        )

        # id_new > current index
        vid_new__vid = dict((v['id_new'], v.index) for v in graph.vs)

        self._log(
            'Copying node attributes from original nodes '
            'to their corresponding orthologues. '
            'Number of nodes without attributes: %u' % (
                len([1 for v in graph.vs if v['ncbi_tax_id'] is None])
            )
        )

        # setting attributes of vertices
        for vid_new_orig, vids_new_all in (
            iteritems(vid_new_orig__vid_new_all)
        ):

            # the first ortholog:
            v_orig = graph.vs[vid_new__vid[vid_new_orig]]
            # now setting its taxon to the target:
            v_orig['ncbi_tax_id'] = target

            for vid_new in vids_new_all:
                # iterating further orthologs:
                v_new = graph.vs[vid_new__vid[vid_new]]

                if v_new == v_orig:

                    continue

                # copying attributes:
                for vattr in v_orig.attributes():

                    if vattr != 'name':

                        v_new[vattr] = copy_mod.deepcopy(v_orig[vattr])

        self._log(
            'Copying node attributes from original nodes finished.'
            'Number of nodes without attributes: %u' % (
                len([1 for v in graph.vs if v['ncbi_tax_id'] is None])
            )
        )

        # removing temporary edge attributes
        del self.graph.es['id_orig']
        del self.graph.vs['id_old']
        del self.graph.vs['id_new']
        del self.graph.es['vids_new__e_orig']
        del self.graph.es['names_new__e_orig']

        self._log('Collapsing any duplicate node or edge.')

        self.collapse_by_name(graph = graph)
        self.genesymbol_labels(remap_all = True)

        if not return_graph:

            self._log(
                'Setting the PyPath object`s default taxon to `%u`.' % (
                    target
                )
            )

            self.ncbi_tax_id = target

            self.update_vname()
            self.update_vindex()

            self.clean_graph()

        self._log(
            'Network successfully translated from organism `%u` to'
            ' `%u`. Nodes before: %u, after: %u. Edges before: %u,'
            ' after %u' % (
                source,
                target,
                vcount_before,
                graph.vcount(),
                ecount_before,
                graph.ecount()
            )
        )

        if return_graph:

            return graph


    homology_translation = orthology_translation


    def random_walk_with_return(self, q, graph=None, c=.5, niter=1000):
        """
        Random walk with return (RWR) starting from one or more query nodes.
        Returns affinity (probability) vector of all nodes in the graph.

        Args
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
        _p = copy_mod.copy(_q)

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

    def nodes_by_resource(self, resource):

        return set(
            v['name']
            for v in self.graph.vs
            if resource in v['sources']
        )

    def nodes_by_resources(self):

        return dict(
            (resource, self.nodes_by_resource(resource))
            for resource in self.sources
        )


    def name_edgelist(self, graph = None):
        """
        Returns an edge list, i.e. a list with tuples of vertex names.
        """

        graph = graph or self.graph

        _sort = (lambda _: _) if graph.is_directed() else sorted
        names = graph.vs['name']

        return [
            tuple(_sort((names[e.source], names[e.target])))
            for e in graph.es
        ]

    def __iter__(self):

        return self.iter_interactions()

    def iter_interactions(
            self,
            signs = True,
            all_undirected = True,
            by_source = False,
            with_references = False,
        ):
        """
        Iterates over edges and yields interaction records.

        :param bool signs:
            Ignoring signs if ``False``. This way each directed interaction
            will yield a single record even if it's ambiguously labeled
            both positive and negative. The default behaviour is to yield
            two records in this case, one with positive and one with negative
            sign.
        :param bool all_undirected:
            Yield records for undirected interactions even if certain sources
            provide direction. If ``False`` only the directed records will
            be provided and the undirected sources and references will be
            added to the directed records.
        :param bool by_source:
            Yield separate records by resources. This way the node pairs
            will be redundant and you need to group later if you want
            unique interacting pairs. By default is ``False`` because for
            most applications unique interactions are preferred.
            If ``False`` the *refrences* field will still be present
            but with ``None`` values.
        :param bool with_references:
            Include the literature references. By default is ``False``
            because you rarely need these and they increase the data size
            significantly.
        """


        source_op = operator.eq if by_source else operator.contains


        def get_references(sources, edge, typ):

            return (
                set(
                    ref.pmid
                    for this_refs in
                    (
                        refs
                        for src, refs in iteritems(edge['refs_by_source'])
                        if source_op(sources, src)
                    )
                    for ref in this_refs & edge['refs_by_type'][typ]
                )

                if with_references else

                None
            )


        def iter_sources(sources, edge, typ):

            sources = (
                sources
                    if by_source else
                (sources,)
                    if sources else
                ()
            )

            for _sources in sources:

                yield (
                    _sources,
                    get_references(_sources, edge, typ)
                )


        for edge in self.graph.es:

            directions = edge['dirs']

            for typ, typ_sources in iteritems(edge['sources_by_type']):

                for direction in (directions.straight, directions.reverse):

                    if not directions.dirs[direction]:
                        # this direction does not exist
                        continue

                    dir_sources = directions.get_dir(
                        direction,
                        sources = True,
                    ) & typ_sources

                    id_a = direction[0]
                    id_b = direction[1]
                    type_a = self.uniprot(id_a)['type']
                    type_b = self.uniprot(id_b)['type']

                    for effect, sign_sources in zip(
                        (1, -1),
                        directions.get_sign(direction, sources = True)
                    ):

                        for sources, references in iter_sources(
                            sign_sources & typ_sources,
                            edge,
                            typ,
                        ):

                            yield network.Interaction(
                                id_a = id_a,
                                id_b = id_b,
                                type_a = type_a,
                                type_b = type_b,
                                type = typ,
                                directed = True,
                                effect = effect,
                                sources = sources,
                                references = references,
                            )

                    sources_with_sign = set.union(
                        *directions.get_sign(direction, sources = True)
                    )
                    sources_without_sign = dir_sources - sources_with_sign

                    for sources, references in iter_sources(
                        sources_without_sign,
                        edge,
                        typ,
                    ):

                        yield network.Interaction(
                            id_a = id_a,
                            id_b = id_b,
                            type_a = type_a,
                            type_b = type_b,
                            type = typ,
                            directed = True,
                            effect = 0,
                            sources = sources,
                            references = references,
                        )

                undirected_sources = (
                    directions.get_dir('undirected', sources = True)
                ) & typ_sources

                if not undirected_sources:

                    continue

                id_a = self.graph.vs[edge.source]['name']
                id_b = self.graph.vs[edge.target]['name']
                type_a = self.graph.vs[edge.source]['type']
                type_b = self.graph.vs[edge.target]['type']

                for sources, references in iter_sources(
                    undirected_sources,
                    edge,
                    typ,
                ):

                    yield network.Interaction(
                        id_a = id_a,
                        id_b = id_b,
                        type_a = type_a,
                        type_b = type_b,
                        type = typ,
                        directed = False,
                        effect = 0,
                        sources = sources,
                        references = references,
                    )

    # shortcuts for the most often used igraph attributes:


    def name_to_label(self, name):

        try:

            return self.nodLab[self.nodDct[name]]

        except (KeyError, IndexError):

            return str(name)


    @property
    def vcount(self):

        return self.graph.vcount()


    @property
    def ecount(self):

        return self.graph.ecount()


    @property
    def vs(self):

        return self.graph.vs


    @property
    def es(self):

        return self.graph.es


    @property
    def vertex_attributes(self):

        return self.graph.vertex_attributes()


    @property
    def edge_attributes(self):

        return self.graph.edge_attributes()


    def reload(self):
        """Reloads the object from the module level."""

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    @staticmethod
    def _disclaimer(self):
        """
        Prints a disclaimer about respecting data licences.
        """

        pypath._disclaimer()


    @staticmethod
    def license(self):
        """
        Prints information about data licences.
        """

        pypath.license()


def init_db(use_omnipath = False, **kwargs):

    pa = PyPath()
    getattr(
        pa,
        'load_omnipath' if use_omnipath else 'init_network'
    )(**kwargs)

    globals()['db'] = pa


def get_db(**kwargs):

    if 'db' not in globals():

        init_db(**kwargs)

    return globals()['db']
