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
Here we define one class, the :py:class:`Interaction`` which provides a
rich API for representing and querying molecular interactions. The
interactions serve as the building elements of the network and the
:py:class:`pypath.network.Network` object largely relies on methods
of the :py:class:`Interaction`` objects.
"""

from __future__ import annotations

from future.utils import iteritems

from typing import Literal

import importlib as imp
import collections
import operator
import itertools
import functools
import json

import pypath.core.evidence as pypath_evidence
import pypath.internals.resource as pypath_resource
import pypath.share.session as session_mod
import pypath.share.common as common
import pypath.utils.mapping as mapping
import pypath.core.entity as entity
import pypath.core.attrs as attrs_mod
import pypath.utils.orthology as orthology

_logger = session_mod.Logger(name = 'interaction')
_log = _logger._log


InteractionKey = collections.namedtuple(
    'InteractionKey',
    [
        'entity_a',
        'entity_b',
    ],
)


InteractionDataFrameRecord = collections.namedtuple(
    'InteractionDataFrameRecord',
    [
        'id_a',
        'id_b',
        'type_a',
        'type_b',
        'directed',
        'effect',
        'type',
        'dmodel',
        'sources',
        'references',
    ],
)
InteractionDataFrameRecord.__new__.__defaults__ = (None,) * 8


class Interaction(attrs_mod.AttributeHandler):
    """
    Represents a unique pair of molecular entities interacting with each
    other. One :py:class:`Interaction` object might represent multiple
    interactions i.e. with different direction or effect or type (e.g.
    transcriptional regulation and post-translational regulation),
    each supported by different evidences.

    :arg str,pypath.entity.Entity a,b:
        The two interacting partners. If an :py:class:`pypath.entity.Entity`
        objects provided the other attributes (entity_type, id_type, taxon)
        will be ignored.
    :arg str id_type_a,id_type_b:
        The identifier types for partner ``a`` and ``b`` e.g. ``'uniprot'``.
    :arg str entity_type_a,entity_type_b:
        The types of the molecular entities ``a`` and ``b``
        e.g. ``'protein'``.
    :arg int taxon_a,taxon_b:
        The NCBI Taxonomy Identifiers of partner ``a`` and ``b``
        e.g. ``9606`` for human.

    :details:
        The arguments ``a`` and ``b`` will be assigned to the attribute ``a``
        and ``b`` in an alphabetical order, hence it's possible that
        argument ``a`` becomes attribute ``b``.
    """

    __slots__ = [
        'a',
        'b',
        'a_b',
        'b_a',
        'nodes',
        'key',
        'evidences',
        'direction',
        'positive',
        'negative',
        'unknown_effect',
    ]

    _get_methods = {
        'datasets',
        'entities',
        'evidences',
        'references',
        'curation_effort',
        'resource_names',
        'resources',
        'data_models',
        'interaction_types',
        'interactions',
        'interactions_0',
        'interactions_directed',
        'interactions_undirected',
        'interactions_non_directed',
        'interactions_undirected_0',
        'interactions_non_directed_0',
        'interactions_signed',
        'interactions_positive',
        'interactions_negative',
        'interactions_mutual',
        'resources_via',
        'resource_names_via',
    }

    _get_methods_autogen = (
        'references',
        'resources',
        'resources_via',
        'resource_names',
        'resource_names_via',
        'data_models',
        'interaction_types',
        'datasets',
    )

    _by_methods = (
        'resource',
        'reference',
        'data_model',
        'interaction_type',
        'interaction_type_and_data_model',
        'interaction_type_and_data_model_and_resource',
    )

    _count_methods = {
        'references',
        'resources',
        'resources_via',
        'resource_names_via',
        'resource_names',
        'curation_effort',
        'entities',
        'proteins',
        'complexes',
        'mirnas',
        'interactions',
        'interactions_0',
        'interactions_directed',
        'interactions_signed',
        'interactions_positive',
        'interactions_negative',
        'data_models',
        'interaction_types',
    }

    _get_method_signature = [
        ('direction', None),
        ('effect', None),
        ('resources', None),
        ('data_model', None),
        ('interaction_type', None),
        ('via', None),
        ('references', None),
    ]

    _degree_modes = (
        'ALL',
        'IN',
        'OUT',
    )

    _degree_directions = {
        'undirected': (None, None),
        'non_directed': (False, None),
        'directed': (True, None),
        'signed': (True, True),
        'positive': (True, 'positive'),
        'negative': (True, 'negative'),
    }

    _entity_types = {
        'protein',
        ('complex', 'complexes'),
        'mirna',
        'lncrna',
        # ugly :(
        (('small_molecule', 'drug', 'metabolite'),),
        None,
    }

    _entity_values = {
        'identifiers',
        'labels',
        None,
    }


    def __init__(
            self,
            a,
            b,
            id_type_a = 'uniprot',
            id_type_b = 'uniprot',
            entity_type_a = 'protein',
            entity_type_b = 'protein',
            taxon_a = 9606,
            taxon_b = 9606,
            attrs = None,
        ):

        a = self._get_entity(
            identifier = a,
            id_type = id_type_a,
            entity_type = entity_type_a,
            taxon = taxon_a,
        )
        b = self._get_entity(
            identifier = b,
            id_type = id_type_b,
            entity_type = entity_type_b,
            taxon = taxon_b,
        )

        self.nodes = tuple(sorted((a, b)))
        self.a = self.nodes[0]
        self.b = self.nodes[1]

        self.key = self._key

        self.a_b = (self.nodes[0], self.nodes[1])
        self.b_a = (self.nodes[1], self.nodes[0])

        self.evidences = pypath_evidence.Evidences()
        self.direction = {
            self.a_b: pypath_evidence.Evidences(),
            self.b_a: pypath_evidence.Evidences(),
            'undirected': pypath_evidence.Evidences(),
        }
        self.positive = {
            self.a_b: pypath_evidence.Evidences(),
            self.b_a: pypath_evidence.Evidences(),
        }
        self.negative = {
            self.a_b: pypath_evidence.Evidences(),
            self.b_a: pypath_evidence.Evidences(),
        }
        self.unknown_effect = {
            self.a_b: pypath_evidence.Evidences(),
            self.b_a: pypath_evidence.Evidences(),
        }

        attrs_mod.AttributeHandler.__init__(self, attrs)


    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        evmodname = self.evidences.__class__.__module__
        enmodname = self.a.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        evmod = __import__(evmodname, fromlist = [evmodname.split('.')[0]])
        enmod = __import__(enmodname, fromlist = [enmodname.split('.')[0]])
        imp.reload(mod)
        imp.reload(evmod)
        imp.reload(enmod)
        new = getattr(mod, self.__class__.__name__)
        evsnew = getattr(evmod, 'Evidences')
        evnew = getattr(evmod, 'Evidence')
        ennew = getattr(enmod, 'Entity')
        setattr(self, '__class__', new)

        for evs in itertools.chain(
            (self.evidences,),
            self.direction.values(),
            self.positive.values(),
            self.negative.values(),
            self.unknown_effect.values(),
        ):

            evs.__class__ = evsnew

            for ev in evs:

                ev.__class__ = evnew

        self.a.__class__ = ennew
        self.b.__class__ = ennew

        self._generate_get_methods()
        self._generate_count_methods()
        self._generate_by_methods()


    def _get_entity(
            self,
            identifier,
            id_type = 'uniprot',
            entity_type = 'protein',
            taxon = 9606,
        ):

        if not isinstance(identifier, entity.Entity):

            identifier = entity.Entity(
                identifier = identifier,
                id_type = id_type,
                entity_type = entity_type,
                taxon = taxon,
            )

        return identifier


    def _check_nodes_key(self, nodes):
        """Checks if *nodes* is contained in the edge.

        :arg list nodes:
             Or [tuple], contains the names of the nodes to be checked.

        :return:
            (*bool*) -- ``True`` if all elements in *nodes* are
            contained in the object :py:attr:`nodes` list.
        """

        return nodes == self.a_b or nodes == self.b_a


    def _check_direction_key(self, direction):
        """
        Checks if *direction* is ``'undirected'`` or contains the nodes of
        the current edge. Used internally to check that *di* is a valid
        key for the object attributes declared on dictionaries.

        :arg tuple di:
            Or [str], key to be tested for validity.

        :return:
            (*bool*) -- ``True`` if *di* is ``'undirected'`` or a tuple
              of node names contained in the edge, ``False`` otherwise.
        """

        return (
            direction == 'undirected' or (
                isinstance(direction, tuple) and
                self._check_nodes_key(direction)
            )
        )


    def id_to_entity(self, identifier):

        return (
            self.a
                if self.a == identifier else
            self.b
                if self.b == identifier else
            None
        )


    def direction_key(self, direction):
        """
        The direction keys are tuples of `Entity` objects; this method
        creates these tuples from a tuple of strings. The two strings
        can be labels or identifiers.
        """

        if direction == 'undirected':

            return direction

        direction = tuple(map(self.id_to_entity, direction))

        return (
            direction
                if direction == self.a_b or direction == self.b_a else
            None
        )


    @staticmethod
    def direction_key_identifiers(direction):

        if direction == 'undirected':

            return direction

        return tuple(ent.identifier for ent in direction)


    @staticmethod
    def _directed_key(direction):

        return direction is not None and direction != 'undirected'


    def add_evidence(
            self,
            evidence,
            direction = 'undirected',
            effect = None,
            references = None,
            attrs = None,
        ):
        """
        Adds directionality information with the corresponding data
        source named. Modifies self attributes :py:attr:`dirs` and
        :py:attr:`sources`.

        :arg resource.NetworkResource,evidence.Evidence evidence:
            Either a ``pypath.evidence.Evidence`` object or a resource as
            ``pypath.resource.NetworkResource`` object. In the latter case
            the references can be provided in a separate argument.
        :arg tuple direction:
            Or [str], the directionality key for which the value on
            :py:attr:`dirs` has to be set ``True``.
        :arg int effect:
            The causal effect of the interaction. 1 or 'stimulation'
            corresponds to a stimulatory, -1 or 'inhibition' to an
            inhibitory while 0 to an unknown or neutral effect.
        :arg set,NoneType references:
            A set of references, used only if the resource have been provided
            as ``NetworkResource`` object.
        :arg dict attrs:
            Custom (resource specific) attributes.
        """

        direction = self.direction_key(direction)

        if direction is None:

            _log(
                'Attempting to add evidence with non matching '
                'interaction partners.'
            )
            return

        evidence = (
            evidence
                if isinstance(
                    evidence,
                    (
                        pypath_evidence.Evidence,
                        pypath_evidence.Evidences,
                    )
                ) else
            pypath_evidence.Evidence(
                resource = evidence,
                references = references,
            )
        )

        if hasattr(evidence, 'update_attrs'):

            evidence.update_attrs(attrs)

        self.evidences += evidence
        self.direction[direction] += evidence

        if direction != 'undirected':

            if effect in {1, 'positive', 'stimulation'}:

                self.positive[direction] += evidence

            elif effect in {-1, 'negative', 'inhibition'}:

                self.negative[direction] += evidence

            elif effect in {0, 'unknown'}:

                self.unknown_effect[direction] += evidence


    def __hash__(self):

        return hash(self.key)


    def __eq__(self, other):

        return hasattr(other, 'key') and self.key == other.key


    @property
    def _key(self):

        return InteractionKey(
            self.a.key,
            self.b.key,
        )


    def __iadd__(self, other):

        if self != other:

            _log(
                'Attempt to merge interactions with '
                'non matching interaction partners.'
            )
            return self

        self._merge_evidences(self, other)
        attrs_mod.AttributeHandler.__iadd__(self, other)

        return self


    def __add__(self, other):

        new = self.__copy__()

        new += other

        return new


    def __copy__(self):

        new = Interaction(*self.key)
        new += self
        new.update_attrs(self)

        return new


    @staticmethod
    def _merge_evidences(one, other):

        one.evidences += other.evidences

        for dir_key in one.direction.keys():

            one.direction[dir_key] += other.direction[dir_key]

        for eff_key in one.positive.keys():

            one.positive[eff_key] += other.positive[eff_key]

        for eff_key in one.negative.keys():

            one.negative[eff_key] += other.negative[eff_key]

        for eff_key in one.unknown_effect.keys():

            one.unknown_effect[eff_key] += other.unknown_effect[eff_key]


    def __repr__(self):

        return '<Interaction: %s [%s]>' % (
            self.__str__(),
            self.evidences.__repr__().strip('<>'),
        )


    def __str__(self):

        return '%s %s=%s=%s=%s %s' % (
            self.a.label or self.a.identifier,
            '<' if self.direction[self.b_a] else '=',
            (
                '(+-)' if (
                    self.positive[self.b_a] and
                    self.negative[self.b_a]
                ) else
                '(+)=' if self.positive[self.b_a] else
                '(-)=' if self.negative[self.b_a] else
                '===='
            ),
            (
                '(+-)' if (
                    self.positive[self.a_b] and
                    self.negative[self.a_b]
                ) else
                '(+)=' if self.positive[self.a_b] else
                '(-)=' if self.negative[self.a_b] else
                '===='
            ),
            '>' if self.direction[self.a_b] else '=',
            self.b.label or self.b.identifier,
        )

    def __contains__(self, other):

        return (
            other == self.a or
            other == self.b or
            self.evidences.__contains__(other) or
            attrs_mod.AttributeHandler.__contains__(self, other)
        )


    def has_data_model(self, data_model):

        return self.evidences.has_data_model(data_model)


    @property
    def data_models(self):

        return {
            ev.resource.data_model
            for ev in self.evidences
        }


    def has_dataset(
            self,
            dataset: str,
            direction: str | tuple = None,
            effect: Literal['positive', 'negative'] = None,
            **kwargs,
        ) -> bool:

        return (
            self.get_evidences(direction = direction, effect = effect).
            has_dataset(dataset, **kwargs)
        )


    def get_direction(
            self,
            direction,
            resources = False,
            evidences = False,
            sources = False,
            resource_names = False,
        ):
        """
        Returns the state (or *resources* if specified) of the given
        *direction*.

        :arg tuple direction:
            Or [str] (if ``'undirected'``). Pair of nodes from which
            direction information is to be retrieved.
        :arg bool resources:
            Optional, ``'False'`` by default. Specifies if the
            :py:attr:`resources` information of the given direction is to
            be retrieved instead.

        :return:
            (*bool* or *set*) -- (if ``resources=True``). Presence/absence
            of the requested direction (or the list of resources if
            specified). Returns ``None`` if *direction* is not valid.
        """

        direction = self.direction_key(direction)

        if direction is not None:

            return self._select_answer_type(
                self.direction[direction],
                resources = resources,
                evidences = evidences,
                resource_names = resource_names,
                sources = sources,
            )

        else:
            return None


    def get_directions(
            self,
            src,
            tgt,
            resources = False,
            evidences = False,
            resource_names = False,
            sources = False,
        ):
        """
        Returns all directions with boolean values or list of sources.

        :arg str src:
            Source node.
        :arg str tgt:
            Target node.
        :arg bool resources:
            Optional, ``False`` by default. Specifies whether to return
            the :py:attr:`resources` attribute instead of :py:attr:`dirs`.

        :return:
            Contains the :py:attr:`dirs` (or :py:attr:`resources` if
            specified) of the given edge.
        """

        query = (src, tgt)

        answer_type_args = {
            'resources': resources,
            'evidences': evidences,
            'resource_names': resource_names,
            'sources': sources,
        }

        query = self.direction_key(query)

        if query is not None:

            return [
                self._select_answer_type(
                    self.direction[query],
                    **answer_type_args
                ),
                self._select_answer_type(
                    self.direction[tuple(reversed(query))],
                    **answer_type_args
                ),
                self._select_answer_type(
                    self.direction['undirected'],
                    **answer_type_args
                ),
            ]

        else:
            return None


    def _select_answer_type(
            self,
            answer,
            resources = False,
            evidences = False,
            resource_names = False,
            sources = False,
        ):

        return (
            answer
                if evidences else
            answer.get_resources()
                if resources else
            answer.get_resource_names()
                if sources or resource_names else
            bool(answer)
        )


    def which_directions(
            self,
            resources = None,
            effect = None,
        ):
        """
        Returns the pair(s) of nodes for which there is information
        about their directionality.

        :arg str effect:
            Either *positive* or *negative*.
        :arg str,set resources:
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
            for _dir, _evidences in iteritems(self.direction)
            if _dir != 'undirected' and
            _evidences and (
                not resources or
                _evidences & resources
            ) and (
                not effect
                or (
                    not resources and
                    getattr(self, effect)[_dir]
                ) or
                getattr(self, effect)[_dir] & resources
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
            for _dir, _evidences in iteritems(getattr(self, _effect))
            if _evidences and (
                not resources or
                _evidences & resources
            )
        )


    @staticmethod
    def _effect_synonyms(effect):

        if not effect or effect == True:

            return effect

        if effect in {'positive', 'stimulation', 'stimulatory'}:

            return 'positive'

        if effect in {'negative', 'inhibition', 'inhibitory'}:

            return 'negative'

        if effect in {'unknown'}:

            return 'unknown'


    def _resources_set(self, resources = None):

        return common.to_set(resources)


    def unset_direction(
            self,
            direction,
            only_sign = False,
            resource = None,
            interaction_type = None,
            via = False,
            source = None,
        ):
        """
        Removes directionality and/or source information of the
        specified *direction*. Modifies attribute :py:attr:`dirs` and
        :py:attr:`sources`.

        :arg tuple direction:
            Or [str] (if ``'undirected'``) the pair of nodes specifying
            the directionality from which the information is to be
            removed.
        :arg set resource:
            Optional, ``None`` by default. If specified, determines
            which specific source(s) is(are) to be removed from
            :py:attr:`sources` attribute in the specified *direction*.
        """

        direction = self.direction_key(direction)

        if direction is not None:

            attrs = (
                (self._effect_synonyms(only_sign),)
                    if only_sign else
                ('direction', 'positive', 'negative', 'unknown_effect')
            )
            resource = resource or source

            for attr in attrs:

                if resource is not None:

                    getattr(self, attr)[direction].remove(
                        resource = resource,
                        interaction_type = interaction_type,
                        via = via,
                    )

                else:
                    getattr(self, attr)[direction] = (
                        pypath_evidence.Evidences()
                    )


    # synonym: old name
    unset_dir = unset_direction


    def unset_sign(
            self,
            direction,
            sign,
            resource = None,
            interaction_type = None,
            via = False,
            source = None,
        ):
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

        self.unset_direction(
            direction = direction,
            only_sign = sign,
            resource = resource,
            interaction_type = interaction_type,
            via = via,
            source = source,
        )


    def unset_interaction_type(self, interaction_type):
        """
        Removes all evidences with a certain ``interaction_type``.
        """

        for ev in tuple(self.evidences):

            if ev.resource.interaction_type == interaction_type:

                self.evidences -= ev

        for attr in ('direction', 'positive', 'negative'):

            for key, evs in getattr(self, attr):

                for ev in tuple(evs):

                    if ev.resource.interaction_type == interaction_type:

                        evs -= ev


    def is_directed(self):
        """
        Checks if edge has any directionality information.

        :return:
            (*bool*) -- Returns ``True`` if any of the :py:attr:`dirs`
            attribute values is ``True`` (except ``'undirected'``),
            ``False`` otherwise.
        """

        return any(
            evs
            for dkey, evs in iteritems(self.direction)
            if dkey != 'undirected'
        )


    def is_directed_by_resources(self, resources = None):
        """
        Checks if edge has any directionality information from some
        resource(s).

        :return:
            (*bool*) -- Returns ``True`` if any of the :py:attr:`dirs`
            attribute values is ``True`` (except ``'undirected'``),
            ``False`` otherwise.
        """

        return self._dir_by_resource(resources, op = operator.or_)


    def is_mutual(self, resources = None):
        """
        Checks if the edge has mutual directions (both A-->B and B-->A).
        """

        return (
            bool(self.direction[self.a_b]) and bool(self.direction[self.b_a])
                if not resources else
            self.is_mutual_by_resources(resources = resources)
        )


    def is_mutual_by_resources(self, resources = None):
        """
        Checks if the edge has mutual directions (both A-->B and B-->A)
        according to some resource(s).
        """

        return self._dir_by_resource(resources, op = operator.and_)


    def is_loop(self):
        """
        :returns:
        ``True`` if the interaction is a loop edge i.e. its endpoints are the
        same node.
        """

        return self.a == self.b


    def _dir_by_resource(self, resources = None, op = operator.or_):

        resources = self._resources_set(resources)

        return op(
            self.direction[self.a_b] & resources,
            self.direction[self.b_a] & resources,
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

        _sign = getattr(self, sign)
        _resources = self._resources_set(resources)

        return (
            any(
                bool(
                    _evidences
                        if not _resources else
                    _evidences & _resources
                )
                for _direction, _evidences in iteritems(_sign)
                if not direction or direction == _direction
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


    def add_sign(
            self,
            direction: tuple,
            sign: str,
            resource: str | set[str] | None = None,
            resource_name: str = None,
            interaction_type: str = 'PPI',
            data_model: str = None,
            attrs: dict = None,
            **kwargs
        ):
        """
        Sets sign and source information on a given direction of the
        edge. Modifies the attributes :py:attr:`positive` and
        :py:attr:`positive_sources` or :py:attr:`negative` and
        :py:attr:`negative_sources` depending on the sign. Direction is
        also updated accordingly, which also modifies the attributes
        :py:attr:`dirs` and :py:attr:`sources`.

        Args
            direction:
                Pair of edge nodes specifying the direction from which the
                information is to be set/updated.
            sign:
                Specifies the type of interaction. Either ``'positive'`` or
                ``'negative'``.
            resource:
                Contains the name(s) of the source(s) from which the
                information was obtained.
            attrs:
                Custom (resource specific) edge attributes.
            kwargs:
                Passed to ``pypath.resource.NetworkResource`` if ``resource``
                is not already a ``NetworkResource`` or ``Evidence``
                instance.
        """

        sign = self._effect_synonyms(sign)

        evidence = (
            resource
                if isinstance(resource, pypath_evidence.Evidence) else
            pypath_evidence.Evidence(
                resource = resource,
                references = references,
            )
                if isinstance(resource, pypath_resource.NetworkResource) else
            pypath_evidence.Evidence(
                resource = pypath_resource.NetworkResource(
                    name = resource_name,
                    interaction_type = interaction_type,
                    data_model = data_model,
                    **kwargs
                )
            )
                if resource_name is not None else
            None
        )

        evidence.update_attrs(attrs)

        direction = self.direction_key(direction)

        if self._directed_key(direction) and evidence is not None:

            ev_attr = getattr(self, sign)
            ev_attr += evidence


    def get_sign(
            self,
            direction: tuple,
            sign: str = None,
            evidences: bool = False,
            resources: bool = False,
            resource_names: bool = False,
            sources: bool = False,
        ) -> list[bool | str]:
        """
        Retrieves the sign information of the edge in the given
        diretion. If specified in *sign*, only that sign's information
        will be retrieved. If specified in *sources*, the sources of
        that information will be retrieved instead.

        Args
            direction:
                Contains the pair of nodes specifying the directionality of
                the edge from which th information is to be retrieved.
            sign:
                Optional, ``None`` by default. Denotes whether to retrieve
                the ``'positive'`` or ``'negative'`` specific information.
            resources:
                Optional, ``False`` by default. Specifies whether to return
                the resources instead of sign.

        Returns
            If ``sign=None`` containing [bool] values
            denoting the presence of positive and negative sign on that
            direction, if ``sources=True`` the [set] of sources for each
            of them will be returned instead. If *sign* is specified,
            returns [bool] or [set] (if ``sources=True``) of that
            specific direction and sign.
        """

        sign = self._effect_synonyms(sign)

        answer_type_args = {
            'resources': resources,
            'evidences': evidences,
            'resource_names': resource_names,
            'sources': sources,
        }

        direction = self.direction_key(direction)

        if self._directed_key(direction):

            return (

                self._select_answer_type(
                    getattr(self, sign)[direction],
                    **answer_type_args
                )

                    if sign else

                [
                    self._select_answer_type(
                        self.positive[direction],
                        **answer_type_args
                    ),
                    self._select_answer_type(
                        self.negative[direction],
                        **answer_type_args
                    )
                ]

            )


    def source(
            self,
            undirected: bool = False,
            resources = None,
            **kwargs
        ):
        """
        Returns the name(s) of the source node(s) for each existing
        direction on the interaction.

        Args
            undirected:
                Optional, ``False`` by default.

        Returns
            (*list*) -- Contains the name(s) for the source node(s).
            This means if the interaction is bidirectional, the list
            will contain both identifiers on the edge. If the
            interaction is undirected, an empty list will be returned.
        """

        return self._partner(
            source_target = 'source',
            undirected = undirected,
            resources = resources,
            **kwargs
        )


    # synonym: old name
    src = source


    def target(
            self,
            undirected: bool = False,
            resources = None,
            **kwargs
        ):
        """
        Returns the name(s) of the target node(s) for each existing
        direction on the interaction.

        Args
            undirected:
                Optional, ``False`` by default.

        Returns
            (*list*) -- Contains the name(s) for the target node(s).
            This means if the interaction is bidirectional, the list
            will contain both identifiers on the edge. If the
            interaction is undirected, an empty list will be returned.
        """

        return self._partner(
            source_target = 'target',
            undirected = undirected,
            resources = resources,
            **kwargs
        )


    # synonym: old name
    tgt = target


    def _partner(
            self,
            source_target,
            undirected = False,
            resources = None,
            **kwargs
        ):

        resources = self._resources_set(resources)
        _slice = slice(0, 1) if source_target == 'source' else slice(1, 2)

        return tuple(itertools.chain(
            (
                _direction[_slice]
                    if _direction != 'undirected' else
                self.nodes
                    if undirected else
                ()
            )
            for _direction, _evidences in iteritems(self.direction)
            if (
                (
                    (
                        not resources and
                        not kwargs and
                        bool(_evidences)
                    ) or
                    (
                        any(
                            ev.match(
                                resource = res,
                                **kwargs
                            )
                            for res in resources or (None,)
                            for ev in _evidences
                        )
                    )
                )
            )
        ))


    def src_by_resource(self, resource):
        """
        Returns the name(s) of the source node(s) for each existing
        direction on the interaction for a specific *resource*.

        :arg str resource:
            Name of the resource according to which the information is to
            be retrieved.

        :return:
            (*list*) -- Contains the name(s) for the source node(s)
            according to the specified *resource*. This means if the
            interaction is bidirectional, the list will contain both
            identifiers on the edge. If the specified *source* is not
            found or invalid, an empty list will be returned.
        """

        return [
            _dir[0]
            for _dir, _evidences in iteritems(self.direction)
            if (
                _dir != 'undirected' and
                resource in _evidences
            )
        ]


    def tgt_by_resource(self, resource):
        """
        Returns the name(s) of the target node(s) for each existing
        direction on the interaction for a specific *resource*.

        :arg str resource:
            Name of the resource according to which the information is to
            be retrieved.

        :return:
            (*list*) -- Contains the name(s) for the target node(s)
            according to the specified *resource*. This means if the
            interaction is bidirectional, the list will contain both
            identifiers on the edge. If the specified *source* is not
            found or invalid, an empty list will be returned.
        """

        return [
            _dir[1]
            for _dir, _evidences in iteritems(self.direction)
            if (
                _dir != 'undirected' and
                resource in _evidences
            )
        ]


    def resources_a_b(
            self,
            resources = False,
            evidences = False,
            resource_names = False,
            sources = False,
        ):
        """
        Retrieves the list of resources for the :py:attr:`a_b`
        direction.

        :return:
            (*set*) -- Contains the names of the sources supporting the
            :py:attr:`a_b` directionality of the edge.
        """

        answer_type_args = {
            'resources': resources,
            'evidences': evidences,
            'resource_names': resource_names,
            'sources': sources,
        }

        return self._select_answer_type(
            self.direction[self.a_b],
            **answer_type_args
        )


    # synonym for old method name
    sources_straight = resources_a_b


    def resources_b_a(
            self,
            resources = False,
            evidences = False,
            resource_names = False,
            sources = False,
        ):
        """
        Retrieves the list of sources for the :py:attr:`b_a` direction.

        :return:
            (*set*) -- Contains the names of the sources supporting the
            :py:attr:`b_a` directionality of the edge.
        """

        answer_type_args = {
            'resources': resources,
            'evidences': evidences,
            'resource_names': resource_names,
            'sources': sources,
        }

        return self._select_answer_type(
            self.direction[self.b_a],
            **answer_type_args
        )


    # synonym for old method name
    sources_reverse = resources_b_a


    def resources_undirected(
            self,
            resources = False,
            evidences = False,
            resource_names = False,
            sources = False,
        ):
        """
        Retrieves the list of resources without directed information.

        :return:
            (*set*) -- Contains the names of the sources supporting the
            edge presence but without specific directionality
            information.
        """

        answer_type_args = {
            'resources': resources,
            'evidences': evidences,
            'resource_names': resource_names,
            'sources': sources,
        }

        return self._select_answer_type(
            self.direction['undirected'],
            **answer_type_args
        )


    sources_undirected = resources_undirected


    def positive_a_b(self):
        """
        Checks if the :py:attr:`a_b` directionality is a positive
        interaction.

        :return:
            (*bool*) -- ``True`` if there is supporting information on
            the :py:attr:`a_b` direction of the edge as activation.
            ``False`` otherwise.
        """

        return bool(self.positive[self.a_b])


    positive_straight = positive_a_b


    def positive_b_a(self):
        """
        Checks if the :py:attr:`b_a` directionality is a positive
        interaction.

        :return:
            (*bool*) -- ``True`` if there is supporting information on
            the :py:attr:`b_a` direction of the edge as activation.
            ``False`` otherwise.
        """

        return bool(self.positive[self.b_a])


    positive_reverse = positive_b_a


    def negative_a_b(self):
        """
        Checks if the :py:attr:`a_b` directionality is a negative
        interaction.

        :return:
            (*bool*) -- ``True`` if there is supporting information on
            the :py:attr:`a_b` direction of the edge as inhibition.
            ``False`` otherwise.
        """

        return bool(self.negative[self.a_b])


    negative_straight = negative_a_b


    def negative_b_a(self):
        """
        Checks if the :py:attr:`b_a` directionality is a negative
        interaction.

        :return:
            (*bool*) -- ``True`` if there is supporting information on
            the :py:attr:`b_a` direction of the edge as inhibition.
            ``False`` otherwise.
        """

        return bool(self.negative[self.b_a])


    negative_reverse = negative_b_a


    def negative_resources_a_b(self, **kwargs):
        """
        Retrieves the list of resources for the :py:attr:`a_b`
        direction and negative sign.

        :return:
            (*set*) -- Contains the names of the resources supporting the
            :py:attr:`a_b` directionality of the edge with a
            negative sign.
        """

        answer_type_args = {
            'resource_names': True
        }
        answer_type_args.update(kwargs)

        return self._select_answer_type(
            self.negative[self.a_b],
            **answer_type_args
        )


    def negative_resources_b_a(self, **kwargs):
        """
        Retrieves the list of resources for the :py:attr:`b_a`
        direction and negative sign.

        :return:
            (*set*) -- Contains the names of the resources supporting the
            :py:attr:`b_a` directionality of the edge with a
            negative sign.
        """

        answer_type_args = {
            'resource_names': True
        }
        answer_type_args.update(kwargs)

        return self._select_answer_type(
            self.negative[self.b_a],
            **answer_type_args
        )


    def positive_resources_a_b(self, **kwargs):
        """
        Retrieves the list of resources for the :py:attr:`a_b`
        direction and positive sign.

        :return:
            (*set*) -- Contains the names of the resources supporting the
            :py:attr:`a_b` directionality of the edge with a
            positive sign.
        """

        answer_type_args = {
            'resource_names': True
        }
        answer_type_args.update(kwargs)

        return self._select_answer_type(
            self.positive[self.a_b],
            **answer_type_args
        )


    def positive_resources_b_a(self, **kwargs):
        """
        Retrieves the list of resources for the :py:attr:`b_a`
        direction and positive sign.

        :return:
            (*set*) -- Contains the names of the resources supporting the
            :py:attr:`b_a` directionality of the edge with a
            positive sign.
        """

        answer_type_args = {
            'resource_names': True
        }
        answer_type_args.update(kwargs)

        return self._select_answer_type(
            self.positive[self.b_a],
            **answer_type_args
        )


    def majority_dir(
            self,
            only_interaction_type = None,
            only_primary = True,
            by_references = False,
            by_reference_resource_pairs = True,
        ):
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


        a_b = self.direction[self.a_b]
        b_a = self.direction[self.b_a]

        if not a_b and not b_a:

            return 'undirected'

        method = (
            'count_references'
                if by_references else
            'count_curation_effort'
                if by_reference_resource_pairs else
            'count_resources'
        )

        n_a_b = getattr(a_b, method)(
            interaction_type = only_interaction_type,
            via = False if only_primary else None,
        )
        n_b_a = getattr(b_a, method)(
            interaction_type = only_interaction_type,
            via = False if only_primary else None,
        )

        return (
            'undirected'
                if n_a_b == 0 and n_b_a == 0 else
            None
                if n_a_b == n_b_a else
            self.a_b
                if n_a_b > n_b_a else
            self.b_a
        )


    def majority_sign(
            self,
            only_interaction_type = None,
            only_primary = True,
            by_references = False,
            by_reference_resource_pairs = True,
        ):
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

        result = {}

        method = (
            'count_references'
                if by_references else
            'count_curation_effort'
                if by_reference_resource_pairs else
            'count_resources'
        )

        for _dir in (self.a_b, self.b_a):

            n_pos = getattr(self.positive[_dir], method)(
                interaction_type = only_interaction_type,
                via = False if only_primary else None,
            )
            n_neg = getattr(self.negative[_dir], method)(
                interaction_type = only_interaction_type,
                via = False if only_primary else None,
            )

            result[_dir] = [
                0 < n_pos >= n_neg,
                0 < n_neg >= n_pos,
            ]

        return result


    def consensus(
            self,
            only_interaction_type = None,
            only_primary = False,
            by_references = False,
            by_reference_resource_pairs = True,
        ):
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

        _dir = self.majority_dir(
            only_interaction_type = only_interaction_type,
            only_primary = only_primary,
            by_references = by_references,
            by_reference_resource_pairs = by_reference_resource_pairs,
        )

        if _dir == 'undirected' or _dir is None:

            _dir = self.majority_dir(
                only_interaction_type = only_interaction_type,
                only_primary = only_primary,
                by_references = False,
                by_reference_resource_pairs = False,
            )

        _effect = self.majority_sign(
            only_interaction_type = only_interaction_type,
            only_primary = only_primary,
            by_references = by_references,
            by_reference_resource_pairs = by_reference_resource_pairs,
        )
        _effect_noref = self.majority_sign(
            only_interaction_type = only_interaction_type,
            only_primary = only_primary,
            by_references = False,
            by_reference_resource_pairs = False,
        )

        if _dir == 'undirected':

            result.append([
                self.a_b[0],
                self.a_b[1],
                'undirected',
                'unknown',
            ])

        else:

            dirs = (self.a_b, self.b_a) if _dir is None else (_dir,)

            for d in dirs:

                d_effect = (
                    _effect[d]
                        if (
                            _effect[d] is not None and
                            _effect[d][0] != _effect[d][1]
                        ) else
                    _effect_noref[d]
                )

                if d_effect is not None:

                    # index #0 is positive
                    if d_effect[0]:

                        result.append([
                            d[0],
                            d[1],
                            'directed',
                            'positive',
                        ])

                    # can not be elif bc of the case of equal weight of
                    # evidences for both positive and negative
                    if d_effect[1]:

                        result.append([
                            d[0],
                            d[1],
                            'directed',
                            'negative',
                        ])

                # directed with unknown effect
                else:

                    result.append([
                        d[0],
                        d[1],
                        'directed',
                        'unknown',
                    ])

        return result


    consensus_edges = consensus


    def merge(self, other):
        """
        Merges current Interaction with another (if and only if they are the
        same class and contain the same nodes). Updates the attributes
        :py:attr:`direction`, :py:attr:`positive` and :py:attr:`negative`.

        :arg pypath.interaction.Interaction other:
            The new Interaction object to be merged with the current one.
        """


        if not self._check_nodes_key(other.nodes):

            _log(
                'Attempting to merge Interaction instances with different '
                'interacting partners.'
            )
            return

        self.evidences += other.evidences

        for attr, _dir in itertools.product(
            ('direction', 'positive', 'negative'),
            (self.a_b, self.b_a, 'undirected')
        ):

            if attr != 'direction' and _dir == 'undirected':

                continue

            getattr(self, attr)[_dir] += getattr(other, attr)[_dir]


    def translate(self, ids, new_attrs = None):
        """
        Translates the node names/identifiers according to the
        dictionary *ids*. Also is able to change attributes like `id_type`,
        `taxon` and `entity_type`.

        :arg dict ids:
            Dictionary containing (at least) the current names of the
            nodes as keys and their translation as values.
        :arg dict new_attrs:
            Dictionary with new IDs as keys and their dicts of their new
            attributes as values. For any attribute not provided here
            the attributes from the original instance will be used.
            E.g. you can provide `{'1956': {'id_type': 'entrez'}}' if the
            new ID type for protein EGFR is Entrez Gene ID.

        :return:
            (*pypath.main.Direction*) -- The copy of current edge object
            with translated node names.
        """

        new_a = ids[self.nodes[0]]
        new_b = ids[self.nodes[1]]
        new_ids = {'a': new_a, 'b': new_b}
        to_old = common.swap_dict_simple(ids)

        all_new_attrs = dict(
            (
                '%s_%s' % (attr, label),
                new_attrs[new_ids[label]][attr]
                    if (
                        new_ids[label] in new_attrs and
                        attr in new_attrs[new_ids[label]]
                    ) else
                getattr(getattr(self, label), attr)
            )
            for attr in ('id_type', 'entity_type', 'taxon')
            for label in ('a', 'b')
        )

        all_new_attrs['attrs'] = new_attrs.get('attrs', self.attrs)

        new = Interaction(
            a = new_a,
            b = new_b,
            **all_new_attrs
        )

        new.evidences += self.evidences

        to_old = dict(
            (
                new_id,
                self.id_to_entity(old_id)
            )
            for new_id, old_id in iteritems(to_old)
        )

        # this is required to handle also loop edges
        new_old_a_b = (
            (
                to_old[new.a.identifier],
                to_old[new.b.identifier],
            )
                if new.a != new.b else
            self.a_b
        )
        new_old_b_a = (
            (
                to_old[new.b.identifier],
                to_old[new.a.identifier],
            )
                if new.a != new.b else
            self.b_a
        )

        for (old_dir, new_dir), attr in itertools.product(
            zip(
                (
                    new_old_a_b,
                    new_old_b_a,
                    'undirected'
                ),
                (
                    new.a_b,
                    new.b_a,
                    'undirected',
                ),
            ),
            ('direction', 'positive', 'negative'),
        ):

            if old_dir == 'undirected' and attr != 'direction':

                continue

            getattr(new, attr)[new_dir] += getattr(self, attr)[old_dir]

        return new


    def orthology_translate_one(self, id_a, id_b, taxon):

        return self.translate(
            ids = {
                self.a: id_a,
                self.b: id_b,
            },
            new_attrs = {
                id_a: {
                    'taxon': (
                        self.a.taxon
                            if id_a == self.a.identifier else
                        taxon
                    ),
                },
                id_b: {
                    'taxon': (
                        self.b.taxon
                            if id_b == self.b.identifier else
                        taxon
                    ),
                },
            },
        )


    def orthology_translate(self, taxon, exclude = None):

        exclude = exclude or set()
        exclude.add(0)

        for new_a, new_b in itertools.product(
            (self.a.identifier,)
                if self.a.taxon in exclude else
            orthology.translate(
                source_id = self.a.identifier,
                target = taxon,
                source = self.a.taxon,
            ),
            (self.b.identifier,)
                if self.b.taxon in exclude else
            orthology.translate(
                source_id = self.b.identifier,
                target = taxon,
                source = self.b.taxon,
            ),
        ):

            yield self.orthology_translate_one(
                id_a = new_a,
                id_b = new_b,
                taxon = taxon,
            )


    def get_evidences(
            self,
            direction = None,
            effect = None,
            resources = None,
            data_model = None,
            interaction_type = None,
            via = None,
            references = None,
            datasets = None,
        ):

        effect = self._effect_synonyms(effect)

        evidences = (

            # any signed
            sum(itertools.chain(
                self.positive.values(),
                self.negative.values(),
            ))

                if effect == True else

            # only positive
            (
                self.positive[direction]
                    if direction in self.positive else
                sum(self.positive.values())
            )

                if effect == 'positive' else

            # only negative
            (
                self.negative[direction]
                    if direction in self.negative else
                sum(self.negative.values())
            )

                if effect == 'negative' else

            # any directed
            sum(self.direction[_dir] for _dir in self.which_dirs())

                if direction == True else

            # one specific direction
            self.direction[direction]

                if direction in self.direction else

            # all evidences (default)
            self.evidences

        )

        return (
            pypath_evidence.Evidences(
                evidences.filter(
                    resource = resources,
                    interaction_type = interaction_type,
                    via = via,
                    data_model = data_model,
                    references = references,
                    datasets = datasets,
                )
            )
        )


    def get_entities(
            self,
            entity_type = None,
            direction = None,
            effect = None,
            resources = None,
            data_model = None,
            interaction_type = None,
            via = None,
            references = None,
            datasets = None,
            return_type = None,
        ):
        """
        Retrieves the entities involved in interactions matching the criteria.
        It either returns both interacting entities in a *set* or an empty
        *set*. This may not sound so useful at the level of this object but
        becomes more useful once we want to collect entities having certain
        kind of interactions across a series of `Interaction` objects.

        :arg str entity_type:
            The type of the molecular entity. Possible values: `protein`,
            `complex`, `mirna`, `small_molecule`.
        :arg str return_type:
            The type of values to return. Default is
            py:class:``pypath.entity.Entity`` objects, alternatives are
            ``labels``  ``identifiers``.
        """

        # TODO: this method could be made slightly more efficient by using
        # not ``get_interactions`` but a simpler logic as here we don't need
        # to handle the directions separately; however this is not very
        # important and for the time being it's good as it is.

        kwargs = locals()
        _ = kwargs.pop('self')
        entity_type = common.to_set(kwargs.pop('entity_type'))
        return_type = kwargs.pop('return_type')

        return_types = {
            'entity': None,
            'entities': None,
            'id': 'identifier',
            'name': 'identifier',
        }

        # allow plurals
        return_type = (
            return_type[:-1]
                if (
                    isinstance(return_type, str) and
                    return_type[-1] == 's'
                ) else
            return_type
        )
        # allow some synonyms
        return_type = (
            return_types[return_type]
                if return_type in return_types else
            return_type
        )

        return (
            set(
                (
                    getattr(en, return_type)
                        if return_type and hasattr(en, return_type) else
                    en
                )
                for en in
                itertools.chain(
                    *self.get_interactions(**kwargs)
                )
                if not entity_type or en.entity_type in entity_type
            )
        )


    @classmethod
    def _generate_entity_methods(cls):

        def _create_entity_method(entity_type, return_type):


            def _entity_method(*args, **kwargs):

                self = args[0]
                kwargs['entity_type'] = entity_type
                kwargs['return_type'] = return_type

                return self.get_entities(*args[1:], **kwargs)


            return _entity_method


        for etype, vtype in itertools.product(
            cls._entity_types,
            cls._entity_values,
        ):

            if etype is None and vtype is None:

                continue

            entity_type = common.sfirst(etype)
            entity_type_label = common.sfirst(entity_type)
            return_type = vtype

            etype_part = (
                ''
                    if not entity_type_label else
                entity_type_label
                    if vtype else
                etype[1]
                    if isinstance(etype, tuple) and len(etype) > 1 else
                '%ss' % entity_type_label
            )
            vtype_part = '%s' % vtype if vtype else ''

            _method_name = '%s%s%s' % (
                etype_part,
                '_' if etype_part and vtype_part else '',
                vtype_part,
            )
            method_name = 'get_%s' % _method_name
            method = _create_entity_method(
                entity_type = entity_type,
                return_type = return_type,
            )

            cls._add_method(
                method_name,
                method,
                signature = (
                    ['self', ('entity_type', None)] +
                    cls._get_method_signature +
                    [('return_type', None)]
                ),
                doc = cls.get_entities.__doc__,
            )
            cls._get_methods.add(_method_name)
            cls._count_methods.add(_method_name)


    def get_interactions(
            self,
            direction = None,
            effect = None,
            resources = None,
            data_model = None,
            interaction_type = None,
            via = None,
            references = None,
            entity_type = None,
            source_entity_type = None,
            target_entity_type = None,
            datasets = None,
        ):
        """
        Returns one or two tuples of the interacting partners: one if only
        one direction, two if both directions match the query criteria.
        The tuple will be empty if no evidence matches the criteria.

        :arg NontType,bool,tuple direction:
            If `None` both undirected and directed, if `True` only directed,
            if a *tuple* of entities only the interactions with that specific
            direction will be considered. Unless you set this parameter to
            `True` this method will return both directions if one or more
            undirected resources present.
            If `False`, only the undirected interactions will be considered,
            and if any resource annotates this interaction as undirected
            both directions will be returned. However the
            ``count_interactions_undirected`` method will return `1`
            in this case.
        :arg NoneType,bool,str effect:
            If `None` also interactions without effect, if `True` only
            the ones with any effect, if a string naming an effect only the
            interactions with that specific effect will be considered.
        :arg NontType,str,set resources:
            Optionally limit the query to one or more resources.
        :arg NontType,str,set data_model:
            Optionally limit the query to one or more data models e.g.
            `activity_flow`.
        :arg NontType,str,set interaction_type:
            Optionally limit the query to one or more interaction types
            e.g. `PPI`.
        :arg NontType,bool,str,set via:
            Optionally limit the query to certain secondary databases or
            if `False` consider only data from primary databases.
        :arg str entity_type:
            Molecule type for both of the entities.
        :arg str source_entity_type:
            Molecule type for the source entity.
        :arg str target_entity_type:
            Molecule type for the target entity.
        """

        effect = self._effect_synonyms(effect)
        direction = (
            self.direction_key(direction)
                if isinstance(direction, tuple) else
            direction
        )

        entity_type = common.to_set(entity_type)
        source_entity_type = common.to_set(source_entity_type) or entity_type
        target_entity_type = common.to_set(target_entity_type) or entity_type

        return tuple(

            # direction key
            _dir

            # possible directions
            for _dir in (self.a_b, self.b_a)

            # conditions by selecting and evaluating evidence collections
            if (
                (
                    not source_entity_type or
                    _dir[0].entity_type in source_entity_type
                ) and
                (
                    not target_entity_type or
                    _dir[1].entity_type in target_entity_type
                )
            )
            and
            (
                self.evaluate_evidences(
                    this_direction = _dir,
                    direction = direction,
                    effect = effect,
                    resources = resources,
                    data_model = data_model,
                    interaction_type = interaction_type,
                    via = via,
                    references = references,
                    datasets = datasets,
                )
            )

        )


    def evaluate_evidences(
            self,
            this_direction,
            direction = None,
            effect = None,
            resources = None,
            data_model = None,
            interaction_type = None,
            via = None,
            references = None,
            datasets = None,
        ):
        """
        Selects the evidence collections matching the direction and effect
        criteria and then evaluates if any of the evidences in these
        collections match the evidence criteria.
        """

        kwargs = locals()
        _ = kwargs.pop('self')

        return any(self.iter_match_evidences(**kwargs))


    def iter_match_evidences(
            self,
            this_direction,
            direction = None,
            effect = None,
            resources = None,
            data_model = None,
            interaction_type = None,
            via = None,
            references = None,
            datasets = None,
        ):
        """
        Selects the evidence collections matching the direction and effect
        criteria and yields collections matching the evidence criteria.
        """

        for evs in self.iter_evidences(
            this_direction = this_direction,
            direction = direction,
            effect = effect,
        ):

            if evs.match(
                resource = resources,
                data_model = data_model,
                interaction_type = interaction_type,
                via = via,
                references = references,
                datasets = datasets,
            ):

                yield evs


    def iter_evidences(
            self,
            this_direction,
            direction = None,
            effect = None,
        ):
        """
        Selects and yields evidence collections matching the direction and
        effect criteria.
        """

        # evidence keys
        for evs_key in ('undirected', this_direction):

            # evidence dicts
            for this_effect in ('direction', 'positive', 'negative'):

                if (
                    # only undirected
                    (
                        direction == False and
                        evs_key == 'undirected' and
                        this_effect == 'direction'
                    ) or
                    # undirected
                    (
                        direction is None and
                        not effect and
                        this_effect == 'direction'
                    ) or
                    # directed
                    (
                        direction != False and
                        evs_key != 'undirected' and
                        this_effect == 'direction' and
                        not effect and (
                            # any direction
                            direction == True or
                            # specific direction
                            direction == this_direction
                        )
                    ) or
                    # with effect
                    (
                        direction != False and
                        evs_key != 'undirected' and
                        this_effect != 'direction' and (
                            # any effect
                            effect == True or
                            # specific effect
                            effect == this_effect
                        )
                    )
                ):

                    # getting the evidence dict and the key from it
                    yield getattr(self, this_effect)[evs_key]


    def get_interactions_0(self, **kwargs):
        """
        Returns unique interacting pairs without being aware of the direction.
        """

        kwargs['direction'] = None
        kwargs['effect'] = None

        result = self.get_interactions(**kwargs)

        return result[:1] if result else ()


    def get_interactions_directed(self, **kwargs):
        """
        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        if 'direction' not in kwargs or kwargs['direction'] is None:

            kwargs['direction'] = True

        return self.get_interactions(**kwargs)


    def get_interactions_undirected(self, **kwargs):
        """
        Only the undirected interactions will be considered, if any resource
        annotates this interaction as undirected both directions will be
        returned, no matter if certain resources provide direction. However
        the ``count_interactions_undirected`` method will return `1` in this
        case.

        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        kwargs['direction'] = False

        return self.get_interactions(**kwargs)


    def get_interactions_undirected_0(self, **kwargs):
        """
        Only the undirected interactions will be considered, if any resource
        annotates this interaction as undirected the interacting pair as
        a sorted tuple will be returned inside a one element tuple.

        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        undir = self.get_interactions_undirected(**kwargs)

        return undir[:1] if undir else ()


    def get_interactions_non_directed(self, **kwargs):
        """
        Only the undirected interactions will be considered, if any resource
        annotates this interaction as undirected both directions will be
        returned, but only if no resource provide direction. However
        the ``count_interactions_non_directed`` method will return `1` in
        this case.

        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        kwargs['direction'] = True

        return (
            self.get_interactions_undirected(**kwargs)
                if not self.get_interactions(**kwargs) else
            ()
        )


    def get_interactions_non_directed_0(self, **kwargs):
        """
        Only the undirected interactions will be considered, if any resource
        annotates this interaction as undirected and none as directed, the
        interacting pair as a sorted tuple will be returned inside a one
        element tuple.

        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        nondir = self.get_interactions_non_directed(**kwargs)

        return nondir[:1] if nondir else ()


    def get_interactions_signed(self, **kwargs):
        """
        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        if 'effect' not in kwargs or kwargs['effect'] is None:

            kwargs['effect'] = True

        return self.get_interactions(**kwargs)


    def get_interactions_positive(self, **kwargs):
        """
        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        kwargs['effect'] = 'positive'

        return self.get_interactions(**kwargs)


    def get_interactions_negative(self, **kwargs):
        """
        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        kwargs['effect'] = 'negative'

        return self.get_interactions(**kwargs)


    def get_interactions_mutual(self, **kwargs):
        """
        Note: undirected interactions does not count as mutual but only
        interactions with explicit direction information for both directions.

        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        if 'direction' not in kwargs or kwargs['direction'] is None:

            kwargs['direction'] = True

        interactions = self.get_interactions(**kwargs)

        return interactions if len(interactions) == 2 else ()


    def is_mutual(self, **kwargs):
        """
        Note: undirected interactions does not count as mutual but only
        interactions with explicit direction information for both directions.

        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        return bool(self.get_interactions_mutual(**kwargs))


    def count_interactions_mutual(self, **kwargs):
        """
        Note: undirected interactions does not count as mutual but only
        interactions with explicit direction information for both directions.

        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        return int(self.is_mutual(**kwargs))


    def count_interactions_undirected(self, **kwargs):
        """
        Returns `True` if any resource annotates this interaction without
        direction.

        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        return bool(self.get_interactions_undirected(**kwargs))


    def count_interactions_non_directed(self, **kwargs):
        """
        Returns `True` if any resource annotates this interaction without
        and no resource with direction.

        Args
            kwargs:
                See the docs of method ``get_interactions``.
        """

        return bool(self.get_interactions_non_directed(**kwargs))


    def get_degrees(
            self,
            mode: str,
            direction = None,
            effect = None,
            resources = None,
            data_model = None,
            interaction_type = None,
            via = None,
            references = None,
        ):
        """
        Returns a *set* of nodes with the connections matching the direction,
        effect and evidence criteria. E.g. if the query concerns the incoming
        degrees with positive effect and the matching evidences show A
        activates B, but not the other way around, only "B" will be returned.

        Args
            mode:
                The type of degrees to be considered. Three possible values are
                ``'IN'``, `'OUT'`` and ``'ALL'`` for incoming, outgoing and all
                connections, respectively. If the ``direction`` is ``False``
                the only possible mode is ``ALL``. If the ``direction`` is
                ``None`` and also directed evidence(s) match the criteria these
                will overwrite the undirected evidences and only the directed
                result will be returned.
        """

        kwargs = locals()
        _ = kwargs.pop('self')
        mode = kwargs.pop('mode')

        idx = {
            'ALL': (0, 2),
            'OUT': (0, 1),
            'IN':  (1, 2),
        }

        if direction == False:

            mode = 'ALL'

        if direction is None and not effect:

            _ = kwargs.pop('direction')

            return (
                self.get_degrees(mode = mode, direction = True, **kwargs) or
                self.get_degrees(mode = mode, direction = False, **kwargs)
            )

        result = set()

        node_pairs = self.get_interactions(**kwargs)

        for pair in node_pairs:

            result.update(
                pair[
                    idx[mode][0]:
                    idx[mode][1]
                ]
            )

        return result


    def get_curation_effort(self, **kwargs):

        return tuple(
            (self.a, self.b, res, ref)
            for res, refs in
            iteritems(self.references_by_resource(**kwargs))
            for ref in refs
        )


    @staticmethod
    def _get(self, method, **kwargs):

        via = kwargs['via'] if 'via' in kwargs else False

        return getattr(
            self.get_evidences(
                **kwargs
            ),
            'get_%s' % method,
        )(via = via)


    @staticmethod
    def _count(method):

        @functools.wraps(method)
        def count_method(*args, **kwargs):

            return len(method(*args, **kwargs))

        return count_method


    @staticmethod
    def _by(method, by = 'resources'):

        by = (by,) if isinstance(by, str) else by

        @functools.wraps(method)
        def by_method(*args, **kwargs):

            self = args[0]
            name_keys = kwargs.pop('name_keys', True)

            for _by in by:

                _ = kwargs.pop(_by, None)

            levels_methods = (
                'get_%s%ss' % (
                    _by[:-1] if _by in {'resources', 'references'} else _by,
                    '_name' if _by == 'resources' and name_keys else ''
                )
                for _by in by
            )

            levels = list(itertools.product(*(
                getattr(self, levels_method)(**kwargs)
                for levels_method in levels_methods
            )))

            result = dict(
                (
                    _levels if len(_levels) > 1 else _levels[0],
                    method(
                        *args,
                        **dict(zip(by, _levels)),
                        **kwargs
                    )
                )
                for _levels in levels
            )

            return dict((k, v) for k, v in iteritems(result) if v)

        return by_method


    @classmethod
    def _by_resource(cls, method):

        return cls._by(method, by = 'resources')


    @classmethod
    def _by_data_model(cls, method):

        return cls._by(method, by = 'data_model')


    @classmethod
    def _by_interaction_type_and_data_model(cls, method):

        return cls._by(method, by = ('interaction_type', 'data_model'))


    @classmethod
    def _by_interaction_type_and_data_model_and_resource(cls, method):

        return cls._by(
            method,
            by = ('interaction_type', 'data_model', 'resources'),
        )


    @classmethod
    def _by_interaction_type(cls, method):

        return cls._by(method, by = 'interaction_type')


    @classmethod
    def _by_reference(cls, method):

        return cls._by(method, by = 'references')


    @classmethod
    def _generate_get_methods(cls):

        def _create_get_method(method):

            @functools.wraps(method)
            def _get_method(*args, **kwargs):

                return cls._get(self = args[0], method = method, **kwargs)

            return _get_method

        for _get in cls._get_methods_autogen:

            method_name = 'get_%s' % _get

            signature = cls._update_get_method_signature(method_name)

            cls._add_method(
                method_name = method_name,
                method = _create_get_method(_get),
                # this is not always correct, to be fixed later
                signature = signature,
                doc = (
                    'Retrieves %s matching the criteria.' % (
                        _get.replace('_', ' ')
                    )
                ),
            )


    @classmethod
    def _generate_degree_methods(cls):

        def _create_degree_method(mode, direction, effect):

            wrap_args = mode, direction, effect

            @functools.wraps(wrap_args)
            def _degree_method(*args, **kwargs):

                mode, direction, effect = wrap_args

                kwargs['direction'] = direction
                kwargs['effect'] = effect

                return cls.get_degrees(self = args[0], mode = mode, **kwargs)

            return _degree_method

        for mode, (dir_label, dir_args) in itertools.product(
            cls._degree_modes,
            iteritems(cls._degree_directions)
        ):

            if dir_label in {'undirected', 'non_directed'} and mode != 'ALL':

                continue

            method_name = 'degrees_%s%s' % (
                dir_label,
                '_%s' % mode.lower() if mode != 'ALL' else ''
            )

            cls._count_methods.add(method_name)
            cls._get_methods.add(method_name)

            method = _create_degree_method(mode, *dir_args)
            _method_name = 'get_%s' % method_name

            cls._add_method(
                method_name = _method_name,
                method = method,
                signature = ['mode'] + cls._get_method_signature,
                doc = cls.get_degrees.__doc__,
            )


    @classmethod
    def _generate_count_methods(cls):

        for _get in cls._count_methods:

            _get_method = getattr(cls, 'get_%s' % _get)

            method_name = 'count_%s' % _get

            signature = cls._update_get_method_signature(method_name)

            cls._add_method(
                method_name = method_name,
                method = cls._count(_get_method),
                signature = cls._get_method_signature,
                doc = _get_method.__doc__,
            )


    @classmethod
    def _update_get_method_signature(cls, method_name):

        signature = cls._get_method_signature

        if 'interaction' in method_name:

            signature = signature + [
                ('entity_type', None),
                ('source_entity_type', None),
                ('target_entity_type', None),
            ]

        return signature


    @classmethod
    def _generate_by_methods(cls):

        for _get, _by in itertools.product(
            cls._get_methods,
            cls._by_methods,
        ):

            _get_method = getattr(cls, 'get_%s' % _get)
            method_name = '%s_by_%s' % (_get, _by)
            method = getattr(cls, '_by_%s' % _by)(_get_method)

            cls._add_method(
                method_name = method_name,
                method = method,
                signature = cls._get_method_signature,
                doc = _get_method.__doc__,
            )


    @classmethod
    def _add_method(cls, method_name, method, signature = None, doc = None):

        common.add_method(
            cls,
            method_name,
            method,
            signature = signature,
            doc = doc,
        )


    def generate_df_records(
            self,
            by_source: bool = False,
            with_references: bool = False,
        ):
        """
        Yields interaction records. It is a generator because one edge can
        be represented by one or more records depending on the signs and
        directions and other parameters

        Args
            by_source:
                Yield separate records by resources. This way the node pairs
                will be redundant and you need to group later if you want
                unique interacting pairs. By default is ``False`` because for
                most applications unique interactions are preferred.
                If ``False`` the *refrences* field will still be present
                but with ``None`` values.
            with_references:
                Include the literature references. By default is ``False``
                because you rarely need these and they increase the data size
                significantly.
        """

        def source_add_via(source, via):

            return '%s%s' % (source, '_%s' % via if via else '')


        def iter_sources(evs):

            sources = evs.get_resource_names_via()

            if by_source:

                for source, via in sources:

                    refs = (
                        {
                            ref.pmid
                            for ref in
                            self.get_references(
                                resources = source,
                                interaction_type = interaction_type,
                                via = via,
                            )
                        }
                            if with_references else
                        None
                    )

                    _source = source_add_via(source, via)

                    yield _source, refs

            else:

                _sources = {
                    source_add_via(source, via)
                    for source, via in sources
                }

                refs = (
                    {
                        ref.pmid
                        for ref in
                        self.get_references(
                            resources = {s[0] for s in sources},
                            interaction_type = interaction_type,
                        )
                    }
                        if with_references else
                    None
                )

                if _sources:

                    yield _sources, refs


        for interaction_type in self.get_interaction_types():

            dmodels = (
                self.get_data_models(interaction_type = interaction_type)
            )
            dmodels = dmodels if by_source else (dmodels,)

            for data_model in dmodels:

                evs_undirected = self.get_evidences(
                    direction = 'undirected',
                    interaction_type = interaction_type,
                    data_model = data_model,
                )

                for _dir in (self.a_b, self.b_a):

                    evs_dir = self.get_evidences(
                        direction = _dir,
                        interaction_type = interaction_type,
                        data_model = data_model,
                    )
                    evs_without_sign = evs_dir.__copy__()

                    for _effect, effect in zip(
                        (1, -1),
                        ('positive', 'negative')
                    ):

                        evs_sign = self.get_evidences(
                            direction = _dir,
                            effect = effect,
                            interaction_type = interaction_type,
                            data_model = data_model,
                        )
                        # to make sure we keep all references:
                        evs_sign += evs_sign.intersection(evs_dir)
                        evs_without_sign -= evs_sign
                        evs_undirected -= evs_sign

                        for sources, refs in iter_sources(evs_sign):

                            yield InteractionDataFrameRecord(
                                id_a = _dir[0].identifier,
                                id_b = _dir[1].identifier,
                                type_a = _dir[0].entity_type,
                                type_b = _dir[1].entity_type,
                                directed = True,
                                effect = _effect,
                                type = interaction_type,
                                dmodel = data_model,
                                sources = sources,
                                references = refs,
                            )

                    if evs_without_sign:

                        evs_undirected -= evs_without_sign

                        for sources, refs in iter_sources(evs_without_sign):

                            yield InteractionDataFrameRecord(
                                id_a = _dir[0].identifier,
                                id_b = _dir[1].identifier,
                                type_a = _dir[0].entity_type,
                                type_b = _dir[1].entity_type,
                                directed = True,
                                effect = 0,
                                type = interaction_type,
                                dmodel = data_model,
                                sources = sources,
                                references = refs,
                            )

                if evs_undirected:

                    for sources, refs in iter_sources(evs_undirected):

                        yield InteractionDataFrameRecord(
                            id_a = self.a.identifier,
                            id_b = self.b.identifier,
                            type_a = self.a.entity_type,
                            type_b = self.b.entity_type,
                            directed = False,
                            effect = 0,
                            type = interaction_type,
                            dmodel = data_model,
                            sources = sources,
                            references = refs,
                        )


    def asdict(self, direction: tuple) -> dict:
        """
        Dictionary representation of the evidences.
        """

        direction = self.a_b if self.a_b == direction else self.b_a

        return {
            'id_a': str(direction[0].identifier),
            'id_b': str(direction[1].identifier),
            'positive': self.positive[direction].asdict(),
            'negative': self.negative[direction].asdict(),
            'directed': self.unknown_effect[direction].asdict(),
            'undirected': self.direction['undirected'].asdict(),
        }


    def serialize_attrs(self, direction: tuple | str | None = None) -> str:
        """
        Serialize the resource specific attributes into a JSON string.
        """

        evs = self.evidences if not direction else self.direction[direction]

        return evs.serialize_attrs()


    def serialize_evidences(self, direction: tuple) -> str:
        """
        Serialize the evidences into a JSON string.
        """

        return self._serialize(self.asdict(direction = direction))


    def _get_attr(self, resource, key, direction):

        if resource in self.direction[direction]:

            return self.direction[direction][resource][key]


    def get_attr(
            self,
            resource: str,
            key: str,
            direction: tuple | str | None = None,
        ):
        """
        Extracts the values of one specific attribute.

        Args
            resource:
                Name of the resource.
            key:
                Name of the attribute.
            direction:
                Direction(s) to consider, either a tuple of entities or
                entity names, or the string `undirected`.

        Returns
            Depends on the arguments. The value of the attribute if direction
            is defined. Otherwise a dict with the value of the attribute for
            each direction. The value of the attribute is `None` if the
            resource or the attribute does not belong to this interaction.
        """

        if direction:

            return self._get_attr(resource, key, direction)

        else:

            return dict(
                (
                    d,
                    self._get_attr(resource, key, d)
                )
                for d in self.direction.keys()
            )


    def dorothea_levels(self, direction: str | tuple | None = None):
        """
        Retrieves the DoRothEA confidence levels.

        Args
            direction:
                Direction(s) to consider, either a tuple of entities or
                entity names, or the string `undirected`.

        Returns
            List of unique single letter strings representing the five
            confidence levels (A-E).
        """

        directions = (direction,) if direction else (self.a_b, self.b_a)

        return sorted(
            {
                level
                for direction in directions
                for level in (
                    self._get_attr('DoRothEA', 'level', direction) or
                    ()
                )
            }
        )


    def dorothea_level(
            self,
            direction: str | tuple,
        ) -> Literal['A', 'B', 'C', 'D', 'E']:
        """
        DoRothEA confidence level for one direction as a single letter.

        Some interactions might have multiple levels due to the ambiguous
        nature of translating gene symbols to UniProt IDs. Here we take the
        highest level and drop the rest. For interactions without DoRothEA
        levels None is returned.
        """

        levels = self.dorothea_levels(direction = direction)

        return common.first(levels)


Interaction._generate_entity_methods()
Interaction._generate_get_methods()
Interaction._generate_degree_methods()
Interaction._generate_count_methods()
Interaction._generate_by_methods()
