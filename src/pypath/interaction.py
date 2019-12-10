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
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import importlib as imp
import collections
import operator

import pypath.evidence as pypath_evidence
import pypath.resource as pypath_resource
import pypath.session_mod as session_mod

_logger = session_mod.Logger(name = 'interaction')
_log = _logger._log


class Interaction(object):
    
    
    _key = collections.namedtuple(
        'InteractionKey',
        [
            'id_a',
            'id_b',
            'id_type_a',
            'id_type_b',
            'entity_type_a',
            'entity_type_b',
            'taxon_a',
            'taxon_b',
        ],
    )
    
    
    def __init__(
            self,
            id_a,
            id_b,
            id_type_a = 'uniprot',
            id_type_b = 'uniprot',
            entity_type_a = 'protein',
            entity_type_b = 'protein',
            taxon_a = 9606,
            taxon_b = 9606,
        ):
        
        self.nodes = tuple(sorted((id_a, id_b)))
        self.id_a = self.nodes[0]
        self.id_b = self.nodes[1]
        
        self.id_type_a, self.id_type_b = (
            (id_type_a, id_type_b)
                if (id_a, id_b) == self.nodes else
            (id_type_b, id_type_a)
        )
        self.entity_type_a, self.entity_type_b = (
            (entity_type_a, entity_type_b)
                if (id_a, id_b) == self.nodes else
            (entity_type_b, entity_type_a)
        )
        self.taxon_a, self.taxon_b = (
            (taxon_a, taxon_b)
                if (id_a, id_b) == self.nodes else
            (taxon_b, taxon_a)
        )
        
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
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


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


    def add_evidence(
            self,
            evidence,
            direction = 'undirected',
            effect = 0,
            references = None,
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
        """

        if not self._check_direction_key(direction):
            
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
        
        self.evidences += evidence
        self.direction[direction] += evidence
        
        if direction != 'undirected':
            
            if effect in {1, 'positive', 'stimulation'}:
                
                self.positive[direction] += evidence
            
            elif effect in {-1, 'negative', 'inhibition'}:
                
                self.negative[direction] += evidence
    
    
    def __hash__(self):
        
        return hash(self.key)
    
    
    def __eq__(self, other):
        
        return self.key == other.key
    
    
    @property
    def key(self):
        
        return self._key(
            self.id_a,
            self.id_b,
            self.id_type_a,
            self.id_type_b,
            self.entity_type_a,
            self.entity_type_b,
            self.taxon_a,
            self.taxon_b,
        )


    def __iadd__(self, other):
        
        if self != other:
            
            _log(
                'Attempt to merge interactions with '
                'non matching interaction partners.'
            )
            return self
        
        self._merge_evidences(self, other)
        
        return self
    
    
    def __add__(self, other):
        
        new = self.__copy__()
        
        new += other
        
        return new
    
    
    def __copy__(self):
        
        new = Interaction(*self.key)
        new += self
        
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
    
    
    def __repr__(self):
        
        return '<Interaction: %s %s=%s=%s=%s %s [%s]>' % (
            self.id_a,
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
            self.id_b,
            self.evidences.__repr__().strip('<>'),
        )
    
    
    def __contains__(self, other):
        
        return self.evidences.__contains__(other)
    
    
    def has_data_model(self, data_model):
        
        return self.evidences.has_data_model(data_model)
    
    
    @property
    def data_models(self):
        
        return {
            ev.data_model
            for ev in self.evidences
        }


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

        if self._check_direction_key(direction):
            
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

        if self.check_nodes(query):
            
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
            bool(self.direction[direction])
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
                    getattr(self, effect)
                ) or
                getattr(self, effect) & resources
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

        if not effect:

            return

        if effect in {'positive', 'stimulation', 'stimulatory'}:

            return 'positive'

        if effect in {'negative', 'inhibition', 'inhibitory'}:

            return 'negative'


    def _resources_set(self, resources = None):

        return common.to_set(resources)


    def unset_direction(
            self,
            direction,
            source = None,
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

        if self._check_direction_key(direction):
            
            resource = resource or source
            
            if resource is not None:

                self.sources[direction].remove(
                    resource = resource,
                    interaction_type = interaction_type,
                    via = via,
                )

            else:
                self.direction[direction] = pypath_evidence.Evidences()


    # synonym: old name
    unset_dir = unset_direction


    def is_directed(self):
        """
        Checks if edge has any directionality information.

        :return:
            (*bool*) -- Returns ``True`` if any of the :py:attr:`dirs`
            attribute values is ``True`` (except ``'undirected'``),
            ``False`` otherwise.
        """

        return any(self.direction.values())


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
            self.direction[self.a_b] and self.direction[self.b_a]
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
            direction,
            sign,
            resource = None,
            resource_name = None,
            interaction_type = 'PPI',
            data_model = None,
            **kwargs
        ):
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
            Specifies the type of interaction. Either ``'positive'`` or
            ``'negative'``.
        :arg set resource:
            Contains the name(s) of the source(s) from which the
            information was obtained.
        :arg **kwargs:
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
                    **kwargs,
                )
            )
                if resource_name is not None else
            None
        )
        
        if self._check_direction_key(direction) and evidence is not None:
            
            getattr(self, sign) += evidence


    def get_sign(
            self,
            direction,
            sign = None,
            evidences = False,
            resources = False,
            resource_names = False,
            sources = False,
        ):
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
        :arg bool resources:
            Optional, ``False`` by default. Specifies whether to return
            the resources instead of sign.

        :return:
            (*list*) -- If ``sign=None`` containing [bool] values
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

        if self._check_direction_key(direction):

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
