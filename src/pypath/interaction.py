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

import pypath.evidence as pypath_evidence
import pypath.session_mod as session_mod

_logger = session_mod.Logger(name = 'interaction')
_log = _logger._log


class InteractionAttributes(object):
    
    
    def __init__(
            self,
            id_a,
            id_b,
            id_type_a,
            id_type_b,
        ):
        
        self.nodes = tuple(sorted((id_a, id_b)))
        
        self.id_type_a, self.id_type_b = (
            (id_type_a, id_type_b)
                if (id_a, id_b) == self.nodes else
            (id_type_b, id_type_a)
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
                if isinstance(evidence, pypath_evidence.Evidence) else
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
