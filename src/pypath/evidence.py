#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  (Planned for) centrally handling cache for all databases/resources.
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

from future.utils import iteritems

import copy

import pypath.refs as refs
import pypath.common as common
import pypath.session_mod as session_mod

_logger = session_mod.Logger(name = 'evidence')
_log = _logger._log


class Evidence(object):
    
    
    def __init__(self, resource, references = None):
        
        self.resource = resource
        self.references = self._process_references(references)
    
    
    @staticmethod
    def _process_references(references):
        
        references = common.to_set(references)
        
        return (
            set(
                (
                    refs.Reference(ref)
                        if not isinstance(ref, refs.Reference) else
                    ref
                )
                for ref in references
            )
        )
    
    
    def __hash__(self):
        
        return self.resource.__hash__()
    
    
    def __eq__(self, other):
        
        return (
            self.resource == other or
            (
                hasattr(other, 'resource') and
                self.resource == other.resource
            )
        )
    
    
    def __iadd__(self, other):
        """
        This will ignore if the other evidence is from different resource:
        still better than attributing wrong references to a resource.
        """
        
        if self == other:
            
            self.references.update(other.references)
            
        else:
            
            _log(
                'Warning: attempt to merge evidences from different '
                'resources. Ignoring the second evidence.'
            )
        
        return self
    
    
    def __add__(self, other):
        
        return Evidence(
            resource = self.resource,
            references = self.references | other.references,
        )
    
    
    @property
    def key(self):
        
        return self.resource.key
    
    
    def merge(self, other):
        """
        Merges two evidences. Returns set of either one or two evidences
        depending on whether the two evidences are from the same resource.
        """
        
        if self == other:
            
            self += other
            return {self}
            
        else:
            
            return {self, other}
    
    
    def __repr__(self):
        
        return '<Evidence %s (%s%u references)>' % (
            self.resource.name,
            'via %s,' % self.resource.via if self.resource.via else '',
            len(self.references),
        )
    
    
    def __copy__(self):
        
        return Evidence(
            resource = self.resource,
            references = copy.copy(self.references),
        )
    
    
    def __contains__(self, other):
        """
        :arg str,tuple,Reference other:
            Either a reference or a database name, or a tuple of a database
            name and an interaction type or a tuple of a database, interaction
            type and a primary database (or None if the query should be
            limited only to primary databases).
        """
        
        return self._contains(self, other)
    
    
    def contains_database(self, database):
        
        return self.resource.name == database
    
    
    def contains_reference(self, reference):
        
        return reference in self.references
    
    
    def has_database_via(self, database, via):
        
        return (
            self.resource.name == database and
            self.via == via
        )
    
    
    def has_interaction_type(
            self,
            interaction_type,
            database = None,
            via = False,
        ):
        """
        If ``via`` is ``False`` then it will be ignored, otherwise if ``None``
        only primary resources are considered.
        """
        
        return (
            self.resource.interaction_type == interaction_type and
            (
                not database or
                self.resource.name == database
            ) and
            (
                via == False or
                self.resource.via == via
            )
        )
    
    
    @staticmethod
    def _contains(obj, other):
        
        if isinstance(other, refs.Reference):
            
            return obj.contains_reference(other)
        
        # this makes possible to accept a NetworkResource or a
        # NetworkResourceKey:
        if (
            hasattr(other, 'name') and
            hasattr(other, 'interaction_type') and
            hasattr(other, 'via')
        ):
            
            other = (other.name, other.interaction_type, other.via)
        
        other = other if isinstance(other, tuple) else (other,)
        
        return (
            obj.contains_database(other[0]) and
            (
                len(other) == 1 or
                obj.has_interaction_type(other[1], other[0])
            ) and
            (
                len(other) <= 2 or
                obj.has_database_via(other[0], other[2])
            )
        )
    
    
    def has_data_model(self, data_model):
        
        return self.resource.data_model == data_model


class Evidences(object):
    
    
    def __init__(self, evidences = ()):
        
        self.evidences = {}
        self.__iadd__(evidences)
    
    
    def __iadd__(self, other):
        
        other = (
            other
                if hasattr(other, '__iter__') else
            (other,)
                if isinstance(other, Evidence) else
            ()
        )
        
        for ev in other:
            
            if ev.key in self.evidences:
                
                self.evidences[ev.key] = self.evidences[ev.key] + ev
                
            else:
                
                self.evidences[ev.key] = ev.__copy__()
        
        return self
    
    
    def __add__(self, other):
        
        return Evidences(
            (
                self.evidences[key].__copy__()
                    if key not in other.evidences else
                other.evidences[key].__copy__()
                    if key not in self.evidences else
                self.evidences[key] + other.evidences[key]
            )
            for key in
            set(self.evidences.keys()) | set(other.evidences.keys())
        )
    
    
    def __iter__(self):
        
        for ev in self.evidences.values():
            
            yield ev
    
    
    def __repr__(self):
        
        return '<Evidences: %s (%u references)>' % (
            ', '.join(sorted(ev.resource.name for ev in self)),
            sum(len(ev.references) for ev in self),
        )
    
    
    def __copy__(self):
        
        return Evidences((ev.__copy__() for ev in self.evidences))
    
    
    def __bool__(self):
        
        return bool(len(self.evidences))
    
    
    def __contains__(self, other):
        """
        :arg str,tuple,Reference other:
            Either a reference or a database name, or a tuple of a database
            name and an interaction type or a tuple of a database, interaction
            type and a primary database (or None if the query should be
            limited only to primary databases).
        """
        
        return Evidence._contains(self, other)
    
    
    def contains_database(self, database):
        
        return any(ev.resource.name == database for ev in self)
    
    
    def contains_reference(self, reference):
        
        return any(reference in ev.references for ev in self)
    
    
    def has_database_via(self, database, via):
        
        return any(
            ev.has_database_via(database, via)
            for ev in self
        )
    
    
    def has_interaction_type(
            self,
            interaction_type,
            database = None,
            via = False,
        ):
        """
        If ``via`` is ``False`` then it will be ignored, otherwise if ``None``
        only primary resources are considered.
        """
        
        return any(
            ev.has_interaction_type(interaction_type, database, via)
            for ev in self
        )
    
    
    def has_data_model(self, data_model):
        
        return any(ev.has_data_model(data_model) for ev in self)
    
    
    def get_resources(self):
        
        return {ev.resource for ev in self}
    
    
    def get_resource_names(self):
        
        return {ev.resource.name for ev in self}
