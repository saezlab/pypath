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


class Evidences(object):
    
    
    def __init__(self, evidences):
        
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
                
                self.evidences[ev.key] += ev
                
            else:
                
                self.evidences[ev.key] = ev
        
        return self
    
    
    def __iter__(self):
        
        for ev in self.evidences.values():
            
            yield ev
    
    
    def __repr__(self):
        
        return '<Evidences: %s (%u references)>' % (
            ', '.join(sorted(ev.resource.name for ev in self)),
            sum(len(ev.references) for ev in self),
        )
