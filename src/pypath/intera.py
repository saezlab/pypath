#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
# this module provides classes to represent and handle
# structural details of protein interactions
# i.e. residues, post-translational modifications,
# short motifs, domains, domain-motif ands
# domain-motif interactions, binding interfaces
#

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import re
import sys
from collections import Counter

# from pypath:
import pypath.common as common

__all__ = [
    'Residue', 'Ptm', 'Motif', 'Domain', 'DomainDomain', 'DomainMotif',
    'Interface'
]

if 'unicode' not in __builtins__:
    unicode = str


class Residue(object):
    
    def __init__(self,
                 number,
                 name,
                 identifier,
                 id_type='uniprot',
                 isoform=1,
                 mutated=False,
                 seq=None):
        non_digit = re.compile(r'[^\d.-]+')
        self.name = name
        self.number = number if type(number) not in [str, unicode] \
            else int(non_digit.sub('', number))
        self.protein = identifier
        self.id_type = id_type
        self.mutated = mutated
        self.seq = seq
        self.isoform = isoform if type(isoform) is int \
            else int(non_digit.sub('', isoform))

    def __hash__(self):
        return hash((self.number, self.name, self.protein))

    def __eq__(self, other):
        if type(other) is Residue and self.protein == other.protein and \
                self.number == other.number and self.name == other.name:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return 'Residue %s-%u in protein %s-%u%s\n' % (
            self.name, self.number, self.protein, self.isoform, ' (mutated)'
            if self.mutated else '')

    def __repr__(self):
        return self.__str__()

    def serialize(self):
        return '%s:%u' % (self.name, self.number)

    def in_isoform(self, isoform, seq=None):
        seq = seq or self.seq
        if seq and seq.has_isoform(isoform):
            if seq.get(self.number, isoform=isoform) == self.name:
                res = Residue(self.number, self.name, self.protein,
                              self.id_type, isoform, self.mutated)
                return res
        return None


class Mutation(object):
    
    def __init__(self, original, mutated, sample, properties={}):
        if original.protein == mutated.protein and \
                original.number == mutated.number:
            self.original = original
            self.mutated = mutated
            self.protein = original.protein
            self.number = original.number
            self.isoform = original.isoform \
                if original.isoform == mutated.isoform else None
            self.sample = sample
            self.prop = properties
        else:
            self.protein = None
            sys.stdour.write(
                '\t:: Warning: Original and mutated residues should'
                'be within the same protein.\n'
                '\t     Here: original: %s-%u %s%u, '
                'mutated: %s-%u %s%u\n' %
                (original.protein, original.isoform, original.name,
                 original.number, mutated.protein, mutated.isoform,
                 mutated.name, mutated.number))

    def __hash__(self):
        return hash((self.original, self.mutated))

    def __str__(self):
        return '%s%u%s' % (self.original.name, self.number, self.mutated.name)

    def __eq__(self, other):
        if self.protein == other.protein and self.original == other.original and \
                self.mutated == other.mutated:
            return True
        return False

    def __contains__(self, other):
        if isinstance(other, Residue):
            if self.protein == other.protein and self.number == other.number and \
                    self.original.name == other.name:
                return True
        elif isinstance(other, Ptm):
            if other.residue == self.original or self.original in other.motif:
                return True
        elif isinstance(other, Domain):
            if self.original in other:
                return True
        elif isinstance(other, DomainMotif):
            if self.original in other.ptm or self.original in other.domain:
                return True
        else:
            return False


class Ptm(object):
    
    def __init__(self,
                 protein,
                 id_type='uniprot',
                 typ='unknown',
                 motif=None,
                 residue=None,
                 source=None,
                 isoform=1,
                 seq=None):
        self.non_digit = re.compile(r'[^\d.-]+')
        self.protein = protein
        self.id_type = id_type
        self.typ = typ
        self.seq = seq
        self.motif = motif
        self.residue = residue
        self.isoform = isoform if type(isoform) is int \
            else int(self.non_digit.sub('', isoform))
        self.sources = set([])
        self.add_source(source)
        self.isoforms = set([])
        self.add_isoform(isoform)

    def __hash__(self):
        return hash((self.residue, self.typ))

    def __str__(self):
        return '%s in protein %s-%u\n    '\
            'Motif: %s\n%s' % (
                'Domain-motif interaction' if self.typ == 'unknown' and self.residue is None
                else 'PTM: %s' % self.typ,
                self.protein, self.isoform,
                'unknown' if self.motif is None else self.motif.__str__(),
                '' if self.residue is None else '\n    Residue: %s' % self.residue.__str__()
            )

    def __eq__(self, other):
        if type(other) == Ptm and \
            self.protein == other.protein and \
            (self.residue == other.residue or
                (
                    (self.residue is None or other.residue is None) and
                    (self.motif is None or other.motif is None or
                        self.motif == other.motif
                     )
                )
             ) and \
            (self.typ == other.typ or
             self.typ is None or other.typ is None):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, other):
        if isinstance(other, Residue):
            if self.residue is not None:
                return other == self.residue
            elif self.motif is not None:
                return other in self.motif
            else:
                return False
        if isinstance(other, Motif):
            return other in self.motif
        elif other == self.protein:
            return True
        elif isinstance(other, Mutation):
            if other.original == self.residue or \
                    other.original in self.motif:
                return True
        else:
            return False
    
    def __deepcopy__(self, memo):
        new = type(self)(self.protein,
                         id_type = self.id_type,
                         typ = self.typ,
                         motif = self.motif,
                         residue = self.residue,
                         isoform = self.isoform)
        new.add_isoform(self.isoforms)
        return new

    def add_source(self, source):
        if source is None:
            return None
        elif type(source) in common.charTypes:
            self._add_source(source)
        else:
            for s in source:
                self._add_source(s)

    def _add_source(self, source):
        self.sources.add(source)

    def serialize(self):
        return '%s-%u:%s:%s:%s:%s:%u' % (
            self.protein, self.isoform, self.typ, ','.join(self.sources),
            ':::0-0' if self.motif is None else self.motif.serialize(), ''
            if self.residue is None else self.residue.name, 0
            if self.residue is None else self.residue.number)

    def print_residue(self):
        return '%s-%u:%s:%u' % (
            self.protein, self.isoform, ''
            if self.residue is None else self.residue.name, 0
            if self.residue is None else self.residue.number)

    def merge(self, other):
        if self == other:
            self.add_source(other.sources)
            self.motif = self.motif if other.motif is None \
                else other.motif if self.motif is None \
                else self.motif.merge(other.motif)
            if (self.typ == 'unknown' or len(self.typ) == 3) and \
                    other.typ != 'unknown':
                self.typ = other.typ
            self.isoform = min(self.isoform, other.isoform)
            self.isoforms = other.isoforms | self.isoforms

    def add_isoform(self, isoform):
        isoform = set([isoform]) if type(isoform) is int \
            else isoform if type(isoform) is set \
            else set([int(self.non_digit.sub('', isoform))])
        self.isoforms = self.isoforms | isoform

    def get_isoforms(self, seq=None):
        result = []
        seq = seq or self.seq
        if seq:
            for isoform in seq.get_isoforms():
                ptm = self.in_isoform(isoform, seq)
                if ptm:
                    result.append(ptm)
        return result

    def in_isoform(self, isoform, seq=None):
        seq = seq or self.seq

        if seq and seq.has_isoform(isoform):

            if seq.get(self.residue.number,
                       isoform=isoform) == self.residue.name:

                res = self.residue.in_isoform(isoform, seq=seq)
                mot = self.motif.in_isoform(isoform, seq=seq)
                ptm = Ptm(self.protein, self.id_type, self.typ, mot, res,
                          self.sources, isoform, seq)

                return ptm
        return None


class Motif(object):
    
    def __init__(self,
                 protein,
                 start,
                 end,
                 id_type='uniprot',
                 regex=None,
                 instance=None,
                 isoform=1,
                 motif_name=None,
                 prob=None,
                 elm=None,
                 description=None,
                 seq=None,
                 source=None):
        non_digit = re.compile(r'[^\d.-]+')
        self.protein = protein
        self.id_type = id_type
        self.seq = seq
        self.isoform = isoform if type(isoform) is int \
            else int(non_digit.sub('', isoform))
        self.start = start if type(start) not in [str, unicode] \
            else int(non_digit.sub('', start))
        self.end = end if type(end) not in [str, unicode] \
            else int(non_digit.sub('', end))
        self.regex = None if regex is None else re.compile(regex)
        self.instance = instance
        self.motif_name = motif_name
        self.prob = prob
        self.elm = elm
        self.description = description
        self.sources = set([])
        if source is not None:
            self.add_source(source)

    def __hash__(self):
        return hash((self.protein, self.start, self.end))

    def __eq__(self, other):
        if other.protein == self.protein and \
                other.start == self.start and \
                other.end == self.end:
            return True
        else:
            return False

    def __contains__(self, other):
        if isinstance(other, Residue):
            if other.protein == self.protein and \
                    other.number >= self.start and \
                    other.number <= self.end:
                return True
        elif isinstance(other, Mutation):
            if other.original in self:
                return True
        elif other == self.protein or \
                other == self.instance or \
                other == self.motif_name:
            return True
        return False

    def add_source(self, source):
        if source is None:
            return None
        elif type(source) in common.charTypes:
            self._add_source(source)
        else:
            for s in source:
                self._add_source(s)

    def _add_source(self, source):
        self.sources.add(source)

    def serialize(self):
        return '%s:%s:%u-%u' % ('unknown' if self.motif_name is None else
                                self.motif_name, self.instance, 0
                                if self.start is None else self.start, 0
                                if self.end is None else self.end)

    def print_residues(self):
        return '%s-%u:%u-%u' % (self.protein, self.isoform, 0
                                if self.start is None else self.start, 0
                                if self.end is None else self.end)

    def merge(self, other):
        if self == other:
            self.instance = self.instance or other.instance
            self.regex = self.regex or other.regex
            self.elm = self.elm or other.elm
            self.prob = self.prob or other.prob
            self.motif_name = self.motif_name or other.motif_name
            self.description = self.description or other.description

    def __str__(self):
        return 'Motif in protein %s-%u:\n'\
            '\tName: %s\n'\
            '\tELM: %s\n'\
            '\tRange: %u-%u\n'\
            '\tRegex: %s\n'\
            '\tInstance: %s\n' % (self.protein, self.isoform,
                                  'unknown' if self.motif_name is None else self.motif_name,
                                  'unknown' if self.elm is None else self.elm,
                                  0 if self.start is None else self.start,
                                  0 if self.end is None else self.end,
                                  'unknown' if self.regex is None else self.regex.pattern,
                                  'unknown' if self.instance is None else self.instance)

    def in_isoform(self, isoform, seq=None):
        seq = seq or self.seq
        if seq and seq.has_isoform(isoform):
            start, end, reg = seq.get_region(self.start, self.start, self.end)
            mot = Motif(self.protein, start, end, self.id_type, self.refex,
                        reg, isoform, self.motif_name, self.prob, self.elm,
                        self.description, seq)
            return mot
        return None


class Domain(object):
    
    def __init__(self,
                 protein,
                 id_type='uniprot',
                 domain=None,
                 domain_id_type='pfam',
                 start=None,
                 end=None,
                 isoform=1,
                 chains={}):
        non_digit = re.compile(r'[^\d.-]+')
        self.protein = protein
        self.id_type = id_type
        self.domain = domain
        self.domain_id_type = domain_id_type
        self.start = start if type(start) not in [str, unicode] \
            else int(non_digit.sub('', start))
        self.end = end if type(end) not in [str, unicode] \
            else int(non_digit.sub('', end))
        self.isoform = isoform if type(isoform) is int \
            else int(non_digit.sub('', isoform))
        self.pdbs = {}
        for pdb, chain in iteritems(chains):
            self.add_chains(pdb, chain)

    def __hash__(self):
        return hash((self.protein, self.domain))

    def __eq__(self, other):
        if (self.start and self.end and other.start and other.end) is None:
            return False
        flk = min(
            max(
                int(
                    min(self.end - self.start, other.end - other.start) * 0.1),
                10), 30)
        if self.protein == other.protein and \
                self.id_type == other.id_type and \
                self.start is not None and self.end is not None and \
                self.start < other.start + flk and \
                self.start > other.start - flk and \
                self.end < other.end + flk and \
                self.end > other.end - flk:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, other):
        if isinstance(other, Residue):
            if other.protein == self.protein and \
                    other.number >= self.start and \
                    other.number <= self.end:
                return True
            else:
                return False
        elif isinstance(other, Motif):
            if other.protein == self.protein and \
                    other.start < self.end and \
                    other.end <= self.start:
                return True
            else:
                return False
        elif isinstance(other, Ptm):
            return other.residue in self or other.motif in self
        elif other == self.protein or \
                other == self.instance or \
                other == self.motif_name:
            return True
        else:
            return False

    def has_position(self):
        return bool(self.start and self.end)

    def get_position(self):
        return (self.start, self.end)

    def add_chains(self, pdb, chain):
        if pdb not in self.pdbs:
            self.pdbs[pdb] = []
        self.pdbs[pdb] = common.addToList(self.pdbs[pdb], chain)

    def serialize(self):
        return '%s-%u:%s:%s:%u-%u' % (
            self.protein, self.isoform, 'unknown'
            if self.domain is None else self.domain, 0
            if self.start is None else self.start, 0
            if self.end is None else self.end, ','.join([
                '%s.%s' % (pdb, '.'.join(chains))
                for pdb, chains in iteritems(self.pdbs)
            ]))

    def __str__(self):
        return 'Domain in protein %s-%u:\n'\
            '\tName: %s\n'\
            '\tRange: %u-%u\n'\
            '\t3D structures: %s\n' % (self.protein, self.isoform,
                                       'unknown' if self.domain is None else self.domain,
                                       0 if self.start is None else self.start,
                                       0 if self.end is None else self.end,
                                       ', '.join(['%s (chains %s)' % (pdb, ', '.join(chains))
                                                  for pdb, chains in iteritems(self.pdbs)])
                                       )

    def merge(self, other):
        if self == other or (self.start and self.end) is None or \
                (other.start and other.end) is None:
            for pdb, chain in iteritems(other.pdbs):
                self.add_chains(pdb, chain)
            self.domain_id_type = self.domain_id_type or other.domain_id_type
            if self.domain_id_type != 'pfam' and other.domain is not None \
                    or (self.domain is None and other.domain is not None):
                self.domain = other.domain


class DomainDomain(object):
    
    def __init__(self,
                 domain_a,
                 domain_b,
                 pdbs=None,
                 sources=None,
                 refs=None,
                 contact_residues=None):
        self.domains = [domain_a, domain_b]
        self.sources = set([])
        self.refs = set([])
        self.pdbs = set([])
        self.add_sources(sources)
        self.add_refs(refs)
        self.add_pdbs(pdbs)
        '''This can be found from 3DComplexes; floating point 
        numbers show the number of residues in contact. Other
        two numbers in the tuple are the length of domain sequences.'''
        self.contact_residues = contact_residues

    def __hash__(self):
        return hash((self.domain_a, self.domain_b))

    def __eq__(self, other):
        if self.__dict__ == other.__dict__:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, other):
        return other in self.domains[0] or other in self.domains[1]

    def add_sources(self, source):
        if source is None:
            return None
        elif type(source) in common.charTypes:
            self._add_source(source)
        else:
            for s in source:
                self._add_source(s)

    def _add_source(self, source):
        self.sources.add(source)

    def add_refs(self, refs):
        self.refs = common.addToSet(self.refs, refs)

    def add_pdbs(self, pdbs):
        self.pdbs = common.addToSet(self.pdbs, pdbs)

    def serialize(self):
        return '|'.join([
            self.domains[0].serialize(), self.domains[1].serialize(),
            ','.join(self.sources), ','.join(self.refs), ','.join(self.pdbs)
        ])

    # domain1|domain2|sources|references|pdb

    def __str__(self):
        return 'Domain-domain interaction:\n'\
            ' %s %s\n'\
            ' Data sources: %s\n'\
            ' References: %s\n'\
            ' 3D structures: %s\n' % (self.domains[0].__str__(),
                                      self.domains[1].__str__(), ', '.join(
                                          self.sources),
                                      ', '.join(self.refs), ', '.join(self.pdbs))


class DomainMotif(object):
    
    def __init__(self, domain, ptm, sources=None, refs=None, pdbs=None):
        self.ptm = ptm
        self.domain = domain
        self.sources = set([])
        self.refs = set([])
        self.pdbs = set([])
        self.add_sources(sources)
        self.add_refs(refs)
        self.add_pdbs(pdbs)
        self.pnetw_score = None

    def __hash__(self):
        return hash((self.domain, self.ptm))

    def __str__(self):
        return 'Domain-motif interaction:\n'\
            '  %s  %s'\
            '  Data sources: %s\n'\
            '  References: %s\n'\
            '  3D structures: %s\n' % (self.domain.__str__(),
                                       self.ptm.__str__(), ', '.join(self.sources),
                                       ', '.join(self.refs), ', '.join(self.pdbs))

    def __eq__(self, other):
        if isinstance(other, DomainMotif) and \
            self.ptm == other.ptm and \
            (self.domain == other.domain or
             (self.domain.start and self.domain.end) is None or
             (other.domain.start and other.domain.end) is None):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __contains__(self, other):
        if other == self.domain or other == self.ptm:
            return True
        elif other == self.domain.protein or other == self.ptm.protein:
            return True
        else:
            return False

    def get_proteins(self):
        return [self.domain.protein, self.ptm.protein]

    def add_sources(self, source):
        if source is None:
            return None
        elif type(source) in common.charTypes:
            self._add_source(source)
        else:
            for s in source:
                self._add_source(s)

    def _add_source(self, source):
        self.sources.add(source)

    def add_refs(self, refs):
        self.refs = common.addToSet(self.refs, refs)

    def add_pdbs(self, pdbs):
        self.pdbs = common.addToSet(self.pdbs, pdbs)

    def serialize(self):
        return '|'.join([
            self.domain.serialize(), self.ptm.serialize(),
            ','.join(self.sources), ','.join(self.refs), ','.join(self.pdbs)
        ])

    def print_residues(self):
        return '%s-%u:%s:%s' % (self.domain.protein, self.domain.isoform,
                                '%s-%u:' % (self.ptm.protein, self.ptm.isoform)
                                if self.ptm.motif is None else
                                self.ptm.motif.print_residues(),
                                self.ptm.print_residue())

    def merge(self, other):
        if self == other:
            self.domain.merge(other.domain)
            self.ptm.merge(other.ptm)
            self.add_sources(other.sources)
            self.add_refs(other.refs)
            self.add_pdbs(other.pdbs)
            self.pnetw_score = self.pnetw_score or other.pnetw_score
    
    def get_line(self):
        return [
            self.domain.protein,
            self.ptm.protein,
            ';'.join(map(lambda i: '%u' % i, sorted(self.ptm.isoforms))),
            self.ptm.residue.name,
            '%u' % self.ptm.residue.number,
            self.ptm.typ,
            ';'.join(sorted(self.sources)),
            ';'.join(sorted(self.refs))
        ]


class Regulation(object):
    def __init__(self, ptm, source, target, effect, sources=None, refs=None):
        self.ptm = ptm if type(ptm) is list else [ptm]
        self.source = source
        self.target = target
        self.effect = effect
        self.sources = set([])
        self.refs = set([])
        self.add_sources(sources)
        self.add_refs(refs)

    def __hash__(self):
        return hash((self.ptm, self.source, self.target, self.effect))

    def __eq__(self, other):
        if isinstance(other, Regulation) and \
                self.ptm == other.ptm and \
                self.source == other.source and \
                self.target == other.target and \
                self.effect == other.effect:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def add_sources(self, source):
        if source is None:
            return None
        elif type(source) in common.charTypes:
            self._add_source(source)
        else:
            for s in source:
                self._add_source(s)

    def _add_source(self, source):
        self.sources.add(source)

    def add_refs(self, refs):
        self.refs = common.addToSet(self.refs, refs)

    def serialize(self):
        return '|'.join([
            self.effect, self.ptm.serialize(), self.target,
            ','.join(self.sources), ','.join(self.refs)
        ])

    def __str__(self):
        return 'Regulation by PTM:\n'\
            ' PTM on %s %s interaction with %s\n'\
            ' %s \n'\
            ' Data sources: %s\n'\
            ' References: %s\n' % (self.source, self.target,
                                   self.ptm.__str__(), ', '.join(self.sources),
                                   ', '.join(self.refs))

    def merge(self, other):
        if self == other:
            self.ptm.merge(other.ptm)
            self.add_sources(other.sources)
            self.add_refs(other.refs)


class Complex(object):
    def __init__(self, proteins, name, long_name, sources):
        self.synonyms = set([])
        sources = sources if type(sources) is list else [sources]
        sources = set(sources)
        proteins = proteins if type(proteins) is list else list(proteins)
        for key, val in iteritems(locals()):
            setattr(self, key, val)
        self.members = sorted(uniqList(self.proteins))
        self.constitution = Counter(proteins)

    def __hash__(self):
        return hash(tuple(sorted(self.proteins)))

    def __contains__(self, item):
        return item in self.proteins

    def __eq__(self, other):
        return self.proteins == other.proteins


class Interface(object):
    def __init__(self,
                 id_a,
                 id_b,
                 source,
                 id_type='uniprot',
                 pdb=None,
                 css=None,
                 stab_en=None,
                 solv_en=None,
                 area=None,
                 isoform_a=1,
                 isoform_b=1):
        '''
        This class is to store residue level information of 
        protein-protein interfaces.
        '''
        self.source = source
        self.isoform_a = isoform_a if type(isoform_a) is int \
            else int(non_digit.sub('', isoform_a))
        self.isoform_b = isoform_b if type(isoform_b) is int \
            else int(non_digit.sub('', isoform_b))
        self.pdb = pdb
        self.id_a = id_a
        self.id_b = id_b
        self.id_type = id_type
        self.types = ['undefined', 'hbonds', 'sbridges', 'ssbonds', 'covbonds']
        for t in self.types:
            self.__dict__[t] = {id_a: [], id_b: []}
        self.area = area
        self.stab_en = stab_en
        self.solv_en = solv_en
        self.css = css

    def add_residues(self, res_a, res_b, typ='undefined'):
        '''
        Adds one pair of residues of type `typ`, 
        where `res_a` and `res_b` are tuples of 
        residue number in sequence and residue type,
        e.g. (124, 'S') -- (means Serine #124)
        `typ` can be undefined, hbonds, sbridges, ssbonds or covbonds
        '''
        if type(res_a) is not tuple or type(res_b) is not tuple \
                or type(res_a[0]) is not int or type(res_b[0]) is not int \
                or (type(res_a[1]) is not unicode and type(res_a[1]) is not str) \
                or (type(res_b[1]) is not unicode and type(res_b[1]) is not str) \
                or typ not in self.__dict__:
            sys.stdout.write(
                '\tWrong parameters for Interface.add_residues()\n')
        else:
            self.__dict__[typ][self.id_a].append(
                Residue(res_a[0], res_a[1], res_a[2], self.id_type))
            self.__dict__[typ][self.id_b].append(
                Residue(res_b[0], res_b[1], res_b[2], self.id_type))

    def numof_residues(self):
        '''
        Returns the number of residue pairs by bound type
        '''
        nbonds = {}
        for t in self.types:
            nbonds[t] = len(self.__dict__[t][id_a])
        return nbonds

    def bond_types(self):
        '''
        Returns the bond types present in this interface
        '''
        types = []
        for t in self.types:
            if len(self.__dict__[t][id_a]) > 0:
                types.append(t)
        return types

    def get_bonds(self, typ=None, mode=None):
        '''
        Gives a generator to iterate throught bonds in 
        this interface. If no type given, bonds of all types 
        returned. 
        '''
        if typ is None:
            typ = self.types
        if type(typ) is str:
            typ = [typ]
        for t in typ:
            if t in self.__dict__:
                for i in range(0, len(self.__dict__[t])):
                    if mode == 'dict':
                        yield {
                            self.id_a: self.__dict__[t][self.id_a][i],
                            self.id_b: self.id_b,
                            'res_b': self.__dict__[t][self.id_b][i],
                            'type': t
                        }
                    else:
                        yield (self.id_a,) + self.__dict__[t][self.id_a][i] + \
                            (self.id_b,) + self.__dict__[t][self.id_b][i] + \
                            (t,)

    def serialize(self):
        res = []
        for t in self.types:
            res.append('%s:%s+%s' % (t, ','.join(self.__dict__[t][self.id_a]),
                                     ','.join(self.__dict__[t][self.id_b])))
        return '%s-%u:%s-%u:%s:%s:%s' % (self.id_a, self.isoform_a, self.id_b,
                                         self.isoform_b, self.source, self.pdb,
                                         ':'.join(res))

    def __str__(self):
        nbonds = self.numof_residues()
        return 'Molecular interface between %s and %s,\n'\
            'as observed in PDB structure %s\n\n'\
            ' Data source: %s\n'\
            ' Number of residues in contact: %u\n'\
            ' Hydrogene bonds: %u\n'\
            ' Covalent bonds: %u\n'\
            ' Saltbridges: %u\n'\
            ' S-S bonds: %u\n'\
            ' Stable energy: %s\n'\
            ' Solvation energy: %s\n'\
            ' Surface area: %s\n'\
            ' Complexation significance score: %s\n' % (self.id_a, self.id_b,
                                                        self.pdb, self.source,
                                                        sum(nbonds.values()), nbonds[
                                                            'hbonds'], nbonds['covbonds'],
                                                        nbonds['sbridges'], nbonds[
                                                            'ssbonds'],
                                                        'n/a' if self.stab_en is None else str(
                                                            self.stab_en),
                                                        'n/a' if self.solv_en is None else str(
                                                            self.solv_en),
                                                        'n/a' if self.area is None else str(
                                                            self.area),
                                                        'n/a' if self.css is None else str(self.css))
