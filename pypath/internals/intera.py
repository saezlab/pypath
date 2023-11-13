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
This module provides classes to represent and handle
structural details of protein interactions
i.e. residues, post-translational modifications,
short motifs, domains, domain-motifs and
domain-motif interactions, binding interfaces.
"""

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import re
import sys
import importlib as imp
import collections
import itertools

# from pypath:
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.utils.mapping as mapping
import pypath.core.evidence as evidence
import pypath.core.entity as entity
import pypath.utils.taxonomy as taxonomy

__all__ = [
    'Residue',
    'Ptm',
    'Motif',
    'Domain',
    'DomainDomain',
    'DomainMotif',
    'Interface',
]

if 'unicode' not in __builtins__:
    unicode = str


COMPLEX_SEP = '_'


class Residue(object):

    def __init__(
            self,
            number,
            name,
            protein,
            id_type = 'uniprot',
            ncbi_tax_id = 9606,
            isoform = 1,
            mutated = False,
            seq = None
        ):

        non_digit = re.compile(r'[^\d.-]+')
        self.name = name
        self.number = (
            number
                if not isinstance(number, str) else
            int(non_digit.sub('', number))
        )

        self.protein = (
            protein
                if hasattr(protein, 'identifier') else
            entity.Entity(
                identifier = protein,
                id_type = id_type,
                taxon = ncbi_tax_id,
            )
        )
        self.mutated = mutated
        self.seq = seq
        self.isoform = (
            isoform
                if type(isoform) is int else
            int(non_digit.sub('', isoform))
        )


    def __hash__(self):

        return hash((self.number, self.name, self.protein))


    def __eq__(self, other):

        return (
            self.protein == other.protein and
            self.number == other.number and
            self.name == other.name
        )

    def __ne__(self, other):

        return not self.__eq__(other)


    def __str__(self):

        return 'Residue %s-%u in protein %s-%u%s\n' % (
            self.name,
            self.number,
            self.protein.identifier,
            self.isoform,
            ' (mutated)' if self.mutated else ''
        )


    def __repr__(self):

        return '<Residue %s-%u:%s%u>' % (
            self.protein.label,
            self.isoform,
            self.name,
            self.number,
        )


    def serialize(self):

        return '%s%u' % (self.name, self.number)


    def in_isoform(self, isoform, seq=None):

        seq = seq or self.seq

        if seq and seq.has_isoform(isoform):

            if seq.get(self.number, isoform=isoform) == self.name:

                res = Residue(
                    number = self.number,
                    name = self.name,
                    protein = self.protein,
                    id_type = self.id_type,
                    isoform = isoform,
                    mutated = self.mutated,
                )
                return res

        return None


class Ptm(object):

    def __init__(
            self,
            protein,
            id_type = 'uniprot',
            ncbi_tax_id = 9606,
            typ = 'unknown',
            motif = None,
            residue = None,
            isoform = 1,
            evidences = None,
            seq = None,
        ):

        self.non_digit = re.compile(r'[^\d.-]+')
        self.protein = (
            protein
                if hasattr(protein, 'identifier') else
            entity.Entity(
                identifier = protein,
                id_type = id_type,
                taxon = ncbi_tax_id,
            )
        )
        self.id_type = id_type
        self.typ = typ.lower()
        self.seq = seq
        self.motif = motif
        self.residue = residue
        self.isoform = (
            isoform
                if type(isoform) is int else
            int(self.non_digit.sub('', isoform))
        )
        self.isoforms = set()
        self.add_isoform(isoform)
        self.evidences = evidence.Evidences()
        self.add_evidences(evidences)


    def __hash__(self):

        return hash((self.residue, self.typ))


    def __str__(self):

        return (
            '%s in protein %s-%u\n    Motif: %s\n%s' % (
                (
                    'Domain-motif interaction'
                        if (
                            self.typ == 'unknown' and
                            self.residue is None
                        ) else
                    'PTM: %s' % self.typ
                ),
                self.protein.label,
                self.isoform,
                (
                    'unknown'
                        if self.motif is None else
                    self.motif.__str__()
                ),
                (
                    ''
                        if self.residue is None else
                    '\n    Residue: %s' % self.residue.__str__()
                ),
            )
        )


    def __repr__(self):

        return '<PTM %s%s>' % (
            (
                self.residue.__repr__().strip('<>')
                    if self.residue else
                self.protein.label
            ),
            ':%s' % self.typ if self.residue else '',
        )


    def __eq__(self, other):

        return (
            isinstance(other, Ptm) and
            self.protein == other.protein and
            (
                self.residue == other.residue or
                (
                    (self.residue is None or other.residue is None) and
                    self.motif == other.motif
                )
            ) and (
                self.typ == other.typ or
                self.typ is None or other.typ is None
            )
        )


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

            return (
                other.original == self.residue or
                other.original in self.motif
            )


    def __deepcopy__(self, memo):

        new = type(self)(
            protein = self.protein,
            id_type = self.id_type,
            typ = self.typ,
            motif = self.motif,
            residue = self.residue,
            isoform = self.isoform,
        )

        new.add_isoform(self.isoforms)

        return new


    def add_evidences(self, evidences):

        self.evidences += evidences


    def serialize(self):

        return '%s-%u:%s:%s:%s:%s:%u' % (
            self.protein,
            self.isoform,
            self.typ,
            ','.join(self.sources),
            ':::0-0' if self.motif is None else self.motif.serialize(),
            '' if self.residue is None else self.residue.name,
            0 if self.residue is None else self.residue.number,
        )

    def print_residue(self):

        return '%s-%u:%s:%u' % (
            self.protein, self.isoform,
            '' if self.residue is None else self.residue.name,
            0 if self.residue is None else self.residue.number,
        )


    def merge(self, other):

        if self == other:

            self.add_evidences(other.evidences)
            self.motif = (
                self.motif
                    if other.motif is None else
                other.motif
                    if self.motif is None else
                self.motif.merge(other.motif)
            )
            if (
                (self.typ == 'unknown' or len(self.typ) == 3) and
                other.typ != 'unknown'
            ):
                self.typ = other.typ

            self.isoform = min(self.isoform, other.isoform)
            self.isoforms = other.isoforms | self.isoforms


    def add_isoform(self, isoform):

        isoform = (
            set([isoform])
                if isinstance(isoform, int) else
            isoform
                if isinstance(isoform, set) else
            {int(self.non_digit.sub('', isoform))}
        )

        self.isoforms = self.isoforms | isoform


    def get_isoforms(self, seq = None):

        result = []
        seq = seq or self.seq

        if seq:

            for isoform in seq.get_isoforms():

                ptm = self.in_isoform(isoform, seq)

                if ptm:

                    result.append(ptm)

        return result


    def in_isoform(self, isoform, seq = None):

        seq = seq or self.seq

        if seq and seq.has_isoform(isoform):

            if (
                seq.get(self.residue.number, isoform = isoform) ==
                self.residue.name
            ):

                res = self.residue.in_isoform(isoform, seq = seq)
                mot = self.motif.in_isoform(isoform, seq = seq)

                ptm = Ptm(
                    protein = self.protein,
                    id_type = self.id_type,
                    typ = self.typ,
                    motif = mot,
                    residue = res,
                    evidences = self.sources,
                    isoform = isoform,
                    seq = seq,
                )

                return ptm


class Motif(object):

    def __init__(
            self,
            protein,
            start,
            end,
            id_type = 'uniprot',
            ncbi_tax_id = 9606,
            regex = None,
            instance = None,
            isoform = 1,
            motif_name = None,
            prob = None,
            elm = None,
            description = None,
            seq = None,
            evidences = None,
        ):

        non_digit = re.compile(r'[^\d.-]+')
        self.protein = (
            protein
                if hasattr(protein, 'identifier') else
            entity.Entity(
                protein,
                id_type = id_type,
                taxon = ncbi_tax_id,
            )
        )
        self.id_type = id_type
        self.seq = seq
        self.isoform = (
            isoform
                if isinstance(isoform, int) else
            int(non_digit.sub('', isoform))
        )
        self.start = (
            start
                if not isinstance(start, str) else
            int(non_digit.sub('', start))
        )
        self.end = (
            end
                if not isinstance(end, str) else
            int(non_digit.sub('', end))
        )
        self.regex = None if regex is None else re.compile(regex)
        self.instance = instance
        self.motif_name = motif_name
        self.prob = prob
        self.elm = elm
        self.description = description
        self.evidences = evidence.Evidences()

        self.add_evidences(evidences)


    def __hash__(self):

        return hash((self.protein, self.start, self.end))


    def __eq__(self, other):

        return (
            other.protein == self.protein and
            other.start == self.start and
            other.end == self.end
        )


    def __contains__(self, other):

        return (
            (
                isinstance(other, Residue) and
                other.protein == self.protein and
                other.number >= self.start and
                other.number <= self.end
            ) or (
                other == self.protein or
                other == self.instance or
                other == self.motif_name
            )
        )


    def add_evidences(self, evidences):

        self.evidences += evidences


    def serialize(self):

        return '%s:%s:%u-%u' % (
            self.motif_name or 'unknown',
            self.instance,
            0 if self.start is None else self.start,
            0 if self.end is None else self.end,
        )


    def print_residues(self):

        return '%s-%u:%u-%u' % (
            self.protein, self.isoform,
            0 if self.start is None else self.start,
            0 if self.end is None else self.end,
        )


    def merge(self, other):

        if self == other:

            self.instance = self.instance or other.instance
            self.regex = self.regex or other.regex
            self.elm = self.elm or other.elm
            self.prob = self.prob or other.prob
            self.motif_name = self.motif_name or other.motif_name
            self.description = self.description or other.description
            self.evidences += other.evidences


    def __str__(self):

        return (
            'Motif in protein %s-%u:\n'
            '\tName: %s\n'
            '\tELM: %s\n'
            '\tRange: %u-%u\n'
            '\tRegex: %s\n'
            '\tInstance: %s\n' % (
                self.protein.label,
                self.isoform,
                self.motif_name or 'unknown',
                self.elm or 'unknown',
                0 if self.start is None else self.start,
                0 if self.end is None else self.end,
                'unknown' if self.regex is None else self.regex.pattern,
                self.instance or 'unknown',
            )
        )


    def __repr__(self):

        rng = self.range_str()

        return '<Motif %sin %s-%u%s>' % (
            '%s ' % self.motif_name if self.motif_name else '',
            self.protein.label,
            self.isoform,
            ' [%s]' % rng if rng else '',
        )


    def range(self):

        return (
            (self.start, self.end)
                if self.start and self.end else
            None
        )


    def range_str(self):

        start_end = self.range()

        return '%s-%s' % start_end if start_end else ''


    def in_isoform(self, isoform, seq = None):

        seq = seq or self.seq

        if seq and seq.has_isoform(isoform):

            start, end, reg = seq.get_region(self.start, self.start, self.end)

            mot = Motif(
                self.protein,
                start,
                end,
                self.id_type,
                self.regex,
                reg,
                isoform,
                self.motif_name,
                self.prob,
                self.elm,
                self.description,
                seq,
            )

            return mot

        return None


class Domain(object):

    def __init__(
        self,
        protein,
        id_type = 'uniprot',
        ncbi_tax_id = 9606,
        domain = None,
        domain_id_type = 'pfam',
        start = None,
        end = None,
        isoform = 1,
        chains = {},
    ):

        non_digit = re.compile(r'[^\d.-]+')
        self.protein = (
            protein
                if hasattr(protein, 'identifier') else
            entity.Entity(
                identifier = protein,
                id_type = id_type,
                taxon = ncbi_tax_id,
            )
        )
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

        if any(
            num is None
            for num in (self.start, self.end, other.start, other.end)
        ):

            return False

        flk = min(
            max(
                int(
                    min(
                        self.end - self.start,
                        other.end - other.start
                    ) * 0.1
                ),
                10
            ),
            30
        )

        return (
            self.protein == other.protein and
            self.id_type == other.id_type and
            self.start is not None and
            self.end is not None and
            self.start < other.start + flk and
            self.start > other.start - flk and
            self.end < other.end + flk and
            self.end > other.end - flk
        )


    def __ne__(self, other):

        return not self.__eq__(other)


    def __contains__(self, other):

        return (
            (
                isinstance(other, Residue) and
                other.protein == self.protein and
                other.number >= self.start and
                other.number <= self.end
            ) or
            (
                isinstance(other, Motif) and
                other.protein == self.protein and
                other.start < self.end and
                other.end <= self.start
            ) or
            (
                isinstance(other, Ptm) and
                (
                    other.residue in self or
                    other.motif in self
                )
            ) or
            (
                other == self.protein or
                other == self.instance or
                other == self.motif_name
            )
        )


    def has_position(self):

        return bool(self.start and self.end)


    def get_position(self):

        return (self.start, self.end)


    def add_chains(self, pdb, chain):

        if pdb not in self.pdbs:

            self.pdbs[pdb] = []

        self.pdbs[pdb] = common.add_to_list(self.pdbs[pdb], chain)


    def serialize(self):

        return '%s-%u:%s:%u-%u:%s' % (
            self.protein,
            self.isoform,
            'unknown' if self.domain is None else self.domain,
            0 if self.start is None else self.start,
            0 if self.end is None else self.end,
            ','.join(
                '%s.%s' % (pdb, '.'.join(chains))
                for pdb, chains in iteritems(self.pdbs)
            )
        )


    def __str__(self):

        return (
            'Domain in protein %s-%u:\n'
            '\tName: %s\n'
            '\tRange: %u-%u\n'
            '\t3D structures: %s\n' % (
                self.protein.label,
                self.isoform,
                self.domain or 'unknown',
                0 if self.start is None else self.start,
                0 if self.end is None else self.end,
                ', '.join(
                    '%s (chains %s)' % (pdb, ', '.join(chains))
                    for pdb, chains in iteritems(self.pdbs)
                )
            )
        )


    def __repr__(self):

        rng = self.range_str()

        return '<Domain %sin %s-%u%s>' % (
            '%s ' % self.domain if self.domain else '',
            self.protein.label,
            self.isoform,
            ' [%s]' % rng if rng else '',
        )


    def range(self):

        return (
            (self.start, self.end)
                if self.start and self.end else
            None
        )


    def range_str(self):

        start_end = self.range()

        return '%s-%s' % start_end if start_end else ''

    def merge(self, other):

        if (
            self == other or
            (self.start and self.end) is None or
            (other.start and other.end) is None
        ):

            for pdb, chain in iteritems(other.pdbs):

                self.add_chains(pdb, chain)

            self.domain_id_type = self.domain_id_type or other.domain_id_type

            if (
                (
                    self.domain_id_type != 'pfam' and
                    other.domain is not None
                ) or
                (
                    self.domain is None and
                    other.domain is not None
                )
            ):

                self.domain = other.domain


class DomainDomain(object):

    def __init__(
            self,
            domain_a,
            domain_b,
            pdbs = None,
            sources = None,
            refs = None,
            contact_residues = None,
        ):

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
        elif type(source) in _const.CHAR_TYPES:
            self._add_source(source)
        else:
            for s in source:
                self._add_source(s)

    def _add_source(self, source):
        self.sources.add(source)

    def add_refs(self, refs):
        self.refs = common.add_to_set(self.refs, refs)

    def add_pdbs(self, pdbs):
        self.pdbs = common.add_to_set(self.pdbs, pdbs)

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
            ' 3D structures: %s\n' % (
                self.domains[0].__str__(),
                self.domains[1].__str__(),
                ', '.join(self.sources),
                ', '.join(self.refs),
                ','.join(self.pdbs)
            )


class DomainMotif(object):


    def __init__(self, domain, ptm, evidences = None, pdbs = None):

        self.ptm = ptm
        self.domain = domain
        self.pdbs = set()
        self.pnetw_score = None

        self.add_pdbs(pdbs)

        self.evidences = evidences or evidence.Evidences()


    def __hash__(self):

        return hash((self.domain, self.ptm))


    def __str__(self):

        return (
            'Domain-motif interaction:\n'
            '  %s  %s'
            '  Data sources: %s\n'
            '  References: %s\n'
            '  3D structures: \n' % (
                self.domain.__str__(),
                self.ptm.__str__(),
                ', '.join(self.evidences.get_resource_names()),
                ', '.join(str(r) for r in self.evidences.get_references()),
            )
        )


    def __repr__(self):

        return '<%s => %s [%s]>' % (
            self.domain.protein.label,
            self.ptm.__repr__().strip('<>').replace('PTM ', ''),
            self.evidences.__repr__().strip('<>')
        )


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


    def key(self):
        """
        Returns a unique key which is a tuple of the proteins, the residue
        and the modification type.
        """

        return (
            self.domain.protein,
            self.ptm.protein,
            self.ptm.residue.name,
            self.ptm.residue.number,
            self.ptm.typ,
        )


    def get_proteins(self):

        return [self.domain.protein, self.ptm.protein]


    def add_pdbs(self, pdbs):

        self.pdbs = common.add_to_set(self.pdbs, pdbs)


    def serialize(self):

        return '|'.join([
            self.domain.serialize(), self.ptm.serialize(),
            ','.join(self.sources), ','.join(self.refs), ','.join(self.pdbs)
        ])


    def print_residues(self):

        return '%s-%u:%s:%s' % (
            self.domain.protein,
            self.domain.isoform,
            '%s-%u:' % (self.ptm.protein, self.ptm.isoform)
                if self.ptm.motif is None else
            self.ptm.motif.print_residues(),
            self.ptm.print_residue(),
        )


    def merge(self, other):

        if self == other:

            self.domain.merge(other.domain)
            self.ptm.merge(other.ptm)
            self.add_evidences(other.evidences)
            self.add_pdbs(other.pdbs)
            self.pnetw_score = self.pnetw_score or other.pnetw_score


    def resources(self, only_primary = False):

        return [
            '%s%s' % (
                res,
                '_%s' % via if via else '',
            )
            for res, via in
            self.evidences.get_resource_names_via(via = None)
            if not only_primary or not via
        ]


    def references(self):

        return self.evidences.get_references()


    def references_by_resource(self, only_primary = True):

        return [
            (
                ev.resource.name,
                ev.resource.via,
                ref,
            )
            for ev in self.evidences
            for ref in ev.references
            if not only_primary or not ev.resource.via
        ]


    def references_by_resource_str(self, only_primary = True):

        return ';'.join(sorted(
            '%s%s:%s' % (
                res,
                '_%s' % via if via else '',
                ref.pmid,
            )
            for res, via, ref
            in self.references_by_resource(only_primary = only_primary)
        ))


    def get_line(self, resources_only_primary = False):
        """
        Returns a list intended to be a row in a data frame of
        enzyme-substrate relationships.

        Elements of the list:
            - enzyme
            - enzyme_genesymbol
            - substrate
            - substrate_genesymbol
            - isoforms
            - residue_type
            - residue_offset
            - modification
            - sources
            - references
            - curation_effort
        """

        return [
            self.domain.protein.identifier,
            self.domain.protein.label,
            self.ptm.protein.identifier,
            self.ptm.protein.label,
            ';'.join(map(lambda i: '%u' % i, sorted(self.ptm.isoforms))),
            self.ptm.residue.name,
            '%u' % self.ptm.residue.number,
            self.ptm.typ,
            ';'.join(sorted(
                self.resources(only_primary = resources_only_primary)
            )),
            self.references_by_resource_str(),
            self.evidences.count_curation_effort(),
        ]


    def add_evidences(self, evidences):

        self.evidences += evidences


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
        elif type(source) in _const.CHAR_TYPES:
            self._add_source(source)
        else:
            for s in source:
                self._add_source(s)

    def _add_source(self, source):
        self.sources.add(source)

    def add_refs(self, refs):
        self.refs = common.add_to_set(self.refs, refs)

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


#TODO this class does not belong here, find a better place
class Complex(object):

    have_stoichiometry = {
        'PDB',
        'Compleat',
        'ComplexPortal',
        'CellPhoneDB',
    }


    def __init__(
            self,
            components,
            ncbi_tax_id = 9606,
            name = None,
            ids = None,
            sources = None,
            interactions = None,
            references = None,
            proteins = None,
            attrs = None,
        ):
        """
        Represents a molecular complex.

        components : list,dict
            Either a list of identifiers or a dict with identifiers as keys
            and stoichiometric coefficients as values. List of identifiers
            also assumed to represent stoichiometry by repetition
            of identifiers.
        ncbi_tax_id : int
            NCBI taxonomy identifier of the complex. It implies all members
            of the complex belong to the same organism. Support for multi-
            organism complexes will be implemented in the future.
        name : str
            A custom name or identifier of the complex.
        ids : dict
            Identifiers. If ``sources`` is a set, list or tuple it should be
            a dict with database names as keys and set of identifiers as
            values. If ``sources`` is a string, it can be a set of
            identifiers or a single identifier.
        sources : set,str
            Database(s) the complex has been defined in.
        interactions : list,dict
            Interactions between the components of the complex. Either
            a list of tuples of component IDs or a dict with tuples as
            keys and custom interaction properties as values.
        proteins : list,dict
            Synonym for `components`, kept for compatibility.
        """

        components = components or proteins

        if not isinstance(components, dict):

            self.components = dict(collections.Counter(components))

        else:

            self.components = components

        self.proteins = self.components
        self.name = name
        self.ids = collections.defaultdict(set)
        self.add_ids(ids, source = sources)
        self.sources = common.to_set(sources)
        self.references = common.to_set(references)
        self.ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(ncbi_tax_id)
        self.attrs = {}
        if isinstance(attrs, dict):
            self.attrs.update(attrs)

        self.interactions = interactions


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        import importlib as imp
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def __str__(self):

        return 'COMPLEX:%s' % (
            COMPLEX_SEP.join(sorted(self.components.keys()))
        )


    def __repr__(self):

        return 'Complex%s: %s' % (
            ' %s' % self.name if self.name else '',
            self.__str__(),
        )


    def __hash__(self):

        return hash(self.__str__())


    def __contains__(self, other):

        return other in self.components


    def __eq__(self, other):

        return self.__hash__() == other.__hash__()


    def __iadd__(self, other):

        self.merge(other)

        return self


    def __lt__(self, other):

        return self.__str__() < other


    def __gt__(self, other):

        return self.__str__() > other


    def __len__(self):

        return len(self.components)


    def merge(self, other):
        """
        Adds the annotations (sources, references, attrs) of the other
        ``Complex`` instance to this one. If the other ``Complex`` has
        different components it does nothing.
        """

        if self != other:

            return

        if (
            set(self.components.values()) == {1} and
            set(other.components.values()) != {1}
        ):
            # this complex has no stoichiometry information
            # but the other has
            self.components = other.components

        self.sources.update(other.sources)
        self.references.update(other.references)

        self.add_ids(other.ids)

        for k, v in iteritems(other.attrs):

            if k not in self.attrs:

                self.attrs[k] = v

            elif isinstance(self.attrs[k], (dict, set)):

                self.attrs[k].update(v)


    def add_ids(self, ids, source = None):

        if not isinstance(ids, dict):

            ids = common.to_set(ids)

        if isinstance(ids, set) and source:

            source = common.to_set(source)

            ids = dict((s, ids) for s in source)

        if isinstance(ids, dict):

            for this_source, this_ids in iteritems(ids):

                this_ids = common.to_set(this_ids)
                self.ids[this_source].update(this_ids)


    def get_interaction(self, component1, component2):

        if self.has_interaction(component1, component2):

            return self.interactions[(component1, component2)]


    def set_interaction(self, component1, component2, interaction):

        key = (component1, component2)

        self.interactions = self.interactions or {}
        self.interactions[key] = interaction


    def has_interaction(self, component1, component2):

        key = (component1, component2)

        return self.interactions and key in self.interactions


    def add_source(self, source):

        self.sources.add(source)


    def iter_proteins(self):

        for protein in self.proteins.keys():

            yield protein


    __iter__ = iter_proteins


    def add_attr(self, source, attr):
        """
        Attributes can store annotations for complexes.
        """

        self.attrs[source] = attr


    @property
    def stoichiometry(self):

        return ':'.join(
            '%u' % (
                cnt
                    if self.sources & self.have_stoichiometry else
                0
            )
            for _id, cnt in
            sorted(
                iteritems(self.components),
                key = lambda id_cnt: id_cnt[0],
            )
        )


    @property
    def stoichiometry_str(self):

        return ';'.join(
            itertools.chain(*(
                (comp,) * cnt
                for comp, cnt in
                sorted(
                    iteritems(self.components),
                    key = lambda comp_cnt: comp_cnt[0],
                )
            ))
        )


    @property
    def stoichiometry_str_genesymbols(self):

        return ';'.join(
            itertools.chain(*(
                (
                    (
                        mapping.map_name0(
                            uniprot,
                            'uniprot',
                            'genesymbol',
                        ) or
                        uniprot
                    ),
                ) * cnt
                for uniprot, cnt in
                sorted(
                    iteritems(self.components),
                    key = lambda comp_cnt: comp_cnt[0],
                )
            ))
        )


    @property
    def genesymbols(self):

        return sorted(
            (
                mapping.map_name0(uniprot, 'uniprot', 'genesymbol') or
                uniprot
            )
            for uniprot in self.components.keys()
        )


    @property
    def genesymbol_str(self):

        return COMPLEX_SEP.join(self.genesymbols)


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

            nbonds[t] = len(self.__dict__[t][self.id_a])

        return nbonds


    def bond_types(self):
        '''
        Returns the bond types present in this interface
        '''
        types = []

        for t in self.types:

            if len(self.__dict__[t][self.id_a]) > 0:

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

                for i in range(0, len(self.__dict__[t][self.id_a])):

                    if mode == 'dict':

                        yield {
                            self.id_a: self.__dict__[t][self.id_a][i],
                            self.id_b: self.__dict__[t][self.id_b][i],
                            'type': t,
                        }

                    else:

                        yield (
                            (self.id_a,) +
                            (self.__dict__[t][self.id_a][i].serialize(),) +
                            (self.id_b,) +
                            (self.__dict__[t][self.id_b][i].serialize(),) +
                            (t,)
                        )


    def serialize(self):

        res = []
        for t in self.types:

            if self.__dict__[t][self.id_a] and self.__dict__[t][self.id_b]:

                res.append(
                    '%s:%s+%s' % (
                        t,
                        ','.join(self.__dict__[t][self.id_a].serialize()),
                        ','.join(self.__dict__[t][self.id_b].serialize()),
                    )
                )

        return (
            '%s-%u:%s-%u:%s:%s:%s' % (
                self.id_a,
                self.isoform_a,
                self.id_b,
                self.isoform_b,
                self.source,
                self.pdb,
                ':'.join(res),
            )
        )


    def __str__(self):

        nbonds = self.numof_residues()

        return (
            'Molecular interface between %s and %s,\n'
            'as observed in PDB structure %s\n\n'
            ' Data source: %s\n'
            ' Number of residues in contact: %u\n'
            ' Hydrogene bonds: %u\n'
            ' Covalent bonds: %u\n'
            ' Saltbridges: %u\n'
            ' S-S bonds: %u\n'
            ' Stable energy: %s\n'
            ' Solvation energy: %s\n'
            ' Surface area: %s\n'
            ' Complexation significance score: %s\n' % (
                self.id_a,
                self.id_b,
                self.pdb,
                self.source,
                sum(nbonds.values()),
                nbonds['hbonds'],
                nbonds['covbonds'],
                nbonds['sbridges'],
                nbonds['ssbonds'],
                'n/a' if self.stab_en is None else str(self.stab_en),
                'n/a' if self.solv_en is None else str(self.solv_en),
                'n/a' if self.area is None else str(self.area),
                'n/a' if self.css is None else str(self.css),
            )
        )


    def __repr__(self):

        nbonds = self.numof_residues()

        return (
            'Interface [%s-%s, %u bonds]' % (
                self.id_a,
                self.id_b,
                sum(nbonds.values()),
            )
        )
