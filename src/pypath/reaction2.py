#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from collections import Counter

import pypath.mapping
import pypath.intera
from pypath.common import *
import pypath.dataio


class Reaction(object):
    '''
    Standard representation of a reaction from reaction networks,
    e.g. ACSN, NCI-PID, PANTHER or Reactome.
    '''

    def __init__(self,
                 source,
                 reactants=[],
                 products=[],
                 references=[],
                 id_type='uniprot',
                 names_reactants=None,
                 names_products=None,
                 ptms_reactants=None,
                 ncbi_tax_id=9606,
                 ptms_products=None,
                 mapper=None,
                 seq=None):
        '''
        source : str
        Source database, e.g. `Reactome`.

        reactants : list
        List of reactants. List of string IDs or list of lists, if
        reactants are complexes.

        products : list
        List of products. List of string IDs or list of lists, if
        products are complexes.

        id_type : str
        Type of reactant and product IDs. Default is UniProt ID.

        mapper : pypath.mapping.Mapper
        If no Mapper instance given, a new one will be initialized.
        '''
        self.mapper = mapper if mapper is not None else mapping.Mapper()
        self.source = source
        self.references = references
        self.ncbi_tax_id = ncbi_tax_id
        self.reactants = {}
        self.products = {}
        reactants = (ids if type(ids) is list else [ids] for ids in reactants)
        products = (ids if type(ids) is list else [ids] for ids in products)
        names_reactants = names_reactants if type(names_reactants) is list \
            else (self.species_name(ids) for ids in self.reactants)
        names_reactants = dict(
            zip((species_id(ids) for ids in reactants), names_reactants))
        names_products = names_products if type(names_products) is list \
            else (self.species_name(ids) for ids in self.products)
        names_products = dict(
            zip((species_id(ids) for ids in products), names_products))
        self.add_species('reactants', reactants, names_reactants, source)
        self.add_species('products', products, names_products, source)
        if source == 'Reactome':
            for i, name in names_reactants.iteritems():
                ptms = self.reactome_ptms(i, name)

    def _add_species(role, ids, name, source):
        '''
        Adds a complex to dict of reactants.
        Single proteins considered a one member "complex".
        '''
        ids = ids if type(ids) is list else [ids]
        ids = self.species_id(ids)
        self.getattr(role)[ids] = Complex(ids, name, name, source)

    def add_species(role, species, names, source):
        for sp, name in zip(species, names):
            sp = sp if type(sp) is list else [sp]
            self._add_species(role, sp, name)

    def species_name(ids):
        return '-'.join([i for i in sorted(ids)])

    def species_id(ids):
        return tuple(sorted(ids))

    def sequences(self, isoforms=True):
        self.seq = dataio.swissprot_seq(self.ncbi_tax_id, isoforms)

    def reactome_ptms(i, name, role):
        renum = re.compile(r'[0-9]+')
        if self.seq is None:
            self.sequences(isoforms=True)
        name = name.split()[0].split(':')
        for protein in name:
            ptms = []
            ptmstr, gs = protein.split('-')
            uniprots = self.mapper.map_name(gs, 'genesymbol', 'uniprot')
            for u in uniprots:
                if u in self.getattr(role)[i] and u in self.seq:
                    s = self.seq[u]
                    for ptm in ptmstr.split(','):
                        num = renum.findall()
                        pass


class Species(object):
    '''
    Represents one species, i.e. reactant or product
    in a reaction.
    '''

    def __init__(self, identifier, ptms=None):
        '''
        identifier : str
        UniProt AC of the protein.

        ptms : list
        List 3 element tuples, containing
        residue name, number and the PTM type (e.g. phosphorylation).
        '''
        self.id = identifier
        self.ptms = []
        if ptms is not None:
            for name, number, typ in ptms:
                res = intera.Residue(number, name, identifier)
                ptm = intera.Ptm(identifier, residue=res, typ=typ)
