#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from future.utils import iteritems
from past.builtins import xrange, range

import sys
import itertools
import imp
import re

import pypath.mapping as mapping
import pypath.dataio as dataio
import pypath.common as common
import pypath.intera as intera
import pypath.urls as urls
import pypath.curl as curl
import pypath.uniprot_input as uniprot_input
import pypath.seq as _se


class SequenceContainer(object):
    
    def __init__(self, preload_seq = [], isoforms = True):
        """
        This is an object to store sequences of multiple
        organisms and select the appropriate one.
        """
        
        self.seq_isoforms = isoforms
        
        for taxon in preload_seq:
            
            self.load_seq(taxon)
    
    def load_seq(self, taxon):
        
        if not hasattr(self, 'seq'):
            self.seq = {}
        
        taxon = taxon or self.ncbi_tax_id
        
        if taxon not in self.seq:
            
            self.seq[taxon] = _se.swissprot_seq(
                organism = taxon,
                isoforms = self.seq_isoforms)
    
    def get_seq(self, protein, taxon = None):
        
        if taxon is not None:
            
            if taxon not in self.seq:
                
                self.load_seq(taxon)
            
            if protein in self.seq[taxon]:
                
                return self.seq[taxon][protein]
        
        else:
            
            for taxon, seq in iteritems(self.seq):
                
                if protein in seq:
                    
                    return seq[protein]


class Proteomes(object):
    
    def __init__(self, preload_prot = [], swissprot_only = True):
        
        if not hasattr(self, '_taxonomy'):
            self._taxonomy = {}
            self._proteomes = {}
        
        for taxon in preload_prot:
            
            self.load_proteome(taxon, swissprot_only)
    
    def load_proteome(self, taxon, swissprot_only = True):
        
        key = (taxon, swissprot_only)
        
        if key not in self._proteomes:
            
            self._proteomes[key] = (
                set(uniprot_input.all_uniprots(*key))
            )
            
            for protein in self._proteomes[key]:
                
                self._taxonomy[protein] = key
            
            if not swissprot_only:
                
                self.load_proteome(taxon, True)
    
    def get_taxon(self, protein, swissprot_only = True):
        
        if not swissprot_only or self.is_swissprot(protein):
            
            return self._taxonomy[protein][0]
    
    def get_taxon_trembl(self, protein):
        
        if self.has_protein(protein):
            
            return self._taxonomy[protein][0]
        
    def has_protein(self, protein):
        
        return protein in self._taxonomy
    
    def is_swissprot(self, protein):
        
        return self.has_protein(protein) and self._taxonomy[protein][1]


class ProteinHomology(Proteomes):
    
    def __init__(self, target, source = None, only_swissprot = True, mapper = None):
        """
        This class translates between homologous UniProt IDs of
        2 organisms based on NCBI HomoloGene data.
        Uses RefSeq and Entrez IDs for translation.
        
        
        :param int target: NCBI Taxonomy ID of the organism
                           to be translated to.
        :param int source: NCBI Taxonomy ID of the default organism
                           to be translated from.
        :param bool only_swissprot: Whether only SwissProt or Trembl IDs
                                    should be used.
        :mapper pypath.mapping.Mapper mapper: A Mapper object.
        """
        
        self.homo = {}
        self.only_swissprot = only_swissprot
        self.target = target
        self.source = source
        self.set_default_source(source)
        self.mapper = mapping.Mapper() if mapper is None else mapper
        
        Proteomes.__init__(self)
        self.load_proteome(self.target, self.only_swissprot)
        
        if source is not None:
            self.homologene_uniprot_dict(source)
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def set_default_source(self, source = None):
        
        self.source = source or self.source
    
    def get_source(self, source = None):
        
        source = source or self.source
        
        if source is None:
            raise ValueError('No source NCBI Taxonomy ID provided.')
        else:
            return source
    
    def translate(self, protein, source = None):
        """
        For one UniProt ID of the source organism returns all orthologues
        from the target organism.
        """
        
        if self.get_taxon(protein) == self.target:
            return [protein]
        
        source = self.get_source(source)
        
        if source not in self.homo:
            
            self.homologene_uniprot_dict(source)
        
        if protein in self.homo[source]:
            
            return self.homo[source][protein]
        
        else:
            
            return []
    
    def homologene_uniprot_dict(self, source):
        """
        Builds orthology translation table as dict from UniProt to Uniprot,
        obtained from NCBI HomoloGene data. Uses RefSeq and Entrez IDs for
        translation.
        """
        
        source = self.get_source(source)
        
        self.homo[source] = {}
        
        hge = dataio.homologene_dict(source, self.target, 'entrez')
        hgr = dataio.homologene_dict(source, self.target, 'refseq')
        
        self.load_proteome(source, self.only_swissprot)
        
        for u in self._proteomes[(source, self.only_swissprot)]:
            
            source_e = self.mapper.map_name(
                u, 'uniprot', 'entrez', source)
            source_r = self.mapper.map_name(
                u, 'uniprot', 'refseqp', source)
            target_u = set([])
            target_r = set([])
            target_e = set([])
            
            for e in source_e:
                if e in hge:
                    target_e.update(hge[e])
            
            for r in source_r:
                if r in hgr:
                    target_r.update(hgr[r])
            
            for e in target_e:
                target_u.update(set(self.mapper.map_name(
                    e, 'entrez', 'uniprot', self.target)))
            
            for r in target_r:
                target_u.update(set(self.mapper.map_name(
                    e, 'refseqp', 'uniprot', self.target)))
            
            target_u = \
                itertools.chain(
                    *map(
                        lambda tu:
                            self.mapper.map_name(
                                tu, 'uniprot', 'uniprot', self.target),
                        target_u
                    )
                )
            
            self.homo[source][u] = sorted(list(target_u))

class PtmHomology(ProteinHomology,SequenceContainer):
    
    def __init__(self, target, source = None, only_swissprot = True,
             mapper = None, strict = True):
        
        ProteinHomology.__init__(self, target,
                                       source,
                                       only_swissprot,
                                       mapper)
        
        SequenceContainer.__init__(self)
        self.load_seq(taxon = self.target)
        
        self.reptm = re.compile(r'([A-Z\d]{6,10})_([A-Z])(\d*)')
        
        self.strict = strict
        
        self.ptm_orthology()
    
    def translate_site(self, protein, res, offset,
                       isoform = 1, typ = 'phosphorylation',
                       source_taxon = None):
        """
        Translates one PTM site.
        """
        
        result = set([])
        
        self.set_default_source(source_taxon)
        
        source = self.get_source(source_taxon)
        
        sourceptm = (protein, isoform, res, offset, source, typ)
        
        if self.get_taxon(protein) == self.target:
            result.add(sourceptm)
            return result
        
        if sourceptm in self.ptmhomo:
            
            if self.target in self.ptmhomo[sourceptm]:
                
                result = self.ptmhomo[sourceptm]
        
        if not result and not self.strict:
            
            tsubs = ProteinHomology.translate(self, protein, source = source)
            
            for tsub in tsubs:
                
                se = self.get_seq(tsub)
                
                if se is None:
                    continue
                
                for toffset in xrange(offset, offset + 3):
                    
                    for i in se.isoforms():
                        
                        tres = se.get(toffset, isoform = i)
                        
                        if tres == res:
                            
                            result.add((
                                tsub,
                                i,
                                tres,
                                toffset,
                                self.target,
                                typ
                            ))
                    
                    if result:
                        break
        
        return result
    
    def translate_domain(self, domain):
        
        return (
            list(
                map(
                    lambda x:
                        intera.Domain(x),
                    ProteinHomology.translate(
                        self,
                        domain.protein,
                        source = self.get_source()
                    )
                )
            )
        )
    
    def translate_ptm(self, ptm):
        
        tptms = self.translate_site(ptm.protein,
                                    ptm.residue.name,
                                    ptm.residue.number,
                                    ptm.residue.isoform,
                                    ptm.typ)
        
        result = []
        
        if self.target in tptms:
            
            for x in tptms[self.target]:
                
                se = self.get_seq(x[0])
                
                if (se is None or x[1] not in se.isof) and self.strict:
                    continue
                
                res = intera.Residue(x[3], x[2], x[0], isoform = x[1])
                start, end, region = (
                    se.get_region(x[3], isoform = x[1])
                    if se is not None and x[1] in se.isof
                    else (None, None, None)
                )
                mot = intera.Motif(x[0], start = start, end = end,
                                instance = region,
                                isoform = x[1])
                
                ptm = intera.Ptm(x[0], motif = mot, residue = res,
                                typ = x[5], isoform = x[1],
                                source = ptm.sources)
                
                result.append(ptm)
        
        return result
    
    def translate_domain_motif(self, dmotif):
        
        ds = self.translate_domain(dmotif.domain)
        ps = self.translate_ptm(dmotif.ptm)
        
        return (
            list(
                map(
                    lambda x:
                        intera.DomainMotif(x[0], x[1],
                                           sources = dmotif.sources,
                                           refs = dmotif.refs),
                        itertools.product(ds, ps)
                )
            )
        )
    
    def translate_residue(self, residue):
        
        return (
            list(
                map(
                    lambda r:
                        intera.Residue(r[3], r[2], r[0], isoform = r[1]),
                    self.translate_site(
                        residue.protein,
                        residue.name,
                        residue.number,
                        residue.isoform
                    )
                )
            )
        )
    
    def translate(self, x, return_strings = False, **kwargs):
        """
        Translates anything:
        
        - one UniProt ID
        - one PTM provided as tuple of (UniProt, amino acid, offest)
        - one PTM provided as string (e.g. `P00533_S231`)
        - instance from pypath.intera: DomainMotif, Domain or Ptm
        
        Additional arguments can be isoform and typ (modification type).
        
        """
        
        result = []
        
        if type(x) is tuple:
            
            result = self.translate_site(*x, **kwargs)
        
        elif type(x) in common.charTypes:
            
            ptm = self.reptm.match(x)
            
            if ptm is not None:
                
                result = self.translate_site(ptm[0], ptm[1],
                                             int(ptm[2]), **kwargs)
        
        if return_strings:
            
            result = list(map(lambda r:
                              '%s_%s%u' % (r[0], r[2], r[3]),
                              result))
            
        elif type(x) is intera.Ptm:
            
            result = self.translate_ptm(x)
            
        elif type(x) is intera.Domain:
            
            result = self.translate_domain(x)
            
        elif type(x) is intera.DomainMotif:
            
            result = self.translate_domain_motif(x)
        
        return result
    
    def ptm_orthology(self):
        """
        Creates an orthology translation dict of phosphosites
        based on phosphorylation sites table from PhosphoSitePlus.
        In the result all PTMs represented by a tuple of the following
        6 elements: UniProt ID, isoform (int), residue one letter code,
        residue number (int), NCBI Taxonomy ID (int), modification type.
        
        """
        
        self.ptmhomo = {}
        
        nondigit = re.compile(r'[^\d]+')
        
        unknown_taxa = set([])
        
        for typ in common.psite_mod_types:
            
            groups = {}
            
            url = urls.urls['psite_%s' % typ[0]]['url']
            c = curl.Curl(url, silent=False, large=True)
            
            data = c.result
            
            for _ in xrange(4):
                null = data.readline()
            
            for r in data:
                
                r = r.decode('utf-8').split('\t')
                
                if len(r) < 10:
                    
                    continue
                
                uniprot = r[2]
                isoform = 1 if '-' not in uniprot else int(uniprot.split('-')[1])
                uniprot = uniprot.split('-')[0]
                aa = r[4][0]
                num = int(nondigit.sub('', r[4]))
                if r[6] not in common.taxa:
                    unknown_taxa.add(r[6])
                    continue
                
                tax = common.taxa[r[6]]
                group = int(r[5])
                
                this_site = (uniprot, isoform, aa, num, tax, typ[1])
                
                if group not in groups:
                    groups[group] = set([])
                
                groups[group].add(this_site)
            
            for group, sites in iteritems(groups):
                
                for site1 in sites:
                    
                    for site2 in sites:
                        
                        if site1[4] == site2[4]:
                            
                            continue
                        
                        if site1 not in self.ptmhomo:
                            
                            self.ptmhomo[site1] = {}
                        
                        if site2[4] not in self.ptmhomo[site1]:
                            
                            self.ptmhomo[site1][site2[4]] = set([])
                        
                        self.ptmhomo[site1][site2[4]].add(site2)
        
        if len(unknown_taxa):
            sys.stdout.write('\t:: Unknown taxa encountered:\n\t   %s\n' %
                             ', '.join(sorted(unknown_taxa)))
