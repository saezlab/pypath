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

class ProteinHomology(object):
    
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
        self.target = target
        self.set_default_source(source)
        self.mapper = mapping.Mapper() if mapper is None else mapper
        
        self.homologene_uniprot_dict()
    
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
        
        source = self.get_source(source)
        
        if source not in self.homo:
            
            self.homologene_uniprot_dict(source)
        
        if protein in self.homo[source]:
            
            return self.homo[source][protein]
    
    def homologene_uniprot_dict(self, source):
        """
        Builds orthology translation table as dict from UniProt to Uniprot,
        obtained from NCBI HomoloGene data. Uses RefSeq and Entrez IDs for
        translation.
        """
        
        self.homo[source] = {}
        
        hge = dataio.homologene_dict(source, self.target, 'entrez')
        hgr = dataio.homologene_dict(source, self.target, 'refseq')
        
        all_source[source] = set(dataio.all_uniprots(
            organism = source,
            swissprot = 'YES'))
        
        if not self.only_swissprot:
            all_source_trembl = dataio.all_uniprots(
                organism = source, swissprot = 'NO')
            all_source[source].update(set(all_source_trembl))
        
        for u in all_source[source]:
            
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
                    e, 'entrez', 'uniprot', target)))
            
            for r in target_r:
                target_u.update(set(self.mapper.map_name(
                    e, 'refseqp', 'uniprot', target)))
            
            target_u = \
                itertools.chain(
                    *map(
                        lambda tu:
                            self.mapper.map_name(
                                tu, 'uniprot', 'uniprot', target),
                        target_u
                    )
                )
            
            self.homo[source][u] = sorted(list(target_u))

class PtmHomology(ProteinHomology):
    
    __init__(self, target, source = None, only_swissprot = True,
             mapper = None, seq = None):
        
        super(PtmHomology, self).__init__(target,
                                          source,
                                          only_swissprot,
                                          mapper)
        
        self.reptm = re.compile(r'([A-Z\d]{6,10})_([A-Z])(\d*)')
        
        self.ptm_orthology()
        self.set_seq(seq)
    
    def set_seq(self, seq = None):
        
        self.seq = (
            seq or
            uniprot_input.swissprot_seq(self.target, isoforms = True)
        )
    
    def translate_site(self, protein, res, offset,
                       source_taxon = None,
                       isoform = 1, typ = 'phosphorylation',
                       strict = False):
        """
        Translates one PTM site.
        """
        
        result = set([])
        
        self.set_default_source(source_taxon)
        
        source = self.get_source(source)
        
        sourceptm = (protein, isoform, res, offset, source, typ)
        
        if sourceptm in self.ptmhomo:
            
            if self.target in self.ptmhomo[sourceptm]:
                
                result = self.ptmhomo[sourceptm]
        
        if not result and not strict:
            
            tsubs = super(PTMHomology, self).translate(protein)
            
            for tsub in tsubs:
                
                if tsub not in self. tseq:
                    continue
                
                for toffset in xrange(offset, offset + 3):
                    
                    for i in self.seq[tsub].isoforms():
                        
                        tres = self.seq[tsub].get(toffset, isoform = i)
                        
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
                    super(PTMHomology, self).translate(domain.protein)
                )
            )
        )
    
    def translate_ptm(self, ptm, strict = False):
        
        tptms = self.translate_site(ptm.protein,
                                    ptm.motif.residue.name,
                                    ptm.motif.residue.number,
                                    ptm.motif.residue.isoform,
                                    ptm.typ,
                                    strict)
        
        result = []
        
        for x in tptms:
            
            res = intera.Residue(x[3], x[2], x[0], isoform = x[1])
            start, end, region = (
                self.tseq[x[0]].get_region(x[3], isoform = x[1])
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
            
            ptm = reptm.match(x)
            
            if ptm is not None:
                
                result = self.translate_site(ptm[0], ptm[1],
                                             int(ptm[2]), **kwargs)
        
        if return_strings:
            
            result = list(map(lambda r:
                              '%s_%s%u' % (r[0], r[2], r[3]),
                              result))
            
        elif type(x) is pypath.intera.Ptm:
            
            result = self.translate_ptm(x)
            
        elif type(x) is pypath.intera.Domain:
            
            result = self.translate_domain(x)
            
        elif type(x) is pypath.intera.DomainMotif:
            
            result = self.translate_domain_motif(x)
        
        return result
    
    def ptm_orthology():
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
