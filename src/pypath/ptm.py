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

import sys
import imp

import pypath.dataio as dataio
import pypath.common as common
import pypath.mapping as mapping
import pypath.homology as homology

class PtmProcessor(object):
    
    methods = {
        'Signor': 'load_signor_ptms',
        'MIMP': 'get_mimp',
        'PhosphoNetworks': 'get_phosphonetworks',
        'phosphoELM': 'get_phosphoelm',
        'dbPTM': 'get_dbptm',
        'PhosphoSite': 'get_psite_phos',
        'HPRD': 'get_hprd_ptms',
        'Li2012': 'li2012_phospho'
    }
    
    organisms_supported = set(['Signor', 'PhosphoSite',
                               'phosphoELM', 'dbPTM'])
    
    enzyme_id_uniprot = set(['PhosphoSite', 'phosphoELM', 'Signor'])
    
    substrate_id_types = {
        'MIMP': [('genesymbol', 'substrate'), ('refseq', 'substrate_refseq')],
        'PhosphoNetworks': ['genesymbol'],
        'phosphoELM': ['uniprot'],
        'Li2012': ['genesymbol'],
        'dbPTM': ['uniprot'],
        'PhosphoSite': ['uniprot'],
        'Signor': ['uniprot'],
        'HPRD': [('refseqp', 'substrate_refseqp')]
    }
    
    __init__(self, source,
             ncbi_tax_id = 9606,
             trace = False,
             mapper = None,
             seq = None,
             enzyme_id_type = 'genesymbol',
             substrate_id_type = 'genesymbol',
             **kwargs):
        """
        Processes enzyme-substrate interaction data from various databases.
        Provedes generators to iterate over these interactions.
        For organisms other than human obtains the organism specific
        interactions from databases whereever it is possible and translates
        the human interactions by homology. By default it does both of them
        and iterates over all the interactions.
        
        :param str source: Either a method name in `dataio` or a database
                           name e.g. `PhosphoSite` or a callable which
                           returns data in list of dicts format.
        :param int ncbi_tax_id: NCBI Taxonomy ID used at the database lookups.
        :param bool trace: Keep data about ambiguous ID mappings and PTM data
                           in mismatch with UniProt sequences.
        :param list map_by_homology_from: Look up by these taxons in database
                                          and map them by homology.
        :param int map_by_homology_to: The target taxon of the homology
                                       translation.
        :param pypath.mapping.Mapper: A `Mapper` instance. If `None` a new
                                      instance will be created.
        :param str enzyme_id_type: The ID type of the enzyme in the database.
        :param str substrate_id_type: The ID type of the substrate in the
                                      database.
        
        :param bool nonhuman_direct_lookup: Use direct lookup at non-human
                                            target species.
        :param **kwargs: Args to be forwarded to the input method.
        
        """
        
        self.nomatch = []
        self.kin_ambig = {}
        self.sub_ambig = {}
        
        self.seq = seq
        self.source = source
        self.trace = trace
        self.ncbi_tax_id = ncbi_tax_id
        self.mapper = mapper
        self.enzyme_id_type = enzyme_id_type
        self.set_method()
        self.set_inputargs(**kwargs)
        self.init_mapper()
        self.load()
    
    def load(self):
        
        self._setup()
        self.load_seq()
        self.load_data()
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def reset_ptmprocessor(self, seq = None, ncbi_tax_id = None):
        
        ncbi_tax_id = ncbi_tax_id or self.ncbi_tax_id
        
        self.set_taxon(ncbi_tax_id)
        self.load_seq(seq)
        self.load_data()
    
    def set_taxon(self, ncbi_tax_id):
        self.ncbi_tax_id = ncbi_tax_id
        self._organism_setup()
    
    def set_method(self):
        """
        Selects the input method.
        """
        
        def f(**kwargs): return []
        
        if hasattr(self.source, '__call__'):
            self.inputm = self.source
        elif hasattr(dataio, self.source):
            self.inputm = getattr(dataio, self.source)
        elif self.source in self.methods:
            self.inputm = self.methods[self.source]
        else:
            self.inputm = f
    
    def set_inputargs(self, **inputargs):
        """
        Sets the arguments to be provided for the input method.
        """
        
        self.inputargs = inputargs
    
    def load_data(self):
        """
        Loads the data by the defined input method.
        """
        
        self.data = self.inputm(**self.inputargs)
    
    def init_mapper(self):
        
        if self.mapper is None:
            
            self.mapper = mapping.Mapper()
    
    def _phosphosite_setup(self):
        
        self.ncbi_tax_id in common.taxids:
        
        if 'strict' not in self.inputargs:
            self.inputargs['strict'] = False
        
        self.inputargs['mapper'] = self.mapper
    
    def _phosphoelm_setup(self):
        
        if self.ncbi_tax_id != 9606 and 'ltp_only' not in self.inputargs:
            
            self.inputargs['ltp_only'] = False
    
    def _setup(self):
        
        setupmethod = '_%s_setup' % self.source.lower()
        
        if hasattr(self, setupmethod):
            
            getattr(self, setupmethod)()
        
        self._organism_setup()
    
    def _organism_setup(self):
        
        if self.source in self.organisms_supported:
            
            self.inputargs['organism'] = self.ncbi_tax_id
    
    def load_seq(self, seq = None):
        
        self.seq = (
            seq or
            uniprot_input.swissprot_seq(organism = self.ncbi_tax_id,
                                        isoforms=True)
        )
    
    def _process(self, p):
        
        for p in self.data:
            
            if p['kinase'] is not None and len(p['kinase']) > 0:
                
                # database specific id conversions
                if source in self.enzyme_id_uniprot:
                    self.enzyme_id_type = 'uniprot'
                
                if not isinstance(p['kinase'], list):
                    p['kinase'] = [p['kinase']]
                
                kinase_ups = self.mapper.map_name(p['kinase'],
                                self.enzyme_id_type,
                                'uniprot',
                                ncbi_tax_id = self.ncbi_tax_id)
                
                if p['substrate'].startswith('HLA'):
                    # human leukocyte antigenes result a result an
                    # extremely high number of combinations
                    continue
                
                self.substrate_ups_all = set([])
                
                for sub_id_type in (
                    self.substrate_id_types[self.source]
                    if self.source in self.substrate_id_types else
                    [self.substrate_id_type]
                ):
                    
                    if type(sub_id_type) is tuple:
                        sub_id_type, sub_id_attr = sub_id_type
                    else:
                        sub_id_attr = 'substrate'
                    
                    self.substrate_ups_all.update(
                        set(
                            self.mapper.map_name(
                                p[sub_id_attr],
                                sub_id_type,
                                'uniprot',
                                self.ncbi_tax_id
                            )
                        )
                    )
                
                # looking up sequences in all isoforms:
                substrate_ups = []
                
                for s in substrate_ups_all:
                    
                    if s in self.seq:
                        
                        for isof in self.seq[s].isoforms():
                            
                            if p['instance'] is not None:
                                
                                if self.seq[s].match(
                                        p['instance'],
                                        p['start'],
                                        p['end'],
                                        isoform=isof):
                                    substrate_ups.append((s, isof))
                            else:
                                
                                if self.seq[s].match(
                                        p['resaa'],
                                        p['resnum'],
                                        isoform=isof):
                                    
                                    substrate_ups.append((s, isof))
                
                if self.trace:
                    
                    if p['substrate'] not in self.sub_ambig:
                        
                        self.sub_ambig[p['substrate']] = substrate_ups
                    
                    for k in p['kinase']:
                        
                        if k not in self.kin_ambig:
                            
                            self.kin_ambig[k] = kinase_ups
                    # generating report on non matching substrates
                    if len(substrate_ups) == 0:
                        
                        for s in substrate_ups_all:
                            
                            if s[0] in self.seq:
                                
                                nomatch.append((s[0], s[1],
                                    ((p['substrate_refseq']
                                    if 'substrate_refseq' in p
                                    else ''),
                                    s, p['instance'],
                                    self.seq[s].get(
                                        p['start'],
                                        p['end'])
                                    )
                                ))
                
                # adding kinase-substrate interactions
                for k in kinase_ups:
                    
                    for s in substrate_ups:
                        
                        res = intera.Residue(
                            p['resnum'],
                            p['resaa'],
                            s[0],
                            isoform=s[1])
                        
                        if p['instance'] is None:
                            
                            reg = self.seq[s[0]].get_region(
                                p['resnum'],
                                p['start'],
                                p['end'],
                                isoform=s[1])
                            
                            if reg is not None:
                                
                                p['instance'] = reg[2]
                                p['start'] = reg[0]
                                p['end'] = reg[1]
                                
                        if 'typ' not in p:
                            p['typ'] = 'phosphorylation'
                            
                        mot = intera.Motif(
                            s[0],
                            p['start'],
                            p['end'],
                            instance=p['instance'],
                            isoform=s[1])
                        
                        ptm = intera.Ptm(s[0],
                                            motif=mot,
                                            residue=res,
                                            typ=p['typ'],
                                            source=[self.source],
                                            isoform=s[1])
                        
                        dom = intera.Domain(protein=k)
                        
                        if 'references' not in p:
                            p['references'] = []
                            
                        dommot = intera.DomainMotif(
                            domain=dom,
                            ptm=ptm,
                            sources=[self.source],
                            refs=p['references'])
                        
                        if self.source == 'MIMP':
                            dommot.mimp_sources = ';'.split(p[
                                'databases'])
                            dommot.npmid = p['npmid']
                            
                        elif source == 'PhosphoNetworks':
                            dommot.pnetw_score = p['score']
                            
                        elif source == 'dbPTM':
                            dommot.dbptm_sources = [p['source']]
                            
                        yield dommot
    
    def __iter__(self):
        """
        Iterates through the enzyme-substrate interactions.
        """
        
        for p in self.data:
            
            for ptm in self._process(p)
                
                yield ptm

    
class PtmHomologyProcessor(PtmProcessor, homology.PtmHomology):
    
    def __init__(self,
        source,
        ncbi_tax_id,
        map_by_homology_from = [9606],
        trace = False,
        mapper = None,
        enzyme_id_type = 'genesymbol',
        substrate_id_type = 'genesymbol',
        homology_only_swissprot = True,
        ptm_homology_strict = False,
        **kwargs):
        """
        Unifies a `pypath.ptm.PtmProcessor` and
        a `pypath.homology.PtmHomology` object to build a set of
        enzyme-substrate interactions from a database and subsequently
        translate them by homology to one different organism.
        Multiple organism can be chosen as the source of the
        enzyme-substrate interactions. For example if you want mouse
        interactions, you can translate them from human and from rat.
        To get the original mouse interactions themselves, use an
        other instance of the `PtmProcessor`.
        To have both the original and the homology translated set,
        and also from multiple databases, whatmore all these merged
        into a single set, use the `PtmAggregator`.
        
        :param str source: Data source for `PtmProcessor`.
        :param int ncbi_tax_id: The NCBI Taxonomy ID the interactions
                                should be translated to.
        :param bool homology_only_swissprot: Use only SwissProt
                                             (i.e. not Trembl) at homology
                                             translation.
        :param bool ptm_homology_strict: Use only those homologous PTM pairs
                                         which are in PhosphoSite data, i.e.
                                         do not look for residues with same
                                         offset in protein sequence.
        
        See further options at `PtmProcessor`.
        
        """
        
        self.map_by_homology_from = map_by_homology_from

        PtmProcessor.__init__(self, source, ncbi_tax_id,
                              trace = trace, mapper = mapper,
                              enzyme_id_type = enzyme_id_type,
                              substrate_id_type = substrate_id_type,
                              **kwargs)
        
        homology.PtmTranslator.__init__(self, target = self.ncbi_tax_id,
                                        only_swissprot = homology_only_swissprot,
                                        strict = ptm_homology_strict)
    
    def __iter__(self):
        """
        Iterates through enzyme-substrate interactions
        translated to another organism by orthology.
        """
        
        for source_taxon in self.map_by_homology_from:
            
            self.set_default_source(source_taxon)
            self.reset_ptmprocessor(ncbi_tax_id = source_taxon)
            
            for ptm in PtmProcessor.__iter__(self):
                
                for tptm in self.translate(ptm):
                    
                    yield tptm
    
class PtmAggregator(object):
    
    def __init__(self,
        sources,
        ncbi_tax_id = 9606,
        map_by_homology_from = [9606],
        trace = False,
        mapper = None,
        enzyme_id_type = 'genesymbol',
        substrate_id_type = 'genesymbol',
        homology_only_swissprot = True,
        ptm_homology_strict = False,
        nonhuman_direct_lookup = True,
        inputargs = {}):
        """
        Docs not written yet.
        """
        
        for k, v in iteritems(locals()):
            setattr(self, k, v)
        
        self.init_mapper()
        
        self.map_by_homology_from = set(self.map_by_homology_from)
        self.map_by_homology_from.discard(self.ncbi_tax_id)
    
    def __iter__(self):
        
        for ptm in self.unique_list:
            
            yield ptm
    
    def build_list(self):
        """
        Builds a full list of enzyme-substrate interactions from
        all the requested sources. This list might contain redundant
        elements which later will be merged by `unique`.
        """
        
        self.full_list = []
        
        for source in sources:
            
            inputargs = (
                self.inputargs[source]
                if source in self.inputargs
                else {}
            )
            
            if ncbi_tax_id == 9606 or self.nonhuman_direct_lookup:
                
                proc = PtmProcessor(source, self.ncbi_tax_id, self.trace,
                                    self.mapper, self.enzyme_id_type,
                                    self.substrate_id_type, **inputargs)
                
                self.full_list.extend(list(proc.__iter__()))
            
            if self.map_by_homology_from:
                
                proc = PtmHomologyProcessor(source, self.ncbi_tax_id,
                                            self.map_by_homology_from,
                                            self.trace, self.mapper,
                                            self.enzyme_id_type,
                                            self.substrate_id_type,
                                            self.homology_only_swissprot,
                                            self.ptm_homology_strict,
                                            **inputargs)
                
                self.full_list.extend(list(proc.__iter__()))
    
    def unique(self):
        """
        Merges the redundant elements of the interaction list.
        Elements are redundant if they agree in all their attributes
        except the sources, references and isoforms.
        """
        
        self.unique_list = set([])
        
        ptms_uniq = []
        for ptm in ptms:
            merged = False
            for i, ptmu in enumerate(ptms_uniq):
                if ptm == ptmu:
                    ptms_uniq[i].merge(ptm)
                    merged = True
            if not merged:
                ptms_uniq.append(ptm)
        return ptms_uniq
    
    def init_mapper(self):
        
        self.mapper = self.mapper or mapping.Mapper()
    
    def assign_to_network(self, pa):
        """
        Assigns enzyme-substrate interactions to edges of a
        network in a `pypath.main.PyPath` instance.
        """
        
        pa.update_vname()
        if 'ptm' not in pa.graph.es.attributes():
            pa.graph.es['ptm'] = [[] for _ in pa.graph.es]
        
        for es in self:
            
            nodes = pa.get_node_pair(es.domain.protein, es.ptm.protein,
                    directed = pa.graph.is_directed())
            
            e = None
            if nodes:
                e = pa.graph.get_eid(
                    nodes[0], nodes[1], error=False)
            
            if isinstance(e, int) and e > 0:
                
                if pa.graph.es[e]['ptm'] is None:
                    pa.graph.es[e]['ptm'] = []
                
                pa.graph.es[e]['ptm'].append(es)
