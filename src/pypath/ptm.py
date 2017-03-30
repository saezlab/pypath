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
import imp
import itertools

import pypath.dataio as dataio
import pypath.common as common
import pypath.mapping as mapping
import pypath.homology as homology
import pypath.uniprot_input as uniprot_input
import pypath.intera as intera
import pypath.progress as progress

class PtmProcessor(homology.Proteomes,homology.SequenceContainer):
    
    methods = {
        'signor': 'load_signor_ptms',
        'mimp': 'get_mimp',
        'phosphonetworks': 'get_phosphonetworks',
        'phosphoelm': 'get_phosphoelm',
        'dbptm': 'get_dbptm',
        'phosphosite': 'get_psite_phos',
        'hprd': 'get_hprd_ptms',
        'li2012': 'li2012_phospho'
    }
    
    organisms_supported = set(['signor', 'phosphosite',
                               'phosphoelm', 'dbptm'])
    
    enzyme_id_uniprot = set(['phosphosite', 'phosphoelm', 'signor'])
    
    substrate_id_types = {
        'mimp': [('genesymbol', 'substrate'), ('refseq', 'substrate_refseq')],
        'phosphonetworks': ['genesymbol'],
        'phosphoelm': ['uniprot'],
        'li2012': ['genesymbol'],
        'dbptm': ['uniprot'],
        'phosphosite': ['uniprot'],
        'signor': ['uniprot'],
        'hprd': [('refseqp', 'substrate_refseqp')]
    }
    
    def __init__(self, input_method,
             ncbi_tax_id = 9606,
             trace = False,
             mapper = None,
             enzyme_id_type = 'genesymbol',
             substrate_id_type = 'genesymbol',
             name = None,
             allow_mixed_organisms = False,
             **kwargs):
        """
        Processes enzyme-substrate interaction data from various databases.
        Provedes generators to iterate over these interactions.
        For organisms other than human obtains the organism specific
        interactions from databases.
        
        :param str input_method: Either a method name in `dataio` or a database
                                 name e.g. `PhosphoSite` or a callable which
                                 returns data in list of dicts format.
        :param int ncbi_tax_id: NCBI Taxonomy ID used at the database lookups.
        :param bool trace: Keep data about ambiguous ID mappings and PTM data
                           in mismatch with UniProt sequences.
        :param pypath.mapping.Mapper: A `Mapper` instance. If `None` a new
                                      instance will be created.
        :param str enzyme_id_type: The ID type of the enzyme in the database.
        :param str substrate_id_type: The ID type of the substrate in the
                                      database.
        
        :param bool nonhuman_direct_lookup: Use direct lookup at non-human
                                            target species.
        :param **kwargs: Args to be forwarded to the input method.
        
        """
        
        self.mammal_taxa = set([9606, 10090, 10116])
        self.nomatch = []
        self.kin_ambig = {}
        self.sub_ambig = {}
        
        self.name = name
        self.allow_mixed_organisms = allow_mixed_organisms
        self.input_method = input_method
        self.trace = trace
        self.ncbi_tax_id = ncbi_tax_id
        
        homology.SequenceContainer.__init__(self)
        self.load_seq(self.ncbi_tax_id)
        
        if self.allow_mixed_organisms:
            for taxon in self.mammal_taxa:
                self.load_seq(taxon = taxon)
        
        homology.Proteomes.__init__(self)
        
        self.mapper = mapper
        self.enzyme_id_type = enzyme_id_type
        self.set_method()
        self.set_inputargs(**kwargs)
        self.init_mapper()
        self.load()
    
    def load(self):
        
        self._setup()
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
        self.load_seq(ncbi_tax_id)
        self.load_data()
    
    def set_taxon(self, ncbi_tax_id):
        self.ncbi_tax_id = ncbi_tax_id
        self._organism_setup()
    
    def set_method(self):
        """
        Selects the input method.
        """
        
        def f(**kwargs): return []
        
        if hasattr(self.input_method, '__call__'):
            self.inputm = self.input_method
            self.name = self.name or self.input_method.__name__
        elif hasattr(dataio, self.input_method):
            self.inputm = getattr(dataio, self.input_method)
            self.name = self.name or self.inputm.__name__
        elif self.input_is(self.methods, '__contains__'):
            self.inputm = getattr(dataio,
                                  self.methods[self.input_method.lower()])
            self.name = self.name or self.input_method
        else:
            self.inputm = f
            self.name = self.name or 'Unknown'
    
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
        
        if 'strict' not in self.inputargs:
            self.inputargs['strict'] = False
        
        if self.inputargs['organism'] in common.taxids:
            self.inputargs['organism'] = (
                common.taxids[self.inputargs['organism']]
            )
        
        self.inputargs['mapper'] = self.mapper
    
    def _phosphoelm_setup(self):
        
        if self.ncbi_tax_id != 9606 and 'ltp_only' not in self.inputargs:
            
            self.inputargs['ltp_only'] = False
    
    def _setup(self):
        
        setupmethod = '_%s_setup' % self.input_method.lower()
        
        self._organism_setup()
        
        if hasattr(self, setupmethod):
            
            getattr(self, setupmethod)()
        
        # database specific id conversions
        if self.input_is(self.enzyme_id_uniprot, '__contains__'):
            self.enzyme_id_type = 'uniprot'
    
    def _organism_setup(self):
        
        if self.input_is(self.organisms_supported, '__contains__'):
            
            if self.ncbi_tax_id in common.taxa:
                self.ncbi_tax_id = common.taxa[self.ncbi_tax_id]
            
            self.inputargs['organism'] = self.ncbi_tax_id
        
        self.load_proteome(self.ncbi_tax_id, False)
    
    def _process(self, p):
        
        # human leukocyte antigenes result a result an
        # extremely high number of combinations
        if not p['kinase'] or p['substrate'].startswith('HLA'):
            return
        
        if not isinstance(p['kinase'], list):
            p['kinase'] = [p['kinase']]
        
        kinase_ups = self.mapper.map_names(p['kinase'],
                        self.enzyme_id_type,
                        'uniprot',
                        ncbi_tax_id = self.ncbi_tax_id)
        
        substrate_ups_all = set([])
        
        for sub_id_type in (
            self.substrate_id_types[self.input_method.lower()]
            if self.input_is(self.substrate_id_types, '__contains__')
            else [self.substrate_id_type]
        ):
            
            if type(sub_id_type) is tuple:
                sub_id_type, sub_id_attr = sub_id_type
            else:
                sub_id_attr = 'substrate'
            
            substrate_ups_all.update(
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
            
            se = self.get_seq(s)
            
            if se is None:
                continue
            
            for isof in se.isoforms():
                
                if p['instance'] is not None:
                    
                    if se.match(
                        p['instance'],
                        p['start'],
                        p['end'],
                        isoform=isof
                    ):
                        
                        substrate_ups.append((s, isof))
                    
                else:
                    
                    if se.match(
                        p['resaa'],
                        p['resnum'],
                        isoform=isof
                    ):
                        
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
                    
                    se = self.get_seq(s[0])
                    
                    if se is None:
                        continue
                    
                    nomatch.append((s[0], s[1],
                        ((p['substrate_refseq']
                        if 'substrate_refseq' in p
                        else ''),
                        s, p['instance'],
                        se.get(
                            p['start'],
                            p['end'])
                        )
                    ))
        
        # adding kinase-substrate interactions
        
        for k in kinase_ups:
            
            for s in substrate_ups:
                
                if (
                    not self.allow_mixed_organisms and (
                        self.get_taxon(k) != self.ncbi_tax_id or
                        self.get_taxon(s[0]) != self.ncbi_tax_id
                    )
                ):
                    continue
                
                se = self.get_seq(s[0])
                
                if se is None:
                    continue
                
                res = intera.Residue(
                    p['resnum'],
                    p['resaa'],
                    s[0],
                    isoform=s[1])
                
                if p['instance'] is None:
                    
                    reg = se.get_region(
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
                                    source=[self.name],
                                    isoform=s[1])
                
                dom = intera.Domain(protein=k)
                
                if 'references' not in p:
                    p['references'] = []
                    
                dommot = intera.DomainMotif(
                    domain=dom,
                    ptm=ptm,
                    sources=[self.name],
                    refs=p['references'])
                
                if self.input_is('mimp'):
                    dommot.mimp_sources = ';'.split(p[
                        'databases'])
                    dommot.npmid = p['npmid']
                    
                elif self.input_is('phosphonetworks'):
                    dommot.pnetw_score = p['score']
                    
                elif self.input_is('dbptm'):
                    dommot.dbptm_sources = [p['source']]
                    
                yield dommot
    
    def input_is(self, i, op = '__eq__'):
        
        return (
            type(self.input_method) in common.charTypes and
            getattr(i, op)(self.input_method.lower())
        )
    
    def __iter__(self):
        """
        Iterates through the enzyme-substrate interactions.
        """
        #prg = progress.Progress(len(self.data), 'Processing PTMs', 1)
        for p in self.data:
            
            #prg.step()
            
            for ptm in self._process(p):
                
                yield ptm
        
        #prg.terminate()


class PtmHomologyProcessor(
        homology.PtmHomology,
        PtmProcessor):
    
    def __init__(self,
        input_method,
        ncbi_tax_id,
        map_by_homology_from = [9606],
        trace = False,
        mapper = None,
        enzyme_id_type = 'genesymbol',
        substrate_id_type = 'genesymbol',
        name = None,
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
        
        :param str input_method: Data source for `PtmProcessor`.
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
        
        self.target_taxon = ncbi_tax_id
        self.input_method = input_method
        self.trace = trace
        self.enzyme_id_type = enzyme_id_type
        self.substrate_id_type = substrate_id_type
        self.name = name
        self.ptmprocargs = kwargs
        
        homology.PtmHomology.__init__(self, target = ncbi_tax_id,
                                        only_swissprot = homology_only_swissprot,
                                        strict = ptm_homology_strict,
                                        mapper = mapper)
    
    def __iter__(self):
        """
        Iterates through enzyme-substrate interactions
        translated to another organism by orthology.
        """
        
        for source_taxon in self.map_by_homology_from:
            
            self.set_default_source(source_taxon)
            
            PtmProcessor.__init__(self, self.input_method, source_taxon,
                              trace = self.trace, mapper = self.mapper,
                              enzyme_id_type = self.enzyme_id_type,
                              substrate_id_type = self.substrate_id_type,
                              name = self.name, allow_mixed_organisms = True,
                              **self.ptmprocargs)
            
            #self.reset_ptmprocessor(ncbi_tax_id = source_taxon)
            
            for ptm in PtmProcessor.__iter__(self):
                
                for tptm in self.translate(ptm):
                    
                    yield tptm
    
class PtmAggregator(object):
    
    def __init__(self,
        input_methods = None,
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
        
        self.builtin_inputs = ['PhosphoSite', 'phosphoELM',
                               'Signor', 'dbPTM', 'HPRD',
                               'Li2012', 'PhosphoNetworks',
                               'MIMP']
        
        for k, v in iteritems(locals()):
            setattr(self, k, v)
        
        self.set_inputs()
        
        self.init_mapper()
        
        self.map_by_homology_from = set(self.map_by_homology_from)
        self.map_by_homology_from.discard(self.ncbi_tax_id)
        
        self.build_list()
        self.unique()
    
    def __iter__(self):
        
        for ptm in itertools.chain(*self.full_list.values()):
            
            yield ptm
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def set_inputs(self):
        
        if self.input_methods is None:
            self.input_methods = self.builtin_inputs
    
    def build_list(self):
        """
        Builds a full list of enzyme-substrate interactions from
        all the requested sources. This list might contain redundant
        elements which later will be merged by `unique`.
        This 'full list' is organised into a dict by pairs of proteins
        in order to make it more efficient to compile a unique set
        for each pair.
        """
        
        def extend_lists(ptms):
            
            for ptm in ptms:
                
                key = (ptm.domain.protein, ptm.ptm.protein)
                
                if key not in self.full_list:
                    
                    self.full_list[key] = []
                
                self.full_list[key].append(ptm)
        
        self.full_list = {}
        
        for input_method in self.input_methods:
            
            inputargs = (
                self.inputargs[input_method]
                if input_method in self.inputargs
                else {}
            )
            
            if self.ncbi_tax_id == 9606 or self.nonhuman_direct_lookup:
                
                proc = PtmProcessor(input_method = input_method,
                                    ncbi_tax_id = self.ncbi_tax_id,
                                    trace = self.trace,
                                    mapper = self.mapper,
                                    enzyme_id_type = self.enzyme_id_type,
                                    substrate_id_type = self.substrate_id_type,
                                    **inputargs)
                
                extend_lists(proc.__iter__())
            
            if self.map_by_homology_from:
                
                proc = PtmHomologyProcessor(
                    input_method = input_method,
                    ncbi_tax_id = self.ncbi_tax_id,
                    map_by_homology_from = self.map_by_homology_from,
                    trace = self.trace,
                    mapper = self.mapper,
                    enzyme_id_type = self.enzyme_id_type,
                    substrate_id_type = self.substrate_id_type,
                    homology_only_swissprot = self.homology_only_swissprot,
                    ptm_homology_strict = self.ptm_homology_strict,
                    **inputargs
                )
                
                extend_lists(proc.__iter__())
    
    def unique(self):
        """
        Merges the redundant elements of the interaction list.
        Elements are redundant if they agree in all their attributes
        except the sources, references and isoforms.
        """
        
        self.unique_list = set([])
        
        for key, ptms in iteritems(self.full_list):
            
            self.full_list[key] = self.uniq_ptms(ptms)
    
    @staticmethod
    def uniq_ptms(ptms):
        
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
    
    def export_table(self, fname):
        
        hdr = ['enzyme', 'substrate', 'isoforms',
               'residue', 'offset', 'modification',
               'sources', 'references']
        
        with open(fname, 'w') as fp:
            
            fp.write('%s\n' % '\t'.join(hdr))
            
            for dm in self:
                
                fp.write('%s\n' % '\t'.join(dm.get_line()))
    
    def assign_to_network(self, pa):
        """
        Assigns enzyme-substrate interactions to edges of a
        network in a `pypath.main.PyPath` instance.
        """
        
        pa.update_vname()
        
        if 'ptm' not in pa.graph.es.attributes():
            pa.graph.es['ptm'] = [[] for _ in pa.graph.es]
        
        for key, ptms in iteritems(self.full_list):
            
            nodes = pa.get_node_pair(key[0], key[1],
                    directed = pa.graph.is_directed())
            
            e = None
            
            if nodes:
                e = pa.graph.get_eid(
                    nodes[0], nodes[1], error = False)
            
            if isinstance(e, int) and e > 0:
                
                if pa.graph.es[e]['ptm'] is None:
                    pa.graph.es[e]['ptm'] = []
                
                pa.graph.es[e]['ptm'].extend(ptms)
