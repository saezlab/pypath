#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
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
from past.builtins import xrange, range

import sys
import importlib as imp
import itertools
import collections
import pickle

import pandas as pd

import pypath.inputs.main as dataio
import pypath.share.common as common
import pypath.utils.mapping as mapping
import pypath.utils.homology as homology
import pypath.inputs.uniprot as uniprot_input
import pypath.internals.intera as intera
import pypath.share.progress as progress
import pypath.share.session as session_mod
import pypath.utils.taxonomy as taxonomy
import pypath.inputs as inputs


builtin_inputs = [
    'PhosphoSite',
    'phosphoELM',
    'SIGNOR',
    'dbPTM',
    'HPRD',
    'Li2012',
    'PhosphoNetworks',
    'MIMP',
    'DEPOD',
    'ProtMapper',
    'KEA',
]


class EnzymeSubstrateProcessor(
        homology.Proteomes,
        homology.SequenceContainer
    ):

    methods = {
        'signor': 'load_signor_ptms',
        'mimp': 'get_mimp',
        'phosphonetworks': 'get_phosphonetworks',
        'phosphoelm': 'get_phosphoelm',
        'dbptm': 'get_dbptm',
        'phosphosite': 'get_psite_phos',
        'hprd': 'get_hprd_ptms',
        'li2012': 'li2012_phospho',
        'depod': 'get_depod',
        'protmapper': 'protmapper_ptms',
        'kea': 'kea.kea_enzyme_substrate',
    }

    organisms_supported = set([
        'signor',
        'phosphosite',
        'phosphoelm',
        'dbptm',
        'depod',
    ])

    enzyme_id_uniprot = set([
        'phosphosite',
        'phosphoelm',
        'signor',
        'depod',
        'protmapper',
        'kea',
    ])

    substrate_id_types = {
        'mimp': [('genesymbol', 'substrate'), ('refseq', 'substrate_refseq')],
        'phosphonetworks': ['genesymbol'],
        'phosphoelm': ['uniprot'],
        'li2012': ['genesymbol'],
        'dbptm': ['uniprot'],
        'phosphosite': ['uniprot'],
        'signor': ['uniprot'],
        'hprd': [('refseqp', 'substrate_refseqp')],
        'depod': ['uniprot'],
        'protmapper': ['uniprot'],
        'kea': ['uniprot'],
    }

    resource_names = dict(
        (
            name.lower(),
            name
        )
        for name in builtin_inputs
    )


    def __init__(
            self,
            input_method,
            ncbi_tax_id = 9606,
            trace = False,
            enzyme_id_type = 'genesymbol',
            substrate_id_type = 'genesymbol',
            name = None,
            allow_mixed_organisms = False,
            **kwargs
        ):
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

        self.mammal_taxa = {9606, 10090, 10116}
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

        self.enzyme_id_type = enzyme_id_type
        self.set_method()
        self.set_inputargs(**kwargs)
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

        def empty_input(*args, **kwargs): return []

        
        # a method provided
        if hasattr(self.input_method, '__call__'):
            
            self.inputm = self.input_method
            self.name = self.name or self.input_method.__name__
            
        # the method is associated to a resource name
        # in the list of built in resources
        elif self.input_is(self.methods, '__contains__'):
            
            self.inputm = inputs.get_method(
                self.methods[self.input_method.lower()]
            )
            self.name = (
                self.name or
                (
                    self.resource_names[self.input_method.lower()]
                    if self.input_method.lower() in self.resource_names else
                    self.input_method
                )
            )
            
        # attempting to look up the method in the inputs module
        else:
            
            self.inputm = inputs.get_method(self.input_method) or empty_input
            self.name = self.name or self.inputm.__name__


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


    def _phosphosite_setup(self):

        if 'strict' not in self.inputargs:
            self.inputargs['strict'] = False

        if self.inputargs['organism'] in taxonomy.taxids:
            self.inputargs['organism'] = (
                taxonomy.taxids[self.inputargs['organism']]
            )


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

            if self.ncbi_tax_id in taxonomy.taxa:
                self.ncbi_tax_id = taxonomy.taxa[self.ncbi_tax_id]

            self.inputargs['organism'] = self.ncbi_tax_id

        self.load_proteome(self.ncbi_tax_id, False)


    def _process(self, p):

        # human leukocyte antigenes result a result an
        # extremely high number of combinations
        if (
            not p['kinase'] or (
                isinstance(p['substrate'], common.basestring) and
                p['substrate'].startswith('HLA')
            )
        ):

            return

        if not isinstance(p['kinase'], list):
            p['kinase'] = [p['kinase']]

        kinase_ups = mapping.map_names(
            p['kinase'],
            self.enzyme_id_type,
            'uniprot',
            ncbi_tax_id = self.ncbi_tax_id,
        )

        substrate_ups_all = set([])

        for sub_id_type in (
            self.substrate_id_types[self.input_method.lower()]
            if self.input_is(self.substrate_id_types, '__contains__')
            else [self.substrate_id_type]
        ):

            if isinstance(sub_id_type, (list, tuple)):
                sub_id_type, sub_id_attr = sub_id_type
            else:
                sub_id_attr = 'substrate'

            substrate_ups_all.update(
                set(
                    mapping.map_name(
                        p[sub_id_attr],
                        sub_id_type,
                        'uniprot',
                        self.ncbi_tax_id,
                    )
                )
            )

        # looking up sequences in all isoforms:
        substrate_ups = []

        for s in substrate_ups_all:

            if 'substrate_isoform' in p and p['substrate_isoform']:

                substrate_ups.append((s, p['substrate_isoform']))

            else:

                se = self.get_seq(s)

                if se is None:
                    continue

                for isof in se.isoforms():

                    if 'instance' in p and p['instance'] is not None:

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
                    isoform=s[1],
                )

                if 'instance' not in p or p['instance'] is None:

                    reg = se.get_region(
                        p['resnum'],
                        p['start'] if 'start' in p else None,
                        p['end'] if 'end' in p else None,
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

                ptm = intera.Ptm(
                    s[0],
                    motif=mot,
                    residue=res,
                    typ=p['typ'],
                    source=[self.name],
                    isoform=s[1],
                )

                dom = intera.Domain(protein=k)

                if 'references' not in p:
                    p['references'] = []

                dommot = intera.DomainMotif(
                    domain=dom,
                    ptm=ptm,
                    sources=[self.name],
                    refs=p['references'],
                )
                

                if self.input_is('mimp') and p['databases']:
                    dommot.mimp_sources = p['databases'].split(';')
                    dommot.add_sources(dommot.mimp_sources)
                    dommot.npmid = p['npmid']
                
                if self.input_is('protmapper') and p['databases']:
                    dommot.protmapper_sources = p['databases']
                    dommot.add_sources(p['databases'])

                elif self.input_is('phosphonetworks'):
                    dommot.pnetw_score = p['score']

                elif self.input_is('dbptm') and p['source']:
                    dommot.dbptm_sources = ['%s_dbPTM' % p['source']]
                    dommot.add_sources(dommot.dbptm_sources)

                yield dommot

    def input_is(self, i, op = '__eq__'):

        return (
            type(self.input_method) in common.char_types and
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
    
    
    def __len__(self):
        
        return len(self.data)
    
    
    def __repr__(self):
        
        return '<Enzyme-substrate processor: %u records>' % len(self)


class EnzymeSubstrateHomologyProcessor(
        homology.EnzymeSubstrateHomology,
        EnzymeSubstrateProcessor
    ):


    def __init__(self,
            input_method,
            ncbi_tax_id,
            map_by_homology_from = [9606],
            trace = False,
            enzyme_id_type = 'genesymbol',
            substrate_id_type = 'genesymbol',
            name = None,
            homology_only_swissprot = True,
            ptm_homology_strict = False,
            **kwargs
        ):
        """
        Unifies a `pypath.core.ptm.EnzymeSubstrateProcessor` and
        a `pypath.utils.homology.EnzymeSubstrateHomology` object to build
        a set of enzyme-substrate interactions from a database and
        subsequently translate them by homology to one different organism.
        Multiple organism can be chosen as the source of the
        enzyme-substrate interactions. For example if you want mouse
        interactions, you can translate them from human and from rat.
        To get the original mouse interactions themselves, use an
        other instance of the `EnzymeSubstrateProcessor`.
        To have both the original and the homology translated set,
        and also from multiple databases, whatmore all these merged
        into a single set, use the `EnzymeSubstrateAggregator`.

        :param str input_method: Data source for `EnzymeSubstrateProcessor`.
        :param int ncbi_tax_id: The NCBI Taxonomy ID the interactions
                                should be translated to.
        :param bool homology_only_swissprot: Use only SwissProt
                                             (i.e. not Trembl) at homology
                                             translation.
        :param bool ptm_homology_strict: Use only those homologous PTM pairs
                                         which are in PhosphoSite data, i.e.
                                         do not look for residues with same
                                         offset in protein sequence.

        See further options at `EnzymeSubstrateProcessor`.

        """

        self.map_by_homology_from = map_by_homology_from

        self.target_taxon = ncbi_tax_id
        self.input_method = input_method
        self.trace = trace
        self.enzyme_id_type = enzyme_id_type
        self.substrate_id_type = substrate_id_type
        self.name = name
        self.ptmprocargs = kwargs

        homology.EnzymeSubstrateHomology.__init__(
            self,
            target = ncbi_tax_id,
            only_swissprot = homology_only_swissprot,
            strict = ptm_homology_strict,
        )


    def __iter__(self):
        """
        Iterates through enzyme-substrate interactions
        translated to another organism by orthology.
        """

        for source_taxon in self.map_by_homology_from:

            self.set_default_source(source_taxon)

            EnzymeSubstrateProcessor.__init__(
                self,
                self.input_method,
                source_taxon,
                trace = self.trace,
                enzyme_id_type = self.enzyme_id_type,
                substrate_id_type = self.substrate_id_type,
                name = self.name, allow_mixed_organisms = True,
                **self.ptmprocargs,
            )

            #self.reset_ptmprocessor(ncbi_tax_id = source_taxon)

            for ptm in EnzymeSubstrateProcessor.__iter__(self):

                for tptm in self.translate(ptm):

                    yield tptm


class EnzymeSubstrateAggregator(session_mod.Logger):
    
    
    def __init__(self,
            input_methods = None,
            ncbi_tax_id = 9606,
            map_by_homology_from = None,
            trace = False,
            enzyme_id_type = 'genesymbol',
            substrate_id_type = 'genesymbol',
            homology_only_swissprot = True,
            ptm_homology_strict = False,
            nonhuman_direct_lookup = True,
            inputargs = None,
            pickle_file = None,
        ):
        """
        Docs not written yet.
        """
        
        session_mod.Logger.__init__(self, name = 'enz_sub')
        
        for k, v in iteritems(locals()):
            setattr(self, k, v)
        
        self.main()
    
    
    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def main(self):
        
        if self.pickle_file:
            
            self.load_from_pickle(pickle_file = self.pickle_file)
            
        else:
            
            self.build()
    
    
    def load_from_pickle(self, pickle_file = None):
        
        self._log('Loading from file `%s`.' % pickle_file)
        
        with open(self.pickle_file, 'rb') as fp:
            
            self.enz_sub, self.references = pickle.load(fp)
        
        self.update_ptm_lookup_dict()
    
    
    def save_to_pickle(self, pickle_file):
        
        self._log('Saving to file file `%s`.' % pickle_file)
        
        with open(pickle_file, 'wb') as fp:
            
            pickle.dump(
                obj = (
                    self.enz_sub,
                    self.references,
                ),
                file = fp,
            )
    
    
    def build(self):
        
        self.inputargs = self.inputargs or {}
        self.map_by_homology_from = set(self.map_by_homology_from or [9606])

        self.set_inputs()

        self.map_by_homology_from = set(self.map_by_homology_from)
        self.map_by_homology_from.discard(self.ncbi_tax_id)

        self.build_list()
        self.unique()


    def __iter__(self):

        for ptm in itertools.chain(*self.enz_sub.values()):

            yield ptm
    
    
    def __len__(self):
        
        return sum([len(esub) for esub in self.enz_sub.values()])
    
    
    def __repr__(self):
        
        return '<Enzyme-substrate database: %s relationships>' % len(self)
    
    
    def set_inputs(self):

        if self.input_methods is None:
            
            self.input_methods = builtin_inputs


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

                if key not in self.enz_sub:

                    self.enz_sub[key] = []

                self.enz_sub[key].append(ptm)
                
                for resource in ptm.sources:
                    
                    self.references[resource][ptm.key()].update(ptm.refs)
        
        
        self.enz_sub = {}
        self.references = collections.defaultdict(
            lambda: collections.defaultdict(set)
        )

        for input_method in self.input_methods:
            
            self._log(
                'Loding enzyme-substrate interactions '
                'from `%s`.' % input_method
            )

            inputargs = (
                self.inputargs[input_method]
                if input_method in self.inputargs
                else {}
            )

            if self.ncbi_tax_id == 9606 or self.nonhuman_direct_lookup:
                
                self._log(
                    'Loading enzyme-substrate interactions '
                    'for taxon `%u`.' % self.ncbi_tax_id
                )
                
                proc = EnzymeSubstrateProcessor(
                    input_method = input_method,
                    ncbi_tax_id = self.ncbi_tax_id,
                    trace = self.trace,
                    enzyme_id_type = self.enzyme_id_type,
                    substrate_id_type = self.substrate_id_type,
                    **inputargs,
                )
                
                extend_lists(proc.__iter__())
            
            if self.map_by_homology_from:
                
                self._log(
                    'Mapping `%s` by homology from taxons %s to %u.' % (
                        input_method,
                        ', '.join(
                            '%u' % tax for tax in self.map_by_homology_from
                        ),
                        self.ncbi_tax_id,
                    )
                )
                
                proc = EnzymeSubstrateHomologyProcessor(
                    input_method = input_method,
                    ncbi_tax_id = self.ncbi_tax_id,
                    map_by_homology_from = self.map_by_homology_from,
                    trace = self.trace,
                    enzyme_id_type = self.enzyme_id_type,
                    substrate_id_type = self.substrate_id_type,
                    homology_only_swissprot = self.homology_only_swissprot,
                    ptm_homology_strict = self.ptm_homology_strict,
                    **inputargs
                )
                
                extend_lists(proc.__iter__())
        
        self.references = dict(self.references)
        self.update_ptm_lookup_dict()
    
    
    def update_ptm_lookup_dict(self):
        
        self.ptm_to_enzyme = collections.defaultdict(set)
        self.ptms = {}
        
        for (enz, sub), ptms in iteritems(self.enz_sub):
            
            for ptm in ptms:
                
                self.ptm_to_enzyme[ptm.ptm].add(enz)
                self.ptms[ptm.ptm] = ptm.ptm
        
        self.ptm_to_enzyme = dict(self.ptm_to_enzyme)


    def unique(self):
        """
        Merges the redundant elements of the interaction list.
        Elements are redundant if they agree in all their attributes
        except the sources, references and isoforms.
        """

        self.unique_list = set([])

        for key, ptms in iteritems(self.enz_sub):

            self.enz_sub[key] = self.uniq_ptms(ptms)


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


    def make_df(self, tax_id = False):
        
        self._log('Creating enzyme-substrate interaction data frame.')
        
        
        hdr = ['enzyme', 'substrate', 'isoforms',
               'residue_type', 'residue_offset', 'modification',
               'sources', 'references']

        self.df = pd.DataFrame(
            [dm.get_line() for dm in self],
            columns = hdr
        ).astype(
            {
                'enzyme': 'category',
                'substrate': 'category',
                'isoforms': 'category',
                'residue_type': 'category',
                'residue_offset': 'int32',
                'modification': 'category',
            }
        )

        self.df['enzyme_genesymbol'] = pd.Series([
            (
                mapping.map_name0(
                    u,
                    id_type = 'uniprot',
                    target_id_type = 'genesymbol',
                    ncbi_tax_id = self.ncbi_tax_id,
                ) or ''
            )
            for u in self.df.enzyme
        ])

        self.df['substrate_genesymbol'] = pd.Series([
            (
                mapping.map_name0(
                    u,
                    id_type = 'uniprot',
                    target_id_type = 'genesymbol',
                    ncbi_tax_id = self.ncbi_tax_id,
                ) or ''
            )
            for u in self.df.substrate
        ])

        hdr.insert(2, 'enzyme_genesymbol')
        hdr.insert(3, 'substrate_genesymbol')

        self.df = self.df.loc[:,hdr]

        if tax_id:

            self.df['ncbi_tax_id'] = [self.ncbi_tax_id] * self.df.shape[0]
        
        self._log(
            'Created enzyme-substrate interaction data frame. '
            'Memory usage: %s.' % common.df_memory_usage(self.df)
        )


    def export_table(self, fname):

        self.make_df()
        self.df.to_csv(fname, sep = '\t', index = False)


    def assign_to_network(self, pa):
        """
        Assigns enzyme-substrate interactions to edges of a
        network in a `pypath.main.PyPath` instance.
        """

        pa.update_vname()

        if 'ptm' not in pa.graph.es.attributes():
            pa.graph.es['ptm'] = [[] for _ in pa.graph.es]

        for key, ptms in iteritems(self.enz_sub):

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
    
    
    @property
    def resources(self):
        
        return set.union(*(es.sources for es in self))
    
    
    def update_summaries(self):
        
        self.summaries = {}
        
        refs_by_resource = dict(
            (
                resource,
                set.union(
                    *itertools.chain(
                        self.references[resource].values()
                    )
                )
            )
            for resource in self.resources
        )
        curation_effort_by_resource = dict(
            (
                resource,
                {
                    key + (ref,)
                    for key, refs in
                    itertools.chain(
                        iteritems(self.references[resource])
                    )
                    for ref in refs
                }
            )
            for resource in self.resources
        )
        
        for resource in sorted(self.resources):
            
            n_total = sum(1 for es in self if resource in es.sources)
            n_unique = sum(
                1 for es in self
                if len(es.sources) == 1 and resource in es.sources
            )
            n_shared = sum(
                1 for es in self
                if len(es.sources) > 1 and resource in es.sources
            )
            
            curation_effort = len(curation_effort_by_resource[resource])
            ce_others = set.union(*(
                ce
                for res, ce in iteritems(curation_effort_by_resource)
                if res != resource
            ))
            curation_effort_shared = len(
                curation_effort_by_resource[resource] &
                ce_others
            )
            curation_effort_unique = len(
                curation_effort_by_resource[resource] -
                ce_others
            )
            
            references = len(refs_by_resource[resource])
            refs_others = set.union(*(
                refs
                for res, refs in iteritems(refs_by_resource)
                if res != resource
            ))
            references_shared = len(refs_by_resource[resource] & refs_others)
            references_unique = len(refs_by_resource[resource] - refs_others)
            
            enzymes = len(set(
                es.domain.protein
                for es in self
                if resource in es.sources
            ))
            substrates = len(set(
                es.ptm.protein
                for es in self
                if resource in es.sources
            ))
            
            modification_types = ', '.join(
                (
                    '%s (%u)' % (typ, cnt)
                    for typ, cnt in
                    sorted(
                        iteritems(collections.Counter(
                            es.ptm.typ
                            for es in self
                            if resource in es.sources
                        )),
                        key = lambda type_cnt: type_cnt[1],
                        reverse = True,
                    )
                    if typ
                )
            )
            
            self.summaries[resource] = {
                'name': resource,
                'n_es_total': n_total,
                'n_es_unique': n_unique,
                'n_es_shared': n_shared,
                'n_enzymes': enzymes,
                'n_substrates': substrates,
                'references': references,
                'references_unique': references_unique,
                'references_shared': references_shared,
                'curation_effort': curation_effort,
                'curation_effort_unique': curation_effort_shared,
                'curation_effort_shared': curation_effort_shared,
                'modification_types': modification_types,
            }
    
    
    def summaries_tab(self, outfile = None, return_table = False):
        
        columns = (
            ('name', 'Resource'),
            ('n_es_total', 'E-S interactions'),
            ('n_es_shared', 'Shared E-S interactions'),
            ('n_es_unique', 'Unique E-S interactions'),
            ('n_enzymes', 'Enzymes'),
            ('n_substrates', 'Substrates'),
            ('references', 'References'),
            ('reference_shared', 'Shared references'),
            ('references_unique', 'Unique references'),
            ('curation_effort', 'Curation effort'),
            ('curation_effort_shared', 'Shared curation effort'),
            ('curation_effort_unique', 'Unique curation effort'),
            ('modification_types', 'Modification types'),
        )
        
        tab = []
        tab.append([f[1] for f in columns])
        
        tab.extend([
            [
                str(self.summaries[src][f[0]])
                for f in columns
            ]
            for src in sorted(self.summaries.keys())
        ])
        
        if outfile:
            
            with open(outfile, 'w') as fp:
                
                fp.write('\n'.join('\t'.join(row) for row in tab))
        
        if return_table:
            
            return tab



def init_db(**kwargs):
    
    globals()['db'] = EnzymeSubstrateAggregator(**kwargs)


def get_db(**kwargs):
    
    if 'db' not in globals():
        
        init_db(**kwargs)
    
    return globals()['db']
