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

from future.utils import iteritems
from past.builtins import xrange, range

import sys
import importlib as imp
import itertools
import collections
import pickle
import traceback

import pandas as pd

import pypath.share.common as common
import pypath_common._constants as _const
import pypath.utils.mapping as mapping
import pypath.utils.orthology as orthology
import pypath.inputs.uniprot as uniprot_input
import pypath.internals.intera as intera
import pypath.share.progress as progress
import pypath.share.session as session_mod
import pypath.utils.taxonomy as taxonomy
import pypath.inputs as inputs
import pypath.core.evidence as evidence
import pypath.core.entity as entity
import pypath.resources as resources


class EnzymeSubstrateProcessor(
        orthology.Proteomes,
        orthology.SequenceContainer
    ):


    def __init__(
            self,
            input_param = None,
            input_method = None,
            ncbi_tax_id = None,
            trace = False,
            id_type_enzyme = None,
            id_type_substrate = None,
            name = None,
            allow_mixed_organisms = None,
            organisms_supported = False,
            **kwargs
        ):
        """
        Processes enzyme-substrate interaction data from various databases.
        Provides generators to iterate over these interactions.
        For organisms other than human obtains the organism specific
        interactions from databases.

        :param str input_method:
            Either a method name in the ``inputs`` module or a database
            name e.g. `PhosphoSite` or a callable which returns data in
            list of dicts format.
        :param int ncbi_tax_id: NCBI Taxonomy ID used at the database lookups.
        :param bool trace: Keep data about ambiguous ID mappings and PTM data
                           in mismatch with UniProt sequences.
        :param pypath.mapping.Mapper: A `Mapper` instance. If `None` a new
                                      instance will be created.
        :param str id_type_enzyme: The ID type of the enzyme in the database.
        :param str id_type_substrate: The ID type of the substrate in the
                                      database.

        :param bool nonhuman_direct_lookup: Use direct lookup at non-human
                                            target species.
        :param **kwargs: Args to be forwarded to the input method.

        """

        if not hasattr(self, '_logger'):

            session_mod.Logger.__init__(self, name = 'enz_sub')

        self.mammal_taxa = {9606, 10090, 10116}
        self.nomatch = []
        self.kin_ambig = {}
        self.sub_ambig = {}

        self.input_param = input_param
        self.name = name
        self.id_type_enzyme = id_type_enzyme
        self.id_type_substrate = id_type_substrate
        self.allow_mixed_organisms = allow_mixed_organisms
        self.input_method = input_method
        self.trace = trace
        self.ncbi_tax_id = ncbi_tax_id
        self.organisms_supported = organisms_supported

        self.setup()

        orthology.SequenceContainer.__init__(self)
        self.load_seq(self.ncbi_tax_id)

        if self.allow_mixed_organisms:

            for taxon in self.mammal_taxa:

                self.load_seq(taxon = taxon)

        orthology.Proteomes.__init__(self)

        self.set_inputargs(**kwargs)
        self.load_enz_sub()


    def setup(self):

        self.name = self._get_param('name')
        self.id_type_enzyme = self._get_param('id_type_enzyme', 'genesymbol')
        self.id_type_substrate = self._get_param(
            'id_type_substrate',
            'genesymbol',
        )
        self.id_type_substrate = common.to_list(self.id_type_substrate)
        self.ncbi_tax_id = self._get_param('ncbi_tax_id', 9606)
        self.organisms_supported = self._get_param(
            'organisms_supported',
            False,
        )
        self.allow_mixed_organisms = self._get_param(
            'allow_mixed_organisms',
            False,
        )
        self.input_method = self._get_param('input_method')
        self.set_method()


    def _get_param(self, label, default = None):

        return (
            getattr(self, label) or (
                getattr(self.input_param, label)
                    if hasattr(self.input_param, label) else
                default
            )
        )


    def load_enz_sub(self):

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


        # attempting to look up the method in the inputs module
        if not hasattr(self.input_method, '__call__'):

            self.input_method = (
                inputs.get_method(self.input_method) or
                empty_input
            )

        self.name = self.name or self.input_method.__name__


    def set_inputargs(self, **inputargs):
        """
        Sets the arguments to be provided for the input method.
        """

        self.inputargs = inputargs


    def load_data(self):
        """
        Loads the data by the defined input method.
        """

        input_method_name = '%s.%s' % (
            self.input_method.__module__,
            self.input_method.__name__,
        )

        self._log(
            'Calling `%s` with arguments %s.' % (
                input_method_name,
                str(self.inputargs)
            )
        )

        self.data = self.input_method(**self.inputargs)

        self._log(
            'Loaded data by `%s`, resulted %u records.' % (
                input_method_name,
                len(self.data),
            )
        )


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

        setupmethod = '_%s_setup' % self.name.lower()

        self._organism_setup()

        if hasattr(self, setupmethod):

            getattr(self, setupmethod)()


    def _organism_setup(self):

        if self.organisms_supported:

            if self.ncbi_tax_id in taxonomy.taxa:

                self.ncbi_tax_id = taxonomy.taxa[self.ncbi_tax_id]

            self.inputargs['organism'] = self.ncbi_tax_id

        self.load_proteome(self.ncbi_tax_id, False)


    def _process(self, p):

        # human leukocyte antigenes result a result an
        # extremely high number of combinations
        if (
            not p['kinase'] or (
                isinstance(p['substrate'], str) and
                p['substrate'].startswith('HLA')
            )
        ):

            return

        if not isinstance(p['kinase'], list):
            p['kinase'] = [p['kinase']]

        kinase_ups = mapping.map_names(
            p['kinase'],
            self.id_type_enzyme,
            'uniprot',
            ncbi_tax_id = self.ncbi_tax_id,
        )

        substrate_ups_all = set()

        for sub_id_type in self.id_type_substrate:

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
                            isoform = isof,
                        ):

                            substrate_ups.append((s, isof))

                    else:

                        if se.match(
                            p['resaa'],
                            p['resnum'],
                            isoform = isof,
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

                    self.nomatch.append(
                        (
                            s[0],
                            s[1],
                            (
                                p['substrate_refseq']
                                    if 'substrate_refseq' in p else
                                '',
                                s,
                                p['instance'],
                                se.get(
                                    p['start'],
                                    p['end']
                                ),
                            ),
                        )
                    )

        # building objects representing the enzyme-substrate interaction(s)

        if 'typ' not in p:
            p['typ'] = 'phosphorylation'

        _resources = tuple(
            (
                self.input_param.get_via(name)
                    if hasattr(self.input_param, 'get_via') else
                name
            )
            for name in (
                p['databases'] if 'databases' in p else ()
            )
        )
        _resources += (
            (self.name,)
                if isinstance(self.input_param, str) else
            (self.input_param,)
        )

        # collecting the evidences
        evidences = evidence.Evidences(
            evidence.Evidence(
                resource = _res,
                references = p['references'] if 'references' in p else None
            )
            for _res in _resources
        )

        for s in substrate_ups:

            # building the objects representing the substrate
            se = self.get_seq(s[0])

            if se is None:
                continue

            res = intera.Residue(
                p['resnum'],
                p['resaa'],
                s[0],
                isoform = s[1],
                ncbi_tax_id = self.ncbi_tax_id,
            )

            if 'instance' not in p or p['instance'] is None:

                reg = se.get_region(
                    p['resnum'],
                    p['start'] if 'start' in p else None,
                    p['end'] if 'end' in p else None,
                    isoform = s[1],
                )

                if reg is not None:

                    p['start'], p['end'], p['instance'] = reg

            mot = intera.Motif(
                    s[0],
                    p['start'],
                    p['end'],
                    instance = p['instance'],
                    isoform = s[1],
                    ncbi_tax_id = self.ncbi_tax_id,
                )

            ptm = intera.Ptm(
                s[0],
                motif = mot,
                residue = res,
                typ = p['typ'],
                evidences = evidences,
                isoform = s[1],
                ncbi_tax_id = self.ncbi_tax_id,
            )

            for k in kinase_ups:

                if (
                    not self.allow_mixed_organisms and (
                        self.get_taxon(k) != self.ncbi_tax_id or
                        self.get_taxon(s[0]) != self.ncbi_tax_id
                    )
                ):
                    continue

                # the enzyme (kinase)
                dom = intera.Domain(
                    protein = k,
                    ncbi_tax_id = self.ncbi_tax_id,
                )

                dommot = intera.DomainMotif(
                    domain = dom,
                    ptm = ptm,
                    evidences = evidences,
                )

                if hasattr(self.input_param, 'extra_attrs'):

                    for attr, key in iteritems(self.input_param.extra_attrs):

                        if key in p:

                            setattr(dommot, attr, p[key])

                yield dommot


    def input_is(self, i, op = '__eq__'):

        return (
            type(self.name) in _const.CHAR_TYPES and
            getattr(i, op)(self.name.lower())
        )


    def __iter__(self):
        """
        Iterates through the enzyme-substrate interactions.
        """

        for p in self.data:

            for enz_sub in self._process(p):

                yield enz_sub


    def __len__(self):

        return len(self.data) if hasattr(self, 'data') else 0


    def __repr__(self):

        return '<Enzyme-substrate processor: %u records>' % len(self)


class EnzymeSubstrateOrthologyProcessor(
        orthology.PtmOrthology,
        EnzymeSubstrateProcessor,
        session_mod.Logger
    ):


    def __init__(
            self,
            ncbi_tax_id,
            input_param = None,
            input_method = None,
            map_by_orthology_from = None,
            trace = False,
            id_type_enzyme = None,
            id_type_substrate = None,
            name = None,
            orthology_only_swissprot = True,
            ptm_orthology_strict = False,
            **kwargs
        ):
        """
        Unifies a `pypath.core.enz_sub.EnzymeSubstrateProcessor` and
        a `pypath.utils.orthology.PtmOrthology` object to build
        a set of enzyme-substrate interactions from a database and
        subsequently translate them by orthology to one different organism.
        Multiple organism can be chosen as the source of the
        enzyme-substrate interactions. For example if you want mouse
        interactions, you can translate them from human and from rat.
        To get the original mouse interactions themselves, use an
        other instance of the `EnzymeSubstrateProcessor`.
        To have both the original and the orthology translated set,
        and also from multiple databases, whatmore all these merged
        into a single set, use the `EnzymeSubstrateAggregator`.

        :param str input_method: Data source for `EnzymeSubstrateProcessor`.
        :param int ncbi_tax_id: The NCBI Taxonomy ID the interactions
                                should be translated to.
        :param bool orthology_only_swissprot: Use only SwissProt
                                             (i.e. not Trembl) at orthology
                                             translation.
        :param bool ptm_orthology_strict: Use only those homologous PTM pairs
                                         which are in PhosphoSite data, i.e.
                                         do not look for residues with same
                                         offset in protein sequence.

        See further options at `EnzymeSubstrateProcessor`.

        """

        if not hasattr(self, '_logger'):

            session_mod.Logger.__init__(self, name = 'enz_sub_orthology')

        self.target_taxon = ncbi_tax_id
        self.map_by_orthology_from = (
            map_by_orthology_from or
            {9606, 10090, 10116}
        )
        self.map_by_orthology_from = common.to_set(self.map_by_orthology_from)
        self.map_by_orthology_from.discard(self.target_taxon)

        self.input_param = input_param
        self.input_method = input_method
        self.trace = trace
        self.id_type_enzyme = id_type_enzyme
        self.id_type_substrate = id_type_substrate
        self.name = name
        self.ptmprocargs = kwargs

        orthology.PtmOrthology.__init__(
            self,
            target = ncbi_tax_id,
            only_swissprot = orthology_only_swissprot,
            strict = ptm_orthology_strict,
        )


    def __iter__(self):
        """
        Iterates through enzyme-substrate interactions
        translated to another organism by orthology.
        """

        for source_taxon in self.map_by_orthology_from:

            self._log(
                'Translating enzyme-substrate interactions '
                'from organism %u to %u.' % (
                    source_taxon,
                    self.target_taxon,
                )
            )

            self.set_default_source(source_taxon)

            EnzymeSubstrateProcessor.__init__(
                self,
                input_param = self.input_param,
                input_method = self.input_method,
                ncbi_tax_id = source_taxon,
                trace = self.trace,
                id_type_enzyme = self.id_type_enzyme,
                id_type_substrate = self.id_type_substrate,
                name = self.name,
                allow_mixed_organisms = True,
                **self.ptmprocargs,
            )

            self._log(
                'Enzyme-substrate interactions loaded from resource `%s` '
                'for organism %s, %u raw records.' % (
                    self.name,
                    source_taxon,
                    len(self),
                )
            )

            for es in EnzymeSubstrateProcessor.__iter__(self):

                for target_es in self.translate(es):

                    yield target_es


    def __repr__(self):

        return (
            '<Enzyme-substrate orthology processor, '
            'target taxon: %u, source taxon(s): %s>' % (
                self.target_taxon,
                ', '.join(str(tax) for tax in self.map_by_orthology_from),
            )
        )


class EnzymeSubstrateAggregator(session_mod.Logger):


    def __init__(self,
            input_param = None,
            exclude = None,
            ncbi_tax_id = 9606,
            map_by_orthology_from = None,
            trace = False,
            orthology_only_swissprot = True,
            ptm_orthology_strict = False,
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
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
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
        self.map_by_orthology_from = (
            (
                {9606, 10090, 10116}
                    if self.ncbi_tax_id != 9606 else
                set()
            )
                if self.map_by_orthology_from is None else
            self.map_by_orthology_from
        )
        self.map_by_orthology_from = set(self.map_by_orthology_from)
        self.map_by_orthology_from.discard(self.ncbi_tax_id)

        self.set_inputs()

        self.build_list()
        self.unique()


    def __iter__(self):

        for ptm in itertools.chain(*self.enz_sub.values()):

            yield ptm


    def __len__(self):

        return sum([len(esub) for esub in self.enz_sub.values()])


    def __repr__(self):

        return '<Enzyme-substrate database: %s relationships>' % len(self)


    def __getitem__(self, *args):

        args = args[0] if isinstance(args[0], tuple) else args

        return self.get_enzyme_substrate(*args)


    def get_enzyme_substrate(self, enzyme, substrate):

        enzyme = entity.Entity(enzyme)
        substrate = entity.Entity(substrate)

        key = (enzyme, substrate)

        if key in self.enz_sub:

            return self.enz_sub[key]


    def set_inputs(self):

        self.input_param = (
            self.input_param or
            resources.get_controller().collect_enzyme_substrate()
        )


    def build_list(self):
        """
        Builds a full list of enzyme-substrate interactions from
        all the requested sources. This list might contain redundant
        elements which later will be merged by `unique`.
        This 'full list' is organised into a dict by pairs of proteins
        in order to make it more efficient to compile a unique set
        for each pair.
        """

        def extend_lists(enz_sub):

            for es in enz_sub:

                key = (es.domain.protein, es.ptm.protein)

                if key not in self.enz_sub:

                    self.enz_sub[key] = []

                self.enz_sub[key].append(es)

                for ev in es.evidences:

                    resource_key = (ev.resource.name, ev.resource.via)

                    self.references[resource_key][es.key()].update(
                        ev.references
                    )

        self._log(
            'Starting to build enzyme-substrate '
            'database for organism `%u`.' % self.ncbi_tax_id
        )

        self.enz_sub = {}
        self.references = collections.defaultdict(
            lambda: collections.defaultdict(set)
        )

        for input_param in self.input_param:

            name = (
                input_param['name']
                    if isinstance(input_param, dict) else
                input_param.name
            )

            try:

                input_method = (
                    input_param['input_method']
                        if isinstance(input_param, dict) else
                    input_param.input_method
                )

                self._log(
                    'Loading enzyme-substrate interactions '
                    'from resource `%s` by method `%s`.' % (
                        name,
                        input_method,
                    )
                )

                args = (
                    input_param
                        if isinstance(input_param, dict) else
                    {'input_param': input_param}
                )

                if (
                    self.ncbi_tax_id == 9606 or (
                        self.nonhuman_direct_lookup and
                        input_param.organisms_supported
                    )
                ):

                    self._log(
                        'Loading enzyme-substrate interactions '
                        'for taxon `%u`.' % self.ncbi_tax_id
                    )

                    proc = EnzymeSubstrateProcessor(
                        ncbi_tax_id = self.ncbi_tax_id,
                        trace = self.trace,
                        **args,
                    )

                    extend_lists(proc.__iter__())

                if self.map_by_orthology_from:

                    source_taxons_str = ', '.join(
                        '%u' % tax for tax in self.map_by_orthology_from
                    )

                    self._log(
                        'Mapping `%s` by orthology from taxons %s to %u.' % (
                            input_method,
                            source_taxons_str,
                            self.ncbi_tax_id,
                        )
                    )

                    proc = EnzymeSubstrateOrthologyProcessor(
                        ncbi_tax_id = self.ncbi_tax_id,
                        map_by_orthology_from = self.map_by_orthology_from,
                        trace = self.trace,
                        orthology_only_swissprot = self.orthology_only_swissprot,
                        ptm_orthology_strict = self.ptm_orthology_strict,
                        **args
                    )

                    extend_lists(proc.__iter__())

                    self._log(
                        'Finished translating `%s` by orthology '
                        'from %s to %u.' % (
                            input_method,
                            source_taxons_str,
                            self.ncbi_tax_id,
                        )
                    )

                self._log(
                    'Finished loading enzyme-substrate data '
                    'from resource `%s`.' % name
                )

            except Exception as e:

                self._log('Failed to load resource `%s`.' % name)
                self._log_traceback()

                try:

                    traceback.print_tb(
                        e.__traceback__,
                        file = sys.stdout,
                    )

                except Exception as e:

                    self._log('Failed handling exception.')
                    self._log_traceback()

        self.references = dict(self.references)
        self.update_ptm_lookup_dict()

        self._log(
            'Finished building enzyme-substrate database '
            'for organism `%u`, resulted %u relationships.' % (
                self.ncbi_tax_id,
                len(self),
            )
        )


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

        self.unique_list = set()

        for key, enz_sub in iteritems(self.enz_sub):

            self.enz_sub[key] = self.uniq_enz_sub(enz_sub)


    @staticmethod
    def uniq_enz_sub(enz_sub):

        enz_sub_uniq = []

        for es in enz_sub:

            merged = False

            for i, es_u in enumerate(enz_sub_uniq):

                if es == es_u:

                    enz_sub_uniq[i].merge(es)
                    merged = True

            if not merged:

                enz_sub_uniq.append(es)

        return enz_sub_uniq


    def make_df(self, tax_id = False, resources_only_primary = False):

        self._log('Creating enzyme-substrate interaction data frame.')


        hdr = [
            'enzyme',
            'enzyme_genesymbol',
            'substrate',
            'substrate_genesymbol',
            'isoforms',
            'residue_type',
            'residue_offset',
            'modification',
            'sources',
            'references',
            'curation_effort',
        ]

        self.df = pd.DataFrame(
            [
                dm.get_line(resources_only_primary = resources_only_primary)
                for dm in self
            ],
            columns = hdr,
        ).astype(
            {
                'enzyme': 'category',
                'substrate': 'category',
                'isoforms': 'category',
                'residue_type': 'category',
                'residue_offset': 'int32',
                'modification': 'category',
                'sources': 'category',
                'references': 'category',
                'curation_effort': 'int32',
            }
        )

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
        Assigns enzyme-substrate interactions to the edges of a
        network in a py:class:``pypath.legacy.main.PyPath`` instance.
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

        return set.union(*(
            es.evidences.get_resource_names_via(via = None)
            for es in self
        ))


    @property
    def resources_sorted(self):

        return sorted(
            self.resources,
            key = lambda res: (res[0], '') if res[1] is None else res
        )


    def update_summaries(self, collect_args = None):

        collect_args = collect_args or {'via': False}

        self.summaries = {}

        resources = [
            res for res in self.resources_sorted
            if (
                res[1] is None or
                'via' not in collect_args or
                collect_args['via'] != False
            )
        ]

        refs_by_resource = dict(
            (
                resource,
                set.union(
                    *itertools.chain(
                        self.references[resource].values()
                    )
                )
            )
            for resource in resources
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
            for resource in resources
        )

        resources_sorted = sorted(resources)

        for resource in resources:

            n_total = sum(
                1
                for es in self
                if resource in es.evidences.get_resource_names(**collect_args)
            )

            n_unique = sum(
                1 for es in self
                if (
                    resource[0] in es.evidences and
                    es.evidences.count_resources(**collect_args) == 1
                )
            )
            n_shared = sum(
                1 for es in self
                if (
                    resource[0] in es.evidences and
                    es.evidences.count_resources(**collect_args) > 1
                )
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
                if resource[0] in es.evidences
            ))
            substrates = len(set(
                es.ptm.protein
                for es in self
                if resource[0] in es.evidences
            ))

            modification_types = ', '.join(
                (
                    '%s (%u)' % (typ, cnt)
                    for typ, cnt in
                    sorted(
                        iteritems(collections.Counter(
                            es.ptm.typ
                            for es in self
                            if resource[0] in es.evidences
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
            ('references_shared', 'Shared references'),
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
            for src in sorted(
                self.summaries.keys(),
                key = lambda res: (res[0], '') if res[1] is None else res,
            )
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
