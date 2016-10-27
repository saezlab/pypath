#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2016 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import os
import sys
import subprocess
import imp
import locale
import numpy as np
import pandas as pd

import pypath.main as main
import pypath.plot as plot
import pypath.data_formats as data_formats
import pypath.dataio as dataio
import pypath.descriptions as descriptions
from pypath.common import *
import pypath.refs as _refs


class Workflow(object):
    def __init__(self,
                 name,
                 network_datasets=[],
                 do_main_table=True,
                 do_compile_main_table=True,
                 do_curation_table=True,
                 do_compile_curation_table=True,
                 do_simgraphs=True,
                 do_multi_barplots=True,
                 do_coverage_groups=True,
                 do_htp_char=True,
                 do_ptms_barplot=True,
                 do_scatterplots=True,
                 do_history_tree=True,
                 do_compile_history_tree=True,
                 do_refs_journals_grid=True,
                 do_refs_years_grid=True,
                 do_dirs_stacked=True,
                 do_refs_composite=True,
                 do_curation_plot=True,
                 do_refs_by_j=True,
                 do_refs_by_db=True,
                 do_refs_by_year=True,
                 do_resource_list=True,
                 do_compile_resource_list=True,
                 do_consistency_dedrogram=True,
                 do_consistency_table=True,
                 title=None,
                 outdir=None,
                 htdata={},
                 inc_raw=None,
                 **kwargs):

        for k, v in iteritems(locals()):
            setattr(self, k, v)

        for k, v in iteritems(kwargs):
            setattr(self, k, v)

        self.title = self.name if self.title is None else self.title

        self.defaults = {
            'ccolors': {
                'p': '#77AADD',
                'm': '#77CCCC',
                'i': '#DDAA77',
                'r': '#CC99BB'
            },
            'ccolors2': {
                'p': '#4477AA',
                'm': '#117777',
                'i': '#774411',
                'r': '#771155'
            },
            'group_colors': ['#4477AA', '#44AAAA', '#DDAA77', '#CC99BB'],
            'group_colors2': ['#77AADD', '#77CCCC', '#DDAA77', '#CC99BB'],
            'table2file': 'curation_tab_%s.tex' % self.name,
            'stable2file': 'curation_tab_stripped_%s.tex' % self.name,
            'latex': '/usr/bin/xelatex',
            'latex_timeout': 10,
            'compile_latex': True,
            'multi_barplots_summary': True,
            'protein_lists': {
                'proteome': 'human proteome',
                'signaling_proteins': 'signaling proteins'
            },
            'intogen_file': None,
            'set_annots': {
                'receptors': 'receptors',
                'tfs': 'transcription factors',
                'kinases': 'kinases',
                'druggability': 'druggable proteins',
                'disease_genes': 'DisGeNet disease related genes'
            },
            'load_annots': {
                'corum': 'CORUM complexes',
                'ptms': 'post-translational modifications',
                'disgenet': 'DisGeNet disease related genes'
            },
            'fiher_file': 'fisher_%s' % self.name,
            'fisher':
            [('dis', 'Disease related genes'), ('rec', 'Receptors'),
             ('tf', 'Transcription factors'), ('kin', 'Kinases'),
             ('dgb', 'Druggable proteins'), ('cdv', 'Cancer drivers'),
             ('sig', 'Signaling proteins')],
            'cat_ordr': [
                'Activity flow', 'Enzyme-substrate', 'Undirected PPI',
                'Process description'
            ],
            'pdf_vcount_order': 'vcount_ordr.pdf',
            'htp_lower': 1,
            'htp_upper': 500,
            'history_tree_fname': 'history_tree.tex',
            'main_table_fname': 'main_table_%s.tex' % self.name,
            'main_table_stripped_fname':
            'main_table_stripped_%s.tex' % self.name,
            'simgraph_vertex_fname':
            'sources_similarity_vertex_%s.pdf' % self.name,
            'simgraph_edge_fname':
            'sources_similarity_edge_%s.pdf' % self.name,
            'simgraph_curation_fname':
            'sources_similarity_curation_%s.pdf' % self.name,
            'refs_journal_grid_fname': 'refs_by_db_journal_%s.pdf' % self.name,
            'refs_year_grid_fname': 'refs_by_db_year_%s.pdf' % self.name,
            'dirs_stacked_fname': 'dirs-signes-by-db-%s_%s.pdf',
            'refs_composite_fname': 'refs-composite_%s.pdf',
            'refs_by_j_file': 'references-by-journal-%u_%s.pdf',
            'refs_by_db_file': 'references-by-db_%s.pdf',
            'refs_by_year_file': 'references-by-year_%s.pdf',
            'curation_plot_fname': 'new-curated-by-year_%s.pdf' % self.name,
            'resource_list_fname': 'resources.tex',
            'resource_list_fname_stripped': 'resources_stripped.tex',
            'consistency_dedrogram_fname':
            'inconsistency-dendrogram_%s.pdf' % self.name,
            'consistency_table_fname': 'inconsistency-table_%s.tex' % self.name
        }

        self.barplots_settings = [
            (
                'Proteins',
                'Number of proteins',
                'proteins',
                lambda gs: gs[0].vcount(),
                lambda gs: len(self.specific(gs[1], gs[0].vs)),
                'vcount'
            ),
            (
                'Interactions',
                'Number of interactions',
                'interactions',
                lambda gs: gs[0].ecount(),
                lambda gs: len(self.specific(gs[1], gs[0].es)),
                'vcount'
            ),
            (
                'Density',
                'Graph density',
                'density',
                lambda gs: gs[0].density(),
                None,
                'vcount'
            ),
            (
                'Transitivity',
                'Graph global transitivity',
                'transitivity',
                lambda gs: gs[0].transitivity_undirected(),
                None,
                'vcount'
            ),
            (
                'Diameter',
                'Graph diameter',
                'diameter',
                lambda gs: gs[0].diameter(),
                None,
                'vcount'
            ),
            (
                'Receptors',
                'Number of receptors',
                'receptors',
                lambda gs: len([v for v in gs[0].vs if v['rec']]),
                lambda gs: len(
                    [v for v in self.specific(gs[1], gs[0].vs) if v['rec']]),
                'vcount'
            ),
            (
                r'Receptors [%]',
                'Percentage of receptors',
                'receptorprop',
                lambda gs: len([v for v in gs[0].vs if v['rec']]) /
                float(gs[0].vcount()) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['rec']]) /
                float(gs[0].vcount()) * 100.0,
                'vcount'
            ),
            (
                r'Receptors [%]',
                'Percentage of all human receptors covered',
                'receptorcov',
                lambda gs: len([v for v in gs[0].vs if v['rec']]) /
                float(len(self.pp.lists['rec'])) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['rec']]) /
                float(len(self.pp.lists['rec'])) * 100.0,
                'vcount'
            ),
            (
                'TFs',
                'Number of transcription factors',
                'tfs',
                lambda gs: len([v for v in gs[0].vs if v['tf']]),
                lambda gs: len(
                    [v for v in self.specific(gs[1], gs[0].vs) if v['tf']]),
                'vcount'
            ),
            (
                r'TFs [%]',
                'Percentage of transcription factors',
                'tfprop',
                lambda gs: len([v for v in gs[0].vs if v['tf']]) /
                float(gs[0].vcount()) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['tf']]) /
                float(gs[0].vcount()) * 100.0,
                'vcount'
            ),
            (
                r'TFs [%]',
                'Percentage of all human TFs covered',
                'tfcov',
                lambda gs: len([v for v in gs[0].vs if v['tf']]) /
                float(len(self.pp.lists['tf'])) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['tf']]) /
                float(len(self.pp.lists['tf'])) * 100.0,
                'vcount'
            ),
            (
                'Kinases',
                'Number of kinases',
                'kinases',
                lambda gs: len([v for v in gs[0].vs if v['kin']]),
                lambda gs: len(
                    [v for v in self.specific(gs[1], gs[0].vs) if v['kin']]),
                'vcount'
            ),
            (
                r'Kinases [%]',
                'Percentage of kinases',
                'kinprop',
                lambda gs: len([v for v in gs[0].vs if v['kin']]) /
                float(gs[0].vcount()) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['kin']]) /
                float(gs[0].vcount()) * 100.0,
                'vcount'
            ),
            (
                r'Kinases [%]',
                'Percentage of all human kinases covered',
                'kincov',
                lambda gs: len([v for v in gs[0].vs if v['kin']]) /
                float(len(self.pp.lists['kin'])) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['kin']]) /
                float(len(self.pp.lists['kin'])) * 100.0,
                'vcount'
            ),
            (
                r'Druggable proteins [%]',
                'Percentage of druggable proteins',
                'dgbprop',
                lambda gs: len([v for v in gs[0].vs if v['dgb']]) /
                float(gs[0].vcount()) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['dgb']]) /
                float(gs[0].vcount()) * 100.0,
                'vcount'
            ),
            (
                r'Druggable proteins [%]',
                'Percentage of all human druggable proteins covered',
                'dgbcov',
                lambda gs: len([v for v in gs[0].vs if v['dgb']]) /
                float(len(self.pp.lists['dgb'])) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['dgb']]) /
                float(len(self.pp.lists['dgb'])) * 100.0,
                'vcount'
            ),
            (
                r'Disease genes [%]',
                'Percentage of disease-gene associations covered',
                'discov',
                lambda gs: len([v for v in gs[0].vs if v['dis']]) /
                float(len(self.pp.lists['dis'])) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['dis']]) /
                float(len(self.pp.lists['dis'])) * 100.0,
                'vcount'
            ),
            (
                'Complexes',
                'Number of complexes covered',
                'complexes',
                lambda gs: len(self.pp.complexes_in_network(graph=gs[0])),
                None,
                'vcount'
            ),
            (
                'E-S interactions',
                'Number of enzyme-substrate interactions covered',
                'ptmnum',
                lambda gs: sum(map(lambda e: len(e['ptm']), gs[0].es)),
                lambda gs: sum(
                    map(lambda e: len(e['ptm']), self.specific(gs[1], gs[0].es))),
                'vcount'
            ),
            (
                'E-S interactions',
                'Number of interactions associated with enzyme-substrate relationship',
                'havingptm',
                lambda gs: sum(map(lambda e: len(e['ptm']) > 0, gs[0].es)),
                lambda gs: sum(
                    map(lambda e: len(e['ptm']) > 0, self.specific(gs[1], gs[0].es))),
                'vcount'
            ),
            (
                'Interactions with \n' + r'E-S relationship [%]',
                'Percentage of interactions associated with enzyme-substrate relationship',
                'ptmprop',
                lambda gs: sum(map(lambda e: len(e['ptm']) > 0, gs[
                               0].es)) / float(gs[0].ecount()) * 100.0,
                lambda gs: sum(map(lambda e: len(e['ptm']) > 0, self.specific(
                    gs[1], gs[0].es))) / float(gs[0].ecount()) * 100.0,
                'vcount'
            ),
            (
                r'CGC genes [%]',
                'Percentage of COSMIC Cancer Gene Census cancer drivers covered',
                'ccgccov',
                lambda gs: len(set(gs[0].vs['name']) & set(self.pp.lists['cgc'])) /
                float(len(self.pp.lists['cgc'])) * 100.0,
                lambda gs: len(
                    set(map(lambda v: v['name'], self.specific(gs[1], gs[0].vs))) &
                    set(self.pp.lists['cgc'])
                ) /
                float(len(self.pp.lists['cgc'])) * 100.0,
                'vcount'
            ),
            (
                r'IntOGen genes [%]',
                'Percentage of IntOGen cancer drivers covered',
                'intocov',
                lambda gs: len(set(gs[0].vs['name']) & set(self.pp.lists['IntOGen'])) /
                float(len(self.pp.lists['IntOGen'])) * 100.0,
                lambda gs: len(
                    set(map(lambda v: v['name'], self.specific(gs[1], gs[0].vs))) &
                    set(self.pp.lists['IntOGen'])
                ) /
                float(len(self.pp.lists['IntOGen'])) * 100.0,
                'vcount'
            ),
            (
                r'Signaling proteins [%]',
                'Percentage of all human signaling proteins covered',
                'sigcov',
                lambda gs: len(set(gs[0].vs['name']) & set(self.pp.lists['sig'])) /
                float(len(self.pp.lists['sig'])) * 100.0,
                lambda gs: len(
                    set(map(lambda v: v['name'], self.specific(gs[1], gs[0].vs))) &
                    set(self.pp.lists['sig'])
                ) /
                float(len(self.pp.lists['sig'])) * 100.0,
                'vcount'
            ),
            (
                r'Signaling proteins [%]',
                'Percentage of signaling proteins',
                'sigpct',
                lambda gs: len(set(gs[0].vs['name']) & set(self.pp.lists['sig'])) /
                float(gs[0].vcount()) * 100.0,
                lambda gs: len(
                    set(map(lambda v: v['name'], self.specific(gs[1], gs[0].vs))) &
                    set(self.pp.lists['sig'])
                ) /
                float(gs[0].vcount()) * 100.0,
                'vcount'
            )
        ]

        self.scatterplots_settings = [
            (
                lambda gs: gs[0].vcount(),
                lambda gs: gs[0].ecount(),
                lambda gs: gs[0].density(),
                {
                    'ylim': [20.0, 200000.0],
                    'xlim': [30.0, 15000.0],
                    'ylog': True,
                    'xlog': True,
                    'xlab': 'Number of proteins',
                    'ylab': 'Number of interacting pairs',
                    'legtitle': 'Density',
                    'title': 'Number of proteins, interactions and graph density'
                },
                'vcount-ecount',
                True
            ),
            (
                lambda gs: len(set(gs[0].vs['name']) &
                               (set(self.pp.lists['dis']))),
                lambda gs: len(set(gs[0].vs['name']) &
                               set(self.pp.lists['cdv'])),
                lambda gs: gs[0].vcount(),
                {
                    'ylim': [0.0, 2000.0],
                    'xlim': [30.0, 5500.0],
                    'ylog': True,
                    'xlog': True,
                    'xlab': 'Number of\ndisease related proteins',
                    'ylab': 'Number of cancer drivers',
                    'legtitle': 'Total number of proteins',
                    'title': 'Number of disease related proteins and cancer drivers',
                    'legstrip': (3, None)
                },
                'dis-cancer',
                False
            ),
            (
                lambda gs: len(set(gs[0].vs['name']) &
                               set(self.pp.lists['rec'])),
                lambda gs: len(set(gs[0].vs['name']) &
                               set(self.pp.lists['tf'])),
                lambda gs: gs[0].vcount(),
                {
                    'ylim': [0.5, 2000.0],
                    'xlim': [5.0, 1200.0],
                    'ylog': True,
                    'xlog': True,
                    'xlab': 'Number of receptors',
                    'ylab': 'Number of transcription factors',
                    'legtitle': 'Total number of proteins',
                    'title': 'Number of receptors and TFs',
                    'legstrip': (3, None)
                },
                'tf-rec',
                False
            ),
            (
                lambda gs: len(self.pp.complexes_in_network(graph=gs[0])),
                lambda gs: sum(map(lambda e: len(e['ptm']) > 0, gs[0].es)),
                lambda gs: gs[0].ecount(),
                {
                    'ylim': [0.5, 5400.0],
                    'xlim': [3.0, 1500.0],
                    'ylog': True,
                    'xlog': True,
                    'xlab': 'Number of complexes',
                    'ylab': 'Number of\nenzyme-substrate relationships',
                    'legtitle': 'Total number of interactions',
                    'title': 'Number of complexes and enzyme-substrate relationships',
                    'legstrip': (3, None)
                },
                'comp-ptm',
                False
            )
        ]

        for k, v in iteritems(self.defaults):
            if not hasattr(self, k):
                setattr(self, k, v)

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def run(self):

        self.load_data()
        self.make_plots()

    def load_data(self):

        self.set_outdir()
        self.init_pypath()
        self.load_protein_lists()
        self.load_annotations()
        self.set_categories()
        self.separate()
        self.barplot_colors()
        if self.do_refs_composite or \
                self.do_refs_journals_grid or \
                self.do_refs_years_grid or \
                self.do_curation_plot or \
                self.do_htp_char:
            self.load_pubmed_data()
        if self.do_dirs_stacked:
            self.get_dirs_data()
        if self.do_consistency_dedrogram or \
                self.do_consistency_table:
            self.inconsistency_data()

    def make_plots(self):

        if self.do_main_table:
            self.main_table()
            if self.do_compile_main_table:
                self.latex_compile(self.main_table_fname)
        if self.do_curation_table:
            self.curation_table()
            if self.do_compile_curation_table:
                self.latex_compile(self.stable2file)
        self.get_multibarplot_ordr()
        if self.do_simgraphs:
            self.make_simgraph_vertex()
            self.make_simgraph_edge()
            self.make_simgraph_curation()
        if self.do_multi_barplots:
            self.make_multi_barplots()
        if self.do_coverage_groups:
            self.get_coverage_groups_data()
            self.make_coverage_groups_plot()
        if self.do_scatterplots:
            self.make_scatterplots()
        if self.do_ptms_barplot:
            self.all_ptms_list()
            self.make_ptms_barplot()
        if self.do_htp_char:
            self.make_htp_characteristics()
        if self.do_history_tree:
            self.make_history_tree()
            if self.do_compile_history_tree:
                self.latex_compile(self.history_tree_fname)
        if self.do_refs_by_j:
            self.make_refs_by_journal()
        if self.do_refs_by_db:
            self.make_refs_by_db()
        if self.do_refs_by_year:
            self.make_refs_by_year()
        if self.do_refs_years_grid:
            self.make_refs_years_grid()
        if self.do_refs_journals_grid:
            self.make_refs_journals_grid()
        if self.do_dirs_stacked:
            self.make_dirs_stacked()
            self.make_dirs_stacked(include_all=False)
        if self.do_refs_composite:
            self.make_refs_composite()
        if self.do_curation_plot:
            self.make_curation_plot()
        if self.do_resource_list:
            self.resource_list_table()
            if self.do_compile_resource_list:
                self.compile_latex(self.resource_list_fname)
        if self.do_consistency_dedrogram:
            self.make_consistency_dendrogram()
        if self.do_consistency_table:
            self.make_consistency_table()

    def set_outdir(self):
        self.outdir = self.name if self.outdir is None else self.outdir
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

    def get_path(self, fname):
        return os.path.join(self.outdir, fname)

    def init_pypath(self):

        self.pp = main.PyPath(9606)
        for netdata in self.network_datasets:
            self.pp.load_resources(getattr(data_formats, netdata))

    def load_protein_lists(self):

        for meth, name in iteritems(self.protein_lists):

            self.console('Loading list of %s' % name)
            getattr(self.pp, '%s_list' % meth)()

        self.console('Loading list of Cancer Gene Census %scancer drivers' %
                     ('and IntOGen ' if self.intogen_file else ''))
        self.pp.cancer_drivers_list(intogen_file=self.intogen_file)

    def load_annotations(self):

        for meth, name in iteritems(self.set_annots):

            self.console('Loading list of %s into protein attribute' % name)
            getattr(self.pp, 'set_%s' % meth)()

        for meth, name in iteritems(self.load_annots):

            self.console('Loading %s into protein attribute' % name)
            getattr(self.pp, 'load_%s' % meth)()

    def set_categories(self):
        self.pp.set_categories()

    def separate(self):
        self.sep = self.pp.separate()
        self.csep = self.pp.separate_by_category()
        self.cats = dict(
            map(lambda c: (c[0], data_formats.catnames[c[1]]),
                iteritems(data_formats.categories)))
        self.cats.update(
            dict(
                map(lambda c: (('All', c[0]), c[1]),
                    filter(lambda c: c[0] in self.pp.has_cats,
                           iteritems(data_formats.catnames)))))
        self.cat_ordr = list(
            filter(lambda c: data_formats.catletters[c] in self.pp.has_cats,
                   self.cat_ordr))

    def fisher_tests(self):

        self.console('Doing Fisher tests, writing results to `%s`' %
                     self.get_path(self.fisher_file))
        with open(self.get_path(self.fisher_file), 'w') as fi:

            for attr, name in self.fisher:
                print(attr, name)
                cont = np.array(
                    [[len(self.pp.lists['proteome']), self.pp.graph.vcount()],
                     [
                         len(self.pp.lists[attr]),
                         len([1 for v in self.pp.graph.vs if v[attr]])
                     ]])
                fi.write('%s:' % name)
                fi.write('\t%s\t%s\n' % stats.fisher_exact(cont))

    def get_data(self, fun, attr):
        if attr is None or not hasattr(self, attr):
            result = \
                list(
                    zip(
                        *[(s, fun((self.sep[s], s))) for s in sorted(self.pp.sources)] +
                        list(
                            map(
                                lambda c:
                                    (('All', c[0]), fun(
                                        (self.csep[c[0]], c[0]))),
                                sorted(self.csep.keys())
                            )
                        )
                    )
                )
            if attr is None:
                return result
            else:
                setattr(self, attr, result)

    def specific(self, s, seq):
        return [
            w for w in seq
            if s in data_formats.categories and s in w['sources'] and len(w[
                'sources'] & getattr(data_formats, data_formats.categories[s]))
            == 1 or s in w['cat'] and len(w['cat']) == 1
        ]

    def make_multi_barplots(self):

        self.multi_barplots = {}

        for par in self.barplots_settings:

            _ = sys.stdout.write('\t:: Plotting %s\n' % par[1])

            data_attr = 'data_%s' % par[2]

            self.get_data(par[3], data_attr)

            if par[4] is not None:
                data_attr2 = 'data_%s' % par[4]
                self.get_data(par[4], data_attr2)

            data = getattr(self, data_attr)
            data2 = getattr(self, data_attr2)

            ordr = self.vcount_ordr if par[5] == 'vcount' else par[5]

            csvname = self.get_path('%s-by-db_%s.csv' % (par[2], self.name))

            with open(csvname, 'w') as csv:
                # print(data)
                # print(list(zip(*data)))
                # print(data2)
                ddata = dict(zip(*data))
                ddata2 = dict(zip(*data2))
                csv.write('Label;%s;%s\n' % (par[0], par[0]))
                csv.write(
                    '\n'.join(
                        map(
                            lambda k:
                                ';'.join(
                                    map(
                                        str,
                                        [k,
                                         ddata[k] if k in ddata else '',
                                         ddata2[k] if k in ddata2 else '']
                                    )
                                ),
                            ordr
                        )
                    )
                )

            self.multi_barplots[par[2]] = \
                plot.MultiBarplot(
                    data[0], data[1],
                    categories=self.cats,
                    color=self.labcol,
                    cat_ordr=self.cat_ordr,
                    ylab=par[0],
                    title=par[1],
                    desc=False,
                    fname=self.get_path('%s-by-db_%s.pdf' %
                                        (par[2], self.name)),
                    order=ordr,
                    y2=None if par[4] is None else data2[1],
                    color2=None if par[4] is None else self.labcol2,
                    summary=self.multi_barplots_summary,
                    do=True
            )

    def barplot_colors(self):

        self.data_protein_counts = \
            list(zip(*[(s,
                        len([v for v in self.pp.graph.vs
                             if s in v['sources']]))
                       for s in sorted(self.pp.sources)] +
                     list(map(lambda c: (('All', c),
                                         self.csep[c].vcount()),
                              sorted(self.csep.keys())))))

        self.labcol = \
            list(
                map(
                    lambda lab:
                        self.ccolors[data_formats.categories[lab]]
                    if lab in data_formats.categories
                    else self.ccolors[lab[1]],
                    self.data_protein_counts[0]
                )
            )

        self.labcol2 = \
            list(
                map(
                    lambda lab:
                        self.ccolors2[data_formats.categories[lab]]
                    if lab in data_formats.categories
                    else self.ccolors2[lab[1]],
                    self.data_protein_counts[0]
                )
            )

    def make_coverage_groups_plot(self):

        self.coverage_groups = \
            plot.MultiBarplot(
                self.labels_coverage_groups,
                self.data_coverage_groups,
                categories=self.cats,
                color=[self.group_colors[0],
                       self.group_colors[2],
                       self.group_colors[1],
                       self.group_colors[3]],
                group_labels=self.coverage_groups_group_labels,
                cat_ordr=self.cat_ordr,
                ylab=r'Coverage [%]',
                title='Coverage of resources on different groups of proteins',
                desc=False,
                grouped=True,
                fname=self.get_path(
                    'interesting-proteins-cov-by-db_%s.pdf' % self.name),
                order=self.vcount_ordr,
                ylim=[0.0, 100.0]
            )

    def get_coverage_groups_data(self):

        covdata = list(
            map(
                lambda fun:
                    dict(zip(*self.get_data(fun, None))),
                [
                    lambda gs: len([v for v in gs[0].vs if v['rec']]) /
                        float(len(self.pp.lists['rec'])) * 100.0,
                    lambda gs: len([v for v in gs[0].vs if v['tf']]) /
                        float(len(self.pp.lists['tf'])) * 100.0,
                    lambda gs: len([v for v in gs[0].vs if v['kin']]) /
                        float(len(self.pp.lists['kin'])) * 100.0,
                    lambda gs: len([v for v in gs[0].vs if v['dgb']]) /
                        float(len(self.pp.lists['dgb'])) * 100.0
                ]
            )
        )

        self.labels_coverage_groups = list(covdata[0].keys())

        self.data_coverage_groups = \
            list(
                map(
                    lambda d:
                        list(
                            map(
                                lambda k:
                                d[k],
                                self.labels_coverage_groups
                            )
                        ),
                    covdata
                )
            )

        self.coverage_groups_group_labels = [
            'Receptors (all: %s)' % locale.format(
                '%d', len(self.pp.lists['rec']), grouping=True),
            'TFs (all: %s)' % locale.format(
                '%d', len(self.pp.lists['tf']), grouping=True),
            'Kinases (all: %s)' % locale.format(
                '%d', len(self.pp.lists['kin']), grouping=True),
            'Druggable proteins (all: %s)' % locale.format(
                '%d', len(self.pp.lists['dgb']), grouping=True)
        ]

    def get_multibarplot_ordr(self):

        self.vcount_ordr_barplot = \
            plot.MultiBarplot(
                self.data_protein_counts[0],
                self.data_protein_counts[1],
                categories=self.cats,
                color=self.labcol,
                cat_ordr=self.cat_ordr,
                ylab='Number of proteins',
                title='Number of proteins',
                desc=True,
                lab_angle=90,
                fname=self.pdf_vcount_order,
                order='y'
            )

        os.remove(self.pdf_vcount_order)

        self.vcount_ordr = self.vcount_ordr_barplot.x

    def make_scatterplots(self):

        for par in self.scatterplots_settings:

            _ = sys.stdout.write('\t:: Plotting %s\n' % par[3]['title'])

            xattr = 'data_%s_x' % par[4]
            yattr = 'data_%s_y' % par[4]
            sattr = 'data_%s_s' % par[4]
            lattr = 'labels_scatterplot_%s' % par[4]

            self.get_data(par[0], xattr)
            self.get_data(par[1], yattr)
            self.get_data(par[2], sattr)

            x = getattr(self, xattr)
            y = getattr(self, yattr)
            s = getattr(self, sattr)

            setattr(self, lattr,
                    list(
                        filter(lambda l: type(l) is not tuple,
                               getattr(self, xattr)[0])))

            labels = getattr(self, lattr)

            x = dict(zip(*x))
            y = dict(zip(*y))
            s = dict(zip(*s))

            x = list(map(lambda l: x[l], labels))
            y = list(map(lambda l: y[l], labels))
            s = np.array(list(map(lambda l: s[l], labels)))

            if par[5]:
                s = s / 1.0

            colors = list(
                map(lambda l: self.ccolors2[data_formats.categories[l]],
                    labels))

            color_labels = []
            for c in ['p', 'm', 'i', 'r']:
                color_labels.append(
                    (data_formats.catnames[c], self.ccolors2[c]))

            for i, l in enumerate(labels):
                if type(l) is tuple:
                    labels[i] = data_formats.catnames[l[1]]

            csvname = self.get_path('%s_%s.csv' % (par[4], self.name))

            with open(csvname, 'w') as csv:

                csv.write('Label;%s;%s;%s\n' %
                          (par[3]['xlab'], par[3]['ylab'], par[3]['legtitle']))
                csv.write('\n'.join(
                    map(lambda l: ';'.join(map(str, l)), zip(labels, x, y,
                                                             s))))

            sp = plot.ScatterPlus(
                x,
                y,
                size=s,
                min_size=30,
                max_size=1000,
                labels=labels,
                color=colors,
                color_labels=color_labels,
                fname=self.get_path('%s_%s.pdf' % (par[4], self.name)),
                **par[3])

    def all_ptms_list(self):

        self.ptms = {
            #'DEPOD': net.load_depod_dmi(return_raw = True),
            'Signor': self.pp.load_signor_ptms(return_raw=True),
            'Li2012': self.pp.load_li2012_ptms(return_raw=True),
            'HPRD': self.pp.load_hprd_ptms(return_raw=True),
            'MIMP': self.pp.load_mimp_dmi(return_raw=True),
            'PhosphoNetworks': self.pp.load_pnetworks_dmi(return_raw=True),
            'PhosphoELM': self.pp.load_phosphoelm(return_raw=True),
            'dbPTM': self.pp.load_dbptm(return_raw=True),
            'PhosphoSite': self.pp.load_psite_phos(return_raw=True)
        }

    def make_ptms_barplot(self):

        self.pp.uniq_ptms()

        self.ptms_all_in_network = set(
            uniqList(flatList(self.pp.graph.es['ptm'])))

        self.data_ptms = list(
            zip(*[(s, len(uniqList(m)), len(
                set(m) & self.ptms_all_in_network))
                  for s, m in self.ptms.items()] + [('All', len(
                      uniqList(flatList(self.ptms.values()))), len(
                          self.ptms_all_in_network))]))

        self.ptms_barplot = \
            plot.MultiBarplot(
                self.data_ptms[0],
                self.data_ptms[1:],
                categories=[0] * len(self.data_ptms[0]),
                color=[self.ccolors['p'], self.ccolors['m']],
                group_labels=['Total', 'In network'],
                cat_names=['E-S resources'],
                ylab=r'E-S interactions',
                title='Enzyme-substrate interactions\nmapped '
                'to combined network of %s resources' % self.title,
                desc=True,
                grouped=True,
                fname=self.get_path(
                    'ptms-by-ptmdb_%s.pdf' % self.name),
                order='y',
                figsize=(6, 4),
                title_font={'size': 'large'},
                legend_font={'size': 'x-large'},
                bar_args={'width': 0.8}
            )

    def make_htp_characteristics(self):

        self.htp_char = \
            plot.HtpCharacteristics(
                self.pp,
                fname=self.get_path('htp-%u_%s.pdf' %
                                    (self.htp_upper, self.name)),
                title='%s resources' % self.title.capitalize(),
                lower=self.htp_lower,
                upper=self.htp_upper,
                htdata=self.htdata
            )
        self.htdata = self.htp_char.htdata

    def make_history_tree(self):

        self.history_tree = plot.HistoryTree(
            fname=self.get_path(self.history_tree_fname),
            compile=False,
            dotlineopacity=1.0)

    def compile_history_tree(self):

        self.console('Running `%s` on `%s`' %
                     (self.latex, self.get_path(self.history_tree_fname)))

        self.history_tree_latex_proc = subprocess.Popen(
            [self.latex, self.get_path(self.history_tree_fname)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        self.history_tree_latex_output, self.history_tree_latex_error = \
            self.history_tree_latex_proc.communicate()
        self.history_tree_latex_return = \
            self.history_tree_latex_proc.returncode

        self.console('LaTeX %s' % ('compiled successfully'
                                   if self.history_tree_latex_return == 0 else
                                   'compilation failed'))

    def make_refs_years_grid(self):

        self.refs_years_grid = \
            plot.BarplotsGrid(
                self.pp,
                'year',
                xmin=1970,
                by='database',
                fname=self.get_path(self.refs_year_grid_fname),
                ylab='References',
                title='Number of references by databases and years',
                data=self.pubmeds,
                uniform_ylim=True
            )

    def make_refs_journals_grid(self):

        self.refs_journals_grid = \
            plot.BarplotsGrid(
                self.pp,
                'journal',
                by='database',
                fname=self.get_path(self.refs_journal_grid_fname),
                ylab='References',
                title='Number of references by databases and journals',
                full_range_x=False,
                xlim=[-0.5, 10.5],
                sort=True,
                desc=True,
                hoffset=10,
                data=self.pubmeds,
                small_xticklabels=True
            )

    def make_refs_by_db(self):
        self.console('Plotting references by database')
        bydb_ordr = reversed(self.pubmeds.database.value_counts().sort_values()
                             .index)
        bydb_ordr = ['All'] + list(bydb_ordr)

        self.refs_by_db = \
            plot.MultiBarplot(
                bydb_ordr,
                [self.pubmeds.pmid.unique().shape[0]] +
                list(self.pubmeds.database.value_counts()),
                fname=self.get_path(self.refs_by_db_file % self.name),
                order=bydb_ordr,
                xlab='Resources',
                ylab='Number of PubMed IDs',
                title='Number of references by resources in %s databases' % self.name.capitalize(),
                color=['#44AA99'] + ['#88CCEE'] *
                len(self.pubmeds.database.unique()),
                figsize=(12, 9),
                desc=False
            )

    def make_refs_by_year(self):
        self.console('Plotting references by year')
        counts = dict(self.pubmeds.year.value_counts())
        ecounts = dict(self.pubmeds_earliest.year.value_counts())
        years = np.arange(1970, max(self.pubmeds.year) + 1)
        #years = years[list(years).index(1970):]
        values = np.array(
            list(map(lambda y: counts[y] if y in counts else 0.0, years)))
        evalues = np.array(
            list(map(lambda y: ecounts[y] if y in counts else 0.0, years)))
        years = list(map(lambda y: '%u' % y if y % 5 == 0 else '', years))

        self.refs_by_year = \
            plot.MultiBarplot(
                years,
                [values, evalues],
                fname=self.get_path(self.refs_by_year_file % self.name),
                order=years,
                grouped=True,
                xlab='Years',
                ylab='Number of PubMed IDs',
                title='Number of references by year in %s databases' % self.name,
                legend_font={'size': 'x-large'},
                color=['#77AADD', '#88CCAA'],
                figsize=(8, 5),
                group_labels=['All references', 'Earliest references'],
                desc=False,
                legloc=2
            )

    def make_refs_by_journal(self):

        self.console('Plotting references by journal')

        for i in [50, 100]:
            byj_ordr = list(
                reversed(self.pubmeds.journal.value_counts().sort_values()
                         .index))[:i]
            byj_vals = list(
                reversed(self.pubmeds.journal.value_counts().sort_values()))[:
                                                                             i]

            size = 7 if i == 100 else 'small'

            self.refs_by_journal = \
                plot.MultiBarplot(
                    byj_ordr,
                    byj_vals,
                    ylog=True,
                    fname=self.get_path(self.refs_by_j_file % (i, self.name)),
                    order=byj_ordr,
                    xlab='Journals',
                    ylab='Number of PubMed IDs',
                    title='Number of references by journals in %s databases' %
                    self.name.capitalize(),
                    color='#4477AA',
                    figsize=(9, 6.75),
                    xticklabel_font={'size': size},
                    desc=False,
                    maketitle=False
                )

    def make_simgraph_vertex(self):
        self.simgraph_vertex = \
            plot.SimilarityGraph(
                self.pp,
                fname=self.get_path(self.simgraph_vertex_fname),
                similarity='vertex',
                size='vertex',
            )

    def make_simgraph_curation(self):
        self.simgraph_curation = \
            plot.SimilarityGraph(
                self.pp,
                fname=self.get_path(self.simgraph_curation_fname),
                similarity='curation',
                size='curation',
            )

    def make_simgraph_edge(self):
        self.simgraph_edge = \
            plot.SimilarityGraph(
                self.pp,
                fname=self.get_path(self.simgraph_edge_fname),
                similarity='edge',
                size='edge',
            )

    def make_dirs_stacked(self, include_all=True):

        names = [
            'Positive (All: {:,g})'.format(self.data_dirs[1][-1]),
            'Negative (All: {:,g})'.format(self.data_dirs[2][-1]),
            'Unknown effect (All: {:,g})'.format(self.data_dirs[3][-1]),
            'Unknown direction (All: {:,g})'.format(self.data_dirs[4][-1]),
        ]
        if include_all:
            x = self.data_dirs[0]
            y = self.data_dirs[1:]
            ordr = ['All'] + \
                list(filter(lambda l: type(l) != tuple, self.vcount_ordr))
        else:
            x = self.data_dirs[0][:-1]
            y = list(map(lambda yy: yy[:-1], self.data_dirs[1:]))
            ordr = list(
                filter(lambda l: l != 'All' and type(l) != tuple,
                       self.vcount_ordr))

        setattr(
            self,
            'dirs_stacked%s' % ('_all' if include_all else ''),
            plot.StackedBarplot(
                x=x,
                y=y,
                names=names,
                colors=list(reversed(self.group_colors2)),
                ylab='Interactions',
                xlab='Resources',
                title='Interactions with direction or sign',
                order=ordr,
                fname=self.get_path(self.dirs_stacked_fname %
                                    ('all'
                                     if include_all else 'wo-all', self.name)),
                desc=False))

    def make_curation_plot(self):

        self.curation_plot = plot.CurationPlot(
            self.pp,
            self.get_path(self.curation_plot_fname),
            colors=self.group_colors,
            pubmeds=self.pubmeds)

    def make_refs_composite(self):
        self.refs_composite = plot.RefsComposite(
            self.pp,
            self.get_path(self.refs_composite_fname % self.name),
            pubmeds=self.pubmeds,
            earliest=self.pubmeds_earliest)

    def get_dirs_data(self):

        self.data_dirs = \
            list(zip(*[
                (
                    s,
                    sum([sum([s in e['dirs'].positive_sources[e['dirs'].straight],
                              s in e['dirs'].positive_sources[e['dirs'].reverse]]) for e in g.es]),
                    sum([sum([s in e['dirs'].negative_sources[e['dirs'].straight],
                              s in e['dirs'].negative_sources[e['dirs'].reverse]]) for e in g.es]),
                    sum([sum([s in ((e['dirs'].sources[e['dirs'].straight] -
                                     e['dirs'].positive_sources[e['dirs'].straight]) -
                                    e['dirs'].negative_sources[e['dirs'].straight]),
                              s in ((e['dirs'].sources[e['dirs'].reverse] -
                                     e['dirs'].positive_sources[e['dirs'].reverse]) -
                                    e['dirs'].negative_sources[e['dirs'].reverse])]) for e in g.es]),
                    sum([s in e['dirs'].sources['undirected'] for e in g.es])
                )
                for s, g in iteritems(self.sep)] +
                [(
                    'All',
                    sum([e['dirs'].is_stimulation()
                         for e in self.pp.graph.es]),
                    sum([e['dirs'].is_inhibition()
                         for e in self.pp.graph.es]),
                    sum([not e['dirs'].is_stimulation() and
                         not e['dirs'].is_inhibition() and
                         e['dirs'].is_directed() for e in self.pp.graph.es]),
                    sum([not e['dirs'].is_stimulation() and
                         not e['dirs'].is_inhibition() and
                         not e['dirs'].is_directed() for e in self.pp.graph.es])
                )]
            ))

    def curation_table(self):

        self.console('Making curation statistics '
                     'LaTeX table, writing to file `%s`' %
                     self.get_path(self.table2file))
        self.pp.curation_tab(
            latex_hdr=True, fname=self.get_path(self.table2file))

        self.console('Making curation statistics '
                     'LaTeX table without header, writing to file `%s`' %
                     self.get_path(self.stable2file))
        self.pp.curation_tab(
            latex_hdr=False, fname=self.get_path(self.stable2file))

    def main_table(self):

        self.pp.basic_stats(
            latex=True, fname=self.get_path(self.main_table_fname))
        self.pp.basic_stats(
            latex=True,
            fname=self.get_path(self.main_table_stripped_fname),
            latex_hdr=False)

    def resource_list_table(self):

        descriptions.resource_list_latex(
            filename=self.get_path(self.resource_list_fname))
        descriptions.resource_list_latex(
            filename=self.get_path(self.resource_list_fname_stripped),
            latex_hdr=False)

    def latex_compile(self, fname):

        self.console('Running `%s` on `%s`' %
                     (self.latex, self.get_path(fname)))

        try:
            self.latex_proc = subprocess.Popen(
                [self.latex, self.get_path(fname)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            self.latex_output, self.latex_error = self.latex_proc.communicate(
                timeout=self.latex_timeout)
            self.latex_return = self.latex_proc.returncode
        except subprocess.TimeoutExpired:
            self.latex_return = 1

        self.console('LaTeX %s' % ('compiled successfully'
                                   if self.latex_return == 0 else
                                   'compilation failed'))

    def inconsistency_data(self):

        if self.inc_raw is None:
            self.inc_raw = self.pp.consistency()

        self.inco =  \
            dict(
                map(
                    lambda m:
                        (m,
                         dict(
                             map(
                                 lambda t:
                                 (t,
                                     dict(
                                         zip(
                                             self.pp.sources,
                                             map(
                                                 lambda n:
                                                 n if type(
                                                     n) is int else len(n),
                                                 map(
                                                     lambda s:
                                                     reduce(
                                                         lambda a, b:
                                                         a +
                                                         b if type(
                                                             a) is int else a | b,
                                                         map(
                                                             lambda ic1:
                                                             ic1[1][
                                                                 t],
                                                             filter(
                                                                 lambda ic1:
                                                                 ic1[0][
                                                                     0] == s,
                                                                 iteritems(
                                                                     self.inc_raw['inconsistency'][m])
                                                             )
                                                         )
                                                     ),
                                                     self.pp.sources
                                                 )
                                             )
                                         )
                                     )
                                  ),
                                 ['total', 'major', 'minor']
                             )
                         )
                         ),
                    ['directions', 'signs', 'directions_edges', 'signs_edges']
                )
            )

        self.cons = \
            dict(
                map(
                    lambda m:
                        (m,
                         dict(
                             map(
                                 lambda t:
                                 (t,
                                     dict(
                                         zip(
                                             self.pp.sources,
                                             map(
                                                 lambda n:
                                                 n if type(
                                                     n) is int else len(n),
                                                 map(
                                                     lambda s:
                                                     reduce(
                                                         lambda a, b:
                                                         a +
                                                         b if type(
                                                             a) is int else a | b,
                                                         map(
                                                             lambda ic1:
                                                             ic1[1][
                                                                 t],
                                                             filter(
                                                                 lambda ic1:
                                                                 ic1[0][
                                                                     0] == s,
                                                                 iteritems(
                                                                     self.inc_raw['consistency'][m])
                                                             )
                                                         )
                                                     ),
                                                     self.pp.sources
                                                 )
                                             )
                                         )
                                     )
                                  ),
                                 ['total', 'major', 'minor']
                             )
                         )
                         ),
                    ['directions', 'signs', 'directions_edges', 'signs_edges']
                )
            )

        self.incid = dict(
            map(lambda s: (s, self.inco['directions_edges']['minor'][s] / float(self.inco['directions_edges']['minor'][s] + self.cons['directions_edges']['total'][s] + 0.0001)),
                self.pp.sources))

        self.incis = dict(
            map(lambda s: (s, self.inco['signs_edges']['minor'][s] / float(self.inco['signs_edges']['minor'][s] + self.cons['signs_edges']['total'][s] + 0.0001)),
                self.pp.sources))

        self.incidp = dict(map(lambda s: (s[0], (len(s[1]['total']) /
                                                 float(len(s[1]['total']) +
                                                       len(self.inc_raw['consistency']['directions_edges'][s[0]]['total']) + 0.0001))
                                          if (len(s[1]['total']) +
                                              len(self.inc_raw['consistency']['directions_edges'][s[0]]['total'])) > 5
                                          else np.nan),
                               iteritems(self.inc_raw['inconsistency']['directions_edges'])))

        self.incidp = dict(map(lambda s: (s[0], (len(s[1]['total']) /
                                                 float(len(s[1]['total']) +
                                                       len(self.inc_raw['consistency']['directions_edges'][s[0]]['total']) + 0.0001))),
                               iteritems(self.inc_raw['inconsistency']['directions_edges'])))

        self.incdf = pd.DataFrame(
            np.vstack([
                np.array(
                    [self.incidp[(s1, s2)] for s1 in sorted(self.pp.sources)])
                for s2 in sorted(self.pp.sources)
            ]),
            index=sorted(self.pp.sources),
            columns=sorted(self.pp.sources))

        self.incdf = self.incdf.loc[(self.incdf.sum(axis=1) != 0), (
            self.incdf.sum(axis=0) != 0)]

        for lab in [
                'HPRD', 'NetPath', 'Reactome', 'ACSN', 'PDZBase', 'NCI-PID'
        ]:
            if lab in self.incdf:
                self.incdf.__delitem__(lab)
                self.incdf = self.incdf.drop(lab)

    def make_consistency_dendrogram(self):

        self.cons_dendro = plot.Dendrogram(
            fname=self.get_path(self.consistency_dedrogram_fname),
            data=self.incdf)

    def make_consistency_table(self):

        undirected = set(
            map(lambda s: s[0],
                filter(lambda s: s[1] == 0, iteritems(self.incid))))
        signed = set(
            map(lambda s: s[0],
                filter(lambda s: s[1] > 0, iteritems(self.incis))))
        directed = set(
            map(lambda s: s[0],
                filter(lambda s: s[1] > 0, iteritems(self.incid)))) - signed

        singles_undirected = dict(map(lambda s:
                                      (s, len(list(filter(lambda e:
                                                          s in e['sources'] and len(
                                                              e['sources']) == 1,
                                                          self.pp.graph.es)))),
                                      undirected))
        singles_directed = dict(map(lambda s:
                                    (s, sum(map(
                                        lambda e: sum([
                                            s in e['dirs'].sources_straight() and len(
                                                e['dirs'].sources_straight()) == 1,
                                            s in e['dirs'].sources_reverse() and len(e['dirs'].sources_reverse()) == 1]),
                                        self.pp.graph.es))), directed))
        singles_signed = dict(map(lambda s:
                                  (s, sum(map(
                                      lambda e: sum([
                                          s in e['dirs'].positive_sources_straight() and len(
                                              e['dirs'].positive_sources_straight()) == 1,
                                          s in e['dirs'].positive_sources_reverse() and len(
                                              e['dirs'].positive_sources_reverse()) == 1,
                                          s in e['dirs'].negative_sources_straight() and len(
                                              e['dirs'].negative_sources_straight()) == 1,
                                          s in e['dirs'].negative_sources_reverse() and len(e['dirs'].negative_sources_reverse()) == 1]),
                                      self.pp.graph.es))), signed))

        scores_undirected = dict(map(lambda s: (s[0], s[1] /
                                                float(max(singles_undirected.values()))), iteritems(singles_undirected)))

        scores_directed = dict(map(lambda s: (s[0], (s[1] /
                                                     float(max(singles_directed.values())) +
                                                     (1 - self.incid[s[0]]) / max(map(lambda x: 1 - x, self.incid.values()))) / 2.0),
                                   iteritems(singles_directed)))

        scores_signed = dict(map(lambda s: (s[0], (s[1] /
                                                   float(max(singles_signed.values())) +
                                                   (1 - self.incis[s[0]]) / max(map(lambda x: 1 - x, self.incis.values())) +
                                                   (1 - self.incid[s[0]]) / max(map(lambda x: 1 - x, self.incid.values()))) / 3.0),
                                 iteritems(singles_signed)))

        tbl = r'''\begin{tabularx}{\textwidth}{p{3.5cm}>{\raggedright\arraybackslash}X>{\raggedright\arraybackslash}X>{\raggedright\arraybackslash}X>{\raggedright\arraybackslash}X}
            \toprule
            Resurce & Specific interactions & Direction inconsistency & Effect inconsistency & Combined score \\
            \midrule
            \multicolumn{5}{l}{Undirected resources} \\
            \midrule
                '''

        for s, sc in reversed(
                sorted(
                    iteritems(scores_undirected), key=lambda s: s[1])):
            tbl += r'''    %s & %u &   &  & %.04f \\
                ''' % (s, singles_undirected[s], scores_undirected[s])

        tbl += r'''    \midrule
            \multicolumn{5}{l}{Directed resources} \\
            \midrule
                '''

        for s, sc in reversed(
                sorted(
                    iteritems(scores_directed), key=lambda s: s[1])):
            tbl += r'''    %s & %u &  %.04f &  & %.04f \\
                ''' % (s, singles_directed[s], self.incid[s],
                       scores_directed[s])

        tbl += r'''    \midrule
            \multicolumn{5}{l}{Signed resources} \\
            \midrule
                '''

        for s, sc in reversed(
                sorted(
                    iteritems(scores_signed), key=lambda s: s[1])):
            tbl += r'''    %s & %u &  %.04f &  %.04f & %.04f \\
                ''' % (s, singles_signed[s], self.incid[s], self.incis[s],
                       scores_signed[s])

        tbl += r'''    \bottomrule
        \end{tabularx}'''

        with open(self.get_path(self.consistency_table_fname), 'w') as f:
            f.write(tbl)

        self.consistency_table = tbl

    def load_pubmed_data(self):
        self.pubmeds, self.pubmeds_earliest = _refs.get_pubmed_data(
            self.pp, htp_threshold=None)

    def console(self, msg):
        _ = sys.stdout.write('\t:: %s\n' % msg)
        sys.stdout.flush()
