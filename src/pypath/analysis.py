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

import pypath.main as main
import pypath.plot as plot
import pypath.data_formats as data_formats
import pypath.dataio as dataio
from pypath.common import *

class Workflow(object):
    
    def __init__(self,
                 name,
                 network_datasest = [],
                 do_multi_barplots = True,
                 do_coverage_groups = True,
                 do_htp_char = True,
                 do_ptms_barplot = True,
                 do_scatterplots = True,
                 do_history_tree = True,
                 do_compile_history_tree = True,
                 do_refs_journals_grid = True,
                 do_refs_years_grid = True,
                 do_dirs_stacked = True,
                 title = None,
                 **kwargs):
        
        for k, v in iteritems(locals()):
            setattr(self, k, v)
        
        self.title = self.name if self.title is None else self.title
        
        self.defaults = {
            'ccolors': {
                'p': '#A9C98B',
                'r': '#E87A62',
                'm': '#FEEA9C',
                'i': '#6FA6A9'
            },
            'ccolors2': {
                'p': '#84B45C',
                'r': '#DE3B34',
                'm': '#FCCC06',
                'i': '#00888B'
            },
            'group_colors': ['#FE5222', '#F5FC8A', '#0EACD3', '#CDEC25'],
            'table2file': 'curation_tab_%s.tex' % self.name,
            'stable2file': 'curation_tab_stripped_%s.tex' % self.name,
            'latex': '/usr/bin/xelatex',
            'compile_latex': True,
            'protein_lists': {
                'proteome': 'human proteome',
                'disease_genes': 'DisGeNet disease related genes',
                'signaling_proteins': 'signaling proteins'
            },
            'intogen_file': None,
            'set_annots': {
                'receptors': 'receptors',
                'tfs': 'transcription factors',
                'kinases': 'kinases',
                'druggability': 'druggable proteins'
            },
            'load_annots' : {
                'corum': 'CORUM complexes',
                'pmts': 'post-translational modifications',
                'disgenet': 'DisGeNet disease related genes'
            },
            'fiher_file': 'fisher_%s' % self.name,
            'fisher': [
                ('dis', 'Disease related genes'),
                ('rec', 'Receptors'),
                ('tf', 'Transcription factors'),
                ('kin', 'Kinases'),
                ('dgb', 'Druggable proteins'),
                ('cdv', 'Cancer drivers')
            ],
            'cat_ordr': ['Activity flow', 'Enzyme-substrate', 'Interaction', 'Process description'],
            'pdf_vcount_order': 'vcount_ordr.pdf',
            'htp_lower': 1,
            'htp_upper': 500,
            'history_tree_fname': 'history_tree.tex',
            'simgraph_vertex_fname': 'sources_similarity_vertex_%s.pdf' % self.name,
            'simgraph_edge_fname': 'sources_similarity_edge_%s.pdf' % self.name,
            'refs_journal_grid_fname': 'refs_by_db_journal_%s.pdf' % self.name,
            'refs_year_grid_fname': 'refs_by_db_year_%s.pdf' % self.name,
            'dirs_stacked_fname': 'dirs-signes-by-db-%s_%s.pdf',
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
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['rec']]),
                'vcount'
            ),
            (
                r'Receptors [%]',
                'Percentage of receptors',
                'receptorprop',
                lambda gs: len([v for v in gs[0].vs if v['rec']]) / \
                    float(gs[0].vcount()) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['rec']]) / \
                    float(gs[0].vcount()) * 100.0,
                'vcount'
            ),
            (
                r'Receptors [%]',
                'Percentage of all human receptors covered',
                'receptorcov',
                lambda gs: len([v for v in gs[0].vs if v['rec']]) / \
                    float(len(self.pp.lists['rec'])) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['rec']]) / \
                    float(len(self.pp.lists['rec'])) * 100.0,
                'vcount'
            ),
            (
                'TFs',
                'Number of transcription factors',
                'tfs',
                lambda gs: len([v for v in gs[0].vs if v['tf']]),
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['tf']]),
                'vcount'
            ),
            (
                r'TFs [%]',
                'Percentage of transcription factors',
                'tfprop',
                lambda gs: len([v for v in gs[0].vs if v['tf']]) / \
                    float(gs[0].vcount()) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['tf']]) / \
                    float(gs[0].vcount()) * 100.0,
                'vcount'
            ),
            (
                r'TFs [%]',
                'Percentage of all human TFs covered',
                'tfcov',
                lambda gs: len([v for v in gs[0].vs if v['tf']]) / \
                    float(len(self.pp.lists['tfs'])) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['tf']]) / \
                    float(len(self.pp.lists['tfs'])) * 100.0,
                'vcount'
            ),
            (
                r'Kinases [%]',
                'Percentage of all human kinases covered',
                'kincov',
                lambda gs: len([v for v in gs[0].vs if v['kin']]) / \
                    float(len(self.pp.lists['kin'])) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['kin']]) / \
                    float(len(self.pp.lists['kin'])) * 100.0,
                'vcount'
            ),
            (
                r'Druggable proteins [%]',
                'Percentage of all human druggable proteins covered',
                'dgbcov',
                lambda gs: len([v for v in gs[0].vs if v['dgb']]) / \
                    float(len(self.pp.lists['dgb'])) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['dgb']]) / \
                    float(len(self.pp.lists['dgb'])) * 100.0,
                'vcount'
            ),
            (
                r'Disease genes [%]',
                'Percentage of disease-gene associations covered',
                'discov',
                lambda gs: len([v for v in gs[0].vs if v['dis']]) / \
                    float(sum(len(x) for x in dis.values())) * 100.0,
                lambda gs: len([v for v in self.specific(gs[1], gs[0].vs) if v['dis']]) / \
                    float(sum(len(x) for x in dis.values())) * 100.0,
                'vcount'
            ),
            (
                'Complexes',
                'Number of complexes covered',
                'complexes',
                lambda gs: len(sens.complexes_in_network(gs[0])),
                None,
                'vcount'
            ),
            (
                'E-S interactions',
                'Number of enzyme-substrate interactions covered',
                'ptms',
                lambda gs: sum(map(lambda e: len(e['ptm']), gs[0].es)),
                lambda gs: sum(map(lambda e: len(e['ptm']), self.specific(gs[1], gs[0].es))),
                'vcount'
            ),
            (
                'E-S interactions',
                'Number of interactions associated with enzyme-substrate relationship',
                'havingptm',
                lambda gs: sum(map(lambda e: len(e['ptm']) > 0, gs[0].es)),
                lambda gs: sum(map(lambda e: len(e['ptm']) > 0, self.specific(gs[1], gs[0].es))),
                'vcount'
            ),
            (
                r'CGC genes [%]',
                'Percentage of COSMIC Cancer Gene Census cancer drivers covered',
                'ccgccov',
                lambda gs: len(set(gs[0].vs['name']) & set(self.pp.lists['cgc'])) / \
                    float(len(self.pp.lists['cgc'])) * 100.0,
                lambda gs: len(
                                set(map(lambda v: v['name'], self.specific(gs[1], gs[0].vs))) & \
                                set(self.pp.lists['cgc'])
                            ) / \
                            float(len(self.pp.lists['cgc'])) * 100.0,
                'vcount'
            ),
            (
                r'IntOGen genes [%]',
                'Percentage of IntOGen cancer drivers covered',
                'intocov',
                lambda gs: len(set(gs[0].vs['name']) & set(self.pp.lists['IntOGen'])) / \
                    float(len(self.pp.lists['IntOGen'])) * 100.0,
                lambda gs: len(
                                set(map(lambda v: v['name'], self.specific(gs[1], gs[0].vs))) & \
                                set(self.pp.lists['IntOGen'])
                            ) / \
                            float(len(self.pp.lists['IntOGen'])) * 100.0,
                'vcount'
            ),
            (
                r'Signaling proteins [%]',
                'Percentage of all human signaling proteins covered',
                'sigcov',
                lambda gs: len(set(gs[0].vs['name']) & upsig) / \
                    float(len(upsig)) * 100.0,
                lambda gs: len(
                                set(map(lambda v: v['name'], self.specific(gs[1], gs[0].vs))) & \
                                set(self.pp.lists['sig'])
                            ) / \
                            float(len(self.pp.lists['sig'])) * 100.0,
                'vcount'
            ),
            (
                r'Signaling proteins [%]',
                'Percentage of signaling proteins',
                'sigpct',
                lambda gs: len(set(gs[0].vs['name']) & upsig) / \
                    float(gs[0].vcount()) * 100.0,
                lambda gs: len(
                                set(map(lambda v: v['name'], self.specific(gs[1], gs[0].vs))) & \
                                set(self.pp.lists['sig'])
                            ) / \
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
                lambda gs: len(set(gs[0].vs['name']) & (set(self.pp.lists['dis']))),
                lambda gs: len(set(gs[0].vs['name']) & set(self.pp.lists['cdv'])),
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
                    'legstrip': (3,None)
                },
                'dis-cancer',
                False
            ),
            (
                lambda gs: len(set(gs[0].vs['name']) & set(self.pp.lists['rec'])),
                lambda gs: len(set(gs[0].vs['name']) & set(self.pp.lists['tfs'])),
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
                    'legstrip': (3,None)
                },
                'tf-rec',
                False
            ),
            (
                lambda gs: len(self.pp.complexes_in_network(grapg = gs[0])),
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
                    'legstrip': (3,None)
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
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def run(self):
        
        self.init_pypath()
        self.load_protein_lists()
        self.load_annotations()
        self.set_categories()
        self.curation_table()
        self.compile_curation_table()
        self.barplot_colors()
        self.get_multibarplot_ordr()
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
                self.compile_history_tree()
        if self.do_refs_years_grid:
            self.make_refs_years_grid()
        if self.do_refs_journals_grid:
            self.make_refs_journals_grid()
        if self.do_dirs_stacked:
            self.make_dirs_stacked()
            self.make_dirs_stacked(include_all = False)
    
    def init_pypath(self):
        
        self.pp = main.PyPath(9606)
        for netdata in self.network_datasest:
            self.pp.load_resources(getattr(data_formats, netdata))
    
    def load_protein_lists(self):
        
        for meth, name in iteritems(self.protein_lists):
            
            self.console('Loading list of %s' % name)
            getattr(self.pp, '%s_list' % meth)
        
        self.console('Loading list of Cancer Gene Census %scancer drivers' % \
            ('and IntOGen ' if self.intogen_file else ''))
        self.pp.cancer_drivers_list(intogen_file = self.intogen_file)
    
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
        del self.csep['t']
        self.cats = dict(map(lambda c: (c[0], data_formats.catnames[c[1]]), iteritems(data_formats.categories)))
        self.cats.update(dict(map(lambda c: (('All', c[0]), c[1]), iteritems(data_formats.catnames))))
    
    def fisher_tests(self):
        
        self.console('Doing Fisher tests, writing results to `%s`' % self.fisher_file)
        with open(self.fisher_file, 'w') as fi:
            
            for attr, name in self.fisher:
                cont = np.array([[len(self.pp.lists['proteome']),
                                self.pp.graph.vcount()],
                                [len(self.pp.lists[attr]),
                                len([1 for v in self.pp.graph.vs if len(v[attr]) > 0])]])
                fi.write('%s:' % name)
                fi.write('\t%s\t%s\n' % stats.fisher_exact(contDisg))
    
    def get_data(self, fun, attr):
        if not hasattr(self, attr):
            setattr(self, attr,
                list(
                    zip(
                        *[(s, fun((self.sep[s], s))) for s in sorted(self.pp.sources)] + \
                        list(
                            map(
                                lambda c:
                                    (('All', c[0]), fun((self.csep[c[0]], c[0]))),
                                sorted(self.csep.keys())
                            )
                        )
                    )
                )
            )
    
    def specific(self, s, seq):
        return [w for w in seq \
                if s in data_formats.categories \
                    and s in w['sources'] \
                    and len(w['sources'] & getattr(data_formats, data_formats.categories[s])) == 1 \
                or s in w['cat'] and len(w['cat']) == 1 \
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
            
            self.multi_barplots[par[2]] = \
                plot.MultiBarplot(
                    data[0], data[1],
                    categories = self.cats,
                    color = self.labcol,
                    cat_ordr = self.cat_ordr,
                    ylab = par[0],
                    title = par[1],
                    desc = False,
                    fname = '%s_%s-by-db.pdf' % (par[2], self.name),
                    order = ordr,
                    y2 = None if par[4] is None else data2[1],
                    color2 = None if par[4] is None else self.labcol2
                )
    
    def barplot_colors(self):
        
        self.data_protein_counts = \
            list(zip(*[(s, len([v for v in net.graph.vs if s in v['sources']])) for s in sorted(net.sources)] + \
            list(map(lambda c: (('All', c), csep[c].vcount()), sorted(csep.keys())))))
        
        self.labcol = \
            list(
                map(
                    lambda lab:
                        self.ccolors[data_formats.categories[lab]] \
                            if lab in data_formats.categories \
                            else self.ccolors[lab[1]],
                    self.data_protein_counts[0]
                )
            )
        
        self.labcol2 = \
            list(
                map(
                    lambda lab:
                        self.ccolors2[data_formats.categories[lab]] \
                            if lab in data_formats.categories \
                            else self.ccolors2[lab[1]],
                    self.data_protein_counts[0]
                )
            )
    
    def make_coverage_groups_plot(self):
        
        self.coverage_groups = \
            plot.MultiBarplot(
                self.labels_coverage_groups,
                self.data_coverage_groups,
                categories = self.cats,
                color = self.group_colors,
                group_labels = self.coverage_groups_group_labels,
                cat_ordr = self.cat_ordr,
                ylab = r'Coverage [%]',
                title = 'Coverage of resources on different groups of proteins',
                desc = False,
                grouped = True,
                fname = 'interesting-proteins-cov_%s-by-db.pdf' % self.name,
                order = self.vcount_ordr,
                ylim = [0.0, 100.0]
            )
    
    def get_coverage_groups_data(self):
        
        covdata = list(
            map(
                lambda fun:
                    dict(zip(*self.get_data(fun))),
                [
                    lambda gs: len([v for v in gs[0].vs if v['rec']]) / \
                                float(len(net.lists['rec'])) * 100.0,
                    lambda gs: len([v for v in gs[0].vs if v['tf']]) / \
                                float(len(net.lists['tfs'])) * 100.0,
                    lambda gs: len([v for v in gs[0].vs if v['kin']]) / \
                                float(len(net.lists['kin'])) * 100.0,
                    lambda gs: len([v for v in gs[0].vs if v['dgb']]) / \
                                float(len(net.lists['dgb'])) * 100.0
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
                                keys
                            )
                        ),
                    self.covdata
                )
            )
        
        self.grouped_parplot_group_labels = [
            'Receptors (all: %s)' % \
                locale.format('%d', len(self.pp.lists['rec']), grouping = True),
            'TFs (all: %s)' % \
                locale.format('%d', len(self.pp.lists['tfs']), grouping = True),
            'Kinases (all: %s)' % \
                locale.format('%d', len(self.pp.lists['kin']), grouping = True),
            'Druggable proteins (all: %s)' % \
                locale.format('%d', len(self.pp.lists['dgb']), grouping = True)
        ]
    
    def get_multibarplot_ordr(self):
        self.vcount_ordr_barplot = \
            plot.MultiBarplot(
                self.data_protein_counts[0],
                self.data_protein_counts[1],
                categories = self.cats,
                color = self.labcol,
                cat_ordr = self.cat_ordr,
                ylab = 'Number of proteins',
                title = 'Number of proteins',
                desc = True,
                lab_angle = 90,
                fname = self.pdf_vcount_order,
                order = 'y'
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
                list(filter(lambda l:
                    type(l) is not tuple, getattr(self, xattr)[0])))
            
            labels = getattr(self, lattr)
            
            x = dict(zip(*x))
            y = dict(zip(*y))
            s = dict(zip(*s))
            
            x = list(map(lambda l: x[l], labels))
            y = list(map(lambda l: y[l], labels))
            s = np.array(list(map(lambda l: s[l], labels)))
            
            if par[5]:
                s = s / 1.0
            
            colors = list(map(lambda l: self.ccolors2[data_formats.categories[l]], labels))
            
            color_labels = []
            for c in ['p', 'm', 'i', 'r']:
                color_labels.append((data_formats.catnames[c], self.ccolors[c]))
            
            for i, l in enumerate(labels):
                if type(l) is tuple:
                    labels[i] = data_formats.catnames[l[1]]
            
            sp = plot.ScatterPlus(
                    x,
                    y,
                    size = s,
                    min_size = 30,
                    max_size = 1000,
                    labels = labels,
                    color = colors,
                    color_labels = color_labels,
                    fname = '%s_%s.pdf' % (par[4], self.name),
                    **par[3]
                )
    
    def all_ptms_list(self):
        
        self.ptms = {
            #'DEPOD': net.load_depod_dmi(return_raw = True),
            'Signor': net.load_signor_ptms(return_raw = True),
            'Li2012': net.load_li2012_ptms(return_raw = True),
            'HPRD': net.load_hprd_ptms(return_raw = True),
            'MIMP': net.load_mimp_dmi(return_raw = True),
            'PhosphoNetworks': net.load_pnetworks_dmi(return_raw = True),
            'PhosphoELM': net.load_phosphoelm(return_raw = True),
            'dbPTM': net.load_dbptm(return_raw = True),
            'PhosphoSite': net.load_psite_phos(return_raw = True)
        }
    
    def make_ptms_barplot(self):
        
        self.pp.uniq_ptms()
        
        self.ptms_all_in_network = set(uniqList(flatList(self.pp.graph.es['ptm'])))
        
        self.data_ptms = list(zip(*[(s,
                    len(uniqList(m)),
                    len(set(m) & self.ptms_all_in_network)
                ) for s, m in self.ptms.items()] + \
                [('All',
                    len(uniqList(flatList(self.ptms.values()))),
                    len(self.ptms_all_in_network))]))
        
        group_colors = ['#FE5222', '#F5FC8A', '#0EACD3', '#CDEC25']
        
        self.ptms_barplot = \
            plot.MultiBarplot(
                self.data_ptms[0],
                self.data_ptms[1:],
                categories = [0] * len(self.data_ptms[0]),
                color = self.group_colors,
                group_labels = ['Total', 'In network'],
                cat_names = ['E-S resources'],
                ylab = r'E-S interactions',
                title = 'Enzyme-substrate interactions\nmapped '\
                    'to combined network of %s resources' % self.title,
                desc = True,
                grouped = True,
                fname = 'ptms-by-ptmdb_%s.pdf' % self.name,
                order = 'y',
                figsize = (12, 9),
                legend_font = {'size': 'large'},
                bar_args = {'width': 0.8}
            )
    
    def make_htp_characteristics(self):
        
        self.htp_char = \
            plot.HtpCharacteristics(
                self.pp,
                fname = 'htp_%s.pdf' % self.name,
                title = '%s resources' % self.title.capitalize(),
                lower = self.htp_lower,
                upper = self.htp_upper
            )
    
    def make_history_tree(self):
        
        self.history_tree = plot.HistoryTree(fname = self.history_tree_fname)
    
    def compile_history_tree(self):
        
        self.console('Running `%s` on `%s`' % (self.latex, self.table2file))
        
        self.curation_tab_latex_proc = subprocess.Popen(
            [self.latex, self.history_tree_fname],
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE)
        self.history_tree_latex_output, self.history_tree_latex_error = \
            self.history_tree_latex_proc.communicate()
        self.history_tree_latex_return = \
            self.history_tree_latex_proc.returncode
        
        self.console('LaTeX ' % (
                'compiled successfully' \
                    if self.history_tree_latex_return == 0 \
                else 'compilation failed'))
    
    def make_refs_years_grid(self):
        
        self.refs_years_grid = \
            plot.BarplotsGrid(
                self.pp,
                'year',
                by = 'database',
                fname = self.refs_year_grid_fname,
                ylab = 'References',
                title = 'Number of references by databases and years'
            )
    
    def make_refs_journals_grid(self):
        
        self.refs_years_grid = \
            plot.BarplotsGrid(
                self.pp,
                'journal',
                by = 'database',
                fname = self.refs_journal_grid_fname,
                ylab = 'References',
                title = 'Number of references by databases and journals'
            )
    
    def make_simgraph_vertex(self):
        self.simgraph_vertex = \
            plot.SimilarityGraph(
                self.pp,
                fname = self.simgraph_vertex_fname,
                similarity = 'vertex',
                size = 'vertex',
            )
    
    def make_simgraph_edge(self):
        self.simgraph_edge = \
            plot.SimilarityGraph(
                self.pp,
                fname = self.simgraph_edge_fname,
                similarity = 'edge',
                size = 'edge',
            )
    
    def make_dirs_stacked(self, include_all = True):
        
        if include_all:
            x = self.data_dirs[0]
            y = self.data_dirs[1:]
            ordr = self.vcount_ordr
        else:
            x = self.data_dirs[0][:-1]
            y = list(map(lambda yy: yy[:-1],  self.data_dirs[1:]))
            ordr = list(filter(lambda l: l != 'All', self.vcount_ordr))
        
        setattr(self, 'dirs_stacked%s' % ('_all' if include_all else ''),
            plot.StackedBarplot(
                x = x,
                y = y,
                names = ['Positive', 'Negative', 'Unknown effect', 'Unknown direction'],
                colors = self.group_colors,
                ylab = 'Interactions',
                xlab = 'Resources',
                title = 'Interactions with direction or sign',
                order = ordr,
                fname = self.dirs_stacked_fname % \
                    ('all' if include_all else 'wo-all', self.name),
                desc = False
            ))
    
    def get_dirs_data(self):
        
        self.data_dirs = \
            list(zip(*[
                (
                    s,
                    sum([sum([s in e['dirs'].positive_sources[e['dirs'].straight],
                        s in e['dirs'].positive_sources[e['dirs'].reverse]]) for e in g.es]),
                    sum([sum([s in e['dirs'].negative_sources[e['dirs'].straight],
                        s in e['dirs'].negative_sources[e['dirs'].reverse]]) for e in g.es]),
                    sum([sum([s in ((e['dirs'].sources[e['dirs'].straight] - \
                        e['dirs'].positive_sources[e['dirs'].straight]) - \
                        e['dirs'].negative_sources[e['dirs'].straight]),
                        s in ((e['dirs'].sources[e['dirs'].reverse] - \
                        e['dirs'].positive_sources[e['dirs'].reverse]) - \
                        e['dirs'].negative_sources[e['dirs'].reverse])]) for e in g.es]),
                    sum([s in e['dirs'].sources['undirected'] for e in g.es])
                ) \
                for s, g in iteritems(self.sep)] + \
                    [(
                        'All',
                        sum([e['dirs'].is_stimulation() for e in self.pp.graph.es]),
                        sum([e['dirs'].is_inhibition() for e in self.pp.graph.es]),
                        sum([not e['dirs'].is_stimulation() and \
                            not e['dirs'].is_inhibition() and \
                            e['dirs'].is_directed() for e in self.pp.graph.es]),
                        sum([not e['dirs'].is_stimulation() and \
                            not e['dirs'].is_inhibition() and \
                            not e['dirs'].is_directed() for e in self.pp.graph.es])
                    )]
            ))
    
    def curation_table(self):
        
        self.console('Making curation statistics '\
            'LaTeX table, writing to file `%s`' % self.table2file)
        self.pp.curation_tab(latex_hdr = True,
                             fname = self.table2file)
        
        self.console('Making curation statistics '\
            'LaTeX table without header, writing to file `%s`' % self.stable2file)
        self.pp.curation_tab(latex_hdr = True,
                             fname = self.table2file)
    
    def compile_curation_table(self):
        
        if self.compile_latex:
            self.console('Running `%s` on `%s`' % (self.latex, self.table2file))
            
            self.curation_tab_latex_proc = subprocess.Popen(
                [self.latex, self.table2file],
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE)
            self.curation_tab_latex_output, self.curation_tab_latex_error = \
                self.curation_tab_latex_proc.communicate()
            self.curation_tab_latex_return = \
                self.curation_tab_latex_proc.returncode
            
            self.console('LaTeX ' % (
                'compiled successfully' \
                    if self.curation_tab_latex_return == 0 \
                else 'compilation failed'))
    
    def console(self, msg):
        _ = sys.stdout.write('\t:: %s\n' % msg)
        sys.stdout.flush()
