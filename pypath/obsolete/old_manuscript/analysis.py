#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
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
#  This code has been used to generate figures for omnipath2/pypath manuscript

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import os
import sys
import subprocess
import imp
import locale
import numpy as np
import pandas as pd
import collections

import pypath.main as main
import pypath.visual.plot as plot
import pypath.resources.data_formats as data_formats
import pypath.inputs.main as dataio
import pypath.resources.descriptions as descriptions
from pypath.share.common import *
import pypath.internals.refs as _refs
import pypath.omnipath.legacy as omnipath
import pypath.share.session as session_mod
import pypath.utils.go as go
import pypath.db_categories as db_categories


# defines a multi-section barplot
MultiBarplotParam = collections.namedtuple(
    'MultiBarplotParam',
    ['ylab', 'title', 'name', 'method', 'smethod', 'order']
)
MultiBarplotParam.__new__.__defaults__ = (None, 'vcount')


# defines a scatterplot
ScatterplotParam = collections.namedtuple(
    'ScatterplotParam',
    [
        'xmethod', # method for variable x
        'ymethod', # method for variable y
        'smethod', # method for variable size
        'gparam',  # graphics params
        'name',    # name included in the file name
        'sdiv',    # have no idea what it is...
    ]
)
ScatterplotParam.__new__.__defaults__ = (False,)


# defines graphics details for a scatterplot
ScatterplotGraphicsParam = collections.namedtuple(
    'ScatterplotGraphicsParam',
    [
        'ylim',
        'xlim',
        'xlog',
        'ylog',
        'xlab',
        'ylab',
        'legtitle',
        'title',
        'legstrip',
    ]
)
ScatterplotGraphicsParam.__new__.__defaults__ = ((None, None),)


class Workflow(omnipath.OmniPath):
    
    def __init__(
            self,
            name,
            network_datasets = [],
            do_main_table = True,
            do_compile_main_table = True,
            do_curation_table = True,
            do_compile_curation_table = True,
            do_simgraphs = True,
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
            do_refs_composite = True,
            do_curation_plot = True,
            do_refs_by_j = True,
            do_refs_by_db = True,
            do_refs_by_year = True,
            do_resource_list = True,
            do_compile_resource_list = True,
            do_consistency_dedrogram = True,
            do_consistency_table = True,
            only_categories = None,
            title = None,
            outdir = None,
            htdata = {},
            inc_raw = None,
            intogen_file = None,
            cosmic_credentials = None,
            network_pickle = None,
            annotation_pickle = None,
            intercell_pickle = None,
            complex_pickle = None,
            enz_sub_pickle = None,
            omnipath_pickle = None,
            load_network = False,
            load_complexes = True,
            load_annotations = True,
            load_intercell = True,
            load_enz_sub = True,
            **kwargs
        ):
        """
        Executes the workflow of comparative analysis of network resources
        by categories. Creates tables and figures.
        
        :arg str name:
            A label included in the name of all files.
        :arg list network_datasets:
            Datasets to use at building the network. Dicts from
            ``pypath.data_formats``.
        :arg bool do_main_table:
            Create a LaTeX table with key numbers.
        :arg bool do_compile_main_table:
            Compile the table by running xelatex.
        :arg bool do_curation_table:
            Create a table of curation effort.
        :arg bool do_compile_curation_table:
            Compile the curation effort table using xelatex.
        :arg bool do_simgraphs:
            Create graph figures of similarities between resources.
        :arg bool do_multi_barplots:
            Create multi section barplots showing size and coverage of the
            resources, categories and their overlaps.
        :arg bool do_coverage_groups:
            ??
        :arg bool do_htp_char:
            High throughput characteristics figure.
        :arg bool do_ptms_barplot:
            Barplot about size of PTM resources.
        :arg bool do_scatterplots:
            Create scatterplots of resource sizes.
        :arg bool do_history_tree:
            Create a TikZ figures about the history of the resources.
        :arg bool do_compile_history_tree:
            Compile the history figure using xelatex.
        :arg bool do_refs_journals_grid:
            Create multifacet plot showing the most often curated journals
            for each resource.
        :arg bool do_refs_years_grid:
            Create multifacet plot showing the publication years of references
            in each resource.
        :arg bool do_dirs_stacked:
            Create a stacked barplot about number of undirected, directed
            and signed interactions.
        :arg bool do_refs_composite:
            Create composite figure about references.
        :arg bool do_curation_plot:
            ??
        :arg bool do_refs_by_j:
            Create histogram of publication years.
        :arg bool do_refs_by_db:
            Create barplot of references per database.
        :arg bool do_refs_by_year:
            Create a barplot with refrences by year.
        :arg bool do_resource_list:
            Create a table with the resources listed.
        :arg bool do_compile_resource_list:
            Compile the table using xelatex.
        :arg bool do_consistency_dedrogram:
            Create a dendrogram using consistencies as a distance metric
            across all resources.
        :arg bool do_consistency_table:
            Create a table with inconsistency statistics across resources.
        :arg str outdir:
            Directory to save the output files.
        """
        
        session_mod.Logger.__init__(self, name = 'analysis')
        
        omnipath.OmniPath.__init__(
            self,
            network_pickle = network_pickle,
            annotation_pickle = annotation_pickle,
            intercell_pickle = intercell_pickle,
            complex_pickle = complex_pickle,
            enz_sub_pickle = enz_sub_pickle,
            load_network = load_network,
            load_complexes = load_complexes,
            load_annotations = load_annotations,
            load_intercell = load_intercell,
            load_enz_sub = load_enz_sub,
        )
        
        for k, v in iteritems(locals()):
            
            if not hasattr(self, k):
                
                setattr(self, k, v)

        for k, v in iteritems(kwargs):
            
            if not hasattr(self, k):
                
                setattr(self, k, v)
        
        self.title = self.title or self.name
        
        self.defaults = {
            # colors for the categories
            'ccolors': {
                'p': '#77AADD',
                'l': '#77CCCC',
                'm': '#DDAA77',
                'i': '#CC99BB',
                'r': '#77AADD',
                'o': '#AAAA44',
            },
            # colors of the shaded parts
            'ccolors2': {
                'p': '#4477AA',
                'l': '#117777',
                'm': '#774411',
                'i': '#771155',
                'r': '#4477AA',
                'o': '#AAAA44',
            },
            'group_colors': [
                '#4477AA',
                '#44AAAA',
                '#DDAA77',
                '#CC99BB',
                '#4477AA',
            ],
            'group_colors2': [
                '#77AADD',
                '#77CCCC',
                '#DDAA77',
                '#CC99BB',
                '#77AADD',
            ],
            'table2file':
                'curation_tab_%s.tex' % self.name,
            'stable2file':
                'curation_tab_stripped_%s.tex' % self.name,
            'latex':
                '/usr/bin/xelatex',
            'latex_timeout': 30,
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
            'fiher_file':
                'fisher_%s' % self.name,
            'fisher': [
                ('dis', 'Disease related genes'),
                ('rec', 'Receptors'),
                ('tf', 'Transcription factors'),
                ('kin', 'Kinases'),
                ('dgb', 'Druggable proteins'),
                ('cdv', 'Cancer drivers'),
                ('sig', 'Signaling proteins'),
            ],
            'cat_ordr_default': [
                'Activity flow',
                'Ligand-receptor',
                'Enzyme-substrate',
                'Undirected PPI',
                'Process description',
            ],
            'pdf_vcount_order':
                'vcount_ordr.pdf',
            'htp_lower': 1,
            'htp_upper': 500,
            'history_tree_fname':
                'history_tree.tex',
            'main_table_fname':
                'main_table_%s.tex' % self.name,
            'main_table_stripped_fname':
                'main_table_stripped_%s.tex' % self.name,
            'simgraph_vertex_fname':
                'sources_similarity_vertex_%s.pdf' % self.name,
            'simgraph_edge_fname':
                'sources_similarity_edge_%s.pdf' % self.name,
            'simgraph_curation_fname':
                'sources_similarity_curation_%s.pdf' % self.name,
            'refs_journal_grid_fname':
                'refs_by_db_journal_%s.pdf' % self.name,
            'refs_year_grid_fname':
                'refs_by_db_year_%s.pdf' % self.name,
            'dirs_stacked_fname':
                'dirs-signes-by-db-%s_%s.pdf',
            'refs_composite_fname':
                'refs-composite_%s.pdf',
            'refs_by_j_file':
                'references-by-journal-%u_%s.pdf',
            'refs_by_db_file':
                'references-by-db_%s.pdf',
            'refs_by_year_file':
                'references-by-year_%s.pdf',
            'curation_plot_fname':
                'new-curated-by-year_%s.pdf' % self.name,
            'resource_list_fname':
                'resources.tex',
            'resource_list_fname_stripped':
                'resources_stripped.tex',
            'consistency_dedrogram_fname':
                'inconsistency-dendrogram_%s.pdf' % self.name,
            'consistency_table_fname':
                'inconsistency-table_%s.tex' % self.name
        }
        
        # settings for barplots
        # each element of the list is for one barplot
        # settings in each tuple:
        #   - Y axis label
        #   - Main title
        #   - Label included in the file name
        #   - Method to calculate the height of each bar
        #     Accepts one argument: a tuple of resource name and
        #     an ``igraph.Graph`` object
        #   - Method to calculate the height of the shaded part of each bar
        #   - 
        self.barplots_settings = [
            MultiBarplotParam(
                'Proteins',
                'Number of proteins',
                'proteins',
                lambda gs: gs[0].vcount(),
                lambda gs: len(self.specific(gs[1], gs[0].vs)),
                'vcount'
            ),
            MultiBarplotParam(
                'Interactions',
                'Number of interactions',
                'interactions',
                lambda gs: gs[0].ecount(),
                lambda gs: len(self.specific(gs[1], gs[0].es)),
                'vcount'
            ),
            MultiBarplotParam(
                'Density',
                'Graph density',
                'density',
                lambda gs: gs[0].density(),
                None,
                'vcount'
            ),
            MultiBarplotParam(
                'Transitivity',
                'Graph global transitivity',
                'transitivity',
                lambda gs: gs[0].transitivity_undirected(),
                None,
                'vcount'
            ),
            MultiBarplotParam(
                'Diameter',
                'Graph diameter',
                'diameter',
                lambda gs: gs[0].diameter(),
                None,
                'vcount'
            ),
            MultiBarplotParam(
                'Receptors',
                'Number of receptors',
                'receptors',
                lambda gs: len([v for v in gs[0].vs if v['rec']]),
                lambda gs: len(
                    [v for v in self.specific(gs[1], gs[0].vs) if v['rec']]
                ),
                'vcount'
            ),
            MultiBarplotParam(
                r'Receptors [%]',
                'Percentage of receptors',
                'receptorprop',
                lambda gs:
                    len([v for v in gs[0].vs if v['rec']]) /
                    float(gs[0].vcount()) * 100.0,
                lambda gs:
                    len([
                        v for v in self.specific(gs[1], gs[0].vs)
                        if v['rec']
                    ]) /
                    float(gs[0].vcount()) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                r'Receptors [%]',
                'Percentage of all human receptors covered',
                'receptorcov',
                lambda gs: len([v for v in gs[0].vs if v['rec']]) /
                    float(len(self.pp.lists['rec'])) * 100.0,
                lambda gs: len([
                        v for v in self.specific(gs[1], gs[0].vs)
                        if v['rec']
                    ]) /
                    float(len(self.pp.lists['rec'])) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                'TFs',
                'Number of transcription factors',
                'tfs',
                lambda gs: len([v for v in gs[0].vs if v['tf']]),
                lambda gs: len(
                    [v for v in self.specific(gs[1], gs[0].vs) if v['tf']]),
                'vcount'
            ),
            MultiBarplotParam(
                r'TFs [%]',
                'Percentage of transcription factors',
                'tfprop',
                lambda gs: len([v for v in gs[0].vs if v['tf']]) /
                    float(gs[0].vcount()) * 100.0,
                lambda gs: len([
                        v for v in self.specific(gs[1], gs[0].vs)
                        if v['tf']
                    ]) /
                    float(gs[0].vcount()) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                r'TFs [%]',
                'Percentage of all human TFs covered',
                'tfcov',
                lambda gs: len([v for v in gs[0].vs if v['tf']]) /
                    float(len(self.pp.lists['tf'])) * 100.0,
                lambda gs: len([
                        v for v in self.specific(gs[1], gs[0].vs)
                        if v['tf']
                    ]) /
                    float(len(self.pp.lists['tf'])) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                'Kinases',
                'Number of kinases',
                'kinases',
                lambda gs: len([v for v in gs[0].vs if v['kin']]),
                lambda gs: len(
                    [v for v in self.specific(gs[1], gs[0].vs) if v['kin']]),
                'vcount'
            ),
            MultiBarplotParam(
                r'Kinases [%]',
                'Percentage of kinases',
                'kinprop',
                lambda gs: len([v for v in gs[0].vs if v['kin']]) /
                    float(gs[0].vcount()) * 100.0,
                lambda gs: len([
                        v for v in self.specific(gs[1], gs[0].vs)
                        if v['kin']
                    ]) /
                    float(gs[0].vcount()) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                r'Kinases [%]',
                'Percentage of all human kinases covered',
                'kincov',
                lambda gs: len([v for v in gs[0].vs if v['kin']]) /
                    float(len(self.pp.lists['kin'])) * 100.0,
                lambda gs: len([
                        v for v in self.specific(gs[1], gs[0].vs)
                        if v['kin']
                    ]) /
                    float(len(self.pp.lists['kin'])) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                r'Druggable proteins [%]',
                'Percentage of druggable proteins',
                'dgbprop',
                lambda gs: len([v for v in gs[0].vs if v['dgb']]) /
                    float(gs[0].vcount()) * 100.0,
                lambda gs: len([
                        v for v in self.specific(gs[1], gs[0].vs)
                        if v['dgb']
                    ]) /
                    float(gs[0].vcount()) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                r'Druggable proteins [%]',
                'Percentage of all human druggable proteins covered',
                'dgbcov',
                lambda gs: len([v for v in gs[0].vs if v['dgb']]) /
                    float(len(self.pp.lists['dgb'])) * 100.0,
                lambda gs: len([
                        v for v in self.specific(gs[1], gs[0].vs)
                        if v['dgb']
                    ]) /
                    float(len(self.pp.lists['dgb'])) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                r'Disease genes [%]',
                'Percentage of disease-gene associations covered',
                'discov',
                lambda gs:  len([v for v in gs[0].vs if v['dis']]) /
                    float(len(self.pp.lists['dis'])) * 100.0,
                lambda gs: len([
                        v for v in self.specific(gs[1], gs[0].vs)
                        if v['dis']
                    ]) /
                    float(len(self.pp.lists['dis'])) * 100.0,
                'vcount'
            ),
            #MultiBarplotParam(
                #'Complexes',
                #'Number of complexes covered',
                #'complexes',
                #lambda gs: len(self.pp.complexes_in_network(graph=gs[0])),
                #None,
                #'vcount'
            #),
            MultiBarplotParam(
                'E-S interactions',
                'Number of enzyme-substrate interactions covered',
                'ptmnum',
                lambda gs: sum(map(lambda e: len(e['ptm']), gs[0].es)),
                lambda gs: sum(
                        map(
                            lambda e: len(e['ptm']),
                            self.specific(gs[1], gs[0].es)
                        )
                    ),
                'vcount'
            ),
            MultiBarplotParam(
                'E-S interactions',
                (
                    'Number of interactions associated with '
                    'enzyme-substrate relationship'
                ),
                'havingptm',
                lambda gs: sum(map(lambda e: len(e['ptm']) > 0, gs[0].es)),
                lambda gs: sum(
                        map(
                            lambda e: len(e['ptm']) > 0,
                            self.specific(gs[1], gs[0].es)
                        )
                    ),
                'vcount'
            ),
            MultiBarplotParam(
                'Interactions with \n' + r'E-S relationship [%]',
                (
                    'Percentage of interactions associated with '
                    'enzyme-substrate relationship'
                ),
                'ptmprop',
                lambda gs: sum(
                        map(
                            lambda e: len(e['ptm']) > 0,
                            gs[0].es
                        )
                    ) /
                    float(gs[0].ecount()) * 100.0,
                lambda gs: sum(
                        map(
                            lambda e: len(e['ptm']) > 0,
                            self.specific(gs[1], gs[0].es)
                        )
                    ) / float(gs[0].ecount()) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                r'CGC genes [%]',
                (
                    'Percentage of COSMIC Cancer Gene Census '
                    'cancer drivers covered'
                ),
                'ccgccov',
                lambda gs: len(
                        set(gs[0].vs['name']) & set(self.pp.lists['cgc'])
                    ) /
                    float(len(self.pp.lists['cgc'])) * 100.0,
                lambda gs: len(
                        set(
                            map(
                                lambda v: v['name'],
                                self.specific(gs[1], gs[0].vs)
                            )
                        ) &
                        set(self.pp.lists['cgc'])
                    ) /
                    float(len(self.pp.lists['cgc'])) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                r'IntOGen genes [%]',
                'Percentage of IntOGen cancer drivers covered',
                'intocov',
                lambda gs: len(
                        set(gs[0].vs['name']) &
                        set(self.pp.lists['IntOGen'])
                    ) /
                    float(len(self.pp.lists['IntOGen'])) * 100.0,
                lambda gs: len(
                        set(
                            map(
                                lambda v: v['name'],
                                self.specific(gs[1], gs[0].vs)
                            )
                        ) &
                        set(self.pp.lists['IntOGen'])
                    ) /
                    float(len(self.pp.lists['IntOGen'])) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                r'Signaling proteins [%]',
                'Percentage of all human signaling proteins covered',
                'sigcov',
                lambda gs: len(
                        set(gs[0].vs['name']) & set(self.pp.lists['sig'])
                    ) / float(len(self.pp.lists['sig'])) * 100.0,
                lambda gs: len(
                        set(
                            map(
                                lambda v: v['name'],
                                self.specific(gs[1], gs[0].vs)
                            )
                        ) &
                        set(self.pp.lists['sig'])
                    ) / float(len(self.pp.lists['sig'])) * 100.0,
                'vcount'
            ),
            MultiBarplotParam(
                r'Signaling proteins [%]',
                'Percentage of signaling proteins',
                'sigpct',
                lambda gs: len(
                        set(gs[0].vs['name']) & set(self.pp.lists['sig'])
                    ) / float(gs[0].vcount()) * 100.0,
                lambda gs: len(
                    set(
                        map(
                            lambda v: v['name'],
                            self.specific(gs[1], gs[0].vs)
                        )
                    ) &
                    set(self.pp.lists['sig'])
                ) / float(gs[0].vcount()) * 100.0,
                'vcount'
            )
        ]

        self.scatterplots_settings = [
            ScatterplotParam(
                lambda gs: gs[0].vcount(),
                lambda gs: gs[0].ecount(),
                lambda gs: gs[0].density(),
                ScatterplotGraphicsParam(
                    ylim = [20.0, 200000.0],
                    xlim = [30.0, 15000.0],
                    ylog = True,
                    xlog = True,
                    xlab = 'Number of proteins',
                    ylab = 'Number of interacting pairs',
                    legtitle = 'Density',
                    title = (
                        'Number of proteins, '
                        'interactions and graph density'
                    )
                ),
                'vcount-ecount',
                True
            ),
            ScatterplotParam(
                lambda gs: len(set(gs[0].vs['name']) &
                    (set(self.pp.lists['dis']))),
                lambda gs: len(set(gs[0].vs['name']) &
                               set(self.pp.lists['cdv'])),
                lambda gs: gs[0].vcount(),
                ScatterplotGraphicsParam(
                    ylim = [0.0, 2000.0],
                    xlim = [30.0, 5500.0],
                    ylog = True,
                    xlog = True,
                    xlab = 'Number of\ndisease related proteins',
                    ylab = 'Number of cancer drivers',
                    legtitle = 'Total number of proteins',
                    title = 'Number of disease related proteins and cancer drivers',
                    legstrip = (3, None)
                ),
                'dis-cancer',
                False
            ),
            ScatterplotParam(
                lambda gs: len(set(gs[0].vs['name']) &
                               set(self.pp.lists['rec'])),
                lambda gs: len(set(gs[0].vs['name']) &
                               set(self.pp.lists['tf'])),
                lambda gs: gs[0].vcount(),
                ScatterplotGraphicsParam(
                    ylim = [0.5, 2000.0],
                    xlim = [5.0, 1200.0],
                    ylog = True,
                    xlog = True,
                    xlab = 'Number of receptors',
                    ylab = 'Number of transcription factors',
                    legtitle = 'Total number of proteins',
                    title = 'Number of receptors and TFs',
                    legstrip = (3, None)
                ),
                'tf-rec',
                False
            ),
            #ScatterplotParam(
                #lambda gs: len(self.pp.complexes_in_network(graph=gs[0])),
                #lambda gs: sum(map(lambda e: len(e['ptm']) > 0, gs[0].es)),
                #lambda gs: gs[0].ecount(),
                #ScatterplotGraphicsParam(
                    #ylim = [0.5, 5400.0],
                    #xlim = [3.0, 1500.0],
                    #ylog = True,
                    #xlog = True,
                    #xlab = 'Number of complexes',
                    #ylab = 'Number of\nenzyme-substrate relationships',
                    #legtitle = 'Total number of interactions',
                    #title = 'Number of complexes and enzyme-substrate relationships',
                    #legstrip = (3, None)
                #),
                #'comp-ptm',
                #False
            #)
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
        """
        Executes the entire workflow. Loads data and creates the figures
        and tables.
        """

        self.load_data()
        self.make_plots()

    def load_data(self):
        """
        Calls methods to load and preprocess the data.
        """
        
        # creating output directory
        self.set_outdir()
        # creating PyPath object
        self.init_pypath()
        
        self.load_omnipath()
        # load list of protein annotations (e.g. kinases, receptors, ...)
        self.load_protein_lists()
        # set the resource categories
        self.set_categories()
        # load the enzyme-substrate interactions
        self.pp.load_ptms2()
        
        if self.omnipath:
            
            self.omnipath.load_ptms2()
        
        # separate the network by resources
        self.separate()
        # assign colors for each resource for the multi-section barplots
        self.barplot_colors()
        # if we create any figure about literature curation
        # need to load the data from PubMed
        if (
            self.do_refs_composite or
            self.do_refs_journals_grid or
            self.do_refs_years_grid or
            self.do_curation_plot or
            self.do_htp_char
        ):
            self.load_pubmed_data()
        # for the directions plot compile the data about directions
        if self.do_dirs_stacked:
            self.get_dirs_data()
        # for consistency plots and tables compile the consistency data
        if (
            self.do_consistency_dedrogram or
            self.do_consistency_table
        ):
            self.inconsistency_data()

    def make_plots(self):
        """
        Calls methods to create figures and tables.
        """

        if self.do_main_table:
            self.main_table()
            if self.do_compile_main_table:
                self.latex_compile(self.main_table_fname)
        
        if self.do_curation_table:
            self.curation_table()
            if self.do_compile_curation_table:
                self.latex_compile(self.table2file)
        
        if self.do_multi_barplots:
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
                self.latex_compile(self.resource_list_fname)
        
        if self.do_consistency_dedrogram:
            self.make_consistency_dendrogram()
        
        if self.do_consistency_table:
            self.make_consistency_table()

    def set_outdir(self):
        """
        Creates the directory to save all output files.
        If no outdir given the ``name`` attribute will be used.
        """
        
        self.outdir = self.name if self.outdir is None else self.outdir
        
        if not os.path.exists(self.outdir):
            
            os.mkdir(self.outdir)
    
    
    def get_path(self, fname):
        """
        Returns the path of an output file by adding the path of the outdir.
        """
        
        return os.path.join(self.outdir, fname)
    
    
    def init_pypath(self):
        """
        Initializes a ``pypath.PyPath`` object using the resources from
        ``network_datasets``. If ``network_datasets`` is already a ``PyPath``
        object it will be used without any change. The ``PyPath`` object
        assigned to the ``pp`` attribute.
        """
        
        if isinstance(self.network_datasets, main.PyPath):
            
            self.pp = self.network_datasets
            
        else:
            
            self.pp = main.PyPath(9606)
            for netdata in self.network_datasets:
                self.pp.load_resources(getattr(data_formats, netdata))
    
    
    def load_omnipath(self):
        """
        Loads OmniPath from a previously compiled and saved pickle file.
        """
        
        self.omnipath = None
        
        if self.omnipath_pickle:
            
            self.omnipath = main.PyPath()
            self.omnipath.init_network(pfile = self.omnipath_pickle)
    
    
    def load_protein_lists(self):
        """
        Loads protein annotations to vertex attributes of the network object.
        By default these are kinases, TFs, receptors, druggability, disease
        relatedness.
        """
        
        # transcription factors
        self.pp.graph.vs['tf'] = [
            v['name'] in self.annot.annots['TFcensus']
            for v in self.pp.graph.vs
        ]
        self.pp.lists['tf'] = set(self.annot.annots['TFcensus'].annot.keys())
        
        # kinases
        self.pp.graph.vs['kin'] = [
            v['name'] in self.annot.annots['Kinases']
            for v in self.pp.graph.vs
        ]
        self.pp.lists['kin'] = set(self.annot.annots['Kinases'].annot.keys())
        
        # disease related genes
        self.pp.graph.vs['dis'] = [
            v['name'] in self.annot.annots['DisGeNet']
            for v in self.pp.graph.vs
        ]
        self.pp.lists['dis'] = set(self.annot.annots['DisGeNet'].annot.keys())
        
        # druggable proteins
        self.pp.graph.vs['dgb'] = [
            v['name'] in self.annot.annots['DGIdb']
            for v in self.pp.graph.vs
        ]
        self.pp.lists['dgb'] = set(self.annot.annots['DGIdb'].annot.keys())
        
        # receptors
        self.pp.graph.vs['rec'] = [
            v['name'] in self.intercell.classes['receptor']
            for v in self.pp.graph.vs
        ]
        self.pp.lists['rec'] = self.intercell.classes['receptor']
        
        # signaling proteins
        goannot = go.get_db()
        sig = goannot.select('signaling')
        
        self.pp.graph.vs['sig'] = [
            v['name'] in sig
            for v in self.pp.graph.vs
        ]
        self.pp.lists['sig'] = sig
        
        # cosmic cancer drivers
        self.pp.graph.vs['cgc'] = [
            v['name'] in self.annot.annots['CancerGeneCensus']
            for v in self.pp.graph.vs
        ]
        self.pp.lists['cgc'] = (
            set(self.annot.annots['CancerGeneCensus'].annot.keys())
        )
        
        self.pp.graph.vs['cgc'] = [
            v['name'] in self.pp.lists['cgc']
            for v in self.pp.graph.vs
        ]
        
        # IntOGen cancer drivers
        self.pp.graph.vs['IntOGen'] = [
            v['name'] in self.annot.annots['IntOGen']
            for v in self.pp.graph.vs
        ]
        self.pp.lists['IntOGen'] = (
            set(self.annot.annots['IntOGen'].annot.keys())
        )
        
        # all cancer drivers
        self.pp.graph.vs['cdv'] = [
            (
                v['name'] in self.annot.annots['IntOGen'] or
                v['name'] in self.pp.lists['cgc']
            )
            for v in self.pp.graph.vs
        ]
        self.pp.lists['cdv'] = self.pp.lists['IntOGen'] | self.pp.lists['cgc']
        
        self.pp.lists['proteome'] = dataio.all_uniprots(swissprot = True)
    
    
    def set_categories(self):
        """
        Sets the resource categories as a vertex and edge attribute on the
        network.
        """
        
        self.pp.set_categories()
    
    
    def separate(self):
        """
        Separates the network by resource category and resource.
        Also defines important variables for category names and their order
        which later will be used by the multi-section barplots.
        """
        
        # separated by resource
        self.sep = self.pp.separate()
        # separated by category
        self.csep = self.pp.separate_by_category()
        
        # dict of resource names to category names
        self.cats = dict(
            (c[0], db_categories.catnames[cc])
            for c in iteritems(db_categories.categories)
            for cc in c[1]
            if not self.only_categories or cc in self.only_categories
        )
        
        # tuples representing the total per each category
        self.cats.update(
            dict(
                (
                    ('All', c[0]),
                    cname
                )
                for c, cname in iteritems(db_categories.catnames)
                if c[0] in self.pp.has_cats and (
                    not self.only_categories or
                    c[0] in self.only_categories
                )
            )
        )
        
        # order of the categories on the multi section barplots
        self.cat_ordr = [
            c
            for c in self.cat_ordr_default # as given in the settings
            # if the category presents in the network
            if db_categories.catletters[c] in self.pp.has_cats and (
                # and it is allowed in the settings
                self.only_categories is None or
                db_categories.catletters[c] in self.only_categories
            )
        ]

    def fisher_tests(self):
        """
        Does Fisher tests for the enrichment of protein categories in
        resources and resource categories.
        """

        self._log('Doing Fisher tests, writing results to `%s`' %
                     self.get_path(self.fisher_file))
        
        with open(self.get_path(self.fisher_file), 'w') as fi:

            for attr, name in self.fisher:
                
                # contingency table
                cont = np.array(
                    [
                        [
                            # size of the proteome
                            len(self.pp.lists['proteome']),
                            # vertices in the network
                            self.pp.graph.vcount()
                        ],
                        [
                            # proteins in the protein category
                            len(self.pp.lists[attr]),
                            # proteins of category in the network
                            len([1 for v in self.pp.graph.vs if v[attr]])
                        ]
                     ]
                )
                
                fi.write('%s:' % name)
                fi.write('\t%s\t%s\n' % stats.fisher_exact(cont))

    def get_data(self, fun, attr = None):
        """
        Executes a method for each network and creates a list of
        tuples of resource names and results from the method.
        The result is either assigned to the attribute ``attr`` or returned
        of ``attr`` is ``None``.
        """
        
        if attr is None or not hasattr(self, attr):
            
            result = \
                list(
                    zip(
                        *[
                            (
                                s,
                                fun((self.sep[s], s))
                            )
                            for s in sorted(self.pp.sources)
                        ] +
                        [
                            (
                                ('All', c[0]), # instead of resource name
                                fun((self.csep[c[0]], c[0]))
                            )
                            for c in
                            sorted(self.csep.keys())
                        ]
                    )
                )
            
            if attr is None:
                
                return result
                
            else:
                
                setattr(self, attr, result)

    def specific(self, s, seq):
        """
        Returns the elements from a sequence of vertices or edges which
        are specific for a resource within a resource category.
        """
        
        return [
            w # vertex or edge
            for w in seq
            if (
                (
                    s in db_categories.categories and # the resource
                                                     # has category
                    s in w['sources'] and # the element belongs
                                          # to the resource
                    len(
                        w['sources'] & # set of resources for the element
                        getattr(data_formats, db_categories.categories[s])
                            # set of resources in the category
                    ) == 1
                ) or (
                    s in w['cat'] and # if not a resource
                                         # but a category
                    len(w['cat']) == 1 # element belongs only to this
                                       # category
                )
            )
        ]

    def make_multi_barplots(self):
        """
        Creates multi-section barplots. A section (subplot) created for
        each resource category. Generates plots according to settings in
        the ``barplot_settings`` attribute.
        """
        
        # keeping the plot objects in dict
        self.multi_barplots = {}
        
        for par in self.barplots_settings:
            
            _ = sys.stdout.write('\t:: Plotting %s\n' % par[1])
            
            # attribute name for y variable
            data_attr = 'data_%s' % par.name
            
            # heights of the bars
            self.get_data(par.method, data_attr)
            
            # data for the shaded area heights
            if par.smethod is not None:
                
                # attribute name for 2nd y variable
                data_attr2 = 'data_%s_2' % par.name
                self.get_data(par.smethod, data_attr2)
            
            data  = getattr(self, data_attr)
            data2 = getattr(self, data_attr2)
            
            # ordering of the columns
            ordr = self.vcount_ordr if par.order == 'vcount' else par.order
            
            # csv file to export the data
            csvname = self.get_path('%s-by-db_%s.csv' % (par.name, self.name))
            
            with open(csvname, 'w') as csv:
                
                ddata  = dict(zip(*data))
                ddata2 = dict(zip(*data2))
                
                csv.write('Label;%s;%s\n' % (par.ylab, par.ylab))
                
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
            
            # creating the barplot object
            self.multi_barplots[par.name] = (
                
                plot.MultiBarplot(
                    data[0], # resource names
                    data[1], # values for each resource
                    categories = self.cats, # resource to category dict
                    color = self.labcol,
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
                
            )

    def barplot_colors(self):
        """
        Creates lists of color codes for multi-section barplots according
        to the current set of resources and categories.
        """
        
        # size of resources and resource categories
        self.data_protein_counts = \
            list(
                zip(
                    *[
                        (
                            s,
                            len([
                                v for v in self.pp.graph.vs
                                if s in v['sources']
                            ])
                        )
                        for s in sorted(self.pp.sources)
                    ] +
                    [
                        (
                            ('All', c),
                            self.csep[c].vcount()
                        )
                        for c in sorted(self.csep.keys())
                    ]
                )
            )
        
        # colors of the bars
        self.labcol = [
            self.ccolors[db_categories.categories[lab]]
                if lab in db_categories.categories else
            self.ccolors[lab[1]] # for resource categories
            for lab in self.data_protein_counts[0]
        ]
        
        # colors of the shaded parts
        self.labcol2 = [
            self.ccolors2[db_categories.categories[lab]]
                if lab in db_categories.categories else
            self.ccolors2[lab[1]] # for resource categories
            for lab in self.data_protein_counts[0]
        ]

    def make_coverage_groups_plot(self):
        """
        Creates multi-section barplot showing the coverage of resources
        and resource categories on a group of proteins.
        """

        self.coverage_groups = (
            
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
            
        )

    def get_coverage_groups_data(self):
        """
        Creates data about the coverage of resources and resource categories
        on various groups of proteins: receptors, kinases, TFs and
        druggable proteins.
        """
        
        covdata = list(
            map(
                lambda fun:
                    # dict of resources and coverage
                    dict(zip(*self.get_data(fun, None))),
                # methods to calculate the coverage in percent
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
        
        # labels for coverage groups: resource names
        # and tuples for resource categories
        self.labels_coverage_groups = list(covdata[0].keys())
        
        # coverage data:
        # list of lists of percentages
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
        
        # labels for protein categories
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
        """
        Creates a ``plot.MultiBarplot`` object which orders by values
        on the y axis, in this case the number of proteins.
        This way obtains the order of resources and resource categories.
        Assigns the order to the ``vcount_ordr`` attribute.
        """
        
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
                order='y',
            )
        
        # we anyways create this figure later
        # now we remove it
        # os.remove(self.pdf_vcount_order)
        
        self.vcount_ordr = self.vcount_ordr_barplot.x
    
    def make_scatterplots(self):
        """
        Creates scatterplots showing pairs of variables.
        Each point represents one resource.
        """
        'xmethod', 'ymethod', 'smethod', 'gparam', 'name', 'sdiv'
        for par in self.scatterplots_settings:

            _ = sys.stdout.write('\t:: Plotting %s\n' % par.gparam.title)

            xattr = 'data_%s_x' % par.name
            yattr = 'data_%s_y' % par.name
            sattr = 'data_%s_s' % par.name
            lattr = 'labels_scatterplot_%s' % par.name

            self.get_data(par.xmethod, xattr)
            self.get_data(par.ymethod, yattr)
            self.get_data(par.smethod, sattr)

            x = getattr(self, xattr)
            y = getattr(self, yattr)
            s = getattr(self, sattr)

            setattr(
                self,
                lattr,
                list(
                    filter(
                        lambda l: not isinstance(l, tuple),
                        getattr(self, xattr)[0]
                    )
                )
            )

            labels = getattr(self, lattr)

            x = dict(zip(*x))
            y = dict(zip(*y))
            s = dict(zip(*s))

            x = list(map(lambda l: x[l], labels))
            y = list(map(lambda l: y[l], labels))
            s = np.array(list(map(lambda l: s[l], labels)))

            if par.sdiv:
                s = s / 1.0

            colors = list(
                map(
                    lambda l:
                        self.ccolors2[db_categories.categories[l]],
                    labels
                )
            )
            
            color_labels = []
            for c in ['p', 'l', 'm', 'i', 'r']:
                
                if (
                    c in self.pp.has_cats and
                    (
                        self.only_categories is None or
                        c in self.only_categories
                    )
                ):
                    
                    color_labels.append(
                        (db_categories.catnames[c], self.ccolors2[c])
                    )

            for i, l in enumerate(labels):
                if type(l) is tuple:
                    labels[i] = db_categories.catnames[l[1]]

            csvname = self.get_path('%s_%s.csv' % (par.name, self.name))

            with open(csvname, 'w') as csv:

                csv.write('Label;%s;%s;%s\n' % (
                    par.gparam.xlab, par.gparam.ylab, par.gparam.legtitle
                ))
                csv.write(
                    '\n'.join(
                        map(
                            lambda l: ';'.join(map(str, l)),
                            zip(labels, x, y, s)
                        )
                    )
                )

            sp = plot.ScatterPlus(
                x,
                y,
                size = s,
                min_size = 30,
                max_size = 1000,
                labels = labels,
                color = colors,
                color_labels = color_labels,
                fname = self.get_path('%s_%s.pdf' % (par.name, self.name)),
                **par.gparam._asdict(),
            )

    def all_ptms_list(self):
        """
        Obtains a list of PTM objects from all databases.
        Creates a dict with resource names as keys and
        lists of PTMs as values.
        Result assigned to the ``ptms`` attribute.
        """

        self.ptms = {
            #'DEPOD': net.load_depod_dmi(return_raw = True),
            'Signor': self.pp.load_signor_ptms(return_raw=True),
            'Li2012': self.pp.load_li2012_ptms(return_raw=True),
            'HPRD': self.pp.load_hprd_ptms(return_raw=True),
            'MIMP': self.pp.load_mimp_dmi(return_raw=True),
            'PhosphoNetworks': self.pp.load_pnetworks_dmi(return_raw=True),
            'PhosphoELM': self.pp.load_phosphoelm(return_raw=True),
            'dbPTM': self.pp.load_dbptm(return_raw=True),
            'PhosphoSite': self.pp.load_psite_phos(return_raw=True),
        }

    def make_ptms_barplot(self):
        """
        Creates a barplot with enzyme-substrate resources on the x axis
        and number of interactions on the y axis.
        Each resource represented by two columns, one showing the total,
        other the specific interactions.
        """

        self.pp.uniq_ptms()

        self.ptms_all_in_network = set(
            uniqList(flatList(self.pp.graph.es['ptm'])))

        self.data_ptms = list(
            zip(*[(s, len(uniqList(m)), len(
                set(m) & self.ptms_all_in_network))
                  for s, m in self.ptms.items()] + [('All', len(
                      uniqList(flatList(self.ptms.values()))), len(
                          self.ptms_all_in_network))]))

        self.ptms_barplot = (
            
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
            
        )

    def make_htp_characteristics(self):
        """
        Creates a 5 plot figures showing the high-throughput characteristics
        of the network.
        The x axis represents the number of interactions curated from one
        single publication. Above this threshold we consider a study HTP.
        The y axis on each plot shows various characteristics of the network
        depending on the HTP threshold.
        """

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
        """
        Creates a TikZ figure of the history of resources.
        """

        self.history_tree = plot.HistoryTree(
            fname=self.get_path(self.history_tree_fname),
            compile=False,
            dotlineopacity=1.0)

    def compile_history_tree(self):
        """
        Compiles the history of resources figure by LaTeX.
        """

        self._log('Running `%s` on `%s`' %
                     (self.latex, self.get_path(self.history_tree_fname)))

        self.history_tree_latex_proc = subprocess.Popen(
            [
                self.latex,
                '-halt-on-error',
                '-output-directory',
                self.outdir,
                self.get_path(self.history_tree_fname),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.history_tree_latex_output, self.history_tree_latex_error = \
            self.history_tree_latex_proc.communicate()
        self.history_tree_latex_return = \
            self.history_tree_latex_proc.returncode

        self._log('LaTeX %s' % (
            'compiled successfully'
                if self.history_tree_latex_return == 0 else
            'compilation failed')
        )

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
        
        self._log('Plotting references by database')
        bydb_ordr = reversed(
            self.pubmeds.database.value_counts().sort_values().index
        )
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
                title = (
                    'Number of references by resources '
                    'in %s databases' % self.title.capitalize()
                ),
                color=['#44AA99'] + ['#88CCEE'] *
                len(self.pubmeds.database.unique()),
                figsize=(12, 9),
                desc=False
            )

    def make_refs_by_year(self):
        
        self._log('Plotting references by year')
        counts = dict(self.pubmeds.year.value_counts())
        ecounts = dict(self.pubmeds_earliest.year.value_counts())
        years = np.arange(1970, max(self.pubmeds.year) + 1)
        #years = years[list(years).index(1970):]
        values = np.array(
            list(map(lambda y: counts[y] if y in counts else 0.0, years)))
        evalues = np.array(
            list(map(lambda y: ecounts[y] if y in ecounts else 0.0, years)))
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
                title= (
                    'Number of references by year '
                    'in %s databases' % self.title
                ),
                legend_font={'size': 'x-large'},
                color=['#77AADD', '#88CCAA'],
                figsize=(8, 5),
                group_labels=['All references', 'Earliest references'],
                desc=False,
                legloc=2
            )

    def make_refs_by_journal(self):

        self._log('Plotting references by journal')

        for i in [50, 100]:
            byj_ordr = list(
                reversed(
                    self.pubmeds.journal.value_counts().sort_values().index
                )
            )[:i]
            byj_vals = list(
                reversed(
                    self.pubmeds.journal.value_counts().sort_values()
                )
            )[:i]

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
                    title = (
                        'Number of references by journals '
                        'in %s databases' % self.title.capitalize()
                    ),
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
                fname = self.get_path(
                    self.dirs_stacked_fname % (
                        'all' if include_all else 'wo-all', self.name
                    )
                ),
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
                    sum([
                        sum([
                            s in
                            e['dirs'].positive_sources[e['dirs'].straight],
                            s in
                            e['dirs'].positive_sources[e['dirs'].reverse]
                        ])
                        for e in g.es
                    ]),
                    sum([
                        sum([
                            s in
                            e['dirs'].negative_sources[e['dirs'].straight],
                            s in
                            e['dirs'].negative_sources[e['dirs'].reverse]
                        ])
                        for e in g.es
                    ]),
                    sum([
                        sum([
                            s in (
                                (
                                    e['dirs'].sources[e['dirs'].straight] -
                                    e['dirs'].positive_sources[
                                        e['dirs'].straight
                                    ]
                                ) -
                                e['dirs'].negative_sources[e['dirs'].straight]
                            ),
                            s in (
                                (
                                    e['dirs'].sources[e['dirs'].reverse] -
                                    e['dirs'].positive_sources[
                                        e['dirs'].reverse
                                    ]
                                ) -
                                e['dirs'].negative_sources[e['dirs'].reverse]
                            )
                        ])
                        for e in g.es
                    ]),
                    sum([s in e['dirs'].sources['undirected'] for e in g.es])
                )
                for s, g in iteritems(self.sep)] +
                [(
                    'All',
                    sum([
                        e['dirs'].is_stimulation()
                        for e in self.pp.graph.es
                    ]),
                    sum([
                        e['dirs'].is_inhibition()
                        for e in self.pp.graph.es
                    ]),
                    sum([
                        not e['dirs'].is_stimulation() and
                        not e['dirs'].is_inhibition() and
                        e['dirs'].is_directed() for e in self.pp.graph.es
                    ]),
                    sum([
                        not e['dirs'].is_stimulation() and
                        not e['dirs'].is_inhibition() and
                        not e['dirs'].is_directed() for e in self.pp.graph.es
                    ])
                )]
            ))

    def curation_table(self):

        self._log('Making curation statistics '
                     'LaTeX table, writing to file `%s`' %
                     self.get_path(self.table2file))
        self.pp.curation_tab(
            latex_hdr=True, fname=self.get_path(self.table2file))

        self._log('Making curation statistics '
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
        """
        Compiles the tex file ``fname`` using the LaTeX compiler defined in
        the ``latex`` attribute.
        """

        self._log(
            'Running `%s` on `%s`' % (self.latex, self.get_path(fname))
        )

        try:
            self.latex_proc = subprocess.Popen(
                [
                    self.latex,
                    '-halt-on-error',
                    '-output-directory',
                    self.outdir,
                    self.get_path(fname)
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.latex_output, self.latex_error = (
                self.latex_proc.communicate(timeout = self.latex_timeout)
            )
            self.latex_return = self.latex_proc.returncode
            
        except subprocess.TimeoutExpired:
            
            self.latex_return = 1
        
        self._log(
            'LaTeX %s' % (
                'compiled successfully'
                    if self.latex_return == 0 else
                'compilation failed'
            )
        )

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
            map(
                lambda s:
                (
                    s,
                    self.inco['directions_edges']['minor'][s] / float(
                        self.inco['directions_edges']['minor'][s] +
                        self.cons['directions_edges']['total'][s] +
                        0.0001
                    )
                ),
                self.pp.sources
            )
        )

        self.incis = dict(
            map(
                lambda s:
                (
                    s,
                    self.inco['signs_edges']['minor'][s] / float(
                        self.inco['signs_edges']['minor'][s] +
                        self.cons['signs_edges']['total'][s] +
                        0.0001
                    )
                ),
                self.pp.sources))
        
        self.incidp = dict(
            map(
                lambda s: (
                    s[0],
                    (
                        len(s[1]['total']) /
                        float(
                            len(s[1]['total']) +
                            len(
                                self.inc_raw['consistency'][
                                    'directions_edges'
                                ][s[0]]['total']
                            ) +
                            0.0001
                        )
                    )
                    if (
                        len(s[1]['total']) +
                        len(
                            self.inc_raw['consistency'][
                                'directions_edges'
                            ][s[0]]['total']
                        )
                    ) > 5
                    else np.nan
                ),
                iteritems(self.inc_raw['inconsistency']['directions_edges'])
            )
        )
        
        self.incidp = dict(
            map(
                lambda s: (
                    s[0],
                    (
                        len(s[1]['total']) /
                        float(
                            len(s[1]['total']) +
                            len(
                                self.inc_raw['consistency'][
                                    'directions_edges'
                                ][s[0]]['total']
                            ) +
                            0.0001
                        )
                    )
                ),
                iteritems(self.inc_raw['inconsistency']['directions_edges'])
            )
        )
        
        self.incdf = pd.DataFrame(
            np.vstack([
                np.array(
                    [
                        self.incidp[(s1, s2)]
                        for s1 in sorted(self.pp.sources)
                    ]
                )
                for s2 in sorted(self.pp.sources)
            ]),
            index=sorted(self.pp.sources),
            columns=sorted(self.pp.sources)
        )
        
        self.incdf = self.incdf.loc[
            (self.incdf.sum(axis=1) != 0),
            (self.incdf.sum(axis=0) != 0)
        ]
        
        for lab in (
            'HPRD', 'NetPath', 'Reactome', 'ACSN', 'PDZBase', 'NCI-PID',
        ):
            
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

        singles_undirected = dict(
            map(
                lambda s:
                (
                    s,
                    len(
                        list(
                            filter(
                                lambda e:
                                    s in e['sources'] and
                                    len(e['sources']) == 1,
                                self.pp.graph.es
                            )
                        )
                    )
                ),
                undirected
            )
        )
        
        singles_directed = dict(
            map(
                lambda s:
                (
                    s,
                    sum(
                        map(
                            lambda e:
                                sum([
                                    s in e['dirs'].sources_straight() and
                                    len(e['dirs'].sources_straight()) == 1,
                                    s in e['dirs'].sources_reverse() and
                                    len(e['dirs'].sources_reverse()) == 1
                                ]),
                            self.pp.graph.es
                        )
                    )
                ),
                directed
            )
        )
        
        singles_signed = dict(
            map(
                lambda s:
                (
                    s,
                    sum(
                        map(
                            lambda e:
                                sum([
                                    s in
                                    e['dirs'].positive_sources_straight() and
                                    len(
                                        e['dirs'].positive_sources_straight()
                                    ) == 1,
                                    s in
                                    e['dirs'].positive_sources_reverse() and
                                    len(
                                        e['dirs'].positive_sources_reverse()
                                    ) == 1,
                                    s in
                                    e['dirs'].negative_sources_straight() and
                                    len(
                                        e['dirs'].negative_sources_straight()
                                    ) == 1,
                                    s in
                                    e['dirs'].negative_sources_reverse() and
                                    len(
                                        e['dirs'].negative_sources_reverse()
                                    ) == 1
                                ]),
                            self.pp.graph.es
                        )
                    )
                ),
                signed
            )
        )
        
        scores_undirected = dict(
            map(
                lambda s:
                (
                    s[0],
                    s[1] / float(max(singles_undirected.values()))
                ),
                iteritems(singles_undirected)
            )
        )
        
        scores_directed = dict(
            map(
                lambda s:
                (
                    s[0],
                    (
                        s[1] / float(max(singles_directed.values())) +
                        (1 - self.incid[s[0]]) /
                        max(
                            map(
                                lambda x: 1 - x,
                                self.incid.values()
                            )
                        )
                    ) / 2.0
                ),
                iteritems(singles_directed)
            )
        )
        
        scores_signed = dict(
            map(
                lambda s: (
                    s[0],
                    (
                        s[1] / float(max(singles_signed.values())) +
                        (1 - self.incis[s[0]]) / max(
                            map(
                                lambda x: 1 - x,
                                self.incis.values()
                            )
                        ) +
                        (1 - self.incid[s[0]]) / max(
                            map(
                                lambda x: 1 - x,
                                self.incid.values()
                            )
                        )
                    ) / 3.0
                ),
                iteritems(singles_signed)
            )
        )
        
        tbl = (
            r'\begin{tabularx}{\textwidth}'
            '\n'
            r'{p{3.5cm}>{\raggedright\arraybackslash}'
            '\n'
            r'X>{\raggedright\arraybackslash}'
            '\n'
            r'X>{\raggedright\arraybackslash}'
            '\n'
            r'X>{\raggedright\arraybackslash}X}'
            '\n'
            r'\toprule'
            '\n'
            r'Resurce & Specific interactions & '
            '\n'
            r'Direction inconsistency & Effect inconsistency & '
            '\n'
            r'Combined score \\'
            '\n'
            r'\midrule'
            '\n'
            r'\multicolumn{5}{l}{Undirected resources} \\'
            '\n'
            r'\midrule'
            '\n'
        )
        
        for s, sc in (
            reversed(
                sorted(
                    iteritems(scores_undirected),
                    key = lambda s: s[1]
                )
            )
        ):
            
            tbl += (
                r'    %s & %u &   &  & %.04f \\' % (
                    s,
                    singles_undirected[s],
                    scores_undirected[s],
                ) + '\n'
            )
        
        tbl += (
            r'    \midrule'
            '\n'
            r'    \multicolumn{5}{l}{Directed resources} \\'
            '\n'
            r'    \midrule'
            '\n'
        )
        
        for s, sc in (
            reversed(
                sorted(
                    iteritems(scores_directed),
                    key = lambda s: s[1]
                )
            )
        ):
            
            tbl += r'    %s & %u &  %.04f &  & %.04f \\' % (
                s,
                singles_directed[s],
                self.incid[s],
                scores_directed[s]
            ) + '\n'
        
        tbl += (
            r'    \midrule'
            '\n'
            r'    \multicolumn{5}{l}{Signed resources} \\'
            '\n'
            r'    \midrule'
            '\n'
        )
        
        for s, sc in reversed(
                sorted(
                    iteritems(scores_signed), key=lambda s: s[1])):
            tbl += r'    %s & %u &  %.04f &  %.04f & %.04f \\' % (
                s,
                singles_signed[s],
                self.incid[s],
                self.incis[s],
                scores_signed[s],
            ) + '\n'
        
        tbl += (
            r'    \bottomrule'
            '\n'
            r'\end{tabularx}'
        )
        
        with open(self.get_path(self.consistency_table_fname), 'w') as fp:
            
            fp.write(tbl)
        
        self.consistency_table = tbl
    
    
    def load_pubmed_data(self):
        """
        Loads data about all references using the web query interface of
        NCBI PubMed (E-utils).
        """
        
        self.pubmeds, self.pubmeds_earliest = _refs.get_pubmed_data(
            self.pp,
            htp_threshold = None,
        )
        
        if self.omnipath:
            
            self.omnipath_pubmeds, self.omnipath_pubmeds_earliest = (
                _refs.get_pubmed_data(
                    self.omnipath,
                    htp_threshold = None,
                )
            )
