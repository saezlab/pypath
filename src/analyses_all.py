#!/usr/bin/python

import imp

from pypath import analysis
import pypath.data_formats as data_formats

causal = analysis.Workflow(
        name = 'causal',
        network_datasets = ['pathway', 'ptm'],
        intogen_file = '../../a2/src/pypath/data/intogene_cancerdrivers.tsv',
        do_curation_table = True,
        do_compile_curation_table = True,
        do_simgraphs = True,
        do_multi_barplots = True,
        do_coverage_groups = True,
        do_htp_char = True,
        do_ptms_barplot = True,
        do_scatterplots = True,
        do_history_tree = False,
        do_compile_history_tree = False,
        do_refs_journals_grid = True,
        do_refs_years_grid = True,
        do_dirs_stacked = True,
        do_refs_composite = True
    )

causal.run()

binary = analysis.Workflow(
        name = 'binary',
        network_datasets = ['pathway', 'ptm', 'interaction_htp'],
        intogen_file = '../../a2/src/pypath/data/intogene_cancerdrivers.tsv',
        do_curation_table = True,
        do_compile_curation_table = True,
        do_multi_barplots = True,
        do_coverage_groups = True,
        do_htp_char = True,
        do_ptms_barplot = True,
        do_scatterplots = True,
        do_history_tree = True,
        do_compile_history_tree = False,
        do_refs_journals_grid = True,
        do_refs_years_grid = True,
        do_dirs_stacked = True,
        do_refs_composite = True
    )

binary.run()

undir_ppi = analysis.Workflow(
        name = 'binary',
        network_datasets = ['interaction_htp'],
        intogen_file = '../../a2/src/pypath/data/intogene_cancerdrivers.tsv',
        do_curation_table = True,
        do_compile_curation_table = True,
        do_multi_barplots = True,
        do_coverage_groups = True,
        do_htp_char = True,
        do_ptms_barplot = True,
        do_scatterplots = True,
        do_history_tree = True,
        do_compile_history_tree = False,
        do_refs_journals_grid = True,
        do_refs_years_grid = True,
        do_dirs_stacked = True,
        do_refs_composite = True,
        multi_barplots_summary = False
    )

undir_ppi.run()

proc_desc = analysis.Workflow(
        name = 'procdesc',
        network_datasets = ['reaction'],
        intogen_file = '../../a2/src/pypath/data/intogene_cancerdrivers.tsv',
        do_curation_table = True,
        do_compile_curation_table = True,
        do_multi_barplots = True,
        do_coverage_groups = True,
        do_htp_char = True,
        do_ptms_barplot = True,
        do_scatterplots = True,
        do_history_tree = True,
        do_compile_history_tree = False,
        do_refs_journals_grid = True,
        do_refs_years_grid = True,
        do_dirs_stacked = True,
        do_refs_composite = True,
        multi_barplots_summary = False
    )

proc_desc.run()

all4 = analysis.Workflow(
        name = 'all',
        network_datasets = ['pathway', 'ptm', 'interaction_htp', 'reaction'],
        intogen_file = '../../a2/src/pypath/data/intogene_cancerdrivers.tsv',
        do_curation_table = True,
        do_compile_curation_table = True,
        do_multi_barplots = True,
        do_coverage_groups = True,
        do_htp_char = False,
        do_ptms_barplot = True,
        do_scatterplots = True,
        do_history_tree = True,
        do_compile_history_tree = False,
        do_refs_journals_grid = True,
        do_refs_years_grid = True,
        do_dirs_stacked = True,
        do_refs_composite = True
    )

all4.run()
