#!/usr/bin/python

import imp

from pypath import analysis
import pypath.data_formats as data_formats

causal = analysis.Workflow(
        name = 'causal',
        network_datasets = ['pathway', 'ptm'],
        intogen_file = '/home/denes/documents/pw/a2/src/pypath/data/intogene_cancerdrivers.tsv',
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

'''
Comparison of OmniPath, ConsensusPathDB and PathwayCommons

This requires 7G of memory!!!
'''
import imp
import pypath
from pypath import data_formats

pa = pypath.PyPath()
pa.load_omnipath()

# nodes: 7988
# edges: 36541

pa.load_resources(data_formats.pathwaycommons1)
pa.load_resources({'CPDB': data_formats.interaction_misc['cpdb']})

pwc_cpdb = set(['PathwayCommons', 'CPDB'])

pwc_edges = sum(map(lambda e: 'PathwayCommons' in e['sources'], pa.graph.es)) # 572502
pwc_nodes = sum(map(lambda v: 'PathwayCommons' in v['sources'], pa.graph.vs)) # 17122

cpdb_edges = sum(map(lambda e: 'CPDB' in e['sources'], pa.graph.es)) # 183108
cpdb_nodes = sum(map(lambda v: 'CPDB' in v['sources'], pa.graph.vs)) # 15638

pwc_op_edges = sum(map(lambda e: 'PathwayCommons' in e['sources'] and len(e['sources'] - pwc_cpdb) > 0, pa.graph.es)) # 31792
pwc_op_nodes = sum(map(lambda v: 'PathwayCommons' in v['sources'] and len(v['sources'] - pwc_cpdb) > 0, pa.graph.vs)) # 7813

cpdb_op_edges = sum(map(lambda e: 'CPDB' in e['sources'] and len(e['sources'] - pwc_cpdb) > 0, pa.graph.es)) # 28961
cpdb_op_nodes = sum(map(lambda v: 'CPDB' in v['sources'] and len(v['sources'] - pwc_cpdb) > 0, pa.graph.vs)) # 7740

pwc_cpdb_edges = sum(map(lambda e: 'CPDB' in e['sources'] and 'PathwayCommons' in e['sources'], pa.graph.es)) # 161631
pwc_cpdb_nodes = sum(map(lambda v: 'CPDB' in v['sources'] and 'PathwayCommons' in v['sources'], pa.graph.vs)) # 15242

pwc_cpdb_op_edges = sum(map(lambda e: 'CPDB' in e['sources'] and 'PathwayCommons' in e['sources'] and len(e['sources'] - pwc_cpdb) > 0, pa.graph.es)) # 27291
pwc_cpdb_op_nodes = sum(map(lambda v: 'CPDB' in v['sources'] and 'PathwayCommons' in v['sources'] and len(v['sources'] - pwc_cpdb) > 0, pa.graph.vs)) # 7710

only_op_edges = sum(map(lambda e: 'CPDB' not in e['sources'] and 'PathwayCommons' not in e['sources'], pa.graph.es)) # 3700
only_op_nodes = sum(map(lambda v: 'CPDB' not in v['sources'] and 'PathwayCommons' not in v['sources'], pa.graph.vs)) # 145

total_edges = pa.graph.ecount() # 597679
total_nodes = pa.graph.vcount() # 17663

dens = pa.graph.density() # 0.003832
trans = pa.graph.transitivity_avglocal_undirected() # 0.3153
diam = pa.graph.diameter() # 7

pwc_directed = sum(map(lambda e: 'PathwayCommons' in e['dirs'].sources_straight(), pa.graph.es)) + \
    sum(map(lambda e: 'PathwayCommons' in e['dirs'].sources_reverse(), pa.graph.es)) # 90339

pwc_op_directed = pwc_directed = sum(map(lambda e: 'PathwayCommons' in e['dirs'].sources_straight() and len(e['dirs'].sources_straight()) > 1, pa.graph.es)) + \
    sum(map(lambda e: 'PathwayCommons' in e['dirs'].sources_reverse() and len(e['dirs'].sources_reverse()) > 1, pa.graph.es)) # 8833

only_op_directed = pwc_op_directed = pwc_directed = sum(map(lambda e: 'PathwayCommons' not in e['dirs'].sources_straight() and len(e['dirs'].sources_straight()) > 0, pa.graph.es)) + \
    sum(map(lambda e: 'PathwayCommons' not in e['dirs'].sources_reverse() and len(e['dirs'].sources_reverse()) > 0, pa.graph.es)) # 37511
