#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import bioigraph
from bioigraph import data_formats

font = 'HelveticaNeueLTStd Med Cn'
net = bioigraph.BioGraph()
net.init_network(pfile = 'cache/default_network.pickle')
#net.init_network({'arn': data_formats.best['arn']})
for pw in ['TGF', 'Notch']:
    for db in ['SignaLink3', 'Signor', 'SPIKE', 'InnateDB', 'BioGRID', 'NetPath', 'MPPI', 'DIP', 'CA1',
        'Macrophage', 'PhosphoSite', 'CancerCellMap', 'HPRD']:
        print pw, db
        nodes = [v.index for v in net.graph.vs if pw in v['slk_pathways']]
        dot = net.export_dot(nodes = nodes, save_graphics = '%s_%s.pdf'%(pw, db), prog = 'dot',
            save_dot = '%s_%s.dot'%(pw, db),
            graph_label = '%s pathway in %s'%(pw, db), return_object = True,
            font = font, 
            vertex_fontsize = 24.0,
            graph_fontsize = 42.0,
            vertex_shape = 'rect', 
            auto_edges = 'DIRECTIONS',
            edge_sources = [db], 
            dir_sources = [db], hide = True)
    print pw, 'all'
    dot = net.export_dot(nodes = nodes, save_graphics = '%s_all.pdf'%pw, prog = 'dot',
            save_dot = '%s_%s.dot'%(pw, 'all'),
            graph_label = '%s pathway in OmniPath'%pw, return_object = True,
            font = font, 
            vertex_fontsize = 24.0,
            graph_fontsize = 42.0,
            vertex_shape = 'rect',
            auto_edges = 'DIRECTIONS')

for pw in ['TGF', 'Notch', 'WNT']:
    for db in ['SignaLink3', 'Signor', 'SPIKE', 'InnateDB', 'BioGRID', 'NetPath', 'MPPI', 'DIP', 'CA1',
        'Macrophage', 'PhosphoSite', 'CancerCellMap', 'HPRD']:
        nodes = [v.index for v in net.graph.vs if pw in v['slk_pathways']]
        dot = net.export_dot(nodes = nodes, save_graphics = '%s_%s.pdf'%(pw, db), prog = 'dot',
            main_title = '%s pathway in %s'%(pw, db), return_object = True,
            save_dot = '%s_%s.dot'%(pw, db),
            label_font = font, 
            edge_sources = [db], 
            dir_sources = [db], hide = True)
    dot = net.export_dot(nodes = nodes, save_graphics = '%s_all.pdf'%pw, prog = 'dot',
            main_title = '%s pathway in OmniPath'%pw, return_object = True,
            label_font = font)
    
pdftk TGF_BioGRID.pdf  TGF_InnateDB.pdf TGF_NetPath.pdf TGF_MPPI.pdf TGF_DIP.pdf TGF_HPRD.pdf TGF_Macrophage.pdf TGF_PhosphoSite.pdf TGF_CA1.pdf TGF_CancerCellMap.pdf TGF_SPIKE.pdf TGF_SignaLink3.pdf  TGF_Signor.pdf TGF_all.pdf cat output TGF.pdf

pdftk WNT_BioGRID.pdf  WNT_InnateDB.pdf WNT_NetPath.pdf  WNT_MPPI.pdf WNT_DIP.pdf WNT_HPRD.pdf WNT_Macrophage.pdf WNT_PhosphoSite.pdf WNT_CA1.pdf WNT_CancerCellMap.pdf WNT_SPIKE.pdf WNT_SignaLink3.pdf  WNT_Signor.pdf WNT_all.pdf cat output WNT.pdf

pdftk Notch_BioGRID.pdf Notch_InnateDB.pdf Notch_NetPath.pdf  Notch_MPPI.pdf Notch_DIP.pdf Notch_HPRD.pdf Notch_Macrophage.pdf Notch_PhosphoSite.pdf Notch_CA1.pdf Notch_CancerCellMap.pdf Notch_SPIKE.pdf Notch_SignaLink3.pdf  Notch_Signor.pdf Notch_all.pdf cat output Notch.pdf