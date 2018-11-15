#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

"""
This is a standalone module with the only purpose of
building the tables for the webservice.
"""

from future.utils import iteritems

import imp
import copy
import pandas as pd

import pypath.ptm as ptm
import pypath.export as export
import pypath.main as main
import pypath.mapping as mapping
import pypath.data_formats as data_formats


class WebserviceTables(object):
    
    def __init__(
            self,
            only_human = False,
            outfile_interactions = 'omnipath_webservice_interactions.tsv',
            outfile_ptms = 'omnipath_webservice_ptms.tsv'
        ):
        
        self.only_human = only_human
        self.outfile_interactions = outfile_interactions
        self.outfile_ptms = outfile_ptms
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def main(self):
        
        self.init_mapper()
        self.interactions()
        self.ptms()
    
    def init_mapper(self):
        
        self.mapper = mapping.Mapper()
    
    def interactions(self):
        
        dataframes = []
        
        tfregulons = copy.deepcopy(data_formats.transcription)
        tfregulons['tfregulons'].inputArgs['levels'] = {
            'A', 'B', 'C', 'D', 'E'
        }
        
        param = (
            ('load_omnipath', {'kinase_substrate_extra': True}),
            ('init_network',  {'lst': tfregulons}),
            ('init_network',  {'lst': data_formats.mirna_target})
        )
        
        for to_call, kwargs in param:
            
            pa = main.PyPath()
            pa.mapper = self.mapper
            getattr(pa, to_call)(**kwargs)
            
            e = export.Export(pa)
            e.webservice_interactions_df()
            dataframes.append(e.df)
            
            if not self.only_human:
                
                graph_human = None
                
                for rodent in (10090, 10116):
                    
                    if pa.ncbi_tax_id == 9606:
                        
                        if pa.graph.ecount() < 100000:
                            
                            graph_human = copy.deepcopy(pa.graph)
                        
                    else:
                        
                        if graph_human:
                            
                            pa.graph = graph_human
                            pa.ncbi_tax_id = 9606
                            pa.genesymbol_labels(remap_all = True)
                            pa.update_vname()
                            
                        else:
                            
                            del e
                            del pa
                            pa = main.PyPath()
                            pa.mapper = self.mapper
                            getattr(pa, to_call)(**kwargs)
                    
                    pa.orthology_translation(rodent)
                    e = export.Export(pa)
                    e.webservice_interactions_df()
                    dataframes.append(e.df)
        
        del e
        del pa
        
        self.df_interactions = pd.concat(dataframes)
        self.df_interactions.to_csv(
            self.outfile_interactions,
            sep = '\t',
            index = False
        )
    
    def ptms(self):
        
        dataframes = []
        
        ptma = ptm.PtmAggregator(mapper = self.mapper)
        ptma.make_df(tax_id = True)
        dataframes.append(ptma.df)
        
        if not self.only_human:
            
            for rodent1, rodent2 in ((10090, 10116), (10116, 10090)):
                
                ptma = ptm.PtmAggregator(
                    ncbi_tax_id = rodent1,
                    map_by_homology_from = (9606, rodent2),
                    mapper = self.mapper
                )
                ptma.make_df(tax_id = True)
                dataframes.append(ptma.df)
        
        self.df_ptms = pd.concat(dataframes)
        self.df_ptms.to_csv(
            self.outfile_ptms,
            sep = '\t',
            index = False
        )
