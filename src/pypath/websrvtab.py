#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

"""
This is a standalone module with the only purpose of
building the tables for the webservice.
"""

from future.utils import iteritems

import imp
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
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def main(self):
        
        self.init_mapper()
        self.ptm()
        self.interactions()
    
    def init_mapper(self):
        
        self.mapper = mapping.Mapper()
    
    def interactions(self):
        
        dataframes = []
        
        param = (
            ('load_omnipath', {'kinase_substrate_extra': True}),
            ('init_network',  {'lst': data_formats.transcription}),
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
                
                graph_human = copy.deepcopy(pa.graph)
                
                # mouse
                pa.orthology_translation(10090)
                e = export.Export(pa)
                e.webservice_interactions_df()
                dataframes.append(e.df)
                
                # return to the human network
                # in order to translate also to rat
                pa.graph = graph_human
                pa.genesymbol_labels(remap_all = True)
                pa.update_vname()
                
                pa.orthology_translation(10116)
                e = export.Export(pa)
                e.webservice_interactions_df()
                dataframes.append(e.df)
        
        self.df_interactions = pd.concat(dataframes)
        self.df_interactions.to_csv(
            self.outfile_interactions,
            sep = '\t',
            index = False
        )
    
    def ptm(self):
        
        dataframes = []
        
        ptma = ptm.PtmAggregator(mapper = self.mapper)
        ptma.make_df()
        dataframes.append(ptma.df)
        
        if not self.only_human:
            
            for rodent1, rodent2 in ((10090, 10116), (10116, 10090)):
                
                ptma = ptm.PtmAggregator(
                    ncbi_tax_id = rodent1,
                    map_by_homology_from = (9606, rodent2),
                    mapper = self.mapper
                )
                ptma.make_df()
                dataframes.append(ptma.df)
        
        self.df_ptms = pd.concat(dataframes)
        self.df_ptms.to_csv(
            self.outfile_ptms,
            sep = '\t',
            index = False
        )
