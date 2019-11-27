#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2019 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
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

import importlib as imp
import copy

import pandas as pd

import pypath.ptm as ptm
import pypath.complex as complex
import pypath.annot as annot
import pypath.intercell as intercell
import pypath.export as export
import pypath.main as main
import pypath.data_formats as data_formats
import pypath.session_mod as session_mod


class WebserviceTables(session_mod.Logger):
    """
    Creates the data frames which the web service uses to serve the data from.
    """
    
    
    def __init__(
            self,
            only_human = False,
            outfile_interactions = 'omnipath_webservice_interactions.tsv',
            outfile_ptms = 'omnipath_webservice_ptms.tsv',
            outfile_complexes = 'omnipath_webservice_complexes.tsv',
            outfile_annotations = 'omnipath_webservice_annotations.tsv',
            outfile_intercell = 'omnipath_webservice_intercell.tsv',
            annot_args = None,
            complex_args = None,
            intercell_args = None,
            rebuild = None,
        ):
        
        session_mod.Logger.__init__(self, name = 'websrvtab')
        self._log('WebserviceTables initialized.')
        
        self.only_human = only_human
        self.outfile_interactions = outfile_interactions
        self.outfile_ptms = outfile_ptms
        self.outfile_complexes = outfile_complexes
        self.outfile_annotations = outfile_annotations
        self.outfile_intercell = outfile_intercell
        self.annot_args = annot_args or {}
        self.complex_args = complex_args or {}
        self.intercell_args = intercell_args or {}
        self.rebuild = rebuild or set()
    
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def main(self):
        
        self.interactions()
        self.ptms()
        self.complexes()
        self.annotations()
        self.intercell()
    
    
    def interactions(self):
        
        self._log('Building `interactions` data frame.')
        dataframes = []
        
        tf_target = copy.deepcopy(data_formats.transcription)
        tf_target['dorothea'].input_args['levels'] = {
            'A', 'B', 'C', 'D',
        }
        tf_target['dorothea'].must_have_references = False
        
        param = {
            'PPI': (
                'load_omnipath',
                {
                    'kinase_substrate_extra': True,
                    'ligand_receptor_extra': True,
                    'pathway_extra': True,
                },
            ),
            'TF-target': (
                'init_network',
                {'lst': tf_target},
            ),
            'miRNA-target': (
                'init_network',
                {'lst': data_formats.mirna_target},
            ),
        }
        
        for name, (to_call, kwargs) in iteritems(param):
            
            self._log('Building %s interactions.' % name)
            
            pa = main.PyPath()
            getattr(pa, to_call)(**kwargs)
            
            e = export.Export(pa)
            e.webservice_interactions_df()
            dataframes.append(e.df)
            
            if not self.only_human:
                
                graph_human = None
                
                for rodent in (10090, 10116):
                    
                    self._log(
                        'Translating %s interactions to organism `%u`' % (
                            name,
                            rodent,
                        )
                    )
                    
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
        self._log('Data frame `interactions` has been exported to `%s`.' % (
            self.outfile_interactions,
        ))
    
    
    def ptms(self):
        
        self._log('Building `ptms` data frame.')
        
        dataframes = []
        
        ptma = ptm.PtmAggregator()
        ptma.make_df(tax_id = True)
        dataframes.append(ptma.df)
        
        if not self.only_human:
            
            for rodent1, rodent2 in ((10090, 10116), (10116, 10090)):
                
                ptma = ptm.PtmAggregator(
                    ncbi_tax_id = rodent1,
                    map_by_homology_from = (9606, rodent2),
                )
                ptma.make_df(tax_id = True)
                dataframes.append(ptma.df)
        
        self.df_ptms = pd.concat(dataframes)
        self.df_ptms.to_csv(
            self.outfile_ptms,
            sep = '\t',
            index = False
        )
        
        self._log('Data frame `ptms` has been exported to `%s`.' % (
            self.outfile_ptms,
        ))
    
    
    def complexes(self):
        
        self._log('Building `complexes` data frame.')
        
        if 'complex' in self.rebuild:
            
            co = complex.ComplexAggregator(**self.complex_args)
            
        else:
            
            co = complex.get_db(**self.complex_args)
        
        co.make_df()
        
        self.df_complexes = co.df
        self.df_complexes.to_csv(
            self.outfile_complexes,
            sep = '\t',
            index = False,
        )
        
        self._log('Data frame `complexes` has been exported to `%s`.' % (
            self.outfile_complexes,
        ))
    
    
    def annotations(self):
        
        self._log('Building `annotations` data frame.')
        
        if 'annot' in self.rebuild:
            
            an = annot.AnnotationTable(**self.annot_args)
            
        else:
            
            an = annot.get_db(**self.annot_args)
        
        an.make_narrow_df()
        
        self.df_annotations = an.narrow_df
        self.df_annotations.to_csv(
            self.outfile_annotations,
            sep = '\t',
            index = False,
        )
        
        self._log('Data frame `annotations` has been exported to `%s`.' % (
            self.outfile_annotations,
        ))
    
    
    def intercell(self):
        
        self._log('Building `intercell` data frame.')
        
        if 'intercell' in self.rebuild:
            
            i = intercell.IntercellAnnotation(**self.intercell_args)
            
        else:
            
            i = intercell.get_db(**self.intercell_args)
        
        i.make_df()
        i.add_classes_to_df()
        
        self.df_intercell = i.df
        self.df_intercell.to_csv(
            self.outfile_intercell,
            sep = '\t',
            index = False,
        )
        
        self._log('Data frame `intercell` has been exported to `%s`.' % (
            self.outfile_intercell,
        ))
