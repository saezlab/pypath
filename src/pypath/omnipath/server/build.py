#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2020 - EMBL
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

import pypath.core.ptm as ptm
import pypath.core.complex as complex
import pypath.core.annot as annot
import pypath.core.intercell as intercell
import pypath.omnipath.export as export
import pypath.legacy.main as main
import pypath.resources.data_formats as data_formats
import pypath.share.session as session_mod
import pypath.omnipath as omnipath


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
            network_datasets = None,
        ):
        
        session_mod.Logger.__init__(self, name = 'websrvtab')
        self._log('WebserviceTables initialized.')
        
        self.only_human = only_human
        self.outfile_interactions = outfile_interactions
        self.outfile_ptms = outfile_ptms
        self.outfile_complexes = outfile_complexes
        self.outfile_annotations = outfile_annotations
        self.outfile_intercell = outfile_intercell
        self.network_datasets = (
            network_datasets or
            (
                'omnipath',
                'tf_target',
                'mirna_mrna',
                'tf_mirna',
                'lncrna_mrna',
            )
        )
    
    
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
        
        for dataset in self.network_datasets:
            
            self._log('Building `%s` interactions.' % dataset)
            
            netw = omnipath.data.get_db(dataset)
            
            exp = export.Export(netw)
            exp.webservice_interactions_df()
            dataframes.append(exp.df)
            
            if dataset not in {'mirna_mrna', 'lncrna_mrna', 'tf_mirna'}:
                
                for rodent in (10090, 10116):
                    
                    self._log(
                        'Translating `%s` interactions to organism `%u`' % (
                            dataset,
                            rodent,
                        )
                    )
                    
                    rodent_netw = netw.homology_translate(rodent)
                    exp = export.Export(rodent_netw)
                    exp.webservice_interactions_df()
                    dataframes.append(exp.df)
                    
                    del rodent_netw
            
            del exp
            del netw
            omnipath.data.remove_db(dataset)
        
        self.df_interactions = pd.concat(dataframes)
        self.df_interactions.to_csv(
            self.outfile_interactions,
            sep = '\t',
            index = False
        )
        self._log('Data frame `interactions` has been exported to `%s`.' % (
            self.outfile_interactions,
        ))
    
    
    def interactions_legacy(self):
        
        self._log(
            'Building `interactions` data frame from '
            '`legacy.main.PyPath` object.'
        )
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
            'lncRNA-target': (
                'init_network',
                {'lst': data_formats.lncrna_target},
            )
        }
        
        for name, (to_call, kwargs) in iteritems(param):
            
            self._log('Building %s interactions.' % name)
            
            pa = main.PyPath()
            getattr(pa, to_call)(**kwargs)
            
            e = export.Export(pa)
            e.webservice_interactions_df()
            dataframes.append(e.df)
            
            if not self.only_human and name != 'lncRNA-target':
                
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
        
        ptma = omnipath.data.get_db('enz_sub')
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
        
        omnipath.data.remove_db('enz_sub')
        del ptma
        
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
        
        co = omnipath.data.get_db('complex')
        
        co.make_df()
        
        self.df_complexes = co.df
        
        del co
        
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
        
        an = omnipath.data.get_db('annotations')
        
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
        
        i = omnipath.data.get_db('intercell')
        
        i.make_df()
        i.add_classes_to_df()
        
        self.df_intercell = i.df
        self.df_intercell.to_csv(
            self.outfile_intercell,
            sep = '\t',
            index = False,
        )
        
        del i
        omnipath.data.remove_db('intercell')
        omnipath.data.remove_db('complex')
        omnipath.data.remove_db('annotations')
        
        self._log('Data frame `intercell` has been exported to `%s`.' % (
            self.outfile_intercell,
        ))
