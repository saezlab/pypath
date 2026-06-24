#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""
This is a standalone module with the only purpose of
building the tables for the webservice.
"""

from future.utils import iteritems

import os
import importlib as imp
import copy

import pandas as pd

import pypath.core.enz_sub as enz_sub
import pypath.core.complex as complex
import pypath.core.annot as annot
import pypath.core.intercell as intercell
import pypath.omnipath.export as export
import pypath.resources.data_formats as data_formats
import pypath.share.session as session_mod
import pypath.omnipath as omnipath


NO_HOMOLOGY_TRANSLATION = {
    'mirna_mrna',
    'lncrna_mrna',
    'tf_mirna',
    'small_molecule',
}


class WebserviceTables(session_mod.Logger):
    """
    Creates the data frames which the web service uses to serve the data from.
    """


    def __init__(
            self,
            only_human = False,
            build_dir = '.',
            outfile_interactions = 'omnipath_webservice_interactions.tsv',
            outfile_enz_sub = 'omnipath_webservice_enz_sub.tsv',
            outfile_complexes = 'omnipath_webservice_complexes.tsv',
            outfile_annotations = 'omnipath_webservice_annotations.tsv',
            outfile_intercell = 'omnipath_webservice_intercell.tsv',
            network_datasets = None,
        ):

        session_mod.Logger.__init__(self, name = 'websrvtab')
        self._log('WebserviceTables initialized.')

        self.only_human = only_human

        databases = (
            'interactions',
            'complexes',
            'enz_sub',
            'annotations',
            'intercell',
        )

        for db in databases:

            attr = 'outfile_%s' % db

            setattr(
                self,
                attr,
                os.path.join(
                    build_dir,
                    locals()[attr]
                )
            )

        os.makedirs(build_dir, exist_ok = True)

        self.network_datasets = (
            network_datasets or
            (
                'omnipath',
                'tf_target',
                'mirna_mrna',
                'tf_mirna',
                'lncrna_mrna',
                'small_molecule',
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
        self.enz_sub()
        self.complexes()
        self.annotations()
        self.intercell()


    def interactions(self):

        self._log('Building `interactions` data frame.')
        dataframes = []

        for dataset in self.network_datasets:

            self._log('Building `%s` interactions.' % dataset)

            netw = omnipath.db.get_db(dataset)

            exp = export.Export(netw)
            exp.webservice_interactions_df()
            dataframes.append(exp.df)

            if dataset not in NO_HOMOLOGY_TRANSLATION:

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
            omnipath.db.remove_db(dataset)

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

        import pypath.legacy.main as main

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


    def enz_sub(self):

        self._log('Building `enz_sub` data frame.')

        dataframes = []

        self._log('Building `enz_sub` data frame for organism `9606`.')
        enz_sub_a = omnipath.db.get_db('enz_sub')
        enz_sub_a.make_df(tax_id = True)
        dataframes.append(enz_sub_a.df)
        self._log(
            'Finished building `enz_sub` data frame for organism `9606`.'
        )
        omnipath.db.remove_db('enz_sub', ncbi_tax_id = 9606)

        if not self.only_human:

            for rodent in (10090, 10116):

                self._log(
                    'Building `enz_sub` data frame for '
                    'organism `%s`.' % rodent
                )

                enz_sub_a = omnipath.db.get_db(
                    'enz_sub',
                    ncbi_tax_id = rodent,
                )
                enz_sub_a.make_df(tax_id = True)
                dataframes.append(enz_sub_a.df)

                self._log(
                    'Finished building `enz_sub` data frame for '
                    'organism `%s`.' % rodent
                )

                omnipath.db.remove_db('enz_sub', ncbi_tax_id = rodent)
                del enz_sub_a

        self.df_enz_sub = pd.concat(dataframes)
        self.df_enz_sub.to_csv(
            self.outfile_enz_sub,
            sep = '\t',
            index = False
        )

        self._log('Data frame `enz_sub` has been exported to `%s`.' % (
            self.outfile_enz_sub,
        ))


    def complexes(self):

        self._log('Building `complexes` data frame.')

        co = omnipath.db.get_db('complex')

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

        an = omnipath.db.get_db('annotations')

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

        i = omnipath.db.get_db('intercell')

        i.make_df()

        self.df_intercell = i.df
        self.df_intercell.to_csv(
            self.outfile_intercell,
            sep = '\t',
            index = False,
        )

        del i
        omnipath.db.remove_db('intercell')
        omnipath.db.remove_db('complex')
        omnipath.db.remove_db('annotations')

        self._log('Data frame `intercell` has been exported to `%s`.' % (
            self.outfile_intercell,
        ))
