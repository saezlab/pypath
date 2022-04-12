#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from typing import Union
import collections
import itertools

import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping


def cytosig_df(long: bool = False) -> Union[pd.DataFrame, pd.Series]:
    """
    CytoSig core data is a matrix of cytokines vs. targets. Here this matrix
    is returned as a pandas data frame, either in a wide or long format.

    Args:
        long: Convert the data frame to long format.

    Returns:
        The original matrix format data frame (rows are signature genes,
        columns are cytokines) by default; if `long` is `True`, a series
        with multi index is returned.
    """

    url = urls.urls['cytosig']['url']
    c = curl.Curl(url, large = True, silent = False)

    # bravo pandas!
    # I would've never expected that anything in pandas works this smooth
    df = pd.read_csv(c.fileobj, sep = '\t')

    if long:

        # well done, pandas.
        # your api and docs though are still an ugly mess
        df = df.stack()

    return df


def cytosig_annotations() -> dict:
    """
    CytoSig is a compendium of expression signatures from cytokine
    perturbation experiments.

    Returns:
        Dict of sets of annotations.
    """

    cytosig = cytosig_df(long = True)

    record = collections.namedtuple(
        'CytosigAnnotation',
        ('cytokine', 'score'),
    )
    result = collections.defaultdict(set)

    for (target, cytokine), score in cytosig.items():

        u_target = mapping.map_name(target, 'genesymbol', 'uniprot')
        u_cytokine = mapping.map_name(cytokine, 'genesymbol', 'uniprot')

        for u_t, u_c in itertools.product(u_target, u_cytokine):

            result[u_t].add(
                record(
                    cytokine = u_c,
                    score = score,
                )
            )

    return dict(result)
