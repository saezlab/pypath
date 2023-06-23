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

from typing import Union
import collections
import itertools

import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.utils.mapping as mapping


_log = session.Logger(name = 'cytosig_input')._log


CUSTOM_MAPPINGS = {
    'Activin A': ('P08476', 'INHBA'),
    'IL12': ('P29459', 'IL12A'),
    'IL36': ('Q9UHA7', 'IL36A'),
    'MCSF': ('P09603', 'CSF1'),
    'TWEAK': ('O43508', 'TNFSF12'),
}


def cytosig_df(long: bool = False) -> Union[pd.DataFrame, pd.Series]:
    """
    CytoSig core data is a matrix of cytokines vs. targets. Here this matrix
    is returned as a pandas data frame, either in a wide or long format.

    Args
        long: Convert the data frame to long format.

    Returns
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

    Returns
        Dict of sets of annotations.
    """

    def map_to_uniprot(genesymbol):

        uniprots = mapping.map_name(genesymbol, 'genesymbol', 'uniprot')

        if not uniprots and genesymbol in CUSTOM_MAPPINGS:

            uniprots = {CUSTOM_MAPPINGS.get(genesymbol)[0]}

        return uniprots


    cytosig = cytosig_df(long = True)

    record = collections.namedtuple(
        'CytosigAnnotation',
        ('cytokine', 'score', 'cytokine_genesymbol', 'target_genesymbol'),
    )
    result = collections.defaultdict(set)
    unmapped = set()

    for (target, cytokine), score in cytosig.items():

        u_target = map_to_uniprot(target)
        u_cytokine = map_to_uniprot(cytokine)

        if not u_cytokine:

            unmapped.add(cytokine)

        for u_t, u_c in itertools.product(u_target, u_cytokine):

            result[u_t].add(
                record(
                    cytokine = u_c,
                    score = score,
                    cytokine_genesymbol = cytokine,
                    target_genesymbol = target,
                )
            )

    if unmapped:

        _log(
            'Could not translate to UniProt IDs the following cytokines: ' +
            ', '.join(sorted(unmapped))
        )

    return dict(result)
