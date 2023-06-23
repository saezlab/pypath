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
Contains auxiliary functions for preparation of arguments for building built
in databases. When users define custom, non built in databases, they can use
the functions from here or define and provide their own functions from their
own code.
"""

import copy

import pypath.share.settings as settings
import pypath.resources.network as netres


def curated_ppi_resources():
    """
    Returns a resource set which more or less corresponds to the literature
    curated activity flow resources. It is an union of the literature curated
    activity flow and enzyme-substrate resources.
    """

    resources = copy.deepcopy(netres.pathway)
    resources.update(copy.deepcopy(netres.enzyme_substrate))

    return resources


def tf_target_resources():
    """
    Returns the resource set for building the TF-target network dataset.
    """

    transcription = (
        netres.dorothea_expand_levels(
            resources = netres.transcription,
            levels = settings.get('tfregulons_levels'),
        )
            if settings.get('dorothea_expand_levels') else
        netres.transcription
    )

    return transcription
