#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
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
Client and parser for iMM1415 Mus musculus genome-scale metabolic model.

iMM1415 (Sigurdsson et al. 2010) is a mouse GEM distributed via BiGG
Models.  It shares the same BiGG JSON format as Recon3D, with mouse NCBI
Entrez Gene IDs in gene-reaction rules.

The :func:`imm1415_network` and :func:`imm1415_transporter_network`
functions yield :class:`GemInteraction` records compatible with
``pypath.inputs.metatlas`` and the ``omnipath-metabo`` GEM processor.

References:
    Sigurdsson MI, Jamshidi N, Steingrimsson E, et al.  A detailed
    genome-wide reconstruction of mouse metabolism based on human Recon 1.
    BMC Syst Biol. 2010;4:140. doi:10.1186/1752-0509-4-140
"""

from ._raw import imm1415_raw
from ._gem import (
    imm1415_genes,
    imm1415_metabolites,
    imm1415_network,
    imm1415_reactions,
    imm1415_transporter_network,
)

__all__ = [
    'imm1415_genes',
    'imm1415_metabolites',
    'imm1415_network',
    'imm1415_raw',
    'imm1415_reactions',
    'imm1415_transporter_network',
]
