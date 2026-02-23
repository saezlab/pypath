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
Client and parser for Recon3D genome-scale metabolic model.

Recon3D (Brunk et al. 2018) is a comprehensive human GEM maintained by
the Virtual Metabolic Human (VMH) project and distributed via BiGG Models.

Two data sources are supported:

- **BiGG JSON** (primary): full model with stoichiometry, gene rules, and
  metabolite/gene annotations including HMDB, ChEBI, and KEGG IDs.
- **VMH MATLAB** (supplementary): ``.mat`` file providing additional HMDB
  annotations parsed via :func:`gem_matlab_extract`.

The :func:`recon3d_network` function yields :class:`GemInteraction` records
compatible with ``pypath.inputs.metatlas`` and the ``omnipath-metabo`` GEM
processor.  Gene IDs are NCBI Entrez Gene IDs (integers as strings).

References:
    Brunk E, Sahoo S, Zielinski DC, et al. Recon3D enables a three-dimensional
    view of gene variation in human metabolism. Nat Biotechnol. 2018;36(3):272-281.
    doi:10.1038/nbt.4072
"""

from ._raw import (
    gem_matlab_extract,
    recon3d_raw,
    recon3d_raw_vmh,
)
from ._gem import (
    recon3d_genes,
    recon3d_metabolites,
    recon3d_network,
    recon3d_reactions,
)

__all__ = [
    'gem_matlab_extract',
    'recon3d_genes',
    'recon3d_metabolites',
    'recon3d_network',
    'recon3d_raw',
    'recon3d_raw_vmh',
    'recon3d_reactions',
]
