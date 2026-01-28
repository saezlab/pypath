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
Client and parser for Metabolic Atlas database.

Metabolic Atlas (https://metabolicatlas.org) provides access to
Genome-Scale Metabolic Models (GEMs) from various organisms and tissues.

This module provides functions to:
- List and download repository models from the Metabolic Atlas API
- Access standard GEMs from GitHub/GitLab repositories
- Parse reaction, metabolite, and gene annotations from GEM TSV files
"""

from ._api import (
    metatlas_models,
    metatlas_integrated_models,
    metatlas_model_files,
)
from ._git import metatlas_gems
from ._gem import (
    metatlas_gem_reactions,
    metatlas_gem_metabolites,
    metatlas_gem_genes,
    metatlas_gem_tsv,
)
from ._records import (
    MetatlasModel,
    MetatlasGem,
    MetatlasReaction,
    MetatlasMetabolite,
    MetatlasGene,
)
