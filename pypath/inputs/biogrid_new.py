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
BioGRID input module using the new YAML-based configuration system.
This demonstrates how to migrate existing input modules.
"""

from typing import List, Optional
from numbers import Number
import collections

from .input_module import InputModule, create_input_function


# For backward compatibility, we can still define the named tuple
BiogridPhysicalInteraction = collections.namedtuple(
    'BiogridPhysicalInteraction',
    (
        'partner_a',
        'partner_b',
        'pmid',
    ),
)


def biogrid_interactions(
    organism: int = 9606,
    htp_limit: Optional[Number] = 1,
    ltp: bool = True,
) -> List[BiogridPhysicalInteraction]:
    """
    Downloads and processes Physical multi-validated BioGRID interactions.
    Keeps only the "low throughput" interactions by default.
    Returns list of interactions.
    
    This is a compatibility wrapper around the new InputModule system.
    
    Args:
        organism: NCBI Taxonomy ID of organism.
        htp_limit: Exclude interactions only from references cited at more
            than this number of interactions.
        ltp: Whether to filter for low throughput only.
    """
    # Create input module with biogrid config
    module = InputModule('biogrid', 'biogrid_physical')
    
    # Update filters based on parameters
    if not ltp:
        # Remove the throughput filter
        module.source_config['filters'] = [
            f for f in module.source_config.get('filters', [])
            if f.get('field') != 'throughput'
        ]
    
    # Add organism filter if not human
    if organism != 9606:
        for f in module.source_config.get('filters', []):
            if f.get('field') in ('taxid_a', 'taxid_b'):
                f['value'] = str(organism)
    
    # Load raw data as DataFrame for post-processing
    df = module.load_processed(raw=True)
    
    # Create interaction records
    interactions = []
    ref_counts = collections.Counter()
    
    for _, row in df.iterrows():
        interaction = BiogridPhysicalInteraction(
            partner_a=row['partner_a'],
            partner_b=row['partner_b'],
            pmid=row['pmid'],
        )
        interactions.append(interaction)
        ref_counts[row['pmid']] += 1
    
    # Apply HTP limit filter if specified
    if htp_limit is not None:
        interactions = [
            i for i in interactions 
            if ref_counts[i.pmid] <= htp_limit
        ]
    
    return interactions


def biogrid_all_interactions(
    organism: int = 9606,
    htp_limit: Optional[Number] = 1,
) -> List[tuple]:
    """
    Downloads and processes all BioGRID interactions including genetic.
    
    Args:
        organism: NCBI Taxonomy ID of organism.
        htp_limit: Exclude interactions only from references cited at more
            than this number of interactions.
    """
    # Use the 'biogrid_all' configuration
    module = InputModule('biogrid', 'biogrid_all')
    
    # Update organism filter
    if organism != 9606:
        for f in module.source_config.get('filters', []):
            if f.get('field') in ('taxid_a', 'taxid_b'):
                f['value'] = str(organism)
    
    # Load and process
    df = module.load_processed(raw=True)
    
    # Define named tuple for all interactions
    BiogridInteraction = collections.namedtuple(
        'BiogridInteraction',
        list(df.columns)
    )
    
    # Create records
    interactions = []
    ref_counts = collections.Counter()
    
    for _, row in df.iterrows():
        interaction = BiogridInteraction(**row.to_dict())
        interactions.append(interaction)
        ref_counts[row['pmid']] += 1
    
    # Apply HTP limit
    if htp_limit is not None:
        interactions = [
            i for i in interactions 
            if ref_counts[i.pmid] <= htp_limit
        ]
    
    return interactions


# Alternative: Create functions using the factory
biogrid_physical_auto = create_input_function('biogrid', 'biogrid_physical')
biogrid_all_auto = create_input_function('biogrid', 'biogrid_all')


# Example of how to use the InputModule directly for more control
class BiogridModule(InputModule):
    """
    Extended BioGRID module with additional methods.
    """
    
    def __init__(self):
        super().__init__('biogrid', 'biogrid_physical')
    
    def get_interaction_counts(self) -> dict:
        """Get publication counts for interactions."""
        df = self.load_processed(raw=True)
        return df['pmid'].value_counts().to_dict()
    
    def get_high_throughput_pmids(self, threshold: int = 10) -> List[str]:
        """Get PMIDs that have more than threshold interactions."""
        counts = self.get_interaction_counts()
        return [pmid for pmid, count in counts.items() if count > threshold]