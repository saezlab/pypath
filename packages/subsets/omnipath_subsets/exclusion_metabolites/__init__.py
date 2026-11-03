"""Exclusion metabolites: chemical entities that mass spectrometry does not measure.

Cofactors, small molecules, ions and atoms stay in the prior knowledge. This
product publishes them as a classed list, so a consumer can leave them out of a
set without losing them from the build.

The list has a curated floor (seed files in the package) and a derived part
(ChEBI descendants and roles). The product reads one loaded release, resolves
no identifier and downloads nothing.
"""

from __future__ import annotations

from omnipath_subsets.exclusion_metabolites.build import build_exclusion_metabolites

__all__ = ["build_exclusion_metabolites"]
