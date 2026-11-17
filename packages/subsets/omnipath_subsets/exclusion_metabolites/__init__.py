"""Exclusion metabolites: chemical entities that mass spectrometry does not measure.

Cofactors, small molecules, ions and atoms stay in the prior knowledge. This
product publishes them as a classed list, so a consumer can leave them out of a
set without losing them from the build.

The list has a curated floor (seed files in the package) and a derived part
(ChEBI descendants and roles). The product reads the published release as Parquet, resolves
no identifier and downloads nothing.
"""

from __future__ import annotations

from omnipath_subsets.exclusion_metabolites.parquet_build import write_exclusion_metabolites

__all__ = ["write_exclusion_metabolites"]
