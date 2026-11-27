"""Metabolite classes: whether a chemical of the release counts as a metabolite.

One table, ``metabolite_class``, built from the published release as Parquet: one row for
each chemical that a source of the rule file claims, and for each member of its
connectivity group. The row names the class, the first source that claims the molecule,
and every resource that asserts it. A network dataset keeps the relations of metabolites
by naming the classes it admits. The product reads the release only; it resolves no
identifier and downloads nothing.
"""

from __future__ import annotations

from omnipath_subsets.metabolite_class.parquet_build import write_metabolite_class
from omnipath_subsets.metabolite_class.rules import RuleError, read_rules

__all__ = ["RuleError", "read_rules", "write_metabolite_class"]
