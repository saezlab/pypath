"""Curated, build-time classification of canonical entities and predicates.

Each classifier reads a human-reviewable rule file (``*.yaml``) and populates a
derived classification column from data already in the database (source
membership, annotations). Run during ``derive`` — no reload required.
"""

from omnipath_build.classify.chemical_class import (
    ChemicalClassStats,
    classify_chemical_class,
)
from omnipath_build.classify.metabolic_domain import (
    MetabolicDomainStats,
    classify_metabolic_domain,
)
from omnipath_build.classify.interaction_class import (
    InteractionClassStats,
    classify_interaction_class,
)

__all__ = [
    'ChemicalClassStats',
    'classify_chemical_class',
    'MetabolicDomainStats',
    'classify_metabolic_domain',
    'InteractionClassStats',
    'classify_interaction_class',
]
