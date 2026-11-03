"""Entity canonicalization: reference library + unique-match resolution.

* :mod:`.policy`      – per entity class: which library, which namespaces vote, how to label
* :mod:`.identifiers` – identifier normalization shared by builder and matcher
* :mod:`.library`     – pin an immutable published reference
* :mod:`.match`       – set-based matching of observations against a library
* :mod:`.label`       – label choice for unmatched observations
* :mod:`.stats`       – resolution counters

See ``README.md`` in this directory for the conceptual description.
"""

from .identifiers import normalize_id, normalize_identifier, normalize_ns
from .label import assign_label
from .match import Match, LibraryMatcher, Node, Vote, votes_for
from .policy import (
    CHEMICAL,
    CHEMICAL_POLICY,
    COMPLEX_POLICY,
    CV_TERM_POLICY,
    GENE_PROTEIN,
    GENE_PROTEIN_POLICY,
    GENERIC_POLICY,
    POLICIES,
    EntityPolicy,
    get_policy,
)
from .stats import ResolutionTracker

__all__ = [
    "CHEMICAL",
    "CHEMICAL_POLICY",
    "COMPLEX_POLICY",
    "CV_TERM_POLICY",
    "GENE_PROTEIN",
    "GENE_PROTEIN_POLICY",
    "GENERIC_POLICY",
    "POLICIES",
    "EntityPolicy",
    "get_policy",
    "normalize_id",
    "normalize_identifier",
    "normalize_ns",
    "assign_label",
    "Match",
    "LibraryMatcher",
    "Node",
    "Vote",
    "votes_for",
    "ResolutionTracker",
]
