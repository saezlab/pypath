"""Reference and curation controlled vocabularies.

This module contains CV terms for describing references (publications, databases)
and curation metadata (evidence sentences, comments).
"""
from .core import CvEnum


class ReferenceTypeCv(CvEnum):
    """Reference source terms from PSI-MI.

    Describes the type or source of a bibliographic reference
    (e.g., PubMed, DOI, preprint servers).
    """

    parent_cv_term = "MI:0444"  # database citation - Database citation list names of databases commonly used to cross reference interaction data

    # PSI-MI standard reference types
    PUBMED = "MI:0446"
    PUBMED_CENTRAL = "MI:1042"
    DOI = "MI:0574"
    BIORXIV = "MI:2347"


class CurationCv(CvEnum):
    """Curation and annotation metadata terms.

    Terms for capturing metadata about the curation process and
    supporting evidence for interactions and entities.
    """

    parent_cv_term = ("OM:0400", "Curation term - Describes the results of a curation process and supporting evidence for interactions and entities.")  # OmniPath-specific term - no PSI-MI parent

    # OmniPath curation terms (OM:0400-0499 range)
    EVIDENCE_SENTENCE = (
        "OM:0401",
        "Sentence or text excerpt from literature supporting the interaction or annotation"
    )

    # PSI-MI standard curation terms
    COMMENT = "MI:0612"
