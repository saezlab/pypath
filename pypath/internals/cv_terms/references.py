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

    # PSI-MI standard reference types
    PUBMED = ("MI:0446", "PubMed database identifier", "https://pubmed.ncbi.nlm.nih.gov")
    PUBMED_CENTRAL = ("MI:1042", "PubMed Central database identifier", "https://www.ncbi.nlm.nih.gov/pmc")
    DOI = ("MI:0574", "Digital Object Identifier", "https://www.doi.org")
    BIORXIV = ("MI:2347", "bioRxiv preprint server", "https://www.biorxiv.org")


class CurationCv(CvEnum):
    """Curation and annotation metadata terms.

    Terms for capturing metadata about the curation process and
    supporting evidence for interactions and entities.
    """

    # OmniPath curation terms (OM:0400-0499 range)
    EVIDENCE_SENTENCE = (
        "OM:0400",
        "Sentence or text excerpt from literature supporting the interaction or annotation"
    )

    # PSI-MI standard curation terms
    COMMENT = (
        "MI:0612",
        "General comment, remark, or additional note from curator"
    )
