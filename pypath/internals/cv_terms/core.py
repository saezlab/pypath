"""Base class for controlled vocabularies with metadata support.

This module provides the CvEnum base class that extends Python's Enum to include
metadata like definitions and source URLs while maintaining string behavior for
easy integration with existing code.
"""
from enum import Enum


class CvEnum(str, Enum):
    """Base class for controlled vocabularies with metadata.

    This class allows CV terms to behave as strings (for their accession values)
    while also carrying additional metadata like definitions and source URLs.

    Each CvEnum subclass can define a parent_cv_term class attribute that represents
    the parent term in the PSI-MI ontology hierarchy.

    Enum members can be provided either as plain strings (for terms where metadata
    comes from an external ontology, e.g. PSI-MI) or as tuples containing the
    accession followed by an optional definition and optional source URL (for
    OmniPath-specific extensions).

    Usage:
        class MyCV(CvEnum):
            parent_cv_term = "MI:XXXX"  # Parent term in PSI-MI ontology

            TERM_NAME = ("ACC:0001", "Definition of the term", "https://source.url")

        # Access as string
        print(MyCV.TERM_NAME)  # "ACC:0001"

        # Access metadata
        print(MyCV.TERM_NAME.definition)  # "Definition of the term"
        print(MyCV.TERM_NAME.source)      # "https://source.url"
        print(MyCV.parent_cv_term)        # "MI:XXXX"
    """

    # Class attribute to store parent CV term (to be overridden by subclasses)
    parent_cv_term: str | None = None

    def __new__(
        cls,
        accession: str,
        definition: str | None = None,
        source: str | None = None
    ):
        """Create a new CV term with metadata.

        Args:
            accession: The CV term accession (e.g., "MI:0326" or "OM:0001")
            definition: Human-readable definition of the term
            source: URL or reference to the term's source
        """
        obj = str.__new__(cls, accession)
        obj._value_ = accession
        obj.definition = definition
        obj.source = source
        return obj

    def __str__(self):
        """Return the accession value as a string."""
        return self.value
