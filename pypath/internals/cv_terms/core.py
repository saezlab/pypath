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

    Usage:
        class MyCV(CvEnum):
            TERM_NAME = ("ACC:0001", "Definition of the term", "https://source.url")

        # Access as string
        print(MyCV.TERM_NAME)  # "ACC:0001"

        # Access metadata
        print(MyCV.TERM_NAME.definition)  # "Definition of the term"
        print(MyCV.TERM_NAME.source)      # "https://source.url"
    """

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
