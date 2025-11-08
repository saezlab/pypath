"""Base class for controlled vocabularies with metadata support.

This module provides the CvEnum base class that extends Python's Enum to include
metadata like definitions and source URLs while maintaining string behavior for
easy integration with existing code.
"""
from enum import Enum, EnumMeta
from typing import Any

_NON_MEMBER_ATTRIBUTES = {'parent_cv_term'}


class CvEnumMeta(EnumMeta):
    """Enum metaclass that preserves metadata attributes without creating members."""

    def __new__(  # type: ignore[override]
        metacls,
        cls: str,
        bases: tuple[type, ...],
        classdict: dict[str, Any],
        **kwds: Any,
    ):
        metadata: dict[str, Any] = {}
        member_names = getattr(classdict, '_member_names', None)

        for name in _NON_MEMBER_ATTRIBUTES:
            if name in classdict:
                metadata[name] = classdict[name]
                del classdict[name]
                if member_names and name in member_names:
                    member_names.pop(name, None)

        enum_cls = super().__new__(metacls, cls, bases, classdict, **kwds)

        for name in _NON_MEMBER_ATTRIBUTES:
            value = metadata.get(name)
            if value is None and hasattr(enum_cls, name):
                value = getattr(enum_cls, name)
            setattr(enum_cls, name, value)

        return enum_cls


class CvEnum(str, Enum, metaclass=CvEnumMeta):
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
