"""Tests for controlled vocabulary terms."""

import pytest
from pypath.internals.cv_terms.core import CvEnum
from pypath.internals import cv_terms


def test_om_accession_uniqueness():
    """Ensure all OM accessions are unique across all CV enums.

    This test validates that no OmniPath accession (OM:XXXX) is assigned
    to multiple enum members, either within the same CV class or across
    different CV classes.
    """
    accession_map: dict[str, list[tuple[str, str]]] = {}

    # Collect all CV enum classes and their OM accessions
    for name in dir(cv_terms):
        obj = getattr(cv_terms, name)

        # Check if it's a CvEnum subclass (but not CvEnum itself)
        if (
            isinstance(obj, type)
            and issubclass(obj, CvEnum)
            and obj is not CvEnum
        ):
            for member in obj:
                accession = member.value

                # Only check OmniPath accessions (OM:XXXX)
                if accession.startswith('OM:'):
                    accession_map.setdefault(accession, []).append(
                        (obj.__name__, member.name)
                    )

    # Find duplicates
    duplicates = {
        accession: occurrences
        for accession, occurrences in accession_map.items()
        if len(occurrences) > 1
    }

    # Format error message if duplicates found
    if duplicates:
        error_lines = ["Found duplicate OmniPath accessions:"]
        for accession, occurrences in sorted(duplicates.items()):
            locations = ", ".join(f"{cls}.{member}" for cls, member in occurrences)
            error_lines.append(f"  {accession}: {locations}")

        pytest.fail("\n".join(error_lines))


def test_om_accession_format():
    """Ensure all OM accessions follow the correct format (OM:XXXX).

    Validates that OmniPath accessions use the proper format with 4 digits.
    """
    import re

    om_pattern = re.compile(r'^OM:\d{4}$')
    invalid_accessions: list[tuple[str, str, str]] = []

    # Collect all CV enum classes and check OM accession format
    for name in dir(cv_terms):
        obj = getattr(cv_terms, name)

        if (
            isinstance(obj, type)
            and issubclass(obj, CvEnum)
            and obj is not CvEnum
        ):
            for member in obj:
                accession = member.value

                # Check OM accessions for proper format
                if accession.startswith('OM:') and not om_pattern.match(accession):
                    invalid_accessions.append(
                        (obj.__name__, member.name, accession)
                    )

    # Format error message if invalid formats found
    if invalid_accessions:
        error_lines = ["Found OmniPath accessions with invalid format:"]
        for cls_name, member_name, accession in invalid_accessions:
            error_lines.append(f"  {cls_name}.{member_name}: '{accession}' (expected OM:XXXX)")

        pytest.fail("\n".join(error_lines))
