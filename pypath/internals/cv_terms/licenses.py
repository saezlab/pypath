"""License and update frequency controlled vocabularies.

This module contains CV terms for describing data source licenses and
update frequency categories.
"""
from .core import CvEnum


class LicenseCV(CvEnum):
    """Common license terms backed by OmniPath CV accessions.

    Describes the licensing terms under which data sources are available.
    """

    # OmniPath license terms (OM:0301-0399 range - note: different from membership roles)
    CC_BY_4_0 = (
        "OM:0501",
        "Creative Commons Attribution 4.0 International license",
        "https://creativecommons.org/licenses/by/4.0/"
    )
    CC0_1_0 = (
        "OM:0502",
        "Creative Commons Zero 1.0 Universal - public domain dedication",
        "https://creativecommons.org/publicdomain/zero/1.0/"
    )
    GPL_3_0 = (
        "OM:0503",
        "GNU General Public License v3.0",
        "https://www.gnu.org/licenses/gpl-3.0.html"
    )
    MIT = (
        "OM:0504",
        "MIT License - permissive free software license",
        "https://opensource.org/licenses/MIT"
    )
    ACADEMIC_FREE = (
        "OM:0505",
        "Free for academic use, restrictions may apply for commercial use"
    )
    UNSPECIFIED = (
        "OM:0599",
        "License not specified or unclear"
    )


class UpdateCategoryCV(CvEnum):
    """Update frequency categories backed by OmniPath CV accessions.

    Describes how frequently a data source is expected to be updated.
    """

    # OmniPath update category terms (OM:0401-0499 range)
    REGULAR = (
        "OM:0401",
        "Regular scheduled updates (e.g., monthly, quarterly, annually)"
    )
    IRREGULAR = (
        "OM:0402",
        "Irregular or occasional updates with no fixed schedule"
    )
    STATIC = (
        "OM:0403",
        "Static resource with no planned future updates"
    )
