"""License and update frequency controlled vocabularies.

This module contains CV terms for describing data source licenses and
update frequency categories.
"""
from enum import StrEnum

from .core import CvEnum


class LicenseUse(StrEnum):
    """Whether a license permits a category of use of the data."""

    ALLOWED = 'allowed'
    REQUIRES_PERMISSION = 'requires_permission'  # a separate authorization is specified
    NOT_ALLOWED = 'not_allowed'
    UNKNOWN = 'unknown'  # the license terms do not settle the question


class LicenseCV(CvEnum):
    """Common license terms backed by OmniPath CV accessions.

    Describes the licensing terms under which data sources are available.

    Each term is (accession, definition, source URL, academic use, commercial use).
    The two uses are browsing classifications, subject to the exact license
    conditions. A software license alone does not classify imported data, and a
    vague commercial restriction does not settle commercial use: both stay unknown.
    A noncommercial license does not make every use by an academic institution
    noncommercial.
    """

    parent_cv_term = ("OM:0500", "License term", "Describes the licensing terms under which data sources are available.")  # OmniPath-specific term - no PSI-MI parent

    def __new__(
        cls,
        accession: str,
        definition: str,
        source: str | None,
        academic: LicenseUse,
        commercial: LicenseUse,
    ):
        obj = CvEnum.__new_member__(cls, accession, definition, source)
        obj.academic = academic
        obj.commercial = commercial
        return obj

    # OmniPath license terms (OM:0301-0399 range - note: different from membership roles)
    CC_BY_4_0 = (
        "OM:0501",
        "Creative Commons Attribution 4.0 International license",
        "https://creativecommons.org/licenses/by/4.0/",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    CC_BY_3_0 = (
        "OM:0514",
        "Creative Commons Attribution 3.0 International license",
        "https://creativecommons.org/licenses/by/3.0/",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    CC0_1_0 = (
        "OM:0502",
        "Creative Commons Zero 1.0 Universal - public domain dedication",
        "https://creativecommons.org/publicdomain/zero/1.0/",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    BSD_3 = (
        "OM:0513",
        "BSD 3-Clause License",
        "https://opensource.org/license/bsd-3-clause",
        LicenseUse.UNKNOWN,
        LicenseUse.UNKNOWN,
    )
    GPL_3_0 = (
        "OM:0503",
        "GNU General Public License v3.0",
        "https://www.gnu.org/licenses/gpl-3.0.html",
        LicenseUse.UNKNOWN,
        LicenseUse.UNKNOWN,
    )
    GPL_2_0 = (
        "OM:0507",
        "GNU General Public License v2.0",
        "https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html",
        LicenseUse.UNKNOWN,
        LicenseUse.UNKNOWN,
    )
    CC_BY_SA_4_0 = (
        "OM:0510",
        "Creative Commons Attribution-ShareAlike 4.0 International license",
        "https://creativecommons.org/licenses/by-sa/4.0/",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    CC_BY_SA_3_0 = (
        "OM:0509",
        "Creative Commons Attribution-ShareAlike 3.0 Unported license",
        "https://creativecommons.org/licenses/by-sa/3.0/",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    CC_BY_NC_4_0 = (
        "OM:0511",
        "Creative Commons Attribution-NonCommercial 4.0 International license",
        "https://creativecommons.org/licenses/by-nc/4.0/",
        LicenseUse.ALLOWED,
        LicenseUse.NOT_ALLOWED,
    )
    MIT = (
        "OM:0504",
        "MIT License - permissive free software license",
        "https://opensource.org/licenses/MIT",
        LicenseUse.UNKNOWN,
        LicenseUse.UNKNOWN,
    )
    ACADEMIC_FREE = (
        "OM:0505",
        "Free for academic use, restrictions may apply for commercial use",
        None,
        LicenseUse.ALLOWED,
        LicenseUse.UNKNOWN,
    )
    BIGG = (
        "OM:0512",
        "BiGG Models license - free for academic and non-profit use",
        "http://bigg.ucsd.edu/license#license",
        LicenseUse.ALLOWED,
        LicenseUse.UNKNOWN,
    )
    HPO = (
        "OM:0506",
        "HPO License - free for academic and non-profit use",
        "https://hpo.jax.org/app/license",
        LicenseUse.ALLOWED,
        LicenseUse.UNKNOWN,
    )
    KEGG_ACADEMIC = (
        "OM:0508",
        "KEGG license - free for academic use, commercial use requires a license",
        "https://www.kegg.jp/kegg/legal.html",
        LicenseUse.ALLOWED,
        LicenseUse.REQUIRES_PERMISSION,
    )
    PUBLIC = (
        "OM:0515",
        "Public domain - copyright waived in favor of the general public",
        "",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    CC_BY_2_0 = (
        "OM:0527",
        "Creative Commons Attribution 2.0 Generic license",
        "https://creativecommons.org/licenses/by/2.0/",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    CC_BY_2_5 = (
        "OM:0516",
        "Creative Commons Attribution 2.5 Generic license",
        "https://creativecommons.org/licenses/by/2.5/",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    CC_BY_ND_3_0 = (
        "OM:0517",
        "Creative Commons Attribution-NoDerivs 3.0 Unported license",
        "https://creativecommons.org/licenses/by-nd/3.0/",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    CC_BY_NC_2_0 = (
        "OM:0518",
        "Creative Commons Attribution-NonCommercial 2.0 Generic license",
        "https://creativecommons.org/licenses/by-nc/2.0/",
        LicenseUse.ALLOWED,
        LicenseUse.NOT_ALLOWED,
    )
    CC_BY_NC_3_0 = (
        "OM:0519",
        "Creative Commons Attribution-NonCommercial 3.0 Unported license",
        "https://creativecommons.org/licenses/by-nc/3.0/",
        LicenseUse.ALLOWED,
        LicenseUse.NOT_ALLOWED,
    )
    CC_BY_NC_SA_3_0 = (
        "OM:0520",
        "Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported license",
        "https://creativecommons.org/licenses/by-nc-sa/3.0/",
        LicenseUse.ALLOWED,
        LicenseUse.NOT_ALLOWED,
    )
    CC_BY_NC_SA_4_0 = (
        "OM:0521",
        "Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International license",
        "https://creativecommons.org/licenses/by-nc-sa/4.0/",
        LicenseUse.ALLOWED,
        LicenseUse.NOT_ALLOWED,
    )
    CC_BY_NC_ND_3_0 = (
        "OM:0522",
        "Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Unported license",
        "https://creativecommons.org/licenses/by-nc-nd/3.0/",
        LicenseUse.ALLOWED,
        LicenseUse.NOT_ALLOWED,
    )
    DSL = (
        "OM:0523",
        "Design Science License",
        "https://www.gnu.org/licenses/dsl.html",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    ARTISTIC_2_0 = (
        "OM:0524",
        "Artistic License 2.0",
        "https://opensource.org/license/artistic-2-0",
        LicenseUse.UNKNOWN,
        LicenseUse.UNKNOWN,
    )
    LGPL_3_0 = (
        "OM:0525",
        "GNU Lesser General Public License v3.0",
        "https://www.gnu.org/licenses/lgpl-3.0.html",
        LicenseUse.UNKNOWN,
        LicenseUse.UNKNOWN,
    )
    COMPOSITE = (
        "OM:0526",
        "Composite - the constituent databases carry their own licenses",
        None,
        LicenseUse.UNKNOWN,
        LicenseUse.UNKNOWN,
    )
    NONCOMMERCIAL = (
        "OM:0528",
        "Non-commercial use only",
        None,
        LicenseUse.ALLOWED,
        LicenseUse.NOT_ALLOWED,
    )
    NAR_OPEN_ACCESS = (
        "OM:0529",
        "Nucleic Acids Research Open Access - databases of the NAR database issue are freely usable",
        None,
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    MPPI = (
        "OM:0530",
        "MPPI license - free use, commercial use included",
        "http://mips.helmholtz-muenchen.de/proj/ppi/",
        LicenseUse.ALLOWED,
        LicenseUse.ALLOWED,
    )
    UNSPECIFIED = (
        "OM:0599",
        "License not specified or unclear",
        None,
        LicenseUse.UNKNOWN,
        LicenseUse.UNKNOWN,
    )


class UpdateCategoryCV(CvEnum):
    """Update frequency categories backed by OmniPath CV accessions.

    Describes how frequently a data source is expected to be updated.
    """

    parent_cv_term = ("OM:0420", "Update frequency term", "Describes how frequently a data source is expected to be updated.")

    # OmniPath update category terms (OM:0421-0429 range)
    REGULAR = (
        "OM:0421",
        "Regular scheduled updates (e.g., monthly, quarterly, annually)"
    )
    IRREGULAR = (
        "OM:0422",
        "Irregular or occasional updates with no fixed schedule"
    )
    STATIC = (
        "OM:0423",
        "Static resource with no planned future updates"
    )
