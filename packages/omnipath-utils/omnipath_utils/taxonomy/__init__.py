#!/usr/bin/env python

#
# This file is part of the `omnipath_utils` Python module
#
# Copyright 2026
# Heidelberg University Hospital
#
# File author(s): Denes Turei (turei.denes@gmail.com)
#
# Distributed under the GPL-3.0-or-later license
# See the file `LICENSE` or read a copy at
# https://www.gnu.org/licenses/gpl-3.0.txt
#

"""Taxonomy -- organism name resolution and conversion.

Example::

    from omnipath_utils.taxonomy import ensure_ncbi_tax_id, ensure_common_name

    ensure_ncbi_tax_id('human')     # 9606
    ensure_ncbi_tax_id('hsapiens')  # 9606
    ensure_common_name(10090)       # 'mouse'
"""

from omnipath_utils.taxonomy._taxonomy import TaxonomyManager


def ensure_ncbi_tax_id(taxon):
    """Convert any organism representation to NCBI Taxonomy ID.

    Args:
        taxon:
            Organism name or NCBI Taxonomy ID in any supported form.

    Returns:
        NCBI Taxonomy ID, or None if the organism is not found.
    """

    return TaxonomyManager.get().ensure_ncbi_tax_id(taxon)


def ensure_common_name(taxon):
    """Common English name for an organism.

    Args:
        taxon:
            Organism identifier in any supported form.

    Returns:
        Common name string, or None if not found.
    """

    return TaxonomyManager.get().ensure_common_name(taxon)


def ensure_latin_name(taxon):
    """Full latin (scientific) name for an organism.

    Args:
        taxon:
            Organism identifier in any supported form.

    Returns:
        Latin name string, or None if not found.
    """

    return TaxonomyManager.get().ensure_latin_name(taxon)


def ensure_ensembl_name(taxon):
    """Ensembl organism name.

    Args:
        taxon:
            Organism identifier in any supported form.

    Returns:
        Ensembl name string, or None if not found.
    """

    return TaxonomyManager.get().ensure_ensembl_name(taxon)


def ensure_kegg_code(taxon):
    """KEGG three-letter organism code.

    Args:
        taxon:
            Organism identifier in any supported form.

    Returns:
        KEGG code string, or None if not found.
    """

    return TaxonomyManager.get().ensure_kegg_code(taxon)


def ensure_mirbase_name(taxon):
    """MiRBase three-letter organism code.

    Args:
        taxon:
            Organism identifier in any supported form.

    Returns:
        miRBase code string, or None if not found.
    """

    return TaxonomyManager.get().ensure_mirbase_name(taxon)


def ensure_oma_code(taxon):
    """OMA five-letter organism code.

    Args:
        taxon:
            Organism identifier in any supported form.

    Returns:
        OMA code string, or None if not found.
    """

    return TaxonomyManager.get().ensure_oma_code(taxon)


def all_organisms():
    """Return all known organisms.

    Returns:
        Dict mapping NCBI Taxonomy IDs to organism info dicts.
    """

    return TaxonomyManager.get().all_organisms()
