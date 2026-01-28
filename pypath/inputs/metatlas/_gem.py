#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""
Standard GEM TSV file parsing functions.

Standard GEMs may contain annotation files in TSV format:
- reactions.tsv - Reaction ID mappings
- metabolites.tsv - Metabolite ID mappings
- genes.tsv - Gene annotations

Note: Column names vary between GEMs. SysBioChalmers GEMs use names like
'rxnKEGGID', 'rxnBiGGID', while others use 'kegg.reaction', 'bigg.reaction'.
The functions return raw dictionaries to accommodate this variation.
"""

from __future__ import annotations

from collections.abc import Generator
import csv
import io

from ._common import _log
from ._git import git_raw_file, _parse_gem_index

__all__ = [
    'metatlas_gem_reactions',
    'metatlas_gem_metabolites',
    'metatlas_gem_genes',
    'metatlas_gem_tsv',
]


def _get_gem_info(gem: str) -> tuple[str, str] | None:
    """
    Retrieves git host and repo for a GEM from the index.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').

    Returns:
        Tuple of (host, repo) or None if not found.
    """

    index = _parse_gem_index()

    for host, repos in index.items():
        for repo in repos:
            if repo.split('/')[-1] == gem:
                return host, repo

    _log(f'GEM {gem} not found in index.')
    return None


def metatlas_gem_tsv(
        gem: str,
        file: str,
        ref: str | None = None,
) -> Generator[dict, None, None]:
    """
    Downloads and parses a TSV file from a standard GEM repository.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        file: Relative path to the TSV file (e.g., 'model/reactions.tsv').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        Dictionaries for each row in the TSV file.
        Column names vary between GEMs.
    """

    gem_info = _get_gem_info(gem)

    if gem_info is None:
        return

    host, repo = gem_info

    _log(f'Downloading {file} from {gem}.')

    content = git_raw_file(host, repo, ref, file)

    if content is None:
        return

    reader = csv.DictReader(io.StringIO(content), delimiter='\t')

    for row in reader:
        yield row


def metatlas_gem_reactions(
        gem: str,
        ref: str | None = None,
) -> Generator[dict, None, None]:
    """
    Parses reaction annotations from a standard GEM.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        Dictionaries with reaction ID mappings.
        Column names vary between GEMs (e.g., 'rxnKEGGID' vs 'kegg.reaction').
    """

    _log(f'Parsing reactions from {gem}.')

    yield from metatlas_gem_tsv(gem, 'model/reactions.tsv', ref)


def metatlas_gem_metabolites(
        gem: str,
        ref: str | None = None,
) -> Generator[dict, None, None]:
    """
    Parses metabolite annotations from a standard GEM.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        Dictionaries with metabolite ID mappings.
        Column names vary between GEMs.
    """

    _log(f'Parsing metabolites from {gem}.')

    yield from metatlas_gem_tsv(gem, 'model/metabolites.tsv', ref)


def metatlas_gem_genes(
        gem: str,
        ref: str | None = None,
) -> Generator[dict, None, None]:
    """
    Parses gene annotations from a standard GEM.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        Dictionaries with gene annotations.
        Column names vary between GEMs.
    """

    _log(f'Parsing genes from {gem}.')

    yield from metatlas_gem_tsv(gem, 'model/genes.tsv', ref)
