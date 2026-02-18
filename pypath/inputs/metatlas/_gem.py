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
Standard GEM file parsing functions.

Standard GEMs may contain annotation files in TSV format:
- reactions.tsv - Reaction ID mappings
- metabolites.tsv - Metabolite ID mappings
- genes.tsv - Gene annotations

Note: Column names vary between GEMs. SysBioChalmers GEMs use names like
'rxnKEGGID', 'rxnBiGGID', while others use 'kegg.reaction', 'bigg.reaction'.
The TSV functions return raw dictionaries to accommodate this variation.

The YAML model file (model/{GEM}.yml) is present in all standard-GEM repos
and contains stoichiometry, flux bounds, gene rules, compartments, and
cross-references. The YAML functions return typed GemReaction and
GemMetabolite named tuples.
"""

from __future__ import annotations

from collections.abc import Generator
import csv
import io
import re

import yaml

from ._common import _log
from ._git import git_raw_file, gem_file_path, _parse_gem_index
from ._records import GemReaction, GemMetabolite

__all__ = [
    'metatlas_gem_reactions',
    'metatlas_gem_metabolites',
    'metatlas_gem_genes',
    'metatlas_gem_tsv',
    'metatlas_gem_yaml',
    'metatlas_gem_yaml_reactions',
    'metatlas_gem_yaml_metabolites',
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
        _log(f'File `{file}` not available in {gem}; '
             'not all GEMs provide TSV annotations.')
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


_RE_BROKEN_QUOTES = re.compile(
    r'^(\s*- \w+: )'        # YAML omap key prefix, e.g. "      - name: "
    r'("(?:[^"\\]|\\.)*")'  # a double-quoted scalar, e.g. '"5"'
    r'(\S.*)$',             # trailing text that shouldn't be there
    re.MULTILINE,
)


def _fix_yaml_quoting(content: str) -> str:
    """
    Fix broken double-quote usage in GEM YAML files.

    Some GEM YAML files contain unescaped double quotes in string values,
    e.g. ``- name: "5"-deoxyadenosine...`` where ``"5"`` is parsed by YAML
    as a complete quoted scalar.  This wraps such values in single quotes.
    """

    def _fix_match(m):
        prefix = m.group(1)
        value = m.group(2) + m.group(3)
        escaped = value.replace("'", "''")
        return f"{prefix}'{escaped}'"

    return _RE_BROKEN_QUOTES.sub(_fix_match, content)


def metatlas_gem_yaml(
        gem: str,
        ref: str | None = None,
) -> dict | None:
    """
    Downloads and parses the YAML model file from a standard GEM repository.

    The YAML file (model/{gem}.yml) is present in all standard-GEM repos and
    contains the full model: reactions with stoichiometry and gene rules,
    metabolites with compartment annotations, and compartment definitions.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Returns:
        Dictionary with keys 'reactions', 'metabolites', 'compartments',
        'genes', and 'metaData'; or None if download failed.
    """

    gem_info = _get_gem_info(gem)

    if gem_info is None:
        return None

    host, repo = gem_info
    path = gem_file_path(gem, 'yml')

    _log(f'Downloading YAML model from {gem}.')

    content = git_raw_file(host, repo, ref, path)

    if content is None:
        _log(f'YAML model file not found for {gem} at {path}.')
        return None

    _log(f'Parsing YAML model for {gem}.')

    content = _fix_yaml_quoting(content)
    data = yaml.safe_load(content)

    # yaml.safe_load with !!omap returns list of (key, value) tuples
    if isinstance(data, list):
        data = dict(data)

    return data


def metatlas_gem_yaml_reactions(
        gem: str,
        ref: str | None = None,
) -> Generator[GemReaction, None, None]:
    """
    Parses reactions from the YAML model of a standard GEM.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        GemReaction named tuples with stoichiometry, bounds, and gene rules.
    """

    data = metatlas_gem_yaml(gem, ref)

    if data is None:
        return

    reactions = data.get('reactions', [])

    _log(f'Processing {len(reactions)} reactions from {gem} YAML.')

    for rxn in reactions:

        if isinstance(rxn, list):
            rxn = dict(rxn)

        mets = rxn.get('metabolites', {})

        if isinstance(mets, list):
            mets = dict(mets)

        subsystem = rxn.get('subsystem', '')

        if isinstance(subsystem, list):
            subsystem = '; '.join(str(s) for s in subsystem)

        eccodes = rxn.get('eccodes', '')

        if isinstance(eccodes, list):
            eccodes = '; '.join(str(e) for e in eccodes)

        yield GemReaction(
            id=rxn.get('id', ''),
            name=rxn.get('name', ''),
            metabolites=mets,
            lower_bound=float(rxn.get('lower_bound', 0)),
            upper_bound=float(rxn.get('upper_bound', 0)),
            gene_reaction_rule=rxn.get('gene_reaction_rule', ''),
            subsystem=subsystem,
            eccodes=eccodes,
        )


def metatlas_gem_yaml_metabolites(
        gem: str,
        ref: str | None = None,
) -> Generator[GemMetabolite, None, None]:
    """
    Parses metabolites from the YAML model of a standard GEM.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit).
            If None, uses the repository's default branch.

    Yields:
        GemMetabolite named tuples with compartment, formula, and charge.
    """

    data = metatlas_gem_yaml(gem, ref)

    if data is None:
        return

    metabolites = data.get('metabolites', [])

    _log(f'Processing {len(metabolites)} metabolites from {gem} YAML.')

    for met in metabolites:

        if isinstance(met, list):
            met = dict(met)

        charge = met.get('charge')

        if charge is not None:
            charge = int(charge)

        yield GemMetabolite(
            id=met.get('id', ''),
            name=met.get('name', ''),
            compartment=met.get('compartment', ''),
            formula=met.get('formula', ''),
            charge=charge,
        )
