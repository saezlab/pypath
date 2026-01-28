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

Standard GEMs contain annotation files in TSV format:
- reactions.tsv - Reaction annotations
- metabolites.tsv - Metabolite annotations
- genes.tsv - Gene annotations
"""

from __future__ import annotations

from collections.abc import Generator
import csv
import io

from ._common import _log
from ._git import git_raw_file, _parse_gem_index
from ._records import MetatlasReaction, MetatlasMetabolite, MetatlasGene

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
            # GEM name is the last part of the repo path
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
        ref: Git reference (branch, tag, or commit). Defaults to 'main'.

    Yields:
        Dictionaries for each row in the TSV file.
    """

    gem_info = _get_gem_info(gem)

    if gem_info is None:
        return

    host, repo = gem_info
    ref = ref or 'main'

    _log(f'Downloading {file} from {gem}@{ref}.')

    content = git_raw_file(host, repo, ref, file)

    if content is None:
        return

    reader = csv.DictReader(io.StringIO(content), delimiter='\t')

    for row in reader:
        yield row


def metatlas_gem_reactions(
        gem: str,
        ref: str | None = None,
) -> Generator[MetatlasReaction, None, None]:
    """
    Parses reaction ID mappings from a standard GEM.

    The reactions.tsv file contains cross-references to external databases
    like KEGG, BiGG, MetaNetX, and Rhea.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit). Defaults to 'main'.

    Yields:
        MetatlasReaction named tuples with ID mappings.
    """

    _log(f'Parsing reactions from {gem}.')

    for row in metatlas_gem_tsv(gem, 'model/reactions.tsv', ref):
        yield MetatlasReaction(
            id=row.get('rxns', ''),
            kegg=row.get('rxnKEGGID', ''),
            bigg=row.get('rxnBiGGID', ''),
            metanetx=row.get('rxnMetaNetXID', ''),
            rhea=row.get('rxnRheaID', ''),
            rhea_master=row.get('rxnRheaMasterID', ''),
            reactome=row.get('rxnREACTOMEID', ''),
            spontaneous=row.get('spontaneous', '') == '1',
        )


def metatlas_gem_metabolites(
        gem: str,
        ref: str | None = None,
) -> Generator[MetatlasMetabolite, None, None]:
    """
    Parses metabolite ID mappings from a standard GEM.

    The metabolites.tsv file contains cross-references to external databases
    like KEGG, HMDB, ChEBI, PubChem, and LIPID MAPS.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit). Defaults to 'main'.

    Yields:
        MetatlasMetabolite named tuples with ID mappings.
    """

    _log(f'Parsing metabolites from {gem}.')

    for row in metatlas_gem_tsv(gem, 'model/metabolites.tsv', ref):
        yield MetatlasMetabolite(
            id=row.get('mets', ''),
            id_no_compartment=row.get('metsNoComp', ''),
            bigg=row.get('metBiGGID', ''),
            kegg=row.get('metKEGGID', ''),
            hmdb=row.get('metHMDBID', ''),
            chebi=row.get('metChEBIID', ''),
            pubchem=row.get('metPubChemID', ''),
            lipidmaps=row.get('metLipidMapsID', ''),
            metanetx=row.get('metMetaNetXID', ''),
        )


def metatlas_gem_genes(
        gem: str,
        ref: str | None = None,
) -> Generator[MetatlasGene, None, None]:
    """
    Parses gene annotations from a standard GEM.

    The genes.tsv file contains gene identifiers and their mappings
    to Ensembl, UniProt, Entrez, and other databases.

    Args:
        gem: Name of the GEM (e.g., 'Human-GEM').
        ref: Git reference (branch, tag, or commit). Defaults to 'main'.

    Yields:
        MetatlasGene named tuples with gene annotations.
    """

    _log(f'Parsing genes from {gem}.')

    for row in metatlas_gem_tsv(gem, 'model/genes.tsv', ref):
        yield MetatlasGene(
            id=row.get('genes', ''),
            ensembl_transcript=row.get('geneENSTID', ''),
            ensembl_protein=row.get('geneENSPID', ''),
            uniprot=row.get('geneUniProtID', ''),
            symbol=row.get('geneSymbols', ''),
            entrez=row.get('geneEntrezID', ''),
            name=row.get('geneNames', ''),
            aliases=row.get('geneAliases', ''),
            compartments=row.get('compartments', ''),
        )
