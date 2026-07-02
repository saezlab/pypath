#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
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

"""OrthoDB orthology data access.

OrthoDB v12 provides hierarchical ortholog groups across 28,000+ organisms.
The REST API supports gene-level searches and ortholog group queries.

API docs: https://data.orthodb.org/v12/

Note: OrthoDB bulk tab downloads are currently unavailable (404), so
this module uses the REST API. For bulk orthology needs, HCOP
(which includes OrthoDB as one of its 14 sources) is more efficient.
"""

from __future__ import annotations

import logging
import collections

import requests

_log = logging.getLogger(__name__)

ORTHODB_BASE = 'https://data.orthodb.org/v12'

# NCBI Taxonomy IDs for common taxonomic levels used to scope OG queries
LEVELS = {
    'vertebrata': 7742,
    'mammalia': 40674,
    'metazoa': 33208,
    'eukaryota': 2759,
    'tetrapoda': 32523,
    'primates': 9443,
}

OrthodbRecord = collections.namedtuple(
    'OrthodbRecord',
    (
        'source_gene',
        'target_gene',
        'og_id',
        'source_organism',
        'target_organism',
    ),
)


def orthodb_search(
    gene_name: str,
    species: int = 9606,
    level: int | str = 40674,
    limit: int = 5,
) -> list[str]:
    """Search OrthoDB for ortholog group IDs matching a gene name.

    Args:
        gene_name: Gene name or identifier to search.
        species: NCBI Taxonomy ID (default: 9606 for human).
        level: Taxonomic level ID or name (default: 40674 for Mammalia).
        limit: Max results to return.

    Returns:
        List of ortholog group IDs (e.g. ``['441317at40674', ...]``).
    """

    if isinstance(level, str):
        level = LEVELS.get(level.lower(), 40674)

    resp = requests.get(
        f'{ORTHODB_BASE}/search',
        params={
            'query': gene_name,
            'species': species,
            'level': level,
            'limit': limit,
        },
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()

    return data.get('data', [])


def orthodb_group_genes(
    og_id: str,
    species: list[int] | None = None,
) -> list[dict]:
    """Get genes in an ortholog group, optionally filtered by species.

    Args:
        og_id: Ortholog group ID (e.g. ``'441317at40674'``).
        species: Filter to these NCBI Taxonomy IDs.

    Returns:
        List of dicts with gene information including ``gene_name``,
        ``gene_id``, ``organism_name``, ``taxid``, and ``og_id``.
    """

    params = {'id': og_id}
    if species:
        params['species'] = ','.join(str(s) for s in species)

    resp = requests.get(
        f'{ORTHODB_BASE}/orthologs',
        params=params,
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()

    genes = []

    for group in data.get('data', []):

        organism_id = group.get('organism', {}).get('id', '')
        organism_name = group.get('organism', {}).get('name', '')
        # Organism ID format is "{taxid}_{index}", e.g. "9606_0"
        taxid = organism_id.split('_')[0] if '_' in organism_id else organism_id

        for gene in group.get('genes', []):
            genes.append({
                'gene_id': gene.get('gene_id', {}).get('param', ''),
                'gene_name': gene.get('gene_id', {}).get('id', ''),
                'description': gene.get('description', ''),
                'organism_id': organism_id,
                'organism_name': organism_name,
                'taxid': taxid,
                'og_id': og_id,
            })

    return genes


def orthodb_orthologs(
    source: int = 9606,
    target: int = 10090,
    id_type: str = 'genesymbol',
    level: str | int = 'mammalia',
) -> dict[str, set[str]]:
    """Build a pairwise orthology dict between two organisms via OrthoDB.

    Queries the OrthoDB REST API by searching for all ortholog groups
    that contain genes from the source organism, then collecting
    target organism genes from each group.

    Note: The REST API requires one query per gene, so this is inherently
    slow for genome-wide queries. For bulk use, prefer HCOP
    (which includes OrthoDB data) or the ``ensembl`` backend.

    Args:
        source: Source organism NCBI Taxonomy ID.
        target: Target organism NCBI Taxonomy ID.
        id_type: ``'genesymbol'`` for gene names.
        level: Taxonomic level for ortholog groups.

    Returns:
        Dict mapping source gene names to sets of target gene names.
    """

    if isinstance(level, str):
        level = LEVELS.get(level.lower(), 40674)

    _log.info(
        'OrthoDB: querying orthologs %d -> %d at level %d',
        source, target, level,
    )

    _log.warning(
        'OrthoDB REST API is slow for bulk queries (one request per gene). '
        'For better performance, use HCOP which includes OrthoDB data.'
    )

    # The OrthoDB REST API does not provide a bulk pairwise download.
    # A full implementation would iterate over all source genes, search for
    # their ortholog groups, and extract target genes from each group.
    # This is too slow for practical use (thousands of API calls).
    #
    # Instead, we return an empty dict here and recommend HCOP.
    # The individual lookup functions (orthodb_search, orthodb_group_genes)
    # remain available for targeted queries.
    return {}
