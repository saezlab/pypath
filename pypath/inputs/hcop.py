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

"""HGNC Comparison of Orthology Predictions (HCOP).

Meta-aggregator of 14 orthology sources: Ensembl Compara, HGNC, HomoloGene,
Inparanoid, OMA, OrthoDB, PANTHER, PhylomeDB, TreeFam, and others.

Data from: http://ftp.ebi.ac.uk/pub/databases/genenames/hcop/
"""

from __future__ import annotations

import gzip
import logging
import collections

from pypath.share.downloads import dm

_log = logging.getLogger(__name__)

HCOP_BASE = 'http://ftp.ebi.ac.uk/pub/databases/genenames/hcop'

# Species name mapping for per-species files
_SPECIES_NAMES = {
    10090: 'mouse',
    10116: 'rat',
    7955: 'zebrafish',
    9031: 'chicken',
    9615: 'dog',
    9913: 'cattle',
    9823: 'pig',
    9796: 'horse',
    9685: 'cat',
    9598: 'chimpanzee',
    9544: 'macaque',
    7227: 'fruitfly',
    6239: 'c.elegans',
    4932: 'yeast',
    8364: 'frog',
    28377: 'anole_lizard',
    9258: 'platypus',
    13616: 'opossum',
    9986: 'rabbit',
    9940: 'sheep',
}

HcopRecord = collections.namedtuple(
    'HcopRecord',
    (
        'ortholog_taxid',
        'human_entrez',
        'human_ensembl',
        'hgnc_id',
        'human_name',
        'human_symbol',
        'human_chr',
        'ortholog_entrez',
        'ortholog_ensembl',
        'ortholog_db_id',
        'ortholog_name',
        'ortholog_symbol',
        'ortholog_chr',
        'support',
        'n_sources',
    ),
)


def hcop_orthologs(
    target_organism: int | None = None,
    min_sources: int = 1,
) -> list[HcopRecord]:
    """Download and parse HCOP orthology data.

    Args:
        target_organism: NCBI Taxonomy ID of the target organism.
            If None, downloads the full file (all organisms).
        min_sources: Minimum number of supporting databases.

    Returns:
        List of HcopRecord named tuples.
    """

    if target_organism and target_organism in _SPECIES_NAMES:
        species = _SPECIES_NAMES[target_organism]
        url = f'{HCOP_BASE}/human_{species}_hcop_fifteen_column.txt.gz'
        has_species_col = False
    else:
        url = f'{HCOP_BASE}/human_all_hcop_sixteen_column.txt.gz'
        has_species_col = True

    _log.info('Downloading HCOP from %s', url)
    path = dm.download(url)

    if not path:
        _log.error('Failed to download HCOP')
        return []

    records = []

    with gzip.open(path, 'rt') as f:
        header = next(f)  # skip header

        for line in f:
            parts = line.rstrip('\n').split('\t')

            if has_species_col:
                if len(parts) < 16:
                    continue
                taxid = int(parts[0]) if parts[0].isdigit() else 0
                # Filter by target organism if specified
                if target_organism and taxid != target_organism:
                    continue
                cols = parts[1:]  # skip species column for uniform indexing
                support = parts[15]
            else:
                if len(parts) < 15:
                    continue
                taxid = target_organism or 0
                cols = parts
                support = parts[14]

            # Parse support: comma-separated list of databases
            sources = [s.strip() for s in support.split(',') if s.strip()]
            n = len(sources)

            if n < min_sources:
                continue

            records.append(HcopRecord(
                ortholog_taxid=taxid,
                human_entrez=cols[0] if cols[0] != '-' else '',
                human_ensembl=cols[1] if cols[1] != '-' else '',
                hgnc_id=cols[2] if cols[2] != '-' else '',
                human_name=cols[3] if cols[3] != '-' else '',
                human_symbol=cols[4] if cols[4] != '-' else '',
                human_chr=cols[5] if cols[5] != '-' else '',
                # skip assert_ids (cols[6])
                ortholog_entrez=cols[7] if cols[7] != '-' else '',
                ortholog_ensembl=cols[8] if cols[8] != '-' else '',
                ortholog_db_id=cols[9] if cols[9] != '-' else '',
                ortholog_name=cols[10] if cols[10] != '-' else '',
                ortholog_symbol=cols[11] if cols[11] != '-' else '',
                ortholog_chr=cols[12] if cols[12] != '-' else '',
                support=support,
                n_sources=n,
            ))

    _log.info('HCOP: %d ortholog pairs for organism %s', len(records), target_organism or 'all')
    return records


def hcop_dict(
    target_organism: int = 10090,
    id_type: str = 'genesymbol',
    min_sources: int = 2,
) -> dict[str, set[str]]:
    """Get a human -> target organism orthology dict.

    Args:
        target_organism: NCBI Taxonomy ID.
        id_type: Identifier type to use. Options:
            'genesymbol', 'entrez', 'ensembl', 'hgnc'
        min_sources: Minimum supporting databases.

    Returns:
        Dict mapping human IDs to sets of ortholog IDs.
    """

    _HUMAN_FIELDS = {
        'genesymbol': 'human_symbol',
        'entrez': 'human_entrez',
        'ensembl': 'human_ensembl',
        'ensg': 'human_ensembl',
        'hgnc': 'hgnc_id',
    }
    _ORTHOLOG_FIELDS = {
        'genesymbol': 'ortholog_symbol',
        'entrez': 'ortholog_entrez',
        'ensembl': 'ortholog_ensembl',
        'ensg': 'ortholog_ensembl',
    }

    human_field = _HUMAN_FIELDS.get(id_type, 'human_symbol')
    ortholog_field = _ORTHOLOG_FIELDS.get(id_type, 'ortholog_symbol')

    records = hcop_orthologs(target_organism, min_sources=min_sources)

    result = collections.defaultdict(set)
    for rec in records:
        human_id = getattr(rec, human_field)
        ortholog_id = getattr(rec, ortholog_field)
        if human_id and ortholog_id:
            result[human_id].add(ortholog_id)

    return dict(result)
