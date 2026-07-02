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

"""Alliance of Genome Resources orthology data.

Curated orthology assertions from 9 model organism databases,
validated via DIOPT benchmarking. Data aggregates predictions from
multiple algorithms (Ensembl Compara, InParanoid, OMA, PANTHER, etc.)
and applies stringent/moderate confidence filters.

Downloads the combined TSV from the FMS (File Management System) API:
https://fms.alliancegenome.org/
"""

from __future__ import annotations

import gzip
import logging
import collections

from pypath.share.downloads import dm

_log = logging.getLogger(__name__)

FMS_COMBINED_URL = (
    'https://fms.alliancegenome.org/download/ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz'
)

AllianceOrtholog = collections.namedtuple(
    'AllianceOrtholog',
    (
        'gene1_id',
        'gene1_symbol',
        'gene1_taxid',
        'gene1_species',
        'gene2_id',
        'gene2_symbol',
        'gene2_taxid',
        'gene2_species',
        'algorithms',
        'n_algorithms',
        'out_of_algorithms',
        'is_best',
        'is_best_reverse',
    ),
)


def alliance_orthologs() -> list[AllianceOrtholog]:
    """Download and parse Alliance of Genome Resources orthology data.

    Downloads the combined TSV file containing all pairwise orthology
    assertions across the 9 Alliance species.

    Returns:
        List of AllianceOrtholog named tuples.
    """

    _log.info('Downloading Alliance orthology from %s', FMS_COMBINED_URL)

    path = dm.download(FMS_COMBINED_URL)

    if not path:
        _log.error('Failed to download Alliance orthology data')
        return []

    records = []

    with gzip.open(path, 'rt') as f:

        for line in f:

            if line.startswith('#') or line.startswith('Gene1ID'):
                continue

            parts = line.rstrip('\n').split('\t')

            if len(parts) < 13:
                continue

            taxid1_raw = parts[2]  # e.g. "NCBITaxon:9606"
            taxid2_raw = parts[6]
            taxid1 = int(taxid1_raw.split(':')[1]) if ':' in taxid1_raw else 0
            taxid2 = int(taxid2_raw.split(':')[1]) if ':' in taxid2_raw else 0

            n_algo = int(parts[9]) if parts[9].isdigit() else 0
            out_of = int(parts[10]) if parts[10].isdigit() else 0

            records.append(AllianceOrtholog(
                gene1_id=parts[0],
                gene1_symbol=parts[1],
                gene1_taxid=taxid1,
                gene1_species=parts[3],
                gene2_id=parts[4],
                gene2_symbol=parts[5],
                gene2_taxid=taxid2,
                gene2_species=parts[7],
                algorithms=parts[8],
                n_algorithms=n_algo,
                out_of_algorithms=out_of,
                is_best=parts[11] == 'Yes',
                is_best_reverse=parts[12] == 'Yes',
            ))

    _log.info('Alliance: parsed %d ortholog records', len(records))
    return records


def alliance_dict(
    source: int = 9606,
    target: int = 10090,
    id_type: str = 'genesymbol',
    best_only: bool = False,
    min_algorithms: int = 1,
) -> dict[str, set[str]]:
    """Get a pairwise orthology dict from Alliance data.

    Args:
        source: Source organism NCBI Taxonomy ID.
        target: Target organism NCBI Taxonomy ID.
        id_type: ``'genesymbol'`` for gene symbols, ``'dbid'`` for
            database-specific IDs (e.g. ``HGNC:5``, ``MGI:2152878``).
        best_only: If True, only include best-score pairs.
        min_algorithms: Minimum number of matching algorithms.

    Returns:
        Dict mapping source gene IDs to sets of target gene IDs.
    """

    records = alliance_orthologs()

    if not records:
        return {}

    data = collections.defaultdict(set)

    for rec in records:

        # Match source -> target direction
        if rec.gene1_taxid == source and rec.gene2_taxid == target:
            src_sym = rec.gene1_symbol
            tgt_sym = rec.gene2_symbol
            src_id = rec.gene1_id
            tgt_id = rec.gene2_id
        elif rec.gene2_taxid == source and rec.gene1_taxid == target:
            src_sym = rec.gene2_symbol
            tgt_sym = rec.gene1_symbol
            src_id = rec.gene2_id
            tgt_id = rec.gene1_id
        else:
            continue

        if best_only and not rec.is_best:
            continue

        if rec.n_algorithms < min_algorithms:
            continue

        if id_type == 'genesymbol':
            src = src_sym
            tgt = tgt_sym
        else:
            src = src_id
            tgt = tgt_id

        if src and tgt:
            data[src].add(tgt)

    _log.info(
        'Alliance: %d source genes with orthologs (%d -> %d)',
        len(data), source, target,
    )
    return dict(data)
