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

"""UniProt FTP ID mapping file access.

Downloads and parses the UniProt idmapping.dat.gz files -- the most
comprehensive source for UniProt cross-references. Available as a
complete dump (all organisms, ~18GB) or per-organism files (~100MB each).
"""

from __future__ import annotations

import gzip
import logging
import os
from collections import defaultdict

from pypath.share.downloads import dm

_log = logging.getLogger(__name__)

FTP_BASES = (
    'https://ftp.uniprot.org/pub/databases/uniprot',
    'https://ftp.ebi.ac.uk/pub/databases/uniprot',
    'https://ftp.expasy.org/databases/uniprot',
)

IDMAPPING_PATH = 'current_release/knowledgebase/idmapping'

# ID type names in the idmapping.dat file -> our canonical names
IDTYPE_MAP = {
    'UniProtKB-ID': 'uniprot_entry',
    'Gene_Name': 'genesymbol',
    'GeneID': 'entrez',
    'RefSeq': 'refseqp',
    'RefSeq_NT': 'refseqn',
    'GI': 'gi',
    'PDB': 'pdb',
    'GO': 'go',
    'UniRef100': 'uniref100',
    'UniRef90': 'uniref90',
    'UniRef50': 'uniref50',
    'UniParc': 'uniparc',
    'EMBL-CDS': 'embl',
    'EMBL': 'embl_id',
    'Ensembl': 'ensg',
    'Ensembl_TRS': 'enst',
    'Ensembl_PRO': 'ensp',
    'HGNC': 'hgnc',
    'KEGG': 'kegg',
    'NCBI_TaxID': '_taxid',  # special: not an ID mapping, it's the organism
    'ChEMBL': 'chembl',
    'DrugBank': 'drugbank',
    'Reactome': 'reactome',
    'STRING': 'string',
}


_CODES = {
    9606: 'HUMAN', 10090: 'MOUSE', 10116: 'RAT',
    559292: 'YEAST', 83333: 'ECOLI', 7227: 'DROME',
    7955: 'DANRE', 6239: 'CAEEL', 9031: 'CHICK',
    3702: 'ARATH', 44689: 'DICDI', 284812: 'SCHPO',
}


def organism_urls(ncbi_tax_id: int) -> list[str]:
    """Get URLs for per-organism idmapping files (primary + mirrors)."""

    code = _CODES.get(ncbi_tax_id)

    if not code:
        return []

    path = f'{IDMAPPING_PATH}/by_organism/{code}_{ncbi_tax_id}_idmapping.dat.gz'

    return [f'{base}/{path}' for base in FTP_BASES]


def idmapping_stream(
    ncbi_tax_id: int | None = None,
    id_types: set[str] | None = None,
):
    """Stream ID mapping records from a UniProt FTP file.

    Args:
        ncbi_tax_id: If given, download the per-organism file. If None,
            download the full idmapping.dat.gz (very large!).
        id_types: If given, only yield rows for these ID type names
            (as they appear in the file, e.g. 'Gene_Name', 'GeneID').
            Pass None to get all types.

    Yields:
        Tuples of (uniprot_ac, id_type_name, id_value).
    """

    if ncbi_tax_id:
        urls = organism_urls(ncbi_tax_id)
    else:
        urls = [f'{base}/{IDMAPPING_PATH}/idmapping.dat.gz' for base in FTP_BASES]

    if not urls:
        _log.warning('No FTP URL for taxid %s', ncbi_tax_id)
        return

    path = None

    for url in urls:
        _log.info('Trying UniProt FTP: %s', url)

        try:
            path = dm.download(url, connecttimeout=10)

            if path and os.path.getsize(path) > 0:
                break

            if path:
                _log.warning('Empty file from %s, trying next mirror', url)
                path = None

        except Exception as e:
            _log.warning('Failed %s: %s', url, e)

    if not path:
        _log.error('All UniProt FTP mirrors failed')
        return

    _log.info('Parsing %s', path)

    with gzip.open(path, 'rt') as f:

        for line in f:
            parts = line.rstrip('\n').split('\t')

            if len(parts) != 3:
                continue

            uniprot_ac, id_type_name, id_value = parts

            if id_types and id_type_name not in id_types:
                continue

            yield uniprot_ac, id_type_name, id_value


def idmapping(
    id_type: str,
    ncbi_tax_id: int = 9606,
) -> dict[str, set[str]]:
    """Load a mapping table from the FTP idmapping file.

    Args:
        id_type: The target ID type name as it appears in the file
            (e.g. 'Gene_Name', 'GeneID', 'Ensembl').
        ncbi_tax_id: Organism.

    Returns:
        Dict mapping UniProt AC -> set of target IDs.
    """

    data = defaultdict(set)

    for uniprot_ac, id_type_name, id_value in idmapping_stream(
        ncbi_tax_id=ncbi_tax_id,
        id_types={id_type},
    ):
        data[uniprot_ac].add(id_value)

    _log.info(
        'FTP idmapping: %d UniProt ACs with %s for organism %d',
        len(data), id_type, ncbi_tax_id,
    )

    return dict(data)


def all_id_types(ncbi_tax_id: int = 9606) -> set[str]:
    """Get all ID type names present in the file for an organism."""

    types = set()

    for _, id_type_name, _ in idmapping_stream(ncbi_tax_id=ncbi_tax_id):
        types.add(id_type_name)

    return types


def full_idmapping_urls() -> list[str]:
    """URLs for the complete idmapping.dat.gz (all organisms)."""
    return [
        f"{base}/{IDMAPPING_PATH}/idmapping.dat.gz"
        for base in FTP_BASES
    ]


def idmapping_full_stream(
    id_types: set[str] | None = None,
    include_taxids: bool = True,
):
    """Stream ALL records from the complete idmapping.dat.gz.

    This file is ~18GB compressed. Lines are streamed from gzip.

    Args:
        id_types: If given, only yield rows for these ID type names.
            NCBI_TaxID rows are always included if include_taxids is True.
        include_taxids: Always yield NCBI_TaxID rows (for organism assignment).

    Yields:
        Tuples of (uniprot_ac, id_type_name, id_value).
    """
    urls = full_idmapping_urls()

    path = None
    for url in urls:
        _log.info("Trying full idmapping: %s", url)
        try:
            path = dm.download(url, connecttimeout=10)
            if path and os.path.getsize(path) > 1000:
                break
            path = None
        except Exception as e:
            _log.warning("Failed %s: %s", url, e)

    if not path:
        _log.error("All mirrors failed for full idmapping.dat.gz")
        return

    _log.info("Streaming full idmapping from %s", path)

    filter_types = set(id_types) if id_types else None
    if filter_types and include_taxids:
        filter_types.add("NCBI_TaxID")

    count = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue

            uniprot_ac, id_type_name, id_value = parts

            if filter_types and id_type_name not in filter_types:
                continue

            yield uniprot_ac, id_type_name, id_value
            count += 1

            if count % 10_000_000 == 0:
                _log.info("Streamed %dM records", count // 1_000_000)

    _log.info("Finished streaming: %d total records", count)
