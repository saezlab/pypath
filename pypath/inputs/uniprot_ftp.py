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
from collections import defaultdict

from pypath.share.downloads import dm

_log = logging.getLogger(__name__)

FTP_BASE = (
    'https://ftp.uniprot.org/pub/databases/uniprot/'
    'current_release/knowledgebase/idmapping'
)

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


def organism_url(ncbi_tax_id: int) -> str | None:
    """Get the URL for a per-organism idmapping file."""

    # The filename pattern is {CODE}_{TAXID}_idmapping.dat.gz
    _CODES = {
        9606: 'HUMAN', 10090: 'MOUSE', 10116: 'RAT',
        559292: 'YEAST', 83333: 'ECOLI', 7227: 'DROME',
        7955: 'DANRE', 6239: 'CAEEL', 9031: 'CHICK',
        3702: 'ARATH', 44689: 'DICDI', 284812: 'SCHPO',
    }

    code = _CODES.get(ncbi_tax_id)

    if not code:
        return None

    return f'{FTP_BASE}/by_organism/{code}_{ncbi_tax_id}_idmapping.dat.gz'


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
        url = organism_url(ncbi_tax_id)

        if not url:
            _log.warning(
                'No per-organism FTP file for taxid %d', ncbi_tax_id,
            )
            return
    else:
        url = f'{FTP_BASE}/idmapping.dat.gz'

    _log.info('Downloading UniProt FTP idmapping from %s', url)

    path = dm.download(url)

    if not path:
        _log.error('Failed to download %s', url)
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
