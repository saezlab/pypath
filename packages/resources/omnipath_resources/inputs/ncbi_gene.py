#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
NCBI Gene — the NCBI gene database (https://www.ncbi.nlm.nih.gov/gene). NCBI is
the authority for Entrez Gene identifiers. This module reads its bulk data files
under ``ftp.ncbi.nlm.nih.gov/gene/DATA/``.

``gene2ensembl`` provides, per taxon, the authoritative Entrez <-> Ensembl gene
(ENSG) and protein (ENSP) correspondence, covering every transcript — unlike the
UniProt idmapping, which cross-references only a single canonical ENSP per protein
and so misses the other-transcript ENSPs resources supply (e.g. STITCH).
"""

from __future__ import annotations

import collections

import pypath.resources.urls as urls
from pypath.share.downloads import dm

__all__ = ['gene2ensembl']


Gene2EnsemblRecord = collections.namedtuple(
    'Gene2EnsemblRecord',
    ['ncbi_tax_id', 'entrez', 'ensembl_gene', 'ensembl_protein'],
)


def gene2ensembl():
    """
    Yield ``Gene2EnsemblRecord`` from NCBI Gene's ``gene2ensembl.gz``.

    Ensembl ids are returned **versionless** (the ``.N`` suffix is stripped) so
    they match the versionless ids resources supply. Missing values (``-`` in the
    file) are returned as ``None``. The file is ~290 MB gzipped and cached by the
    download manager, so re-runs do not re-fetch.

    File columns: ``tax_id, GeneID, Ensembl_gene_identifier,
    RNA_nucleotide_accession.version, Ensembl_rna_identifier,
    protein_accession.version, Ensembl_protein_identifier``.
    """

    import gzip

    url = urls.urls['gene2ensembl']['url']
    path = dm.download(url)

    if not path:

        return

    def _strip(value):

        if not value or value == '-':

            return None

        return value.split('.', 1)[0]

    with gzip.open(path, 'rt') as fp:

        for line in fp:

            if line.startswith('#'):

                continue

            fields = line.rstrip('\n').split('\t')

            if len(fields) < 7:

                continue

            yield Gene2EnsemblRecord(
                ncbi_tax_id = int(fields[0]) if fields[0].isdigit() else 0,
                entrez = fields[1] if fields[1] != '-' else None,
                ensembl_gene = _strip(fields[2]),
                ensembl_protein = _strip(fields[6]),
            )
