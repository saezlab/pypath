#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
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

from __future__ import annotations

from typing import NamedTuple

__all__ = [
    'UniProtRecord',
]


class UniProtRecord(NamedTuple):
    """UniProt protein record."""
    accession: str | None = None
    entry_name: str | None = None
    protein_name: str | None = None
    length: int | None = None
    mass: int | None = None
    sequence: str | None = None
    gene_primary: str | None = None
    gene_synonym: str | None = None
    organism_id: int | None = None
    cc_disease: str | None = None
    ft_mutagen: str | None = None
    cc_subcellular_location: str | None = None
    cc_ptm: str | None = None
    lit_pubmed_id: str | None = None
    cc_function: str | None = None
    xref_ensembl: str | None = None
    xref_kegg: str | None = None
    cc_pathway: str | None = None
    cc_activity_regulation: str | None = None
    keyword: list[str] | None = None
    ec: str | None = None
    go: list[str] | None = None
    ft_transmem: str | None = None
    protein_families: str | None = None
    xref_refseq: str | None = None
    xref_alphafolddb: str | None = None
    xref_pdb: str | None = None
    xref_chembl: str | None = None
    xref_phosphositeplus: str | None = None
    xref_signor: str | None = None
    xref_pathwaycommons: str | None = None
    xref_intact: str | None = None
    xref_biogrid: str | None = None
    xref_complexportal: str | None = None
