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

from __future__ import annotations

from typing import NamedTuple


class DiseaseClass(NamedTuple):
    id: str
    name: str


class ProteinClass(NamedTuple):
    id: str
    name: str


class Disease(NamedTuple):
    id: str
    name: str
    classes: tuple[DiseaseClass]
    ngenes: int
    nvariants: int
    type: str
    semantic_type: str


class Gene(NamedTuple):
    uniprot: str
    genesymbol: str
    entrez: str
    dsi: float
    dpi: float
    pli: float
    protein_class: ProteinClass


class Variant(NamedTuple):
    id: str
    genesymbol: str
    dsi: float
    dpi: float
    consequence_type: str


class DiseaseDiseaseAssociation(NamedTuple):
    disease1: Disease
    disease2: Disease
    jaccard_genes: float
    pvalue_jaccard_genes: float
    jaccard_variants: float
    pvalue_jaccard_variants: float
    source: str
    ngenes: int


class GeneDiseaseAssociation(NamedTuple):
    gene: Gene
    disease: Disease
    score: float
    ei: float
    el: str
    year_initial: int
    year_final: int
    source: str


class VariantDiseaseAssociation(NamedTuple):
    variant: Variant
    disease: Disease
    score: float
    ei: float
    year_initial: int
    year_final: int
    source: str


class VariantGeneMapping(NamedTuple):
    entrez: str
    genesymbol: str
    sources: tuple[str]


class IdType(NamedTuple):
    name: str
    code: str
    label: str


class DiseaseIdMapping(NamedTuple):
    name: str
    id_types: tuple[IdType]


class DisGeNetAnnotation(NamedTuple):
    disease: str
    type: str
    score: float
    dsi: float
    dpi: float
    nof_pmids: int
    nof_snps: int
    source: str
