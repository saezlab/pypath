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

import enum


class Get(enum.StrEnum):

    SOURCE = 'source'
    PVALUE = 'pvalue'
    LIMIT = 'limit'
    DISEASE = 'disease'
    VARIANT = 'variant'
    GENE = 'gene'
    MIN_SCORE = 'min_score'
    MAX_SCORE = 'max_score'
    MIN_EI = 'min_ei'
    MAX_EI = 'max_ei'
    TYPE = 'type'
    DISEASE_CLASS = 'disease_class'
    MIN_DSI = 'min_dsi'
    MAX_DSI = 'max_dsi'
    MIN_DPI = 'min_dpi'
    MAX_DPI = 'max_dpi'
    MIN_PLI = 'min_pli'
    MAX_PLI = 'max_pli'
    MIN_YEAR = 'min_year'
    MAX_YEAR = 'max_year'
    OFFSET = 'offset'
    FORMAT = 'format'


class Querytype(enum.StrEnum):

    DDA = 'dda'
    GDA = 'gda'
    VDA = 'vda'


class By(enum.StrEnum):

    SOURCE = 'source'
    DISEASE = 'disease'
    VARIANT = 'variant'
    GENE = 'gene'


class Idtype(enum.StrEnum):

    UNIPROT = 'uniprot'
    ICD9CM = 'icd9cm'
    ICD10 = 'icd10'
    MESH = 'mesh'
    OMIM = 'omim'
    DO = 'do'
    EFO = 'efo'
    NCI = 'nci'
    HPO = 'hpo'
    MONDO = 'mondo'
    ORDO = 'ordo'


class Diseaseclass(enum.StrEnum):

    C01 = 'C01'
    C04 = 'C04'
    C05 = 'C05'
    C06 = 'C06'
    C07 = 'C07'
    C08 = 'C08'
    C09 = 'C09'
    C10 = 'C10'
    C11 = 'C11'
    C12 = 'C12'
    C13 = 'C13'
    C14 = 'C14'
    C15 = 'C15'
    C16 = 'C16'
    C17 = 'C17'
    C18 = 'C18'
    C19 = 'C19'
    C20 = 'C20'
    C21 = 'C21'
    C22 = 'C22'
    C23 = 'C23'
    C24 = 'C24'
    C25 = 'C25'
    C26 = 'C26'
    F01 = 'F01'
    F02 = 'F02'
    F03 = 'F03'


class Diseasetype(enum.StrEnum):

    DISEASE = 'disease'
    PHENOTYPE = 'phenotype'
    GROUP = 'group'


class Source(enum.StrEnum):

    CURATED = 'CURATED'
    BEFREE = 'BEFREE'
    ALL = 'ALL'
    CLINVAR = 'CLINVAR'
    GWASCAT = 'GWASCAT'
    GWASDB = 'GWASDB'
    UNIPROT = 'UNIPROT'
    INFERRED = 'INFERRED'
    ANIMAL_MODELS = 'ANIMAL_MODELS'
    CGI = 'CGI'
    CLINGEN = 'CLINGEN'
    CTD_HUMAN = 'CTD_human'
    CTD_MOUSE = 'CTD_mouse'
    CTD_RAT = 'CTD_rat'
    GENOMICS_ENGLAND = 'GENOMICS_ENGLAND'
    HPO = 'HPO'
    LHGDN = 'LHGDN'
    MGD = 'MGD'
    ORPHANET = 'ORPHANET'
    PSYGENET = 'PSYGENET'
    RGD = 'RGD'
