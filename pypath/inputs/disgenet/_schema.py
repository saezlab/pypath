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

from pypath.inputs.disgenet._field import (
    DiseaseClasses,
    Float,
    Int,
    ProteinClass,
    Str,
    Tuple,
)


SCHEMA = {
    'entrez': (Str, 'geneid'),
    'genesymbol': (Str, 'gene_symbol'),
    'uniprot': (Str, 'uniprotid'),
    'dsi': (Float, '{entity_type}_dsi'),
    'dpi': (Float, '{entity_type}_dpi'),
    'pli': (Float, 'gene_pli'),
    'protein_class': (ProteinClass, ('protein_class', 'protein_class_name')),
    'classes': (DiseaseClasses, ('disease{i}_class', 'disease{i}_class_name')),
    'jaccard_genes': Float,
    'jaccard_variants': Float,
    'pvalue_jaccard_genes': Float,
    'pvalue_jaccard_variants': Float,
    'type': (Str, 'disease_type'),
    'id': (Str, '{entity_type}id{i}'),
    'name': (Str, 'disease{i}_name'),
    'semantic_type': (Str, 'disease_semantic_type'),
    'ngenes': (Int, 'disease{i}_ngenes'),
    'nvariants': (Int, 'disease{i}_nvariants'),
    'consequence_type': (Str, 'variant_consequence_type'),
    'score': Float,
    'ei': Float,
    'el': Str,
    'year_initial': Int,
    'year_final': Int,
    'source': Str,
}
