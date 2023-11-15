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

import os

import pypath.share.common as common
import pypath_common._constants as _const
import pypath.internals.input_formats as input_formats

__all__ = ['uniprot', 'mirbase', 'basic']


uniprot = {
    ('embl', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'embl',
        ),
    ('genesymbol', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'genesymbol',
        ),
    ('genesymbol-syn', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'genesymbol-syn'
        ),
    ('entrez', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'entrez',
        ),
    ('hgnc', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'hgnc',
        ),
    # turned off this as the uploadlists also provides ENST IDs
    #('enst', 'uniprot'):
        #input_formats.UniprotMapping(
            #id_type_a = 'enst',
        #),
    ('refseqp', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'refseqp',
        ),
    ('uniprot-entry', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'uniprot-entry',
        ),
    ('protein-name', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'protein-name',
        ),
    ('protein-name-all', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'protein-name',
            swissprot = None,
        ),
    ('ec', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'ec',
        ),
}


mirbase = {
    ('mir-mat-name', 'mirbase'):
        input_formats.FileMapping(
            id_type_a = 'mir-mat-name',
            id_type_b = 'mirbase',
            input_ = 'mirbase.mirbase_mature',
            col_a = 1,
            col_b = 0,
            separator = None,
            header = 0,
        ),
    ('mir-name', 'mir-pre'):
        input_formats.FileMapping(
            id_type_a = 'mir-name',
            id_type_b = 'mir-pre',
            input_ = 'mirbase.mirbase_precursor',
            col_a = 1,
            col_b = 0,
            separator = None,
            header = 0,
        ),
    ('mir-pre', 'mirbase'):
        input_formats.FileMapping(
            id_type_a = 'mir-pre',
            id_type_b = 'mirbase',
            input_ = 'mirbase.mirbase_ids',
            col_a = 1,
            col_b = 0,
            separator = None,
            header = 0,
        ),
    ('mir-name', 'mirbase'):
        input_formats.FileMapping(
            id_type_a = 'mir-name',
            id_type_b = 'mirbase',
            input_ = 'mirbase.mirbase_precursor_to_mature',
            col_a = 0,
            col_b = 1,
            separator = None,
            header = 0,
        ),
}


ipi = {
    ('ipi', 'uniprot'):
        input_formats.FileMapping(
            id_type_a = 'ipi',
            id_type_b = 'uniprot',
            input_ = 'ipi._ipi_uniprot_pairs',
            col_a = 0,
            col_b = 1,
            separator = None,
            header = 0,
        )
}


basic = {
    ('uniprot-sec', 'uniprot-pri'):
        input_formats.FileMapping(
            id_type_a = 'uniprot-sec',
            id_type_b = 'uniprot-pri',
            input_ = 'uniprot.get_uniprot_sec',
            col_a = 0,
            col_b = 1,
            separator = None,
            header = 0,
            ncbi_tax_id = _const.NOT_ORGANISM_SPECIFIC,
        ),
    ('genesymbol', 'trembl'):
        input_formats.UniprotMapping(
            id_type_a = 'genesymbol',
            id_type_b = 'trembl',
            swissprot = 'false',
        ),
    ('genesymbol', 'swissprot'):
        input_formats.UniprotMapping(
            id_type_a = 'genesymbol',
            id_type_b = 'swissprot',
        ),
    ('genesymbol-syn', 'swissprot'):
        input_formats.UniprotMapping(
            id_type_a = 'genesymbol-syn',
            id_type_b = 'swissprot',
        ),
    ('genesymbol', 'uniprot'):
        input_formats.UniprotMapping(
            id_type_a = 'genesymbol',
            swissprot = None,
        ),
}
