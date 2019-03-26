#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import os

import pypath.common as common
import pypath.input_formats as input_formats

__all__ = ['misc', 'uniprot', 'mirbase', 'basic', 'uniprot_old']


misc = {
    ('uniprot-sec', 'uniprot-pri'):
        input_formats.FileMapping(
            id_type_a = 'uniprot-sec',
            id_type_b = 'uniprot-pri',
            entity_type = 'protein',
            fname = os.path.join(common.ROOT, 'data', 'sec_ac.txt'),
            col_a = 0,
            col_b = 1,
            separator = None,
            header = 0,
        ),
    ('trembl', 'genesymbol'):
        input_formats.FileMapping(
            id_type_a = 'trembl',
            id_type_b = 'genesymbol',
            entity_type = 'protein',
            fname = os.path.join(common.ROOT, 'data', 'trembl3.tab'),
            col_a = 0,
            col_b = 1,
            separator = '\t',
            header = 0,
        ),
    ('genesymbol', 'swissprot'):
        input_formats.FileMapping(
            id_type_a = 'genesymbol',
            id_type_b = 'swissprot',
            entity_type = 'protein',
            fname = os.path.join(common.ROOT, 'data', 'swissprot3.tab'),
            col_a = 1,
            col_b = 0,
            separator = '\t',
            header = 0,
        )
    ('genesymbol', 'uniprot'):
        input_formats.FileMapping(
            id_type_a = 'genesymbol',
            id_type_b = 'uniprot',
            entity_type = 'protein',
            fname = os.path.join(common.ROOT, 'data', 'uniprot3.tab'),
            col_a = 1,
            col_b = 0,
            separator = '\t',
            header = 0,
        ),
    ('genesymbol-fallback', 'uniprot'):
        input_formats.FileMapping(
            id_type_a = 'genesymbol-fallback',
            id_type_b = 'uniprot',
            entity_type = 'protein',
            fname = os.path.join(
                common.ROOT, 'data', 'human-genesymbol-all.tab'
            ),
            col_a = 1,
            col_b = 0,
            separator = '\t',
            header=0,
        ),
    ('refseq', 'uniprot'):
        input_formats.FileMapping(
            id_type_a = 'refseq',
            id_type_b = 'uniprot',
            entity_type = 'protein',
            fname = os.path.join(
                common.ROOT, 'data', 'uniprot-refseq-human-1.tab'
            ),
            col_a = 1,
            col_b = 0,
            separator = '\t',
            header=0,
            ncbi_tax_id = 9606,
        ),
    ('entrez', 'uniprot'):
        input_formats.FileMapping(
            id_type_a = 'entrez',
            id_type_b = 'uniprot',
            entity_type = 'protein',
            fname = os.path.join(common.ROOT, 'data', 'entrez_uniprot.csv'),
            col_a = 1,
            col_b = 0,
            separator = ';',
            header = 0,
        ),
}


uniprot = {
    ('embl', 'uniprot'):
        input_formats.UniprotMapping(
            id_type = 'embl',
        ),
    ('genesymbol', 'uniprot'): 
        input_formats.UniprotMapping(
            id_type = 'genesymbol',
        ),
    ('genesymbol-syn', 'uniprot'):
        input_formats.UniprotMapping(
            id_type = 'genesymbol-syn'
        ),
    ('entrez', 'uniprot'):
        input_formats.UniprotMapping(
            id_typye = 'entrez',
        ),
    ('hgnc', 'uniprot'):
        input_formats.UniprotMapping(
            id_type = 'hgnc',
        ),
    ('enst', 'uniprot'):
        input_formats.UniprotMapping(
            id_type = 'enst',
        ),
    ('refseqp', 'uniprot'):
        input_formats.UniprotMapping(
            id_type = 'refseqp',
        ),
    ('uniprot-entry', 'uniprot'):
        input_formats.UniprotMapping(
            id_type = 'uniprot-entry',
        ),
    ('protein-name', 'uniprot'):
        input_formats.UniprotMapping(
            id_type = 'protein-name',
        ),
    ('protein-name-all', 'uniprot'):
        input_formats.UniprotMapping(
            id_type = 'protein-name',
            swissprot = None,
        )
}


mirbase = {
    ('mir-mat-name', 'mirbase'):
        input_formats.FileMapping(
            id_type_a = 'mir-mat-name',
            id_type_b = 'mirbase',
            fname = 'mirbase_mature'
            col_a = 1,
            col_b = 0,
            separator = None,
            header = 0,
        ),
    ('mir-name', 'mir-pre'):
        input_formats.FileMapping(
            id_type_a = 'mir-name',
            id_type_b = 'mir-pre',
            fname = 'mirbase_precursor',
            col_a = 1,
            col_b = 0,
            separator = None,
            header = 0,
        ),
    ('mir-pre', 'mirbase'):
        input_formats.FileMapping(
            id_type_a = 'mir-pre',
            id_type_b = 'mirbase',
            fname = 'mirbase_ids',
            col_a = 1,
            col_b = 0,
            separator = None,
            header = 0,
        ),
    ('mir-name', 'mirbase'):
        input_formats.FileMapping(
            id_type_a = 'mir-name',
            id_type_b = 'mirbase',
            fname = 'mirbase_precursor_to_mature',
            col_a = 1,
            col_b = 0,
            separator = None,
            header = 0,
        ),
}


basic = {
    ('uniprot-sec', 'uniprot-pri'): input_formats.FileMapping(
        'get_uniprot_sec', 0, 1, None, header=0, ncbi_tax_id = 0),
    ('genesymbol', 'trembl'): input_formats.UniprotMapping(
        'genesymbol', swissprot='no', bi_directional = True),
    ('genesymbol', 'swissprot'): input_formats.UniprotMapping('genesymbol'),
    ('genesymbol-syn', 'swissprot'):
    input_formats.UniprotMapping('genesymbol-syn'),
    ('genesymbol', 'uniprot'): input_formats.UniprotMapping(
        'genesymbol', bi_directional = True, swissprot=None)
}


# the part below should be removed
#
# this is all what is needed for corrections of unirpot ids
# i.e. to get primary swissprot id for all proteins
uniprot_old = [
    {
        id_type_a = 'uniprot-sec',
        id_type_b = 'uniprot-pri',
        entity_type = 'protein',
        'src': 'file',
        'par': input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'sec_ac.txt'),
            0, 1, None, header = 0
        )
    },
    {
        id_type_a = 'trembl',
        id_type_b = 'genesymbol',
        entity_type = 'protein',
        'src': 'file',
        'par': input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'trembl3.tab'),
            0, 1, '\t', header = 0
        )
    },
    {
        id_type_a = 'genesymbol',
        id_type_b = 'swissprot',
        entity_type = 'protein',
        'src': 'file',
        'par': input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'swissprot3.tab'),
            1,
            0,
            '\t',
            header = 0
        )
    },
    {
        id_type_a = 'genesymbol-fallback',
        id_type_b = 'uniprot',
        entity_type = 'protein',
        'src': 'file',
        'par': input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'human-genesymbol-all.tab'),
            1,
            0,
            '\t',
            header = 0
        ),
    }
]
