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


misc = [
    {
        "one": "uniprot-sec",
        "two": "uniprot-pri",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'sec_ac.txt'),
            0,
            1,
            None,
            header=0)
    },
    {
        "one": "trembl",
        "two": "genesymbol",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'trembl3.tab'),
            0,
            1,
            "\t",
            header=0)
    },
    {
        "one": "genesymbol",
        "two": "swissprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'swissprot3.tab'),
            1,
            0,
            "\t",
            header=0)
    },
    {
        "one": "genesymbol",
        "two": "uniprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'uniprot3.tab'),
            1,
            0,
            "\t",
            header=0,
            bi_directional = True)
    },
    {
        "one": "genesymbol-fallback",
        "two": "uniprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'human-genesymbol-all.tab'),
            1,
            0,
            "\t",
            header=0)
    },
    {
        "one": "refseq",
        "two": "uniprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'uniprot-refseq-human-1.tab'),
            1,
            0,
            "\t",
            header=0)
    },
    {
        "one": "entrez",
        "two": "uniprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'entrez_uniprot.csv'),
            1,
            0,
            ";",
            header=0)
    },
]

uniprot = {
    ('embl', 'uniprot'): input_formats.UniprotMapping('embl'),
    ('genesymbol', 'uniprot'): input_formats.UniprotMapping(
        'genesymbol', bi_directional = True),
    ('genesymbol-syn', 'uniprot'):
    input_formats.UniprotMapping('genesymbol-syn'),
    ('entrez', 'uniprot'): input_formats.UniprotMapping('entrez'),
    ('hgnc', 'uniprot'): input_formats.UniprotMapping('hgnc'),
    ('enst', 'uniprot'): input_formats.UniprotMapping('enst'),
    ('refseqp', 'uniprot'): input_formats.UniprotMapping('refseqp'),
    ('uniprot-entry', 'uniprot'):
    input_formats.UniprotMapping('uniprot-entry'),
    ('protein-name', 'uniprot'): input_formats.UniprotMapping('protein-name'),
    ('protein-name-all', 'uniprot'): input_formats.UniprotMapping(
        'protein-name', swissprot=None)
}

mirbase = {
    ('mir-mat-name', 'mirbase'): input_formats.FileMapping(
        'mirbase_mature', 1, 0, None, header = 0),
    ('mir-name', 'mir-pre'): input_formats.FileMapping(
        'mirbase_precursor', 1, 0, None, header = 0),
    ('mir-pre', 'mirbase'): input_formats.FileMapping(
        'mirbase_ids', 1, 0, None, header = 0),
    ('mir-name', 'mirbase'): input_formats.FileMapping(
        'mirbase_precursor_to_mature',
        1, 0, None, header = 0),
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

# this is all what is needed for corrections of unirpot ids
# i.e. to get primary swissprot id for all proteins
uniprot_old = [
    {
        "one": "uniprot-sec",
        "two": "uniprot-pri",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'sec_ac.txt'),
            0, 1, None, header = 0
        )
    },
    {
        "one": "trembl",
        "two": "genesymbol",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'trembl3.tab'),
            0, 1, "\t", header = 0
        )
    },
    {
        "one": "genesymbol",
        "two": "swissprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'swissprot3.tab'),
            1,
            0,
            "\t",
            header = 0
        )
    },
    {
        "one": "genesymbol-fallback",
        "two": "uniprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(
            os.path.join(common.ROOT, 'data', 'human-genesymbol-all.tab'),
            1,
            0,
            "\t",
            header = 0
        ),
    }
]
