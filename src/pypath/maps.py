#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2016 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import os

import pypath.common as common
import pypath.input_formats as input_formats

mapList = [
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
            bi=True)
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
    {
        "one": "hgnc",
        "two": "uniprot",
        "typ": "protein",
        "src": "mysql",
        "par": input_formats.MysqlMapping(
            "hgnc_new", "gsy", "u", "mapping", None, bi=True)
    },
    {
        "one": "uniprot",
        "two": "hgncapprov",
        "typ": "protein",
        "src": "mysql",
        "par": input_formats.MysqlMapping(
            "hgnc_prim", "u", "gsy", "mapping", None, bi=True)
    },
    {
        "one": "uniprot",
        "two": "hgnc",
        "typ": "protein",
        "src": "mysql",
        "par":
        input_formats.MysqlMapping("hgnc_names", "u", "gsy", "mapping", None)
    }
    #,
    #{
    #"one": "uniprot",
    #"two": "hgnc",
    #"typ": "protein",
    #"src": "mysql",
    #"par": input_formats.UniprotMapping("hgnc_names","u","gsy","mapping",None)
    #}
]

mapListUniprot = {
    ('embl', 'uniprot'): input_formats.UniprotMapping('embl'),
    ('genesymbol', 'uniprot'): input_formats.UniprotMapping(
        'genesymbol', bi=True),
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

mapListMirbase = {
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

mapListBasic = {
    ('uniprot-sec', 'uniprot-pri'): input_formats.FileMapping(
        'get_uniprot_sec', 0, 1, None, header=0, ncbi_tax_id = 0),
    ('genesymbol', 'trembl'): input_formats.UniprotMapping(
        'genesymbol', swissprot='no', bi=True),
    ('genesymbol', 'swissprot'): input_formats.UniprotMapping('genesymbol'),
    ('genesymbol-syn', 'swissprot'):
    input_formats.UniprotMapping('genesymbol-syn'),
    ('genesymbol', 'uniprot'): input_formats.UniprotMapping(
        'genesymbol', bi=True, swissprot=None)
}

# this is all what is needed for corrections of unirpot ids
# i.e. to get primary swissprot id for all proteins
mapListUniprotOld = [{
    "one": "uniprot-sec",
    "two": "uniprot-pri",
    "typ": "protein",
    "src": "file",
    "par": input_formats.FileMapping(
        os.path.join(common.ROOT, 'data', 'sec_ac.txt'), 0, 1, None, header=0)
}, {
    "one": "trembl",
    "two": "genesymbol",
    "typ": "protein",
    "src": "file",
    "par": input_formats.FileMapping(
        os.path.join(common.ROOT, 'data', 'trembl3.tab'), 0, 1, "\t", header=0)
}, {
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
}, {
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
}]

otherMappings = [{
    "one": "entrez",
    "two": "uniprot",
    "typ": "protein",
    "src": "mysql",
    "par":
    input_formats.MysqlMapping("geneid", "geneid", "u", "mapping", "ncbi")
}, {
    "one": "uniprot",
    "two": "genesymbol",
    "typ": "protein",
    "src": "mysql",
    "par":
    input_formats.MysqlMapping("uniprot_gs", "u", "gs", "mapping", "ncbi")
}]
