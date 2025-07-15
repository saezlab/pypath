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
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.cache as cache
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping


def celltalkdb_download(organism = 9606):
    """
    Downloads ligand-receptor pairs from CellTalkDB GitHub repository.

    :param int,str organism:
        Human and mouse supported, in case of incomprehensible value will
        fall back to human.

    :return:
        Generator yielding dataset specific records as named tuples.
    """

    taxid = taxonomy.ensure_common_name(organism)
    organism_name = taxid.lower() if taxid in {'Mouse', 'Human'} else 'human'

    cache_fname = 'celltalkdb_%s_lr_pair.rds' % organism_name
    cache_path = os.path.join(cache.get_cachedir(), cache_fname)

    url = urls.urls['celltalkdb']['url'] % organism_name

    c = curl.Curl(
        url = url,
        cache = cache_path,
        silent = False,
        large = True,
    )

    # Read the R data file using pyreadr
    import pyreadr
    rdata = pyreadr.read_r(c.fileobj.name)
    
    # Get the data frame (typically there's only one object in the .rds file)
    df = list(rdata.values())[0]
    
    # Convert to named tuples
    header = df.columns.tolist()
    record = collections.namedtuple('CellTalkDbRecord', header)

    for _, row in df.iterrows():
        yield record(*row.values)


def celltalkdb_interactions(organism = 9606):
    """
    Retrieves ligand-receptor interactions from CellTalkDB
    http://tcm.zju.edu.cn/celltalkdb/index.php

    :param int,str organism:
        Human and mouse supported, in case of incomprehensible value will
        fall back to human.

    :return:
        List of interactions as named tuples.
    """

    CellTalkDBInteraction = collections.namedtuple(
        'CellTalkDBInteraction',
        [
            'ligand_genesymbol',
            'receptor_genesymbol',
            'reference',
        ]
    )

    return [
        CellTalkDBInteraction(
            ligand_genesymbol = rec.ligand_gene_symbol,
            receptor_genesymbol = rec.receptor_gene_symbol,
            reference = rec.evidence,
        )
        for rec in celltalkdb_download(organism = organism)
    ]


def celltalkdb_annotations(organism = 9606):
    """
    Retrieves annotation of protein ligand and receptor roles from CellTalkDB
    http://tcm.zju.edu.cn/celltalkdb/index.php

    :param int,str organism:
        Human and mouse supported, in case of incomprehensible value will
        fall back to human.

    :return:
        Dictionary of annotations with UniProt IDs as keys.
    """

    CellTalkDBAnnotation = collections.namedtuple(
        'CellTalkDBAnnotation',
        [
            'role',
            'pmid',
        ]
    )

    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)
    ncbi_tax_id = ncbi_tax_id if ncbi_tax_id in {9606, 10090} else 9606

    annot = collections.defaultdict(set)

    for rec in celltalkdb_download(organism = ncbi_tax_id):

        for role in ('ligand', 'receptor'):

            uniprots = mapping.map_name(
                getattr(rec, '%s_gene_symbol' % role),
                'genesymbol',
                'uniprot',
                ncbi_tax_id = ncbi_tax_id,
            )

            for uniprot in uniprots:

                annot[uniprot].add(
                    CellTalkDBAnnotation(
                        role = role,
                        pmid = rec.evidence,
                    )
                )

    return dict(annot)
