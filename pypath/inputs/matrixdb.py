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

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.inputs.uniprot_db as uniprot_db


def matrixdb_interactions(organism = 9606):

    url = urls.urls['matrixdb']['url']
    c = curl.Curl(url, silent = False, large = True)
    f = c.result
    i = []
    lnum = 0

    for l in f:
        if lnum == 0:
            lnum += 1

            continue

        l = l.replace('\n', '').replace('\r', '')
        l = l.split('\t')
        specA = 0 if l[9] == '-' else int(l[9].split(':')[1].split('(')[0])
        specB = 0 if l[10] == '-' else int(l[10].split(':')[1].split('(')[0])

        if organism is None or (specA == organism and specB == organism):
            pm = [
                p.replace('pubmed:', '') for p in l[8].split('|')
                if p.startswith('pubmed:')
            ]
            met = [
                m.split('(')[1].replace(')', '').strip('"')
                for m in l[6].split('|')
                if '(' in m
            ]
            l = [l[0], l[1]]
            interaction = ()

            for ll in l:
                ll = ll.split('|')
                uniprot = ''

                for lll in ll:
                    nm = lll.split(':')

                    if nm[0] == 'uniprotkb' and len(nm[1]) == 6:
                        uniprot = nm[1]

                interaction += (uniprot, )

            interaction += ('|'.join(pm), '|'.join(met))

            if len(interaction[0]) > 5 and len(interaction[1]) > 5:
                i.append(list(interaction))

        lnum += 1

    f.close()

    return i


def _matrixdb_protein_list(category, organism = 9606):
    """
    Returns a set of proteins annotated by MatrixDB.

    :arg str category:
        The protein annotation category. Possible values: `ecm`, `membrane`
        or `secreted`.
    """

    url = urls.urls['matrixdb']['%s_proteins' % category]
    c = curl.Curl(url, silent = False, large = True)

    proteins = set()

    # header row
    _ = next(c.result)

    for l in c.result:
        if not l:
            continue

        proteins.add(
            l.strip().replace('"', '').split('\t')[0]
        )

    proteins = mapping.map_names(proteins, 'uniprot', 'uniprot')

    if organism:

        uniprots = uniprot_db.all_uniprots(
            organism = organism,
            swissprot = True,
        )
        proteins = proteins & set(uniprots)

    return proteins


def matrixdb_membrane_proteins(organism = 9606):
    """
    Returns a set of membrane protein UniProt IDs retrieved from MatrixDB.
    """

    return _matrixdb_protein_list('membrane', organism = organism)


def matrixdb_secreted_proteins(organism = 9606):
    """
    Returns a set of secreted protein UniProt IDs retrieved from MatrixDB.
    """

    return _matrixdb_protein_list('secreted', organism = organism)


def matrixdb_ecm_proteins(organism = 9606):
    """
    Returns a set of ECM (extracellular matrix) protein UniProt IDs
    retrieved from MatrixDB.
    """

    return _matrixdb_protein_list('ecm', organism = organism)


def matrixdb_annotations(organism = 9606):

    MatrixdbAnnotation = collections.namedtuple(
        'MatrixdbAnnotation',
        ('mainclass',),
    )
    annot = collections.defaultdict(set)

    for cls in ('membrane', 'secreted', 'ecm'):
        cls_annot = MatrixdbAnnotation(mainclass = cls)

        method = globals()['matrixdb_%s_proteins' % cls]

        for uniprot in method(organism = organism):
            annot[uniprot].add(cls_annot)

    return dict(annot)
