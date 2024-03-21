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
import itertools

import pypath.utils.mapping as mapping
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.uniprot_db as uniprot_db
import pypath.share.session as session

_logger = session.Logger(name = 'homologene_input')
_log = _logger._log


def get_homologene():
    """
    Downloads the latest release of the NCBI HomoloGene database.
    Returns file pointer.
    """

    url = urls.urls['homologene']['url_rescued']

    c = curl.Curl(
        url = url,
        silent = False,
        large = True,
        timeout = 1800,
        ignore_content_length = True,
    )

    return c.result


def homologene_dict(source, target, id_type):
    """
    Returns orthology translation table as dict, obtained
    from NCBI HomoloGene data.

    Args
        source (int): NCBI Taxonomy ID of the source species (keys).
        target (int): NCBI Taxonomy ID of the target species (values).
        id_type (str): ID type to be used in the dict. Possible values:
            'RefSeq', 'Entrez', 'GI', 'GeneSymbol'.

    Returns
        Dict of sets: keys are IDs of the source organism, values are sets
        of IDs of the target organism.
    """
    ids = {
        'refseq': 5,
        'refseqp': 5,
        'genesymbol': 3,
        'gi': 4,
        'entrez': 2
    }

    try:
        id_col = ids[id_type.lower()]
    except KeyError:
        _log(
            'Unknown ID type: `%s`. Please use RefSeq, '
            'Entrez, GI or GeneSymbol.' % id_type
        )
        raise

    hg = get_homologene()
    hgroup = None
    result = collections.defaultdict(set)

    for l in hg:

        l = l.strip().split('\t')
        this_hgroup = l[0].strip()

        if this_hgroup != hgroup:

            this_source = None
            this_target = None
            hgroup = this_hgroup

        this_taxon = int(l[1].strip())

        if this_taxon == source:

            this_source = l[id_col]

        elif this_taxon == target:

            this_target = l[id_col]

        if (
            this_source and
            this_target
        ):

            result[this_source].add(this_target)

    return dict(result)


def homologene_uniprot_dict(source, target, only_swissprot = True):
    """
    Returns orthology translation table as dict from UniProt to Uniprot,
    obtained from NCBI HomoloGene data. Uses RefSeq and Entrez IDs for
    translation.

    Args
        source (int): NCBI Taxonomy ID of the source species (keys).
        target(int): NCBI Taxonomy ID of the target species (values).
        only_swissprot (bool): Use only SwissProt IDs.

    Returns
        Dict of sets: keys are UniProt IDs of the source organism, values
        are sets of UniProt IDs of the target organism.
    """

    result = {}

    hge = homologene_dict(source, target, 'entrez')
    hgr = homologene_dict(source, target, 'refseq')

    all_source = set(uniprot_db.all_uniprots(
        organism = source,
        swissprot = 'YES',
    ))

    if not only_swissprot:

        all_source_trembl = uniprot_db.all_uniprots(
            organism = source,
            swissprot = 'NO',
        )
        all_source.update(set(all_source_trembl))

    for u in all_source:

        source_e = mapping.map_name(u, 'uniprot', 'entrez', source)
        source_r = mapping.map_name(u, 'uniprot', 'refseqp', source)
        target_u = set()
        target_r = set()
        target_e = set()

        for e in source_e:

            if e in hge:

                target_e.update(hge[e])

        for r in source_r:

            if r in hgr:

                target_r.update(hgr[r])

        for e in target_e:

            target_u.update(
                mapping.map_name(e, 'entrez', 'uniprot', target)
            )

        for r in target_r:

            target_u.update(
                mapping.map_name(e, 'refseqp', 'uniprot', target)
            )


        target_u = (
            itertools.chain(
                *map(
                    lambda tu:
                        mapping.map_name(tu, 'uniprot', 'uniprot', target),
                    target_u
                )
            )
        )

        result[u] = sorted(list(target_u))

    return result
