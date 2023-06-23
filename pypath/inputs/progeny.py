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

import pypath.share.session as session
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping
import pypath.inputs.rdata as rdata

_logger = session.Logger(name = 'progeny_input')
_log = _logger._log


def progeny_raw(organism = 9606):
    """
    Pathway responsive genes: signatures based on transcriptomics data
    from PROGENy (https://github.com/saezlab/progeny).

    Args
        organism (int,str): Name or NCBI Taxonomy ID of the organism. Human
            and mouse are supported.

    Returns
        (pandas.DataFrame): A data frame of genes, pathways, weights and
            p-values for each association.
    """

    _organism = taxonomy.ensure_common_name(organism)

    if _organism not in ('Human', 'Mouse'):

        msg = (
            'Wrong organism: `%s`; '
            'only human and mouse are available.' % organism
        )
        _log(msg)
        raise ValueError(msg)

    _organism = _organism.lower()

    url = urls.urls['progeny']['url'] % _organism
    c = curl.Curl(url, large = True, silent = False)

    rdata_path = c.fileobj.name
    c.fileobj.close()

    rdata_parsed = rdata.rdata.parser.parse_file(rdata_path)
    rdata_converted = rdata.rdata.conversion.convert(rdata_parsed)

    key = 'model_%s_full' % _organism

    return rdata_converted[key]


def progeny_annotations(organism = 9606):
    """
    Pathway responsive genes: signatures based on transcriptomics data
    from PROGENy (https://github.com/saezlab/progeny).

    Args
        organism (int,str): Name or NCBI Taxonomy ID of the organism. Human
            and mouse are supported.

    Returns
        (dict): Dict of sets, keys are UniProt IDs, values are pathway
            association records, each with a weight and p-value.
    """

    record = collections.namedtuple(
        'ProgenyAnnotation',
        (
            'pathway',
            'weight',
            'p_value',
        )
    )

    raw = progeny_raw(organism = organism)
    result = collections.defaultdict(set)

    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)

    for rec in raw.itertuples():

        uniprots = mapping.map_name(
            rec.gene,
            'genesymbol',
            'uniprot',
            ncbi_tax_id = ncbi_tax_id,
        )

        annot = record(
            pathway = rec.pathway,
            weight = rec.weight,
            p_value = rec[4], # omg, stupid pandas
        )

        for uniprot in uniprots:

            result[uniprot].add(annot)

    return dict(result)
