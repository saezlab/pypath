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

import pyreadr

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def talklr_raw(putative = True):

    url = urls.urls['talklr']['url']
    c = curl.Curl(url, large = True, silent = False)
    rdata_path = c.fileobj.name
    c.fileobj.close()

    rdata = pyreadr.read_r(rdata_path)['receptor_ligand']
    rdata.columns = [col.replace('.', '_') for col in rdata.columns]

    if not putative:

        rdata = rdata[rdata.Pair_Evidence != 'putative']

    return rdata


def talklr_interactions(putative = True):

    TalklrInteraction = collections.namedtuple(
        'TalklrInteraction',
        ['ligand', 'receptor', 'pmids', 'resources', 'putative'],
    )

    resource_columns = {
        'DLRP': ('DLRP',),
        'Guide2Pharma': ('IUPHAR',),
        'HPRD': ('HPRD',),
        'STRING': ('STRING_binding', 'STRING_experiment'),
        'HPMR': ('HPMR',),
    }

    raw = talklr_raw(putative = putative)

    for rec in raw.itertuples():

        resources = tuple(
            resource
            for resource, labels in resource_columns.items()
            if any(
                common.is_str(getattr(rec, lab))
                for lab in labels
            )
        )
        pmids = (
            tuple(pmid.strip() for pmid in rec.PMID_Manual.split(','))
                if common.is_str(rec.PMID_Manual) else
            ()
        )

        yield(
            TalklrInteraction(
                ligand = rec.Ligand_ApprovedSymbol,
                receptor = rec.Receptor_ApprovedSymbol,
                pmids = pmids,
                resources = resources,
                putative = rec.Pair_Evidence == 'putative',
            )
        )


def talklr_annotations(putative = True):

    TalklrAnnotation = collections.namedtuple(
        'TalklrAnnotation',
        ['role', 'pmid', 'putative']
    )


    result = collections.defaultdict(set)

    for ia in talklr_interactions(putative = putative):

        for role in ('ligand', 'receptor'):

            uniprots = mapping.map_name(
                getattr(ia, role),
                'genesymbol',
                'uniprot',
                ncbi_tax_id = 9606,
            )

            for uniprot in uniprots:

                for pmid in (ia.pmids or (None,)):

                    result[uniprot].add(
                        TalklrAnnotation(
                            role = role,
                            pmid = pmid,
                            putative = ia.putative,
                        )
                    )

    return dict(result)
