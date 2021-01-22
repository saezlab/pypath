#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import collections

import pyreadr

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls


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
        ['ligand', 'receptor', 'pmids', 'resources'],
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
            )
        )

