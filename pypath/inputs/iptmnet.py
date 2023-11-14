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

import re
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath_common._constants as _const
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping
import pypath.inputs.common as inputs_common


resite = re.compile(r'([A-Z])(\d+)')


IptmnetInteraction = collections.namedtuple(
    'IptmnetInteraction',
    [
        'enzyme',
        'substrate',
        'enzyme_isoform',
        'substrate_isoform',
        'ptm_type',
        'resaa',
        'resnum',
        'score',
        'references',
    ]
)


def iptmnet_interactions(organism = 9606):

    ptm_url = urls.urls['iptmnet']['ptms']
    score_url = urls.urls['iptmnet']['scores']

    c = curl.Curl(score_url, large = True, silent = False)

    scores = {}

    for line in c.result:

        line = line.strip('\n\r').split('\t')

        if not line[2]:

            continue

        site = resite.match(line[1])

        if not site:

            continue

        resaa, resnum = site.groups()

        resnum = int(resnum)
        score = int(line[4])
        substrate, isoform = inputs_common._try_isoform(line[0])
        enzyme = line[2]

        key = (
            enzyme,
            substrate,
            isoform,
            line[3].lower(), # PTM type
            resaa,
            resnum,
        )

        scores[key] = score

    c = curl.Curl(ptm_url, large = True, silent = False)

    for line in c.result:

        line = line.strip('\n\r').split('\t')

        if not line or not line[6]:

            continue

        organism_ = line[4].strip()
        ncbi_tax_id = (
            taxonomy.ensure_ncbi_tax_id(organism_)
                if organism_ else
            _const.NOT_ORGANISM_SPECIFIC
        )

        if organism and ncbi_tax_id != organism:

            continue

        substrate, s_isoform = inputs_common._try_isoform(line[2])
        ptm_type = line[0].lower()

        enzyme, e_isoform = inputs_common._try_isoform(line[6])

        enzyme_ids = (
            mapping.map_name(
                line[6],
                'pro',
                'uniprot',
                ncbi_tax_id = organism,
            )
                if line[6].startswith('PR:') else
            (enzyme,)
        )

        refs = line[9].split(',')
        resnum, resaa = resite.match(line[5]).groups()

        key = (
            line[6],
            substrate,
            isoform,
            ptm_type,
            resaa,
            resnum,
        )

        score = scores[key] if key in scores else None

        for _enzyme in enzyme_ids:

            yield IptmnetInteraction(
                enzyme = _enzyme,
                substrate = substrate,
                enzyme_isoform = e_isoform,
                substrate_isoform = s_isoform,
                ptm_type = ptm_type,
                resaa = resaa,
                resnum = resnum,
                score = score,
                references = refs,
            )
