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

"""
Web scraper for the CancerDrugs_DB database.
"""

import csv
import collections
import itertools
import re

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.share.session as session

_logger = session.Logger(name = 'cancerdrugsdb_input')
_log = _logger._log


def cancerdrugsdb_download():
    """
    Downloads a curated set of interactions of cancer drugs licensed
    in most parts of the world with gene targets (where available).
    From https://www.anticancerfund.org/en/cancerdrugs-db.
    This function downloads a single dataset.

    Args
        None.

    Returns
        A list of dicts, each is a record as it provided by the database.
    """

    url = urls.urls['cancerdrugs_db']['url']

    c = curl.Curl(url, large = True, silent = False)

    return list(csv.DictReader(c.result, delimiter = '\t'))


def cancerdrugsdb_interactions():
    """
    Returns drug-gene interactions from Cancer Drugs Database
    (https://www.anticancerfund.org/en/cancerdrugs-db).

    Args
        None.

    Returns
        List of named tuples, each describing a drug-gene interaction.
        Identifiers of type PubChem and UniProt.
    """

    # note: drug-target interactions are from GDIdb. if they're in OnmiPath
    # already, do we skip the interactions altogether, or is it better to
    # keep them here for consistency? the dataset is very small anyways.
    # if we keep them, do we test for consistency between the interactions?

    reid = re.compile(r'>(\w+)<')

    def strip_id(field):

        match = reid.search(field)

        return match.group(1) if match else None


    def yes_no(field):

        return field == 'Y'


    CancerDrugsInteraction = collections.namedtuple(
        'CancerdrugsdbInteraction',
        (
            'drug_pubchem',
            'drug_chembl',
            'drug_drugbank',
            'drug_label',
            'target_uniprot',
            'ema_approved',
            'fda_approved',
            'european_national_approved',
            'who_approved',
            'generic',
            'approval_year',
            'indications',
        ),
    )

    result = []
    unmapped_drug = []
    no_targets = []

    data = cancerdrugsdb_download()

    for rec in data:

        chembl = strip_id(rec.get('ChEMBL'))
        drugbank = strip_id(rec.get('DrugBank ID'))

        if chembl is None:

            unmapped_drug.append(rec.get('Product'))
            continue

        pubchems = mapping.map_name(chembl, 'chembl', 'pubchem')

        targets = rec.get('Targets')

        if not targets:

            no_targets.append(rec.get('Product'))
            continue

        target_uniprots = mapping.map_names(
            (tar.strip() for tar in targets.split(';')),
            'genesymbol',
            'uniprot',
        )

        for pubchem, uniprot in itertools.product(pubchems, target_uniprots):

            result.append(
                CancerDrugsInteraction(
                    drug_pubchem = pubchem,
                    drug_chembl = chembl,
                    drug_drugbank = drugbank,
                    drug_label = rec.get('Product'),
                    target_uniprot = uniprot,
                    ema_approved = yes_no(rec.get('EMA')),
                    fda_approved = yes_no(rec.get('FDA')),
                    european_national_approved = yes_no(rec.get('EN')),
                    who_approved = yes_no(rec.get('WHO')),
                    generic = yes_no(rec.get('Generic')),
                    approval_year = int(rec.get('Year')) if rec.get('Year') else None,
                    indications = tuple(
                        i.strip()
                        for i in rec.get('Indications').split(';')
                    ),
                )
            )


    _log(
        'Could not find CHEMBL IDs for %u '
        'CancerDrugs_DB Products.' % len(unmapped_drug)
    )

    _log(
        '%u CancerDrugs_DB Products had no targets.'
        % len(no_targets)
    )

    return result


def cancerdrugsdb_annotations():
    """
    Returns drug annotations from CancerDrugs_DB.

    Args
        None.

    Returns
        (dict): Keys are PubChem IDs, values are sets of annotations.
    """

    record = collections.namedtuple(
        'CancerdrugsdbAnnotation',
        (
            'drug_label',
            'ema_approved',
            'fda_approved',
            'european_national_approved',
            'who_approved',
            'generic',
            'approval_year',
            'indications',
        ),
    )

    result = collections.defaultdict(set)
    data = cancerdrugsdb_interactions()

    for rec in data:

        result[rec.drug_pubchem].add(
            record(
                **dict(
                    i
                    for i in rec._asdict().items()
                    if i[0] in record._fields
                )
            )
        )

    return dict(result)
