#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from __future__ import annotations

import gzip
import collections
import pandas as pd

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common


def adrecs_drug_identifiers() -> list[tuple]:
    """
    Drug identifiers from the AdReCS database.

    IUPAC name, synonyms, DrugBank, MeSH, KEGG and TDD IDs of drugs.
    http://www.bio-add.org/ADReCS/index.jsp

    Returns:
        Dict of sets of tuples of drug identifiers. Top level keys are
        PubChem CIDs.
    """

    AdrecsDrug = collections.namedtuple(
        'AdrecsDrug',
        ('drug', 'synonyms', 'drugbank', 'pubchem_cid', 'mesh', 'kegg', 'tdd'),
    )

    url = urls.urls['adrecs']['drug_information']
    c = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(c.outfile, cell_range = 'A1:H2527')
    result = []

    notavail = lambda x: None if x == 'Not Available' else x

    for line in contents[1:]:

        line = [notavail(x) for x in line]

        line[2] = (
            tuple(sorted(x.strip() for x in line[2].split('|')))
                if line[2] else
            ()
        )

        result.append(AdrecsDrug(*line[1:]))

    return result


def adrecs_terminology() -> list[tuple]:
    """
    Retrieves ADR terminology and hierachy.
    Returns:
        list of namedtuples with all attributes
        (adrecs_id, adr_id, adr_term, adr_synonyms, meddra_code)
    """

    url = urls.urls['adrecs']['url_terminology']
    path = curl.Curl(
        url,
        silent = False,
        large = True
    ).fileobj.name

    with gzip.open(path, 'rb') as f_in:

        df = pd.read_excel(f_in)
        df.replace({'Not Available': None}, inplace=True)
        df['ADR_SYNONYMS'] = df['ADR_SYNONYMS'].str.split('|')
        df['ADR_SYNONYMS'] = df['ADR_SYNONYMS'].apply(lambda x: tuple([element.strip() for element in x]) if x is not None else None)
        df.columns= df.columns.str.lower()

        return list(df.itertuples(name='AdrecsOntology', index=False))


def adrecs_drugs() -> list[tuple]:
    """
    Retrieves Drug-ADR pairs
    Returns:
        list of namedtuples with drug_id and adr_id
    """

    url = urls.urls['adrecs']['url_adrecs_drugs']
    c = curl.Curl(url, large = True, silent = False)
    drugs = c.result

    fields = [
        'drug_id',
        'adr_id'
    ]

    result = set()
    record = collections.namedtuple('AdrecsDrug', fields)

    # requisite fields from all_fields
    indices = [0, 1]

    for line in drugs:

        if not line.strip():
            continue

        line = line.strip().split('\t')
        line = dict(zip(fields, (line[i] for i in indices)))

        result.add(
            record(**dict(zip(fields, (line.get(f, None) for f in fields))))
        )

    return list(result)
