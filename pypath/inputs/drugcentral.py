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

from typing import Optional, Union

import csv
import collections

import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.common as common
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy

_logger = session.Logger(name = 'drugcentral_input')
_log = _logger._log


def drugcentral_drugs() -> list[tuple]:
    """
    Drug names and structures from Drug Central.

    Returns
        List of drugs, each represented by a named tuple.
    """

    DrugcentralDrug = collections.namedtuple(
        'DrugcentralDrug',
        (
            'drugcentral',
            'inn',
            'cas',
            'smiles',
            'inchikey',
            'inchi',
        )
    )

    url = urls.urls['drugcentral']['SMILES_InChI']
    c = curl.Curl(url, large = True, silent = False)
    drugs = list(csv.DictReader(c.result, delimiter = '\t'))

    result = [
        DrugcentralDrug(
            drugcentral = drug['ID'],
            inn = drug['INN'],
            cas = drug['CAS_RN'],
            smiles = drug['SMILES'],
            inchikey = drug['InChIKey'],
            inchi = drug['InChI'],
        )
        for drug in drugs
    ]

    return result


def drugcentral_interactions(
        organism: Optional[Union[str, int]] = None,
        comments: bool = False,
    ) -> list[tuple]:
    """
    Retrieves drug-target interactions from Drug Central.

    Args
        organism:
            Organism name or NCBI Taxonomy ID. If not provided,
            all organisms will be retained.
        comments:
            Include comments in the result.

    Returns
        List of drug-target relationships, represented as named tuples.
    """

    DrugcentralInteraction = collections.namedtuple(
        'DrugcentralInteraction',
        (
            'drug',
            'drug_name',
            'uniprot',
            'target_type',
            'canonical',
            'act_value',
            'act_type',
            'relation',
            'effect',
            'tdl',
            'organism',
            'comment',
        ),
    )

    url = urls.urls['drugcentral']['interactions']
    c = curl.Curl(url, large = True, silent = False)
    interactions = list(csv.DictReader(c.result, delimiter = '\t'))

    organism_latin = taxonomy.ensure_latin_name(organism)

    if organism and not organism_latin:

        msg = f'Could not find latin name for organism: `{organism}`.'
        _log(msg)

    drugs = dict(
        (d.drugcentral, d)
        for d in drugcentral_drugs()
    )

    result = [
        DrugcentralInteraction(
            drug = drugs.get(i['STRUCT_ID'], None),
            drug_name = i['DRUG_NAME'],
            uniprot = uniprot,
            target_type = i['TARGET_CLASS'],
            canonical = i['MOA'] == '1',
            act_value = common.try_float(i['ACT_VALUE']) or None,
            act_type = i['ACT_TYPE'],
            relation = i['RELATION'] or None, # what is relation??
            effect = i['ACTION_TYPE'] or None,
            tdl = i['TDL'],
            organism = i['ORGANISM'],
            comment = i['ACT_COMMENT'] if comments else None,
        )
        for i in interactions
        for uniprot in i['ACCESSION'].split('|')
        if not organism_latin or i['ORGANISM'] == organism_latin
    ]

    return result


def drugcentral_mapping(
        id_type: str,
        target_id_type: str,
    ) -> dict[str, set[str]]:
    """
    Identifier translation table from Drug Central.

    Available ID types: drugcentral, inn, cas, smiles, inchikey, inchi.

    Args
        id_type:
            The identifier type to be used as keys.
        target_id_type:
            The identifier type that will be collected into the values.

    Returns
        An identifier translation table.
    """

    drugs = drugcentral_drugs()

    result = collections.defaultdict(set)

    for d in drugs:

        the_id = getattr(d, id_type)
        target_id = getattr(d, target_id_type)

        if the_id and target_id:

            result[the_id].add(target_id)

    return dict(result)
