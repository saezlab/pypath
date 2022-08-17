#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Tennur Kılıç
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import csv
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common

def drug_central(
        organism: str = "Homo sapiens",
        SMILES: bool = False,
        InChI: bool = False,
        CAS_RN: bool = False,
    ) -> list[tuple]:
    """
    Retrieves drug-target interactions datasets from Drug Central.

    Args:
        organism (str): Which organism to use for processing.
        SMILES (bool): Whether to include SMILES structures from Drug Central.
        InChI (bool): Whether to include InChI formats and InChI keys from Drug Central.
        CAS_RN (bool): Whether to include CAS Registry Number from Drug Central.

    Returns:
        namedtuple.
    """

    fields = (
        'DRUG_NAME',
        'TARGET_NAME',
        'TARGET_CLASS',
        'TARGET_ACCESSION',
        'GENE',
        'ACT_VALUE',
        'ACT_TYPE',
        'ACTION_TYPE',
        'TDL',
        'ORGANISM',
        'SMILES',
        'InChI',
        'InChIKey',
        'CAS_RN',
    )

    url = urls.urls['drugcentral']['interactions']
    c = curl.Curl(url, large = True, silent = False)
    interactions = list(csv.DictReader(c.result, delimiter = '\t'))

    interactions = common.unique_list(interactions)

    result = []

    if SMILES == True or InChI == True or CAS_RN == True:

        url = urls.urls['drugcentral']['SMILES_InChI']
        c = curl.Curl(url, large = True, silent = False)
        structures = list(csv.DictReader(c.result, delimiter = '\t'))

        temp_struct = []

        for rec in structures:

            if rec not in temp_struct:

                temp_struct.append(rec)

        structures = temp_struct

        if SMILES == True and InChI == True and CAS_RN == True:

            InteractionsandStructures = collections.namedtuple('InteractionsandStructures', fields[0:], defaults = (None,) * len(fields))

        elif SMILES == True and InChI == True and CAS_RN == False:

            fields = fields[0:13]
            InteractionsandStructures = collections.namedtuple('InteractionsandStructures', fields, defaults = (None,) * len(fields))

        elif SMILES == True and InChI == False and CAS_RN == True:

            fields = fields[0:11] + fields[13:]
            InteractionsandStructures = collections.namedtuple('InteractionsandStructures', fields, defaults = (None,) * len(fields))

        elif SMILES == True and InChI == False and CAS_RN == False:

            fields = fields[0:11]
            InteractionsandStructures = collections.namedtuple('InteractionsandStructures', fields, defaults = (None,) * len(fields))

        elif SMILES == False and InChI == True and CAS_RN == True:

            fields = fields[0:10] + fields[11:]
            InteractionsandStructures = collections.namedtuple('InteractionsandStructures', fields, defaults = (None,) * len(fields))

        elif SMILES == False and InChI == False and CAS_RN == True:

            fields = fields[13:]
            InteractionsandStructures = collections.namedtuple('InteractionsandStructures', fields, defaults = (None,) * len(fields))

        elif SMILES == False and InChI == True and CAS_RN == False:

            fields = fields[0:10] + fields[11:13]
            InteractionsandStructures = collections.namedtuple('InteractionsandStructures', fields, defaults = (None,) * len(fields))

        for inter_attr in interactions:

            if organism == inter_attr['ORGANISM']:

                result.append(
                    InteractionsandStructures(
                        DRUG_NAME = inter_attr['DRUG_NAME'],
                        TARGET_NAME = inter_attr['TARGET_NAME'],
                        TARGET_CLASS = inter_attr['TARGET_CLASS'],
                        TARGET_ACCESSION = inter_attr['ACCESSION'],
                        GENE = inter_attr['GENE'],
                        ACT_VALUE = inter_attr['ACT_VALUE'],
                        ACT_TYPE = inter_attr['ACT_TYPE'],
                        ACTION_TYPE = inter_attr['ACTION_TYPE'],
                        TDL = inter_attr['TDL'],
                        ORGANISM = inter_attr['ORGANISM'],
                        )
                    )

                for struct_attr in structures:

                    if inter_attr['STRUCT_ID'] == struct_attr['ID']:

                        if SMILES == True and InChI == True and CAS_RN == True:

                            result[-1] = result[-1]._replace(
                                SMILES = struct_attr['SMILES'],
                                InChI = struct_attr['InChI'],
                                InChIKey = struct_attr['InChIKey'],
                                CAS_RN = struct_attr['CAS_RN'],
                            )

                        elif SMILES == True and InChI == True and CAS_RN == False:

                            result[-1] = result[-1]._replace(
                                SMILES = struct_attr['SMILES'],
                                InChI = struct_attr['InChI'],
                                InChIKey = struct_attr['InChIKey'],
                            )

                        elif SMILES == True and InChI == False and CAS_RN == True:

                            result[-1] = result[-1]._replace(
                                SMILES = struct_attr['SMILES'],
                                CAS_RN = struct_attr['CAS_RN'],
                            )

                        elif SMILES == True and InChI == False and CAS_RN == False:

                            result[-1] = result[-1]._replace(
                                SMILES = struct_attr['SMILES'],
                            )

                        elif SMILES == False and InChI == True and CAS_RN == True:

                            result[-1] = result[-1]._replace(
                                InChI = struct_attr['InChI'],
                                InChIKey = struct_attr['InChIKey'],
                                CAS_RN = struct_attr['CAS_RN'],
                            )

                        elif SMILES == False and InChI == False and CAS_RN == True:

                            result[-1] = result[-1]._replace(
                                CAS_RN = struct_attr['CAS_RN'],
                            )

                        elif SMILES == False and InChI == True and CAS_RN == False:

                            result[-1] = result[-1]._replace(
                                InChI = struct_attr['InChI'],
                                InChIKey = struct_attr['InChIKey'],
                            )

    else:

        DrugTargetInteractions = collections.namedtuple('DrugTargetInteractions', fields[0:10])

        for inter_attr in interactions:

            if organism == inter_attr['ORGANISM']:

                result.append(
                    DrugTargetInteractions(
                        DRUG_NAME = inter_attr['DRUG_NAME'],
                        TARGET_NAME = inter_attr['TARGET_NAME'],
                        TARGET_CLASS = inter_attr['TARGET_CLASS'],
                        TARGET_ACCESSION = inter_attr['ACCESSION'],
                        GENE = inter_attr['GENE'],
                        ACT_VALUE = inter_attr['ACT_VALUE'],
                        ACT_TYPE = inter_attr['ACT_TYPE'],
                        ACTION_TYPE = inter_attr['ACTION_TYPE'],
                        TDL = inter_attr['TDL'],
                        ORGANISM = inter_attr['ORGANISM'],
                    )
                )

    return result
