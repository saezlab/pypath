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
import pypath.share.session as session
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy

_logger = session.Logger(name = 'drugcentral_input')
_log = _logger._log



def drugcentral_drugs() -> list[tuple]:
    """
    Drug names and structures from Drug Central.

    Returns:
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


def drug_central(
        organism: Union[str, int] = 'Homo sapiens',
    ) -> list[tuple]:
    """
    Retrieves drug-target interactions from Drug Central.

    Args:
        organism:
            Organism name or NCBI Taxonomy ID.

    Returns:
        List of drug-target relationships, represented as named tuples.
    """

    fields = (
        'drug',
        'target',
        'target_class',
        'target_accession',
        'gene',
        'act_value',
        'act_type',
        'action_type',
        'tdl',
        'organism',
    )

    url = urls.urls['drugcentral']['interactions']
    c = curl.Curl(url, large = True, silent = False)
    interactions = list(csv.DictReader(c.result, delimiter = '\t'))

    organism_latin = taxonomy.ensure_latin_name(organism)

    if not organism_latin:

        msg = f'Could not find latin name for organism: `{organism}`.'
        _log(msg)

        raise ValueError(msg)

    result = []



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
