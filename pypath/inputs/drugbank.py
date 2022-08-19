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
import base64

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.settings as settings

_logger = session.Logger(name = 'drugbank')
_log = _logger._log


def _drugbank_download(user: str, passwd: str, *args, **kwargs):

    defaults = {
        'large': True,
        'silent': False,
        'compr': 'zip',
    }

    defaults.update(kwargs)

    auth_str = base64.b64encode(f"{user}:{passwd}".encode())

    defaults['req_headers'] = [
        f'Authorization: Basic {auth.decode()}',
        settings.get('user_agent'),
    ]

    return curl.Curl(*args, **defaults)


def drugbank_raw_interactions(
        user: str,
        passwd: str,
        pharma_active: bool = False,
    ) -> list[tuple] :
    """
    Retrieves protein identifiers from Drugbank.

    Args:
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.
        pharma_active:
            Only pharmacologically active relations.

    Returns:
        List of drug-protein relations.
    """

    csv_name = 'pharmacologically_active.csv' if pharma_active else 'all.csv'

    fields = (
        'drugbank_id',
        'uniprot_id',
        'relation',
    )

    DrugbankProtein = collections.namedtuple(
        'DrugbankProtein',
        fields,
        defaults = (None,) * len(fields),
    )

    result = []

    for rel in ('carrier', 'enzyme', 'target', 'transporter'):

        url = urls.urls['drugbank'][f'drug_{rel}_identifiers']

        c = _drugbank_download(
            user = user,
            passwd = passwd,
            files_needed = (csv_name,),
        )

        _ = next(c.result[csv_name])

        for l in c.result[csv_name]:

            drugs, uniprot = l.strip().split(',')

            result.extend(
                DrugbankProtein(
                    drugbank_id = drug,
                    uniprot_id = uniprot,
                    relation = rel,
                )
                for drug in drugs
            )

    return result


def drugbank(
        user: str,
        passwd: str,
        addprotid: bool = True,
        pharma_active: bool = False,
    ) -> list[tuple] :
    """
    Retrieves structures, external links and protein identifiers from Drugbank.

    Args:
        user (str): E-mail address for login to DrugBank.
        passwd (str): Password for login to DrugBank.
        addprotid (bool): Wheter to include protein identifiers from DrugBank.
        pharma_active (bool): Wheter to include pharmacologically active identifiers.

    Returns:
        namedtuple.
    """

    fields = ('DrugBank_ID','Name','CAS_Number','Drug_Groups','InChIKey','InChI','SMILES','Formula',
                'KEGG_Compound_ID','KEGG_Drug_ID','PubChem_Compound_ID','PubChem_Substance_ID','ChEBI_ID',
                'ChEMBL_ID','Drug_Type','PharmGKB_ID','HET_ID','Target_UniProt_ID','Transporter_UniProt_ID',
                'Enzym_UniProt_ID','Carrier_UniProt_ID')

    credentials = {'user': user, 'passwd': passwd}

    auth_str = base64.b64encode(
        ('%s:%s' % (credentials['user'], credentials['passwd'])).encode()
    ).decode()

    decoded = 'Basic %s' % auth_str

    req_hdrs = ['Authorization: %s' % decoded]
    req_hdrs.extend([settings.get('user_agent')])

    url = urls.urls['drugbank']['all_structures']
    c = curl.Curl(
        url,
        large = True,
        silent = False,
        req_headers = req_hdrs,
        compr = 'zip',
        files_needed = ('structure links.csv',),
    )

    structure_links = list(
        csv.DictReader(
            c.result['structure links.csv'],
            delimiter = ',',
        )
    )

    url = urls.urls['drugbank']['all_drug']
    c = curl.Curl(
        url,
        large = True,
        silent = False,
        req_headers = req_hdrs,
        compr = 'zip',
        files_needed = ('drug links.csv',),
    )

    drug_links = list(
        csv.DictReader(
            c.result['drug links.csv'],
            delimiter = ',',
        )
    )

    if addprotid:

        Combine = collections.namedtuple('Combine', fields,defaults = ("",) * len(fields))

    else:
        Combine = collections.namedtuple('Combine', fields[:17],defaults = ("",) * len(fields[:17]))

    result = []

    for struct_attr in structure_links:

        for drug_attr in drug_links:

            if struct_attr['DrugBank ID'] == drug_attr['DrugBank ID']:

                result.append(
                    Combine(
                        DrugBank_ID = struct_attr['DrugBank ID'],
                        Name = struct_attr['Name'],
                        CAS_Number = struct_attr['CAS Number'],
                        Drug_Groups = struct_attr['Drug Groups'],
                        InChIKey = struct_attr['InChIKey'],
                        InChI = struct_attr['InChI'],
                        SMILES = struct_attr['SMILES'],
                        Formula = struct_attr['Formula'],
                        KEGG_Compound_ID = struct_attr['KEGG Compound ID'],
                        KEGG_Drug_ID = struct_attr['KEGG Drug ID'],
                        PubChem_Compound_ID = struct_attr['PubChem Compound ID'],
                        PubChem_Substance_ID = struct_attr['PubChem Substance ID'],
                        ChEBI_ID = struct_attr['ChEBI ID'],
                        ChEMBL_ID = struct_attr['ChEMBL ID'],
                        Drug_Type = drug_attr['Drug Type'],
                        PharmGKB_ID = drug_attr['PharmGKB ID'],
                        HET_ID = drug_attr['HET ID'],
                    )
                )

    if addprotid:

        identifiers_list = add_prot_id(user, passwd, pharma_active)
        index = 0

        for res_attr in result:

            for iden_attr in identifiers_list:

                if res_attr.DrugBank_ID == iden_attr.DrugBank_ID:

                    result[index] = result[index]._replace(
                        Target_UniProt_ID = iden_attr.Target_UniProt_ID,
                        Transporter_UniProt_ID = iden_attr.Transporter_UniProt_ID,
                        Enzym_UniProt_ID = iden_attr.Enzym_UniProt_ID,
                        Carrier_UniProt_ID = iden_attr.Carrier_UniProt_ID,
                    )

                    break

            index += 1

    return result
