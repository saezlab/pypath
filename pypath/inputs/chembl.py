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

import json
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls


def chembl_targets() -> list[tuple]:
    """
    Retrieves targets data from ChEMBL.

    Returns:
        namedtuple.
    """

    fields_target = ('accession','target_chembl_id')

    Target = collections.namedtuple('Target', fields_target,defaults = ("None",) * len(fields_target))

    trgtlst = []

    flag = 0

    while True:

        if flag == 0:

            url = urls.urls['chembl']['url'] + urls.urls['chembl']['target']
            c = curl.Curl(url, large=True, silent=False)
            flag = 1

        else:

            if lst['page_meta']['next']:

                url = urls.urls['chembl']['url'] + lst['page_meta']['next']
                c = curl.Curl(url, large=True, silent=False)

            else:

                break

        fileObject = open(c.fileobj.name)
        lst = json.loads(fileObject.read())

        for trgt_attr in lst['targets']:

            if trgt_attr['target_components']:

                trgtlst.append(
                    Target(
                        accession = trgt_attr['target_components'][0]['accession'],
                        target_chembl_id = trgt_attr['target_chembl_id'],
                        )
                    )

            else:

                trgtlst.append(
                    Target(
                        target_chembl_id = trgt_attr['target_chembl_id'],
                        )
                    )

    return trgtlst

def chembl_assays() -> List[tuple] :
    """
    Retrieves assays data from ChEMBL.

    Returns:
        namedtuple.
    """

    fields_assay = ('assay_chembl_id','assay_organism','assay_type','confidence_score','target_chembl_id')

    Assay = collections.namedtuple('Assay', fields_assay,defaults = ("None",) * len(fields_assay))

    assylst = []

    flag = 0

    while True:

        if flag == 0:

            url = urls.urls['chembl']['url'] + urls.urls['chembl']['assay']
            c = curl.Curl(url, large=True, silent=False)
            flag = 1

        else:

            if lst['page_meta']['next']:

                url = urls.urls['chembl']['url'] + lst['page_meta']['next']
                c = curl.Curl(url, large=True, silent=False)

            else:

                break

        fileObject = open(c.fileobj.name)
        lst = json.loads(fileObject.read())

        for assy_attr in lst['assays']:

            assylst.append(
                Assay(
                    assay_chembl_id = assy_attr['assay_chembl_id'],
                    assay_organism = assy_attr['assay_organism'],
                    assay_type = assy_attr['assay_type'],
                    confidence_score = assy_attr['confidence_score'],
                    target_chembl_id = assy_attr['target_chembl_id'],
                    )
                )

    return assylst

def chembl_molecules() -> List[tuple] :
    """
    Retrieves molecules data from ChEMBL.

    Returns:
        namedtuple.
    """

    fields_molecule = ('alogp','conanicle_smiles','chirality','full_mwt','heavy_atoms','standard_inchi_key','molecular_species',
                        'molecul_type','molecule_chembl_id','parent_chembl_id','prodrug','standard_inchi', 'xrefs')

    Molecule = collections.namedtuple('Molecule', fields_molecule,defaults = ("None",) * len(fields_molecule))

    mlcllst = []

    flag = 0

    while True:

        if flag == 0:

            url = urls.urls['chembl']['url'] + urls.urls['chembl']['molecule']
            c = curl.Curl(url, large=True, silent=False)
            flag = 1

        else:

            if lst['page_meta']['next']:

                url = urls.urls['chembl']['url'] + lst['page_meta']['next']
                c = curl.Curl(url, large=True, silent=False)

            else:

                break

        fileObject = open(c.fileobj.name)
        lst = json.loads(fileObject.read())

        for mlcl_attr in lst['molecules']:

            xrefs = []
            mlcllst.append(
                Molecule(
                    chirality = mlcl_attr['chirality'],
                    molecul_type = mlcl_attr['molecule_type'],
                    prodrug = mlcl_attr['prodrug'],
                    )
                )

            if mlcl_attr['molecule_hierarchy'] != None:
                mlcllst[-1] = mlcllst[-1]._replace(
                    molecule_chembl_id = mlcl_attr['molecule_hierarchy']['molecule_chembl_id'],
                    parent_chembl_id = mlcl_attr['molecule_hierarchy']['parent_chembl_id'],
                )

            if mlcl_attr['molecule_properties'] != None:
                mlcllst[-1] = mlcllst[-1]._replace(
                    alogp = mlcl_attr['molecule_properties']['alogp'],
                    full_mwt = mlcl_attr['molecule_properties']['full_mwt'],
                    heavy_atoms = mlcl_attr['molecule_properties']['heavy_atoms'],
                    molecular_species = mlcl_attr['molecule_properties']['molecular_species'],
                )

            if mlcl_attr['molecule_structures'] != None:
                mlcllst[-1] = mlcllst[-1]._replace(
                    conanicle_smiles = mlcl_attr['molecule_structures']['canonical_smiles'],
                    standard_inchi_key = mlcl_attr['molecule_structures']['standard_inchi_key'],
                    standard_inchi = mlcl_attr['molecule_structures']['standard_inchi'],
                )

            if mlcl_attr['cross_references'] != None:

                for rec in mlcl_attr['cross_references']:

                    xrefs.append({'xref_id' : rec['xref_id'], 'xref_src': rec['xref_src']})

                mlcllst[-1] = mlcllst[-1]._replace(
                    xrefs = xrefs
                )


    return mlcllst

def chembl_activities(
        pchembl_value_none: bool = False,
        standard_relation: bool = '=',
    ) -> List[tuple] :
    """
    Retrieves activities data from ChEMBL.

    Args:
        pchembl_value_none (bool): Whether the pchembl value should be none or not.
        standard_relation (str): Which standard relation in needed.

    Returns:
        namedtuple.
            standard_flag and standard_units attributes are not included in the returned namedtuple.
            Only records returned are the ones where data_validity_comment is none.
    """

    fields_activity = ('assay_chembl_id','data_validity_comment','molecule_chembl_id','pchembl_value',
                        'standard_relation','standard_value','target_chembl_id')

    Activity = collections.namedtuple('Activity', fields_activity,defaults = ("None",) * len(fields_activity))

    actvtylst = []

    flag = 0

    while True:

        if flag == 0:

            if pchembl_value_none == True:

                url = urls.urls['chembl']['url'] + urls.urls['chembl']['activity']+'&pchembl_value__isnull=true'

            else:

                url = urls.urls['chembl']['url'] + urls.urls['chembl']['activity']+'&pchembl_value__isnull=false'

            url = url + '&standard_relation__exact='+standard_relation
            c = curl.Curl(url, large=True, silent=False)
            flag = 1

        else:

            if lst['page_meta']['next']:

                url = urls.urls['chembl']['url'] + lst['page_meta']['next']
                c = curl.Curl(url, large=True, silent=False)

            else:

                break

        fileObject = open(c.fileobj.name)
        lst = json.loads(fileObject.read())


        for actvty_attr in lst['activities']:

            if actvty_attr['data_validity_comment'] == None:

                actvtylst.append(
                    Activity(
                        assay_chembl_id = actvty_attr['assay_chembl_id'],
                        data_validity_comment = actvty_attr['data_validity_comment'],
                        molecule_chembl_id = actvty_attr['molecule_chembl_id'],
                        pchembl_value = actvty_attr['pchembl_value'],
                        standard_relation = actvty_attr['standard_relation'],
                        standard_value = actvty_attr['standard_value'],
                        target_chembl_id = actvty_attr['target_chembl_id'],
                        )
                    )


    return actvtylst
