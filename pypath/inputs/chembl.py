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

from typing import Literal

import json
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls


def chembl_targets() -> list[tuple]:
    """
    Retrieves targets data from ChEMBL.

    Returns:
        List of drug target records as named tuples.
    """

    fields_target = (
        'accession',
        'target_chembl_id',
    )

    ChemblTarget = collections.namedtuple(
        'ChemblTarget',
        fields_target,
        defaults = (None,) * len(fields_target),
    )

    tgt_lst = []
    page_dct = {}

    while True:

        if not page_dct:

            url = (
                f"{urls.urls['chembl']['url']}"
                f"{urls.urls['chembl']['target']}"
            )

        elif page_dct['page_meta']['next']:

            url = (
                f"{urls.urls['chembl']['url']}"
                f"{page_dct['page_meta']['next']}"
            )

        else:

            break

        c = curl.Curl(url, large=True, silent=False)
        fileobj = open(c.fileobj.name)
        page_dct = json.loads(fileobj.read())

        tgt_lst.extend(
            ChemblTarget(
                accession = (
                    tgt['target_components'][0]['accession']
                        if tgt['target_components'] else
                    None
                ),
                target_chembl_id = tgt['target_chembl_id'],
            )
            for tgt in page_dct['targets']
        )

    return tgt_lst


def chembl_assays() -> list[tuple] :
    """
    Retrieves assays data from ChEMBL.

    Returns:
        List of assay records as named tuples.
    """

    fields_assay = (
        'assay_chembl_id',
        'assay_organism',
        'assay_type',
        'confidence_score',
        'target_chembl_id',
    )

    ChemblAssay = collections.namedtuple(
        'ChemblAssay',
        fields_assay,
        defaults = (None,) * len(fields_assay),
    )

    assay_lst = []
    page_dct = {}

    while True:

        if not page_dct:

            url = (
                f"{urls.urls['chembl']['url']}"
                f"{urls.urls['chembl']['assay']}"
            )

        elif page_dct['page_meta']['next']:

            url = (
                f"{urls.urls['chembl']['url']}"
                f"{page_dct['page_meta']['next']}"
            )

        else:

            break

        c = curl.Curl(url, large=True, silent=False)
        fileobj = open(c.fileobj.name)
        page_dct = json.loads(fileobj.read())

        assay_lst.extend(
            ChemblAssay(
                assay_chembl_id = assy_attr['assay_chembl_id'],
                assay_organism = assy_attr['assay_organism'],
                assay_type = assy_attr['assay_type'],
                confidence_score = assy_attr['confidence_score'],
                target_chembl_id = assy_attr['target_chembl_id'],
            )
            for assy_attr in page_dct['assays']
        )

    return assay_lst


def chembl_molecules() -> list[tuple]:
    """
    Retrieves molecules data from ChEMBL.

    Returns:
        Molecule records as named tuples.
    """

    def _get(mol, key0, key1):
    
        molecule_properties = mol.get(f'molecule_{key0}', {})
        
        if molecule_properties:
        
            return molecule_properties.get(key1, None)
            
        else:
        
            return None


    fields_molecule = (
        'alogp',
        'canonical_smiles',
        'chirality',
        'full_mwt',
        'heavy_atoms',
        'std_inchi_key',
        'species',
        'type',
        'chembl',
        'parent_chembl',
        'prodrug',
        'std_inchi',
        'xrefs',
    )

    ChemblMolecule = collections.namedtuple(
        'ChemblMolecule',
        fields_molecule,
        defaults = (None,) * len(fields_molecule),
    )

    mol_lst = []
    page_dct = {}

    while True:

        if not page_dct:

            url = urls.urls['chembl']['url'] + urls.urls['chembl']['molecule']
            c = curl.Curl(url, large=True, silent=False)

        elif page_dct['page_meta']['next']:

            url = (
                f"{urls.urls['chembl']['url']}"
                f"{page_dct['page_meta']['next']}"
            )

        else:

            break

        c = curl.Curl(url, large=True, silent=False)
        fileobj = open(c.fileobj.name)
        page_dct = json.loads(fileobj.read())

        mol_lst.extend(
            ChemblMolecule(
                chirality = mol['chirality'],
                type = mol['molecule_type'],
                prodrug = mol['prodrug'],

                chembl = _get(mol, 'hierarchy', 'molecule_chembl_id'),
                parent_chembl = _get(mol, 'hierarchy', 'parent_chembl_id'),

                alogp = _get(mol, 'properties', 'alogp'),
                full_mwt = _get(mol, 'properties', 'full_mwt'),
                heavy_atoms = _get(mol, 'properties', 'heavy_atoms'),
                species = _get(mol, 'properties', 'molecular_species'),

                canonical_smiles = _get(mol, 'structures', 'canonical_smiles'),
                std_inchi_key = _get(mol, 'structures', 'standard_inchi_key'),
                std_inchi = _get(mol, 'structures', 'standard_inchi'),

                xrefs = (
                    [
                        {
                            'xref_id': rec['xref_id'],
                            'xref_src': rec['xref_src'],
                        }
                        for rec in mol['cross_references']
                    ]
                        if mol['cross_references'] else
                    None
                )
            )
            for mol in page_dct['molecules']
        )

    return mol_lst


def chembl_activities(
        #TODO: are these below all the allowed values?
        standard_relation: Literal['=', '>', '<', '>=', '<='],
        pchembl_value_none: bool = False,
    ) -> list[tuple] :
    """
    Retrieves activities data from ChEMBL.

    Args:
        pchembl_value_none:
            # TODO: it is allowed to be None or must be None?
            Whether the pchembl value should be none or not.
        standard_relation:
            Which standard relation in needed.

    Returns:
        List of activity records as named tuples. `standard_flag` and
        `standard_units` attributes are not included in the returned records.
        # TODO: then why the data_validity_comment is part of the records?
        Only records without `data_validity_comment` are returned.
    """

    fields_activity = (
        'assay_chembl',
        'data_validity_comment',
        'chembl',
        'pchembl',
        'standard_relation',
        'standard_value',
        'target_chembl',
    )

    ChemblActivity = collections.namedtuple(
        'ChemblActivity',
        fields_activity,
        defaults = (None,) * len(fields_activity),
    )

    activity_lst = []
    page_dct = {}

    while True:

        if not page_dct:


            url = (
                f"{urls.urls['chembl']['url']}"
                f"{urls.urls['chembl']['activity']}"
                f"&pchembl_value__isnull={str(pchembl_value_none).lower()}"
                f"&standard_relation__exact={standard_relation}"
            )

        elif page_dct['page_meta']['next']:

            url = (
                f"{urls.urls['chembl']['url']}"
                f"{page_dct['page_meta']['next']}"
            )

        else:

            break

        c = curl.Curl(url, large=True, silent=False)
        fileobj = open(c.fileobj.name)
        page_dct = json.loads(fileobj.read())


        activity_lst.extend(
            ChemblActivity(
                assay_chembl = act['assay_chembl_id'],
                data_validity_comment = act['data_validity_comment'],
                chembl = act['molecule_chembl_id'],
                pchembl = act['pchembl_value'],
                standard_relation = act['standard_relation'],
                standard_value = act['standard_value'],
                target_chembl = act['target_chembl_id'],
            )
            for act in page_dct['activities']
            if act['data_validity_comment'] is None
        )

    return activity_lst
