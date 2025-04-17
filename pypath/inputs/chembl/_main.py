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

from typing import Literal
from collections.abc import Generator

import json
import collections

from pypath.share import curl
from pypath.resources import urls

DATA = Literal[
    "target",
    "assay",
    "molecule",
    "activity",
    "document",
    "drug_indication",
    "mechanism",
]

ChemblTarget = collections.namedtuple(
    "ChemblTarget",
    [
        "chembl_id",
        "target_type",
        "preferred_name",
        "ncbi_taxa_id",
        "organism",
        "components",
        "num_components",
    ]
)

ChemblComponent = collections.namedtuple(
    "ChemblComponent",
    [
        "uniprot_accession",
        "component_type",
        "component_description",
        "component_id",
        "component_relationship",
        "component_number",
    ]
)

ChemblAssay = collections.namedtuple(
    "ChemblAssay",
    [
        "assay_chembl_id",
        "assay_type",
        "assay_type_description",
        "assay_category",
        "target_chembl_id",
        "organism",
        "tax_id",
        "tissue",
        "cell_type",
        "subcellular_fraction",
        "parameters",
        "source_id",
        "confidence_score",
        "confidence_description",
        "document_chembl_id",
        "description",
    ]
)

ChemblParam = collections.namedtuple(
    "ChemblParam",
    [
        "standard_type",
        "standard_value",
        "standard_units",
        "standard_relation",
    ]
)

ChemblActivity = collections.namedtuple(
    "ChemblActivity",
    [
        "activity_id",
        "action_type",
        "standard_relation",
        "standard_value",
        "standard_upper_value",
        "standard_units",
        "standard_type",
        "ligand_efficiency",
        "validity_comment",
        "potential_duplicate",
        "pchembl_value",
        "source_id",
        "molecule_chembl_id",
        "target_chembl_id",
        "target_taxa_id",
        "assay_id",
        "document_chembl_id",
    ]
)

ChemblDocument = collections.namedtuple(
    "ChemblDocument",
    [
        "document_chembl_id",
        'pubmed_id',
        "patent_id",
        "doc_type",
        "journal",
        "year",
        "doi",
    ]
)

def chembl_general(data_type: DATA,
                   max_pages: int) -> Generator[dict]:
    """
    Retrieves data from ChEMBL.

    This function will retrieve the specified data from ChEMBL
    and yield the data as a JSON object.

    Args:
        data_type (DATA): The type of data to retrieve.
        max_pages (int): The maximum number of pages to retrieve.

    Yields:
        dict: The JSON object of the retrieved data.
    """

    page_dict: dict = {}
    page_count = 0

    while True:

        if not page_dict:

            url_base = urls.urls['chembl']['url']
            url_path = urls.urls['chembl'][data_type]
            url = f"{url_base}{url_path}"

        elif page_dict['page_meta']['next']:

            url_base = urls.urls['chembl']['url']
            url_path = page_dict['page_meta']['next']
            url = f"{url_base}{url_path}"

        else:
            break
        
        c = curl.Curl(url, large=True, silent=False)
        with open(c.fileobj.name, mode="r", encoding='utf-8') as read_file:
            page_dict = json.load(read_file) # load the json object

            # remove unwanted page_meta
            for key in page_dict.keys():
                if key == 'page_meta':
                    continue
                else:
                    data_dict = page_dict[key]

                    # extract data from each page
                    for data in data_dict:
                        yield data
        
        # allows user to specify maximum number of pages
        page_count += 1
        if page_count >= max_pages:
            break

def get_targets(max_pages: int) -> Generator[ChemblTarget]:
    """
    Retrieves target data from ChEMBL.

    This generator function retrieves the target data from ChEMBL and
    yields the data as named tuples of the type `ChemblTarget`.

    The function uses the `chembl_general` function to retrieve the
    data from ChEMBL.

    Args:
        max_pages (int): The maximum number of pages to retrieve.
    Yields:
        ChemblTarget: The named tuple of the retrieved data.
    """

    targets= chembl_general(data_type="target", max_pages=max_pages)

    # loop through the pages and yield the target data
    for target in targets:

        # components
        components = tuple(target_components(target))
        yield ChemblTarget(
            chembl_id = target['target_chembl_id'],
            target_type = target['target_type'],
            preferred_name = target["pref_name"],
            ncbi_taxa_id = target["tax_id"],
            organism = target["organism"],
            components = components,
            num_components = len(components),
        )

def target_components(target: dict) -> Generator[ChemblComponent]:
    """
    Retrieves target component data from ChEMBL.
    """
    comp_count = 0
    for component in target['target_components']:
        comp_count += 1
        yield ChemblComponent(
            uniprot_accession = component['accession'],
            component_type = component['component_type'],
            component_description = component['component_description'],
            component_id = component['component_id'],
            component_relationship = component['relationship'],
            component_number = comp_count,
        )

def get_assays(max_pages: int) -> Generator[ChemblAssay]:
    """
    Retrieves assay data from ChEMBL.

    This generator function retrieves the assay data from ChEMBL and
    yields the data as named tuples of the type `ChemblAssay`.

    The function uses the `chembl_general` function to retrieve the
    data from ChEMBL.

    Args:
        max_pages (int): The maximum number of pages to retrieve.

    Yields:
        ChemblAssay: The named tuple of the retrieved data.
    """

    assays = chembl_general(data_type="assay", max_pages=max_pages)

    for assay in assays:
        
        # Get the assay parameters
        parameters = tuple(param_assay(assay['assay_parameters']))
        
        # Create the ChemblAssay named tuple
        yield ChemblAssay(
            assay_chembl_id = assay['assay_chembl_id'],
            assay_type = assay['assay_type'],
            assay_type_description = assay['assay_type_description'],
            assay_category = assay['assay_category'],
            target_chembl_id= assay['target_chembl_id'],
            organism = assay['assay_organism'],
            tax_id = assay['assay_tax_id'],
            tissue = assay['assay_tissue'],
            cell_type = assay['assay_cell_type'],
            subcellular_fraction = assay['assay_subcellular_fraction'],
            parameters = parameters,
            source_id = assay['src_id'],
            confidence_score = assay['confidence_score'],
            confidence_description = assay['confidence_description'],
            document_chembl_id = assay['document_chembl_id'],
            description = assay['description'],
        )

def param_assay(parameters: dict) -> Generator[ChemblParam]:
    """
    Retrieves assay parameter data from ChEMBL.

    Args:
        parameters (dict): The dictionary of parameters.

    Yields:
        ChemblParam: The named tuple of the retrieved parameter data.
    """
    if parameters:
        yield from (ChemblParam
        (
            standard_type=parameter["standard_type"],
            standard_value=parameter["standard_value"],
            standard_units=parameter["standard_units"],
            standard_relation=parameter["standard_relation"]
        )
        for parameter in parameters
    )

def get_activity(max_pages: int) -> Generator[ChemblActivity]:
    """
    Retrieves activity data from Chembl.
    This generator function retrieves the assay data from ChEMBL and
    yields the data as named tuples of the type `ChemblActivity.

    The function uses the `chembl_general` function to retrieve the
    data from ChEMBL.

    Args:
        max_pages (int): The maximum number of pages to retrieve.

    Yields:
        ChemblActivity: The named tuple of the retrieved data.   
    """
    activities = chembl_general(data_type="activity", max_pages=max_pages)

    yield from (ChemblActivity
        (
            activity_id=activity["activity_id"],
            action_type=activity["action_type"],
            standard_relation=activity["standard_relation"],
            standard_value=activity["standard_value"],
            standard_upper_value=activity['standard_upper_value'],
            standard_units=activity["standard_units"],
            standard_type=activity["standard_type"],
            ligand_efficiency=activity["ligand_efficiency"],
            validity_comment=activity["data_validity_comment"],
            potential_duplicate=activity["potential_duplicate"],
            pchembl_value=activity["pchembl_value"],
            source_id=activity["src_id"],
            molecule_chembl_id=activity["molecule_chembl_id"],
            target_chembl_id=activity["target_chembl_id"],
            target_taxa_id=activity["target_tax_id"],
            assay_id=activity["assay_chembl_id"],
            document_chembl_id=activity["document_chembl_id"],
        )
        for activity in activities
    )

def get_documents(max_pages: int) -> Generator[ChemblDocument]:
    """
    Retrieves the Chembl document information.

    This function is a wrapper around the `chembl_general` function.
    It retrieves the document information from ChEMBL and
    yields the data as named tuples of the type `ChemblDocument`.

    Args:
        max_pages (int): The maximum number of pages to retrieve.

    Yields:
        ChemblDocument: The named tuple of the retrieved data.
    """
    documents = chembl_general(data_type="document", max_pages=max_pages)

    yield from (ChemblDocument
        (
            document_chembl_id=document["document_chembl_id"],
            pubmed_id=document["pubmed_id"],
            patent_id=document["patent_id"],
            doc_type=document["doc_type"],
            journal=document["journal"],
            year=document["year"],
            doi=document["doi"],
        )
        for document in documents
    )

def chembl_molecules() -> list[tuple]:
    """
    Retrieves molecules data from ChEMBL.

    Returns
        Molecule records as named tuples.
    """

    def _get(mol, key0, key1):
    
        molecule_properties = mol.get(f'molecule_{key0}', {})
        
        if molecule_properties:
        
            return molecule_properties.get(key1, None)
            
        else:
        
            return None


    fields_molecule = (
        'name',
        'alogp',
        'canonical_smiles',
        'chirality',
        'full_mwt',
        'heavy_atoms',
        'species',
        'qed_weighted',
        'type',
        'structure_type',
        'chembl',
        'parent_chembl',
        'prodrug',
        'std_inchi_key',
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
        fileobj = open(c.fileobj.name, encoding='utf-8')
        page_dct = json.loads(fileobj.read())

        mol_lst.extend(
            ChemblMolecule(
                name = mol['pref_name'],
                chirality = mol['chirality'],
                type = mol['molecule_type'],
                prodrug = mol['prodrug'],
                structure_type = mol['structure_type'],

                chembl = _get(mol, 'hierarchy', 'molecule_chembl_id'),
                parent_chembl = _get(mol, 'hierarchy', 'parent_chembl_id'),

                alogp = _get(mol, 'properties', 'alogp'),
                full_mwt = _get(mol, 'properties', 'full_mwt'),
                heavy_atoms = _get(mol, 'properties', 'heavy_atoms'),
                species = _get(mol, 'properties', 'molecular_species'),
                qed_weighted = _get(mol, 'properties', 'qed_weighted'),

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
        standard_relation: Literal['=', '>', '<', '>=', '<=', None] = None,
        pchembl_value_none: bool = False,
        limit: int = 1000
    ) -> list[tuple] :
    """
    Retrieves activities data from ChEMBL.

    Args
        pchembl_value_none:
            # TODO: it is allowed to be None or must be None?
            Whether the pchembl value should be none or not.
        standard_relation:
            Which standard relation in needed.

    Returns
        List of activity records as named tuples. 
        `standard_units` attribute is not included in the returned records.
        # TODO: then why the data_validity_comment is part of the records?
        Only records without `data_validity_comment` are returned.
    """

    fields_activity = (
        'assay_chembl',
        'data_validity_comment',
        'chembl',
        'pchembl',
        'standard_flag',
        'standard_relation',
        'standard_value',
        'standard_type',
        'target_chembl',
        'document'
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
                #f"&pchembl_value__isnull={str(pchembl_value_none).lower()}"
                #f"&standard_relation__exact={standard_relation}"
            )

        elif page_dct['page_meta']['next'] and len(activity_lst) < limit:

            url = (
                f"{urls.urls['chembl']['url']}"
                f"{page_dct['page_meta']['next']}"
            )

        else:

            break

        c = curl.Curl(url, large=True, silent=False)
        fileobj = open(c.fileobj.name, encoding='utf-8')
        page_dct = json.loads(fileobj.read())


        activity_lst.extend(
            ChemblActivity(
                assay_chembl = act['assay_chembl_id'],
                data_validity_comment = act['data_validity_comment'],
                chembl = act['molecule_chembl_id'],
                pchembl = act['pchembl_value'],
                standard_flag = True if act['standard_flag'] == 1 else False,
                standard_relation = act['standard_relation'],
                standard_value = act['standard_value'],
                standard_type = act['standard_type'],
                target_chembl = act['target_chembl_id'],
                document = act['document_chembl_id'],
            )
            for act in page_dct['activities']
            if act['data_validity_comment'] is None
        )

    return activity_lst


def chembl_documents() -> dict[str, str] :
    """
    Retrieves ChEMBL document ID to PubMed ID conversion.

    Returns
        Dictionary of ChEMBL document IDs as keys and PubMed IDs as values.   
    """

    page_dct = {}
    document_dict = {}

    while True:
        if not page_dct:
            url = (
                f"{urls.urls['chembl']['url']}"
                f"{urls.urls['chembl']['document']}"
            )

        elif page_dct['page_meta']['next']:

            url = (
                    f"{urls.urls['chembl']['url']}"
                    f"{page_dct['page_meta']['next']}"
                )
            
        else:
                
            break

        c = curl.Curl(url, large=True, silent=False)
        fileobj = open(c.fileobj.name, encoding='utf-8')
        page_dct = json.loads(fileobj.read())

        for doc in page_dct['documents']:
            if doc['pubmed_id']:
                document_dict[doc['document_chembl_id']]= doc['pubmed_id']
    
    return document_dict


def chembl_drug_indications(
    max_phase_threshold: int = 0,
    ) -> list[tuple]:
    """
    Retrieves drug indications data from ChEMBL.

    Args
        max_phase_threshold:
            The threshold for maximum phase of the drug 
            for which the indication is valid.
    Returns
        List of drug indications as namedtuples.
    """

    fields_indication = (
        'efo_id',
        'efo_term',
        'max_phase',
        'mesh_heading',
        'mesh_id',
        'molecule_chembl',
    )

    ChemblIndication = collections.namedtuple(
        'ChemblIndication',
        fields_indication,
        defaults = (None,) * len(fields_indication),
    )

    indication_lst = []
    page_dct = {}

    while True:

        if not page_dct:

            url = (
                f"{urls.urls['chembl']['url']}"
                f"{urls.urls['chembl']['drug_indication']}"
            )

        elif page_dct['page_meta']['next']:
            url = (
                f"{urls.urls['chembl']['url']}"
                f"{page_dct['page_meta']['next']}"
            )
        
        else:
            break
        
        c = curl.Curl(url, large=True, silent=False)
        fileobj = open(c.fileobj.name, encoding='utf-8')
        page_dct = json.loads(fileobj.read())

        indication_lst.extend(
            ChemblIndication(
                efo_id = ind['efo_id'],
                efo_term = ind['efo_term'],
                max_phase = float(ind['max_phase_for_ind']),
                mesh_heading = ind['mesh_heading'],
                mesh_id = ind['mesh_id'],
                molecule_chembl = ind['molecule_chembl_id'],
            )
            for ind in page_dct['drug_indications']
            if float(ind['max_phase_for_ind']) > max_phase_threshold and max_phase_threshold != 0 \
                or max_phase_threshold == 0
        )
    
    return indication_lst


def chembl_mechanisms() -> list[tuple]:
    """
    Retrieves mechanism data from ChEMBL.

    Returns
        List of mechanisms as namedtuples.
    """

    fields_mechanism = (
        'action_type',
        'direct_interaction',
        'disease_efficacy',
        'mechanism_of_action',
        'chembl',
        'target_chembl',
    )

    ChemblMechanism= collections.namedtuple(
        'ChemblMechanism',
        fields_mechanism,
        defaults = (None,) * len(fields_mechanism),
    )

    mechanism_lst = []
    page_dct = {}

    while True:

        if not page_dct:
            url = (
                f"{urls.urls['chembl']['url']}"
                f"{urls.urls['chembl']['mechanism']}"
            )

        elif page_dct['page_meta']['next']:
            url = (
                f"{urls.urls['chembl']['url']}"
                f"{page_dct['page_meta']['next']}"
            )

        else:
            break

        c = curl.Curl(url, large=True, silent=False)
        fileobj = open(c.fileobj.name, encoding='utf-8')
        page_dct = json.loads(fileobj.read())

        mechanism_lst.extend(
            ChemblMechanism(
                action_type = mech['action_type'],
                direct_interaction = True if mech['direct_interaction'] == 1 else False,
                disease_efficacy = True if mech['disease_efficacy'] == 1 else False,
                mechanism_of_action = mech['mechanism_of_action'],
                chembl = mech['molecule_chembl_id'],
                target_chembl = mech['target_chembl_id'],
            )
            for mech in page_dct['mechanisms']
        )
    
    return mechanism_lst