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

ChemblMolecule = collections.namedtuple(
    "ChemblMolecule",
    [
        "molecule_chembl_id",
        "preferred_name",
        "molecule_type",
        "structure_type",
        "chirality",
        "biotherapeutic",
        "inorganic_flag",
        "natural_flag",
        "polymer_flag",
        "helm_notation",
        "molecule_properties",
        "structure",
    ]
)

ChemblMolProps = collections.namedtuple(
    "ChemblMolProps",
    [
        "mol_formula",
        "full_mwt",
        "monoisotopic_mwt",
        "molecular_species",
        "logp",
        "logd",
        "alogp",
    ]
)

ChemblMolStruct = collections.namedtuple(
    "ChemblMolStruct",
    [
        "canonical_smiles",
        "inchi",
        "inchi_key",
    ]
)

ChemblMechanism = collections.namedtuple(
    "ChemblMechanism",
    [
        "action_type",
        "molecule_chembl_id",
        "target_chembl_id",
        "mechanism_id",
        "drug_rec_id",
        "direct_interaction",
        "variant_sequence",
        "molecular_mechanism",
        "mechanism_refs",
    ]
)
def chembl_general(data_type: DATA,
                   max_pages: int | None = None) -> Generator[dict]:
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

            # remove unwanted page meta data
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
        if max_pages and page_count >= max_pages:
            break

def get_targets(max_pages: int | None = None) -> Generator[ChemblTarget]:
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

def get_assays(max_pages: int | None = None) -> Generator[ChemblAssay]:
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
            standard_type = parameter["standard_type"],
            standard_value = parameter["standard_value"],
            standard_units = parameter["standard_units"],
            standard_relation = parameter["standard_relation"]
        )
        for parameter in parameters
    )

def get_activity(max_pages: int | None = None) -> Generator[ChemblActivity]:
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
            activity_id = activity["activity_id"],
            action_type = activity["action_type"],
            standard_relation = activity["standard_relation"],
            standard_value = activity["standard_value"],
            standard_upper_value = activity['standard_upper_value'],
            standard_units = activity["standard_units"],
            standard_type = activity["standard_type"],
            ligand_efficiency = activity["ligand_efficiency"],
            validity_comment = activity["data_validity_comment"],
            potential_duplicate = activity["potential_duplicate"],
            pchembl_value = activity["pchembl_value"],
            source_id = activity["src_id"],
            molecule_chembl_id = activity["molecule_chembl_id"],
            target_chembl_id = activity["target_chembl_id"],
            target_taxa_id = activity["target_tax_id"],
            assay_id = activity["assay_chembl_id"],
            document_chembl_id = activity["document_chembl_id"],
        )
        for activity in activities
    )

def get_documents(max_pages: int | None = None) -> Generator[ChemblDocument]:
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
            document_chembl_id = document["document_chembl_id"],
            pubmed_id = document["pubmed_id"],
            patent_id = document["patent_id"],
            doc_type = document["doc_type"],
            journal = document["journal"],
            year = document["year"],
            doi = document["doi"],
        )
        for document in documents
    )

def get_molecules(max_pages: int | None = None) -> Generator[ChemblMolecule]:
    """
    Retrieves molecule information from ChEMBL
    """

    molecules = chembl_general(data_type="molecule", max_pages=max_pages)

    for molecule in molecules:

        # Get the molecule properties
        molecule_properties = molecule_props(molecule["molecule_properties"])

        # Get the molecule structures
        structures = molecule_strucs(molecule["molecule_structures"])

        yield ChemblMolecule(
            molecule_chembl_id = molecule['molecule_chembl_id'],
            preferred_name = molecule['pref_name'],
            molecule_type = molecule['molecule_type'],
            structure_type = molecule['structure_type'],
            chirality = molecule['chirality'],
            biotherapeutic = molecule['biotherapeutic'],
            inorganic_flag = molecule['inorganic_flag'],
            natural_flag = molecule['natural_product'],
            polymer_flag = molecule['polymer_flag'],
            helm_notation = molecule['helm_notation'],
            molecule_properties = molecule_properties,
            structure = structures,
        )
def molecule_props(properties: dict) -> ChemblMolProps | None:
    """
    Retrieves molecule properties from ChEMBL.

    Args:
        properties (dict): The dictionary of molecule properties.

    Returns:
        ChemblMolProps: The named tuple of the molecule properties.
    """
    return ChemblMolProps(
        mol_formula=properties.get("full_molformula"),
        full_mwt=properties.get("full_mwt"),
        monoisotopic_mwt=properties.get("mw_monoisotopic"),
        molecular_species=properties.get("molecular_species"),
        logp=properties.get("cx_logp"),
        logd=properties.get("cx_logd"),
        alogp=properties.get("alogp"),
    ) if properties else None

def molecule_strucs(structure: dict) -> ChemblMolStruct | None:
    """
    Retrieves molecule structure from ChEMBL.

    Args:
        structure (dict): The dictionary of molecule structure.

    Returns:
        ChemblMolStruct: The named tuple of the molecule structure.
    """
    return ChemblMolStruct(
        canonical_smiles=structure.get("canonical_smiles"),
        inchi=structure.get("standard_inchi"),
        inchi_key=structure.get("standard_inchi_key"),
    ) if structure else None

def get_mechanisms(max_pages: int | None = None) -> Generator[ChemblMechanism]:
    """
    Retrieves mechanism data from ChEMBL.
    """

    mechanisms = chembl_general(data_type="mechanism", max_pages=max_pages)

    for mechanism in mechanisms:

        # Get the mechanism references
        references = [ref["ref_url"] for ref in mechanism["mechanism_refs"]]

        yield ChemblMechanism(
            action_type = mechanism["action_type"],
            molecule_chembl_id = mechanism["molecule_chembl_id"],
            target_chembl_id = mechanism["target_chembl_id"],
            mechanism_id = mechanism["mec_id"],
            drug_rec_id = mechanism["record_id"],
            direct_interaction = mechanism["direct_interaction"],
            variant_sequence = mechanism["variant_sequence"],
            molecular_mechanism = mechanism["molecular_mechanism"],
            mechanism_refs = references,
        )

# TODO: check if drug indications are required
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