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

from __future__ import annotations

from collections import namedtuple
import csv
import re
import pandas as pd
import bs4

import pypath.share.curl as curl
import pypath.resources.urls as urls


def ddinter_n_drugs() -> int:
    """
    Retrieves number of unique drug records from DDInter2 CSV files.
    """
    
    drug_ids = set()
    
    # DDInter2 has files with codes A, B, D, H, L, P, R, V
    for code in ['A', 'B', 'D', 'H', 'L', 'P', 'R', 'V']:
        url = urls.urls['ddinter']['interaction'] % code
        c = curl.Curl(url, silent=True, large=True)
        
        reader = csv.DictReader(c.result)
        for row in reader:
            drug_ids.add(row['DDInterID_A'])
            drug_ids.add(row['DDInterID_B'])
    
    return len(drug_ids)


def ddinter_identifiers(drug: str) -> list[tuple]:
    """
    DrugBank, ChEMBL and PubChem identifiers of a drug in DDInter2.

    Args:
        drug:
            DDInter drug identifier.

    Returns:
        list of mapping namedtuple (single)
    """
    
    url = urls.urls['ddinter']['mapping'] % drug
    c = curl.Curl(url, silent=True, large=True)
    
    soup = bs4.BeautifulSoup(c.fileobj, 'html.parser')
    
    # Initialize result dictionary
    mapping_dict = {'ddinter': drug}
    
    # Find all links in the "Useful Links" section
    links = soup.find_all('a', href=True)
    
    for link in links:
        href = link.get('href', '')
        
        # Extract DrugBank ID
        if 'drugbank.com/drugs/' in href:
            drugbank_match = re.search(r'/drugs/(DB\d+)', href)
            if drugbank_match:
                mapping_dict['drugbank'] = drugbank_match.group(1)
        
        # Extract PubChem SID
        elif 'pubchem.ncbi.nlm.nih.gov' in href and 'sid=' in href:
            pubchem_match = re.search(r'sid=(\d+)', href)
            if pubchem_match:
                mapping_dict['pubchem'] = pubchem_match.group(1)
        
        # Extract ChEMBL ID  
        elif 'chembl' in href and 'compound_report_card' in href:
            chembl_match = re.search(r'compound_report_card/(CHEMBL\d+)', href)
            if chembl_match:
                mapping_dict['chembl'] = chembl_match.group(1)
    
    # Create namedtuple
    fields = ('ddinter', 'drugbank', 'chembl', 'pubchem')
    record = namedtuple('DdinterIdentifiers', fields, defaults=(None,) * len(fields))
    
    result = [record(
        ddinter=mapping_dict.get('ddinter'),
        drugbank=mapping_dict.get('drugbank'),
        chembl=mapping_dict.get('chembl'), 
        pubchem=mapping_dict.get('pubchem')
    )]
    
    return result


def ddinter_mappings(return_df: bool = False) -> list[tuple] | pd.DataFrame:
    """
    DrugBank, ChEMBL and PubChem identifiers of all drugs in DDInter2.

    Args:
        return_df:
            Return a pandas data frame.

    Returns:
        list of mapping namedtuples
    """

    result = []
    drug_ids = set()
    
    # Collect all unique drug IDs from CSV files
    for code in ['A', 'B', 'D', 'H', 'L', 'P', 'R', 'V']:
        url = urls.urls['ddinter']['interaction'] % code
        c = curl.Curl(url, silent=True, large=True)
        
        reader = csv.DictReader(c.result)
        for row in reader:
            drug_ids.add(row['DDInterID_A'])
            drug_ids.add(row['DDInterID_B'])
    
    # Get identifiers for each drug
    for drug_id in drug_ids:
        result.extend(ddinter_identifiers(drug_id))

    return pd.DataFrame(result) if return_df else result


def ddinter_drug_interactions(
        drug: str,
        return_df: bool = False,
    ) -> list[tuple] | pd.DataFrame:
    """
    Interactions of one single drug from the DDInter2 database.

    Args:
        drug:
            A DDInter drug identifier.
        return_df:
            Return a pandas data frame.

    Returns:
        Drug-drug interaction tuples with drug ids, drug names, interaction
        level and actions.
    """

    result = set()
    record = namedtuple(
        'DdinterInteraction',
        (
            'drug1_id',
            'drug1_name', 
            'drug2_id',
            'drug2_name',
            'level',
        ),
    )

    # Search through all CSV files for interactions involving this drug
    for code in ['A', 'B', 'D', 'H', 'L', 'P', 'R', 'V']:
        url = urls.urls['ddinter']['interaction'] % code
        c = curl.Curl(url, silent=True, large=True)
        
        reader = csv.DictReader(c.result)
        for row in reader:
            if row['DDInterID_A'] == drug or row['DDInterID_B'] == drug:
                interaction = record(
                    drug1_id = row['DDInterID_A'],
                    drug1_name = row['Drug_A'],
                    drug2_id = row['DDInterID_B'],
                    drug2_name = row['Drug_B'],
                    level = row['Level'],
                )
                result.add(interaction)

    return pd.DataFrame(result) if return_df else list(result)


def ddinter_interactions(return_df: bool = False) -> list[tuple] | pd.DataFrame:
    """
    All drug-drug interactions from the DDInter2 database.

    Args:
        return_df:
            Return a pandas data frame.

    Returns:
        list of interaction namedtuples with drug ids, drug names, and interaction levels
    """
    
    result = set()
    record = namedtuple(
        'DdinterInteraction',
        (
            'drug1_id',
            'drug1_name',
            'drug2_id', 
            'drug2_name',
            'level',
        ),
    )

    # Read all interactions from all CSV files
    for code in ['A', 'B', 'D', 'H', 'L', 'P', 'R', 'V']:
        url = urls.urls['ddinter']['interaction'] % code
        c = curl.Curl(url, silent=True, large=True)
        
        reader = csv.DictReader(c.result)
        for row in reader:
            interaction = record(
                drug1_id = row['DDInterID_A'],
                drug1_name = row['Drug_A'],
                drug2_id = row['DDInterID_B'],
                drug2_name = row['Drug_B'],
                level = row['Level'],
            )
            result.add(interaction)

    return pd.DataFrame(result) if return_df else list(result)
