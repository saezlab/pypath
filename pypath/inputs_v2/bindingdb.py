#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
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

"""
Parse BindingDB data and emit Entity records.

This module converts BindingDB binding affinity data into Entity records using
the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv

from pypath.internals.silver_schema import Entity as SilverEntity, Resource
from pypath.internals.cv_terms import EntityTypeCv, IdentifierNamespaceCv, LicenseCV, UpdateCategoryCV, InteractionParameterCv, CurationCv
from pypath.share.downloads import download_and_open
from ..internals.tabular_builder import (
    Annotations,
    Column,
    Entity,
    Identifiers,
    Member,
    Members,
)


def get_resource() -> Resource:
    """
    Define the resource metadata.

    Returns:
        Resource object containing BindingDB metadata.
    """
    return Resource(
        id='bindingdb',
        name='BindingDB',
        license=LicenseCV.CC_BY_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        publication='PMID:26481362',  # BindingDB in 2015 paper
        url='https://www.bindingdb.org/',
        description=(
            'BindingDB is a public, web-accessible database of measured binding '
            'affinities, focusing chiefly on the interactions of proteins considered '
            'to be drug-targets with small, drug-like molecules. It contains binding '
            'data for over 2 million protein-ligand complexes with experimental '
            'measurements including Ki, Kd, IC50, and EC50 values.'
        ),
    )


def bindingdb_interactions(
    dataset: str = 'All',
    max_lines: int | None = None,
    force_refresh: bool = False,
) -> Generator[SilverEntity, None, None]:
    """
    Download and parse BindingDB protein-ligand binding data as Entity records.

    Args:
        dataset: BindingDB dataset to read ('All', 'Articles', 'ChEMBL', 'PubChem', 'Patents', etc.)
        max_lines: Maximum number of lines to read. If None, reads the entire file.
        force_refresh: If True, force redownload

    Yields:
        Entity records with type INTERACTION, representing protein-ligand binding interactions
    """
    # Download and open file using download manager
    # For zip files, opener.result is a dict with extracted file handles
    url = f'https://bindingdb.org/rwd/bind/downloads/BindingDB_{dataset}_202506_tsv.zip'
    opener = download_and_open(
        url,
        filename=f'{dataset}_table.zip',
        subfolder='bindingdb',
        large=True,
        encoding='utf-8',
        ext='zip',
    )
    # Processing pattern to extract taxonomy ID
    tax_processing = {'extract_value': r'(\d+)'}

    # Define the schema mapping
    schema = Entity(
        entity_type=EntityTypeCv.INTERACTION,
        identifiers=Identifiers(
            # Use BindingDB Reactant_set_id as the interaction identifier
            Column('BindingDB Reactant_set_id', cv=IdentifierNamespaceCv.BINDINGDB),
        ),
        annotations=Annotations(
            # Binding affinity measurements
            Column('Ki (nM)', cv=InteractionParameterCv.KI),
            Column('Kd (nM)', cv=InteractionParameterCv.KD),
            Column('IC50 (nM)', cv=InteractionParameterCv.IC50),
            Column('EC50 (nM)', cv=InteractionParameterCv.EC50),
            Column('kon (M-1-s-1)', cv=InteractionParameterCv.KON),
            Column('koff (s-1)', cv=InteractionParameterCv.KOFF),
            # Experimental conditions
            Column('pH', cv=InteractionParameterCv.PH),
            Column('Temp (C)', cv=InteractionParameterCv.TEMPERATURE_CELSIUS),
            # References and metadata
            Column('PMID', cv=IdentifierNamespaceCv.PUBMED),
            Column('Article DOI', cv=IdentifierNamespaceCv.DOI),
            Column('Patent Number', cv=IdentifierNamespaceCv.PATENT_NUMBER),
            Column('Curation/DataSource', cv=CurationCv.COMMENT),
        ),
        members=Members(
            # Interactor A: Ligand (small molecule)
            Member(
                entity=Entity(
                    entity_type=EntityTypeCv.SMALL_MOLECULE,
                    identifiers=Identifiers(
                        Column('BindingDB MonomerID', cv=IdentifierNamespaceCv.BINDINGDB),
                        Column('BindingDB Ligand Name', cv=IdentifierNamespaceCv.NAME),
                        Column('Ligand InChI Key', cv=IdentifierNamespaceCv.STANDARD_INCHI_KEY),
                        Column('Ligand InChI', cv=IdentifierNamespaceCv.STANDARD_INCHI),
                        Column('Ligand SMILES', cv=IdentifierNamespaceCv.SMILES),
                        Column('PubChem CID', cv=IdentifierNamespaceCv.PUBCHEM_COMPOUND),
                        Column('PubChem SID', cv=IdentifierNamespaceCv.PUBCHEM),
                        Column('ChEBI ID of Ligand', cv=IdentifierNamespaceCv.CHEBI),
                        Column('ChEMBL ID of Ligand', cv=IdentifierNamespaceCv.CHEMBL_COMPOUND),
                        Column('DrugBank ID of Ligand', cv=IdentifierNamespaceCv.DRUGBANK),
                        Column('KEGG ID of Ligand', cv=IdentifierNamespaceCv.KEGG_COMPOUND),
                        Column('ZINC ID of Ligand', cv=IdentifierNamespaceCv.NAME),
                        Column('Ligand HET ID in PDB', cv=IdentifierNamespaceCv.PDB),
                    ),
                ),
            ),
            # Interactor B: Target (protein)
            Member(
                entity=Entity(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=Identifiers(
                        Column('Target Name', cv=IdentifierNamespaceCv.NAME),
                        Column('UniProt (SwissProt) Primary ID of Target Chain', cv=IdentifierNamespaceCv.NAME),
                        Column('UniProt (SwissProt) Recommended Name of Target Chain', cv=IdentifierNamespaceCv.NAME),
                        Column('UniProt (TrEMBL) Primary ID of Target Chain', cv=IdentifierNamespaceCv.UNIPROT),
                        Column('UniProt (TrEMBL) Submitted Name of Target Chain', cv=IdentifierNamespaceCv.NAME),
                    ),
                    annotations=Annotations(
                        Column('Target Source Organism According to Curator or DataSource', processing=tax_processing, cv=IdentifierNamespaceCv.NCBI_TAX_ID),
                        Column('BindingDB Target Chain Sequence', cv=IdentifierNamespaceCv.NAME),
                        Column('Number of Protein Chains in Target (>1 implies a multichain complex)', cv=IdentifierNamespaceCv.NAME),
                    ),
                ),
            ),
        ),
    )

    # Parse and yield entities
    # For zip files, opener.result is a dict with filename -> file handle
    if opener and opener.result:
        for file_handle in opener.result.values():
            # BindingDB has 638 columns but only 50 unique column names.
            # Columns 38-637 contain 50 repetitions of 12 target chain columns.
            # We keep only columns 0-48 to capture 94.54%+ to keep it simple for now.
            # TODO: In the future we should decide whether to capture the complex and then have complex - small molecule interaction. 
            # or whether we create a hyperedge to connect all chains to the ligand.
            
            header_line = file_handle.readline().strip()
            header = header_line.split('\t')

            # Keep first 49 columns (indices 0-48)
            columns_to_keep = min(49, len(header))
            filtered_header = header[:columns_to_keep]

            # Create a generator that yields filtered lines (as strings, not lists)
            def filtered_rows():
                for line in file_handle:
                    # Split the line, keep only the desired columns, and rejoin
                    columns = line.strip().split('\t')
                    yield '\t'.join(columns[:columns_to_keep])

            # Create DictReader with filtered columns
            reader = csv.DictReader(filtered_rows(), fieldnames=filtered_header, delimiter='\t')

            if max_lines:
                for i, row in enumerate(reader):
                    if i >= max_lines:
                        break
                    yield schema(row)
            else:
                for row in reader:
                    yield schema(row)
            break  # Only process the first file in the zip