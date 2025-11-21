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
import re

from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    InteractionParameterCv,
    AffinityUnitCv,
    CurationCv,
    ResourceAnnotationCv,
    ResourceCv,
)
from pypath.share.downloads import download_and_open
from ..internals.tabular_builder import (
    AnnotationsBuilder,
    Column,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
    Map,
)


# Patterns to identify different identifier types in the ligand name field
_LIGAND_NAME_PATTERNS = {
    IdentifierNamespaceCv.CHEMBL_COMPOUND: re.compile(r'^CHEMBL\d+$'),
    IdentifierNamespaceCv.ZINC: re.compile(r'^ZINC\d+$'),
    IdentifierNamespaceCv.CAS: re.compile(r'^\d{2,7}-\d{2}-\d$'),
    IdentifierNamespaceCv.CHEBI: re.compile(r'^CHEBI[:\s]?\d+$', re.IGNORECASE),
    IdentifierNamespaceCv.PUBCHEM_COMPOUND: re.compile(r'^CID[:\s]?\d+$', re.IGNORECASE),
    IdentifierNamespaceCv.KEGG_COMPOUND: re.compile(r'^C\d{5}$'),
}

# Pattern for patent references (we'll skip these as they're not useful identifiers)
_PATENT_PATTERN = re.compile(r'^(US|EP|WO|JP|CN|KR|AU|CA)\d+')


def _parse_ligand_name(row: dict) -> list[tuple[IdentifierNamespaceCv, str]]:
    """
    Parse the BindingDB Ligand Name field which contains multiple identifiers
    separated by '::'.

    Returns a list of (namespace, value) tuples for each identified component.
    """
    raw_value = row.get('BindingDB Ligand Name', '')
    if not raw_value:
        return []

    results: list[tuple[IdentifierNamespaceCv, str]] = []
    parts = raw_value.split('::')

    for part in parts:
        part = part.strip()
        if not part:
            continue

        # Skip patent references
        if _PATENT_PATTERN.match(part):
            continue

        # Try to match known identifier patterns
        matched = False
        for namespace, pattern in _LIGAND_NAME_PATTERNS.items():
            if pattern.match(part):
                # Extract numeric part for some identifiers
                if namespace == IdentifierNamespaceCv.CHEBI:
                    # Normalize CHEBI:123 -> 123
                    value = re.sub(r'^CHEBI[:\s]?', '', part, flags=re.IGNORECASE)
                elif namespace == IdentifierNamespaceCv.PUBCHEM_COMPOUND:
                    # Normalize CID123 -> 123
                    value = re.sub(r'^CID[:\s]?', '', part, flags=re.IGNORECASE)
                else:
                    value = part
                results.append((namespace, value))
                matched = True
                break

        # If no pattern matched, treat as a name
        if not matched:
            results.append((IdentifierNamespaceCv.NAME, part))

    return results


def _ligand_names(row: dict) -> list[str]:
    """Extract all NAME-type values from the ligand name field."""
    return [v for ns, v in _parse_ligand_name(row) if ns == IdentifierNamespaceCv.NAME]


def _ligand_chembl(row: dict) -> list[str]:
    """Extract CHEMBL IDs from the ligand name field."""
    return [v for ns, v in _parse_ligand_name(row) if ns == IdentifierNamespaceCv.CHEMBL_COMPOUND]


def _ligand_zinc(row: dict) -> list[str]:
    """Extract ZINC IDs from the ligand name field."""
    return [v for ns, v in _parse_ligand_name(row) if ns == IdentifierNamespaceCv.ZINC]


def _ligand_cas(row: dict) -> list[str]:
    """Extract CAS numbers from the ligand name field."""
    return [v for ns, v in _parse_ligand_name(row) if ns == IdentifierNamespaceCv.CAS]


def _ligand_chebi(row: dict) -> list[str]:
    """Extract CHEBI IDs from the ligand name field."""
    return [v for ns, v in _parse_ligand_name(row) if ns == IdentifierNamespaceCv.CHEBI]


def _ligand_pubchem_cid(row: dict) -> list[str]:
    """Extract PubChem CIDs from the ligand name field."""
    return [v for ns, v in _parse_ligand_name(row) if ns == IdentifierNamespaceCv.PUBCHEM_COMPOUND]


def _ligand_kegg(row: dict) -> list[str]:
    """Extract KEGG compound IDs from the ligand name field."""
    return [v for ns, v in _parse_ligand_name(row) if ns == IdentifierNamespaceCv.KEGG_COMPOUND]


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing BindingDB metadata.
    """
    yield Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.BINDINGDB),
            Identifier(type=IdentifierNamespaceCv.NAME, value='BindingDB'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='26481362'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://www.bindingdb.org/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'BindingDB is a public, web-accessible database of measured binding '
                'affinities, focusing chiefly on the interactions of proteins considered '
                'to be drug-targets with small, drug-like molecules. It contains binding '
                'data for over 2 million protein-ligand complexes with experimental '
                'measurements including Ki, Kd, IC50, and EC50 values.'
            )),
        ],
    )


def bindingdb_interactions(
    dataset: str = 'All',
    max_lines: int | None = None,
    force_refresh: bool = False,
) -> Generator[Entity, None, None]:
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
    tax_column = Column('Target Source Organism According to Curator or DataSource')
    tax_value = Map(col=tax_column, extract=[r'(\d+)'])

    # Define the schema mapping
    schema = EntityBuilder(
        entity_type=EntityTypeCv.INTERACTION,
        identifiers=IdentifiersBuilder(
            # Use BindingDB Reactant_set_id as the interaction identifier
            CV(term=IdentifierNamespaceCv.BINDINGDB, value=Column('BindingDB Reactant_set_id')),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            # Binding affinity measurements
            CV(
                term=InteractionParameterCv.KI,
                value=Column('Ki (nM)'),
                unit=AffinityUnitCv.NANOMOLAR,
            ),
            CV(
                term=InteractionParameterCv.KD,
                value=Column('Kd (nM)'),
                unit=AffinityUnitCv.NANOMOLAR,
            ),
            CV(
                term=InteractionParameterCv.IC50,
                value=Column('IC50 (nM)'),
                unit=AffinityUnitCv.NANOMOLAR,
            ),
            CV(
                term=InteractionParameterCv.EC50,
                value=Column('EC50 (nM)'),
                unit=AffinityUnitCv.NANOMOLAR,
            ),
            CV(
                term=InteractionParameterCv.KON,
                value=Column('kon (M-1-s-1)'),
                unit=AffinityUnitCv.PER_MOLAR_PER_SECOND,
            ),
            CV(
                term=InteractionParameterCv.KOFF,
                value=Column('koff (s-1)'),
                unit=AffinityUnitCv.PER_SECOND,
            ),
            # Experimental conditions
            CV(term=InteractionParameterCv.PH, value=Column('pH')),
            CV(
                term=InteractionParameterCv.TEMPERATURE_CELSIUS,
                value=Column('Temp (C)'),
                unit=AffinityUnitCv.DEGREE_CELSIUS,
            ),
            # References and metadata
            CV(term=IdentifierNamespaceCv.PUBMED, value=Column('PMID')),
            CV(term=IdentifierNamespaceCv.DOI, value=Column('Article DOI')),
            CV(term=IdentifierNamespaceCv.PATENT_NUMBER, value=Column('Patent Number')),
            CV(term=CurationCv.COMMENT, value=Column('Curation/DataSource')),
        ),
        membership=MembershipBuilder(
            # Interactor A: Ligand (small molecule)
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.SMALL_MOLECULE,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.BINDINGDB, value=Column('BindingDB MonomerID')),
                        # Parse ligand names from the compound field (contains :: separated values)
                        CV(term=IdentifierNamespaceCv.NAME, value=_ligand_names),
                        # Extract cross-references embedded in the ligand name field
                        CV(term=IdentifierNamespaceCv.CHEMBL_COMPOUND, value=_ligand_chembl),
                        CV(term=IdentifierNamespaceCv.ZINC, value=_ligand_zinc),
                        CV(term=IdentifierNamespaceCv.CAS, value=_ligand_cas),
                        CV(term=IdentifierNamespaceCv.CHEBI, value=_ligand_chebi),
                        CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=_ligand_pubchem_cid),
                        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=_ligand_kegg),
                        # Standard identifiers from dedicated columns
                        CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=Column('Ligand InChI Key')),
                        CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=Column('Ligand InChI')),
                        CV(term=IdentifierNamespaceCv.SMILES, value=Column('Ligand SMILES')),
                        CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=Column('PubChem CID')),
                        CV(term=IdentifierNamespaceCv.PUBCHEM, value=Column('PubChem SID')),
                        CV(term=IdentifierNamespaceCv.CHEBI, value=Column('ChEBI ID of Ligand')),
                        CV(term=IdentifierNamespaceCv.CHEMBL_COMPOUND, value=Column('ChEMBL ID of Ligand')),
                        CV(term=IdentifierNamespaceCv.DRUGBANK, value=Column('DrugBank ID of Ligand')),
                        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=Column('KEGG ID of Ligand')),
                        CV(term=IdentifierNamespaceCv.ZINC, value=Column('ZINC ID of Ligand')),
                        CV(term=IdentifierNamespaceCv.PDB, value=Column('Ligand HET ID in PDB')),
                    ),
                ),
            ),
            # Interactor B: Target (protein)
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.NAME, value=Column('Target Name')),
                        CV(term=IdentifierNamespaceCv.UNIPROT, value=Column('UniProt (SwissProt) Primary ID of Target Chain')),
                        CV(term=IdentifierNamespaceCv.NAME, value=Column('UniProt (SwissProt) Recommended Name of Target Chain')),
                        CV(term=IdentifierNamespaceCv.UNIPROT_TREMBL, value=Column('UniProt (TrEMBL) Primary ID of Target Chain')),
                        CV(term=IdentifierNamespaceCv.NAME, value=Column('UniProt (TrEMBL) Submitted Name of Target Chain')),
                    ),
                    annotations=AnnotationsBuilder(
                        CV(
                            term=IdentifierNamespaceCv.NCBI_TAX_ID,
                            value=tax_value,
                        ),
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
