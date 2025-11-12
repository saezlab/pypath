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
Parse HMDB (Human Metabolome Database) data and emit Entity records.

This module downloads HMDB metabolite XML data and converts it into Entity
records using the declarative schema pattern from tabular_builder.
"""

from __future__ import annotations

from collections.abc import Generator

import lxml.etree as etree

from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    MoleculeAnnotationsCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceAnnotationCv,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    Annotations,
    Column,
    EntityBuilder,
    Identifiers,
)
from pypath.share.downloads import download_and_open


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing HMDB metadata.
    """
    yield Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.HMDB),
            Identifier(type=IdentifierNamespaceCv.NAME, value='Human Metabolome Database'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='37953221'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://hmdb.ca/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'The Human Metabolome Database (HMDB) is a comprehensive database '
                'containing detailed information about small molecule metabolites '
                'found in the human body. It includes chemical, clinical, and '
                'biochemical/molecular biology data, with over 220,000 metabolite '
                'entries including both water-soluble and lipid-soluble metabolites.'
            )),
        ],
    )


def _iterate_metabolites(opener, max_records: int | None = None) -> Generator[dict[str, str], None, None]:
    """
    Iterate through HMDB XML metabolites and yield dictionaries.

    This function handles the XML parsing and yields flat dictionaries
    suitable for processing with the tabular_builder pattern.

    Args:
        opener: The file opener result from download_and_open
        max_records: Maximum number of records to parse

    Yields:
        Dictionary representations of metabolites
    """
    if not opener or not opener.result:
        return

    # The zip file contains hmdb_metabolites.xml
    xml_file = None
    for filename, file_handle in opener.result.items():
        if 'metabolites.xml' in filename.lower():
            xml_file = file_handle
            break

    if not xml_file:
        return

    # XML namespace used in HMDB
    xmlns = '{http://www.hmdb.ca}'

    def get_text(elem: etree._Element, tag: str, parent: etree._Element = None) -> str:
        """Helper to safely extract text from an XML element."""
        parent = parent if parent is not None else elem
        child = parent.find(f'{xmlns}{tag}')
        return child.text.strip() if child is not None and child.text else ''

    def get_list(elem: etree._Element, container_tag: str, item_tag: str) -> str:
        """Helper to extract a list of text values as semicolon-delimited string."""
        container = elem.find(f'{xmlns}{container_tag}')
        if container is None:
            return ''
        items = [
            item.text.strip() for item in container.findall(f'{xmlns}{item_tag}')
            if item.text and item.text.strip()
        ]
        return ';'.join(items)

    # Use iterparse for memory-efficient parsing of large XML
    context = etree.iterparse(xml_file, events=('end',), tag=f'{xmlns}metabolite')

    count = 0
    for _, elem in context:
        if max_records is not None and count >= max_records:
            break

        # Extract PubMed IDs from general_references
        pubmed_ids = []
        general_refs = elem.find(f'{xmlns}general_references')
        if general_refs is not None:
            for ref in general_refs.findall(f'{xmlns}reference'):
                pubmed_id = ref.find(f'{xmlns}pubmed_id')
                if pubmed_id is not None and pubmed_id.text:
                    pubmed_ids.append(pubmed_id.text.strip())

        yield {
            'accession': get_text(elem, 'accession'),
            'name': get_text(elem, 'name'),
            'traditional_iupac': get_text(elem, 'traditional_iupac'),
            'iupac_name': get_text(elem, 'iupac_name'),
            'synonyms': get_list(elem, 'synonyms', 'synonym'),
            'inchikey': get_text(elem, 'inchikey'),
            'inchi': get_text(elem, 'inchi'),
            'smiles': get_text(elem, 'smiles'),
            'chebi_id': get_text(elem, 'chebi_id'),
            'pubchem_compound_id': get_text(elem, 'pubchem_compound_id'),
            'kegg_id': get_text(elem, 'kegg_id'),
            'drugbank_id': get_text(elem, 'drugbank_id'),
            'cas_registry_number': get_text(elem, 'cas_registry_number'),
            'description': get_text(elem, 'description'),
            'pubmed_ids': ';'.join(pubmed_ids) if pubmed_ids else '',
        }

        count += 1

        # Clear the element to free memory
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

    del context


def hmdb_metabolites(
    max_records: int | None = None,
    force_refresh: bool = False,
) -> Generator[Entity, None, None]:
    """
    Download and parse HMDB metabolite data as Entity records.

    This function downloads the HMDB metabolites XML file, converts each
    XML element to a dictionary, and uses a declarative schema to build
    Entity records.

    Args:
        max_records: Maximum number of records to parse. If None, parses all records.
        force_refresh: If True, force redownload of the data.

    Yields:
        Entity records with type SMALL_MOLECULE, representing metabolites
        with their identifiers (HMDB ID, InChI, SMILES, etc.), chemical
        properties (molecular weight, formula), and cross-references to other
        databases (ChEBI, PubChem, KEGG, DrugBank, CAS).
    """
    # Download HMDB metabolites XML
    url = 'https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip'
    opener = download_and_open(
        url,
        filename='hmdb_metabolites.zip',
        subfolder='hmdb',
        large=True,
        ext='zip',
        default_mode='rb',  # Binary mode required for lxml
        force=force_refresh,
    )

    # Define declarative schema for HMDB metabolites
    schema = EntityBuilder(
        entity_type=EntityTypeCv.SMALL_MOLECULE,
        identifiers=Identifiers(
            # Primary HMDB accession
            Column('accession', cv=IdentifierNamespaceCv.HMDB),
            # Primary name
            Column('name', cv=IdentifierNamespaceCv.NAME),
            # Traditional IUPAC name
            Column('traditional_iupac', cv=IdentifierNamespaceCv.IUPAC_TRADITIONAL_NAME),
            # IUPAC name
            Column('iupac_name', cv=IdentifierNamespaceCv.IUPAC_NAME),
            # Synonyms (semicolon-delimited)
            Column('synonyms', delimiter=';', cv=IdentifierNamespaceCv.SYNONYM),
            # Chemical structure identifiers
            Column('inchikey', cv=IdentifierNamespaceCv.STANDARD_INCHI_KEY),
            Column('inchi', cv=IdentifierNamespaceCv.STANDARD_INCHI),
            Column('smiles', cv=IdentifierNamespaceCv.SMILES),
            # Cross-references to other databases
            Column(
                'chebi_id',
                cv=IdentifierNamespaceCv.CHEBI,
                processing={
                    'extract_prefix': r'^(CHEBI:)?',
                    'extract_value': r'^(?:CHEBI:)?(.+)',
                },
            ),
            Column('pubchem_compound_id', cv=IdentifierNamespaceCv.PUBCHEM_COMPOUND),
            Column('kegg_id', cv=IdentifierNamespaceCv.KEGG_COMPOUND),
            Column('drugbank_id', cv=IdentifierNamespaceCv.DRUGBANK),
            Column('cas_registry_number', cv=IdentifierNamespaceCv.CAS),
        ),
        annotations=Annotations(
            # Source annotation
            # Description
            Column('description', cv=MoleculeAnnotationsCv.DESCRIPTION),
            # PubMed references (semicolon-delimited)
            Column('pubmed_ids', delimiter=';', cv=IdentifierNamespaceCv.PUBMED),
        ),
    )

    # Parse and yield entities
    if opener and opener.result:
        for row in _iterate_metabolites(opener, max_records):
            yield schema(row)
