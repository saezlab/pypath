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
Parse SIGNOR interaction data and emit Entity records.

This module converts SIGNOR MITAB (causalTab) format interactions into
Entity records with type INTERACTION.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.internals.silver_schema import Entity, Identifier, Annotation, Membership
from pypath.internals.cv_terms import IdentifierNamespaceCv
from ..mitab import (
    mitab_field_uniprot,
    mitab_parse_identifiers,
    mitab_parse_mi_terms,
)

__all__ = [
    'signor_interactions',
]

def download_interactions(organism: int = 9606) -> Generator[str]:
    """
    Download SIGNOR interaction data in causalTab (MITAB) format.

    Args:
        organism: NCBI taxonomy ID (9606 for human, 10090 for mouse, 10116 for rat)

    Yields:
        Lines from the causalTab file
    """
    if isinstance(organism, int):
        if organism in taxonomy.taxids:
            _organism = taxonomy.taxids[organism]
        else:
            raise ValueError(f'Unknown organism: {organism}')
    else:
        _organism = organism

    if _organism not in {'human', 'rat', 'mouse'}:
        raise ValueError(f'Organism {_organism} not supported by SIGNOR')

    url = urls.urls['signor']['all_url_new']

    # Download file with POST form data
    file_path = dm.download(
        url,
        filename=f'signor_{_organism}_causalTab.txt',
        subfolder='signor',
        query={
            'organism': _organism,
            'format': 'causalTab',
            'submit': 'Download',
        },
        post=True,
    )

    # Read and yield lines from the downloaded file
    if file_path:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:  # Skip empty lines
                    yield line


def _parse_entity_type(mitab_type: str) -> str:
    """Parse MITAB interactor type field and return MI accession or default."""
    if not mitab_type or mitab_type == '-':
        return 'MI:0326'  # Default to protein

    # Try to extract MI ID from format like psi-mi:"MI:0326"(protein)
    terms = mitab_parse_mi_terms(mitab_type)
    if terms and terms[0].get('mi_id'):
        return terms[0]['mi_id']

    return 'MI:0326'  # Default to protein


def signor_interactions(organism: int = 9606) -> Generator[Entity]:
    """
    Download and parse SIGNOR interactions as Entity records.

    Args:
        organism: NCBI taxonomy ID (9606 for human, 10090 for mouse, 10116 for rat)

    Yields:
        Entity records with type INTERACTION, containing participant entities as members
    """
    lines = download_interactions(organism=organism)

    # Skip header line
    header_skipped = False

    for line in lines:
        # Skip header
        if not header_skipped:
            header_skipped = True
            continue

        if not line or line.startswith('#'):
            continue

        # Split MITAB line (tab-delimited)
        fields = line.split('\t')

        # MITAB 2.8 format has 46 columns
        # Pad with empty strings if we have fewer
        while len(fields) < 46:
            fields.append('-')

        # Extract key fields (using 0-based indexing)
        id_a = fields[0]  # Column 1: ID A
        id_b = fields[1]  # Column 2: ID B
        alt_ids_a = fields[2]  # Column 3: Alt IDs A
        alt_ids_b = fields[3]  # Column 4: Alt IDs B
        detection_methods = fields[6]  # Column 7
        pmids = fields[8]  # Column 9
        interaction_types = fields[11]  # Column 12
        source_dbs = fields[12]  # Column 13
        interaction_ids = fields[13]  # Column 14
        biological_role_a = fields[16]  # Column 17
        biological_role_b = fields[17]  # Column 18
        experimental_role_a = fields[18]  # Column 19
        experimental_role_b = fields[19]  # Column 20
        interactor_type_a = fields[20]  # Column 21
        interactor_type_b = fields[21]  # Column 22
        causal_regulatory_mechanism = fields[44]  # Column 45
        causal_statement = fields[45]  # Column 46

        # Build interaction identifiers
        identifiers = []

        # Extract SIGNOR interaction ID
        interaction_id_list = mitab_parse_identifiers(interaction_ids)
        for id_info in interaction_id_list:
            if id_info['database'].lower() in ('signor', 'signor-interaction'):
                identifiers.append(
                    Identifier(
                        type=IdentifierNamespaceCv.SIGNOR,
                        value=id_info['id']
                    )
                )

        # If no identifiers found, skip this interaction
        if not identifiers:
            continue

        # Build interaction annotations
        annotations = []

        # Add interaction types (use MI accessions directly)
        interaction_type_terms = mitab_parse_mi_terms(interaction_types)
        for term in interaction_type_terms:
            mi_id = term.get('mi_id')
            if mi_id:
                annotations.append(Annotation(term=mi_id))

        # Add detection methods (use MI accessions directly)
        detection_method_terms = mitab_parse_mi_terms(detection_methods)
        for term in detection_method_terms:
            mi_id = term.get('mi_id')
            if mi_id:
                annotations.append(Annotation(term=mi_id))

        # Add causal statement (use MI accessions directly)
        causal_stmt_terms = mitab_parse_mi_terms(causal_statement)
        for term in causal_stmt_terms:
            mi_id = term.get('mi_id')
            if mi_id:
                annotations.append(Annotation(term=mi_id))

        # Add causal mechanism (use MI accessions directly)
        causal_mech_terms = mitab_parse_mi_terms(causal_regulatory_mechanism)
        for term in causal_mech_terms:
            mi_id = term.get('mi_id')
            if mi_id:
                annotations.append(Annotation(term=mi_id))

        # Add PubMed references
        if pmids and pmids != '-':
            for pmid_entry in pmids.split('|'):
                if ':' in pmid_entry:
                    db, pmid_value = pmid_entry.split(':', 1)
                    if db.lower() == 'pubmed':
                        annotations.append(
                            Annotation(
                                term='MI:0446',  # PUBMED
                                value=pmid_value.replace('"', '')
                            )
                        )

        # Build participant entities as members
        members = []

        # Participant A
        uniprot_a = mitab_field_uniprot(id_a) or mitab_field_uniprot(alt_ids_a)
        if uniprot_a:
            entity_type_a = _parse_entity_type(interactor_type_a)

            # Build annotations for participant A
            participant_a_annotations = []

            # Biological role (use MI accessions directly)
            bio_role_terms_a = mitab_parse_mi_terms(biological_role_a)
            for term in bio_role_terms_a:
                mi_id = term.get('mi_id')
                if mi_id:
                    participant_a_annotations.append(Annotation(term=mi_id))

            # Experimental role (use MI accessions directly)
            exp_role_terms_a = mitab_parse_mi_terms(experimental_role_a)
            for term in exp_role_terms_a:
                mi_id = term.get('mi_id')
                if mi_id:
                    participant_a_annotations.append(Annotation(term=mi_id))

            member_a = Entity(
                source='signor',
                type=entity_type_a,
                identifiers=[
                    Identifier(
                        type=IdentifierNamespaceCv.UNIPROT,
                        value=uniprot_a
                    )
                ],
            )

            members.append(
                Membership(
                    member=member_a,
                    annotations=participant_a_annotations if participant_a_annotations else None
                )
            )

        # Participant B
        uniprot_b = mitab_field_uniprot(id_b) or mitab_field_uniprot(alt_ids_b)
        if uniprot_b:
            entity_type_b = _parse_entity_type(interactor_type_b)

            # Build annotations for participant B
            participant_b_annotations = []

            # Biological role (use MI accessions directly)
            bio_role_terms_b = mitab_parse_mi_terms(biological_role_b)
            for term in bio_role_terms_b:
                mi_id = term.get('mi_id')
                if mi_id:
                    participant_b_annotations.append(Annotation(term=mi_id))

            # Experimental role (use MI accessions directly)
            exp_role_terms_b = mitab_parse_mi_terms(experimental_role_b)
            for term in exp_role_terms_b:
                mi_id = term.get('mi_id')
                if mi_id:
                    participant_b_annotations.append(Annotation(term=mi_id))

            member_b = Entity(
                source='signor',
                type=entity_type_b,
                identifiers=[
                    Identifier(
                        type=IdentifierNamespaceCv.UNIPROT,
                        value=uniprot_b
                    )
                ],
            )

            members.append(
                Membership(
                    member=member_b,
                    annotations=participant_b_annotations if participant_b_annotations else None
                )
            )

        # Create the interaction entity
        yield Entity(
            source='signor',
            type='OM:0013',  # EntityTypeCv.INTERACTION
            identifiers=identifiers,
            annotations=annotations if annotations else None,
            members=members if members else None,
        )
