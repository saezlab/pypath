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
Shared utilities for parsing OBO (Open Biomedical Ontologies) format files.

This module provides common functionality for processing OBO files
into Entity records with CV_TERM type.
"""

from __future__ import annotations

from typing import Dict, Set, Any
from pypath.formats.obo import Obo, OboValue


def process_obo_term(term) -> Dict[str, str]:
    """
    Process a single OBO term record into a flat dictionary suitable for tabular_builder.

    This is the shared processing function that extracts fields from OBO term records
    and returns them as semicolon-delimited strings where needed (for multi-value fields).

    Args:
        term: An OBO term record from the Obo parser

    Returns:
        Flat dictionary with string values, suitable for Column processing:
        - accession: Term identifier (e.g., "GO:0005737")
        - name: Term name
        - definition: Term definition text (if present)
        - namespace: Ontology namespace (if present)
        - synonyms: Semicolon-delimited synonym strings
        - alt_ids: Semicolon-delimited alternative IDs
        - is_a: Semicolon-delimited parent term IDs
        - xrefs: Semicolon-delimited cross-references
        - comment: Optional comment text
        - relationship_*: Semicolon-delimited targets for each relationship type
        - is_obsolete: 'true' or empty string
    """
    # Extract single-value fields
    term_id = term.id.value if term.id else ''
    name = term.name.value if term.name else ''
    definition = term.definition.value if term.definition else ''
    namespace = term.namespace.value if term.namespace else ''

    # Process multi-value attributes
    synonyms = []
    is_a = []
    relationships = {}
    xrefs = []
    alt_ids = []
    is_obsolete = ''
    comment = ''

    if hasattr(term, 'attrs') and term.attrs:
        # Extract synonyms
        if 'synonym' in term.attrs:
            synonyms = [syn.value for syn in term.attrs['synonym']]

        # Extract is_a relationships (parent terms)
        if 'is_a' in term.attrs:
            is_a = [parent.value.split('!')[0].strip() for parent in term.attrs['is_a']]

        # Extract other relationships (part_of, regulates, etc.)
        if 'relationship' in term.attrs:
            for rel in term.attrs['relationship']:
                # Format: "relationship_type TERM_ID ! comment"
                parts = rel.value.split()
                if len(parts) >= 2:
                    rel_type = parts[0]
                    rel_target = parts[1].split('!')[0].strip()
                    if rel_type not in relationships:
                        relationships[rel_type] = []
                    relationships[rel_type].append(rel_target)

        # Extract cross-references
        if 'xref' in term.attrs:
            xrefs = [xref.value for xref in term.attrs['xref']]

        # Extract alternative IDs
        if 'alt_id' in term.attrs:
            alt_ids = [alt_id.value for alt_id in term.attrs['alt_id']]

        # Check if obsolete
        if 'is_obsolete' in term.attrs:
            if any(obs.value.lower() == 'true' for obs in term.attrs['is_obsolete']):
                is_obsolete = 'true'

        # Extract comment
        if 'comment' in term.attrs:
            comments = list(term.attrs['comment'])
            if comments:
                # OboValue has value and modifiers fields - concatenate them
                obo_val = comments[0]
                parts = [obo_val.value]
                if obo_val.modifiers:
                    parts.append(obo_val.modifiers)
                comment = ' '.join(parts)
            else:
                comment = ''

    # Build flat dictionary with semicolon-delimited multi-values
    result = {
        'accession': term_id,
        'name': name,
        'definition': definition,
        'namespace': namespace,
        'synonyms': ';'.join(synonyms),
        'alt_ids': ';'.join(alt_ids),
        'is_a': ';'.join(is_a),
        'xrefs': ';'.join(xrefs),
        'comment': comment,
        'is_obsolete': is_obsolete,
    }

    # Add relationship columns dynamically
    for rel_type, targets in relationships.items():
        result[f'relationship_{rel_type}'] = ';'.join(targets)

    return result
