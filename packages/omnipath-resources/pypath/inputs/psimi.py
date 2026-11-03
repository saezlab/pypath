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

from typing import List, Optional, NamedTuple

import pypath.resources.urls as urls
import pypath.formats.obo as obo


class PsimiRecord(NamedTuple):
    """PSI-MI ontology term record."""
    id: str
    name: Optional[str]
    definition: Optional[str]
    definition_refs: Optional[str]  # Format: "[PMID:12345,PMID:67890]"
    parent_ids: Optional[List[str]]
    parent_names: Optional[List[str]]
    synonyms: Optional[List[str]]
    alt_ids: Optional[List[str]]


def psimi_ontology(flatten=True):
    """
    Downloads and parses the PSI-MI (Molecular Interaction) ontology.

    Args:
        flatten: If True, returns flattened PsimiRecord named tuples.
            If False, returns raw OBO records.

    Returns:
        An iterator of PsimiRecord named tuples representing PSI-MI ontology terms.

        When flatten=True, each record contains:
        - id: The PSI-MI identifier (e.g., 'MI:0001')
        - name: The term name
        - definition: The term definition text
        - definition_refs: References for the definition (e.g., "[PMID:12345,PMID:67890]")
        - parent_ids: List of parent term IDs from is_a relationships
        - parent_names: List of parent term names
        - synonyms: List of synonym texts
        - alt_ids: List of alternative IDs
    """

    url = urls.urls['psimi']['url']
    reader = obo.Obo(url)

    if not flatten:
        return reader

    for record in reader:
        yield _flatten_obo_record(record)


def _flatten_obo_record(record):
    """
    Flatten OBO record for easier consumption.
    Returns structured data with lists instead of pipe-separated strings.
    """

    # Extract basic fields
    term_id = record.id.value if record.id else None
    name = record.name.value if record.name else None

    # Extract definition and references
    definition = None
    definition_refs = None
    if record.definition:
        definition = record.definition.value
        definition_refs = record.definition.modifiers

    # Extract parent relationships (is_a)
    parent_ids = None
    parent_names = None
    if 'is_a' in record.attrs:
        parents = list(record.attrs['is_a'])
        parent_ids = [p.value for p in parents]
        parent_names = [p.comment for p in parents if p.comment]

    # Extract synonyms
    synonyms = None
    if 'synonym' in record.attrs:
        synonym_list = list(record.attrs['synonym'])
        synonyms = [s.value for s in synonym_list]

    # Extract alternative IDs
    alt_ids = None
    if 'alt_id' in record.attrs:
        alt_ids = [a.value for a in record.attrs['alt_id']]

    return PsimiRecord(
        id=term_id,
        name=name,
        definition=definition,
        definition_refs=definition_refs,
        parent_ids=parent_ids,
        parent_names=parent_names,
        synonyms=synonyms,
        alt_ids=alt_ids,
    )
