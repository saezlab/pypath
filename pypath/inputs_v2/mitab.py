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

import re
from typing import Optional, List, Set, Generator, NamedTuple

import pypath.share.common as common


class MitabInteraction(NamedTuple):
    """MITAB 2.8 format interaction record (46 columns)."""
    id_a: str  # 1. Unique identifier for interactor A
    id_b: str  # 2. Unique identifier for interactor B
    alt_ids_a: str  # 3. Alternative identifier for interactor A
    alt_ids_b: str  # 4. Alternative identifier for interactor B
    aliases_a: str  # 5. Aliases for A
    aliases_b: str  # 6. Aliases for B
    detection_methods: str  # 7. Interaction detection methods
    first_author: str  # 8. First author
    pmids: str  # 9. Publication identifier(s)
    taxid_a: str  # 10. NCBI Taxonomy identifier for interactor A
    taxid_b: str  # 11. NCBI Taxonomy identifier for interactor B
    interaction_types: str  # 12. Interaction types
    source_dbs: str  # 13. Source databases
    interaction_ids: str  # 14. Interaction identifier(s)
    confidence_scores: str  # 15. Confidence score
    complex_expansion: str  # 16. Complex expansion
    biological_role_a: str  # 17. Biological role A
    biological_role_b: str  # 18. Biological role B
    experimental_role_a: str  # 19. Experimental role A
    experimental_role_b: str  # 20. Experimental role B
    interactor_type_a: str  # 21. Interactor type A
    interactor_type_b: str  # 22. Interactor type B
    xrefs_a: str  # 23. Xref for interactor A
    xrefs_b: str  # 24. Xref for interactor B
    xrefs_interaction: str  # 25. Xref for interaction
    annotations_a: str  # 26. Annotations for interactor A
    annotations_b: str  # 27. Annotations for interactor B
    annotations_interaction: str  # 28. Annotations for interaction
    host_organism: str  # 29. Host organism
    parameters: str  # 30. Parameters
    creation_date: str  # 31. Creation date
    update_date: str  # 32. Update date
    checksum_a: str  # 33. Checksum for interactor A
    checksum_b: str  # 34. Checksum for interactor B
    checksum_interaction: str  # 35. Checksum for interaction
    negative: str  # 36. Negative
    features_a: str  # 37. Features for interactor A
    features_b: str  # 38. Features for interactor B
    stoichiometry_a: str  # 39. Stoichiometry for interactor A
    stoichiometry_b: str  # 40. Stoichiometry for interactor B
    identification_method_a: str  # 41. Participant identification method for interactor A
    identification_method_b: str  # 42. Participant identification method for interactor B
    biological_effect_a: str  # 43. Biological effect of interactor A
    biological_effect_b: str  # 44. Biological effect of interactor B
    causal_regulatory_mechanism: str  # 45. Causal regulatory mechanism
    causal_statement: str  # 46. Causal statement


def mitab_field_list(field):
    """
    Extract list of names from a MITAB field with format: db:id(name)|db:id(name)
    Returns the names in parentheses.
    """
    if not field or field == '-':
        return []

    return common.unique_list(
        map(
            lambda x: x.split('(')[1][:-1] if '(' in x else x,
            field.split('|')
        )
    )


def mitab_field_uniprot(field):
    """
    Extract UniProt ID from a MITAB identifier field.
    """
    if not field or field == '-':
        return None

    uniprots = list(
        filter(
            lambda x: len(x) == 2 and x[0] == 'uniprotkb',
            map(
                lambda x: x.split(':'),
                field.split('|')
            )
        )
    )

    return uniprots[0][1].replace('"', '') if uniprots else None


def mitab_parse_identifier(field: str) -> Optional[dict]:
    """
    Parse a MITAB identifier field (format: databaseName:identifier).
    Returns a dict with 'database' and 'id' keys.
    """
    if not field or field == '-':
        return None

    # Take first identifier if multiple
    first_id = field.split('|')[0]

    if ':' in first_id:
        parts = first_id.split(':', 1)
        return {
            'database': parts[0],
            'id': parts[1].replace('"', '')
        }

    return {'database': 'unknown', 'id': first_id}


def mitab_parse_identifiers(field: str) -> List[dict]:
    """
    Parse multiple MITAB identifiers (format: db:id|db:id|...).
    Returns a list of dicts with 'database' and 'id' keys.
    """
    if not field or field == '-':
        return []

    identifiers = []
    for item in field.split('|'):
        if ':' in item:
            parts = item.split(':', 1)
            identifiers.append({
                'database': parts[0],
                'id': parts[1].replace('"', '')
            })
        else:
            identifiers.append({
                'database': 'unknown',
                'id': item
            })

    return identifiers


def mitab_parse_mi_term(field: str) -> Optional[dict]:
    """
    Parse PSI-MI term field (format: psi-mi:"MI:xxxx"(term name)).
    Returns a dict with 'mi_id', 'name', and optionally 'database' keys.
    """
    if not field or field == '-':
        return None

    # Take first term if multiple
    first_term = field.split('|')[0]

    match = re.match(r'([^:]+):"?([^"(]+)"?\(([^)]+)\)', first_term)
    if match:
        return {
            'database': match.group(1),
            'mi_id': match.group(2),
            'name': match.group(3)
        }

    return None


def mitab_parse_mi_terms(field: str) -> List[dict]:
    """
    Parse multiple PSI-MI term fields (format: psi-mi:"MI:xxxx"(name)|...).
    Returns a list of dicts with 'mi_id', 'name', and 'database' keys.
    """
    if not field or field == '-':
        return []

    terms = []
    for item in field.split('|'):
        match = re.match(r'([^:]+):"?([^"(]+)"?\(([^)]+)\)', item)
        if match:
            terms.append({
                'database': match.group(1),
                'mi_id': match.group(2),
                'name': match.group(3)
            })

    return terms


def mitab_parse_taxid(field: str) -> Optional[dict]:
    """
    Parse taxonomy ID field (format: taxid:9606(Homo sapiens)).
    Returns a dict with 'taxid' and 'organism' keys.
    """
    if not field or field == '-':
        return None

    # Take first taxid if multiple
    first_taxid = field.split('|')[0]

    match = re.match(r'taxid:(-?\d+)(?:\(([^)]+)\))?', first_taxid)
    if match:
        return {
            'taxid': int(match.group(1)),
            'organism': match.group(2) if match.group(2) else None
        }

    return None


def mitab_parse_pubmeds(field: str) -> Set[str]:
    """
    Extract PubMed IDs from publication identifier field.
    Returns a set of PubMed IDs.
    """
    if not field or field == '-':
        return set()

    pubmeds = set()
    for item in field.split('|'):
        if ':' in item:
            db, pub_id = item.split(':', 1)
            if db.lower() == 'pubmed':
                pubmeds.add(pub_id)

    return pubmeds


def mitab_parse_aliases(field: str) -> List[dict]:
    """
    Parse alias field (format: db:name(alias type)|db:name(alias type)).
    Returns a list of dicts with 'database', 'name', and 'type' keys.
    """
    if not field or field == '-':
        return []

    aliases = []
    for item in field.split('|'):
        match = re.match(r'([^:]+):([^(]+)(?:\(([^)]+)\))?', item)
        if match:
            aliases.append({
                'database': match.group(1),
                'name': match.group(2),
                'type': match.group(3) if match.group(3) else None
            })

    return aliases


def mitab_parse_confidence(field: str) -> List[dict]:
    """
    Parse confidence score field (format: scoreType:value|scoreType:value).
    Returns a list of dicts with 'type' and 'value' keys.
    """
    if not field or field == '-':
        return []

    scores = []
    for item in field.split('|'):
        if ':' in item:
            score_type, value = item.split(':', 1)
            scores.append({
                'type': score_type,
                'value': value
            })

    return scores


def mitab_parse_features(field: str) -> List[dict]:
    """
    Parse feature field (format: type:range(text)|type:range(text)).
    Returns a list of dicts with 'type', 'range', and 'text' keys.
    """
    if not field or field == '-':
        return []

    features = []
    for item in field.split('|'):
        match = re.match(r'([^:]+):([^(]+)(?:\(([^)]+)\))?', item)
        if match:
            features.append({
                'type': match.group(1),
                'range': match.group(2),
                'text': match.group(3) if match.group(3) else None
            })

    return features


def mitab_parse_parameters(field: str) -> List[dict]:
    """
    Parse parameter field (format: type:value(text)|type:value(text)).
    Returns a list of dicts with 'type', 'value', and 'text' keys.
    """
    if not field or field == '-':
        return []

    parameters = []
    for item in field.split('|'):
        match = re.match(r'([^:]+):([^(]+)(?:\(([^)]+)\))?', item)
        if match:
            parameters.append({
                'type': match.group(1),
                'value': match.group(2),
                'text': match.group(3) if match.group(3) else None
            })

    return parameters




def mitab_interactions(
        data: Generator[str, None, None],
        organism: Optional[int] = None,
        only_proteins: bool = False,
        skip_header: bool = True,
        columns: int = 15,
    ) -> Generator[MitabInteraction, None, None]:
    """
    Parse a MITAB format file.

    Data is streamed line-by-line for low latency and memory usage.

    Args:
        data: Generator or iterator of MITAB lines
        organism: NCBI taxonomy ID to filter interactions (None = no filter)
        only_proteins: Only keep protein-protein interactions
        skip_header: Skip the first line (header)
        columns: Minimum number of columns expected (15, 25, 27, 42, or 46)

    Yields:
        MitabInteraction namedtuples with all 46 columns
    """
    # Parse lines
    for line_num, line in enumerate(data):

        if skip_header and line_num == 0:
            continue

        if isinstance(line, bytes):
            line = line.decode('utf-8')

        line = line.strip()
        if not line:
            continue

        fields = line.split('\t')

        # Check minimum number of columns
        if len(fields) < columns:
            continue

        # Parse taxonomy IDs for filtering
        taxid_a = mitab_parse_taxid(fields[9] if len(fields) > 9 else '-')
        taxid_b = mitab_parse_taxid(fields[10] if len(fields) > 10 else '-')

        # Apply organism filter
        if organism is not None:
            if not (
                (taxid_a and taxid_a['taxid'] == organism) or
                (taxid_b and taxid_b['taxid'] == organism)
            ):
                continue

        # Apply protein-only filter
        if only_proteins:
            id_a = mitab_parse_identifier(fields[0])
            id_b = mitab_parse_identifier(fields[1])
            if not (
                id_a and id_a['database'] == 'uniprotkb' and
                id_b and id_b['database'] == 'uniprotkb'
            ):
                continue

        # Create interaction record with all 46 fields
        # Use '-' as default for missing fields
        interaction = MitabInteraction(
            id_a=fields[0] if len(fields) > 0 else '-',
            id_b=fields[1] if len(fields) > 1 else '-',
            alt_ids_a=fields[2] if len(fields) > 2 else '-',
            alt_ids_b=fields[3] if len(fields) > 3 else '-',
            aliases_a=fields[4] if len(fields) > 4 else '-',
            aliases_b=fields[5] if len(fields) > 5 else '-',
            detection_methods=fields[6] if len(fields) > 6 else '-',
            first_author=fields[7] if len(fields) > 7 else '-',
            pmids=fields[8] if len(fields) > 8 else '-',
            taxid_a=fields[9] if len(fields) > 9 else '-',
            taxid_b=fields[10] if len(fields) > 10 else '-',
            interaction_types=fields[11] if len(fields) > 11 else '-',
            source_dbs=fields[12] if len(fields) > 12 else '-',
            interaction_ids=fields[13] if len(fields) > 13 else '-',
            confidence_scores=fields[14] if len(fields) > 14 else '-',
            complex_expansion=fields[15] if len(fields) > 15 else '-',
            biological_role_a=fields[16] if len(fields) > 16 else '-',
            biological_role_b=fields[17] if len(fields) > 17 else '-',
            experimental_role_a=fields[18] if len(fields) > 18 else '-',
            experimental_role_b=fields[19] if len(fields) > 19 else '-',
            interactor_type_a=fields[20] if len(fields) > 20 else '-',
            interactor_type_b=fields[21] if len(fields) > 21 else '-',
            xrefs_a=fields[22] if len(fields) > 22 else '-',
            xrefs_b=fields[23] if len(fields) > 23 else '-',
            xrefs_interaction=fields[24] if len(fields) > 24 else '-',
            annotations_a=fields[25] if len(fields) > 25 else '-',
            annotations_b=fields[26] if len(fields) > 26 else '-',
            annotations_interaction=fields[27] if len(fields) > 27 else '-',
            host_organism=fields[28] if len(fields) > 28 else '-',
            parameters=fields[29] if len(fields) > 29 else '-',
            creation_date=fields[30] if len(fields) > 30 else '-',
            update_date=fields[31] if len(fields) > 31 else '-',
            checksum_a=fields[32] if len(fields) > 32 else '-',
            checksum_b=fields[33] if len(fields) > 33 else '-',
            checksum_interaction=fields[34] if len(fields) > 34 else '-',
            negative=fields[35] if len(fields) > 35 else '-',
            features_a=fields[36] if len(fields) > 36 else '-',
            features_b=fields[37] if len(fields) > 37 else '-',
            stoichiometry_a=fields[38] if len(fields) > 38 else '-',
            stoichiometry_b=fields[39] if len(fields) > 39 else '-',
            identification_method_a=fields[40] if len(fields) > 40 else '-',
            identification_method_b=fields[41] if len(fields) > 41 else '-',
            biological_effect_a=fields[42] if len(fields) > 42 else '-',
            biological_effect_b=fields[43] if len(fields) > 43 else '-',
            causal_regulatory_mechanism=fields[44] if len(fields) > 44 else '-',
            causal_statement=fields[45] if len(fields) > 45 else '-',
        )

        yield interaction