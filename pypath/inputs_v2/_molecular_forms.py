"""Preserve asserted MITAB participant forms before identifier normalization.

Feature syntax follows HUPO-PSI MITAB 2.7, columns 37/38:
https://github.com/HUPO-PSI/miTab/blob/master/PSI-MITAB27Format.md
Fuzzy ranges remain descriptions; no exact position is invented for them.
"""

from __future__ import annotations

import re
from typing import Any

from omnipath_core.molecular_forms import (
    molecular_form_from_identifiers,
    normalize_molecular_form,
)

# Exact source/vocabulary labels. Unknown features stay in source descriptions.
# Includes spellings present in the SIGNOR causalTab export, without repairing
# their labels or asserting ontology terms that the export does not supply.
_MODIFICATION_LABELS = {
    'phosphorylated residue', 'de-phosphorylated residue',
    'ubiquitinylated lysine', 'acetylated residue', 'de-acetylated residue',
    'polyubiquitinated residue', 'de-methylated residue', 'sumoylated lysine',
    'monoubiquitinated residue', 'methylated residue', 'carboxylated residue',
    'de-ubiquitinated residue', 'glycosylated residue', 'hydroxylated residue',
    'palmitoylated residue', 'chemical modificated residue',
    'de-glycosylated residue', 'post translatedal modificated residue',
    'de-sumoylated residue', 'neddylated lysine', 'trimethylated residue',
}
_VARIANT_LABELS = {
    'mutation', 'mutation decreasing', 'mutation increasing',
    'mutation disrupting', 'mutation disrupting strength', 'sequence variant',
}


def mitab_participant_form(
        row: dict[str, Any],
        suffix: str,
        *,
        entity_type: str,
    ) -> dict | None:
    """Capture the primary participant ID, not its alternative cross-references."""
    if entity_type not in {'protein', 'rna_product', 'transcript', 'noncoding_rna_product'}:
        return None
    primary = next((row.get(column) for column in (
        f'\ufeff#ID(s) interactor {suffix}',
        f'#ID(s) interactor {suffix}',
        f'ID(s) interactor {suffix}',
    ) if row.get(column)), None)
    identifiers = []
    for token in str(primary or '').split('|'):
        prefix, sep, value = token.strip().strip('"').partition(':')
        if not sep:
            continue
        value = value.strip().strip('"').split('(', 1)[0].strip().strip('"')
        if prefix.lower() in {'uniprot', 'uniprotkb', 'refseq', 'ensembl'}:
            identifiers.append({'ns': prefix.lower(), 'id': value})
    form = molecular_form_from_identifiers(identifiers) or {}
    sequences = form.get('sequence_identifiers') or []
    coordinate = {
        'identifier': sequences[0] if len(sequences) == 1 else None,
        'coordinate_system': 'protein' if entity_type == 'protein' else 'transcript',
        'position_base': 1,
    }
    modifications, variants = [], []
    for feature in str(row.get(f'Feature(s) interactor {suffix}') or '').split('|'):
        feature = feature.strip()
        if not feature or feature == '-':
            continue
        # Controlled accessions contain a colon, unlike feature label names.
        match = re.fullmatch(r'((?:MOD|MI):\d+|[^:]+):(.+)', feature)
        if not match:
            continue
        term, ranges = match.groups()
        lower = term.lower()
        is_variant = (lower in _VARIANT_LABELS
                      or term in {'MI:0118', 'MI:0119', 'MI:0120', 'MI:0121', 'MI:1128'})
        is_modification = term.startswith('MOD:') or lower in _MODIFICATION_LABELS
        if not (is_variant or is_modification):
            continue
        # Descriptions retain fuzzy/unknown/n/c ranges and source feature labels.
        for range_value in ranges.split(','):
            position = re.fullmatch(r'(\d+)-(\d+)(?:\((.*)\))?', range_value.strip())
            start, end = (int(position[1]), int(position[2])) if position else (None, None)
            item = {'position': start, 'end_position': end,
                    'coordinate_reference': coordinate, 'description': feature}
            if is_variant:
                substitution = re.search(r'(?:p\.)?([A-Z*])([0-9]+)([A-Z*])',
                                         position[3] or '') if position else None
                if substitution and start == end == int(substitution[2]):
                    item.update(reference=substitution[1], alternate=substitution[3])
                variants.append(item)
            else:
                item['term'] = term
                modifications.append(item)
    form.update(modifications=modifications, variants=variants)
    return normalize_molecular_form(form, allow_resolved=False)
