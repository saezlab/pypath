"""Preserve asserted MITAB participant forms before identifier normalization.

Feature syntax follows HUPO-PSI MITAB 2.7, columns 37/38:
https://github.com/HUPO-PSI/miTab/blob/master/PSI-MITAB27Format.md
Fuzzy ranges remain descriptions; no exact position is invented for them.
"""

from __future__ import annotations

import hashlib
import re
from typing import Any

from omnipath_core.molecular_forms import (
    molecular_form_from_identifiers,
    normalize_molecular_form,
)

# Source/vocabulary feature labels; unknown features stay in source
# descriptions. Labels are not repaired and no ontology accession is asserted
# that the export does not supply.
# IntAct and SIGNOR mutation labels all start with "mutation" ("mutation with
# no effect", "mutation decreasing interaction strength", ...).
_VARIANT_LABELS = {
    'variant',
    'sequence variant',
    'disease causing amino-acid variant',
    'trapping mutant',
}
# PSI-MI binding region feature types; the label is the region type.
_REGION_LABELS = {
    'binding-associated region',
    'necessary binding region',
    'sufficient binding region',
    'direct binding region',
}
# Modified residues are named by PSI-MOD terms ("O4'-phospho-L-tyrosine",
# "N6,N6-dimethyl-L-lysine", "acylated residue", "protein modification") or,
# in the SIGNOR causalTab export, similar labels ("phosphorylated residue",
# "sumoylated lysine", "post translatedal modificated residue").
_MODIFICATION_NAME_RE = re.compile(
    r'\b(?:alanine|arginine|asparagine|aspartic acid|cysteine|cystine|'
    r'glutamic acid|glutamine|glycine|histidine|isoleucine|leucine|lysine|'
    r'methionine|selenomethionine|phenylalanine|proline|serine|threonine|'
    r'tryptophan|tyrosine|valine|citrulline|residues?|modification)\b'
)


def sequence_form(sequence: object, *, system: str = 'protein') -> dict | None:
    """Identify an explicitly supplied sequence without assigning an accession.

    Hash the uppercase, whitespace-free sequence. RNA and DNA stay distinct;
    no T/U conversion or reference-sequence alignment is performed. The source
    sequence remains in the raw record.
    """
    alphabets = {
        'protein': 'ACDEFGHIKLMNPQRSTVWYBXZJUO*',
        'transcript': 'ACGUNRYSWKMBDHV',
        'genomic': 'ACGTNRYSWKMBDHV',
    }
    if system not in alphabets:
        raise ValueError(f'Unsupported sequence system: {system}')
    text = str(sequence or '').strip()
    if text.lower() in {
        'sequence unavailable',
        'not available',
        'none',
        'null',
        'n/a',
        '-',
    }:
        return None
    value = re.sub(r'\s+', '', text).upper()
    if not value or any(letter not in alphabets[system] for letter in value):
        return None
    identifier = hashlib.sha256(value.encode('ascii')).hexdigest()
    return normalize_molecular_form(
        {
            'sequence_identifiers': [
                {'ns': f'{system}_sequence_sha256', 'id': identifier}
            ],
        },
        allow_resolved=False,
    )


def combine_forms(*forms: dict | None) -> dict | None:
    """Combine descriptions of one occurrence while retaining ambiguous IDs."""
    values = [
        normalize_molecular_form(form, allow_resolved=False)
        for form in forms
        if form
    ]
    result: dict = {}
    for key in (
        'sequence_identifiers',
        'modifications',
        'variants',
        'regions',
    ):
        items = []
        for form in values:
            for item in (form or {}).get(key) or []:
                if key != 'sequence_identifiers' or item not in items:
                    items.append(item)
        result[key] = items or None
    isoforms = []
    for form in values:
        isoform = (form or {}).get('isoform_identifier')
        if isoform and isoform not in isoforms:
            isoforms.append(isoform)
    result['isoform_identifier'] = isoforms[0] if len(isoforms) == 1 else None
    if len(isoforms) > 1:
        sequences = result['sequence_identifiers'] or []
        result['sequence_identifiers'] = sequences + [
            i for i in isoforms if i not in sequences
        ]
    return normalize_molecular_form(result, allow_resolved=False)


def protein_variants(
    description: object,
    *,
    identifier: dict | None = None,
    coordinate_reference: dict | None = None,
) -> list[dict]:
    """Parse explicit single-letter substitutions, retaining other source syntax.

    Each source token is retained. Unknown/deletion/insertion descriptions do
    not acquire guessed positions or alleles, and multiple substitutions are
    individual features of the same reported form.
    """
    text = str(description or '').strip()
    if not text and identifier is None:
        return []
    tokens = re.split(r'\s*[,;/]\s*', text) if text else ['']
    variants = []
    for token in tokens:
        item = {
            'identifier': identifier,
            'description': token or text or None,
            'coordinate_reference': coordinate_reference,
        }
        substitution = re.fullmatch(
            r'(?:p\.)?([ACDEFGHIKLMNPQRSTVWY*])([1-9]\d*)([ACDEFGHIKLMNPQRSTVWY*])',
            token,
        )
        if substitution:
            item.update(
                reference=substitution[1],
                alternate=substitution[3],
                position=int(substitution[2]),
                end_position=int(substitution[2]),
            )
        variants.append(item)
    # Conflicting substitutions at one position are alternatives or an
    # ambiguous source description, not a construct carrying both alleles.
    by_position = {}
    for item in variants:
        position = item.get('position')
        if position is not None:
            alternate = item.get('alternate')
            if position in by_position and by_position[position] != alternate:
                variants = [
                    {
                        'identifier': identifier,
                        'description': text,
                        'coordinate_reference': coordinate_reference,
                    }
                ]
                break
            by_position[position] = alternate
    return normalize_molecular_form(
        {'variants': variants}, allow_resolved=False
    )['variants']


def mitab_participant_form(
    row: dict[str, Any],
    suffix: str,
    *,
    entity_type: str,
) -> dict | None:
    """Capture the primary participant ID, not its alternative cross-references."""
    if entity_type not in {
        'protein',
        'rna_product',
        'transcript',
        'noncoding_rna_product',
    }:
        return None
    primary = next(
        (
            row.get(column)
            for column in (
                f'\ufeff#ID(s) interactor {suffix}',
                f'#ID(s) interactor {suffix}',
                f'ID(s) interactor {suffix}',
            )
            if row.get(column)
        ),
        None,
    )
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
        'coordinate_system': 'protein'
        if entity_type == 'protein'
        else 'transcript',
        'position_base': 1,
    }
    modifications, variants, regions = [], [], []
    for feature in str(row.get(f'Feature(s) interactor {suffix}') or '').split(
        '|'
    ):
        feature = feature.strip()
        if not feature or feature == '-':
            continue
        # Controlled accessions contain a colon, unlike feature label names.
        match = re.fullmatch(r'((?:MOD|MI):\d+|[^:]+):(.+)', feature)
        if not match:
            continue
        term, ranges = match.groups()
        lower = term.strip('"').lower()
        is_variant = (
            lower.startswith('mutation')
            or lower in _VARIANT_LABELS
            or term in {'MI:0118', 'MI:0119', 'MI:0120', 'MI:0121', 'MI:1128'}
        )
        is_region = lower in _REGION_LABELS
        is_modification = not (is_variant or is_region) and (
            term.startswith('MOD:') or bool(_MODIFICATION_NAME_RE.search(lower))
        )
        if not (is_variant or is_region or is_modification):
            continue
        # Descriptions retain fuzzy/unknown/n/c ranges and source feature labels.
        for range_value in ranges.split(','):
            position = re.fullmatch(
                r'(\d+)-(\d+)(?:\((.*)\))?', range_value.strip()
            )
            start, end = (
                (int(position[1]), int(position[2]))
                if position and 0 < int(position[1]) <= int(position[2])
                else (None, None)
            )
            item = {
                'position': start,
                'end_position': end,
                'coordinate_reference': coordinate,
                'description': feature,
            }
            if is_variant:
                substitution = (
                    re.search(
                        r'(?:p\.)?([A-Z*])([0-9]+)([A-Z*])', position[3] or ''
                    )
                    if position
                    else None
                )
                if substitution and start == end == int(substitution[2]):
                    item.update(
                        reference=substitution[1], alternate=substitution[3]
                    )
                variants.append(item)
            elif is_region:
                regions.append({'type': lower, **item})
            else:
                item['term'] = term.strip('"')
                modifications.append(item)
    form.update(modifications=modifications, variants=variants, regions=regions)
    return normalize_molecular_form(form, allow_resolved=False)
