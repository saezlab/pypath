"""Special-case ID translation routing.

Handles ID types that need non-standard translation paths:
- Gene symbols: case normalization, synonym lookup
- RefSeq: version number handling
- Chain translation: via intermediate ID type
- miRNA: reciprocal name type fallback
- CURIE prefix stripping
"""

from __future__ import annotations

import logging

_log = logging.getLogger(__name__)


def map_genesymbol_fallbacks(
    name: str,
    target_id_type: str,
    ncbi_tax_id: int,
    mapper,
    strict: bool = False,
) -> set[str]:
    """Try various gene symbol normalization strategies.

    Order of attempts:
    1. Direct lookup (exact match)
    2. UPPER case
    3. Capitalized (for rodent symbols like Trp53)
    4. Synonym lookup (genesymbol-syn)
    5. If not strict: append "1" (isoform suffixes)
    """

    # 1. Direct lookup (already tried by caller)

    # 2. Try uppercase
    result = mapper._direct_lookup(
        name.upper(), 'genesymbol', target_id_type, ncbi_tax_id
    )
    if result:
        return result

    # 3. Try capitalized (first letter upper, rest lower)
    result = mapper._direct_lookup(
        name.capitalize(), 'genesymbol', target_id_type, ncbi_tax_id
    )
    if result:
        return result

    # 4. Try gene symbol synonyms
    result = mapper._direct_lookup(
        name, 'genesymbol-syn', target_id_type, ncbi_tax_id
    )
    if result:
        return result
    result = mapper._direct_lookup(
        name.upper(), 'genesymbol-syn', target_id_type, ncbi_tax_id
    )
    if result:
        return result

    if strict:
        return set()

    # 5. Try with "1" appended (for isoform suffixes)
    result = mapper._direct_lookup(
        f'{name}1', 'genesymbol', target_id_type, ncbi_tax_id
    )
    if result:
        return result

    return set()


def map_refseq(
    name: str,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int,
    mapper,
    strict: bool = False,
) -> set[str]:
    """Handle RefSeq IDs which may have version suffixes (.1, .2, etc.)."""

    # Try exact match first
    result = mapper._direct_lookup(name, id_type, target_id_type, ncbi_tax_id)
    if result:
        return result

    # Try without version number
    if '.' in name:
        base = name.rsplit('.', 1)[0]
        result = mapper._direct_lookup(
            base, id_type, target_id_type, ncbi_tax_id
        )
        if result:
            return result

    if strict:
        return set()

    # Try common version numbers
    base = name.split('.')[0] if '.' in name else name
    for ver in range(1, 20):
        result = mapper._direct_lookup(
            f'{base}.{ver}', id_type, target_id_type, ncbi_tax_id
        )
        if result:
            return result

    return set()


def map_ensembl_strip_version(
    name: str,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int,
    mapper,
) -> set[str]:
    """Ensembl IDs sometimes have version suffixes (.1, .2, etc.)."""

    if '.' in name:
        base = name.rsplit('.', 1)[0]
        return mapper._direct_lookup(base, id_type, target_id_type, ncbi_tax_id)

    return set()


def chain_map(
    name: str,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int,
    mapper,
    via: str = 'uniprot',
) -> set[str]:
    """Two-step translation: source -> intermediate -> target."""

    intermediates = mapper.map_name(name, id_type, via, ncbi_tax_id)

    if not intermediates:
        return set()

    result = set()
    for intermediate in intermediates:
        targets = mapper.map_name(
            intermediate, via, target_id_type, ncbi_tax_id
        )
        result.update(targets)

    return result


def map_mirna_fallback(
    name: str,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int,
    mapper,
) -> set[str]:
    """MiRNA name fallback: try reciprocal name types.

    miRNA names in data sources often confuse mature and precursor forms.
    If direct mapping fails, try the other form:
    - mir-mat-name (mature, e.g. hsa-miR-21-5p) <-> mir-name (precursor, e.g. hsa-mir-21)

    Strategy: swap the name type assumption, map to miRBase accession
    as intermediate, then to the target.
    """

    # Define reciprocal pairs: (assumed_type, try_type, intermediate, other_intermediate)
    RECIPROCALS = [
        ('mir-name', 'mir-mat-name', 'mir-pre', 'mirbase'),
        ('mir-mat-name', 'mir-name', 'mirbase', 'mir-pre'),
    ]

    for assumed, try_as, _inter_a, inter_b in RECIPROCALS:
        if id_type != assumed:
            continue

        # Try as the other name type -> get accession
        accessions = mapper._direct_lookup(name, try_as, inter_b, ncbi_tax_id)

        if not accessions:
            continue

        if target_id_type == inter_b:
            return accessions

        # Chain: accession -> target
        result = set()
        for acc in accessions:
            targets = mapper._direct_lookup(
                acc, inter_b, target_id_type, ncbi_tax_id
            )
            result.update(targets)

        if result:
            return result

    return set()


def strip_prefix(
    name: str,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int,
    mapper,
) -> set[str]:
    """Try removing a CURIE-style prefix (e.g. CHEBI:12345 -> 12345)."""
    if ':' not in name:
        return set()
    stripped = name.split(':', 1)[1]
    return mapper._direct_lookup(stripped, id_type, target_id_type, ncbi_tax_id)


# ------------------------------------------------------------------
# HMDB ID normalisation
# ------------------------------------------------------------------

import re

_HMDB_RE = re.compile(r"^HMDB(\d+)$", re.IGNORECASE)


def normalise_hmdb(hmdb_id: str) -> str:
    """Normalise an HMDB ID to the current 7-digit format.

    Converts old 5-digit format (``HMDB00001``) to the current 7-digit
    format (``HMDB0000001``).  IDs already in 7-digit format are
    returned unchanged.  Non-HMDB strings are returned as-is.

    Args:
        hmdb_id: HMDB identifier string.

    Returns:
        Normalised HMDB ID string.

    Examples:
        >>> normalise_hmdb("HMDB00001")
        "HMDB0000001"
        >>> normalise_hmdb("HMDB0000001")
        "HMDB0000001"
    """

    m = _HMDB_RE.match(hmdb_id)

    if m is None:
        return hmdb_id

    digits = m.group(1)

    if len(digits) < 7:
        digits = digits.zfill(7)

    return f"HMDB{digits}"
