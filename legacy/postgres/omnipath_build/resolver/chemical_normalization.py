"""Canonical form for a chemical identifier value, by build identifier type
(spec 011 T040).

WP1's resolver lookup joins on the normalized key, not the raw one -- a
mention citing ``chebi:15377`` and one citing ``CHEBI:15377`` must produce
the identical join key. The regexes here mirror
``omnipath_utils.mapping._id_types._NORMALIZERS`` (T024-T026) exactly, one
per chemical namespace. They are duplicated rather than imported: the two
repos share no runtime dependency, the same choice already made for
``resolver_chemical.sql``'s ``source_pattern`` CTE and
``omnipath_build/db/schema.py``'s ``_CHEMICAL_VALUE_PATTERNS``. Keep all
three in sync with any change to the omnipath-utils source of truth.

This module is a hot path -- called once per identifier of every entity a
build ingests -- so :func:`normalize_chemical_identifier` does one dict
lookup before ever touching a regex, and returns ``None`` immediately for
every non-chemical identifier type (the overwhelming majority of calls).
"""

from __future__ import annotations

from collections.abc import Callable
import re


def _normalize_chebi(value: str) -> str | None:
    m = re.fullmatch(r'(?:chebi:)?\s*(\d+)', value.strip(), re.IGNORECASE)
    return f'CHEBI:{m.group(1)}' if m else None


def _normalize_hmdb(value: str) -> str | None:
    m = re.fullmatch(r'hmdb(\d+)', value.strip(), re.IGNORECASE)
    return f'HMDB{m.group(1).zfill(7)}' if m else None


def _normalize_swisslipids(value: str) -> str | None:
    m = re.fullmatch(r'(?:slm:)?\s*(\d+)', value.strip(), re.IGNORECASE)
    return f'SLM:{m.group(1)}' if m else None


def _normalize_lipidmaps(value: str) -> str | None:
    v = value.strip()
    return v.upper() if re.fullmatch(r'lm[a-z0-9]+', v, re.IGNORECASE) else None


def _normalize_kegg(value: str) -> str | None:
    m = re.fullmatch(r'c(\d+)', value.strip(), re.IGNORECASE)
    return f'C{m.group(1).zfill(5)}' if m else None


def _normalize_pubchem(value: str) -> str | None:
    m = re.fullmatch(r'0*(\d+)', value.strip())
    return m.group(1) if m else None


def _normalize_chembl(value: str) -> str | None:
    v = value.strip().replace(' ', '')
    m = re.fullmatch(r'chembl(\d+)', v, re.IGNORECASE)
    return f'CHEMBL{m.group(1)}' if m else None


def _normalize_inchikey(value: str) -> str | None:
    m = re.fullmatch(
        r'(?:inchikey=)?([a-z]{14}-[a-z]{10}-[a-z])', value.strip(), re.IGNORECASE,
    )
    return m.group(1).upper() if m else None


def _build_normalizers() -> dict[str, Callable[[str], str | None]]:
    from omnipath_build.duckdb_load import (
        RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE,
        STANDARD_INCHI_KEY_TYPE,
    )

    by_slug: dict[str, Callable[[str], str | None]] = {
        'chebi': _normalize_chebi,
        'hmdb': _normalize_hmdb,
        'swisslipids': _normalize_swisslipids,
        'lipidmaps': _normalize_lipidmaps,
        'kegg': _normalize_kegg,
        'pubchem': _normalize_pubchem,
        'chembl': _normalize_chembl,
    }
    normalizers = {
        type_name: by_slug[slug]
        for slug, type_name in RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE.items()
        if slug in by_slug
    }
    normalizers[STANDARD_INCHI_KEY_TYPE] = _normalize_inchikey
    return normalizers


_NORMALIZERS: dict[str, Callable[[str], str | None]] | None = None


def normalize_chemical_identifier(
    identifier_type: str,
    value: str | None,
) -> str | None:
    """The normalized form of ``value`` for a chemical build identifier
    type (a ``vocab_identifier_type.name`` such as ``'Chebi:MI:0474'``).

    Returns ``None`` for a non-chemical identifier type, an empty value, or
    a value that does not match that namespace's raw form at all -- exactly
    the cases where the raw ``value`` must be kept as the join key falls
    back to it.
    """

    global _NORMALIZERS
    if not value:
        return None
    if _NORMALIZERS is None:
        _NORMALIZERS = _build_normalizers()
    normalizer = _NORMALIZERS.get(identifier_type)
    return normalizer(value) if normalizer else None


def is_chemical_identifier_type(identifier_type: str | None) -> bool:
    """Whether ``identifier_type`` is one of the chemical namespaces this
    module normalizes -- the structure key or one of the seven
    cross-reference namespaces (spec 011 T024-T026/T040).
    """

    global _NORMALIZERS
    if identifier_type is None:
        return False
    if _NORMALIZERS is None:
        _NORMALIZERS = _build_normalizers()
    return identifier_type in _NORMALIZERS
