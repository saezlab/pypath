"""Unified translation core -- single fallback chain for all APIs."""

from __future__ import annotations

import logging

from omnipath_utils._constants import DEFAULT_ORGANISM
from omnipath_utils.mapping._id_types import IdTypeRegistry

_log = logging.getLogger(__name__)



_db_engine = None


def _get_db_session():
    """Get a DB session if OMNIPATH_UTILS_DB_URL is set.

    Caches the SQLAlchemy engine to avoid re-creating it on every call.
    """
    global _db_engine

    import os
    db_url = os.environ.get("OMNIPATH_UTILS_DB_URL")
    if not db_url:
        return None
    try:
        from sqlalchemy import create_engine
        from sqlalchemy.orm import Session

        if _db_engine is None:
            _db_engine = create_engine(db_url)

        return Session(_db_engine)
    except Exception:
        return None


def _translate_via_db(
    identifiers, id_type, target_id_type, ncbi_tax_id,
    raw=False, uniprot_cleanup=True, full_uniprot='fallback',
):
    """Translate via PostgreSQL if available. Returns None if DB not configured."""
    session = _get_db_session()
    if session is None:
        return None
    try:
        from omnipath_utils.db._query import translate_ids
        result, _backends = translate_ids(
            session, identifiers, id_type, target_id_type, ncbi_tax_id,
            full_uniprot=full_uniprot,
        )
        if not raw and target_id_type == "uniprot" and uniprot_cleanup:
            from omnipath_utils.mapping._cleanup import uniprot_cleanup_batch
            result = uniprot_cleanup_batch(result, ncbi_tax_id, session=session)
        return result
    except Exception as exc:
        _log.debug("DB translate failed: %s", exc)
        return None
    finally:
        session.close()


def translate_core(
    identifiers: list[str],
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
    strict: bool = False,
    raw: bool = False,
    backend: str | None = None,
    uniprot_cleanup: bool = True,
    full_uniprot: str = 'fallback',
) -> dict[str, set[str]]:
    """Unified translation with fallback chain.

    This function is the single entry point for ALL translation APIs.
    It works in memory mode (dict lookup) and can be extended for DB mode.

    Args:
        identifiers: List of source identifiers.
        id_type: Source ID type.
        target_id_type: Target ID type.
        ncbi_tax_id: Organism (default: 9606).
        strict: Skip fuzzy fallbacks (gene symbol case, append "1", etc.).
        raw: Skip ALL special-case handling. Direct table lookup only.
        backend: Force a specific backend. None = auto-select.
        uniprot_cleanup: Apply UniProt cleanup when target is uniprot.
        full_uniprot: Use of the comprehensive full-UniProt table (DB mode):
            ``'fallback'`` (default), ``'never'``, ``'both'``, ``'only'``.

    Returns:
        Dict mapping each source ID to a set of target IDs.
    """
    # Try DB mode first (fast SQL queries)
    db_result = _translate_via_db(
        identifiers, id_type, target_id_type, ncbi_tax_id,
        raw=raw, uniprot_cleanup=uniprot_cleanup, full_uniprot=full_uniprot,
    )
    if db_result is not None:
        return db_result

    # Fall back to in-memory mode
    from omnipath_utils.mapping._mapper import Mapper

    ncbi_tax_id = ncbi_tax_id or DEFAULT_ORGANISM
    reg = IdTypeRegistry.get()
    id_type = reg.resolve(id_type) or id_type

    # Normalise HMDB IDs (old 5-digit → 7-digit format)
    if id_type == "hmdb":
        from omnipath_utils.mapping._special import normalise_hmdb
        identifiers = [normalise_hmdb(i) for i in identifiers]
    target_id_type = reg.resolve(target_id_type) or target_id_type

    mapper = Mapper.get()

    if raw:
        # Raw mode: direct table lookup only, no fallbacks
        table = mapper.which_table(
            id_type,
            target_id_type,
            ncbi_tax_id,
            backend=backend,
        )
        if table:
            return {name: table[name] for name in identifiers}
        return {name: set() for name in identifiers}

    # Full mode: use map_name per ID (includes all special cases)
    # For vectorized performance, try table lookup first for the batch
    table = mapper.which_table(
        id_type,
        target_id_type,
        ncbi_tax_id,
        backend=backend,
    )

    result = {}
    needs_fallback = []

    if table:
        for name in identifiers:
            hits = table[name]
            if hits:
                result[name] = hits
            else:
                needs_fallback.append(name)
    else:
        needs_fallback = list(identifiers)

    # For IDs that did not get a direct hit, use full map_name
    if needs_fallback:
        for name in needs_fallback:
            result[name] = mapper.map_name(
                name,
                id_type,
                target_id_type,
                ncbi_tax_id,
                strict=strict,
                uniprot_cleanup_flag=uniprot_cleanup,
            )

    # Apply cleanup to ALL results if target is uniprot
    # (the direct table hits did not go through map_name cleanup)
    if uniprot_cleanup and target_id_type == 'uniprot' and table:
        from omnipath_utils.mapping._cleanup import (
            uniprot_cleanup as _cleanup,
        )

        for name in identifiers:
            if name in result and result[name] and name not in needs_fallback:
                result[name] = _cleanup(
                    result[name],
                    ncbi_tax_id,
                    mapper=mapper,
                )

    return result
