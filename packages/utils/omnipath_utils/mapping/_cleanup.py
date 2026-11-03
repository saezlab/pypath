"""UniProt ID cleanup — validate, normalize, prefer SwissProt."""

from __future__ import annotations

import re
import logging

_log = logging.getLogger(__name__)

RE_UNIPROT = re.compile(
    r'^([OPQ][0-9][A-Z0-9]{3}[0-9]'
    r'|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$'
)


def is_uniprot_ac(name: str) -> bool:
    """Check if a string matches UniProt AC format."""
    return bool(RE_UNIPROT.match(name))


def uniprot_cleanup(
    uniprots: set[str],
    ncbi_tax_id: int,
    mapper=None,
) -> set[str]:
    """Clean up UniProt IDs: prefer SwissProt, validate, normalize.

    Pipeline:
    1. Translate secondary -> primary ACs
    2. For TrEMBL IDs: try to find SwissProt via gene symbol
    3. Filter to organism proteome (if reflist available)
    4. Validate AC format
    """
    if not uniprots or not mapper:
        return uniprots or set()

    # Step 1: secondary -> primary
    result = _primary_uniprot(uniprots, ncbi_tax_id, mapper)

    # Step 2: TrEMBL -> SwissProt via gene symbol
    result = _trembl_to_swissprot(result, ncbi_tax_id, mapper)

    # Step 3: filter to organism proteome
    result = _filter_organism(result, ncbi_tax_id)

    # Step 4: format validation
    result = {u for u in result if is_uniprot_ac(u)}

    return result


def _primary_uniprot(
    uniprots: set[str],
    ncbi_tax_id: int,
    mapper,
) -> set[str]:
    """Translate secondary UniProt ACs to primary."""
    result = set()

    for up in uniprots:
        primary = mapper._direct_lookup(
            up, 'uniprot-sec', 'uniprot-pri', ncbi_tax_id
        )
        if primary:
            result.update(primary)
        else:
            result.add(up)  # probably already primary

    return result


def _trembl_to_swissprot(
    uniprots: set[str],
    ncbi_tax_id: int,
    mapper,
) -> set[str]:
    """For TrEMBL IDs, try to find SwissProt via gene symbol.

    For each ID: if it is a SwissProt, keep it. If it is a TrEMBL,
    look up its gene symbol and find the SwissProt with the same symbol.
    If no SwissProt found, keep the TrEMBL.
    """
    try:
        from omnipath_utils.reflists import all_swissprots

        swissprot_set = all_swissprots(ncbi_tax_id)
    except Exception:
        # If reflists not available, can not filter
        return uniprots

    if not swissprot_set:
        return uniprots

    result = set()

    for up in uniprots:
        if up in swissprot_set:
            # Already SwissProt, keep it
            result.add(up)
            continue

        # It is a TrEMBL -- try to find SwissProt via gene symbol
        genesymbols = mapper._direct_lookup(
            up, 'trembl', 'genesymbol', ncbi_tax_id
        )
        if not genesymbols:
            genesymbols = mapper._direct_lookup(
                up, 'uniprot', 'genesymbol', ncbi_tax_id
            )

        if genesymbols:
            # Map gene symbol -> SwissProt
            for gs in genesymbols:
                swissprots = mapper._direct_lookup(
                    gs, 'genesymbol', 'swissprot', ncbi_tax_id
                )
                if swissprots:
                    result.update(swissprots)
                    break
            else:
                # No SwissProt found via gene symbol, keep TrEMBL
                result.add(up)
        else:
            # Can not resolve, keep as-is
            result.add(up)

    return result


def _filter_organism(uniprots: set[str], ncbi_tax_id: int) -> set[str]:
    """Filter to IDs present in the organism proteome."""
    try:
        from omnipath_utils.reflists import all_uniprots

        proteome = all_uniprots(ncbi_tax_id)
    except Exception:
        return uniprots

    if not proteome:
        return uniprots

    filtered = uniprots & proteome

    # If filtering removed everything, return originals
    # (the proteome list might be incomplete)
    return filtered if filtered else uniprots


# ------------------------------------------------------------------
# Batch cleanup — vectorized version for server and DB-backed API
# ------------------------------------------------------------------




def _reflist_from_db(session, list_name: str, ncbi_tax_id: int) -> set[str]:
    """Query a reference list from the database."""
    try:
        from sqlalchemy import text as sql_text
        from omnipath_utils.db._connection import SCHEMA

        rows = session.execute(
            sql_text(f"""
                SELECT identifier
                FROM {SCHEMA}.reflist
                WHERE list_name = :name AND ncbi_tax_id = :tax
            """),
            {"name": list_name, "tax": ncbi_tax_id},
        )
        return {row[0] for row in rows}
    except Exception:
        return set()


def uniprot_cleanup_batch(
    results: dict[str, set[str]],
    ncbi_tax_id: int,
    session=None,
) -> dict[str, set[str]]:
    """Clean up UniProt IDs in batch: prefer SwissProt, validate, normalize.

    Vectorized version of ``uniprot_cleanup()`` — does O(1) SQL queries
    instead of O(n) per-ID lookups.  Operates on the full results dict
    from a translation batch.

    Args:
        results: Dict mapping source IDs to sets of UniProt ACs.
        ncbi_tax_id: NCBI Taxonomy ID.
        session: SQLAlchemy session (required for DB queries).

    Returns:
        Same dict structure with cleaned UniProt ACs.
    """

    if not results or not session:
        return results

    from omnipath_utils.db._query import translate_ids

    # Collect all UniProt ACs across all results
    all_acs: set[str] = set()
    for targets in results.values():
        all_acs.update(targets)

    if not all_acs:
        return results

    ac_list = list(all_acs)

    # Step 1: secondary -> primary (batch)
    sec_pri, _ = translate_ids(session, ac_list, "uniprot-sec", "uniprot-pri", ncbi_tax_id)

    # Step 2: SwissProt membership check (from DB reflist table)
    swissprot_set = _reflist_from_db(session, "swissprot", ncbi_tax_id)

    # Step 3: For TrEMBL IDs, batch lookup gene symbols + swissprot
    trembl_acs = [ac for ac in ac_list if swissprot_set and ac not in swissprot_set]
    trembl_to_gs: dict[str, set[str]] = {}
    gs_to_sp: dict[str, set[str]] = {}

    if trembl_acs:
        trembl_to_gs, _ = translate_ids(
            session, trembl_acs, "trembl", "genesymbol", ncbi_tax_id,
        )
        # Also try uniprot -> genesymbol for IDs not found via trembl
        missing_gs = [ac for ac in trembl_acs if not trembl_to_gs.get(ac)]
        if missing_gs:
            uniprot_gs, _ = translate_ids(
                session, missing_gs, "uniprot", "genesymbol", ncbi_tax_id,
            )
            for ac, gs in uniprot_gs.items():
                if gs:
                    trembl_to_gs[ac] = gs

        # Collect all gene symbols and batch lookup swissprot
        all_gs: set[str] = set()
        for gs_set in trembl_to_gs.values():
            all_gs.update(gs_set)

        if all_gs:
            gs_to_sp, _ = translate_ids(
                session, list(all_gs), "genesymbol", "swissprot", ncbi_tax_id,
            )

    # Step 4: Proteome filter (from DB reflist table)
    proteome = (
        _reflist_from_db(session, "swissprot", ncbi_tax_id)
        | _reflist_from_db(session, "trembl", ncbi_tax_id)
    )

    # Apply all steps to each result
    cleaned: dict[str, set[str]] = {}

    for src_id, targets in results.items():
        clean_set: set[str] = set()

        for ac in targets:
            # 1. Secondary -> primary
            primaries = sec_pri.get(ac)
            resolved = primaries if primaries else {ac}

            for rid in resolved:
                # 2. TrEMBL -> SwissProt
                if swissprot_set and rid not in swissprot_set:
                    gs = trembl_to_gs.get(rid, set())
                    found_sp = False
                    for g in gs:
                        sp = gs_to_sp.get(g, set())
                        if sp:
                            clean_set.update(sp)
                            found_sp = True
                            break
                    if not found_sp:
                        clean_set.add(rid)
                else:
                    clean_set.add(rid)

        # 3. Proteome filter
        if proteome:
            filtered = clean_set & proteome
            if filtered:
                clean_set = filtered

        # 4. Format validation
        clean_set = {u for u in clean_set if is_uniprot_ac(u)}

        cleaned[src_id] = clean_set if clean_set else targets

    return cleaned
