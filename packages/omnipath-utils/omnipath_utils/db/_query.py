"""Database query helpers."""

from __future__ import annotations

import logging
from functools import cache, lru_cache
from collections import defaultdict

from sqlalchemy import text
from sqlalchemy.orm import Session

from omnipath_utils.db._connection import SCHEMA

_log = logging.getLogger(__name__)


_FTP_TABLE = 'id_mapping_ftp'
_LONG_TABLE = 'id_mapping_long'
_ftp_exists_cache: dict = {}
_ftp_types_cache: dict = {}

# Long-value id_types -- names and structures -- live in ``id_mapping_long``
# (text values), not the ``varchar(64)`` database-ID ``id_mapping`` table (R2).
NAME_TYPES: frozenset[str] = frozenset(
    {'name', 'synonym', 'iupac', 'traditional_iupac'}
)
STRUCTURE_TYPES: frozenset[str] = frozenset({'inchi', 'smiles'})
LONG_VALUE_TYPES: frozenset[str] = NAME_TYPES | STRUCTURE_TYPES


def _is_long_query(source_type: str, target_type: str) -> bool:
    """A query is long-valued if either side is a name or structure type."""
    return source_type in LONG_VALUE_TYPES or target_type in LONG_VALUE_TYPES


# Database-ID id_types whose stored canonical value keeps an uppercase CURIE
# prefix (the "banana" problem -- ChEBI stores ``CHEBI:15377``, not ``15377``).
_BANANA_PREFIXES: dict[str, str] = {'chebi': 'CHEBI'}


@cache
def _accepted_prefixes(id_type: str) -> frozenset[str]:
    """Lowercase CURIE prefixes/aliases that may decorate an id of this type."""
    try:
        from omnipath_utils.mapping._id_types import IdTypeRegistry

        info = IdTypeRegistry.get().info(id_type) or {}
    except Exception:  # registry unavailable -> only the type name
        info = {}
    prefixes = {id_type.lower()}
    cp = info.get('curie_prefix')
    if cp:
        prefixes.add(str(cp).lower())
    for alias in info.get('aliases', []):
        prefixes.add(str(alias).lower())
    return frozenset(prefixes)


def strip_curie(id_type: str, identifier: str) -> str:
    """Normalise a CURIE-decorated identifier to the stored form for this type.

    Removes a leading ``<prefix>:`` when ``<prefix>`` (case-insensitive) is a
    known CURIE prefix / alias of ``id_type`` (so ``chebi:17612``, ``CHEBI:17612``
    and ``17612`` all collapse), then re-applies the canonical uppercase prefix
    for "banana" id_types (ChEBI). Identifiers without a matching prefix are
    returned unchanged.
    """
    s = str(identifier).strip()
    if ':' in s:
        head, _, rest = s.partition(':')
        if rest and head.lower() in _accepted_prefixes(id_type):
            s = rest.strip()
    banana = _BANANA_PREFIXES.get(id_type.lower())
    if banana and s and not s.upper().startswith(banana + ':'):
        s = f'{banana}:{s}'
    return s


def uniprot_organism(
    session: Session,
    accessions: list[str],
) -> dict[str, int | None]:
    """Map UniProt accessions to their source organism (NCBI taxon id).

    Reads the resolver's accession -> organism map, so it answers even for
    accessions with no gene/protein resolution (the TrEMBL long tail of
    less-studied organisms). Each accession is normalised (a CURIE prefix is
    stripped and an isoform ``-<n>`` suffix removed) before lookup. Returns a
    dict keyed by the ORIGINAL input; the value is the taxon id, or ``None`` when
    the accession is unknown.
    """
    import re

    from sqlalchemy import text

    normalised = {
        ac: re.sub(r'-[0-9]+$', '', strip_curie('uniprot', ac))
        for ac in accessions
    }
    keys = sorted({k for k in normalised.values() if k})

    found: dict[str, int] = {}
    if keys:
        rows = session.execute(
            text(
                f'SELECT uniprot, ncbi_tax_id FROM {SCHEMA}.resolver_uniprot_taxon '
                'WHERE uniprot = ANY(:keys)'
            ),
            {'keys': keys},
        ).all()
        found = {u: t for u, t in rows}

    return {ac: found.get(normalised[ac]) for ac in accessions}


def _lookup_key(source_type: str, identifier: str) -> str:
    """The normalised storage key for one query identifier.

    Names lowercased, structures verbatim, database IDs CURIE-stripped (+ HMDB
    digit-padded), InChIKey case- and prefix-normalised.
    """
    if source_type in NAME_TYPES:
        from omnipath_utils.mapping._id_types import normalize_name

        return normalize_name(identifier)
    if source_type == 'inchikey':
        # C8: InChIKey's own prefix convention is 'InChIKey=' (an '=', not
        # a ':') -- strip_curie's colon-only split never catches it, so this
        # goes through the real per-namespace normalizer instead. Falls back
        # to a plain strip+uppercase for a value that doesn't fully match
        # (e.g. malformed input) rather than dropping the query outright.
        from omnipath_utils.mapping._id_types import normalize_identifier

        normalized = normalize_identifier('inchikey', identifier)
        return normalized if normalized else str(identifier).strip().upper()
    if source_type == 'smiles':
        # C7: two equivalent SMILES spellings must match one stored key.
        # Degrades to a verbatim (byte) comparison without the chemistry
        # toolkit -- identifier-to-identifier translation stays unaffected
        # (research R14), only cross-spelling SMILES matching is lost.
        from omnipath_utils.mapping._chemistry import (
            canonicalize_smiles,
            chemistry_available,
        )

        raw = str(identifier).strip()
        if chemistry_available():
            canonical = canonicalize_smiles(raw)
            if canonical:
                return canonical
        return raw
    if source_type in STRUCTURE_TYPES:
        return str(identifier).strip()
    s = strip_curie(source_type, identifier)
    if source_type == 'hmdb':
        from omnipath_utils.mapping._special import normalise_hmdb
        s = normalise_hmdb(s)
    return s


def _rekey(
    result: dict,
    identifiers: list[str] | None,
    source_type: str,
) -> dict[str, set[str]]:
    """Map a result keyed by normalised lookup keys back to the original input.

    So the response echoes the identifier the user actually sent.
    """
    if identifiers is None:
        return dict(result)
    out: dict[str, set[str]] = {}
    for orig in identifiers:
        hits = result.get(_lookup_key(source_type, orig))
        if hits:
            out.setdefault(orig, set()).update(hits)
    return out


def _ftp_table_exists(session: Session) -> bool:
    """Whether the comprehensive FTP table is present (cached per process)."""
    if 'v' not in _ftp_exists_cache:
        present = session.execute(
            text(f"SELECT to_regclass('{SCHEMA}.{_FTP_TABLE}')")
        ).scalar()
        _ftp_exists_cache['v'] = present is not None
    return _ftp_exists_cache['v']


def _ftp_relevant(
    source_type: str,
    target_type: str,
    ftp_types: frozenset[str],
) -> bool:
    """FR-018 gate predicate for the full-UniProt fallback.

    ``id_mapping_ftp`` can only contain a mapping when **both** sides are
    id_types present in that table (it is 100% UniProt-family). A pair with a
    non-FTP type on either side has zero rows there, so consulting it is
    fruitless -- skipping is result-identical (R8).
    """
    return source_type in ftp_types and target_type in ftp_types


def _ftp_types(session: Session) -> frozenset[str]:
    """The set of id_types present in ``id_mapping_ftp`` (cached per process).

    Read cheaply from a precomputed ``build_info`` row (``table_name =
    'ftp_types'``); if that is absent (a build predating FR-018) it is derived
    once from the table's distinct type ids and cached. Empty when the FTP
    table is not present.
    """
    if 'v' in _ftp_types_cache:
        return _ftp_types_cache['v']

    if not _ftp_table_exists(session):
        _ftp_types_cache['v'] = frozenset()
        return _ftp_types_cache['v']

    # Precomputed by DatabaseBuilder.record_ftp_types as one row per type
    # (build_info.source_type is varchar(64), so a row per name, not a joined
    # string).
    precomputed = session.execute(
        text(
            f"SELECT source_type FROM {SCHEMA}.build_info "
            "WHERE table_name = 'ftp_types' AND source_type IS NOT NULL"
        )
    ).scalars().all()
    if precomputed:
        types = frozenset(t for t in precomputed if t)
    else:
        _log.warning(
            'ftp_types not precomputed in build_info; deriving once from '
            '%s.%s (run DatabaseBuilder.record_ftp_types to cache it)',
            SCHEMA,
            _FTP_TABLE,
        )
        rows = session.execute(
            text(
                f'SELECT name FROM {SCHEMA}.id_type WHERE id IN ('
                f' SELECT DISTINCT source_type_id FROM {SCHEMA}.{_FTP_TABLE}'
                f' UNION'
                f' SELECT DISTINCT target_type_id FROM {SCHEMA}.{_FTP_TABLE})'
            )
        ).scalars().all()
        types = frozenset(rows)

    _ftp_types_cache['v'] = types
    return types


def _query_table(
    session: Session,
    table: str,
    identifiers: list[str] | None,
    source_type: str,
    target_type: str,
    ncbi_tax_id: int,
) -> tuple[defaultdict, set]:
    """Forward + reverse-for-missing lookup against a single mapping table."""
    result = defaultdict(set)
    backends_used = set()
    params: dict = {
        'src_type': source_type,
        'tgt_type': target_type,
        'tax': ncbi_tax_id,
    }
    id_filter = ''
    if identifiers is not None:
        id_filter = 'AND m.source_id = ANY(:ids)'
        params['ids'] = identifiers

    rows = session.execute(
        text(f"""
            SELECT m.source_id, m.target_id, b.name
            FROM {table} m
            JOIN {SCHEMA}.id_type st ON m.source_type_id = st.id
            JOIN {SCHEMA}.id_type tt ON m.target_type_id = tt.id
            JOIN {SCHEMA}.backend b ON m.backend_id = b.id
            WHERE st.name = :src_type
            AND tt.name = :tgt_type
            AND m.ncbi_tax_id = :tax
            {id_filter}
        """),
        params,
    )
    for row in rows:
        result[row[0]].add(row[1])
        backends_used.add(row[2])

    if identifiers is not None:
        missing = [i for i in identifiers if i not in result]
        if missing:
            rev_rows = session.execute(
                text(f"""
                    SELECT m.target_id, m.source_id, b.name
                    FROM {table} m
                    JOIN {SCHEMA}.id_type st ON m.source_type_id = st.id
                    JOIN {SCHEMA}.id_type tt ON m.target_type_id = tt.id
                    JOIN {SCHEMA}.backend b ON m.backend_id = b.id
                    WHERE st.name = :tgt_type
                    AND tt.name = :src_type
                    AND (m.ncbi_tax_id = :tax OR m.ncbi_tax_id = 0)
                    AND m.target_id = ANY(:ids)
                """),
                {
                    'src_type': source_type,
                    'tgt_type': target_type,
                    'tax': ncbi_tax_id,
                    'ids': missing,
                },
            )
            for row in rev_rows:
                result[row[0]].add(row[1])
                backends_used.add(f'{row[2]}(rev)')

    return result, backends_used


# ---------------------------------------------------------------------------
# The generic 'name' axis (spec 011 FR-016..021): 'name', 'synonym', 'iupac'
# and 'traditional_iupac' stop being unrelated types the caller must each try
# in turn, and become one graded lookup. A name is what the minting resource
# recommends; every other name-like value it records -- including a
# systematic IUPAC name -- is a synonym for this purpose (FR-016).
# ---------------------------------------------------------------------------

_NAME_PREFERRED_TYPES: frozenset[str] = frozenset({'name'})
_NAME_SYNONYM_TYPES: frozenset[str] = frozenset(
    {'synonym', 'iupac', 'traditional_iupac'}
)


def _name_axis_buckets(
    name_scope: str, forward: bool,
) -> list[tuple[str, frozenset[str]]]:
    """The ordered (label, real id_types) buckets to consult for one query.

    ``forward`` is True when the 'name' side is the *source* (translating
    FROM a name, FR-017: preferred first, fall back to synonyms); False when
    it is the *target* (translating TO a name, FR-018: preferred only by
    default, no fallback).
    """
    if name_scope == 'preferred':
        return [('preferred', _NAME_PREFERRED_TYPES)]
    if name_scope == 'synonym':
        return [('synonym', _NAME_SYNONYM_TYPES)]
    if name_scope == 'any':
        return [('any', _NAME_PREFERRED_TYPES | _NAME_SYNONYM_TYPES)]
    if forward:
        return [
            ('preferred', _NAME_PREFERRED_TYPES),
            ('synonym', _NAME_SYNONYM_TYPES),
        ]
    return [('preferred', _NAME_PREFERRED_TYPES)]


def _query_name_axis(
    session: Session,
    table: str,
    identifiers: list[str] | None,
    forward: bool,
    other_type: str,
    ncbi_tax_id: int,
    name_scope: str,
) -> tuple[defaultdict[str, set[str]], set[str], dict[str, str]]:
    """Query the generic 'name' axis: union name/synonym/iupac/traditional_iupac,
    honoring ``name_scope``, and report which bucket answered each identifier.

    Each real type is queried through :func:`_query_table` (one call per
    type), reusing its existing forward + reverse-for-missing lookup rather
    than reimplementing it -- a name-like type is stored in only one
    direction for most cross-referenced ids (R3), and the reverse fallback is
    what makes those still answer.
    """
    result: defaultdict[str, set[str]] = defaultdict(set)
    matched_as: dict[str, str] = {}
    backends_used: set[str] = set()
    remaining = identifiers

    for bucket_label, types in _name_axis_buckets(name_scope, forward):
        if remaining is not None and not remaining:
            break
        for real_type in sorted(types):
            src_type, tgt_type = (
                (real_type, other_type) if forward else (other_type, real_type)
            )
            rows, backends = _query_table(
                session, table, remaining, src_type, tgt_type, ncbi_tax_id,
            )
            for ident, targets in rows.items():
                if not targets:
                    continue
                result[ident] |= targets
                matched_as.setdefault(ident, bucket_label)
            backends_used |= backends
        if remaining is not None:
            remaining = [i for i in remaining if i not in result]

    return result, backends_used, matched_as


# ---------------------------------------------------------------------------
# Deprecated / aged-ID recovery (007 US2, FR-006..009). A post-primary fallback:
# only still-unmapped ids of a recoverable source_type consult the recovery tables,
# so the valid-id common path is untouched (FR-007).
# ---------------------------------------------------------------------------

# query source_type -> [(recovery source_type, recovery target_type, label, tax0)].
# The recovered current id is in the recovery target_type namespace; a
# self-referential row (target_id == source_id) means the id is DELETED (no
# successor). ``tax0`` marks organism-agnostic recovery tables (loaded at tax 0).
_RECOVERY: dict[str, list[tuple[str, str, str, bool]]] = {
    'uniprot': [
        ('uniprot-sec', 'uniprot-pri', 'uniprot_sec', True),
        ('uniprot-deleted', 'uniprot-deleted', 'uniprot_deleted', True),
    ],
    'uniprot-pri': [
        ('uniprot-sec', 'uniprot-pri', 'uniprot_sec', True),
        ('uniprot-deleted', 'uniprot-deleted', 'uniprot_deleted', True),
    ],
    'entrez': [
        ('entrez-history', 'entrez', 'gene_history', False),
        # self-referential rows = discontinued with no successor -> deleted
        ('entrez-history', 'entrez-history', 'gene_history', False),
    ],
    'ensg': [
        ('ensembl-history', 'ensg', 'ensembl_history', False),
        # self-referential rows = retired with no successor -> deleted
        ('ensembl-history', 'ensembl-history', 'ensembl_history', False),
    ],
}

# Collapse an id_type to the namespace used for "same-namespace" recovery returns.
_RECOVERY_NS = {
    'uniprot-pri': 'uniprot',
    'uniprot-sec': 'uniprot',
    'uniprot-deleted': 'uniprot',
}


def _ns(id_type: str) -> str:
    return _RECOVERY_NS.get(id_type, id_type)


def _recover_query(
    session: Session,
    ids: list[str],
    source_type: str,
    target_type: str,
    ncbi_tax_id: int,
) -> defaultdict:
    """Forward-only deprecated -> current lookup (no reverse, unlike _query_table)."""
    result = defaultdict(set)
    rows = session.execute(
        text(f"""
            SELECT m.source_id, m.target_id
            FROM {SCHEMA}.id_mapping m
            JOIN {SCHEMA}.id_type st ON m.source_type_id = st.id
            JOIN {SCHEMA}.id_type tt ON m.target_type_id = tt.id
            WHERE st.name = :src AND tt.name = :tgt
            AND (m.ncbi_tax_id = :tax OR m.ncbi_tax_id = 0)
            AND m.source_id = ANY(:ids)
        """),
        {'src': source_type, 'tgt': target_type, 'tax': ncbi_tax_id, 'ids': ids},
    )
    for src, tgt in rows:
        result[src].add(tgt)
    return result


def _recover(
    session: Session,
    missing: list[str],
    source_type: str,
    target_type: str,
    ncbi_tax_id: int,
    result: defaultdict,
    meta: dict,
) -> set[str]:
    """Consult the recovery tables for still-missing ids of a recoverable type.

    Folds recovered current targets into ``result`` and records per-id flags
    (``recovered`` / ``recovery_source`` / ``ambiguous`` / ``deleted``) in ``meta``
    (keyed by the normalised lookup id). Returns the backend labels used. Deprecated
    ids that recover to >1 current identifier return **all** candidates flagged
    ``ambiguous`` (clarification); a deleted id with no successor is flagged
    ``deleted`` rather than left empty (FR-008/009).
    """
    specs = _RECOVERY.get(source_type)
    backends: set[str] = set()
    if not specs:
        return backends
    still = list(missing)
    for rec_src, rec_tgt, label, tax0 in specs:
        if not still:
            break
        tax = 0 if tax0 else ncbi_tax_id
        rec_map = _recover_query(session, still, rec_src, rec_tgt, tax)
        if not rec_map:
            continue
        backends.add(label)
        for src_id, currents in rec_map.items():
            # Self-referential row -> the id is explicitly deleted (no successor).
            if rec_src == rec_tgt and currents == {src_id}:
                meta[src_id] = {
                    'recovered': False, 'deleted': True,
                    'recovery_source': label, 'ambiguous': False,
                }
                continue
            cur_ns = _ns(rec_tgt)
            if _ns(target_type) == cur_ns:
                targets = set(currents)
            else:
                # Re-translate the recovered current id to the requested target
                # via the PRIMARY route (no nested recovery).
                sub, sub_b = translate_ids(
                    session, sorted(currents), cur_ns, target_type,
                    ncbi_tax_id, recover=False,
                )
                backends |= sub_b
                targets = set().union(*sub.values()) if sub else set()
            if targets:
                result[src_id] |= targets
                meta[src_id] = {
                    'recovered': True, 'deleted': False,
                    'recovery_source': label, 'ambiguous': len(currents) > 1,
                }
        still = [i for i in still if i not in result and i not in meta]
    return backends


def translate_ids(
    session: Session,
    identifiers: list[str] | None,
    source_type: str,
    target_type: str,
    ncbi_tax_id: int,
    full_uniprot: str = 'fallback',
    recover: bool = False,
    recovery_meta: dict | None = None,
    name_scope: str = 'default',
    name_meta: dict | None = None,
) -> tuple[dict[str, set[str]], set[str]]:
    """Translate IDs via the database.

    Queries the curated ``id_mapping`` first and, per ``full_uniprot``, the
    comprehensive full-UniProt table ``id_mapping_ftp`` (the complete UniProt
    idmapping, all organisms). Results are deduplicated (set union).

    Args:
        session: SQLAlchemy session.
        identifiers: Source IDs to translate, or ``None`` for the full table.
        source_type / target_type: ID type names.
        ncbi_tax_id: NCBI Taxonomy ID.
        full_uniprot: How to use the comprehensive full-UniProt table:
            ``'fallback'`` (default — curated, then full-UniProt only for
            still-unresolved identifiers), ``'never'`` (curated only),
            ``'both'`` (curated + full-UniProt), ``'only'`` (full-UniProt only).
            The full table is only consulted on explicit request
            (``both``/``only``) for a whole-table query, never as a blanket
            fallback.
        name_scope: Only consulted when ``source_type`` or ``target_type`` is
            the generic ``'name'`` axis (FR-016..018): ``'default'``
            (preferred first, fall back to synonyms when translating from a
            name; preferred only when translating to one), ``'preferred'``,
            ``'synonym'`` or ``'any'``.
        name_meta: When given and the name axis is queried, populated with
            ``{original_identifier: 'preferred' | 'synonym' | 'any'}`` naming
            the bucket that answered each input (FR-020's ``matched_as``).

    Returns:
        Tuple of (results dict, set of backend names used).
    """
    # Normalise every query identifier to its storage key once: names
    # lowercased, structures verbatim, database IDs CURIE-stripped (chebi:17612,
    # CHEBI:17612, 17612 all collapse) + HMDB digit-padded. The response is
    # re-keyed back to the caller's original inputs at the end.
    norm_ids = (
        [_lookup_key(source_type, i) for i in identifiers]
        if identifiers is not None else None
    )

    # The generic 'name' axis (FR-016..018): union name/synonym/iupac/
    # traditional_iupac behind one lookup instead of four unrelated types.
    if source_type == 'name' or target_type == 'name':
        forward = source_type == 'name'
        other_type = target_type if forward else source_type
        # Chemicals are organism-agnostic (rows at tax 0), matching the
        # plain long-value path below.
        result, backends_used, matched_as = _query_name_axis(
            session, f'{SCHEMA}.{_LONG_TABLE}', norm_ids, forward,
            other_type, 0, name_scope,
        )
        if name_meta is not None:
            for orig in identifiers or []:
                bucket = matched_as.get(_lookup_key(source_type, orig))
                if bucket:
                    name_meta[orig] = bucket
        return _rekey(result, identifiers, source_type), backends_used

    # Long-value (name / structure) queries route to the separate
    # ``id_mapping_long`` table (R2). Chemicals are organism-agnostic (rows at
    # tax 0), and this path never touches the UniProt FTP table (FR-018 holds
    # trivially).
    if _is_long_query(source_type, target_type):
        result, backends_used = _query_table(
            session, f'{SCHEMA}.{_LONG_TABLE}',
            norm_ids, source_type, target_type, 0,
        )
        # C6: an InChI/SMILES with no stored match still translates to its
        # structure key -- computed live rather than requiring it to
        # already be one of our own stored source values. Only whatever
        # the stored lookup missed; a no-op without the chemistry toolkit
        # (compute_inchikey degrades to None itself, research R14).
        if target_type == 'inchikey' and source_type in STRUCTURE_TYPES:
            missing = [i for i in (norm_ids or []) if i not in result]
            if missing:
                from omnipath_utils.mapping._chemistry import compute_inchikey

                kwarg = 'inchi' if source_type == 'inchi' else 'smiles'
                for ident in missing:
                    key = compute_inchikey(**{kwarg: ident})
                    if key:
                        result[ident] = {key}
                        backends_used = backends_used | {'chemistry_toolkit'}
        return _rekey(result, identifiers, source_type), backends_used

    if full_uniprot == 'only':
        result, backends_used = defaultdict(set), set()
    else:
        result, backends_used = _query_table(
            session, f'{SCHEMA}.id_mapping',
            norm_ids, source_type, target_type, ncbi_tax_id,
        )

    # Decide whether (and for which identifiers) to consult the full table.
    ftp_ids = norm_ids
    want_ftp = full_uniprot in ('both', 'only')
    if full_uniprot == 'fallback' and norm_ids is not None:
        missing = [i for i in norm_ids if i not in result]
        if missing:
            want_ftp = True
            ftp_ids = missing

    # FR-018 gate: the ~744 M-row FTP table is 100% UniProt-family, so a pair
    # with a non-FTP type on either side (e.g. ``name → chebi``, ``chebi →
    # hmdb``) has zero rows there. Skip the fruitless scan -- result-identical
    # to the prior unconditional fallback, a pure latency win (R8).
    if want_ftp and not _ftp_relevant(
        source_type, target_type, _ftp_types(session)
    ):
        want_ftp = False

    if want_ftp and _ftp_table_exists(session):
        ftp_result, ftp_backends = _query_table(
            session, f'{SCHEMA}.{_FTP_TABLE}',
            ftp_ids, source_type, target_type, ncbi_tax_id,
        )
        for key, vals in ftp_result.items():
            result[key] |= vals  # set union -> deduplicated
        backends_used |= ftp_backends

    # Post-primary deprecated-ID recovery (US2): only still-unmapped ids of a
    # recoverable source_type consult the recovery tables (FR-007 — the valid path
    # is untouched). Off by default; the /mapping/translate route opts in.
    if recover and norm_ids is not None and source_type in _RECOVERY:
        still_missing = [i for i in norm_ids if not result.get(i)]
        if still_missing:
            _meta: dict = {}
            rec_backends = _recover(
                session, still_missing, source_type, target_type,
                ncbi_tax_id, result, _meta,
            )
            backends_used |= rec_backends
            if recovery_meta is not None and _meta:
                # Re-key recovery flags from normalised keys back to the caller's
                # original inputs (mirrors _rekey for the result).
                for orig in identifiers or []:
                    m = _meta.get(_lookup_key(source_type, orig))
                    if m:
                        recovery_meta[orig] = m

    return _rekey(result, identifiers, source_type), backends_used


def get_full_table(
    session: Session,
    source_type: str,
    target_type: str,
    ncbi_tax_id: int,
) -> dict[str, set[str]]:
    """Get a full mapping table from the database."""

    result, _ = translate_ids(session, None, source_type, target_type, ncbi_tax_id)

    return result


@lru_cache(maxsize=1)
def _all_known_prefixes() -> frozenset[str]:
    """Union of accepted CURIE prefixes across every id_type (cached)."""
    try:
        from omnipath_utils.mapping._id_types import IdTypeRegistry

        reg = IdTypeRegistry.get()
        out: set[str] = set()
        for name in reg.all_names():
            out |= _accepted_prefixes(name)
        return frozenset(out)
    except Exception:
        return frozenset()


def _identify_candidates(identifier: str) -> list[str]:
    """Candidate stored forms to search for an identifier of unknown type.

    The raw value, its lowercase (names), and -- when it carries a known CURIE
    prefix -- the prefix-stripped value (and the ChEBI banana form).
    """
    s = str(identifier).strip()
    cands = {s, s.lower()}
    if ':' in s:
        head, _, rest = s.partition(':')
        rest = rest.strip()
        if rest and head.lower() in _all_known_prefixes():
            cands.add(rest)
            if head.lower() == 'chebi':
                cands.add(f'CHEBI:{rest}')
    return [c for c in cands if c]


def identify_ids(
    session: Session,
    identifiers: list[str],
    ncbi_tax_id: int,
) -> dict[str, list[dict]]:
    """Identify what type(s) each identifier belongs to.

    Searches the id_mapping and id_mapping_long tables as both source and
    target to find which ID types contain each identifier. CURIE-decorated
    inputs (``chebi:17612``) are matched against their stored form.

    Returns:
        Dict mapping each identifier to a list of dicts with
        "id_type", "role" ("source" or "target"), and "count".
    """
    result = {}

    for identifier in identifiers:
        matches = []
        cands = _identify_candidates(identifier)

        # Search as source
        rows = session.execute(
            text(f"""
                SELECT st.name, count(DISTINCT m.target_id)
                FROM {SCHEMA}.id_mapping m
                JOIN {SCHEMA}.id_type st ON m.source_type_id = st.id
                WHERE m.source_id = ANY(:ids)
                AND (m.ncbi_tax_id = :tax OR m.ncbi_tax_id = 0)
                GROUP BY st.name
            """),
            {'ids': cands, 'tax': ncbi_tax_id},
        ).fetchall()

        for row in rows:
            matches.append({
                'id_type': row[0],
                'role': 'source',
                'mappings_count': row[1],
            })

        # Search as target
        rows = session.execute(
            text(f"""
                SELECT tt.name, count(DISTINCT m.source_id)
                FROM {SCHEMA}.id_mapping m
                JOIN {SCHEMA}.id_type tt ON m.target_type_id = tt.id
                WHERE m.target_id = ANY(:ids)
                AND (m.ncbi_tax_id = :tax OR m.ncbi_tax_id = 0)
                GROUP BY tt.name
            """),
            {'ids': cands, 'tax': ncbi_tax_id},
        ).fetchall()

        for row in rows:
            matches.append({
                'id_type': row[0],
                'role': 'target',
                'mappings_count': row[1],
            })

        # Also consult the long-value table (names / structures).
        long_rows = session.execute(
            text(f"""
                SELECT st.name, 'source' AS role, count(DISTINCT m.target_id)
                FROM {SCHEMA}.{_LONG_TABLE} m
                JOIN {SCHEMA}.id_type st ON m.source_type_id = st.id
                WHERE m.source_id = ANY(:ids)
                GROUP BY st.name
                UNION ALL
                SELECT tt.name, 'target' AS role, count(DISTINCT m.source_id)
                FROM {SCHEMA}.{_LONG_TABLE} m
                JOIN {SCHEMA}.id_type tt ON m.target_type_id = tt.id
                WHERE m.target_id = ANY(:ids)
                GROUP BY tt.name
            """),
            {'ids': cands},
        ).fetchall()

        for row in long_rows:
            matches.append({
                'id_type': row[0],
                'role': row[1],
                'mappings_count': row[2],
            })

        result[identifier] = matches

    return result


def get_all_mappings(
    session: Session,
    identifiers: list[str],
    source_type: str,
    ncbi_tax_id: int,
) -> dict[str, dict[str, list[str]]]:
    """Get all mappings for identifiers across all target types.

    Returns:
        Dict mapping each identifier to a dict of {target_type: [target_ids]}.
    """
    result: dict[str, dict[str, list[str]]] = {}

    # Normalise each identifier to its storage key (CURIE-stripped for db-ids,
    # lowercased for names) and map it back to the caller's original input.
    norm: dict[str, str] = {}
    for ident in identifiers:
        norm.setdefault(_lookup_key(source_type, ident), ident)
    norm_ids = list(norm.keys())

    rows = session.execute(
        text(f"""
            SELECT m.source_id, tt.name, m.target_id
            FROM {SCHEMA}.id_mapping m
            JOIN {SCHEMA}.id_type st ON m.source_type_id = st.id
            JOIN {SCHEMA}.id_type tt ON m.target_type_id = tt.id
            WHERE st.name = :src_type
            AND (m.ncbi_tax_id = :tax OR m.ncbi_tax_id = 0)
            AND m.source_id = ANY(:ids)
            ORDER BY m.source_id, tt.name
        """),
        {'src_type': source_type, 'tax': ncbi_tax_id, 'ids': norm_ids},
    ).fetchall()

    for src_id, tgt_type, tgt_id in rows:
        orig = norm.get(src_id, src_id)
        result.setdefault(orig, {}).setdefault(tgt_type, []).append(tgt_id)

    # Union the long-value table (names / structures) for long source types.
    if source_type in LONG_VALUE_TYPES:
        long_rows = session.execute(
            text(f"""
                SELECT m.source_id, tt.name, m.target_id
                FROM {SCHEMA}.{_LONG_TABLE} m
                JOIN {SCHEMA}.id_type st ON m.source_type_id = st.id
                JOIN {SCHEMA}.id_type tt ON m.target_type_id = tt.id
                WHERE st.name = :src_type
                AND m.source_id = ANY(:ids)
                ORDER BY m.source_id, tt.name
            """),
            {'src_type': source_type, 'ids': norm_ids},
        ).fetchall()
        for src_id, tgt_type, tgt_id in long_rows:
            orig = norm.get(src_id, src_id)
            bucket = result.setdefault(orig, {}).setdefault(tgt_type, [])
            if tgt_id not in bucket:
                bucket.append(tgt_id)

    return result


def chain_translate(
    session: Session,
    identifiers: list[str],
    source_type: str,
    target_type: str,
    ncbi_tax_id: int,
    via: str = 'uniprot',
) -> tuple[dict[str, set[str]], set[str]]:
    """Chain translation: source -> via -> target, all in the DB.

    Uses two JOINs to find: source_id -> via_id -> target_id.
    Returns (result_dict, backends_used).
    """
    from collections import defaultdict

    result = defaultdict(set)
    backends_used = set()

    rows = session.execute(
        text(f"""
            SELECT m1.source_id, m2.target_id, b1.name, b2.name
            FROM {SCHEMA}.id_mapping m1
            JOIN {SCHEMA}.id_type st1 ON m1.source_type_id = st1.id
            JOIN {SCHEMA}.id_type vt ON m1.target_type_id = vt.id
            JOIN {SCHEMA}.id_mapping m2 ON m1.target_id = m2.source_id
                AND m2.ncbi_tax_id = m1.ncbi_tax_id
            JOIN {SCHEMA}.id_type st2 ON m2.source_type_id = st2.id
            JOIN {SCHEMA}.id_type tt2 ON m2.target_type_id = tt2.id
            JOIN {SCHEMA}.backend b1 ON m1.backend_id = b1.id
            JOIN {SCHEMA}.backend b2 ON m2.backend_id = b2.id
            WHERE st1.name = :src_type
            AND vt.name = :via_type
            AND st2.name = :via_type
            AND tt2.name = :tgt_type
            AND m1.source_id = ANY(:ids)
            AND (m1.ncbi_tax_id = :tax OR m1.ncbi_tax_id = 0)
            AND (m2.ncbi_tax_id = :tax OR m2.ncbi_tax_id = 0)
        """),
        {
            'src_type': source_type,
            'via_type': via,
            'tgt_type': target_type,
            'tax': ncbi_tax_id,
            'ids': identifiers,
        },
    ).fetchall()

    for row in rows:
        result[row[0]].add(row[1])
        backends_used.add(f'{row[2]}+{row[3]}')

    return dict(result), backends_used
