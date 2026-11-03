"""Mapping (ID translation) endpoints."""

from __future__ import annotations

import logging

from litestar import Controller, get, post
from sqlalchemy.orm import Session
from litestar.params import Parameter

from omnipath_utils.db._query import (
    identify_ids,
    translate_ids,
    chain_translate,
    get_all_mappings,
    uniprot_organism,
)
from omnipath_utils.mapping._chemistry import chemistry_available

_log = logging.getLogger(__name__)


def _resolve_and_cleanup(
    result: dict[str, set[str]],
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int,
    session=None,
) -> dict[str, set[str]]:
    """Apply alias resolution and UniProt cleanup to results."""

    from omnipath_utils.mapping._id_types import IdTypeRegistry

    reg = IdTypeRegistry.get()
    target_resolved = reg.resolve(target_id_type) or target_id_type

    if target_resolved == 'uniprot':
        try:
            if session:
                # Batch cleanup via DB (vectorized — O(1) queries)
                from omnipath_utils.mapping._cleanup import uniprot_cleanup_batch
                result = uniprot_cleanup_batch(result, ncbi_tax_id, session=session)
            else:
                # Per-ID cleanup via in-memory mapper (fallback)
                from omnipath_utils.mapping._mapper import Mapper
                from omnipath_utils.mapping._cleanup import uniprot_cleanup

                mapper = Mapper.get()
                for src_id in result:
                    if result[src_id]:
                        result[src_id] = uniprot_cleanup(
                            result[src_id], ncbi_tax_id, mapper=mapper,
                        )
        except Exception:
            _log.debug(
                'UniProt cleanup failed for target %s, returning raw results.',
                target_id_type,
                exc_info=True,
            )

    return result


def _apply_fallbacks(
    session: Session,
    id_list: list[str],
    id_type_resolved: str,
    target_resolved: str,
    ncbi_tax_id: int,
    result: dict[str, set[str]],
) -> dict[str, set[str]]:
    """Apply the same fallback chain as the Python API via DB queries.

    Fallbacks:
    1. Gene symbol: uppercase, capitalize, synonyms
    2. Chain translation via uniprot
    """

    # Gene symbol fallbacks
    if id_type_resolved in ('genesymbol', 'genesymbol-syn'):
        unmapped = [i for i in id_list if not result.get(i)]

        if unmapped:
            # Try uppercase
            upper_map = {}
            upper_ids = []
            for name in unmapped:
                u = name.upper()
                if u != name:
                    upper_map[u] = name
                    upper_ids.append(u)

            if upper_ids:
                upper_result, _ub = translate_ids(
                    session,
                    upper_ids,
                    id_type_resolved,
                    target_resolved,
                    ncbi_tax_id,
                )
                for upper, targets in upper_result.items():
                    if targets and upper in upper_map:
                        result[upper_map[upper]] = targets

            # Try capitalized (for rodent symbols)
            still_unmapped = [i for i in unmapped if not result.get(i)]
            if still_unmapped:
                cap_map = {}
                cap_ids = []
                for name in still_unmapped:
                    c = name.capitalize()
                    if c != name:
                        cap_map[c] = name
                        cap_ids.append(c)

                if cap_ids:
                    cap_result, _ = translate_ids(
                        session,
                        cap_ids,
                        id_type_resolved,
                        target_resolved,
                        ncbi_tax_id,
                    )
                    for cap, targets in cap_result.items():
                        if targets and cap in cap_map:
                            result[cap_map[cap]] = targets

            # Try synonyms
            still_unmapped = [i for i in unmapped if not result.get(i)]
            if still_unmapped:
                syn_result, _sb = translate_ids(
                    session,
                    still_unmapped,
                    'genesymbol-syn',
                    target_resolved,
                    ncbi_tax_id,
                )
                for src, targets in syn_result.items():
                    if targets:
                        result[src] = targets

                # Also try uppercase synonyms
                still_unmapped = [
                    i for i in still_unmapped if not result.get(i)
                ]
                if still_unmapped:
                    upper_syns = [n.upper() for n in still_unmapped]
                    syn_upper, _ = translate_ids(
                        session,
                        upper_syns,
                        'genesymbol-syn',
                        target_resolved,
                        ncbi_tax_id,
                    )
                    for idx, name in enumerate(still_unmapped):
                        upper = upper_syns[idx]
                        if syn_upper.get(upper):
                            result[name] = syn_upper[upper]

    # Chain via uniprot if neither side is uniprot
    still_unmapped = [i for i in id_list if not result.get(i)]
    if (
        still_unmapped
        and id_type_resolved != 'uniprot'
        and target_resolved != 'uniprot'
    ):
        chain_result, chain_backends = chain_translate(
            session,
            still_unmapped,
            id_type_resolved,
            target_resolved,
            ncbi_tax_id,
            via='uniprot',
        )
        for src, targets in chain_result.items():
            if targets:
                result[src] = targets

    return result


#: FR-021: a preferred-name match ranks above a synonym match, and — within
#: the same bucket — a name mapping to one structure ranks above one mapping
#: to several. ``name_scope=any`` is documented as unranked (contracts/
#: translation-api.md), so this ordering applies to every other scope.
_NAME_BUCKET_RANK = {'preferred': 0, 'synonym': 1}


def _rank_name_matches(
    mapped: dict[str, list[str]], matched_as: dict[str, str],
) -> dict[str, list[str]]:
    ranked_keys = sorted(
        mapped,
        key=lambda k: (
            _NAME_BUCKET_RANK.get(matched_as.get(k, 'synonym'), 1),
            len(mapped[k]),
            k,
        ),
    )
    return {k: mapped[k] for k in ranked_keys}


def _needs_chemistry_toolkit(id_type: str, target_id_type: str) -> bool:
    """Whether this (id_type, target_id_type) pair depends on the live
    chemistry toolkit (spec 011 R14) -- SMILES canonicalization on either
    side (C7), or the InChI/SMILES -> InChIKey live-compute fallback (C6).
    A plain InChIKey<->id lookup is stored-data only and does not depend
    on it (research: "identifier-to-identifier translation is
    unaffected").
    """

    return 'smiles' in (id_type, target_id_type) or (
        id_type in ('inchi', 'smiles') and target_id_type == 'inchikey'
    )


def _structure_key_capability_meta() -> dict:
    available = chemistry_available()
    return {
        'available': available,
        'reason': None if available else (
            'the chemistry toolkit (omnipath-utils[chem]) is not '
            'installed in this server process'
        ),
    }


def _build_translate_response(
    id_list: list[str],
    result: dict[str, set[str]],
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int,
    raw: bool,
    backends_used: set[str],
    backend: str | None,
    loading: bool,
    recovery: dict | None = None,
    name_scope: str | None = None,
    matched_as: dict | None = None,
    structure_capability: dict | None = None,
) -> dict:
    """Build the response dict for translate endpoints."""

    mapped = {k: sorted(v) for k, v in result.items() if v}
    if matched_as and name_scope != 'any':
        mapped = _rank_name_matches(mapped, matched_as)
    unmapped = [i for i in id_list if i not in result or not result[i]]

    meta = {
        'id_type': id_type,
        'target_id_type': target_id_type,
        'ncbi_tax_id': ncbi_tax_id,
        'total_input': len(id_list),
        'total_mapped': len(mapped),
        'raw': raw,
        'backend': sorted(backends_used) if backends_used else backend,
        'loading': loading,
    }

    if name_scope is not None:
        meta['name_scope'] = name_scope
        meta['matched_as'] = {
            k: v for k, v in (matched_as or {}).items() if k in mapped
        }

    if structure_capability is not None:
        meta['structure_key_computation'] = structure_capability

    if loading:
        meta['loading_note'] = (
            f'Table {id_type} -> {target_id_type} is being loaded. '
            'Try again in a few minutes.'
        )

    response = {
        'results': mapped,
        'unmapped': unmapped,
        'meta': meta,
    }

    # Deprecated-ID recovery flags (007 US2), a parallel map keyed by input id:
    # {recovered, recovery_source, ambiguous, deleted}. Present only when the
    # recovery stage flagged at least one input (an id the primary route missed).
    if recovery:
        response['recovery'] = recovery
        meta['total_recovered'] = sum(
            1 for m in recovery.values() if m.get('recovered')
        )
        meta['total_deleted'] = sum(
            1 for m in recovery.values() if m.get('deleted')
        )

    return response


def _maybe_trigger_load(
    result: dict,
    raw: bool,
    id_type_resolved: str,
    target_resolved: str,
    ncbi_tax_id: int,
) -> bool:
    """Trigger background loading if no results found. Returns loading flag."""

    if result or raw:
        return False

    from omnipath_utils.db._loader import is_pending, request_table
    from omnipath_utils.db._connection import get_db_url

    if is_pending(id_type_resolved, target_resolved, ncbi_tax_id):
        return True

    queued = request_table(
        id_type_resolved,
        target_resolved,
        ncbi_tax_id,
        db_url=get_db_url(),
    )
    return queued


class MappingController(Controller):
    """ID translation endpoints."""

    path = '/mapping'

    @get('/translate')
    async def translate(  # noqa: D417
        self,
        session: Session,
        identifiers: str = Parameter(
            description='Comma-separated identifiers',
        ),
        id_type: str = Parameter(description='Source ID type'),
        target_id_type: str = Parameter(description='Target ID type'),
        ncbi_tax_id: int = Parameter(
            default=9606,
            description='NCBI Taxonomy ID',
        ),
        raw: bool = Parameter(
            default=False,
            description='Skip special-case handling, return raw DB results',
        ),
        backend: str | None = Parameter(
            default=None,
            required=False,
            description='Force specific backend',
        ),
        full_uniprot: str = Parameter(
            default='fallback',
            description=(
                'Use of the comprehensive full-UniProt table (all organisms): '
                "'fallback' (default — curated first, full-UniProt only for "
                "unresolved IDs), 'never', 'both', 'only'."
            ),
        ),
        name_scope: str = Parameter(
            default='default',
            required=False,
            description=(
                "Only consulted when id_type or target_id_type is 'name': "
                "'default' (preferred names, then synonyms, when "
                "translating from a name; preferred only when translating "
                "to one), 'preferred', 'synonym' or 'any'."
            ),
        ),
    ) -> dict:
        """Translate identifiers from one type to another.

        Args:
            identifiers: Comma-separated identifiers.
            id_type: Source ID type.
            target_id_type: Target ID type.
            ncbi_tax_id: NCBI Taxonomy ID.
            raw: Skip special-case handling.
            backend: Force specific backend (ignored for DB mode).
            full_uniprot: How to use the comprehensive full-UniProt table.
            name_scope: Restrict a 'name'-axis query to a name bucket.
        """

        from omnipath_utils.mapping._id_types import IdTypeRegistry

        reg = IdTypeRegistry.get()
        id_type_resolved = reg.resolve(id_type) or id_type
        target_resolved = reg.resolve(target_id_type) or target_id_type
        is_name_axis = 'name' in (id_type_resolved, target_resolved)
        is_structure_query = _needs_chemistry_toolkit(
            id_type_resolved, target_resolved,
        )

        id_list = [i.strip() for i in identifiers.split(',') if i.strip()]

        # US2: deprecated-ID recovery is a default API fallback (recover=True),
        # consulted only for ids the primary route misses; raw mode stays literal.
        recovery: dict = {}
        matched_as: dict = {}
        result, backends_used = translate_ids(
            session,
            id_list,
            id_type_resolved,
            target_resolved,
            ncbi_tax_id,
            full_uniprot=full_uniprot,
            recover=not raw,
            recovery_meta=recovery,
            name_scope=name_scope,
            name_meta=matched_as,
        )

        # A name-axis miss is never queued for on-demand loading -- the
        # public service performs none, and reports an honest miss instead
        # of a retry note (FR: "Honest misses").
        loading = False if is_name_axis else _maybe_trigger_load(
            result,
            raw,
            id_type_resolved,
            target_resolved,
            ncbi_tax_id,
        )

        if not raw:
            result = _apply_fallbacks(
                session,
                id_list,
                id_type_resolved,
                target_resolved,
                ncbi_tax_id,
                result,
            )
            result = _resolve_and_cleanup(
                result,
                id_type_resolved,
                target_resolved,
                ncbi_tax_id,
                session=session,
            )

        return _build_translate_response(
            id_list,
            result,
            id_type,
            target_id_type,
            ncbi_tax_id,
            raw,
            backends_used,
            backend,
            loading,
            recovery=recovery,
            name_scope=name_scope if is_name_axis else None,
            matched_as=matched_as if is_name_axis else None,
            structure_capability=(
                _structure_key_capability_meta() if is_structure_query else None
            ),
        )

    @post('/translate')
    async def translate_post(
        self,
        session: Session,
        data: dict,
    ) -> dict:
        """Translate identifiers via POST with JSON body.

        Body fields:
            identifiers: List of identifiers.
            id_type: Source ID type.
            target_id_type: Target ID type.
            ncbi_tax_id: NCBI Taxonomy ID (default: 9606).
            raw: Skip special-case handling (default: false).
            backend: Force specific backend (default: null).
            full_uniprot: Use of the comprehensive full-UniProt table
                ('fallback' (default), 'never', 'both', 'only').
            name_scope: Restrict a 'name'-axis query to a name bucket
                ('default' (default), 'preferred', 'synonym', 'any').
        """

        from omnipath_utils.mapping._id_types import IdTypeRegistry

        reg = IdTypeRegistry.get()

        id_list = data.get('identifiers', [])
        id_type = data.get('id_type', '')
        target_id_type = data.get('target_id_type', '')
        ncbi_tax_id = data.get('ncbi_tax_id', 9606)
        raw = data.get('raw', False)
        backend = data.get('backend', None)
        full_uniprot = data.get('full_uniprot', 'fallback')
        name_scope = data.get('name_scope', 'default')

        id_type_resolved = reg.resolve(id_type) or id_type
        target_resolved = reg.resolve(target_id_type) or target_id_type
        is_name_axis = 'name' in (id_type_resolved, target_resolved)
        is_structure_query = _needs_chemistry_toolkit(
            id_type_resolved, target_resolved,
        )

        recovery: dict = {}
        matched_as: dict = {}
        result, backends_used = translate_ids(
            session,
            id_list,
            id_type_resolved,
            target_resolved,
            ncbi_tax_id,
            full_uniprot=full_uniprot,
            recover=not raw,
            recovery_meta=recovery,
            name_scope=name_scope,
            name_meta=matched_as,
        )

        loading = False if is_name_axis else _maybe_trigger_load(
            result,
            raw,
            id_type_resolved,
            target_resolved,
            ncbi_tax_id,
        )

        if not raw:
            result = _apply_fallbacks(
                session,
                id_list,
                id_type_resolved,
                target_resolved,
                ncbi_tax_id,
                result,
            )
            result = _resolve_and_cleanup(
                result,
                id_type_resolved,
                target_resolved,
                ncbi_tax_id,
                session=session,
            )

        return _build_translate_response(
            id_list,
            result,
            id_type,
            target_id_type,
            ncbi_tax_id,
            raw,
            backends_used,
            backend,
            loading,
            recovery=recovery,
            name_scope=name_scope if is_name_axis else None,
            matched_as=matched_as if is_name_axis else None,
            structure_capability=(
                _structure_key_capability_meta() if is_structure_query else None
            ),
        )

    @get('/identify')
    async def identify(
        self,
        session: Session,
        identifiers: str = Parameter(
            description='Comma-separated identifiers',
        ),
        ncbi_tax_id: int = Parameter(
            default=9606,
            description='NCBI Taxonomy ID',
        ),
    ) -> dict:
        """Try to identify the type of given identifiers.

        Searches all mapping tables in the database to find which ID types
        contain the given identifiers.
        """

        id_list = [i.strip() for i in identifiers.split(',') if i.strip()]
        result = identify_ids(session, id_list, ncbi_tax_id)

        return {
            'results': result,
            'meta': {
                'ncbi_tax_id': ncbi_tax_id,
                'total_input': len(id_list),
            },
        }

    @get('/all')
    async def all_mappings(
        self,
        session: Session,
        identifiers: str = Parameter(
            description='Comma-separated identifiers',
        ),
        id_type: str = Parameter(description='Source ID type'),
        ncbi_tax_id: int = Parameter(
            default=9606,
            description='NCBI Taxonomy ID',
        ),
    ) -> dict:
        """Return all known mappings for identifiers across all target types."""

        from omnipath_utils.mapping._id_types import IdTypeRegistry

        reg = IdTypeRegistry.get()
        id_type_resolved = reg.resolve(id_type) or id_type

        id_list = [i.strip() for i in identifiers.split(',') if i.strip()]
        result = get_all_mappings(
            session, id_list, id_type_resolved, ncbi_tax_id
        )

        return {
            'results': result,
            'meta': {
                'id_type': id_type,
                'ncbi_tax_id': ncbi_tax_id,
                'total_input': len(id_list),
            },
        }

    @get('/organism')
    async def organism(
        self,
        session: Session,
        uniprot: str = Parameter(
            description='Comma-separated UniProt accessions',
        ),
    ) -> dict:
        """Return the source organism (NCBI taxon id) of UniProt accessions.

        Backed by the resolver's accession -> organism map, so it answers even
        for accessions with no gene/protein resolution (the TrEMBL long tail).
        Accessions may carry a CURIE prefix or an isoform suffix. The response
        is keyed by the accession as sent; an unknown accession maps to null.
        """

        accessions = [a.strip() for a in uniprot.split(',') if a.strip()]
        result = uniprot_organism(session, accessions)

        return {
            'results': result,
            'meta': {'total_input': len(accessions)},
        }

    @get('/id-types')
    async def id_types(self) -> list[dict]:
        """List all supported ID types."""

        from omnipath_utils.mapping._id_types import IdTypeRegistry

        reg = IdTypeRegistry.get()

        return [
            {
                'name': name,
                'label': (reg.info(name) or {}).get('label'),
                'entity_type': (reg.info(name) or {}).get('entity_type'),
                'curie_prefix': (reg.info(name) or {}).get('curie_prefix'),
                'url_pattern': reg.url_pattern(name),
                'id_pattern': reg.id_pattern(name),
            }
            for name in reg.all_names()
        ]

    @get('/loading')
    async def loading_status(self) -> dict:
        """Check which tables are currently being loaded."""

        from omnipath_utils.db._loader import _pending

        return {
            'loading': [
                {
                    'source_type': s,
                    'target_type': t,
                    'ncbi_tax_id': n,
                }
                for s, t, n in sorted(_pending)
            ],
            'count': len(_pending),
        }

    @get('/table')
    async def get_table(
        self,
        session: Session,
        id_type: str = Parameter(description='Source ID type'),
        target_id_type: str = Parameter(description='Target ID type'),
        ncbi_tax_id: int = Parameter(default=9606),
    ) -> dict:
        """Get a full translation table."""

        from omnipath_utils.mapping._id_types import IdTypeRegistry

        reg = IdTypeRegistry.get()
        id_type_r = reg.resolve(id_type) or id_type
        target_r = reg.resolve(target_id_type) or target_id_type

        data, _ = translate_ids(session, None, id_type_r, target_r, ncbi_tax_id)

        # Convert sets to sorted lists for JSON
        result = {k: sorted(v) for k, v in data.items()}

        return {
            'table': result,
            'meta': {
                'id_type': id_type,
                'target_id_type': target_id_type,
                'ncbi_tax_id': ncbi_tax_id,
                'count': len(result),
            },
        }

