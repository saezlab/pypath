"""Canonical-target resolver export (Milestone I).

Projects the symmetric ``id_mapping`` model into canonical-target rows
``(source_type, source_id, canonical_target, taxonomy_id)`` for an entity family,
honouring a shared policy file. Builds on the DB translation layer (not a raw
pair dump). ``taxonomy_id`` is text (cast from the Integer ``ncbi_tax_id``;
chemicals use the organism-agnostic ``0`` convention).

- protein → per-taxon canonical **UniProt** (reuses ``uniprot_cleanup_batch``).
- chemical → organism-agnostic **Standard InChI Key**.
"""

from __future__ import annotations

import datetime
import json
import logging
import os
from dataclasses import dataclass, field
from typing import Any

import yaml

_log = logging.getLogger(__name__)

MANIFEST_FILENAME = 'manifest.json'

# Canonical target id_type per family (overridable via --canonical-type).
_DEFAULT_CANONICAL = {'protein': 'uniprot', 'chemical': 'inchikey'}
# Canonical-type display names (CLI / contract) → internal id_types.yaml keys.
_CANONICAL_ALIASES = {
    'uniprot': 'uniprot',
    'standard inchi key': 'inchikey',
    'inchikey': 'inchikey',
}
_CHEMICAL_TAXON = 0  # organism-agnostic convention for chemicals


@dataclass(frozen=True)
class PolicyRule:
    key_type: str
    action: str = 'accept'  # accept | candidate_only | ignore
    requires_taxonomy: bool = False


@dataclass
class ExportStats:
    family: str
    canonical_type: str
    files: list[str] = field(default_factory=list)
    rows: int = 0


def load_policy(path: str) -> dict[str, list[PolicyRule]]:
    """Read the shared policy YAML: ``entity_family → [{key_type, action, …}]``."""
    with open(path, encoding='utf-8') as handle:
        raw = yaml.safe_load(handle) or {}
    policy: dict[str, list[PolicyRule]] = {}
    for family, rules in raw.items():
        policy[family] = [
            PolicyRule(
                key_type=rule['key_type'],
                action=rule.get('action', 'accept'),
                requires_taxonomy=bool(rule.get('requires_taxonomy', False)),
            )
            for rule in (rules or [])
        ]
    return policy


def _canonical_type(family: str, requested: str | None) -> str:
    if requested:
        key = _CANONICAL_ALIASES.get(requested.strip().lower())
        if not key:
            raise ValueError(f'Unknown canonical type: {requested!r}')
        return key
    return _DEFAULT_CANONICAL[family]


def _accepted_source_types(
    policy: dict[str, list[PolicyRule]],
    family: str,
    canonical_type: str,
) -> list[str]:
    """Source key types to project (policy ``accept``), excluding the target itself."""
    rules = policy.get(family)
    if not rules:
        raise ValueError(
            f'No policy rules for entity family {family!r}; nothing to export.'
        )
    return [
        rule.key_type
        for rule in rules
        if rule.action == 'accept' and rule.key_type != canonical_type
    ]


def _taxa_for_family(family: str, taxa: list[int] | None) -> list[int]:
    if family == 'chemical':
        return [_CHEMICAL_TAXON]
    return taxa or [9606]


def export_resolver(
    *,
    family: str,
    policy_path: str,
    output_dir: str,
    canonical_type: str | None = None,
    db_url: str | None = None,
    taxa: list[int] | None = None,
    max_records: int | None = None,
    export_id: str | None = None,
) -> ExportStats:
    """Project ``id_mapping`` to canonical-target parquet for ``family``.

    Writes the parquet tables plus a ``manifest.json`` recording the export
    identity (``export_id``), the omnipath-utils version, the per-family row
    counts, and a fingerprint of the backing DB. omnipath-build records the
    consumed ``export_id`` in its own build manifest (Principle V).
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    from omnipath_utils.db._connection import get_engine

    canonical = _canonical_type(family, canonical_type)
    policy = load_policy(policy_path)
    source_types = _accepted_source_types(policy, family, canonical)
    os.makedirs(output_dir, exist_ok=True)
    engine = get_engine(db_url)
    stats = ExportStats(family=family, canonical_type=canonical)
    per_taxon: dict[str, int] = {}

    for taxon in _taxa_for_family(family, taxa):
        source_col, source_id_col, target_col, taxon_col = _project_rows(
            engine, source_types, canonical, taxon, family, max_records
        )
        if not source_id_col:
            _log.warning('no %s rows for taxon %s', family, taxon)
            continue
        table = pa.table(
            {
                'source_type': source_col,
                'source_id': source_id_col,
                'canonical_target': target_col,
                'taxonomy_id': taxon_col,
            }
        )
        if family == 'chemical':
            path = os.path.join(output_dir, 'chemicals.parquet')
        else:
            os.makedirs(os.path.join(output_dir, 'proteins'), exist_ok=True)
            path = os.path.join(output_dir, 'proteins', f'{taxon}.parquet')
        pq.write_table(table, path, compression='snappy')
        stats.files.append(path)
        stats.rows += table.num_rows
        per_taxon[str(taxon)] = table.num_rows
        _log.info('wrote %s (%d rows)', path, table.num_rows)

    if not stats.files:
        raise ValueError(
            f'No data to export for family {family!r} → {canonical!r}; the '
            f'canonical target has no backing mappings in this build.'
        )

    _write_manifest(
        output_dir=output_dir,
        family=family,
        canonical=canonical,
        source_types=source_types,
        per_taxon=per_taxon,
        rows=stats.rows,
        files=stats.files,
        engine=engine,
        export_id=export_id,
    )
    return stats


def _utils_version() -> str:
    try:
        from importlib.metadata import version

        return version('omnipath-utils')
    except Exception:  # pragma: no cover - best-effort provenance
        return 'unknown'


def _db_fingerprint(engine) -> dict[str, Any]:
    """Cheap provenance: curated row count + full-UniProt reltuples estimate."""
    from sqlalchemy import text

    from omnipath_utils.db._connection import SCHEMA

    fp: dict[str, Any] = {}
    try:
        with engine.connect() as conn:
            fp['curated_rows'] = conn.execute(
                text(f'SELECT count(*) FROM {SCHEMA}.id_mapping')
            ).scalar()
            fp['full_uniprot_rows'] = conn.execute(
                text(
                    'SELECT reltuples::bigint FROM pg_class '
                    f"WHERE oid = to_regclass('{SCHEMA}.id_mapping_ftp')"
                )
            ).scalar()
    except Exception as err:  # pragma: no cover - best-effort provenance
        fp['error'] = str(err)
    return fp


def _write_manifest(
    *,
    output_dir: str,
    family: str,
    canonical: str,
    source_types: list[str],
    per_taxon: dict[str, int],
    rows: int,
    files: list[str],
    engine,
    export_id: str | None,
) -> str:
    """Write/merge ``manifest.json`` (one entry per exported family)."""
    path = os.path.join(output_dir, MANIFEST_FILENAME)
    manifest: dict[str, Any] = {}
    if os.path.exists(path):
        try:
            with open(path, encoding='utf-8') as handle:
                manifest = json.load(handle) or {}
        except Exception:  # pragma: no cover - regenerate on corruption
            manifest = {}

    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    # ``export_id`` is stable per output dir; set on first write.
    manifest.setdefault(
        'export_id', export_id or now.replace(':', '').replace('-', '')
    )
    manifest['omnipath_utils_version'] = _utils_version()
    manifest['db_fingerprint'] = _db_fingerprint(engine)
    families = manifest.setdefault('families', {})
    families[family] = {
        'canonical_type': canonical,
        'source_types': source_types,
        'rows': rows,
        'per_taxon': per_taxon,
        'files': [os.path.relpath(f, output_dir) for f in files],
        'created_at': now,
    }
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
    _log.info('wrote %s (export_id=%s)', path, manifest['export_id'])
    return manifest['export_id']


def _project_rows(
    engine,
    source_types: list[str],
    canonical: str,
    taxon: int,
    family: str,
    max_records: int | None,
):
    """One clean pass: emit aligned (source_type, source_id, target, taxon) columns."""
    from collections import defaultdict

    from sqlalchemy.orm import Session

    from omnipath_utils.db._query import get_full_table, translate_ids

    s_type: list[str] = []
    s_id: list[str] = []
    target: list[str] = []
    taxon_text: list[str] = []
    with Session(engine) as session:
        for source_type in source_types:
            if family == 'protein':
                # The comprehensive full-UniProt idmapping is stored in the
                # native ``uniprot -> X`` direction. Query that direction over
                # both the curated and full-UniProt tables, then invert to the
                # ``X -> uniprot`` projection the resolver consumes. (A direct
                # ``X -> uniprot`` whole-table query would only ever see the
                # curated table — the full table is never a blanket fallback.)
                native, _ = translate_ids(
                    session, None, canonical, source_type, taxon,
                    full_uniprot='both',
                )
                inverted: defaultdict[str, set[str]] = defaultdict(set)
                for uniprot_ac, xs in native.items():
                    for x in xs:
                        inverted[x].add(uniprot_ac)
                from omnipath_utils.mapping._cleanup import uniprot_cleanup_batch

                table = uniprot_cleanup_batch(
                    dict(inverted), taxon, session=session
                )
            else:
                table = get_full_table(session, source_type, canonical, taxon)
            if not table:
                continue
            for source_id, targets in table.items():
                for tgt in targets:
                    s_type.append(source_type)
                    s_id.append(str(source_id))
                    target.append(str(tgt))
                    taxon_text.append(str(taxon))
                    if max_records and len(s_id) >= max_records:
                        return s_type, s_id, target, taxon_text
    return s_type, s_id, target, taxon_text
