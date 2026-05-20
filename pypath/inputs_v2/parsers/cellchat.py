"""
Parser for CellChatDB .rda files.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any
import warnings

import rdata


def _first_handle(opener) -> Any | None:
    """Extract the first file handle from an opener result."""
    if not opener or not opener.result:
        return None
    if isinstance(opener.result, dict):
        return next(iter(opener.result.values()), None)
    if hasattr(opener.result, '__next__'):
        return opener.result
    return opener.result


def _read_bytes(opener) -> bytes:
    handle = _first_handle(opener)
    if handle is None:
        return b''
    if hasattr(handle, 'read'):
        data = handle.read()
        return data if isinstance(data, bytes) else data.encode('utf-8')
    return b''.join(
        chunk if isinstance(chunk, bytes) else str(chunk).encode('utf-8')
        for chunk in handle
    )


def _load_db(opener) -> dict[str, Any]:
    data = _read_bytes(opener)
    if not data:
        return {}
    parsed = rdata.parser.parse_data(data)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message='Missing constructor for R class "tbl.*',
            category=UserWarning,
        )
        converted = rdata.conversion.convert(parsed)
    if not converted:
        return {}
    db = next(iter(converted.values()))
    return db if isinstance(db, dict) else {}


def _clean(value: Any) -> str:
    if value is None:
        return ''
    if hasattr(value, 'item'):
        value = value.item()
    text = str(value).strip()
    return '' if text.lower() in {'nan', 'none', 'na'} else text


def _split_symbols(value: Any) -> list[str]:
    text = _clean(value)
    if not text:
        return []
    return [
        item.strip()
        for item in text.replace(';', ',').split(',')
        if item.strip()
    ]


def _id_token(value: Any) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', _clean(value)).strip('_')


def _row_dict(row: Any) -> dict[str, Any]:
    return {
        str(key).replace('.', '_'): _clean(value)
        for key, value in row.items()
    }


def _complex_members(complexes: Any, name: str) -> list[str]:
    if complexes is None or not name:
        return []
    try:
        if name not in complexes.index:
            return []
        row = complexes.loc[name]
    except Exception:
        return []
    return [
        gene
        for key, value in row.items()
        if str(key).startswith('subunit_') and (gene := _clean(value))
    ]


def _cofactor_members(cofactors: Any, name: str) -> list[str]:
    if cofactors is None or not name:
        return []
    try:
        if name not in cofactors.index:
            return []
        row = cofactors.loc[name]
    except Exception:
        return []
    return [
        gene
        for key, value in row.items()
        if str(key).startswith('cofactor') and (gene := _clean(value))
    ]


def _participant_genes(db: dict[str, Any], row: dict[str, Any], role: str) -> list[str]:
    name = row.get(role) or ''
    complex_genes = _complex_members(db.get('complex'), name)
    if complex_genes:
        return complex_genes
    return _split_symbols(row.get(f'{role}_symbol'))


def _iter_interaction_rows(
    db: dict[str, Any],
    taxon_id: str,
) -> Generator[dict[str, Any], None, None]:
    interactions = db.get('interaction')
    if interactions is None:
        return

    for _, raw_row in interactions.iterrows():
        row = _row_dict(raw_row)
        row.update(
            {
                'taxon_id': str(taxon_id),
                'ligand_name': row.get('ligand', ''),
                'ligand_genes': _participant_genes(db, row, 'ligand'),
                'receptor_name': row.get('receptor', ''),
                'receptor_genes': _participant_genes(db, row, 'receptor'),
            }
        )
        yield row


def iter_cellchat_interactions(
    opener,
    *,
    taxon_id: str,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield normalized CellChatDB ligand-receptor interaction records."""
    db = _load_db(opener)
    yield from _iter_interaction_rows(db, taxon_id)


def iter_cellchat_cofactor_interactions(
    opener,
    *,
    taxon_id: str,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield normalized CellChatDB cofactor-target interaction records."""
    db = _load_db(opener)
    cofactors = db.get('cofactor')

    for row in _iter_interaction_rows(db, taxon_id):
        for column, role, effect, target_role in (
            ('agonist', 'agonist', 'stimulation', 'ligand'),
            ('antagonist', 'antagonist', 'inhibition', 'ligand'),
            ('co_A_receptor', 'co-activating receptor', 'stimulation', 'receptor'),
            ('co_I_receptor', 'co-inhibitory receptor', 'inhibition', 'receptor'),
        ):
            cofactor_name = row.get(column, '')
            cofactor_genes = _cofactor_members(cofactors, cofactor_name)
            if not cofactor_genes:
                continue

            target_name = row.get(f'{target_role}_name', '')
            for cofactor_gene in cofactor_genes:
                yield {
                    **row,
                    'cofactor_role': role,
                    'effect': effect,
                    'target_role': target_role,
                    'cofactor_group': cofactor_name,
                    'cofactor_gene': cofactor_gene,
                    'target_name': target_name,
                    'target_genes': row.get(f'{target_role}_genes', []),
                    'target_family': row.get(f'{target_role}_family', ''),
                    'target_location': row.get(f'{target_role}_location', ''),
                    'target_keyword': row.get(f'{target_role}_keyword', ''),
                    'target_secreted_type': row.get(f'{target_role}_secreted_type', ''),
                    'target_transmembrane': row.get(f'{target_role}_transmembrane', ''),
                    'cofactor_interaction_name': (
                        f'{_id_token(cofactor_gene)}_{effect}_{_id_token(target_name)}_'
                        f'{_id_token(row.get("interaction_name", ""))}'
                    ),
                }


def iter_cellchat_complexes(
    opener,
    *,
    taxon_id: str,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield normalized CellChatDB protein complex records."""
    db = _load_db(opener)
    complexes = db.get('complex')
    if complexes is None:
        return

    for name, row in complexes.iterrows():
        genes = [
            gene
            for key, value in row.items()
            if str(key).startswith('subunit_') and (gene := _clean(value))
        ]
        if genes:
            yield {
                'name': _clean(name),
                'genes': genes,
                'taxon_id': str(taxon_id),
            }


def iter_cellchat_cofactors(
    opener,
    *,
    taxon_id: str,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield normalized CellChatDB cofactor group records."""
    db = _load_db(opener)
    cofactors = db.get('cofactor')
    if cofactors is None:
        return

    for name, _row in cofactors.iterrows():
        genes = _cofactor_members(cofactors, _clean(name))
        if genes:
            yield {
                'name': _clean(name),
                'genes': genes,
                'taxon_id': str(taxon_id),
            }
