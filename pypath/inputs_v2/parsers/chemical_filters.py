"""Shared chemical resolver filters for large ligand datasets."""

from __future__ import annotations

from functools import lru_cache
import os
from pathlib import Path
from typing import Any

try:
    import polars as pl
except ImportError:  # pragma: no cover - optional dependency
    pl = None


DEFAULT_CHEMICAL_LOOKUP_PATHS = (
    Path('minimal/data/chemicals/chemical_identifier_lookup.parquet'),
    Path('data/chemicals/chemical_identifier_lookup.parquet'),
)
DEFAULT_ALLOWED_SOURCES = frozenset({'chebi', 'hmdb'})


def chemical_resolver_inchikey_filter_enabled(kwargs: dict[str, Any]) -> bool:
    return bool(
        kwargs.get(
            'filter_chemical_resolver_inchikeys',
            os.environ.get(
                'OMNIPATH_FILTER_CHEMICAL_RESOLVER_INCHIKEYS',
                '1',
            ).lower()
            not in {'0', 'false', 'no'},
        )
    )


def row_has_allowed_chemical_inchikey(
    row: dict[str, Any],
    key: str,
    *,
    kwargs: dict[str, Any],
) -> bool:
    if not chemical_resolver_inchikey_filter_enabled(kwargs):
        return True
    inchikey = clean_inchikey(row.get(key))
    allowed = allowed_chemical_inchikeys(kwargs)
    if not allowed:
        return True
    return bool(inchikey and inchikey in allowed)


def clean_inchikey(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {'none', 'inchikey=none'}:
        return None
    text = text.removeprefix('InChIKey=').strip()
    return text or None


def allowed_chemical_inchikeys(kwargs: dict[str, Any]) -> frozenset[str]:
    sources = kwargs.get('chemical_resolver_sources')
    if sources is None:
        source_values = DEFAULT_ALLOWED_SOURCES
    elif isinstance(sources, str):
        source_values = frozenset(
            item.strip().lower()
            for item in sources.replace(',', ' ').split()
            if item.strip()
        )
    else:
        source_values = frozenset(str(item).strip().lower() for item in sources)
    return _allowed_chemical_inchikeys(_lookup_path(kwargs), source_values)


def chemical_resolver_lookup_path(kwargs: dict[str, Any]) -> Path:
    return Path(_lookup_path(kwargs))


def chemical_resolver_sources(kwargs: dict[str, Any]) -> frozenset[str]:
    sources = kwargs.get('chemical_resolver_sources')
    if sources is None:
        return DEFAULT_ALLOWED_SOURCES
    if isinstance(sources, str):
        return frozenset(
            item.strip().lower()
            for item in sources.replace(',', ' ').split()
            if item.strip()
        )
    return frozenset(str(item).strip().lower() for item in sources)


def chemical_resolver_filter_sql(
    *,
    lookup_path: Path,
    inchikey_expr: str,
    sources: frozenset[str] = DEFAULT_ALLOWED_SOURCES,
) -> str:
    escaped_path = str(lookup_path).replace("'", "''")
    source_list = ', '.join(
        "'" + source.replace("'", "''") + "'"
        for source in sorted(sources)
    )
    if pl is not None and lookup_path.exists():
        columns = set(pl.scan_parquet(lookup_path).collect_schema().names())
    else:
        columns = {'source', 'standard_inchi_key'}
    if {'key_identifier_type_id', 'canonical_identifier'} <= columns:
        identifier_type_path = lookup_path.with_name('identifier_type.parquet')
        escaped_type_path = str(identifier_type_path).replace("'", "''")
        allowed_select = f"""
          SELECT DISTINCT l.canonical_identifier
          FROM read_parquet('{escaped_path}') l
          JOIN read_parquet('{escaped_type_path}') t
            ON t.identifier_type_id = l.key_identifier_type_id
          WHERE lower(split_part(t.name, ':', 1)) IN ({source_list})
            AND l.canonical_identifier IS NOT NULL
            AND l.canonical_identifier <> ''
        """
    else:
        allowed_select = f"""
          SELECT DISTINCT standard_inchi_key
          FROM read_parquet('{escaped_path}')
          WHERE lower(source) IN ({source_list})
            AND standard_inchi_key IS NOT NULL
            AND standard_inchi_key <> ''
        """
    return f"""
        {inchikey_expr} IS NOT NULL
        AND {inchikey_expr} <> ''
        AND {inchikey_expr} IN ({allowed_select})
    """


def _lookup_path(kwargs: dict[str, Any]) -> str:
    explicit = (
        kwargs.get('chemical_resolver_lookup_path')
        or os.environ.get('OMNIPATH_CHEMICAL_RESOLVER_LOOKUP')
    )
    if explicit:
        return str(explicit)
    for path in DEFAULT_CHEMICAL_LOOKUP_PATHS:
        if path.exists():
            return str(path)
    return str(DEFAULT_CHEMICAL_LOOKUP_PATHS[0])


@lru_cache(maxsize=16)
def _allowed_chemical_inchikeys(
    lookup_path: str,
    sources: frozenset[str],
) -> frozenset[str]:
    path = Path(lookup_path)
    if not path.exists():
        print(
            f'Chemical resolver lookup not found; ligand filter disabled: {path}',
            flush=True,
        )
        return frozenset()
    if pl is None:
        raise ImportError('polars is required for chemical resolver filtering.')

    scan = pl.scan_parquet(path)
    columns = set(scan.collect_schema().names())
    if {'source', 'standard_inchi_key'} <= columns:
        frame = scan.filter(
            pl.col('source').str.to_lowercase().is_in(sorted(sources))
            & pl.col('standard_inchi_key').is_not_null()
            & (pl.col('standard_inchi_key') != '')
        ).select('standard_inchi_key')
    elif {'key_identifier_type_id', 'canonical_identifier'} <= columns:
        identifier_type_path = path.with_name('identifier_type.parquet')
        if not identifier_type_path.exists():
            raise FileNotFoundError(
                f'Missing identifier type parquet for normalized chemical '
                f'lookup: {identifier_type_path}'
            )
        type_ids = (
            pl.scan_parquet(identifier_type_path)
            .filter(
                pl.col('name')
                .str.to_lowercase()
                .str.split(':')
                .list.first()
                .is_in(sorted(sources))
            )
            .select('identifier_type_id')
            .collect()
            .get_column('identifier_type_id')
            .to_list()
        )
        frame = scan.filter(
            pl.col('key_identifier_type_id').is_in(type_ids)
            & pl.col('canonical_identifier').is_not_null()
            & (pl.col('canonical_identifier') != '')
        ).select(pl.col('canonical_identifier').alias('standard_inchi_key'))
    else:
        raise ValueError(
            f'Unsupported chemical resolver lookup schema in {path}: '
            f'{sorted(columns)}'
        )

    values = frame.unique().collect().get_column('standard_inchi_key').to_list()
    return frozenset(clean_inchikey(value) for value in values if clean_inchikey(value))
