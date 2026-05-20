"""BindingDB TSV parser (2026-03+ schema)."""

from __future__ import annotations

from collections.abc import Generator
import csv
import math
import os
from pathlib import Path
import re
import zipfile

try:
    import duckdb
except ImportError:  # pragma: no cover - optional dependency
    duckdb = None

from pypath.inputs_v2.parsers.chemical_filters import (
    chemical_resolver_filter_sql,
    chemical_resolver_lookup_path,
    chemical_resolver_sources,
    row_has_allowed_chemical_inchikey,
)
from pypath.inputs_v2.parsers.base import iter_parquet


_CHEMBL_RUN_RE = re.compile(r'(CHEMBL\d+)(?=CHEMBL\d+)')
_AFFINITY_VALUE_RE = re.compile(r'([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)')
_BINDINGDB_AFFINITY_COLUMNS = (
    'Ki (nM)',
    'Kd (nM)',
    'IC50 (nM)',
    'EC50 (nM)',
)
BINDINGDB_MIN_PCHEMBL = 5.0

# Keep the raw parser schema narrow. These are the only BindingDB columns used by
# pypath.inputs_v2.bindingdb.interactions_schema.
_BINDINGDB_COLUMNS = [
    'BindingDB Reactant_set_id',
    'Ki (nM)',
    'Kd (nM)',
    'IC50 (nM)',
    'EC50 (nM)',
    'kon (M-1-s-1)',
    'koff (s-1)',
    'pH',
    'Temp (C)',
    'PMID',
    'Article DOI',
    'Patent Number',
    'Curation/DataSource',
    'BindingDB MonomerID',
    'BindingDB Ligand Name',
    'Ligand InChI Key',
    'Ligand InChI',
    'Ligand SMILES',
    'PubChem CID',
    'PubChem SID',
    'ChEBI ID of Ligand',
    'ChEMBL ID of Ligand',
    'DrugBank ID of Ligand',
    'KEGG ID of Ligand',
    'ZINC ID of Ligand',
    'Target Name',
    'Target Source Organism According to Curator or DataSource',
    'UniProt (SwissProt) Primary ID of Target Chain 1',
    'UniProt (SwissProt) Recommended Name of Target Chain 1',
    'UniProt (TrEMBL) Primary ID of Target Chain 1',
    'UniProt (TrEMBL) Submitted Name of Target Chain 1',
]
BINDINGDB_FILTERED_PARQUET_CACHE_VERSION = 2


def _normalize_row(row: dict[str, str | None]) -> dict[str, str]:
    normalized = {key: '' if value is None else str(value) for key, value in row.items()}
    value = normalized.get('ChEMBL ID of Ligand', '').strip()
    if value and value.count('CHEMBL') > 1 and '::' not in value and ';' not in value and '|' not in value:
        normalized['ChEMBL ID of Ligand'] = _CHEMBL_RUN_RE.sub(r'\1::', value)
    pchembl_value = _bindingdb_pchembl_value(normalized)
    normalized['pchembl_value'] = (
        f'{pchembl_value:.6g}'
        if pchembl_value is not None
        else ''
    )
    return normalized


def _pchembl_from_nm(value: object) -> float | None:
    """Convert a BindingDB nM affinity cell to a pChEMBL-like value.

    Values with a leading ``>`` are weak-activity upper bounds in pChEMBL space,
    so they are not used to pass the potency filter.
    """
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.startswith('>'):
        return None
    match = _AFFINITY_VALUE_RE.search(text)
    if not match:
        return None
    try:
        nm_value = float(match.group(1))
    except ValueError:
        return None
    if nm_value <= 0:
        return None
    return 9.0 - math.log10(nm_value)


def _bindingdb_pchembl_value(row: dict[str, object]) -> float | None:
    values = [
        value
        for value in (
            _pchembl_from_nm(row.get(column))
            for column in _BINDINGDB_AFFINITY_COLUMNS
        )
        if value is not None
    ]
    return max(values) if values else None


def _bindingdb_pchembl_filter_enabled(kwargs: dict[str, object]) -> bool:
    return bool(kwargs.get('filter_bindingdb_pchembl', True))


def _bindingdb_chemical_resolver_filter_enabled(kwargs: dict[str, object]) -> bool:
    return bool(kwargs.get('filter_chemical_resolver_inchikeys', False))


def _bindingdb_min_pchembl(kwargs: dict[str, object]) -> float:
    return float(kwargs.get('bindingdb_min_pchembl', BINDINGDB_MIN_PCHEMBL))


def _row_has_allowed_pchembl(
    row: dict[str, str],
    *,
    kwargs: dict[str, object],
) -> bool:
    if not _bindingdb_pchembl_filter_enabled(kwargs):
        return True
    value = row.get('pchembl_value')
    if value in (None, ''):
        return False
    try:
        return float(value) > _bindingdb_min_pchembl(kwargs)
    except ValueError:
        return False


def _row_passes_filters(
    row: dict[str, str],
    *,
    kwargs: dict[str, object],
) -> bool:
    if not _row_has_allowed_pchembl(row, kwargs=kwargs):
        return False
    if _bindingdb_chemical_resolver_filter_enabled(kwargs):
        return row_has_allowed_chemical_inchikey(
            row,
            'Ligand InChI Key',
            kwargs=kwargs,
        )
    return True


def _quote_identifier(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _quote_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _read_header(tsv_path: Path) -> set[str]:
    with tsv_path.open('r', encoding='utf-8', errors='replace', newline='') as handle:
        header = handle.readline().rstrip('\r\n')
    return set(header.split('\t')) if header else set()


def _bindingdb_tsv_path(opener, *, extract: bool = True) -> Path | None:
    """Return an on-disk TSV path, extracting the zip member if necessary.

    DuckDB's CSV reader operates on paths. The downloaded BindingDB archive is a
    zip containing a single very large TSV, so we extract it in streaming chunks
    next to the zip. This costs disk space but keeps memory bounded and is
    reused by later preparses.
    """
    archive_path_raw = getattr(opener, 'path', None)
    if not archive_path_raw:
        return None

    archive_path = Path(archive_path_raw)
    if archive_path.suffix.lower() != '.zip' or not archive_path.exists():
        return archive_path if archive_path.exists() else None

    tsv_path = archive_path.with_suffix('.tsv')
    if tsv_path.exists() and tsv_path.stat().st_mtime_ns >= archive_path.stat().st_mtime_ns:
        return tsv_path
    if not extract:
        return None

    tmp_path = tsv_path.with_suffix('.tsv.tmp')
    if tmp_path.exists():
        tmp_path.unlink()

    with zipfile.ZipFile(archive_path) as archive:
        member_name = next((name for name in archive.namelist() if name.lower().endswith('.tsv')), None)
        if member_name is None:
            return None
        print(f'Extracting BindingDB TSV for DuckDB: {tsv_path}', flush=True)
        with archive.open(member_name) as src, tmp_path.open('wb') as dst:
            while chunk := src.read(1024 * 1024):
                dst.write(chunk)

    tmp_path.replace(tsv_path)
    return tsv_path


def _iter_duckdb_tsv(
    tsv_path: Path,
    max_lines: int | None = None,
    batch_size: int = 50_000,
) -> Generator[dict[str, str], None, None]:
    if duckdb is None:
        raise ImportError('duckdb is required for DuckDB-based BindingDB parsing.')

    available_columns = _read_header(tsv_path)
    select_exprs = [
        (
            f'{_quote_identifier(column)} AS {_quote_identifier(column)}'
            if column in available_columns
            else f'NULL AS {_quote_identifier(column)}'
        )
        for column in _BINDINGDB_COLUMNS
    ]
    limit_clause = f' LIMIT {int(max_lines)}' if max_lines is not None else ''
    query = f"""
        SELECT {', '.join(select_exprs)}
        FROM read_csv(
            {_quote_string(tsv_path)},
            delim='\t',
            header=true,
            all_varchar=true,
            strict_mode=false,
            null_padding=true,
            ignore_errors=true
        )
        {limit_clause}
    """

    connection = duckdb.connect(':memory:')
    try:
        cursor = connection.execute(query)
        columns = [desc[0] for desc in cursor.description]
        while rows := cursor.fetchmany(batch_size):
            for row in rows:
                yield _normalize_row(dict(zip(columns, row)))
    finally:
        connection.close()


def _filtered_parquet_path(tsv_path: Path) -> Path:
    return tsv_path.with_name(
        f'{tsv_path.stem}_chemical_filter_v'
        f'{BINDINGDB_FILTERED_PARQUET_CACHE_VERSION}.parquet'
    )


def _ensure_filtered_parquet(
    tsv_path: Path,
    *,
    kwargs: dict[str, object],
) -> Path | None:
    if duckdb is None:
        raise ImportError('duckdb is required for BindingDB filtered Parquet.')
    lookup_path = chemical_resolver_lookup_path(kwargs)
    if not lookup_path.exists():
        print(
            f'Chemical resolver lookup not found; BindingDB filter disabled: '
            f'{lookup_path}',
            flush=True,
        )
        return None

    parquet_path = _filtered_parquet_path(tsv_path)
    if parquet_path.exists() and parquet_path.stat().st_mtime_ns >= max(
        tsv_path.stat().st_mtime_ns,
        lookup_path.stat().st_mtime_ns,
    ):
        return parquet_path

    tmp_path = parquet_path.with_suffix('.parquet.tmp')
    if tmp_path.exists():
        tmp_path.unlink()

    available_columns = _read_header(tsv_path)
    select_exprs = [
        (
            f'{_quote_identifier(column)} AS {_quote_identifier(column)}'
            if column in available_columns
            else f'NULL AS {_quote_identifier(column)}'
        )
        for column in _BINDINGDB_COLUMNS
    ]
    filter_sql = chemical_resolver_filter_sql(
        lookup_path=lookup_path,
        inchikey_expr=_quote_identifier('Ligand InChI Key'),
        sources=chemical_resolver_sources(kwargs),
    )
    print(f'Creating filtered BindingDB Parquet dataset: {parquet_path}', flush=True)
    connection = duckdb.connect(':memory:')
    try:
        connection.execute(
            f"""
            COPY (
              SELECT {', '.join(select_exprs)}
              FROM read_csv(
                {_quote_string(tsv_path)},
                delim='\t',
                header=true,
                all_varchar=true,
                strict_mode=false,
                null_padding=true,
                ignore_errors=true
              )
              WHERE {filter_sql}
            ) TO {_quote_string(tmp_path)}
            (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
    finally:
        connection.close()

    tmp_path.replace(parquet_path)
    return parquet_path


def _iter_csv_fallback(opener, max_lines: int | None = None) -> Generator[dict[str, str], None, None]:
    if not opener or not opener.result:
        return

    for file_handle in opener.result.values():
        reader = csv.DictReader(file_handle, delimiter='\t')
        for i, row in enumerate(reader):
            if max_lines is not None and i >= max_lines:
                break
            yield _normalize_row({column: row.get(column, '') for column in _BINDINGDB_COLUMNS})
        break


def _raw(
    opener,
    max_lines: int | None = None,
    use_duckdb: bool | None = None,
    batch_size: int = 50_000,
    **kwargs: object,
) -> Generator[dict[str, str], None, None]:
    """Parse BindingDB TSV rows.

    The default path streams rows from the downloaded archive with
    ``csv.DictReader``. Set ``use_duckdb=True`` or
    ``OMNIPATH_BINDINGDB_USE_DUCKDB=1`` to extract an on-disk TSV and let DuckDB
    stream a projected subset of columns into the bronze/preparse writer.
    """
    if use_duckdb is None:
        use_duckdb = _bindingdb_chemical_resolver_filter_enabled(kwargs) or (
            os.environ.get('OMNIPATH_BINDINGDB_USE_DUCKDB', '').lower()
            in {'1', 'true', 'yes'}
        )

    if use_duckdb:
        try:
            # Avoid extracting the full ~8 GB TSV for small max_lines smoke tests;
            # the CSV fallback can stream those directly from the zip handle.
            tsv_path = _bindingdb_tsv_path(opener, extract=max_lines is None)
            if tsv_path is not None:
                if (
                    max_lines is None
                    and _bindingdb_chemical_resolver_filter_enabled(kwargs)
                ):
                    parquet_path = _ensure_filtered_parquet(
                        tsv_path,
                        kwargs=kwargs,
                    )
                    if parquet_path is not None:
                        for row in iter_parquet(
                            path=parquet_path,
                            batch_size=batch_size,
                        ):
                            row = _normalize_row(row)
                            if _row_passes_filters(row, kwargs=kwargs):
                                yield row
                        return
                for row in _iter_duckdb_tsv(
                    tsv_path,
                    max_lines=max_lines,
                    batch_size=batch_size,
                ):
                    if _row_passes_filters(row, kwargs=kwargs):
                        yield row
                return
        except Exception as error:
            print(f'BindingDB DuckDB parser unavailable; falling back to csv. Reason: {error}', flush=True)

    for row in _iter_csv_fallback(opener, max_lines=max_lines):
        if _row_passes_filters(row, kwargs=kwargs):
            yield row
