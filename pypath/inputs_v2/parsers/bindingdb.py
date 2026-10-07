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
    'Number of Protein Chains in Target (>1 implies a multichain complex)',
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
    'BindingDB Target Chain 1 Sequence',
    'UniProt (SwissProt) Primary ID of Target Chain 1',
    'UniProt (SwissProt) Recommended Name of Target Chain 1',
    'UniProt (TrEMBL) Primary ID of Target Chain 1',
    'UniProt (TrEMBL) Submitted Name of Target Chain 1',
]


# Preserve explicitly reported multichain targets instead of assigning their
# measurements to chain 1 alone.
_BINDINGDB_COLUMNS.extend(
    stem + str(index)
    for index in range(2, 51)
    for stem in (
        'UniProt (SwissProt) Primary ID of Target Chain ',
        'UniProt (TrEMBL) Primary ID of Target Chain ',
        'UniProt (SwissProt) Recommended Name of Target Chain ',
        'UniProt (TrEMBL) Submitted Name of Target Chain ',
    )
)

def _normalize_row(row: dict[str, str | None]) -> dict[str, str]:
    normalized = {key: '' if value is None else str(value) for key, value in row.items()}
    for key, sequence in row.items():
        match = re.fullmatch(r'BindingDB Target Chain\s*(\d+)?\s*Sequence(?:\s*(\d+))?', key)
        if match and sequence:
            normalized[f'BindingDB Target Chain {match[1] or match[2] or "1"} Sequence'] = str(sequence)
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
    return _row_has_allowed_pchembl(row, kwargs=kwargs)


def _quote_identifier(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _quote_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _read_header(tsv_path: Path) -> set[str]:
    with tsv_path.open('r', encoding='utf-8', errors='replace', newline='') as handle:
        header = handle.readline().rstrip('\r\n')
    return set(header.split('\t')) if header else set()


def _selected_columns(header: object) -> list[str]:
    """Keep all reported chains, including sequence-column spelling variants."""
    columns = list(_BINDINGDB_COLUMNS)
    for column in header or []:
        if re.fullmatch(
            r'BindingDB Target Chain\s*(\d+)?\s*Sequence(?:\s*(\d+))?', column
        ) or re.fullmatch(
            r'UniProt \((?:SwissProt|TrEMBL)\) '
            r'(?:Primary ID|Recommended Name|Submitted Name) '
            r'of Target Chain \d+',
            column,
        ):
            if column not in columns:
                columns.append(column)
    return columns


def _bindingdb_tsv_path(opener, *, extract: bool = True) -> Path | None:
    """Return an on-disk TSV path, extracting the zip member if necessary.

    DuckDB's CSV reader operates on paths. The downloaded BindingDB archive is a
    zip containing a single very large TSV, so we extract it in streaming chunks
    next to the zip. This costs disk space but keeps memory bounded and is
    reused by later parser runs.
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
        for column in _selected_columns(sorted(available_columns))
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


def _iter_csv_fallback(opener, max_lines: int | None = None) -> Generator[dict[str, str], None, None]:
    if not opener or not opener.result:
        return

    for file_handle in opener.result.values():
        reader = csv.DictReader(file_handle, delimiter='\t')
        for i, row in enumerate(reader):
            if max_lines is not None and i >= max_lines:
                break
            yield _normalize_row({column: row.get(column, '')
                                  for column in _selected_columns(reader.fieldnames)})
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
    stream a projected subset of columns.
    """
    if use_duckdb is None:
        use_duckdb = (
            os.environ.get('OMNIPATH_BINDINGDB_USE_DUCKDB', '').lower()
            in {'1', 'true', 'yes'}
        )

    if use_duckdb:
        try:
            # Avoid extracting the full ~8 GB TSV for small max_lines smoke tests;
            # the CSV fallback can stream those directly from the zip handle.
            tsv_path = _bindingdb_tsv_path(opener, extract=max_lines is None)
            if tsv_path is not None:
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


_SEQUENCE_COLUMN_RE = re.compile(r'BindingDB Target Chain\s*(\d+)?\s*Sequence(?:\s*(\d+))?')
_AFFINITY_SQL_RE = r'([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)'


def _chembl_runs(value: str | None) -> str | None:
    """``_normalize_row``'s ChEMBL fix (a lookahead regex, which DuckDB lacks)."""
    if value is None:
        return None
    stripped = value.strip()
    if (
        stripped and stripped.count('CHEMBL') > 1
        and '::' not in stripped and ';' not in stripped and '|' not in stripped
    ):
        return _CHEMBL_RUN_RE.sub(r'\1::', stripped)
    return value


def raw_table(db, name: str, opener, max_lines: int | None = None, **kwargs: object) -> int:
    """The rows of :func:`_raw` as DuckDB table ``name``, built in SQL.

    Same columns, values, filter and order as the ``csv`` row parser: selected
    columns ('' when empty or missing), chain sequences under their normalized
    name, ChEMBL runs split, ``pchembl_value`` formatted as ``%.6g`` and filtered.
    """
    tsv_path = _bindingdb_tsv_path(opener, extract=True)
    with tsv_path.open('r', encoding='utf-8', errors='replace', newline='') as handle:
        header = handle.readline().rstrip('\r\n').split('\t')
    columns = _selected_columns(header)
    q = _quote_identifier
    present = set(header)
    select = [f"coalesce({q(c)}, '') AS {q(c)}" if c in present else f"'' AS {q(c)}" for c in columns]
    # Chain sequences: the last non-empty source column wins, as in _normalize_row.
    targets: dict[str, list[str]] = {}
    for column in columns:
        if match := _SEQUENCE_COLUMN_RE.fullmatch(column):
            target = f'BindingDB Target Chain {match[1] or match[2] or "1"} Sequence'
            targets.setdefault(target, []).append(column)
    fallback = {t: ("''" if t in columns else 'NULL') for t in targets}
    sequence = {
        t: f"coalesce({', '.join(f'nullif({q(s)}, {chr(39)}{chr(39)})' for s in reversed(sources))}, {fallback[t]})"
        for t, sources in targets.items()
    }
    affinities = [
        f"""CASE WHEN trim({q(c)}) = '' OR starts_with(trim({q(c)}), '>') THEN NULL
            ELSE 9 - log10(nullif(greatest(TRY_CAST(nullif(regexp_extract(trim({q(c)}), '{_AFFINITY_SQL_RE}', 1), '') AS DOUBLE), 0), 0)) END"""
        for c in _BINDINGDB_AFFINITY_COLUMNS
    ]
    function = f'{name}_chembl_runs'
    db.create_function(function, _chembl_runs, ['VARCHAR'], 'VARCHAR', null_handling='special')
    limit = f'LIMIT {int(max_lines)}' if max_lines is not None else ''
    db.execute('SET preserve_insertion_order=true')
    db.execute(f"""CREATE OR REPLACE TEMP TABLE {name}_read AS SELECT {', '.join(select)}
        FROM read_csv({_quote_string(tsv_path)}, delim='\t', header=true, all_varchar=true, quote='"',
                      escape='"', strict_mode=false, null_padding=true) {limit}""")
    chembl = q('ChEMBL ID of Ligand')
    replaced = {chembl: f"CASE WHEN length({chembl}) - length(replace({chembl}, 'CHEMBL', '')) > 6 THEN {function}({chembl}) ELSE {chembl} END"}
    outputs = [c for c in columns if c not in targets] + list(targets)
    expressions = [sequence.get(c) or replaced.get(q(c)) or q(c) for c in outputs]
    pchembl = f"greatest({', '.join(affinities)})"
    min_pchembl = _bindingdb_min_pchembl(kwargs)
    keep = (
        f"pchembl_value <> '' AND pchembl_value::DOUBLE > {min_pchembl!r}"
        if _bindingdb_pchembl_filter_enabled(kwargs) else 'true'
    )
    db.execute(f"""CREATE OR REPLACE TABLE {name} AS
        SELECT row_number() OVER (ORDER BY read_order) - 1 AS rid, * EXCLUDE (read_order) FROM (
            SELECT rowid AS read_order, {', '.join(f'{e} AS {q(c)}' for e, c in zip(expressions, outputs))},
                   coalesce(printf('%.6g', {pchembl}), '') AS pchembl_value
            FROM {name}_read) WHERE {keep}""")
    db.execute(f'DROP TABLE {name}_read')
    db.remove_function(function)
    db.execute('SET preserve_insertion_order=false')
    return db.execute(f'SELECT count(*) FROM {name}').fetchone()[0]
