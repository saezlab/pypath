"""Stream PubChem SDF records into chemical resolver lookup rows.

PubChem is too large to treat as an in-memory table. This module parses SDF
records from gzip shards as text streams, keeps only CID and standard InChI key
fields, and writes normalized resolver rows in bounded parquet batches.
"""

from __future__ import annotations

import re
import gzip
import os
import time
import multiprocessing
from typing import BinaryIO
from pathlib import Path
from urllib.parse import urljoin, urlparse
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor, as_completed

import polars as pl
import pyarrow.parquet as pq

from pypath.share.downloads import download_and_open
from pypath.internals.cv_terms import (
    IdentifierNamespaceCv,
    cv_term_label_accession,
)
from omnipath_build.resolver.paths import ensure_chemicals_data_dir
from omnipath_build.resolver.parquet import write_parquet_from_dict_rows
from omnipath_build.resolver.identifier_types import (
    IDENTIFIER_TYPE_SCHEMA,
    identifier_type_id,
    identifier_type_rows,
)

PUBCHEM_CURRENT_SDF_BASE_URL = (
    'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/'
)
PUBCHEM_FIRST_COMPOUND_SDF_URL = (
    PUBCHEM_CURRENT_SDF_BASE_URL + 'Compound_000000001_000500000.sdf.gz'
)
PUBCHEM_COMPOUND_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.PUBCHEM_COMPOUND
)
STANDARD_INCHI_KEY_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.STANDARD_INCHI_KEY
)
PUBCHEM_IDENTIFIER_LOOKUP_SCHEMA: dict[str, pl.DataType] = {
    'key_identifier_type_id': pl.UInt32,
    'key_value': pl.Utf8,
    'canonical_identifier_type_id': pl.UInt32,
    'canonical_identifier': pl.Utf8,
}
PUBCHEM_IDENTIFIER_LOOKUP_OUTPUT_FILENAME = 'chemical_identifier_lookup.parquet'
IDENTIFIER_TYPE_OUTPUT_FILENAME = 'identifier_type.parquet'

_FIELD_RE = re.compile(r'^>\s*<([^>]+)>')
_PUBCHEM_SDF_FILENAME_RE = re.compile(r'Compound_\d{9}_\d{9}\.sdf\.gz')
_PUBCHEM_SUBFOLDER = 'pubchem'
_PUBCHEM_CURRENT_INDEX_FILENAME = 'current-full-sdf-index.html'
_PUBCHEM_SHARD_DIR = '_pubchem_shards'
_TARGET_FIELDS = {
    'PUBCHEM_COMPOUND_CID': 'pubchem_cid',
    'PUBCHEM_IUPAC_INCHIKEY': 'standard_inchi_key',
    'PUBCHEM_OPENEYE_INCHIKEY': 'standard_inchi_key',
}


def _clean(value: str | None) -> str | None:
    if value is None:
        return None
    text = value.strip()
    return text or None


def _row_from_record(record: dict[str, str]) -> dict | None:
    pubchem_cid = _clean(record.get('pubchem_cid'))
    standard_inchi_key = _clean(record.get('standard_inchi_key'))
    if not pubchem_cid or not standard_inchi_key:
        return None
    return {
        'key_type': PUBCHEM_COMPOUND_TYPE,
        'key_value': pubchem_cid,
        'standard_inchi_key': standard_inchi_key,
    }


def iter_pubchem_sdf_rows(lines: Iterable[str]) -> Iterable[dict]:
    """Yield PubChem resolver rows from SDF lines without retaining records."""

    record: dict[str, str] = {}
    active_field: str | None = None
    active_values: list[str] = []

    def flush_field() -> None:
        nonlocal active_field, active_values
        if (
            active_field is not None
            and active_values
            and active_field not in record
        ):
            record[active_field] = '\n'.join(active_values).strip()
        active_field = None
        active_values = []

    for raw_line in lines:
        line = raw_line.rstrip('\n\r')
        if line == '$$$$':
            flush_field()
            row = _row_from_record(record)
            if row is not None:
                yield row
            record = {}
            continue

        field_match = _FIELD_RE.match(line)
        if field_match:
            flush_field()
            active_field = _TARGET_FIELDS.get(field_match.group(1))
            continue

        if active_field is not None:
            if line == '':
                flush_field()
            else:
                active_values.append(line)

    flush_field()
    row = _row_from_record(record)
    if row is not None:
        yield row


def _iter_text_from_gzip_fileobj(fileobj: BinaryIO) -> Iterable[str]:
    with gzip.open(
        fileobj, mode='rt', encoding='utf-8', errors='replace'
    ) as handle:
        yield from handle


def _close_opener(opener: object) -> None:
    close = getattr(opener, 'close', None)
    if close is not None:
        close()


def _remove_bad_pubchem_cache(cache_path: Path, error: Exception) -> None:
    if cache_path.exists():
        cache_path.unlink()
        print(
            f'[pubchem] removed invalid cached gzip path={cache_path} '
            f'error={type(error).__name__}: {error}',
            flush=True,
        )


def iter_pubchem_sdf_gz_locations(
    source: str | Path | None = None,
    *,
    shard_count: int | None = None,
) -> Iterable[str | Path]:
    """Yield PubChem SDF gzip file locations.

    With no source, discover every current PubChem full-SDF shard from the NCBI
    listing. A source can be a single local/remote `.sdf.gz` file. A shard
    count limits the discovered PubChem list before any shard is streamed.
    """

    if shard_count is not None and shard_count < 1:
        raise ValueError('PubChem shard count must be at least 1.')

    if source is None:
        emitted = 0
        for location in _iter_current_pubchem_sdf_urls():
            yield location
            emitted += 1
            if shard_count is not None and emitted >= shard_count:
                break
        return

    if shard_count not in (None, 1):
        raise ValueError(
            'PubChem shard count only applies to discovered PubChem shards. '
            'Use pubchem_url for an explicit single shard.'
        )

    location = str(source)
    if location.startswith(('http://', 'https://', 'ftp://')):
        _ensure_sdf_gz(location)
        yield location
        return

    path = Path(location)
    _ensure_sdf_gz(path.name)
    yield path


def _iter_current_pubchem_sdf_urls() -> Iterable[str]:
    opener = download_and_open(
        url=PUBCHEM_CURRENT_SDF_BASE_URL,
        filename=_PUBCHEM_CURRENT_INDEX_FILENAME,
        subfolder=_PUBCHEM_SUBFOLDER,
        large=False,
        encoding='utf-8',
        default_mode='r',
        force_download=True,
    )

    try:
        html = opener.result
    finally:
        opener.close()

    filenames = sorted(set(_PUBCHEM_SDF_FILENAME_RE.findall(html)))
    if not filenames:
        raise ValueError(
            f'No PubChem SDF shards found at {PUBCHEM_CURRENT_SDF_BASE_URL}'
        )
    print(
        f'[pubchem] discovered_shards={len(filenames)} '
        f'base={PUBCHEM_CURRENT_SDF_BASE_URL}',
        flush=True,
    )
    for filename in filenames:
        yield urljoin(PUBCHEM_CURRENT_SDF_BASE_URL, filename)


def _ensure_sdf_gz(location: str) -> None:
    if not location.endswith('.sdf.gz'):
        raise ValueError(
            f'Expected a single PubChem SDF .gz file, got {location!r}'
        )


def iter_pubchem_sdf_gz_rows(path_or_url: str | Path) -> Iterable[dict]:
    """Stream PubChem SDF gzip rows from a local path or URL."""

    location = str(path_or_url)
    if location.startswith(('http://', 'https://', 'ftp://')):
        yield from _iter_downloaded_pubchem_sdf_gz_rows(location)

        return

    with Path(location).open('rb') as handle:
        yield from iter_pubchem_sdf_rows(_iter_text_from_gzip_fileobj(handle))


def _iter_downloaded_pubchem_sdf_gz_rows(
    location: str,
    *,
    retry_bad_cache: bool = True,
) -> Iterable[dict]:
    opener = download_and_open(
        url=location,
        filename=Path(location).name,
        subfolder=_PUBCHEM_SUBFOLDER,
        large=True,
        encoding='utf-8',
        default_mode='r',
        ext='gz',
    )
    cache_path = Path(opener.path)

    try:
        yield from iter_pubchem_sdf_rows(opener.result)
    except (gzip.BadGzipFile, EOFError) as error:
        _close_opener(opener)
        opener = None
        _remove_bad_pubchem_cache(cache_path, error)
        if retry_bad_cache:
            yield from _iter_downloaded_pubchem_sdf_gz_rows(
                location,
                retry_bad_cache=False,
            )
            return
        raise
    finally:
        if opener is not None:
            _close_opener(opener)


def iter_pubchem_compound_rows(
    source: str | Path | None = None,
    *,
    filter_inchikeys: frozenset[str] | None = None,
    shard_count: int | None = None,
    progress_every: int = 100_000,
) -> Iterable[dict]:
    """Stream PubChem resolver rows from all selected SDF gzip shards."""

    locations = iter_pubchem_sdf_gz_locations(
        source,
        shard_count=shard_count,
    )
    for shard_index, location in enumerate(locations, start=1):
        yield from iter_pubchem_location_rows(
            location,
            shard_index=shard_index,
            filter_inchikeys=filter_inchikeys,
            progress_every=progress_every,
        )


def iter_pubchem_location_rows(
    location: str | Path,
    *,
    shard_index: int,
    filter_inchikeys: frozenset[str] | None = None,
    progress_every: int = 100_000,
) -> Iterable[dict]:
    """Stream PubChem resolver rows from one SDF gzip shard."""

    filter_label = 'enabled' if filter_inchikeys is not None else 'disabled'
    filter_key_count = (
        len(filter_inchikeys) if filter_inchikeys is not None else 0
    )
    started = time.monotonic()
    rows_seen = 0
    matched_rows = 0
    completed = False
    print(
        f'[pubchem] shard={shard_index} start '
        f'location={location} filter={filter_label} '
        f'filter_keys={filter_key_count}',
        flush=True,
    )
    try:
        for row in iter_pubchem_sdf_gz_rows(location):
            rows_seen += 1
            if (
                filter_inchikeys is not None
                and row.get('standard_inchi_key') not in filter_inchikeys
            ):
                if progress_every > 0 and rows_seen % progress_every == 0:
                    _log_pubchem_shard_progress(
                        shard_index,
                        rows_seen,
                        matched_rows,
                        started,
                    )
                continue
            matched_rows += 1
            if progress_every > 0 and rows_seen % progress_every == 0:
                _log_pubchem_shard_progress(
                    shard_index,
                    rows_seen,
                    matched_rows,
                    started,
                )
            yield row
        completed = True
    finally:
        status = 'done' if completed else 'stopped'
        print(
            f'[pubchem] shard={shard_index} {status} '
            f'rows_seen={rows_seen} matched_rows={matched_rows} '
            f'elapsed={time.monotonic() - started:.1f}s '
            f'location={location}',
            flush=True,
        )


def _log_pubchem_shard_progress(
    shard_index: int,
    rows_seen: int,
    matched_rows: int,
    started: float,
) -> None:
    print(
        f'[pubchem] shard={shard_index} progress '
        f'rows_seen={rows_seen} matched_rows={matched_rows} '
        f'elapsed={time.monotonic() - started:.1f}s',
        flush=True,
    )


def materialize_pubchem_compound_sdf(
    output_dir: str | Path | None = None,
    *,
    source: str | Path | None = None,
    max_records: int | None = None,
    filter_inchikeys: frozenset[str] | None = None,
    pubchem_shards: int | None = None,
) -> dict[str, int]:
    """Write resolver parquet rows from PubChem SDF shard URLs or paths."""

    output_dir = (
        Path(output_dir)
        if output_dir is not None
        else ensure_chemicals_data_dir()
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = iter_pubchem_compound_rows(
        source,
        filter_inchikeys=filter_inchikeys,
        shard_count=pubchem_shards,
    )
    if max_records is not None:
        rows = _take(rows, max_records)
    rows = _normalized_pubchem_rows(rows)
    row_count = write_parquet_from_dict_rows(
        rows,
        PUBCHEM_IDENTIFIER_LOOKUP_SCHEMA,
        output_dir / PUBCHEM_IDENTIFIER_LOOKUP_OUTPUT_FILENAME,
    )
    type_count = write_parquet_from_dict_rows(
        identifier_type_rows({PUBCHEM_COMPOUND_TYPE, STANDARD_INCHI_KEY_TYPE}),
        IDENTIFIER_TYPE_SCHEMA,
        output_dir / IDENTIFIER_TYPE_OUTPUT_FILENAME,
    )
    return {
        'chemical_identifier_lookup_rows': row_count,
        'identifier_type_rows': type_count,
    }


def materialize_pubchem_compound_shards(
    output_dir: str | Path,
    *,
    source: str | Path | None = None,
    filter_inchikeys: frozenset[str] | None = None,
    pubchem_shards: int | None = None,
    jobs: int = 1,
) -> list[tuple[Path, int]]:
    """Write one normalized lookup parquet per PubChem shard.

    Completed shard files are reused on reruns. A shard is considered complete
    only if both the parquet file and its ``.done`` sidecar exist.
    """

    if jobs < 1:
        raise ValueError('PubChem jobs must be at least 1.')

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    shard_dir = output_dir / _PUBCHEM_SHARD_DIR
    shard_dir.mkdir(parents=True, exist_ok=True)

    locations = list(
        iter_pubchem_sdf_gz_locations(
            source,
            shard_count=pubchem_shards,
        )
    )
    tasks = [
        _pubchem_shard_task(
            shard_dir,
            shard_index,
            location,
            filter_inchikeys=filter_inchikeys,
        )
        for shard_index, location in enumerate(locations, start=1)
    ]
    if not tasks:
        return []

    workers = min(jobs, len(tasks))
    print(
        f'[pubchem] shard_jobs={workers} shard_count={len(tasks)} '
        f'shard_dir={shard_dir}',
        flush=True,
    )
    if workers == 1:
        results = [_write_pubchem_shard_parquet(task) for task in tasks]
    else:
        results = []
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=_process_context(),
        ) as executor:
            futures = {
                executor.submit(_write_pubchem_shard_parquet, task): task
                for task in tasks
            }
            for future in as_completed(futures):
                results.append(future.result())

    results.sort(key=lambda item: item[0])
    return [(path, row_count) for _, path, row_count in results]


def iter_pubchem_lookup_parquet_rows(
    paths: Iterable[str | Path],
) -> Iterable[dict]:
    """Yield normalized PubChem lookup rows from shard parquet files."""

    for path in paths:
        parquet_file = pq.ParquetFile(path)
        for batch in parquet_file.iter_batches():
            for row in batch.to_pylist():
                yield row


def materialize_pubchem_first_compound_sdf(
    output_dir: str | Path | None = None,
    *,
    url: str = PUBCHEM_FIRST_COMPOUND_SDF_URL,
    max_records: int | None = None,
) -> dict[str, int]:
    """Write resolver rows from the first current PubChem compound shard."""

    return materialize_pubchem_compound_sdf(
        output_dir,
        source=url,
        max_records=max_records,
    )


def _take(rows: Iterable[dict], max_records: int) -> Iterable[dict]:
    emitted = 0
    for row in rows:
        if emitted >= max_records:
            break
        yield row
        emitted += 1


def _normalized_pubchem_rows(rows: Iterable[dict]) -> Iterable[dict]:
    key_identifier_type_id = identifier_type_id(PUBCHEM_COMPOUND_TYPE)
    canonical_identifier_type_id = identifier_type_id(STANDARD_INCHI_KEY_TYPE)
    for row in rows:
        standard_inchi_key = row.get('standard_inchi_key')
        yield {
            'key_identifier_type_id': key_identifier_type_id,
            'key_value': row.get('key_value'),
            'canonical_identifier_type_id': canonical_identifier_type_id,
            'canonical_identifier': standard_inchi_key,
        }


def _pubchem_shard_task(
    shard_dir: Path,
    shard_index: int,
    location: str | Path,
    *,
    filter_inchikeys: frozenset[str] | None,
) -> dict:
    return {
        'shard_dir': str(shard_dir),
        'shard_index': shard_index,
        'location': str(location),
        'filter_inchikeys': filter_inchikeys,
    }


def _process_context() -> multiprocessing.context.BaseContext | None:
    methods = multiprocessing.get_all_start_methods()
    if 'fork' in methods:
        return multiprocessing.get_context('fork')
    return None


def _write_pubchem_shard_parquet(task: dict) -> tuple[int, Path, int]:
    shard_dir = Path(task['shard_dir'])
    shard_index = int(task['shard_index'])
    location = task['location']
    filter_inchikeys = task.get('filter_inchikeys')
    output_path = _pubchem_shard_path(shard_dir, shard_index, location)
    done_path = output_path.with_suffix(output_path.suffix + '.done')
    existing_count = _read_done_count(done_path)

    if output_path.exists() and existing_count is not None:
        print(
            f'[pubchem] shard={shard_index} skip existing '
            f'rows={existing_count} path={output_path}',
            flush=True,
        )
        return shard_index, output_path, existing_count

    tmp_path = output_path.with_suffix(
        output_path.suffix + f'.tmp.{os.getpid()}'
    )
    done_tmp_path = done_path.with_suffix(
        done_path.suffix + f'.tmp.{os.getpid()}'
    )
    tmp_path.unlink(missing_ok=True)
    done_tmp_path.unlink(missing_ok=True)

    try:
        rows = _normalized_pubchem_rows(
            iter_pubchem_location_rows(
                location,
                shard_index=shard_index,
                filter_inchikeys=filter_inchikeys,
            )
        )
        row_count = write_parquet_from_dict_rows(
            rows,
            PUBCHEM_IDENTIFIER_LOOKUP_SCHEMA,
            tmp_path,
        )
        done_tmp_path.write_text(
            f'{row_count}\nlocation={location}\n',
            encoding='utf-8',
        )
        tmp_path.replace(output_path)
        done_tmp_path.replace(done_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        done_tmp_path.unlink(missing_ok=True)
        raise

    return shard_index, output_path, row_count


def _pubchem_shard_path(
    shard_dir: Path,
    shard_index: int,
    location: str | Path,
) -> Path:
    location_text = str(location)
    parsed_path = urlparse(location_text).path
    name = Path(parsed_path or location_text).name
    safe_name = re.sub(r'[^A-Za-z0-9_.-]+', '_', name)
    if safe_name.endswith('.sdf.gz'):
        safe_name = safe_name[:-7]
    return shard_dir / f'{shard_index:05d}_{safe_name}.parquet'


def _read_done_count(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        first_line = path.read_text(encoding='utf-8').splitlines()[0]
        return int(first_line)
    except (IndexError, OSError, ValueError):
        return None
