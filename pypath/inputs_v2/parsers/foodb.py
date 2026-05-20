"""Parse FooDB data files.

This parser downloads and processes FooDB CSV files from a tar archive:
- Food.csv: Food information
- Compound.csv: Compound/metabolite information
- Content.csv: Food-compound relationships with concentrations
- CompoundSynonym.csv: Compound synonyms
- CompoundExternalDescriptor.csv: External database identifiers (ChEBI, KEGG)

It produces flat rows with member-list fields for use with MembersFromList.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Generator
import csv
from pathlib import Path
import tarfile
import tempfile
from typing import Protocol, TextIO

import polars as pl

from pypath.share.downloads import DATA_DIR

# Delimiter for member lists (using || as per reactome pattern)
MEMBER_DELIMITER = '||'

# Cache version retained only to avoid accidentally reusing older pickle files.
# The streaming parser below intentionally does not read or write pickle caches:
# historical FooDB caches were ~1.5 GB and can exceed worker memory before
# parsing even starts.
_CACHE_VERSION = 4


def _get_cache_path() -> Path:
    cache_dir = DATA_DIR / 'foodb'
    return cache_dir / f'foodb_data_v{_CACHE_VERSION}.pkl'


class _CsvWriter(Protocol):
    def writerow(self, row: list[str]) -> object: ...


def _load_cached_data(force_refresh: bool = False) -> list[dict[str, object]] | None:
    _ = force_refresh
    return None


def _save_cached_data(data: list[dict[str, object]]) -> None:
    _ = data


def _prepared_cache_available(
    *,
    force_refresh: bool = False,
    **_kwargs: object,
) -> bool:
    _ = force_refresh
    return False


def _csv_member_name(files: dict[str, object], csv_name: str) -> str | None:
    """Return archive member name for a CSV from opener.result keys."""
    for filename in files.keys():
        if filename.endswith(csv_name) or filename.endswith(f'/{csv_name}'):
            return filename
    return None


def _csv_path_from_archive(
    opener: object,
    files: dict[str, object],
    csv_name: str,
) -> Path | None:
    """Return an on-disk path for a CSV, extracting the tar member if necessary.

    Polars lazy `scan_csv` needs a path (not a tar member file object). We keep
    extraction selective via Download.needed and materialize only the required
    CSVs next to the downloaded archive.
    """
    archive_path_raw = getattr(opener, 'path', None)
    if not archive_path_raw:
        return None

    archive_path = Path(archive_path_raw)
    member_name = _csv_member_name(files, csv_name)
    if member_name is None:
        return None

    output_path = archive_path.parent / member_name
    if output_path.exists():
        return output_path

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, mode='r:') as archive:
        member = archive.getmember(member_name)
        src = archive.extractfile(member)
        if src is None:
            return None
        with output_path.open('wb') as dst:
            while chunk := src.read(1024 * 1024):
                dst.write(chunk)

    return output_path


def _clean_pl(col: str) -> pl.Expr:
    """Clean a string column lazily."""
    return pl.col(col).fill_null('').cast(pl.Utf8).str.strip_chars()


def _blank_to_dash_pl(col: str, alias: str | None = None) -> pl.Expr:
    """Polars expression equivalent of `_blank_to_dash`."""
    cleaned = _clean_pl(col)
    return (
        pl.when(cleaned == '')
        .then(pl.lit('-'))
        .otherwise(cleaned)
        .alias(alias or col)
    )


def _load_compounds_pl(path: Path) -> pl.DataFrame:
    """Load Compound.csv with Polars, selecting only columns used downstream."""
    fieldnames = [
        'id', 'public_id', 'name', 'state', 'annotation_quality',
        'description', 'cas_number', 'moldb_smiles', 'moldb_inchi',
        'moldb_mono_mass', 'moldb_inchikey', 'moldb_iupac',
        'kingdom', 'superklass', 'klass', 'subklass'
    ]
    needed_cols = [
        'id', 'public_id', 'name', 'cas_number', 'moldb_smiles',
        'moldb_mono_mass', 'moldb_inchikey', 'moldb_iupac',
        'kingdom', 'superklass', 'klass', 'subklass',
    ]
    return (
        pl.scan_csv(
            path,
            has_header=False,
            skip_rows=1,
            new_columns=fieldnames,
            infer_schema_length=1000,
        )
        .select(needed_cols)
        .with_columns(pl.col('id').cast(pl.Int64, strict=False))
        .filter(pl.col('id').is_not_null())
        .collect(engine='streaming')
    )


def _load_synonyms_pl(path: Path) -> pl.DataFrame:
    """Load CompoundSynonym.csv and aggregate the first five synonyms."""
    return (
        pl.scan_csv(path, infer_schema_length=1000)
        .select('source_id', 'source_type', 'synonym')
        .filter(pl.col('source_type') == 'Compound')
        .with_columns(
            pl.col('source_id').cast(pl.Int64, strict=False),
            _clean_pl('synonym').alias('synonym'),
        )
        .filter(pl.col('source_id').is_not_null() & (pl.col('synonym') != ''))
        .with_columns(pl.int_range(pl.len()).over('source_id').alias('_rank'))
        .filter(pl.col('_rank') < 5)
        .group_by('source_id')
        .agg(pl.col('synonym').str.join(delimiter=';'))
        .rename({'source_id': 'id', 'synonym': 'synonyms'})
        .collect(engine='streaming')
    )


def _load_external_ids_pl(path: Path) -> pl.DataFrame:
    """Load external descriptors and extract first ChEBI/KEGG IDs per compound."""
    df = (
        pl.scan_csv(path, infer_schema_length=1000)
        .select('compound_id', 'external_id')
        .with_columns(
            pl.col('compound_id').cast(pl.Int64, strict=False),
            _clean_pl('external_id').alias('external_id'),
        )
        .filter(pl.col('compound_id').is_not_null() & (pl.col('external_id') != ''))
        .collect(engine='streaming')
    )
    chebis = (
        df.filter(pl.col('external_id').str.starts_with('CHEBI:'))
        .unique(subset=['compound_id'], keep='first')
        .select(pl.col('compound_id').alias('id'), pl.col('external_id').alias('chebi'))
    )
    keggs = (
        df.filter(pl.col('external_id').str.contains(r'^C\d{5}'))
        .unique(subset=['compound_id'], keep='first')
        .select(pl.col('compound_id').alias('id'), pl.col('external_id').alias('kegg'))
    )
    return chebis.join(keggs, on='id', how='full', coalesce=True)


def _clean_scalar(value: object) -> str:
    """Return a stripped string, treating null-like values as empty."""
    if value is None:
        return ''
    return str(value).strip()


def _member_scalar(value: object) -> str:
    """Return the placeholder used for empty aligned member fields."""
    cleaned = _clean_scalar(value)
    return cleaned if cleaned else '-'


def _load_food_records(path: Path, food_cols: list[str]) -> list[dict[str, object]]:
    """Load the small Food.csv table as cleaned dictionaries."""
    df_food = (
        pl.scan_csv(path, infer_schema_length=1000)
        .select(food_cols)
        .with_columns(pl.col('id').cast(pl.Int64, strict=False))
        .filter(pl.col('id').is_not_null())
        .collect(engine='streaming')
    )
    records: list[dict[str, object]] = []
    for row in df_food.iter_rows(named=True):
        record = {'id': int(row['id'])}
        for col in food_cols:
            if col != 'id':
                record[col] = _clean_scalar(row.get(col))
        records.append(record)
    return records


def _load_compound_lookup(paths: dict[str, Path]) -> dict[int, dict[str, object]]:
    """Load compound metadata once; it is small compared with Content.csv."""
    df_full_compounds = (
        _load_compounds_pl(paths['Compound.csv'])
        .join(_load_synonyms_pl(paths['CompoundSynonym.csv']), on='id', how='left')
        .join(_load_external_ids_pl(paths['CompoundExternalDescriptor.csv']), on='id', how='left')
    )
    lookup: dict[int, dict[str, object]] = {}
    for row in df_full_compounds.iter_rows(named=True):
        compound_id = row.get('id')
        if compound_id is not None:
            lookup[int(compound_id)] = row
    return lookup


class _FoodMemberSpool:
    """Append member rows to per-food temporary CSV files with bounded handles."""

    def __init__(self, directory: Path, *, max_open: int | None = None) -> None:
        self.directory = directory
        self.max_open = max_open or _default_max_open_files()
        self._open: OrderedDict[int, tuple[TextIO, _CsvWriter]] = OrderedDict()

    def write(self, food_id: int, values: list[str]) -> None:
        writer = self._writer(food_id)
        writer.writerow(values)

    def path(self, food_id: int) -> Path:
        return self.directory / f'{food_id}.tsv'

    def close(self) -> None:
        for handle, _writer in self._open.values():
            handle.close()
        self._open.clear()

    def _writer(self, food_id: int) -> _CsvWriter:
        if food_id in self._open:
            handle, writer = self._open.pop(food_id)
            self._open[food_id] = (handle, writer)
            return writer

        if len(self._open) >= self.max_open:
            _old_food_id, (old_handle, _old_writer) = self._open.popitem(last=False)
            old_handle.close()

        handle = self.path(food_id).open('a', encoding='utf-8', newline='')
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        self._open[food_id] = (handle, writer)
        return writer


def _default_max_open_files() -> int:
    """Return a conservative per-food spool handle budget."""
    try:
        import resource

        soft_limit, _hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
    except (ImportError, OSError, ValueError):
        return 128

    if soft_limit == getattr(resource, 'RLIM_INFINITY', object()):
        return 2048
    return max(32, min(2048, int(soft_limit) - 64))


def _iter_content_member_rows(
    path: Path,
    compound_lookup: dict[int, dict[str, object]],
) -> Generator[tuple[int, list[str]], None, None]:
    """Stream Content.csv and yield prepared member-field rows."""
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get('source_type') != 'Compound':
                continue
            try:
                food_id = int(row.get('food_id') or '')
                source_id = int(row.get('source_id') or '')
            except ValueError:
                continue

            compound = compound_lookup.get(source_id, {})
            yield food_id, [
                str(source_id),
                _member_scalar(compound.get('public_id')),
                _member_scalar(compound.get('name')),
                _member_scalar(compound.get('cas_number')),
                _member_scalar(compound.get('moldb_inchikey')),
                _member_scalar(compound.get('moldb_smiles')),
                _member_scalar(compound.get('moldb_mono_mass')),
                _member_scalar(compound.get('moldb_iupac')),
                _member_scalar(compound.get('chebi')),
                _member_scalar(compound.get('kegg')),
                _member_scalar(compound.get('synonyms')),
                _member_scalar(compound.get('kingdom')),
                _member_scalar(compound.get('superklass')),
                _member_scalar(compound.get('klass')),
                _member_scalar(compound.get('subklass')),
                _member_scalar(row.get('orig_content')),
                _member_scalar(row.get('orig_min')),
                _member_scalar(row.get('orig_max')),
                _member_scalar(row.get('orig_unit')),
                _member_scalar(row.get('citation')),
                _member_scalar(row.get('orig_method')),
                _member_scalar(row.get('orig_food_part')),
                _member_scalar(row.get('preparation_type')),
            ]


def _read_member_lists(path: Path, member_cols: list[str]) -> dict[str, list[str]]:
    """Read one food's spooled members into list-valued mapper fields."""
    members = {col: [] for col in member_cols}
    if not path.exists():
        return members

    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            if len(row) < len(member_cols):
                row = [*row, *(['-'] * (len(member_cols) - len(row)))]
            for col, value in zip(member_cols, row, strict=False):
                members[col].append(value)
    return members


def _raw(
    opener: object,
    force_refresh: bool = False,
    **_kwargs: object,
) -> Generator[dict[str, object], None, None]:
    """Parse FooDB CSV files from tar archive and produce flat rows.

    Args:
        opener: The opener from Download.open() containing the extracted tar archive files
        force_refresh: If True, ignore cached data and reparse from source
    """
    # 0. Check cache
    cached_data = _load_cached_data(force_refresh=force_refresh)
    if cached_data is not None:
        yield from cached_data
        return

    # Validate opener
    opener_result = getattr(opener, 'result', None)
    if not opener or not opener_result:
        raise ValueError('FooDB download failed - no data available from opener')

    if not isinstance(opener_result, dict):
        raise ValueError(f'Expected dict of files from tar archive, got {type(opener_result)}')

    files = opener_result
    paths = {
        name: _csv_path_from_archive(opener, files, name)
        for name in (
            'Food.csv',
            'Content.csv',
            'Compound.csv',
            'CompoundSynonym.csv',
            'CompoundExternalDescriptor.csv',
        )
    }
    missing = [name for name, path in paths.items() if path is None]
    if missing:
        raise ValueError(f"FooDB required CSV(s) not found in archive: {', '.join(missing)}")

    member_cols = [
        'member_compound_id', 'member_compound_public_id', 'member_compound_name',
        'member_cas', 'member_inchikey', 'member_smiles', 'member_mass',
        'member_iupac', 'member_chebi', 'member_kegg', 'member_synonyms',
        'member_kingdom', 'member_superklass', 'member_klass', 'member_subklass',
        'member_content', 'member_min', 'member_max', 'member_unit',
        'member_citation', 'member_method', 'member_food_part', 'member_preparation'
    ]

    food_cols = [
        'id', 'public_id', 'name', 'name_scientific', 'description',
        'food_group', 'food_subgroup', 'ncbi_taxonomy_id'
    ]
    food_records = _load_food_records(paths['Food.csv'], food_cols)
    compound_lookup = _load_compound_lookup(paths)

    # Content.csv is hundreds of MB and is not ordered by food_id. Spool member
    # rows to disk first, then read and yield one food at a time. This avoids the
    # previous all-food group-by and the large delimiter-joined string columns.
    with tempfile.TemporaryDirectory(prefix='foodb_members_') as spool_root:
        spool = _FoodMemberSpool(Path(spool_root))
        try:
            for food_id, values in _iter_content_member_rows(
                paths['Content.csv'],
                compound_lookup,
            ):
                spool.write(food_id, values)
        finally:
            spool.close()

        for food in food_records:
            yield {
                **food,
                **_read_member_lists(spool.path(int(food['id'])), member_cols),
            }


_raw.prepared_cache_available = _prepared_cache_available
