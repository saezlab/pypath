"""
Parse FooDB data files.

This parser downloads and processes FooDB CSV files from a tar archive:
- Food.csv: Food information
- Compound.csv: Compound/metabolite information  
- Content.csv: Food-compound relationships with concentrations
- CompoundSynonym.csv: Compound synonyms
- CompoundExternalDescriptor.csv: External database identifiers (ChEBI, KEGG)

It produces flat rows with delimited member fields for use with MembersFromList.
Data is cached in a pickle file to speed up subsequent access.
"""

from __future__ import annotations

import io
import os
import pickle
import tarfile
from collections.abc import Generator
from pathlib import Path
from typing import Any

import pandas as pd
import polars as pl
import numpy as np

from pypath.share.downloads import DATA_DIR


# Delimiter for member lists (using || as per reactome pattern)
MEMBER_DELIMITER = '||'

# Module-level cache for parsed data
_DATA_CACHE: list[dict] | None = None

# Cache version to invalidate older pickled formats. The streaming parser below
# no longer writes this cache; bumping avoids loading old ~1.5 GB pickle files.
_CACHE_VERSION = 3


def _get_cache_path() -> Path:
    cache_dir = DATA_DIR / 'foodb'
    return cache_dir / f'foodb_data_v{_CACHE_VERSION}.pkl'


def _load_cached_data(force_refresh: bool = False) -> list[dict] | None:
    global _DATA_CACHE
    if force_refresh:
        return None
    if _DATA_CACHE is not None:
        return _DATA_CACHE
    
    pickle_path = _get_cache_path()
    if pickle_path.exists():
        try:
            with open(pickle_path, 'rb') as f:
                data = pickle.load(f)
            _DATA_CACHE = data
            return data
        except Exception:
            pass
    return None


def _save_cached_data(data: list[dict]) -> None:
    global _DATA_CACHE
    _DATA_CACHE = data
    pickle_path = _get_cache_path()
    try:
        pickle_path.parent.mkdir(parents=True, exist_ok=True)
        with open(pickle_path, 'wb') as f:
            pickle.dump(data, f)
    except Exception:
        pass


def _prepared_cache_available(
    *,
    force_refresh: bool = False,
    **_kwargs: Any,
) -> bool:
    return not force_refresh and _get_cache_path().exists()


def _clean_str(s: pd.Series) -> pd.Series:
    """Clean string series: replace NaN with empty string, strip whitespace."""
    return s.fillna('').astype(str).str.strip()


def _blank_to_dash(s: pd.Series) -> pd.Series:
    """
    Convert empty strings to '-' to preserve alignment in delimited lists.
    The tabular_builder filters out '-' as empty.
    """
    # Ensure it's cleaned first
    s = _clean_str(s)
    # Replace empty with '-'
    return s.replace('', '-')


def _find_csv_file(files: dict[str, Any], csv_name: str) -> Any | None:
    """
    Find a CSV file in the extracted archive files dict.

    Keep archive member handles streaming. The previous implementation read each
    member into a StringIO; for FooDB's 744 MB Content.csv this alone can push a
    machine over memory limits before pandas starts parsing.
    """
    for filename, handle in files.items():
        # Match by basename (handles nested paths in tar archives)
        if filename.endswith(csv_name) or filename.endswith(f'/{csv_name}'):
            if isinstance(handle, bytes):
                return io.BytesIO(handle)
            if isinstance(handle, str):
                return io.StringIO(handle)
            return handle
    return None


def _csv_member_name(files: dict[str, Any], csv_name: str) -> str | None:
    """Return archive member name for a CSV from opener.result keys."""
    for filename in files.keys():
        if filename.endswith(csv_name) or filename.endswith(f'/{csv_name}'):
            return filename
    return None


def _csv_path_from_archive(opener: Any, files: dict[str, Any], csv_name: str) -> Path | None:
    """
    Return an on-disk path for a CSV, extracting the tar member if necessary.

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


def _load_compounds(files: dict[str, Any]) -> pd.DataFrame:
    """Load Compound.csv with correct headers."""
    handle = _find_csv_file(files, 'Compound.csv')
    if handle is None:
        raise ValueError("Compound.csv not found in archive")
    
    # Correct column order (header in file is mismatched with data in the original file)
    # We skip the header row (skiprows=1) and provide our own names.
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

    # Use C engine for speed, handle potential bad lines if any
    df = pd.read_csv(
        handle,
        names=fieldnames,
        skiprows=1,
        usecols=needed_cols,
        dtype={'id': 'Int64'}, # nullable int
        low_memory=False
    )
    
    # Filter valid IDs
    df = df.dropna(subset=['id'])
    
    return df


def _load_synonyms(files: dict[str, Any]) -> pd.DataFrame:
    """Load CompoundSynonym.csv and aggregate synonyms."""
    handle = _find_csv_file(files, 'CompoundSynonym.csv')
    if handle is None:
        raise ValueError("CompoundSynonym.csv not found in archive")
    
    df = pd.read_csv(handle, low_memory=False)
    
    # Filter for compounds
    df = df[df['source_type'] == 'Compound'].copy()
    
    # Clean synonyms
    df['synonym'] = _clean_str(df['synonym'])
    df = df[df['synonym'] != '']
    
    # Aggregate: take top 5, join with ';'
    df_head = df.groupby('source_id').head(5)
    
    synonyms = df_head.groupby('source_id')['synonym'].agg(';'.join).reset_index()
    synonyms.rename(columns={'source_id': 'id', 'synonym': 'synonyms'}, inplace=True)
    
    return synonyms


def _load_external_ids(files: dict[str, Any]) -> pd.DataFrame:
    """Load External Descriptors and extract ChEBI/KEGG."""
    handle = _find_csv_file(files, 'CompoundExternalDescriptor.csv')
    if handle is None:
        raise ValueError("CompoundExternalDescriptor.csv not found in archive")
    
    df = pd.read_csv(handle, low_memory=False)
    
    # Clean
    df['external_id'] = _clean_str(df['external_id'])
    df = df[df['external_id'] != '']
    df['compound_id'] = pd.to_numeric(df['compound_id'], errors='coerce')
    df = df.dropna(subset=['compound_id'])
    
    # Identify types
    # ChEBI
    chebi_mask = df['external_id'].str.startswith('CHEBI:')
    chebis = df[chebi_mask].copy()
    # Deduplicate: take first per compound
    chebis = chebis.drop_duplicates(subset=['compound_id'])
    chebis = chebis[['compound_id', 'external_id']].rename(columns={'external_id': 'chebi'})
    
    # KEGG: C followed by 5 digits
    kegg_mask = df['external_id'].str.match(r'^C\d{5}')
    keggs = df[kegg_mask].copy()
    keggs = keggs.drop_duplicates(subset=['compound_id'])
    keggs = keggs[['compound_id', 'external_id']].rename(columns={'external_id': 'kegg'})
    
    # Merge
    exts = pd.merge(chebis, keggs, on='compound_id', how='outer')
    exts.rename(columns={'compound_id': 'id'}, inplace=True)
    
    return exts


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


def _raw(
    opener,
    force_refresh: bool = False,
    **_kwargs,
) -> Generator[dict[str, Any], None, None]:
    """
    Parse FooDB CSV files from tar archive and produce flat rows with delimited member fields.
    
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
    if not opener or not opener.result:
        raise ValueError("FooDB download failed - no data available from opener")
    
    if not isinstance(opener.result, dict):
        raise ValueError(f"Expected dict of files from tar archive, got {type(opener.result)}")
    
    files = opener.result
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
    df_food = (
        pl.scan_csv(paths['Food.csv'], infer_schema_length=1000)
        .select(food_cols)
        .with_columns(pl.col('id').cast(pl.Int64, strict=False))
        .filter(pl.col('id').is_not_null())
        .with_columns(_clean_pl(col).alias(col) for col in food_cols if col != 'id')
        .collect(engine='streaming')
    )

    df_full_compounds = (
        _load_compounds_pl(paths['Compound.csv'])
        .join(_load_synonyms_pl(paths['CompoundSynonym.csv']), on='id', how='left')
        .join(_load_external_ids_pl(paths['CompoundExternalDescriptor.csv']), on='id', how='left')
    )

    content_cols = [
        'food_id', 'source_id', 'source_type', 'orig_content', 'orig_min',
        'orig_max', 'orig_unit', 'citation', 'orig_method', 'orig_food_part',
        'preparation_type',
    ]
    lf_merged = (
        pl.scan_csv(paths['Content.csv'], infer_schema_length=1000)
        .select(content_cols)
        .filter(pl.col('source_type') == 'Compound')
        .with_columns(
            pl.col('food_id').cast(pl.Int64, strict=False),
            pl.col('source_id').cast(pl.Int64, strict=False),
        )
        .filter(pl.col('food_id').is_not_null() & pl.col('source_id').is_not_null())
        .join(df_full_compounds.lazy(), left_on='source_id', right_on='id', how='left')
        .with_columns(
            pl.col('source_id').cast(pl.Int64).cast(pl.Utf8).alias('member_compound_id'),
            _blank_to_dash_pl('public_id', 'member_compound_public_id'),
            _blank_to_dash_pl('name', 'member_compound_name'),
            _blank_to_dash_pl('cas_number', 'member_cas'),
            _blank_to_dash_pl('moldb_inchikey', 'member_inchikey'),
            _blank_to_dash_pl('moldb_smiles', 'member_smiles'),
            _blank_to_dash_pl('moldb_mono_mass', 'member_mass'),
            _blank_to_dash_pl('moldb_iupac', 'member_iupac'),
            _blank_to_dash_pl('chebi', 'member_chebi'),
            _blank_to_dash_pl('kegg', 'member_kegg'),
            _blank_to_dash_pl('synonyms', 'member_synonyms'),
            _blank_to_dash_pl('kingdom', 'member_kingdom'),
            _blank_to_dash_pl('superklass', 'member_superklass'),
            _blank_to_dash_pl('klass', 'member_klass'),
            _blank_to_dash_pl('subklass', 'member_subklass'),
            _blank_to_dash_pl('orig_content', 'member_content'),
            _blank_to_dash_pl('orig_min', 'member_min'),
            _blank_to_dash_pl('orig_max', 'member_max'),
            _blank_to_dash_pl('orig_unit', 'member_unit'),
            _blank_to_dash_pl('citation', 'member_citation'),
            _blank_to_dash_pl('orig_method', 'member_method'),
            _blank_to_dash_pl('orig_food_part', 'member_food_part'),
            _blank_to_dash_pl('preparation_type', 'member_preparation'),
        )
        .select('food_id', *member_cols)
    )

    df_members = (
        lf_merged
        .group_by('food_id')
        .agg(pl.col(col).str.join(delimiter=MEMBER_DELIMITER) for col in member_cols)
        .collect(engine='streaming')
    )

    df_final = df_food.join(df_members, left_on='id', right_on='food_id', how='left')
    df_final = df_final.with_columns(pl.col(col).fill_null('') for col in member_cols)

    # Yield records directly; do not write/load the old ~1.5 GB pickle cache.
    yield from df_final.iter_rows(named=True)


_raw.prepared_cache_available = _prepared_cache_available
