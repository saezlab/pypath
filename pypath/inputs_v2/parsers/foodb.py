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
import pickle
from collections.abc import Generator
from pathlib import Path
from typing import Any

import pandas as pd
import numpy as np

from pypath.share.downloads import DATA_DIR


# Delimiter for member lists (using || as per reactome pattern)
MEMBER_DELIMITER = '||'

# Module-level cache for parsed data
_DATA_CACHE: list[dict] | None = None

# Cache version to invalidate older pickled formats
_CACHE_VERSION = 2  # Bumped version since we changed the parsing approach


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


def _find_csv_file(files: dict[str, Any], csv_name: str) -> io.StringIO | None:
    """
    Find a CSV file in the extracted archive files dict.
    
    Args:
        files: Dict mapping filename to file handle from opener.result
        csv_name: Name of the CSV file to find (e.g., 'Food.csv')
        
    Returns:
        A StringIO or file handle for the CSV file, or None if not found
    """
    for filename, handle in files.items():
        # Match by basename (handles nested paths in tar archives)
        if filename.endswith(csv_name) or filename.endswith(f'/{csv_name}'):
            # If handle is bytes, decode it
            if hasattr(handle, 'read'):
                content = handle.read()
                if isinstance(content, bytes):
                    return io.StringIO(content.decode('utf-8'))
                elif isinstance(content, str):
                    return io.StringIO(content)
                # Reset handle if readable
                handle.seek(0)
                return handle
    return None


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
    
    # Use C engine for speed, handle potential bad lines if any
    df = pd.read_csv(
        handle,
        names=fieldnames,
        skiprows=1,
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

    # 1. Load Food
    food_handle = _find_csv_file(files, 'Food.csv')
    if food_handle is None:
        raise ValueError("Food.csv not found in archive")
    
    df_food = pd.read_csv(food_handle, dtype={'id': 'Int64'}, low_memory=False)
    df_food = df_food.dropna(subset=['id'])
    
    # Keep only necessary food columns and clean them
    food_cols = [
        'id', 'public_id', 'name', 'name_scientific', 'description',
        'food_group', 'food_subgroup', 'ncbi_taxonomy_id'
    ]
    # Ensure columns exist (Food.csv should have them)
    existing_food_cols = [c for c in food_cols if c in df_food.columns]
    df_food = df_food[existing_food_cols]
    
    for col in existing_food_cols:
        if col != 'id':
            df_food[col] = _clean_str(df_food[col])
            
    # 2. Load Content (The link between Food and Compound)
    content_handle = _find_csv_file(files, 'Content.csv')
    if content_handle is None:
        raise ValueError("Content.csv not found in archive")
    
    df_content = pd.read_csv(content_handle, low_memory=False)
    
    # Filter for source_type == 'Compound'
    df_content = df_content[df_content['source_type'] == 'Compound']
    # Ensure IDs are numeric
    df_content['food_id'] = pd.to_numeric(df_content['food_id'], errors='coerce')
    df_content['source_id'] = pd.to_numeric(df_content['source_id'], errors='coerce') # this is compound_id
    df_content = df_content.dropna(subset=['food_id', 'source_id'])
    
    # 3. Load Compound Info
    df_compounds = _load_compounds(files)
    df_synonyms = _load_synonyms(files)
    df_ext = _load_external_ids(files)
    
    # Merge Compound details
    # Compound base + Synonyms + Ext IDs
    # df_compounds has 'id'
    df_full_compounds = df_compounds.merge(df_synonyms, on='id', how='left')
    df_full_compounds = df_full_compounds.merge(df_ext, on='id', how='left')
    
    # 4. Merge Content with Compounds
    df_merged = df_content.merge(
        df_full_compounds, 
        left_on='source_id', 
        right_on='id', 
        how='left',
        suffixes=('_content', '_compound')
    )
    
    # 5. Prepare data for aggregation
    # Helper to safe get column
    def get_col(df, col):
        return df[col] if col in df.columns else pd.Series([''] * len(df))

    # Compound fields
    df_merged['member_compound_id'] = df_merged['source_id'].astype(str)
    df_merged['member_compound_public_id'] = get_col(df_merged, 'public_id')
    df_merged['member_compound_name'] = get_col(df_merged, 'name')
    df_merged['member_cas'] = get_col(df_merged, 'cas_number')
    df_merged['member_inchikey'] = get_col(df_merged, 'moldb_inchikey')
    df_merged['member_smiles'] = get_col(df_merged, 'moldb_smiles')
    df_merged['member_mass'] = get_col(df_merged, 'moldb_mono_mass')
    df_merged['member_iupac'] = get_col(df_merged, 'moldb_iupac')
    df_merged['member_chebi'] = get_col(df_merged, 'chebi')
    df_merged['member_kegg'] = get_col(df_merged, 'kegg')
    df_merged['member_synonyms'] = get_col(df_merged, 'synonyms')
    df_merged['member_kingdom'] = get_col(df_merged, 'kingdom')
    df_merged['member_superklass'] = get_col(df_merged, 'superklass')
    df_merged['member_klass'] = get_col(df_merged, 'klass')
    df_merged['member_subklass'] = get_col(df_merged, 'subklass')
    
    # Content fields
    df_merged['member_content'] = get_col(df_merged, 'orig_content')
    df_merged['member_min'] = get_col(df_merged, 'orig_min')
    df_merged['member_max'] = get_col(df_merged, 'orig_max')
    df_merged['member_unit'] = get_col(df_merged, 'orig_unit')
    df_merged['member_citation'] = get_col(df_merged, 'citation')
    df_merged['member_method'] = get_col(df_merged, 'orig_method')
    df_merged['member_food_part'] = get_col(df_merged, 'orig_food_part')
    df_merged['member_preparation'] = get_col(df_merged, 'preparation_type')
    
    # List of all member columns to aggregate
    member_cols = [
        'member_compound_id', 'member_compound_public_id', 'member_compound_name',
        'member_cas', 'member_inchikey', 'member_smiles', 'member_mass',
        'member_iupac', 'member_chebi', 'member_kegg', 'member_synonyms',
        'member_kingdom', 'member_superklass', 'member_klass', 'member_subklass',
        'member_content', 'member_min', 'member_max', 'member_unit',
        'member_citation', 'member_method', 'member_food_part', 'member_preparation'
    ]
    
    # Apply _blank_to_dash to all member columns to ensure alignment
    for col in member_cols:
        df_merged[col] = _blank_to_dash(df_merged[col])
        
    # 6. Aggregate by food_id
    # We group by food_id and join all member columns with delimiter
    agg_funcs = {col: lambda x: MEMBER_DELIMITER.join(x) for col in member_cols}
    
    df_members = df_merged.groupby('food_id').agg(agg_funcs).reset_index()
    
    # 7. Merge back with Food info
    df_final = df_food.merge(df_members, left_on='id', right_on='food_id', how='left')
    
    # For foods with no members, the member columns will be NaN. Fill with empty string.
    for col in member_cols:
        df_final[col] = df_final[col].fillna('')
        
    # 8. Save to cache and yield
    # Materialize list of dicts
    final_data = df_final.to_dict('records')
    
    _save_cached_data(final_data)
    
    yield from final_data
