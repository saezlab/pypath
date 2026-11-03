"""
Parse Phenol-Explorer data files.

This parser downloads and merges multiple Phenol-Explorer data files:
- compounds.csv: Base compound information
- compounds-structures.csv: SMILES structures
- foods.csv: Food information
- composition-data.xlsx: Food-compound relationships with concentrations

It produces flat rows with delimited member fields for use with MembersFromList.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Generator
import csv
import io
from typing import Any

from pypath.share.downloads import download_and_open


# Delimiter for member lists (using || as per reactome pattern)
MEMBER_DELIMITER = '||'


def _is_valid_value(v: Any) -> bool:
    """Check if value is valid (not None, not NaN, not empty string)."""
    if v is None:
        return False
    if isinstance(v, float):
        import math
        if math.isnan(v):
            return False
    if isinstance(v, str) and (not v.strip() or v.lower() == 'nan'):
        return False
    return True


def _clean_value(v: Any) -> str:
    """Clean a value for use in delimited string, returning empty string if invalid."""
    if not _is_valid_value(v):
        return ''
    return str(v)


def _format_chebi(v: str) -> str:
    """Format CHEBI ID with prefix if needed."""
    if not v:
        return ''
    v = str(v).strip()
    if not v:
        return ''
    if v.startswith('CHEBI:'):
        return v
    return f'CHEBI:{v}'


def _clean_pubmed(pubmed_str: Any) -> str:
    """Clean pubmed IDs, removing NaN values."""
    if not _is_valid_value(pubmed_str):
        return ''
    s = str(pubmed_str)
    # Filter out 'nan' values and clean
    pmids = [p.strip() for p in s.split(';') if p.strip() and p.strip().lower() != 'nan']
    return ';'.join(pmids)


def _build_compound_lookup(
    download_compounds,
    download_compounds_structures,
) -> dict[str, dict[str, Any]]:
    """Build a lookup dict from compound name to compound details."""
    compounds_by_name: dict[str, dict[str, Any]] = {}
    
    compounds_opener = download_compounds.open()
    if compounds_opener and compounds_opener.result:
        handle = next(iter(compounds_opener.result.values()), None) if isinstance(compounds_opener.result, dict) else compounds_opener.result
        if handle:
            for row in csv.DictReader(handle):
                name = row.get('name', '')
                if name:
                    compounds_by_name[name] = {
                        'id': row.get('id', ''),
                        'chebi_id': row.get('chebi_id', ''),
                        'pubchem_compound_id': row.get('pubchem_compound_id', ''),
                        'cas_number': row.get('cas_number', ''),
                        'molecular_weight': row.get('molecular_weight', ''),
                        'formula': row.get('formula', ''),
                        'synonyms': row.get('synonyms', ''),
                        'aglycones': row.get('aglycones', ''),
                        'compound_class': row.get('compound_class', ''),
                        'compound_subclass': row.get('compound_subclass', ''),
                    }
    
    # Enrich with SMILES from structures file
    structures_opener = download_compounds_structures.open()
    if structures_opener and structures_opener.result:
        handle = next(iter(structures_opener.result.values()), None) if isinstance(structures_opener.result, dict) else structures_opener.result
        if handle:
            for row in csv.DictReader(handle):
                name = row.get('name', '')
                if name and name in compounds_by_name:
                    compounds_by_name[name]['smiles'] = row.get('smiles', '')
    
    return compounds_by_name


def _load_composition_data(
    download_composition,
    compound_lookup: dict[str, dict[str, Any]],
) -> dict[str, list[dict]]:
    """Load composition data and enrich with compound details."""
    import pandas as pd
    
    composition_by_food: dict[str, list[dict]] = defaultdict(list)
    
    composition_opener = download_composition.open()
    if composition_opener and composition_opener.result:
        comp_handle = next(iter(composition_opener.result.values()), None) if isinstance(composition_opener.result, dict) else composition_opener.result
        if comp_handle:
            content = comp_handle.read()
            df = pd.read_excel(io.BytesIO(content))
            for _, row in df.iterrows():
                food_name = row.get('food', '')
                compound_name = row.get('compound', '')
                if food_name and compound_name:
                    compound_details = compound_lookup.get(compound_name, {})
                    
                    composition_by_food[food_name].append({
                        'compound_name': compound_name,
                        'compound_id': compound_details.get('id', ''),
                        'chebi_id': _format_chebi(compound_details.get('chebi_id', '')),
                        'pubchem_compound_id': compound_details.get('pubchem_compound_id', ''),
                        'cas_number': compound_details.get('cas_number', ''),
                        'smiles': compound_details.get('smiles', ''),
                        'molecular_weight': compound_details.get('molecular_weight', ''),
                        'formula': compound_details.get('formula', ''),
                        'synonyms': compound_details.get('synonyms', ''),
                        'aglycones': compound_details.get('aglycones', ''),
                        'compound_class': row.get('compound_group', '') or compound_details.get('compound_class', ''),
                        'compound_subclass': row.get('compound_sub_group', '') or compound_details.get('compound_subclass', ''),
                        'mean': _clean_value(row.get('mean')),
                        'min': _clean_value(row.get('min')),
                        'max': _clean_value(row.get('max')),
                        'sd': _clean_value(row.get('sd')),
                        'units': _clean_value(row.get('units', '')),
                        'n': _clean_value(row.get('n')),
                        'N': _clean_value(row.get('N')),
                        'experimental_method': _clean_value(row.get('experimental_method_group', '')),
                        'pubmed_ids': _clean_pubmed(row.get('pubmed_ids', '')),
                    })
    
    return composition_by_food


def _build_member_fields(compositions: list[dict], delimiter: str = MEMBER_DELIMITER) -> dict[str, str]:
    """
    Build delimited member fields from composition data.
    
    IMPORTANT: We use '-' for empty values to preserve index alignment.
    The Column._normalize_token in tabular_builder filters out '-' as empty.
    """
    EMPTY = '-'
    
    if not compositions:
        return {
            'member_compound_id': '',
            'member_compound_name': '',
            'member_chebi': '',
            'member_pubchem': '',
            'member_cas': '',
            'member_smiles': '',
            'member_formula': '',
            'member_synonyms': '',
            'member_compound_class': '',
            'member_compound_subclass': '',
            'member_molecular_weight': '',
            'member_aglycones': '',
            'member_mean': '',
            'member_min': '',
            'member_max': '',
            'member_sd': '',
            'member_units': '',
            'member_n': '',
            'member_N': '',
            'member_experimental_method': '',
            'member_pubmed': '',
        }
    
    def _v(val: str) -> str:
        """Return value or placeholder if empty."""
        return val if val else EMPTY
    
    return {
        'member_compound_id': delimiter.join(_v(c.get('compound_id', '')) for c in compositions),
        'member_compound_name': delimiter.join(_v(c.get('compound_name', '')) for c in compositions),
        'member_chebi': delimiter.join(_v(c.get('chebi_id', '')) for c in compositions),
        'member_pubchem': delimiter.join(_v(c.get('pubchem_compound_id', '')) for c in compositions),
        'member_cas': delimiter.join(_v(c.get('cas_number', '')) for c in compositions),
        'member_smiles': delimiter.join(_v(c.get('smiles', '')) for c in compositions),
        'member_formula': delimiter.join(_v(c.get('formula', '')) for c in compositions),
        'member_synonyms': delimiter.join(_v(c.get('synonyms', '')) for c in compositions),
        'member_compound_class': delimiter.join(_v(c.get('compound_class', '')) for c in compositions),
        'member_compound_subclass': delimiter.join(_v(c.get('compound_subclass', '')) for c in compositions),
        'member_molecular_weight': delimiter.join(_v(c.get('molecular_weight', '')) for c in compositions),
        'member_aglycones': delimiter.join(_v(c.get('aglycones', '')) for c in compositions),
        'member_mean': delimiter.join(_v(c.get('mean', '')) for c in compositions),
        'member_min': delimiter.join(_v(c.get('min', '')) for c in compositions),
        'member_max': delimiter.join(_v(c.get('max', '')) for c in compositions),
        'member_sd': delimiter.join(_v(c.get('sd', '')) for c in compositions),
        'member_units': delimiter.join(_v(c.get('units', '')) for c in compositions),
        'member_n': delimiter.join(_v(c.get('n', '')) for c in compositions),
        'member_N': delimiter.join(_v(c.get('N', '')) for c in compositions),
        'member_experimental_method': delimiter.join(_v(c.get('experimental_method', '')) for c in compositions),
        'member_pubmed': delimiter.join(_v(c.get('pubmed_ids', '')) for c in compositions),
    }


def _raw(
    opener,
    download_compounds,
    download_compounds_structures,
    download_composition,
    **_kwargs,
) -> Generator[dict[str, Any], None, None]:
    """
    Parse foods CSV and produce flat rows with delimited member fields.
    
    Output row structure:
    - Food fields (id, name, food_group, food_subgroup, etc.)
    - Member fields as ||-delimited strings:
      - member_compound_id, member_compound_name, member_chebi, member_smiles, etc.
      - membership annotation fields: member_mean, member_min, member_max, etc.
    """
    if not opener or not opener.result:
        return
    
    handle = next(iter(opener.result.values()), None) if isinstance(opener.result, dict) else opener.result
    if not handle:
        return
    
    # Build compound lookup for enriching members
    compound_lookup = _build_compound_lookup(download_compounds, download_compounds_structures)
    
    # Load composition data (Excel file)
    composition_by_food = _load_composition_data(download_composition, compound_lookup)
    
    for row in csv.DictReader(handle):
        food_name = row.get('name', '')
        compositions = composition_by_food.get(food_name, [])
        
        # Build delimited member fields
        member_fields = _build_member_fields(compositions)
        
        # Merge food row with member fields
        output_row = {
            'id': row.get('id', ''),
            'name': food_name,
            'food_group': row.get('food_group', ''),
            'food_subgroup': row.get('food_subgroup', ''),
            'scientific_name': row.get('food_source_scientific_name', ''),
            'botanical_family': row.get('food_source_botanical_family', ''),
            **member_fields,
        }
        yield output_row
