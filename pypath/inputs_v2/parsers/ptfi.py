"""
Parse PTFI (Protein-Food Interaction Database) data.

This parser processes PTFI JSON files downloaded from ptfidiscover.markerlab.com.
It extracts:
- Food specimen metadata (FOODON ID, name, taxonomy)
- Compound composition data with chemical formulas
- Statistical measurements (mean, median, min, max concentrations)

The data is structured as food entities with compound members, similar to FooDB.
"""

from __future__ import annotations

import json
import re
from collections.abc import Generator
from pathlib import Path
from typing import Any

# Delimiter for member lists (using || as per pypath standard)
MEMBER_DELIMITER = '||'


def parse_compound_name(name: str) -> dict[str, str]:
    """
    Parse compound name to extract chemical formula and PTFI ID.

    Example: "C28H34O4-470 (PTF28874)" -> {'formula': 'C28H34O4', 'ptfi_id': 'PTF28874'}
    """
    result = {}

    # Extract PTFI ID
    ptfi_match = re.search(r'PTF\d+', name)
    if ptfi_match:
        result['ptfi_id'] = ptfi_match.group()

    # Extract chemical formula
    formula_match = re.search(r'^([A-Z][a-z]?\d*)+', name)
    if formula_match:
        result['formula'] = formula_match.group()

    return result


def flatten_compounds(compounds: list, parent_category: str = None) -> Generator[dict[str, Any], None, None]:
    """
    Recursively flatten the hierarchical compound structure.

    Yields compound data with their parent category information.
    """
    for compound in compounds:
        level = compound.get('level', 0)

        # Level 0 = category, skip but process children
        if level == 0:
            category_name = compound.get('name', 'Unknown')
            if 'children' in compound and compound['children']:
                yield from flatten_compounds(compound['children'], category_name)
        else:
            # Level 1+ = actual compound with data
            analyte = compound.get('analyte', {})

            # Parse compound name for identifiers
            name = analyte.get('name', compound.get('name', ''))
            parsed = parse_compound_name(name)

            compound_data = {
                'name': name,
                'ptfi_id': parsed.get('ptfi_id', ''),
                'formula': parsed.get('formula', ''),
                'category': parent_category or compound.get('parent', ''),
                'value': str(compound.get('value', '')),
                'units': analyte.get('units', ''),
                # Statistical data
                'mean': str(analyte.get('group_mean', '')),
                'median': str(analyte.get('group_median', '')),
                'min': str(analyte.get('group_min', '')),
                'max': str(analyte.get('group_max', '')),
                'n_samples': str(analyte.get('group_n', '')),
            }

            yield compound_data

            # Process nested children if any
            if 'children' in compound and compound['children']:
                yield from flatten_compounds(compound['children'], parent_category)


def _blank_to_dash(s: str) -> str:
    """
    Convert empty strings to '-' to preserve alignment in delimited lists.
    The tabular_builder filters out '-' as empty.
    """
    return '-' if s == '' or s is None else str(s)


def _raw(
    data_dir: str,
    **_kwargs,
) -> Generator[dict[str, Any], None, None]:
    """
    Parse PTFI JSON files and produce flat rows with delimited member fields.

    Args:
        data_dir: Path to directory containing PTFI JSON files

    Yields:
        Dictionary records with food metadata and delimited compound member lists
    """
    data_path = Path(data_dir)

    if not data_path.exists():
        raise FileNotFoundError(f"PTFI data directory not found: {data_dir}")

    # Process each JSON file in the directory
    json_files = sorted(data_path.glob('*.json'))

    if not json_files:
        raise ValueError(f"No JSON files found in {data_dir}")

    for json_file in json_files:
        with open(json_file) as f:
            data = json.load(f)

        # Extract specimen metadata
        specimen_id = data.get('id', '')
        name = data.get('name', '')
        num_samples = data.get('num_samples', 0)

        # Extract FOODON ID from parents (most specific one)
        parents = data.get('parents', [])
        foodon_id = ''
        for parent in parents:
            if parent.get('id', '').startswith('FOODON_'):
                foodon_id = parent['id']

        # Extract all compounds from both arrays
        all_compounds = []
        all_compounds.extend(flatten_compounds(data.get('compounds', [])))
        all_compounds.extend(flatten_compounds(data.get('otherCompounds', [])))

        # Build delimited member fields
        if all_compounds:
            member_fields = {
                'member_ptfi_id': [],
                'member_name': [],
                'member_formula': [],
                'member_category': [],
                'member_value': [],
                'member_units': [],
                'member_mean': [],
                'member_median': [],
                'member_min': [],
                'member_max': [],
                'member_n_samples': [],
            }

            for compound in all_compounds:
                member_fields['member_ptfi_id'].append(_blank_to_dash(compound['ptfi_id']))
                member_fields['member_name'].append(_blank_to_dash(compound['name']))
                member_fields['member_formula'].append(_blank_to_dash(compound['formula']))
                member_fields['member_category'].append(_blank_to_dash(compound['category']))
                member_fields['member_value'].append(_blank_to_dash(compound['value']))
                member_fields['member_units'].append(_blank_to_dash(compound['units']))
                member_fields['member_mean'].append(_blank_to_dash(compound['mean']))
                member_fields['member_median'].append(_blank_to_dash(compound['median']))
                member_fields['member_min'].append(_blank_to_dash(compound['min']))
                member_fields['member_max'].append(_blank_to_dash(compound['max']))
                member_fields['member_n_samples'].append(_blank_to_dash(compound['n_samples']))

            # Join all member fields with delimiter
            member_data = {
                key: MEMBER_DELIMITER.join(values)
                for key, values in member_fields.items()
            }
        else:
            # No compounds - empty member fields
            member_data = {
                'member_ptfi_id': '',
                'member_name': '',
                'member_formula': '',
                'member_category': '',
                'member_value': '',
                'member_units': '',
                'member_mean': '',
                'member_median': '',
                'member_min': '',
                'member_max': '',
                'member_n_samples': '',
            }

        # Yield combined record
        yield {
            'specimen_id': specimen_id,
            'name': name,
            'num_samples': str(num_samples),
            'foodon_id': foodon_id,
            **member_data,
        }
