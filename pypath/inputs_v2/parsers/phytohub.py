"""
Parse PhytoHub data pages.

This parser reads and merges the PhytoHub foods and compounds HTML tables.

It produces two kinds of rows:
- _raw: one row per food, with ||-delimited member fields for its compounds
- _raw_compounds: one row per compound (including those without a food source)
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_html_tables, read_text

# Delimiter for member lists (using || as per reactome pattern)
MEMBER_DELIMITER = '||'

# Output key -> column name in the PhytoHub compounds table
COMPOUND_COLUMNS = {
    'compound_id': 'Identifier',
    'primary_key': 'Primary key',
    'compound_name': 'Chemical Name',
    'systematic_name': 'Systematic name',
    'iupac_name': 'IUPAC Name',
    'synonyms': 'Synonyms',
    'cas_number': 'CAS Number',
    'formula': 'Formula',
    'monoisotopic_mass': 'Monoisotopic mass',
    'average_mass': 'Average mass',
    'inchi_key': 'InChI Key',
    'inchi': 'Inchi',
    'smiles': 'Smiles',
    'phytochemical_family': 'Phytochemical family',
    'phytochemical_class': 'Phytochemical class',
    'phytochemical_subclass': 'Phytochemical subclass',
    'metabolite_family': 'Metabolite family',
    'metabolite_class': 'Metabolite class',
    'metabolite_subclass': 'Metabolite subclass',
    'solubility': 'Solubility',
    'predicted_logp': 'Predicted alogps logp',
    'predicted_logs': 'Predicted alogps logs',
    'pka_basic': 'Pka strongest basic',
    'pka_acidic': 'Pka strongest acidic',
    'metabolites_precursor': 'Metabolites / Precursor',
    'intake_biomarkers': 'Intake biomarkers',
}

# Output key -> member field name (cas_number is exposed as member_cas)
MEMBER_KEYS = {
    key: ('member_cas' if key == 'cas_number' else f'member_{key}')
    for key in COMPOUND_COLUMNS
}


def _clean_value(v: Any) -> str:
    """Clean a value, returning empty string if missing or 'Not Available'."""
    if v is None:
        return ''
    text = str(v).strip()
    if not text or text.lower() in ('not available', 'nan', 'n/a'):
        return ''
    return text


def _get_handle(opener):
    """Return the file handle from an opener result, or None."""
    if not opener or not opener.result:
        return None
    return (
        next(iter(opener.result.values()), None)
        if isinstance(opener.result, dict)
        else opener.result
    )


def _table_rows(
    opener: Any,
    table_id: str,
    required: tuple[str, ...],
) -> Generator[dict[str, str], None, None]:
    """Read a named PhytoHub table, validating its headers and row shapes."""
    if not _get_handle(opener):
        return

    for table in iter_html_tables(read_text(opener)):
        if table['attrs'].get('id') != table_id:
            continue

        rows = table['rows']
        if not rows:
            break
        headers = [cell['text'] for cell in rows[0]]
        missing = [name for name in required if name not in headers]
        if missing:
            raise ValueError(
                f'PhytoHub table {table_id!r} lacks expected columns: '
                f'{", ".join(missing)}'
            )

        for row in rows[1:]:
            if len(row) != len(headers):
                raise ValueError(
                    f'PhytoHub table {table_id!r} row has {len(row)} cells; '
                    f'expected {len(headers)}'
                )
            record = {
                header: cell['text']
                for header, cell in zip(headers, row, strict=True)
                if header
            }
            if any(record.values()):
                yield record
        return

    raise ValueError(f'PhytoHub page has no table with id {table_id!r}')


def _load_compounds_by_food(
    download_compounds: Any,
) -> dict[str, list[dict[str, str]]]:
    """Read the compounds table and group compounds by their food sources."""
    compounds_by_food: dict[str, list[dict[str, str]]] = defaultdict(list)

    for row in _table_rows(
        download_compounds.open(),
        'entries',
        ('Identifier', 'Chemical Name', 'Food sources'),
    ):
        food_sources_raw = _clean_value(row.get('Food sources', ''))
        if not food_sources_raw:
            continue

        compound_data = {
            key: _clean_value(row.get(column, ''))
            for key, column in COMPOUND_COLUMNS.items()
        }

        # Associate this compound with each food listed in 'Food sources'
        for food in food_sources_raw.split(' / '):
            food_name = food.strip()
            if food_name:
                compounds_by_food[food_name.lower()].append(compound_data)

    return compounds_by_food


def _build_member_fields(
    compounds: list[dict[str, str]],
    delimiter: str = MEMBER_DELIMITER,
) -> dict[str, str]:
    """
    Build delimited member fields from compound data.

    Empty values are replaced by '-' to preserve index alignment.
    The Column._normalize_token in tabular_builder filters out '-' as empty.
    """
    EMPTY = '-'

    if not compounds:
        return {member_key: '' for member_key in MEMBER_KEYS.values()}

    return {
        member_key: delimiter.join(c.get(key, '') or EMPTY for c in compounds)
        for key, member_key in MEMBER_KEYS.items()
    }


def _raw(
    opener,
    download_compounds,
    **_kwargs,
) -> Generator[dict[str, Any], None, None]:
    """
    Parse the PhytoHub foods table into rows with delimited member fields.
    """
    handle = _get_handle(opener)
    if not handle:
        return

    # Load compounds grouped by food name
    compounds_by_food = _load_compounds_by_food(download_compounds)

    for row in _table_rows(
        opener,
        'foods-index',
        ('Name', 'Food id', 'Group'),
    ):
        food_name = _clean_value(row.get('Name', ''))
        compounds = compounds_by_food.get(food_name.lower(), [])
        food_group, _, food_subgroup = _clean_value(
            row.get('Group', '')
        ).partition(', ')

        yield {
            'id': _clean_value(row.get('Food id', '')),
            'name': food_name,
            'synonyms': _clean_value(row.get('Synonyms', '')),
            'food_group': food_group,
            'food_subgroup': food_subgroup,
            'food_phytochemicals_count': _clean_value(
                row.get('Food phytochemicals', '')
            ),
            'intake_biomarkers': _clean_value(row.get('Intake biomarkers', '')),
            **_build_member_fields(compounds),
        }


def _raw_compounds(
    opener,
    **_kwargs,
) -> Generator[dict[str, Any], None, None]:
    """
    Parse the PhytoHub compounds table into one flat row per compound.

    Includes compounds without food sources (e.g. human/animal metabolites).
    Keys are the original table headers, with values cleaned.
    """
    for row in _table_rows(
        opener,
        'entries',
        ('Identifier', 'Chemical Name'),
    ):
        yield {k: _clean_value(v) for k, v in row.items() if k}
