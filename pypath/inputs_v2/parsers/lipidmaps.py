"""
LIPID MAPS Structure Database SDF parser.

Parses SDF format lipid structure files and yields flat dictionaries.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.formats.sdf import SdfReader


def _raw(opener, **_kwargs: object) -> Generator[dict[str, str], None, None]:
    """
    Parse LIPID MAPS SDF file.

    Args:
        opener: File opener from download_and_open

    Yields:
        Dictionary for each lipid with flattened structure and annotation fields
    """
    if not opener or not opener.result:
        return

    for file_handle in opener.result.values():
        sdf_reader = SdfReader(
            file_handle,
            names={
                'HMDB_ID': 'HMDB_ID',
                'PUBCHEM_CID': 'PUBCHEM_CID',
                'SWISSLIPIDS_ID': 'SWISSLIPIDS_ID',
                'LM_ID': 'LM_ID',
                'ABBREVIATION': 'ABBREVIATION',
                'CHEBI_ID': 'CHEBI_ID',
                'SYNONYMS': 'SYNONYMS',
                'INCHI': 'INCHI',
                'INCHI_KEY': 'INCHI_KEY',
                'COMMON_NAME': 'COMMON_NAME',
                'SYSTEMATIC_NAME': 'SYSTEMATIC_NAME',
                'SMILES': 'SMILES',
                'FORMULA': 'FORMULA',
            },
            fields={
                'EXACT_MASS',
                'CATEGORY',
                'MAIN_CLASS',
                'SUB_CLASS',
                'NAME',
            },
        )

        for record in sdf_reader:
            flat_record = {}
            if 'name' in record:
                flat_record.update(record['name'])
            if 'annot' in record:
                flat_record.update(record['annot'])
            yield flat_record
        break
