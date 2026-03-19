"""
ChEMBL-specific data parsers.
"""

from __future__ import annotations

from typing import Any, Generator

from pypath.inputs_v2.parsers.base import iter_sqlite


def molecules_parser(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Parses ChEMBL molecules by joining dictionary, structures, and properties.

    Args:
        opener: The opener object for handling database connections.
        **kwargs: Additional arguments passed to iter_sqlite (e.g., sqlite_path, db_rel_path).

    Yields:
        Dictionary representing a joined molecule record.
    """
    query = """
    SELECT
        md.*,
        cs.canonical_smiles,
        cs.standard_inchi,
        cs.standard_inchi_key,
        cp.full_mwt,
        cp.full_molformula,
        cp.mw_monoisotopic,
        cp.alogp,
        cp.cx_logp,
        cp.cx_logd,
        cp.molecular_species
    FROM molecule_dictionary md
    LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
    LEFT JOIN compound_properties cp ON md.molregno = cp.molregno
    """
    yield from iter_sqlite(opener, query=query, **kwargs)
