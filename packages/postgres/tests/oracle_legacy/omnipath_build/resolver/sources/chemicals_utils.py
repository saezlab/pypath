"""utils-backed chemical resolver export (the chemical repoint).

The native chemical resolver streams a capped/sharded PubChem SDF; this module
instead pulls the **authoritative** chemical structure mappings from the
omnipath-utils Postgres (`omnipath_utils.resolver_chemical`: full PubChem
cid→InChIKey + UniChem/ChEBI-bridged cross-refs) and writes them into the same
``chemicals/lookup/`` parquet partition the build already globs. Read via DuckDB
ATTACH and streamed straight to parquet (``COPY … TO``) so the ~124M rows never
materialise in Python — mirrors the global UniProt slice (``proteins.py``).

Requires ``OMNIPATH_BUILD_UTILS_PG_URL`` (the utils Postgres). When unset the
export is skipped (the native PubChem source still provides a capped fallback).
"""

from __future__ import annotations

from pathlib import Path

from omnipath_build.resolver.paths import ensure_chemicals_data_dir
from omnipath_build.resolver.identifier_types import identifier_type_id
from omnipath_build.resolver.sources.chemicals import (
    CHEBI_TYPE,
    CHEMBL_COMPOUND_TYPE,
    CHEMICAL_IDENTIFIER_LOOKUP_PARTITION_DIRNAME,
    HMDB_TYPE,
    KEGG_COMPOUND_TYPE,
    LIPIDMAPS_TYPE,
    PUBCHEM_COMPOUND_TYPE,
    STANDARD_INCHI_KEY_TYPE,
    SWISSLIPIDS_TYPE,
)
from omnipath_build.resolver.sources.proteins import _utils_pg_url

UTILS_LOOKUP_FILENAME = 'utils.parquet'

#: utils ``resolver_chemical.source_type`` -> build identifier type label. Only
#: source types with a build identifier type are exported (others are dropped).
_SOURCE_TYPE_TO_BUILD_TYPE: dict[str, str] = {
    'pubchem': PUBCHEM_COMPOUND_TYPE,
    'chembl': CHEMBL_COMPOUND_TYPE,
    'chebi': CHEBI_TYPE,
    'hmdb': HMDB_TYPE,
    'lipidmaps': LIPIDMAPS_TYPE,
    'swisslipids': SWISSLIPIDS_TYPE,
    'kegg': KEGG_COMPOUND_TYPE,
}


def materialize_utils_chemicals(
    output_dir: str | Path | None = None,
    *,
    skip_existing: bool = True,
) -> dict[str, int]:
    """Export ``omnipath_utils.resolver_chemical`` → ``lookup/utils.parquet``.

    Returns ``{'utils_chemical_rows': n}``. Raises if the utils Postgres URL is
    unset — callers wrap this and warn+continue (the native sources remain).
    """

    import duckdb

    base = ensure_chemicals_data_dir() if output_dir is None else Path(output_dir)
    base.mkdir(parents=True, exist_ok=True)
    lookup_dir = base / CHEMICAL_IDENTIFIER_LOOKUP_PARTITION_DIRNAME
    lookup_dir.mkdir(parents=True, exist_ok=True)
    parquet_path = lookup_dir / UTILS_LOOKUP_FILENAME
    if skip_existing and parquet_path.exists():
        return {'utils_chemical_rows': -1}  # -1 = pre-existing, not recomputed

    type_map = {
        source_type: identifier_type_id(build_type)
        for source_type, build_type in _SOURCE_TYPE_TO_BUILD_TYPE.items()
    }
    inchikey_type_id = identifier_type_id(STANDARD_INCHI_KEY_TYPE)
    case_sql = 'CASE source_type ' + ' '.join(
        f"WHEN '{source_type}' THEN {type_id}"
        for source_type, type_id in type_map.items()
    ) + ' END'
    in_list = ', '.join(f"'{source_type}'" for source_type in type_map)

    url = _utils_pg_url()
    literal = "'" + url.replace("'", "''") + "'"
    target = "'" + str(parquet_path).replace("'", "''") + "'"
    con = duckdb.connect()
    try:
        con.execute('INSTALL postgres; LOAD postgres;')
        con.execute(f'ATTACH {literal} AS up (TYPE postgres, READ_ONLY)')
        con.execute(
            f"""
            COPY (
              SELECT
                ({case_sql})::UINTEGER AS key_identifier_type_id,
                source_id AS key_value,
                {inchikey_type_id}::UINTEGER AS canonical_identifier_type_id,
                inchikey AS canonical_identifier
              FROM up.omnipath_utils.resolver_chemical
              WHERE source_type IN ({in_list})
                AND source_id IS NOT NULL
                AND inchikey IS NOT NULL
            ) TO {target} (FORMAT PARQUET)
            """
        )
        rows = con.execute(
            f'SELECT count(*) FROM read_parquet({target})'
        ).fetchone()[0]
    finally:
        con.close()
    return {'utils_chemical_rows': int(rows)}
